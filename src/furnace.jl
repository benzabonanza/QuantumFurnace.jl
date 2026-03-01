function construct_lindbladian(jumps::Vector{JumpOp}, config::Config{Lindbladian}, hamiltonian::HamHam;
    trotter::Union{TrottTrott, Nothing}=nothing)

    domain_name = replace(string(typeof(config.domain)), "Domain" => "")
    println("Constructing Liouvillian ($(domain_name))")

    ham_or_trott = if config.domain isa TrotterDomain
        trotter === nothing && error("A Trotter object must be provided for the TrotterDomain")
        trotter
    else # For Bohr, Energy, Time domains
        hamiltonian
    end

    precomputed_data = _precompute_data(config, ham_or_trott)

    dim = size(hamiltonian.data, 1)
    T = eltype(hamiltonian.eigvals)
    CT = Complex{T}
    total_lindbladian = zeros(CT, dim^2, dim^2)

    # Build Workspace{Lindbladian} with LiouvillianScratch
    sc = LiouvillianScratch(CT, dim)
    Id = Matrix{CT}(I, dim, dim)
    ws = Workspace{Lindbladian, typeof(config.domain), typeof(config.construction), T}(
        nothing, nothing, nothing, nothing,  # physics data
        nothing, nothing, nothing, nothing,  # G fields
        nothing, nothing, nothing, nothing, nothing,  # channel fields
        nothing, nothing, nothing, nothing,  # domain precomputed (transition, gnf, energy_labels, oft_domain_prefactor)
        nothing, nothing, nothing, nothing, nothing, nothing, nothing,  # domain-specific (oft_nufft_prefactors, bohr_alpha, bohr_keys, bohr_is, bohr_js, b_minus, b_plus)
        nothing,  # coherent_unitaries
        nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing,  # trajectory fields
        Id,       # Id
        sc,       # scratch
    )

    # Precompute all B's for the A's if for KMS DB and with_coherent.
    Btot = _precompute_coherent_B(jumps, ham_or_trott, config, precomputed_data)
    if Btot !== nothing
        _vectorize_liouvillian_coherent!(total_lindbladian, Btot, ws)
    end

    # Jumps arrive in the correct basis (trotter.eigvecs for TrotterDomain,
    # hamiltonian.eigvecs for other domains -- basis selection is at the source).

    # Accumulate Liouvillian in-place (no per-jump dim^2 x dim^2 allocations).
    for (k, jump) in pairs(jumps)
        _jump_contribution!(total_lindbladian, jump, ham_or_trott, config, precomputed_data, ws;
            coherent_term=nothing)
    end

    return total_lindbladian
end

# ============================================================================
# Public entry points
# ============================================================================

"""
    run_lindblad(jumps, config, hamiltonian, trotter=nothing) -> LindbladResults

Dense Liouvillian spectral analysis via Arpack shift-invert.

Constructs the full Lindbladian superoperator, finds the two eigenvalues nearest zero
(steady state and gap mode), and returns spectral data in a `LindbladResults` struct.

# Arguments
- `jumps::Vector{JumpOp}`: Jump operators
- `config::Config{Lindbladian}`: Lindbladian configuration
- `hamiltonian::HamHam`: Hamiltonian with eigenbasis data
- `trotter::Union{TrottTrott, Nothing}=nothing`: Trotter object (required for TrotterDomain)

# Returns
`LindbladResults` with eigenvalues, fixed point, gap mode, spectral gap, and metadata.
"""
function run_lindblad(
    jumps::Vector{JumpOp},
    config::Config{Lindbladian,D,C,T},
    hamiltonian::HamHam{T},
    trotter::Union{TrottTrott, Nothing}=nothing,
) where {D, C, T<:AbstractFloat}

    t_start = time()

    validate_config!(config)
    _print_press(config)

    liouv = construct_lindbladian(jumps, config, hamiltonian, trotter=trotter)
    @printf("Done.\n")

    # Arpack shift-invert eigensolver
    shift = 1e-9 * (1 + 1im)
    eigvals_near_zero, eigvecs_near_zero = eigs(liouv, nev=2, sigma=shift, tol=1e-12)
    sorted_permutation_eigen = sortperm(abs.(real.(eigvals_near_zero)))

    ss_index = sorted_permutation_eigen[1]
    gap_index = sorted_permutation_eigen[2]
    spectral_gap = eigvals_near_zero[gap_index]

    steady_state_vec = eigvecs_near_zero[:, ss_index]
    steady_state_dm = reshape(steady_state_vec, size(hamiltonian.data))
    hermitianize!(steady_state_dm)
    steady_state_dm ./= tr(steady_state_dm)

    gap_vec = eigvecs_near_zero[:, gap_index]
    gap_mode_op = reshape(gap_vec, size(hamiltonian.data))

    wall_time = time() - t_start
    metadata = _capture_metadata(wall_time_seconds=wall_time)

    return LindbladResults{T}(
        config,
        Complex{T}.(eigvals_near_zero[sorted_permutation_eigen]),
        Complex{T}.(steady_state_dm),
        Complex{T}.(gap_mode_op),
        Complex{T}(spectral_gap),
        metadata,
    )
end

"""
    run_thermalize(jumps, config, hamiltonian, trotter=nothing; initial_dm=nothing, rng, rescale_by_inv_prob, save_every) -> ThermalizeResults

Density-matrix Kraus evolution toward the Gibbs state.

Evolves an initial density matrix via random jump channels, recording trace distance
to the Gibbs state at configurable intervals. Returns the final state and convergence
history in a `ThermalizeResults` struct.

# Arguments
- `jumps::Vector{JumpOp}`: Jump operators
- `config::Config{Thermalize}`: Thermalization configuration (provides mixing_time, delta)
- `hamiltonian::HamHam`: Hamiltonian with eigenbasis data
- `trotter::Union{TrottTrott, Nothing}=nothing`: Trotter object (required for TrotterDomain)

# Keyword Arguments
- `initial_dm::Union{Nothing, Matrix{<:Complex}}=nothing`: Initial density matrix (defaults to maximally mixed I/d)
- `rng::AbstractRNG=Random.default_rng()`: Random number generator
- `rescale_by_inv_prob::Bool=true`: Rescale delta by 1/p_jump for physical mixing time
- `save_every::Int=1`: Record trace distance every `save_every` steps. Default 1 preserves per-step recording. Convergence cutoff is only checked at save points.

# Returns
`ThermalizeResults` with final density matrix, trace distances, time steps, and metadata.
"""
function run_thermalize(
    jumps::Vector{JumpOp},
    config::Config{Thermalize,D,C,T},
    hamiltonian::HamHam{T},
    trotter::Union{TrottTrott, Nothing}=nothing;
    initial_dm::Union{Nothing, Matrix{<:Complex}}=nothing,
    rng::AbstractRNG = Random.default_rng(),
    rescale_by_inv_prob::Bool = true,
    save_every::Int = 1,
) where {D, C, T<:AbstractFloat}

    dim = size(hamiltonian.data, 1)

    # Default initial_dm: maximally mixed state I/d
    evolving_dm = if initial_dm === nothing
        Matrix{Complex{T}}(I(dim) / dim)
    else
        copy(initial_dm)
    end

    t_start = time()

    validate_config!(config)
    @assert save_every >= 1 "save_every must be >= 1"
    _print_press(config)

    if config.domain isa TrotterDomain
        @assert trotter !== nothing
        ham_or_trott = trotter
        gibbs = Hermitian(trotter.eigvecs' * hamiltonian.eigvecs * hamiltonian.gibbs *
                          hamiltonian.eigvecs' * trotter.eigvecs)
    else
        ham_or_trott = hamiltonian
        gibbs = hamiltonian.gibbs
    end

    precomputed_data = _precompute_data(config, ham_or_trott)

    p_jump = 1.0 / length(jumps)
    coherent_unitaries = _precompute_coherent_unitary(jumps, hamiltonian, config, precomputed_data;
        trotter=trotter, delta_scale = rescale_by_inv_prob ? (1.0 / p_jump) : 1.0)

    CT = eltype(evolving_dm)
    scratch = ThermalizeScratch(CT, dim)

    # Precompute per-jump CPTP channels (K0, U_residual) -- eliminates eigen() from hot loop
    (; K0s, U_residuals) = _precompute_per_jump_channels(
        jumps, ham_or_trott, config, precomputed_data;
        rescale_by_inv_prob = rescale_by_inv_prob,
    )

    # Precompute jump_weight_scaling for _accumulate_rho_jump!
    jump_weight_scaling = rescale_by_inv_prob ? (precomputed_data.gamma_norm_factor / p_jump) : precomputed_data.gamma_norm_factor

    num_steps = Int(ceil(config.mixing_time / config.delta))

    convergence_cutoff = 1e-5
    trace_distances = [trace_distance_h(Hermitian(evolving_dm), gibbs)]
    recorded_steps = Int[0]

    for step in 1:num_steps
        idx = rand(rng, 1:length(jumps))
        jump = jumps[idx]

        # Coherent evolution (already precomputed, same as before)
        _apply_coherent_unitary!(
            evolving_dm,
            coherent_unitaries === nothing ? nothing : coherent_unitaries[idx],
            scratch,
        )

        # Accumulate rho_jump (rho-dependent omega loop, no R/eigen)
        _accumulate_rho_jump!(
            scratch, evolving_dm, jump, ham_or_trott, config, precomputed_data;
            jump_weight_scaling = jump_weight_scaling,
        )

        # Apply precomputed channel (no eigendecomposition)
        _apply_precomputed_channel!(
            evolving_dm, K0s[idx], U_residuals[idx], scratch,
        )

        if step % save_every == 0
            dist = trace_distance_h(Hermitian(evolving_dm), gibbs)
            push!(trace_distances, dist)
            push!(recorded_steps, step)
            @printf("Dist to Gibbs: %s\n", dist)
            if dist < convergence_cutoff
                break
            end
        end
    end

    time_steps = T.(recorded_steps .* config.delta)

    wall_time = time() - t_start
    metadata = _capture_metadata(wall_time_seconds=wall_time)
    metadata[:save_every] = save_every

    return ThermalizeResults{T}(
        config,
        Complex{T}.(evolving_dm),
        T.(trace_distances),
        T.(time_steps),
        metadata,
    )
end
