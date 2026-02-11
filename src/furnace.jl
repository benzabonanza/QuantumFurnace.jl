function run_lindbladian(jumps::Vector{JumpOp}, config::AbstractLiouvConfig, hamiltonian::HamHam; 
    trotter::Union{TrottTrott, Nothing}=nothing)

    validate_config!(config)
    print_press(config)

    liouv = construct_lindbladian(jumps, config, hamiltonian, trotter=trotter)
    @printf("Done.\n")

    # Arpack eigs
    # eigvals_near_zero, eigvecs_near_zero = eigs(liouv, nev=2, which=:SM, tol=1e-12)
    shift = 1e-9 * (1 + 1im)
    eigvals_near_zero, eigvecs_near_zero = eigs(liouv, nev=2, sigma=shift, tol=1e-12)
    sorted_permutation_eigen = sortperm(abs.(real.(eigvals_near_zero)))
    
    ss_index = sorted_permutation_eigen[1]   # Smallest
    gap_index = sorted_permutation_eigen[2]  # Second smallest
    spectral_gap = eigvals_near_zero[gap_index] # Spectral gap

    steady_state_vec = eigvecs_near_zero[:, ss_index]
    steady_state_dm = reshape(steady_state_vec, size(hamiltonian.data))
    steady_state_dm = (steady_state_dm + steady_state_dm') / 2
    steady_state_dm ./= tr(steady_state_dm) # Normalize

    # Gap mode used for LSI
    gap_vec = eigvecs_near_zero[:, gap_index]
    gap_mode_op = reshape(gap_vec, size(hamiltonian.data))

    result = HotSpectralResults(
        data = liouv,
        fixed_point = steady_state_dm,
        gap_mode = gap_mode_op,
        spectral_gap = spectral_gap,
        hamiltonian = hamiltonian,
        trotter = trotter, 
        config = config
    )
    return result
end

function construct_lindbladian(jumps::Vector{JumpOp}, config::AbstractLiouvConfig, hamiltonian::HamHam;
    trotter::Union{TrottTrott, Nothing}=nothing)
    
    domain_name = replace(string(typeof(config.domain)), "Domain" => "")
    println("Constructing Liouvillian ($(domain_name))")

    ham_or_trott = if config.domain isa TrotterDomain
        trotter === nothing && error("A Trotter object must be provided for the TrotterDomain")
        trotter
    else # For Bohr, Energy, Time domains
        hamiltonian
    end

    precomputed_data = precompute_data(config.domain, config, ham_or_trott)

    #! uncomment for multi-threads
    # total_lindbladian = @distributed (+) for jump in jumps
    #     jump_contribution(config.domain, jump, ham_or_trott, config, precomputed_data)
    # end

    dim = size(hamiltonian.data, 1)
    total_lindbladian = zeros(ComplexF64, dim^2, dim^2)
    ws = LindbladianWorkspace(dim)

    # Precompute all B's for the A's if for KMS DB and with_coherent.
    Btot = precompute_coherent_total_B(jumps, hamiltonian, config, precomputed_data; trotter=trotter)
    if Btot !== nothing
        vectorize_liouvillian_coherent!(total_lindbladian, Btot, ws)
    end

    # Accumulate Liouvillian in-place (no per-jump dim^2×dim^2 allocations).
    for (k, jump) in pairs(jumps)
        jump_contribution!(total_lindbladian, config.domain, jump, ham_or_trott, config, precomputed_data, ws; 
            coherent_term=nothing)
    end

    return total_lindbladian
end

function run_thermalization(
    jumps::Vector{JumpOp},
    config::AbstractThermalizeConfig,
    evolving_dm::Matrix{ComplexF64},
    hamiltonian::HamHam;
    trotter::Union{TrottTrott, Nothing}=nothing,
    rng::AbstractRNG = Random.default_rng(),
    rescale_by_inv_prob::Bool = false,
    )

    dim = size(hamiltonian.data, 1)
    validate_config!(config)
    print_press(config)

    if config.domain isa TrotterDomain
        @assert trotter !== nothing
        ham_or_trott = trotter
        gibbs = Hermitian(trotter.eigvecs' * hamiltonian.eigvecs * hamiltonian.gibbs *
                          hamiltonian.eigvecs' * trotter.eigvecs)
    else
        ham_or_trott = hamiltonian
        gibbs = hamiltonian.gibbs
    end

    precomputed_data = precompute_data(config.domain, config, ham_or_trott)

    # precompute coherent U_B = exp(-i delta B(jump)) per jump to avoid allocations
    p_jump = 1.0 / length(jumps)
    coherent_unitaries = precompute_coherent_unitary_terms(jumps, hamiltonian, config, precomputed_data; 
        trotter=trotter, delta_scale = rescale_by_inv_prob ? (1.0 / p_jump) : 1.0)

    scratch = KrausScratch(ComplexF64, dim)

    num_steps = Int(ceil(config.mixing_time / config.delta))

    convergence_cutoff = 1e-5
    distances_to_gibbs = [trace_distance_h(Hermitian(evolving_dm), gibbs)]

    for step in 1:num_steps
        idx = rand(rng, 1:length(jumps))
        jump = jumps[idx]

        jump_contribution!(config.domain,
            evolving_dm,
            jump,
            ham_or_trott,
            config,
            precomputed_data,
            scratch;
            coherent_unitary_cache = (coherent_unitaries === nothing ? nothing : coherent_unitaries[idx]),
            jump_prob = p_jump,
            rescale_by_inv_prob = rescale_by_inv_prob,
            )

        dist = trace_distance_h(Hermitian(evolving_dm), gibbs)
        push!(distances_to_gibbs, dist)
        @printf("Dist to Gibbs: %s\n", dist)
        if dist < convergence_cutoff
            num_steps = step
            break
        end
    end

    time_steps = [0.0:config.delta:(num_steps * config.delta);]
    return HotAlgorithmResults(evolving_dm, distances_to_gibbs, time_steps, hamiltonian, trotter, config)
end

