"""
    KrylovWorkspace{T, PD}

Pre-allocated workspace for matrix-free Lindbladian action (matvec).

Stores precomputed data (transition function, energy labels, gamma normalization),
the coherent correction matrix B (if applicable), jump operator references,
and scratch matrices for zero-allocation dissipator accumulation.

Tied to a specific (config, hamiltonian) pair at construction time.

# Type Parameters
- `T <: Complex`: element type of all matrices (e.g. `ComplexF64`)
- `PD <: NamedTuple`: concrete type of the precomputed data tuple (varies by domain)

# Fields
- `precomputed_data::PD`: from `_precompute_data(config, ham_or_trott)`
- `B_total::Union{Nothing, Matrix{T}}`: precomputed coherent B (nothing for GNS / with_coherent=false)
- `jumps::Vector{JumpOp}`: reference to jump operators (kept for external access)
- `jump_eigenbases::Vector{Matrix{T}}`: concrete-typed eigenbasis matrices (avoids JumpOp abstract field boxing)
- `jump_hermitian::Vector{Bool}`: hermitian flags for each jump operator

Scratch matrices (all dim x dim, zeroed or overwritten each matvec call):
- `jump_oft::Matrix{T}`: A(omega) buffer (written by `oft!`)
- `tmp1::Matrix{T}`: scratch for `mul!` results
- `tmp2::Matrix{T}`: scratch for `mul!` results
- `LdagL::Matrix{T}`: scratch for L'*L product
- `rho_out::Matrix{T}`: output accumulator (zeroed at start of each matvec)

Channel fields (populated only for ThermalizeConfig, nothing for LiouvConfig):
- `channel_K0::Union{Nothing, Matrix{T}}`: I - alpha * R_total (Chen Eq. 3.2)
- `channel_U_residual::Union{Nothing, Matrix{T}}`: sqrt_psd(S) residual TP fix
- `channel_U_coherent::Union{Nothing, Matrix{T}}`: exp(-i*delta*B_total) coherent unitary
- `channel_rho_jump::Union{Nothing, Matrix{T}}`: scratch for jump sandwich accumulation
- `channel_delta::Union{Nothing, Float64}`: delta from ThermalizeConfig
"""
struct KrylovWorkspace{T<:Complex, PD<:NamedTuple}
    # Precomputed data (immutable after construction)
    precomputed_data::PD
    B_total::Union{Nothing, Matrix{T}}
    jumps::Vector{JumpOp}

    # Concrete-typed jump data for zero-allocation hot path
    # (JumpOp.in_eigenbasis is Matrix{<:Complex} -- abstract element type causes boxing)
    jump_eigenbases::Vector{Matrix{T}}
    jump_hermitian::Vector{Bool}

    # Scratch matrices for dissipator accumulation (dim x dim)
    jump_oft::Matrix{T}
    tmp1::Matrix{T}
    tmp2::Matrix{T}
    LdagL::Matrix{T}
    rho_out::Matrix{T}

    # Channel fields (populated only for ThermalizeConfig constructor)
    channel_K0::Union{Nothing, Matrix{T}}
    channel_U_residual::Union{Nothing, Matrix{T}}
    channel_U_coherent::Union{Nothing, Matrix{T}}
    channel_rho_jump::Union{Nothing, Matrix{T}}
    channel_delta::Union{Nothing, Float64}
end

"""
    KrylovWorkspace(config, hamiltonian, jumps; trotter=nothing)

Construct a `KrylovWorkspace` pre-allocating all scratch matrices for the given
(config, hamiltonian) pair. Mirrors `construct_lindbladian` setup in `furnace.jl`.

# Arguments
- `config::AbstractLiouvConfig`: Lindbladian configuration (EnergyDomain, TimeDomain, etc.)
- `hamiltonian::HamHam`: Hamiltonian with eigenbasis data
- `jumps::Vector{JumpOp}`: Jump operators (stored by reference)
- `trotter::Union{TrottTrott, Nothing}=nothing`: Trotter object (required for TrotterDomain)
"""
function KrylovWorkspace(
    config::AbstractLiouvConfig,
    hamiltonian::HamHam,
    jumps::Vector{JumpOp};
    trotter::Union{TrottTrott, Nothing}=nothing,
)
    # Determine ham_or_trott (mirrors construct_lindbladian in furnace.jl)
    ham_or_trott = if config.domain isa TrotterDomain
        trotter === nothing && error("A Trotter object must be provided for the TrotterDomain")
        trotter
    else
        hamiltonian
    end

    # Precompute domain-specific data (transition, energy_labels, gamma_norm_factor, ...)
    precomputed_data = _precompute_data(config, ham_or_trott)

    # Precompute coherent B_total (returns nothing for GNS / with_coherent=false)
    B_total = _precompute_coherent_total_B(jumps, ham_or_trott, config, precomputed_data)

    # Determine dimensions and element type
    dim = size(hamiltonian.data, 1)
    CT = Complex{eltype(hamiltonian.eigvals)}

    # Extract concrete-typed eigenbasis matrices and hermitian flags
    # (JumpOp.in_eigenbasis is Matrix{<:Complex} -- abstract type parameter
    #  causes boxing allocations in the hot path)
    jump_eigenbases = [Matrix{CT}(j.in_eigenbasis) for j in jumps]
    jump_hermitian  = [j.hermitian for j in jumps]

    # Allocate scratch matrices
    jump_oft = zeros(CT, dim, dim)
    tmp1     = zeros(CT, dim, dim)
    tmp2     = zeros(CT, dim, dim)
    LdagL    = zeros(CT, dim, dim)
    rho_out  = zeros(CT, dim, dim)

    return KrylovWorkspace{CT, typeof(precomputed_data)}(
        precomputed_data, B_total, jumps,
        jump_eigenbases, jump_hermitian,
        jump_oft, tmp1, tmp2, LdagL, rho_out,
        nothing, nothing, nothing, nothing, nothing,  # channel fields
    )
end

# ---------------------------------------------------------------------------
# ThermalizeConfig -> LiouvConfig conversion (moved from krylov_eigsolve.jl)
# ---------------------------------------------------------------------------

"""
    _thermalize_to_liouv_config(tc::ThermalizeConfig) -> LiouvConfig

Build a `LiouvConfig` from a `ThermalizeConfig` by copying all shared Lindbladian
parameters and stripping `mixing_time` and `delta`.

Required because `apply_lindbladian!` dispatches on `AbstractLiouvConfig`, not
`AbstractThermalizeConfig`.
"""
function _thermalize_to_liouv_config(tc::ThermalizeConfig)
    LiouvConfig(
        num_qubits = tc.num_qubits,
        with_coherent = tc.with_coherent,
        with_linear_combination = tc.with_linear_combination,
        domain = tc.domain,
        beta = tc.beta,
        sigma = tc.sigma,
        gaussian_parameters = tc.gaussian_parameters,
        a = tc.a,
        b = tc.b,
        num_energy_bits = tc.num_energy_bits,
        t0 = tc.t0,
        w0 = tc.w0,
        eta = tc.eta,
        num_trotter_steps_per_t0 = tc.num_trotter_steps_per_t0,
    )
end

"""
    _thermalize_to_liouv_config(tc::ThermalizeConfigGNS) -> LiouvConfigGNS

Build a `LiouvConfigGNS` from a `ThermalizeConfigGNS`. GNS configs always have
`with_coherent = false`.
"""
function _thermalize_to_liouv_config(tc::ThermalizeConfigGNS)
    LiouvConfigGNS(
        num_qubits = tc.num_qubits,
        with_coherent = false,
        with_linear_combination = tc.with_linear_combination,
        domain = tc.domain,
        beta = tc.beta,
        sigma = tc.sigma,
        gaussian_parameters = tc.gaussian_parameters,
        a = tc.a,
        b = tc.b,
        num_energy_bits = tc.num_energy_bits,
        t0 = tc.t0,
        w0 = tc.w0,
        eta = tc.eta,
        num_trotter_steps_per_t0 = tc.num_trotter_steps_per_t0,
    )
end

# ---------------------------------------------------------------------------
# R_total accumulation helpers (physics convention: R = sum rate^2 * L' * L)
# ---------------------------------------------------------------------------

"""
    _accumulate_R_total!(R, ws, config, hamiltonian) -> nothing

Accumulate R_total = sum over all jumps and frequencies of rate^2 * (L' * L)
in physics convention. Used at workspace construction time (not per-matvec).

Matches the R accumulation in `_jump_contribution!` for EnergyDomain
(jump_workers.jl) but summed over all jumps.
"""
function _accumulate_R_total!(
    R::Matrix{T},
    ws_eigenbases::Vector{Matrix{T}},
    ws_hermitian::Vector{Bool},
    precomputed_data,
    config::AbstractLiouvConfig{EnergyDomain},
    hamiltonian::HamHam,
) where {T<:Complex}
    (; transition, gamma_norm_factor, energy_labels) = precomputed_data
    bohr_freqs = hamiltonian.bohr_freqs
    inv_4sigma2 = 1.0 / (4 * config.sigma^2)
    prefactor = (config.w0 / (config.sigma * sqrt(2 * pi))) * gamma_norm_factor

    dim = size(R, 1)
    jump_oft = zeros(T, dim, dim)
    LdagL = zeros(T, dim, dim)

    for (k, eigenbasis) in enumerate(ws_eigenbases)
        is_herm = ws_hermitian[k]
        if is_herm
            for w_raw in energy_labels
                w_raw > 1e-12 && continue
                w = abs(w_raw)
                @. jump_oft = eigenbasis * exp(-(w - bohr_freqs)^2 * inv_4sigma2)
                rate2 = prefactor * transition(w)
                # R += rate^2 * (L' * L)  [physics convention]
                mul!(LdagL, jump_oft', jump_oft)
                @. R += rate2 * LdagL
                if w > 1e-12
                    rate2_neg = prefactor * transition(-w)
                    # Negative freq: L_neg = L', so L_neg'*L_neg = L*L'
                    mul!(LdagL, jump_oft, jump_oft')
                    @. R += rate2_neg * LdagL
                end
            end
        else
            for w in energy_labels
                @. jump_oft = eigenbasis * exp(-(w - bohr_freqs)^2 * inv_4sigma2)
                rate2 = prefactor * transition(w)
                mul!(LdagL, jump_oft', jump_oft)
                @. R += rate2 * LdagL
            end
        end
    end
    return nothing
end

function _accumulate_R_total!(
    R::Matrix{T},
    ws_eigenbases::Vector{Matrix{T}},
    ws_hermitian::Vector{Bool},
    precomputed_data,
    config::AbstractLiouvConfig{D},
    ham_or_trott::Union{HamHam, TrottTrott},
) where {T<:Complex, D<:Union{TimeDomain, TrotterDomain}}
    (; transition, gamma_norm_factor, energy_labels, oft_nufft_prefactors) = precomputed_data
    prefactor = config.w0 * config.t0^2 * (config.sigma * sqrt(2 / pi)) / (2 * pi) * gamma_norm_factor

    dim = size(R, 1)
    jump_oft = zeros(T, dim, dim)
    LdagL = zeros(T, dim, dim)

    for (k, eigenbasis) in enumerate(ws_eigenbases)
        is_herm = ws_hermitian[k]
        if is_herm
            for w_raw in energy_labels
                w_raw > 1e-12 && continue
                w = abs(w_raw)
                nufft_pf = _prefactor_view(oft_nufft_prefactors, w)
                @. jump_oft = eigenbasis * nufft_pf
                rate2 = prefactor * transition(w)
                mul!(LdagL, jump_oft', jump_oft)
                @. R += rate2 * LdagL
                if w > 1e-12
                    rate2_neg = prefactor * transition(-w)
                    mul!(LdagL, jump_oft, jump_oft')
                    @. R += rate2_neg * LdagL
                end
            end
        else
            for w in energy_labels
                nufft_pf = _prefactor_view(oft_nufft_prefactors, w)
                @. jump_oft = eigenbasis * nufft_pf
                rate2 = prefactor * transition(w)
                mul!(LdagL, jump_oft', jump_oft)
                @. R += rate2 * LdagL
            end
        end
    end
    return nothing
end

function _accumulate_R_total!(
    R::Matrix{T},
    ws_eigenbases::Vector{Matrix{T}},
    ws_hermitian::Vector{Bool},
    precomputed_data,
    config::AbstractLiouvConfig{BohrDomain},
    hamiltonian::HamHam,
) where {T<:Complex}
    (; alpha, gamma_norm_factor) = precomputed_data
    dim = size(R, 1)
    jump_oft = zeros(T, dim, dim)
    A_nu2_dag = zeros(T, dim, dim)

    for (k, eigenbasis) in enumerate(ws_eigenbases)
        for nu_2 in keys(hamiltonian.bohr_dict)
            # alpha_A = B_nu2
            @. jump_oft = alpha(hamiltonian.bohr_freqs, nu_2) * eigenbasis

            # Build A_nu2_dag via index scatter
            fill!(A_nu2_dag, 0)
            indices = hamiltonian.bohr_dict[nu_2]
            @inbounds for idx in indices
                i = idx[1]; j = idx[2]
                A_nu2_dag[j, i] = conj(eigenbasis[i, j])
            end

            # R += gamma_norm_factor * (A_nu2_dag * alpha_A)
            # This matches the R accumulation in _jump_contribution! for BohrDomain
            mul!(R, A_nu2_dag, jump_oft, gamma_norm_factor, 1.0)
        end
    end
    return nothing
end

# ---------------------------------------------------------------------------
# ThermalizeConfig workspace constructor
# ---------------------------------------------------------------------------

"""
    KrylovWorkspace(config::AbstractThermalizeConfig, hamiltonian, jumps; trotter=nothing)

Construct a `KrylovWorkspace` with precomputed CPTP channel matrices for
the faithful Chen channel (Eq. 3.2).

Precomputes R_total, K0, U_residual, U_coherent at construction time so the
per-matvec cost is only the rho-dependent sandwich terms.

# Arguments
- `config::AbstractThermalizeConfig`: Thermalization configuration (provides delta)
- `hamiltonian::HamHam`: Hamiltonian with eigenbasis data
- `jumps::Vector{JumpOp}`: Jump operators (stored by reference)
- `trotter::Union{TrottTrott, Nothing}=nothing`: Trotter object (required for TrotterDomain)
"""
function KrylovWorkspace(
    config::AbstractThermalizeConfig,
    hamiltonian::HamHam,
    jumps::Vector{JumpOp};
    trotter::Union{TrottTrott, Nothing}=nothing,
)
    # Convert to LiouvConfig for precomputation dispatch
    config_liouv = _thermalize_to_liouv_config(config)

    # Determine ham_or_trott
    ham_or_trott = if config.domain isa TrotterDomain
        trotter === nothing && error("A Trotter object must be provided for the TrotterDomain")
        trotter
    else
        hamiltonian
    end

    # Precompute domain-specific data (uses LiouvConfig for dispatch)
    precomputed_data = _precompute_data(config_liouv, ham_or_trott)

    # Precompute coherent B_total (returns nothing for GNS / with_coherent=false)
    B_total = _precompute_coherent_total_B(jumps, ham_or_trott, config_liouv, precomputed_data)

    # Determine dimensions and element type
    dim = size(hamiltonian.data, 1)
    CT = Complex{eltype(hamiltonian.eigvals)}

    # Extract concrete-typed eigenbasis matrices and hermitian flags
    jump_eigenbases = [Matrix{CT}(j.in_eigenbasis) for j in jumps]
    jump_hermitian  = [j.hermitian for j in jumps]

    # Allocate scratch matrices
    jump_oft = zeros(CT, dim, dim)
    tmp1     = zeros(CT, dim, dim)
    tmp2     = zeros(CT, dim, dim)
    LdagL    = zeros(CT, dim, dim)
    rho_out  = zeros(CT, dim, dim)

    # --- Precompute channel matrices ---
    delta = config.delta

    # 1. Compute R_total (physics convention)
    R_total = zeros(CT, dim, dim)
    _accumulate_R_total!(R_total, jump_eigenbases, jump_hermitian,
                         precomputed_data, config_liouv, ham_or_trott)
    hermitianize!(R_total)

    # 2. Compute K0, S, U_residual (Chen Eq. 3.2)
    alpha_chen = 1 - sqrt(1 - delta)
    K0 = Matrix{CT}(I, dim, dim) .- alpha_chen .* R_total

    # S = (2*alpha - delta)*R - alpha^2 * R^2
    R2 = R_total * R_total
    S = (2 * alpha_chen - delta) .* R_total .- (alpha_chen^2) .* R2
    hermitianize!(S)

    # PSD guard: clamp negative eigenvalues to zero
    eig = eigen(Hermitian(S))
    eig.values .= max.(eig.values, 0.0)
    U_residual = Matrix{CT}(Diagonal(sqrt.(eig.values)) * eig.vectors')

    # 3. Compute U_coherent = exp(-i*delta*B_total) if coherent
    U_coherent = if B_total !== nothing
        Matrix{CT}(exp(-1im * delta * Hermitian(B_total)))
    else
        nothing
    end

    # 4. Allocate channel scratch
    channel_rho_jump = zeros(CT, dim, dim)

    return KrylovWorkspace{CT, typeof(precomputed_data)}(
        precomputed_data, B_total, jumps,
        jump_eigenbases, jump_hermitian,
        jump_oft, tmp1, tmp2, LdagL, rho_out,
        K0, U_residual, U_coherent, channel_rho_jump, Float64(delta),
    )
end
