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
- `B_total::Union{Nothing, Matrix{T}}`: precomputed coherent B (nothing for GNS / with_coherent(construction)=false)
- `jumps::Vector{JumpOp}`: reference to jump operators (kept for external access)
- `jump_eigenbases::Vector{Matrix{T}}`: concrete-typed eigenbasis matrices (avoids JumpOp abstract field boxing)
- `jump_hermitian::Vector{Bool}`: hermitian flags for each jump operator

Scratch matrices (all dim x dim, zeroed or overwritten each matvec call):
- `jump_oft::Matrix{T}`: A(omega) buffer (written by `oft!`)
- `tmp1::Matrix{T}`: scratch for `mul!` results
- `tmp2::Matrix{T}`: scratch for `mul!` results
- `LdagL::Matrix{T}`: scratch for L'*L product
- `rho_out::Matrix{T}`: output accumulator (zeroed at start of each matvec)

Channel fields (populated only for Config{Thermalize}, nothing for Config{Lindbladian}):
- `channel_K0::Union{Nothing, Matrix{T}}`: I - alpha * R_total (Chen Eq. 3.2)
- `channel_U_residual::Union{Nothing, Matrix{T}}`: sqrt_psd(S) residual TP fix
- `channel_U_coherent::Union{Nothing, Matrix{T}}`: exp(-i*delta*B_total) coherent unitary
- `channel_rho_jump::Union{Nothing, Matrix{T}}`: scratch for jump sandwich accumulation
- `channel_delta::Union{Nothing, Float64}`: delta from Config{Thermalize}

Precomputed effective Hamiltonian for optimized Lindbladian matvec (Phase 32):
- `G_left::Union{Nothing, Matrix{T}}`: i*B^T - 0.5*R_total^T (left action in L(rho) = G_left*rho + rho*G_right + sandwiches)
- `G_right::Union{Nothing, Matrix{T}}`: -i*B^T - 0.5*R_total^T (right action)
- `G_left_adj::Union{Nothing, Matrix{T}}`: adjoint left action (= G_right for Hermitian R_total, differs for BohrDomain)
- `G_right_adj::Union{Nothing, Matrix{T}}`: adjoint right action (= G_left for Hermitian R_total, differs for BohrDomain)
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

    # Channel fields (populated only for Config{Thermalize} constructor)
    channel_K0::Union{Nothing, Matrix{T}}
    channel_U_residual::Union{Nothing, Matrix{T}}
    channel_U_coherent::Union{Nothing, Matrix{T}}
    channel_rho_jump::Union{Nothing, Matrix{T}}
    channel_delta::Union{Nothing, Float64}

    # Precomputed effective Hamiltonian for optimized Lindbladian matvec (Phase 32)
    # G_left = i*B^T - 0.5*R_total, G_right = -i*B^T - 0.5*R_total
    # Used directly in gemm!('N','N',...) for zero-allocation hot path.
    G_left::Union{Nothing, Matrix{T}}
    G_right::Union{Nothing, Matrix{T}}
    G_left_adj::Union{Nothing, Matrix{T}}
    G_right_adj::Union{Nothing, Matrix{T}}
end

"""
    KrylovWorkspace(config::Config{Lindbladian}, hamiltonian, jumps; trotter=nothing)

Construct a `KrylovWorkspace` pre-allocating all scratch matrices for the given
(config, hamiltonian) pair. Mirrors `construct_lindbladian` setup in `furnace.jl`.

# Arguments
- `config::Config{Lindbladian}`: Lindbladian configuration (EnergyDomain, TimeDomain, etc.)
- `hamiltonian::HamHam`: Hamiltonian with eigenbasis data
- `jumps::Vector{JumpOp}`: Jump operators (stored by reference)
- `trotter::Union{TrottTrott, Nothing}=nothing`: Trotter object (required for TrotterDomain)
"""
function KrylovWorkspace(
    config::Config{Lindbladian},
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

    # Precompute G_left/G_right for optimized Lindbladian matvec (Phase 32)
    # R_total = sum_i scalar_i * L_i'L_i (physics convention, reuse existing helper)
    R_total = zeros(CT, dim, dim)
    _accumulate_R_total!(R_total, jump_eigenbases, jump_hermitian,
                         precomputed_data, config, ham_or_trott)
    hermitianize!(R_total)

    # Build G_left/G_right (stored as the actual matrices used in gemm!('N','N',...))
    # Convention (L*rho*L'): L(rho) = G_left*rho + rho*G_right + sandwiches
    # Non-sandwich part: i*B^T*rho - i*rho*B^T  (from _vectorize_liouvillian_coherent!)
    #                  - 0.5*R*rho - 0.5*rho*R   (new anticommutator with kron(I,R) and kron(R^T,I))
    # => G_left = i*B^T - 0.5*R,  G_right = -i*B^T - 0.5*R
    if B_total !== nothing
        B_T = Matrix{CT}(transpose(B_total))
        G_left  = 1im .* B_T .- 0.5 .* R_total
        G_right = -1im .* B_T .- 0.5 .* R_total
    else
        G_left  = -0.5 .* R_total
        G_right = -0.5 .* R_total
    end
    G_left  = Matrix{CT}(G_left)
    G_right = Matrix{CT}(G_right)

    # Adjoint: HS adjoint of rho -> M*rho is rho -> M'*rho, of rho -> rho*N is rho -> rho*N'.
    # G_left_adj  = G_left'  = (i*B^T - 0.5*R)' = -i*B^T - 0.5*R = G_right  (B Hermitian, R Hermitian after hermitianize!)
    # G_right_adj = G_right' = (-i*B^T - 0.5*R)' = i*B^T - 0.5*R = G_left
    # This holds for ALL domains because R_total is Hermitianized before G construction.
    G_left_adj  = G_right
    G_right_adj = G_left

    return KrylovWorkspace{CT, typeof(precomputed_data)}(
        precomputed_data, B_total, jumps,
        jump_eigenbases, jump_hermitian,
        jump_oft, tmp1, tmp2, LdagL, rho_out,
        nothing, nothing, nothing, nothing, nothing,  # channel fields
        G_left, G_right, G_left_adj, G_right_adj,     # Phase 32 precomputed effective Hamiltonian
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
    config::Config{<:Any, EnergyDomain},
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
    config::Config{<:Any, D},
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
    config::Config{<:Any, BohrDomain},
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
# Config{Thermalize} workspace constructor
# ---------------------------------------------------------------------------

"""
    KrylovWorkspace(config::Config{Thermalize}, hamiltonian, jumps; trotter=nothing)

Construct a `KrylovWorkspace` with precomputed CPTP channel matrices for
the faithful Chen channel (Eq. 3.2).

Precomputes R_total, K0, U_residual, U_coherent at construction time so the
per-matvec cost is only the rho-dependent sandwich terms.

# Arguments
- `config::Config{Thermalize}`: Thermalization configuration (provides delta)
- `hamiltonian::HamHam`: Hamiltonian with eigenbasis data
- `jumps::Vector{JumpOp}`: Jump operators (stored by reference)
- `trotter::Union{TrottTrott, Nothing}=nothing`: Trotter object (required for TrotterDomain)
"""
function KrylovWorkspace(
    config::Config{Thermalize},
    hamiltonian::HamHam,
    jumps::Vector{JumpOp};
    trotter::Union{TrottTrott, Nothing}=nothing,
)
    # Determine ham_or_trott
    ham_or_trott = if config.domain isa TrotterDomain
        trotter === nothing && error("A Trotter object must be provided for the TrotterDomain")
        trotter
    else
        hamiltonian
    end

    # Precompute domain-specific data (Config accepts both Lindbladian and Thermalize)
    precomputed_data = _precompute_data(config, ham_or_trott)

    # Precompute coherent B_total (returns nothing for GNS / with_coherent=false)
    B_total = _precompute_coherent_total_B(jumps, ham_or_trott, config, precomputed_data)

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
                         precomputed_data, config, ham_or_trott)
    hermitianize!(R_total)

    # Precompute G_left/G_right for optimized Lindbladian matvec (Phase 32)
    # Convention (L*rho*L'): G_left = i*B^T - 0.5*R,  G_right = -i*B^T - 0.5*R
    if B_total !== nothing
        B_T = Matrix{CT}(transpose(B_total))
        G_left  = 1im .* B_T .- 0.5 .* R_total
        G_right = -1im .* B_T .- 0.5 .* R_total
    else
        G_left  = -0.5 .* R_total
        G_right = -0.5 .* R_total
    end
    G_left  = Matrix{CT}(G_left)
    G_right = Matrix{CT}(G_right)

    # Adjoint G matrices: G_left_adj = G_left' = G_right, G_right_adj = G_right' = G_left
    # (holds for all domains since R_total is Hermitianized)
    G_left_adj  = G_right
    G_right_adj = G_left

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
        G_left, G_right, G_left_adj, G_right_adj,
    )
end
