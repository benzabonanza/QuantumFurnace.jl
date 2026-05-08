# Single-jump B_bohr variant removed in Phase 35; callers use [jump] wrapper.
function B_bohr(hamiltonian::HamHam{T}, jumps::AbstractVector{<:JumpOp}, config::Config{<:Any, <:Any, KMS}) where {T<:AbstractFloat}

    dim = size(hamiltonian.data, 1)
    CT = Complex{T}
    unique_freqs = collect(keys(hamiltonian.bohr_dict))

    f = _pick_f(config)  # Picks rates for B in Bohr domain

    n_jumps = length(jumps)
    n_freqs = length(unique_freqs)
    if Threads.nthreads() > 1 && n_jumps * n_freqs >= OMEGA_THREAD_THRESHOLD
        return _B_bohr_threaded(hamiltonian, jumps, f, unique_freqs, dim, n_jumps, CT)
    end

    B = zeros(CT, dim, dim)
    for jump in jumps
        f_A_nu_1 = zeros(CT, dim, dim)
        for nu_2 in unique_freqs
            indices = hamiltonian.bohr_dict[nu_2]
            @. f_A_nu_1 = f(hamiltonian.bohr_freqs, nu_2) * jump.in_eigenbasis
            # B += A_nu_2' * f_A_nu_1 expanded per-index:
            # A_nu_2'[j,i] = conj(jump.in_eigenbasis[i,j]) for (i,j) in indices
            @inbounds for idx in indices
                i, j = idx[1], idx[2]
                val = conj(jump.in_eigenbasis[i, j])
                @inbounds for col in 1:dim
                    B[j, col] += val * f_A_nu_1[i, col]
                end
            end
        end
    end
    return B
end

# qf-6af.5: thread the (jump_idx, freq_idx) outer×inner double loop. Each
# task accumulates a private dim×dim B partial and uses its own f_A_nu_1
# scratch; reduce after `@sync`. Used by Workspace setups that hit B_bohr
# from `_precompute_coherent_B` (BohrDomain + EnergyDomain Lindbladian KMS).
function _B_bohr_threaded(
    hamiltonian::HamHam{T},
    jumps::AbstractVector{<:JumpOp},
    f,
    unique_freqs::Vector,
    dim::Int,
    n_jumps::Int,
    ::Type{CT},
) where {T<:AbstractFloat, CT}
    n_freqs = length(unique_freqs)
    n_work  = n_jumps * n_freqs
    nt = min(Threads.nthreads(), n_work)
    chunks = _partition_range(1:n_work, nt)
    n_chunks = length(chunks)

    B_partials  = [zeros(CT, dim, dim) for _ in 1:n_chunks]
    f_A_buffers = [Matrix{CT}(undef, dim, dim) for _ in 1:n_chunks]

    old_blas = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    try
        @sync for (idx, chunk) in enumerate(chunks)
            Threads.@spawn _B_bohr_chunk!(
                B_partials[idx], f_A_buffers[idx],
                hamiltonian, jumps, f, unique_freqs,
                dim, n_freqs, chunk)
        end
    finally
        BLAS.set_num_threads(old_blas)
    end

    B = zeros(CT, dim, dim)
    @inbounds for idx in 1:n_chunks
        B .+= B_partials[idx]
    end
    return B
end

function _B_bohr_chunk!(
    B_partial::Matrix{CT},
    f_A_nu_1::Matrix{CT},
    hamiltonian::HamHam,
    jumps::AbstractVector{<:JumpOp},
    f,
    unique_freqs::Vector,
    dim::Int,
    n_freqs::Int,
    chunk::UnitRange{Int},
) where {CT}
    @inbounds for w_idx in chunk
        # Linear (jump_idx, freq_idx) decoding: outer jump, inner freq.
        jump_idx = ((w_idx - 1) ÷ n_freqs) + 1
        freq_idx = ((w_idx - 1) % n_freqs) + 1
        jump = jumps[jump_idx]
        nu_2 = unique_freqs[freq_idx]
        indices = hamiltonian.bohr_dict[nu_2]

        @. f_A_nu_1 = f(hamiltonian.bohr_freqs, nu_2) * jump.in_eigenbasis
        @inbounds for idx in indices
            i, j = idx[1], idx[2]
            val = conj(jump.in_eigenbasis[i, j])
            @inbounds for col in 1:dim
                B_partial[j, col] += val * f_A_nu_1[i, col]
            end
        end
    end
    return nothing
end

function _pick_f(config::Config{<:Any, <:Any, KMS})

    beta = config.beta
    sigma = config.sigma
    if config.with_linear_combination
        a = config.a
        s = config.s
        return (nu_1, nu_2) -> create_f(nu_1, nu_2, beta, sigma, a, s)
    else
        gaussian_parameters = config.gaussian_parameters
        return (nu_1, nu_2) -> create_f_gauss(nu_1, nu_2, beta, sigma, gaussian_parameters)
    end
end

function create_f(nu_1::Real, nu_2::Real, beta::Real, sigma::Real, a::Real, s::Real)
    alpha = create_alpha(nu_1, nu_2, beta, sigma, a, s)
    return tanh(-beta * (nu_1 - nu_2) / 4) * alpha / (2im)
end

function create_f_gauss(nu_1::Real, nu_2::Real, beta::Real, sigma::Real,
    gaussian_parameters::Union{Tuple{<:Real, <:Real}, Tuple{Nothing, Nothing}})
    """Tanh * alpha."""
    alpha_nu1_nu2 = create_alpha_gauss(nu_1, nu_2, sigma, gaussian_parameters)
    return tanh(-beta * (nu_1 - nu_2) / 4) * alpha_nu1_nu2 / (2im)
end

_pick_alpha(config::Config{<:Any, <:Any, KMS}) = _pick_alpha_kms(config)
_pick_alpha(config::Config{<:Any, <:Any, GNS}) = _pick_alpha_gns(config)

# 2-arg forms: compute alpha directly via dispatch (zero allocation on hot path)
function _pick_alpha(config::Config{<:Any, <:Any, KMS}, nu_1::Real, nu_2::Real)
    if config.with_linear_combination
        return create_alpha(nu_1, nu_2, config.beta, config.sigma, config.a, config.s)
    else
        return create_alpha_gauss(nu_1, nu_2, config.sigma, config.gaussian_parameters)
    end
end

function _pick_alpha(config::Config{<:Any, <:Any, GNS}, nu_1::Real, nu_2::Real)
    if config.with_linear_combination
        return create_alpha_gns(nu_1, nu_2, config.beta, config.sigma, config.a, config.s)
    else
        return create_alpha_gauss(nu_1, nu_2, config.sigma, config.gaussian_parameters)
    end
end

function _pick_alpha_kms(config::Config{<:Any, <:Any, KMS})

    sigma = config.sigma
    if config.with_linear_combination
        beta = config.beta
        a = config.a
        s = config.s
        return (nu_1, nu_2) -> create_alpha(nu_1, nu_2, beta, sigma, a, s)
    else
        gaussian_parameters = config.gaussian_parameters
        return (nu_1, nu_2) -> create_alpha_gauss(nu_1, nu_2, sigma, gaussian_parameters)
    end
end

function create_alpha(nu_1::Real, nu_2::Real, beta::Real, sigma::Real, a::Real, s::Real)

    sqrtA = sqrt(beta * (4 * a + 1) / 4)
    sqrtB = sqrt(beta / 16) * abs(nu_1 + nu_2)
    C = beta * (nu_1 + nu_2) / 4
    prefactor = exp(a * beta^2 * sigma^2 / 2) / 2
    u_min = sqrt(beta * sigma^2 * (1 + s) / 2)
    z_plus = sqrtA * u_min + sqrtB / u_min
    z_minus = sqrtA * u_min - sqrtB / u_min

    alpha_nu_1 = (prefactor * exp(-C) * exp(-(nu_1 - nu_2)^2 / (8 * sigma^2)) * exp(- 2 * sqrtA * sqrtB) *
                    (erfc(z_minus) + exp(4 * sqrtA * sqrtB) * erfc(z_plus)))

    return alpha_nu_1
end

function _pick_alpha_gns(config::Config{<:Any, <:Any, GNS})

    sigma = config.sigma
    if config.with_linear_combination
        beta = config.beta
        a = config.a
        s = config.s
        return (nu_1, nu_2) -> create_alpha_gns(nu_1, nu_2, beta, sigma, a, s)
    else
        gaussian_parameters = config.gaussian_parameters  # but now β = 2ω_γ / σ_γ^2
        return (nu_1, nu_2) -> create_alpha_gauss(nu_1, nu_2, sigma, gaussian_parameters)
    end
end

"""
Coming from unshifted γ(ω) in the energy domain, leads to a partially shifted Kossakowski matrix α.
Difference vs the KMS DB case: |ν1 + ν2| → |ν1 + ν2 + β σ^2 / 2| (and thus not skew symmetric due to the Gaussian filters)
Also: No f and no b-funcs will be constructed because in this setup there is no fine-tuned B in the Lindbladian.
"""
function create_alpha_gns(nu_1::Real, nu_2::Real, beta::Real, sigma::Real, a::Real, s::Real)
    sqrtA = sqrt(beta * (4 * a + 1) / 4)
    sqrtB = sqrt(beta / 16) * abs(nu_1 + nu_2 + beta * sigma^2 / 2)
    C = beta * (nu_1 + nu_2) / 4
    prefactor = exp(a * beta^2 * sigma^2 / 2) / 2
    u_min = sqrt(beta * sigma^2 * (1 + s) / 2)
    z_plus = sqrtA * u_min + sqrtB / u_min
    z_minus = sqrtA * u_min - sqrtB / u_min

    alpha_nu_1 = (prefactor * exp(-C) * exp(-(nu_1 - nu_2)^2 / (8 * sigma^2)) * exp(- 2 * sqrtA * sqrtB) *
                    (erfc(z_minus) + exp(4 * sqrtA * sqrtB) * erfc(z_plus)))

    return alpha_nu_1
end

function create_alpha_gauss(
    nu_1::Real,
    nu_2::Real,
    sigma::Real,
    gaussian_parameters::Union{Tuple{<:Real, <:Real}, Tuple{Nothing, Nothing}})

    (w_gamma, sigma_gamma) = gaussian_parameters
    combined_sigma = sigma^2 + sigma_gamma^2
    prefactor = sigma_gamma / sqrt(combined_sigma)
    alpha_fn(nu_1) = prefactor * (exp(-(nu_1 + nu_2 + 2 * w_gamma)^2 / (8 * combined_sigma))
                                    * exp(-(nu_1 - nu_2)^2 / (8 * sigma^2)))
    return alpha_fn(nu_1)
end

#* TOOLS --------------------------------------------------------------------------------------------------------------------
function create_bohr_dict(bohr_freqs::Matrix{T}) where {T<:AbstractFloat}
    """Creates a dictionary, where the keys are the Bohr frequencies, and the values are a list of their sparse indices
    in the Bohr matrix. (With special care on the diagonal elements, that are identically 0.)"""

    bohr_dict = DefaultDict{T, Vector{CartesianIndex{2}}}(() -> CartesianIndex{2}[])
    dim = size(bohr_freqs, 1)
    bohr_dict[zero(T)] = CartesianIndex{2}.(1:dim, 1:dim) # nu = 0.0 is the diagonal and might be other offdiags
    for j in 1:dim
        for i in 1:(j - 1)
            push!(bohr_dict[bohr_freqs[i, j]], CartesianIndex{2}(i, j))
            push!(bohr_dict[-bohr_freqs[i, j]], CartesianIndex{2}(j, i))
        end
    end
    return bohr_dict
end

function check_alpha_skew_symmetry(alpha::Function, nu_1::Real, nu_2::Real, beta::Real)
    @assert norm(alpha(nu_1, nu_2) - alpha(-nu_2, -nu_1) * exp(-beta * (nu_1 + nu_2) / 2)) < 1e-14
end
