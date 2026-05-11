# Single-jump B_bohr variant removed in Phase 35; callers use [jump] wrapper.
#
# qf-sta: outer-ν₂ / inner-jump restructure with sparse row-aware f-eval.
# The previous loop computed `f(bohr_freqs, ν₂) * in_eb` over the full
# d×d matrix for every (jump, ν₂), wasting work on rows that never appear
# in `bohr_dict[ν₂]`. The new loop walks `bohr_dict[ν₂]` directly and
# evaluates `f(bohr_freqs[i, :], ν₂)` once per (ν₂, i), reusing across
# all jumps. For non-degenerate spectra |bohr_dict[ν]|≈1, so this matches
# the optimal `O(d³)` f-evaluation count (vs. `O(n_jumps · n_freqs · d²)`
# in the old loop). Measured 720× speedup at n=7 (48.6s → 67ms) with
# zero precomputed sparsity (per-call allocation stays at 16·nt·d²
# bytes — 540 MB at n=11, vs ~64 GB if a full per-ν CSR layout were
# precomputed).
function B_bohr(hamiltonian::HamHam{T}, jumps::AbstractVector{<:JumpOp}, config::Config{<:Any, <:Any, KMS}) where {T<:AbstractFloat}

    dim = size(hamiltonian.data, 1)
    CT = Complex{T}
    unique_freqs = collect(keys(hamiltonian.bohr_dict))

    f = _pick_f(config)  # Picks rates for B in Bohr domain

    n_jumps = length(jumps)
    n_freqs = length(unique_freqs)
    bohr_freqs = hamiltonian.bohr_freqs
    # qf-qmi.1: hoist concrete-typed views once so the hot loop is type-stable.
    in_ebs = [(jump.in_eigenbasis::Matrix{CT}) for jump in jumps]

    if Threads.nthreads() > 1 && n_freqs >= OMEGA_THREAD_THRESHOLD
        return _B_bohr_threaded(hamiltonian, in_ebs, f, unique_freqs, bohr_freqs, dim, n_jumps, CT)
    end

    B = zeros(CT, dim, dim)
    f_row = Vector{CT}(undef, dim)

    for nu_2 in unique_freqs
        indices = hamiltonian.bohr_dict[nu_2]
        last_i = 0
        @inbounds for idx in indices
            i = idx[1]; j = idx[2]
            if i != last_i
                for col in 1:dim
                    f_row[col] = f(bohr_freqs[i, col], nu_2)
                end
                last_i = i
            end
            for jump_idx in 1:n_jumps
                in_eb = in_ebs[jump_idx]
                val = conj(in_eb[i, j])
                @inbounds for col in 1:dim
                    B[j, col] += val * f_row[col] * in_eb[i, col]
                end
            end
        end
    end
    return B
end

# qf-6af.5 / qf-sta: thread the outer ν₂ loop. Each task accumulates a
# private dim×dim B partial and uses its own f_row scratch; reduce after
# `@sync`. Used by Workspace setups that hit B_bohr from
# `_precompute_coherent_B` (BohrDomain + EnergyDomain Lindbladian KMS).
function _B_bohr_threaded(
    hamiltonian::HamHam{T},
    in_ebs::Vector{Matrix{CT}},
    f,
    unique_freqs::Vector,
    bohr_freqs::AbstractMatrix{<:Real},
    dim::Int,
    n_jumps::Int,
    ::Type{CT},
) where {T<:AbstractFloat, CT}
    n_freqs = length(unique_freqs)
    nt = min(Threads.nthreads(), n_freqs)
    chunks = _partition_range(1:n_freqs, nt)
    n_chunks = length(chunks)

    B_partials = [zeros(CT, dim, dim) for _ in 1:n_chunks]

    old_blas = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    try
        @sync for (cidx, chunk) in enumerate(chunks)
            Threads.@spawn _B_bohr_chunk!(
                B_partials[cidx], hamiltonian, in_ebs, f, unique_freqs,
                bohr_freqs, dim, n_jumps, chunk, CT)
        end
    finally
        BLAS.set_num_threads(old_blas)
    end

    B = zeros(CT, dim, dim)
    @inbounds for cidx in 1:n_chunks
        B .+= B_partials[cidx]
    end
    return B
end

function _B_bohr_chunk!(
    B_partial::Matrix{CT},
    hamiltonian::HamHam,
    in_ebs::Vector{Matrix{CT}},
    f,
    unique_freqs::Vector,
    bohr_freqs::AbstractMatrix{<:Real},
    dim::Int,
    n_jumps::Int,
    chunk::UnitRange{Int},
    ::Type{CT},
) where {CT}
    f_row = Vector{CT}(undef, dim)
    @inbounds for k in chunk
        nu_2 = unique_freqs[k]
        indices = hamiltonian.bohr_dict[nu_2]
        last_i = 0
        for idx in indices
            i = idx[1]; j = idx[2]
            if i != last_i
                for col in 1:dim
                    f_row[col] = f(bohr_freqs[i, col], nu_2)
                end
                last_i = i
            end
            for jump_idx in 1:n_jumps
                in_eb = in_ebs[jump_idx]
                val = conj(in_eb[i, j])
                @inbounds for col in 1:dim
                    B_partial[j, col] += val * f_row[col] * in_eb[i, col]
                end
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

"""
    default_smooth_s(beta, sigma) -> Real

Default smooth-Metropolis regularisation `s` that preserves a *constant
absolute smoothing width* `σ·√s = SMOOTH_S_REF_WIDTH = 0.05` across
β-sweeps along the σ = c/β line (c = O(1)) — i.e. the kinky-Metropolis
γ(ω) is mollified by the same Δω in absolute units, independent of β.

Domain assumption: callers use σ in 1/β units (σ = c/β with c ∈ {0.25, …, 3}
typical). The formula uses σ directly, so it remains well-defined off this
line, but the kink-depth-relative interpretation only holds on it. β is in
the signature both for the dispatch contract and to flag the domain to
readers (and to enable a sanity assertion later if needed).

Calibration: at β = 10, σ = 1/β = 0.1, thesis-numerics default `s = 0.25`,
giving σ·√s = 0.05. Holding that absolute width fixed,

    s_default(σ) = (SMOOTH_S_REF_WIDTH / σ)²    = β²/400   at σ = 1/β

| β  |  σ = 1/β  |  s_default |
|----|-----------|------------|
|  5 | 0.20      |  0.0625    |
| 10 | 0.10      |  0.25      |  ← calibration point (qf-3il)
| 20 | 0.05      |  1.0       |
| 50 | 0.02      |  6.25      |

Rationale: the σ-Gaussian filter f_σ narrows linearly with σ along σ = c/β,
so the same fractional smoothing of the kinky-Metro kink at ω = −βσ²/2
requires σ·√s ∝ const. Larger β (narrower σ) needs *more* `s` to keep the
energy-kernel quadrature error — and hence the energy-register size r_D —
uniform across the β-sweep. The smoothing parameter affects only γ in the
dissipator; the leading b_+ quadrature error is independent of s.

This is opt-in: callers pass `s = default_smooth_s(β, σ)` to Config
explicitly. Existing fixtures locked at `s = 0.25` (β = 10 regression tests)
stay numerically identical because the formula evaluates to 0.25 at the
calibration point.
"""
const SMOOTH_S_REF_WIDTH = 0.05  # σ·√s calibration: (β=10, σ=0.1, s=0.25)
default_smooth_s(beta::Real, sigma::Real) = (SMOOTH_S_REF_WIDTH / sigma)^2

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
