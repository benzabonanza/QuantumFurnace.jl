"""
    oft_domain_prefactor(::EnergyDomain, w0, sigma)
    oft_domain_prefactor(::TimeDomain, w0, sigma, t0)
    oft_domain_prefactor(::TrotterDomain, w0, sigma, t0)

Compute the domain-dependent scalar prefactor for OFT-based rate calculations.

EnergyDomain:  w0 / (sigma * sqrt(2*pi))
Time/Trotter:  w0 * t0^2 * sigma * sqrt(2/pi) / (2*pi)

No method for BohrDomain -- callers use `gamma_norm_factor` directly.
"""
oft_domain_prefactor(::EnergyDomain, w0::Real, sigma::Real) = w0 / (sigma * sqrt(2 * pi))
oft_domain_prefactor(::TimeDomain, w0::Real, sigma::Real, t0::Real) = w0 * t0^2 * (sigma * sqrt(2 / pi)) / (2 * pi)
oft_domain_prefactor(::TrotterDomain, w0::Real, sigma::Real, t0::Real) = w0 * t0^2 * (sigma * sqrt(2 / pi)) / (2 * pi)

function _precompute_labels(config::Config{<:Any, D}) where {D<:Union{BohrDomain, EnergyDomain}}
    # Dissipative register only (qf-9z0): EnergyDomain has no t-grid.
    energy_labels = _create_energy_labels(register_r_D(config), register_w0_D(config))
    truncated_energy_labels = _truncate_energy_labels(energy_labels, config)
    return (truncated_energy_labels,)  # Energy labels
end

function _precompute_labels(config::Config{<:Any, D}) where {D<:Union{TimeDomain, TrotterDomain}}
    # Dissipative register only (qf-9z0). Coherent labels are built per-term
    # in `_precompute_data` from the `b_minus / b_plus` registers.
    r_D = register_r_D(config)
    w0_D = register_w0_D(config)
    t0_D = register_t0_D(config)
    energy_labels = _create_energy_labels(r_D, w0_D)
    truncated_energy_labels = _truncate_energy_labels(energy_labels, config)
    time_labels = energy_labels .* (t0_D / w0_D)
    return (truncated_energy_labels, time_labels) # Energy and time labels (dissipative grid)
end

function _precompute_data(
    config::Config{Lindbladian, BohrDomain},
    ham_or_trott::Union{HamHam, TrottTrott}
)

    alpha = _pick_alpha(config)
    # Grid-independent normalisation: `pick_gamma_sup(config)` is the
    # closed-form continuum sup of γ — 1.0 for every standard family.
    # Replaces the prior `1.0 / maximum(transition.(energy_labels))`
    # which sampled the sup on a discrete `(r_D, w0_D)`-dependent grid.
    gamma_norm_factor = 1.0 / pick_gamma_sup(config)
    return (
        alpha = alpha,
        gamma_norm_factor = gamma_norm_factor
    )
end

function _precompute_data(
    config::Config{Thermalize, BohrDomain},
    hamiltonian::HamHam
)
    alpha = _pick_alpha(config)

    # Grid-independent normalisation (qf-etx); see Config{Lindbladian, BohrDomain} branch above.
    gamma_norm_factor = 1.0 / pick_gamma_sup(config)

    # Cache the Bohr buckets as plain Int index pairs to avoid CartesianIndex overhead
    # and avoid rebuilding any per-frequency index lists inside jump_contribution!.
    bohr_keys = collect(keys(hamiltonian.bohr_dict))
    bohr_is = Vector{Vector{Int}}(undef, length(bohr_keys))
    bohr_js = Vector{Vector{Int}}(undef, length(bohr_keys))
    @inbounds for (k, nu) in pairs(bohr_keys)
        idxs = hamiltonian.bohr_dict[nu]
        is = Vector{Int}(undef, length(idxs))
        js = Vector{Int}(undef, length(idxs))
        @inbounds for t in eachindex(idxs)
            is[t] = idxs[t][1]
            js[t] = idxs[t][2]
        end
        bohr_is[k] = is
        bohr_js[k] = js
    end

    return (
        alpha = alpha,
        gamma_norm_factor = gamma_norm_factor,
        bohr_keys = bohr_keys,
        bohr_is = bohr_is,
        bohr_js = bohr_js,
    )

end

function _precompute_data(
    config::Config{<:Any, EnergyDomain},
    ham_or_trott::Union{HamHam, TrottTrott}
)
    energy_labels, = _precompute_labels(config)
    transition = pick_transition(config)
    # Grid-independent normalisation (qf-etx); see Config{Lindbladian, BohrDomain} branch.
    gamma_norm_factor = 1.0 / pick_gamma_sup(config)
    # EnergyDomain dissipator only consults `w0_D` (no time grid).
    dp = oft_domain_prefactor(config.domain, register_w0_D(config), config.sigma)

    return (
        transition = transition,
        gamma_norm_factor = gamma_norm_factor,
        energy_labels = energy_labels,
        oft_domain_prefactor = dp,
    )
end

"""
    _precompute_data(config::Config{Lindbladian, BohrDomain, DLL}, hamiltonian::HamHam)

DLL Bohr-domain precompute (Ding–Li–Lin 2024, Eq. 3.4 first form). DLL has no
γ(ω) transition rates and no γ-norm factor — the filter `f̂(ν)` already encodes
the KMS weight. Returns just the resolved filter.
"""
function _precompute_data(
    config::Config{Lindbladian, BohrDomain, DLL},
    hamiltonian::HamHam,
)
    return (filter = _resolve_filter(config),)
end

"""
    _precompute_data(config::Config{Lindbladian, TimeDomain, DLL}, hamiltonian::HamHam)

DLL Time-domain precompute (Ding–Li–Lin 2024, Eq. 3.4 third form). Builds the
trapezoidal time grid `t_m = m · t0` for `m ∈ [-N/2, N/2 − 1]`, `N =
2^num_energy_bits`, then truncates by `filter_time_cutoff`. Stores the uniform
spacing `t0` as the integration weight together with a single-slice NUFFT
prefactor evaluated at `ω = 0`:

    pf[i, j] = Σ_m time_kernel(filter, t_m) · cis((λ_i − λ_j) · t_m)

This is exactly the DFT in `dll_lindblad_op_time`'s explicit triple loop;
computing it via FINUFFT amortises the per-jump cost from `O(Nt · n²)` to
`O(Nt log Nt + n² log(1/ε))` (FINUFFT precision `ε = 1e-12`). The Riemann sum
itself — and hence the Eq. 3.15 quadrature error structure — is unchanged.

No `gamma_norm_factor`, `transition`, or `w0` — the DLL Lindblad operator is
the OFT at `ω = 0` of the coupling, with the filter providing the only
weighting.
"""
function _precompute_data(
    config::Config{Lindbladian, TimeDomain, DLL},
    hamiltonian::HamHam{T},
) where {T<:AbstractFloat}
    filter = _resolve_filter(config)
    # Build the time grid directly from `(r_D, t0_D)` — the DLL construction
    # has no ω-grid, so `w0_D` plays no role here. DLL has no `b_-/b_+` split
    # either; G shares the same dissipative grid (qf-9z0).
    r_D = register_r_D(config)
    t0_D = register_t0_D(config)
    N = 2^r_D
    raw_time_labels = collect((-N÷2):(N÷2 - 1)) .* t0_D
    oft_time_labels = _truncate_time_labels_for_oft(raw_time_labels, config.sigma; filter=filter)

    # Single-slice NUFFT at ω = 0 per channel; replaces the per-jump explicit
    # `cis()` triple loop in `dll_lindblad_op_time` with a single FINUFFT eval.
    # For single-channel filters this is a length-1 list — the consumer loop
    # in `_jump_contribution!` (DLL TimeDomain) iterates uniformly.
    sub_filters = _filter_channels_for_dll_oft(filter)
    oft_nufft_at_zero_list = Matrix{Complex{T}}[]
    for sub in sub_filters
        nufft = _prepare_oft_nufft_prefactors(
            hamiltonian.bohr_freqs, oft_time_labels, T[zero(T)], sub; eps=1e-12,
        )
        push!(oft_nufft_at_zero_list, Matrix(@view nufft.data[:, :, 1]))
    end

    return (
        filter = filter,
        time_labels = oft_time_labels,
        t0 = t0_D,
        oft_nufft_at_zero_list = oft_nufft_at_zero_list,
    )
end

# Helper: enumerate the per-channel filters used to build the DLL TimeDomain
# OFT prefactor stack. For single-channel DLL filters the list has length 1
# (the filter itself); the `DLLMultiChannelFilter` overload lives in
# `src/dll_multichannel.jl`.
@inline _filter_channels_for_dll_oft(filter::AbstractFilter) = (filter,)

function _precompute_data(
    config::Config{<:Any, D},
    ham_or_trott::Union{HamHam, TrottTrott}
) where {D<:Union{TimeDomain, TrotterDomain}}
    energy_labels, time_labels = _precompute_labels(config)  # dissipative grid (qf-9z0)
    oft_time_labels = _truncate_time_labels_for_oft(time_labels, config.sigma; filter=_resolve_filter(config))

    transition = pick_transition(config)
    # Grid-independent normalisation (qf-etx); see Config{Lindbladian, BohrDomain} branch.
    gamma_norm_factor = 1.0 / pick_gamma_sup(config)

    # Coherent term B (KMS-only). Each leg is built on its own register —
    # `b_-(t)` on the outer triple, `b_+(τ)` on the inner triple. qf-9z0.3
    # plumbs these grids through `B_time / B_trotter`; here we only need the
    # truncated per-leg dictionaries.
    b_minus, b_plus = if with_coherent(config.construction)
        time_labels_b_minus = _create_energy_labels(register_r_b_minus(config),
            register_w0_b_minus(config)) .* (register_t0_b_minus(config) / register_w0_b_minus(config))
        time_labels_b_plus  = _create_energy_labels(register_r_b_plus(config),
            register_w0_b_plus(config))  .* (register_t0_b_plus(config)  / register_w0_b_plus(config))
        _b_minus = _compute_truncated_func(_compute_b_minus, time_labels_b_minus, config.beta, config.sigma)
        chosen_b_plus, b_plus_args = _select_b_plus_calculator(config)
        _b_plus = _compute_truncated_func(chosen_b_plus, time_labels_b_plus, b_plus_args...)
        (_b_minus, _b_plus)
    else
        (nothing, nothing)
    end

    # OFT NUFFT prefactors (dissipative path; same call for TimeDomain/TrotterDomain).
    # `nthreads=1` preserves byte-deterministic FINUFFT output across runs
    # (load-bearing for the byte-identity tests in `test_dll_filter.jl` and
    # the regression checks in `test_regression.jl`). FINUFFT is not the
    # construction-time hot spot per the qf-6af motivation; enabling its
    # OpenMP threading gives a single-digit-percent improvement at most while
    # introducing ULP-level non-determinism that breaks the byte-identity
    # contract above.
    oft_nufft_prefactors = _prepare_oft_nufft_prefactors(
        ham_or_trott.bohr_freqs,
        oft_time_labels,
        energy_labels,
        _resolve_filter(config);
        eps=1e-12,
        nthreads=1,
    )

    dp = oft_domain_prefactor(config.domain, register_w0_D(config), config.sigma, register_t0_D(config))

    return (
        transition = transition,
        gamma_norm_factor = gamma_norm_factor,
        energy_labels = energy_labels,
        oft_nufft_prefactors = oft_nufft_prefactors,
        b_minus = b_minus,
        b_plus = b_plus,
        oft_domain_prefactor = dp,
    )
end

function _select_b_plus_calculator(config::Config{<:Any, <:Any, KMS})
    if !config.with_linear_combination
        # Gaussian
        return (_compute_b_plus, (config.beta, config.gaussian_parameters[1], config.gaussian_parameters[2]))
    else
        if config.a != 0.0
            # a-regularized smooth (covers Metro a>0,s=0 and Glauberish a>0,s>0)
            return (_compute_b_plus_smooth, (config.beta, config.sigma, config.a, config.s))
        else
            # eta-regularized smooth Metro (covers Chen plain s=0 and thesis-main s>0)
            s_val = something(config.s, 0.0)
            return (_compute_b_plus_metro, (config.beta, config.sigma, config.eta, s_val))
        end
    end
end

# ---------------------------------------------------------------------------
# CPTP weak-measurement channel construction (Chen et al. Eq. 3.2)
# ---------------------------------------------------------------------------

"""
    _build_cptp_channel(R, delta) -> (; K0, U_residual, alpha)

Construct the CPTP weak-measurement channel matrices from the accumulated R matrix
(Chen et al. Eq. 3.2).

Returns a NamedTuple with:
- `K0 = I - alpha * R` (no-event Kraus operator)
- `U_residual = sqrt_psd(S)` where `S = (2*alpha - delta)*R - alpha^2 * R^2`
- `alpha = 1 - sqrt(1 - delta)` (the step-size parameter)

The PSD guard clamps negative eigenvalues of S to zero before taking the square root,
handling numerical noise robustly.

# Arguments
- `R::Matrix{<:Complex}`: accumulated R matrix (Hermitianized), either per-operator R^a or summed R_total
- `delta::Real`: step size parameter

# Returns
NamedTuple `(; K0, U_residual, alpha)` with `K0::Matrix` and `U_residual::Matrix`.
"""
function _build_cptp_channel(R::Matrix{T}, delta::Real) where {T<:Complex}
    dim = size(R, 1)
    alpha = 1 - sqrt(1 - delta)

    # K0 = I - alpha * R
    K0 = Matrix{T}(I, dim, dim) .- alpha .* R

    # S = (2*alpha - delta)*R - alpha^2 * R^2
    R2 = R * R
    S = (2 * alpha - delta) .* R .- (alpha^2) .* R2
    hermitianize!(S)

    # PSD guard: clamp negative eigenvalues to zero
    eig = eigen(Hermitian(S))
    eig.values .= max.(eig.values, 0.0)
    U_residual = Matrix{T}(Diagonal(sqrt.(eig.values)) * eig.vectors')

    return (; K0, U_residual, alpha)
end

"""
    _apply_precomputed_channel!(evolving_dm, K0, U_residual, scratch)

Apply precomputed CPTP channel to evolving density matrix. Uses pre-stored K0
and U_residual (no eigendecomposition). `scratch.rho_jump` must be pre-filled
by `_accumulate_rho_jump!`.

Computes: rho_next = K0 * rho * K0' + rho_jump + U_residual * rho * U_residual'
"""
function _apply_precomputed_channel!(
    evolving_dm::Matrix{<:Complex},
    K0::Matrix{<:Complex},
    U_residual::Matrix{<:Complex},
    scratch::ThermalizeScratch{<:Complex},
)
    # rho_next = K0 * rho * K0'
    mul!(scratch.sandwich_tmp, K0, evolving_dm)
    mul!(scratch.rho_next, scratch.sandwich_tmp, K0')

    # + rho_jump (pre-filled by _accumulate_rho_jump!)
    scratch.rho_next .+= scratch.rho_jump

    # + U_residual * rho * U_residual'
    mul!(scratch.sandwich_tmp, U_residual, evolving_dm)
    mul!(scratch.rho_next, scratch.sandwich_tmp, U_residual', 1.0, 1.0)

    # Keep it a density matrix numerically
    hermitianize!(scratch.rho_next)
    copyto!(evolving_dm, scratch.rho_next)
    return evolving_dm
end

"""
    _precompute_per_jump_channels(jumps, ham_or_trott, config, precomputed_data; rescale_by_inv_prob=true)

Precompute per-jump CPTP channel matrices (K0, U_residual) for all jumps.
Eliminates per-step eigendecomposition in the `run_thermalize` hot loop.

For each jump operator, computes R^a via `_precompute_R`, optionally scales by `1/p_jump`
(when `rescale_by_inv_prob=true`), and builds the CPTP channel via `_build_cptp_channel`.

Returns a NamedTuple `(; K0s, U_residuals)` where each is a `Vector{Matrix{CT}}` of length `n_jumps`.
"""
function _precompute_per_jump_channels(
    jumps::AbstractVector{<:JumpOp},
    ham_or_trott::Union{HamHam, TrottTrott},
    config::Config{Thermalize},
    precomputed_data;
    rescale_by_inv_prob::Bool = true,
)
    CT = if ham_or_trott isa HamHam
        Complex{eltype(ham_or_trott.eigvals)}
    else
        eltype(ham_or_trott.eigvecs)  # already Complex{T}
    end
    dim = size(ham_or_trott isa HamHam ? ham_or_trott.data : ham_or_trott.eigvecs, 1)
    n_jumps = length(jumps)
    p_jump = 1.0 / n_jumps

    K0s = Vector{Matrix{CT}}(undef, n_jumps)
    U_residuals = Vector{Matrix{CT}}(undef, n_jumps)

    # Create temporary scratch for R construction (construction-time only)
    builder_scratch = ThermalizeScratch(CT, dim)

    @inbounds for a in 1:n_jumps
        _precompute_R([jumps[a]], ham_or_trott, config, precomputed_data, builder_scratch)
        R_a = copy(builder_scratch.R)
        if rescale_by_inv_prob
            R_a .*= (1.0 / p_jump)
        end
        (; K0, U_residual) = _build_cptp_channel(R_a, config.delta)
        K0s[a] = K0
        U_residuals[a] = U_residual
    end

    return (; K0s, U_residuals)
end

"""
    _apply_one_dm_substep!(evolving_dm, scratch, jump, U_coherent, K0, U_residual,
                           ham_or_trott, config, precomputed_data, jump_weight_scaling)

Apply one per-jump substep `e^{δ𝓛_a}` to the evolving density matrix:
coherent unitary, accumulate ρ_jump, then the precomputed CPTP channel.
Used by `run_thermalize` for both `:sweep` (called S times in order) and
`:random` (called once with a random jump index per outer δ-step).
"""
@inline function _apply_one_dm_substep!(
    evolving_dm::Matrix{<:Complex},
    scratch::ThermalizeScratch{<:Complex},
    jump::JumpOp,
    U_coherent::Union{Nothing, Matrix{<:Complex}},
    K0::Matrix{<:Complex},
    U_residual::Matrix{<:Complex},
    ham_or_trott::Union{HamHam, TrottTrott},
    config::Config{Thermalize},
    precomputed_data,
    jump_weight_scaling::Real,
)
    _apply_coherent_unitary!(evolving_dm, U_coherent, scratch)
    _accumulate_rho_jump!(
        scratch, evolving_dm, jump, ham_or_trott, config, precomputed_data;
        jump_weight_scaling = jump_weight_scaling,
    )
    _apply_precomputed_channel!(evolving_dm, K0, U_residual, scratch)
    return nothing
end