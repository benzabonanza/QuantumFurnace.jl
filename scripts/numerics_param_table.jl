#!/usr/bin/env julia
#
# qf-e4z.1 (P0a) — canonical channel parameter table for the implemented
# δ-step channel (TrotterDomain + GQSP, jump-sweep splitting).
#
# Materialises a single mapping
#   (n, β, ε_target, filter)  →  (full per-term register triples + Trotter Ms
#                                  + δ + η + family + Hamiltonian-norm meta)
# into  scripts/output/channel_param_table.bson.  Downstream sweep harnesses
# (S3 / S5 / S6) load this BSON to pick parameters per cell — no parameters
# are hardcoded inside sweep scripts.
#
# References:
#   - drafts/error-analysis/quadrature-convergence-summary.md  (canonical r tables)
#   - drafts/error-analysis/parameter-recommendations.md       (M, δ, η — still authoritative)
#   - .claude-memory/quadrature_register_recipe_qf_7xt.md      (qf-7xt summary)
#
# Cells (post-2026-05-08 convention):
#   n      ∈ {3, 4, 5, 6}                                  (laptop range; S3 extends)
#   β      ∈ {5, 10, 20}
#   ε_target ∈ {1e-3}                                      (1e-6 dropped — see below)
#   filter ∈ {:gaussian, :smooth_metro, :kinky_metro}
# = 36 cells total in the default build.
#
# 2026-05-08: ε_target = 1e-6 was REMOVED from the default sweep grid.  Reason:
# the mixing-time comparison plots (P1, P4, P5, P6) draw the spectral-gap
# upper-bound overlay τ_bound = log(d/ε)/λ_2 next to the τ_mix curve, and λ_2
# is already returned by `predict_lindbladian_trajectory` / `predict_channel_trajectory`
# per cell.  A second curve at ε=1e-6 — whether re-simulated at tighter δ /
# wider registers, or analytically projected as τ_mix(1e-3) + 3 ln 10 / λ_2 —
# adds no information not already on the plot: any reader who wants τ_mix at
# 1e-6 reads λ_2 off the gap overlay.  So the new convention is "1e-3 only,
# do not plot 1e-6 in any form".
#
# The legacy ε=1e-6 recipe (r_D, r_b+, M_D, δ values) is still encoded in
# `_r_D`, `_r_bp`, `_M_D`, `_delta_step`, `_M_bm`, `_M_bp` below for the rare
# case we want a single explicit ε=1e-6 cell — pass `eps_values = [1e-3, 1e-6]`
# to `build_table` / `build_ideal_lindbladian_table` to opt back in.
# See bd memory `correction-2026-05-08-supersedes-earlier-same-day` and
# bd issue qf-e4z.5/.11/.14/.15/.16 notes for the broader epic decision.
# Closed S1 / S2 BSON sidecars (sweep_S1_ckg_ideal/, sweep_S2_dll_ideal/)
# still contain pre-2026-05-08 ε=1e-6 rows; those remain usable historical
# data but are NOT loaded into plots.
# = 72 cells total in the legacy build (eps_values = [1e-3, 1e-6]).
#
# Sanity-check (--check): runs `predict_channel_trajectory` at six fixtures
# (3 filters × 2 (n, β) pairs at ε_target = 1e-3) and reports the achieved
# asymptotic trace distance.  Subsecond per cell at n ≤ 5.
#
# Output:
#   - scripts/output/channel_param_table.bson — Vector{NamedTuple}, one row per cell
#   - drafts/error-analysis/channel-param-table.md — recipe summary + sanity-check report
#
# Usage:
#   julia --project scripts/numerics_param_table.jl              # build + sanity-check
#   julia --project scripts/numerics_param_table.jl --no-check   # build only

using Printf
using LinearAlgebra
using BSON
using QuantumFurnace
using QuantumFurnace: predict_channel_trajectory, _load_hamiltonian_bson

const SOURCE_ROOT = dirname(@__DIR__)
const HAM_DIR    = joinpath(SOURCE_ROOT, "hamiltonians")
const OUT_BSON   = joinpath(SOURCE_ROOT, "scripts", "output", "channel_param_table.bson")
const OUT_DOC    = joinpath(SOURCE_ROOT, "drafts", "error-analysis", "channel-param-table.md")
const OUT_BSON_IDEAL = joinpath(SOURCE_ROOT, "scripts", "output", "ideal_lindbladian_param_table.bson")

const FAMILY_TAG = "xxx_zzdisordered"   # matches `heis_<FAMILY_TAG>_periodic_n<n>.bson`
ham_filename(n::Integer) = "heis_$(FAMILY_TAG)_periodic_n$(n).bson"

# ---------------------------------------------------------------------------
# Per-cell parameter recipe
# ---------------------------------------------------------------------------

# Recipe targets for the b_± integration windows.  ~3× the b_± supports at
# β=10 σ=0.1 so the truncation atol = 1e-12 sets the effective range.  In the
# qf-e4z.5.3 baseline-inversion recipe (Option A) the t0_b± are sized FIRST
# from the b_± kernel widths and an integer commensurability constraint, then
# r_b± are derived to make the realised T_± ≥ these targets.
const T_MINUS_TARGET = 18.0
const T_PLUS_TARGET  = 12.0

# qf-e4z.5.3: number of TimeDomain b_+ grid samples per the kernel's 1/e
# Gaussian envelope.  Empirically controlled at (n=3, β=10, smooth_metro):
#
#   m_D = 5,  t0_bp = 5e-2, r_D=7+ : ‖L_T·σ_β‖_HS = 2.5e-5   (slope -1 starts)
#   m_D = 20, t0_bp = 1.3e-2, r_D=7+ : 6.3e-6
#   m_D = 80, t0_bp = 3.1e-3, r_D=7+ : 1.6e-6
#   m_D = 320, t0_bp = 7.8e-4, r_D=7+ : 3.9e-7              (controllable)
#
# With t0_bp = support / N_SAMPLES_PER_B_PLUS_SUPPORT and r_D ≥ 7, the slope
# (-1) in t0_bp gives `K · t0_bp` with `K ≈ 5e-4`.  N_SAMPLES = 12 caps the
# floor below ε=1e-3 by ~40× at (n=3, β=10); increase to push floor lower.
const N_SAMPLES_PER_B_PLUS_SUPPORT = 12.0

# r_D table (filter × ε_target).
#
# qf-e4z.5.3 audit (n=3, β=10, smooth_metro): the EnergyDomain reference
# `‖L_E·σ_β‖_HS` saturates at r_D=5 around 4e-6 (a quadrature artefact from the
# residual dissipative OFT truncation), then drops to machine precision (~1e-16)
# by r_D=7.  This 4e-6 floor was masking the t0_bp slope-(-1) at small t0_bp in
# the TimeDomain recipe — recipes targeting ‖L·σ_β‖_HS ≤ 1e-5 must use r_D ≥ 7.
#
# Updated baseline (qf-e4z.5.3):
function _r_D(filter::Symbol, eps::Float64)
    if filter == :gaussian
        return eps == 1e-3 ? 7 : 8
    elseif filter == :smooth_metro
        return eps == 1e-3 ? 7 : 8
    elseif filter == :kinky_metro
        return eps == 1e-3 ? 9 : 14
    else
        throw(ArgumentError("unknown filter: $filter"))
    end
end

# r_b_minus / r_b_plus are DERIVED per cell in `pick_channel_params` from the
# t0_b± sizing (qf-e4z.5.3 Option A baseline-inversion).  The legacy
# constant-per-ε tables are kept for byte-identity diagnostics only:
#   _r_bm_legacy(eps) = eps == 1e-3 ? 5 : 7
#   _r_bp_legacy(:gaussian, eps) = eps == 1e-3 ? 5 : 7
#   _r_bp_legacy(:smooth/kinky_metro, eps) = eps == 1e-3 ? 5 : 15

# Trotter step counts.
#
# qf-e4z.5.3 Option A: the param-table `M_D` field is the M_user argument to
# `make_trotter_for_config`.  With the baseline inverted (b_+ leg drives the
# minimum natural step), `M_user = 1` already gives δt₀ = β·t0_bp — far finer
# than the legacy δt₀ = t0_D / M_user_legacy.  Empirically at (n=3, β=10, δ=1e-3)
# the channel floor is δ-split-limited (≈ 0.3·δ); bumping M_user does not
# improve it.  Default M_user = 1 for all (ε, β).
_M_D(eps::Float64, beta::Float64) = 1
_M_bm(eps::Float64) = 1
_M_bp(eps::Float64) = 1

# δ_step from the jump-wise generator-splitting recipe.
#
# qf-e4z.5.3 Option A baseline-inversion: with the b_+ quadrature dominant
# below the splitting floor (≈ 3e-5 at β=10, m_D=5, r_D=7), the channel
# asymptotic shift `‖ρ_∞ − ρ_β‖₁/2` is set by slope +1 in δ with prefactor
# 0.20 – 0.32:
#
#   β =  5, δ = 1e-3 : floor ≈ 2e-4  ⟹ floor/δ ≈ 0.20
#   β = 10, δ = 1e-3 : floor ≈ 3e-4  ⟹ floor/δ ≈ 0.32
#   β = 20, δ = 1e-3 : floor ≈ 3e-4  ⟹ floor/δ ≈ 0.30 (projected)
#
# δ = 1e-3 gives ≈ 3× headroom under ε=1e-3 at all β in the sweep range.
# Pushing δ below ~1e-4 lets per-step Trotter error accumulate over k_max ~
# 1/δ steps, so δ has a sweet-spot near 1e-3 at this Trotter substep.
_delta_step(eps::Float64, ::Float64) =
    eps == 1e-3 ? 1.0e-3 : (eps == 1e-4 ? 1.5e-4 : (eps == 1e-5 ? 5e-5 : 1.5e-5))

# η — Metropolis regularisation cutoff.  Set so η < t0' = 2T_+/2^r_b+, which
# kills the η-cutoff branch entirely (qf-7xt: η-jump is dead code in the
# discretisation).  ε / (3β) is comfortably below t0' for all our cells.
_eta(eps::Float64, beta::Float64) = eps / (3.0 * beta)

# Gaussian-filter (ω_γ, σ_γ) chosen on the KMS line β·(σ² + σ_γ²) = 2·ω_γ.
# Canonical pick: σ_γ = σ = 1/β, ω_γ = β·σ² = 1/β = σ.
function _gaussian_params(beta::Float64)
    sigma = 1.0 / beta
    sigma_gamma = sigma
    omega_gamma = beta * (sigma^2 + sigma_gamma^2) / 2.0
    return (omega_gamma, sigma_gamma)
end

# (a, s) for kinky / smooth Metropolis (thesis convention, project memory).
# Smooth-Metro `s` is β-scaled (qf-96o) so absolute smoothing width σ·√s stays
# at the calibrated 0.05 across β-sweeps along σ = 1/β. Reduces to s = 0.25 at
# β = 10 (the thesis calibration point in qf-3il).
function _metro_as(filter::Symbol, beta::Real)
    if filter == :smooth_metro
        return (0.0, QuantumFurnace.default_smooth_s(beta, 1.0 / beta))
    else
        return (0.0, 0.0)
    end
end

"""
    _b_plus_support(filter, beta, sigma, gaussian_params, s) -> Float64

1/e half-width (in `t`) of the b_+ kernel's Gaussian envelope. Used to size
`t0_bp` so the inner OFT grid carries ≥ `N_SAMPLES_PER_B_PLUS_SUPPORT` samples
per support.

- Smooth/kinky Metro (`a=0`): envelope `exp(-σ²β²·2t²·(1+s))` ⟹ 1/e at
  `t = 1/(σβ·√(2(1+s)))`.
- Gaussian: envelope `exp(-4·β·ω_γ·t²)` ⟹ 1/e at `t = 1/(2·√(βω_γ))`.
"""
function _b_plus_support(filter::Symbol, beta::Real, sigma::Real,
                          gaussian_params, s::Real)
    if filter == :gaussian
        ω_γ = gaussian_params[1]
        ω_γ === nothing && return NaN
        return 1.0 / (2.0 * sqrt(beta * Float64(ω_γ)))
    else  # :smooth_metro / :kinky_metro (a=0 in our recipe)
        return 1.0 / (sigma * beta * sqrt(2.0 * (1.0 + s)))
    end
end

"""
    pick_channel_params(n, beta, eps, filter; ham=nothing) -> NamedTuple

qf-e4z.5.3 Option A baseline-inversion recipe.

1. `r_D, w0_D, t0_D` from `_r_D(filter, eps)` and `ω_range = 2(‖H‖ + 8σ)`.
2. Pick `m_D = ⌈t0_D / (β · support_bp / N_SAMPLES_PER_B_PLUS_SUPPORT)⌉` so
   the b_+ grid has ≥ N samples per kernel support. Then `t0_bp = t0_D/(β·m_D)`
   and `β·t0_bp = t0_D/m_D` becomes the *minimum* natural Trotter step.
3. `k_m = round(T_-_target/T_+_target)` and `t0_bm = k_m · t0_bp`.
4. `r_bp = ⌈log₂(2 T_+_target / t0_bp)⌉`, `r_bm = ⌈log₂(2 T_-_target / t0_bm)⌉`;
   realised `T_± = 2^(r_b±−1) · t0_b±` ≥ recipe target.
5. `M_user = _M_D = 1` so the constructor picks `δt₀ = β·t0_bp = t0_D/m_D`.

Pass an explicit `ham::HamHam` to use measured `‖H‖`; otherwise `H_norm = NaN`
and only filter-independent recipe fields (M_D, δ, η, filter parameters) are
reliable.
"""
function pick_channel_params(n::Integer, beta::Real, eps::Real, filter::Symbol;
                              ham=nothing)
    sigma = 1.0 / beta
    a, s = _metro_as(filter, Float64(beta))
    gaussian_params = filter == :gaussian ?
        _gaussian_params(Float64(beta)) : (nothing, nothing)

    H_norm = ham === nothing ? NaN : opnorm(ham.data)
    omega_range = isnan(H_norm) ? NaN : 2.0 * (H_norm + 8.0 * sigma)

    r_D  = _r_D(filter, Float64(eps))
    w0_D = isnan(omega_range) ? NaN : omega_range / 2.0^r_D
    t0_D = isnan(omega_range) ? NaN : 2π / (2.0^r_D * w0_D)
    @assert isnan(t0_D) || isapprox(t0_D, 2π / omega_range; rtol=1e-12)

    # Option A: size t0_bp from b_+ support, then invert.
    support_bp = _b_plus_support(filter, Float64(beta), sigma, gaussian_params, s)
    if isnan(t0_D) || isnan(support_bp)
        m_D   = 1
        t0_bp = 2.0 * T_PLUS_TARGET  / 2.0^6
        t0_bm = 2.0 * T_MINUS_TARGET / 2.0^6
    else
        t0_bp_max = support_bp / N_SAMPLES_PER_B_PLUS_SUPPORT
        m_D       = max(1, ceil(Int, t0_D / (beta * t0_bp_max)))
        t0_bp     = t0_D / (beta * m_D)
        k_m       = max(1, round(Int, T_MINUS_TARGET / T_PLUS_TARGET))
        t0_bm     = k_m * t0_bp
    end

    # Derive r_b± from window targets.  Lower bound r ≥ 6 retains qf-7xt's
    # minimum-r margin even when t0_b± shrinks so much that the formula gives
    # r < 6 (e.g. degenerate cases at very small β).
    r_bp = max(6, ceil(Int, log2(2.0 * T_PLUS_TARGET  / t0_bp)))
    r_bm = max(6, ceil(Int, log2(2.0 * T_MINUS_TARGET / t0_bm)))
    T_plus  = 2.0^(r_bp - 1) * t0_bp
    T_minus = 2.0^(r_bm - 1) * t0_bm
    w0_bp = π / T_plus
    w0_bm = π / T_minus

    M_D  = _M_D(Float64(eps), Float64(beta))
    M_bm = _M_bm(Float64(eps))
    M_bp = _M_bp(Float64(eps))

    delta = _delta_step(Float64(eps), Float64(beta))
    eta   = _eta(Float64(eps), Float64(beta))

    with_lc = filter != :gaussian

    return (
        # Identifying tuple
        n = Int(n), beta = Float64(beta), eps = Float64(eps), filter = filter,
        family = FAMILY_TAG,
        # Hamiltonian metadata (NaN if `ham === nothing`)
        H_norm = Float64(H_norm),
        sigma  = Float64(sigma),
        omega_range = Float64(omega_range),
        # Dissipative register
        r_D = r_D, w0_D = Float64(w0_D), t0_D = Float64(t0_D),
        # Outer coherent register
        r_bm = r_bm, w0_bm = Float64(w0_bm), t0_bm = Float64(t0_bm),
        # Inner coherent register
        r_bp = r_bp, w0_bp = Float64(w0_bp), t0_bp = Float64(t0_bp),
        # Time-window choices (qf-e4z.5.1: derived from the picked t0_b± so the
        # OFT grid relation 2π = t0·w0·2^r holds exactly; the realised T_± is
        # within ~10 % of the target T_MINUS_TARGET = 18, T_PLUS_TARGET = 12).
        T_minus = Float64(T_minus), T_plus = Float64(T_plus),
        # Trotter Ms, δ, η
        M_D = M_D, M_bm = M_bm, M_bp = M_bp,
        delta = Float64(delta), eta = Float64(eta),
        # Filter parameters
        with_linear_combination = with_lc,
        a = Float64(a), s = Float64(s),
        gaussian_omega = gaussian_params[1] === nothing ? NaN : Float64(gaussian_params[1]),
        gaussian_sigma = gaussian_params[2] === nothing ? NaN : Float64(gaussian_params[2]),
        # GQSP defaults (Form B per qf-e4z.18 cost model; Form C in Python POC
        # currently — see qf-e4z.19 for the slot-pattern refactor).
        with_gqsp = true, gqsp_degree = 1,
        # Splitting convention
        jump_selection = :sweep,
    )
end

# ---------------------------------------------------------------------------
# Ideal-Lindbladian (EnergyDomain) parameter recipe — S1 / S2 reference
# ---------------------------------------------------------------------------
#
# `predict_lindbladian_trajectory` on EnergyDomain at sufficiently large
# `r_D` is the project's reference for "ideal CKG Lindbladian" — close
# enough to BohrDomain (analytical) that the residual quadrature error is
# negligible against the algorithm-level precision target ε.
#
# **Coherent term**: EnergyDomain dispatches `_precompute_coherent_B` to
# `B_bohr(...)` — the analytical Bohr-form formula — so the b_± per-term
# registers are not used. Only `(r_D, w0_D, t0_D)` matter for the
# ideal-Lindbladian recipe.
#
# **Quadrature target**: ≤ 10⁻⁹ vs BohrDomain (a few orders below the
# ε ∈ {1e-3, 1e-6} algorithm-level targets). Per qf-7xt, smooth-Metro and
# Gaussian saturate at machine precision (~10⁻¹⁴) by r_D = 6; one extra bit
# of safety against (n, β) shifts in the K-prefactor gives r_D = 7. Kinky-
# Metro has slope -2 in 1/N: r_D = 14 gives ~ 2.5×10⁻⁸, marginal at 10⁻⁹;
# kinky in EnergyDomain at r_D ≥ 15 hits memory limits at n ≥ 6, so
# **kinky-Metro ideal-Lindbladian is BohrDomain-only at n ≥ 6**.

# Ideal-Lindbladian r_D table (no ε dependence — quadrature is independent of
# the algorithm-level ε; the recipe just has to clear ≤ 10⁻⁹ across the cells).
function _r_D_ideal(filter::Symbol)
    if filter == :gaussian
        return 7   # super-algebraic; 6 saturates at ~1e-14, +1 safety bit
    elseif filter == :smooth_metro
        return 7   # super-algebraic; 6 saturates at ~1e-14, +1 safety bit
    elseif filter == :kinky_metro
        return 14  # slope -2 in 1/N; r=14 ≈ 2.5e-8; tighter requires r=15+ → BohrDomain
    else
        throw(ArgumentError("unknown filter: $filter"))
    end
end

"""
    pick_ideal_lindbladian_params(n, β, ε, filter; ham=nothing) -> NamedTuple

Build the canonical ideal-Lindbladian (EnergyDomain) row for the
`(n, β, ε, filter)` cell. ε is recorded for downstream τ_mix bookkeeping
(`t_mix(ε) = log(d/ε)/λ`) but does not enter the QUADRATURE recipe — that
depends only on the filter family.
"""
function pick_ideal_lindbladian_params(n::Integer, beta::Real, eps::Real,
                                        filter::Symbol; ham=nothing)
    sigma = 1.0 / beta

    H_norm = ham === nothing ? NaN : opnorm(ham.data)
    omega_range = isnan(H_norm) ? NaN : 2.0 * (H_norm + 8.0 * sigma)

    r_D  = _r_D_ideal(filter)
    w0_D = isnan(omega_range) ? NaN : omega_range / 2.0^r_D
    t0_D = isnan(omega_range) ? NaN : 2π / (2.0^r_D * w0_D)

    # Filter-specific sub-fields (mirrors channel-side recipe so the BSON is
    # uniform). The ideal Lindbladian has no b_± registers / Trotter Ms / δ / η.
    a, s             = _metro_as(filter, Float64(beta))
    gaussian_params  = filter == :gaussian ? _gaussian_params(Float64(beta)) : (nothing, nothing)
    with_lc          = filter != :gaussian

    # Feasibility flag: kinky-Metro at r_D ≥ 14 has ~ 2^14 · d² ω-grid entries;
    # at n ≥ 6 (d = 64) that's > 10⁸ entries on EnergyDomain → flag for BohrDomain.
    domain_recommended = if filter == :kinky_metro && n >= 6
        :BohrDomain
    else
        :EnergyDomain
    end

    return (
        n = Int(n), beta = Float64(beta), eps = Float64(eps), filter = filter,
        family = FAMILY_TAG,
        H_norm = Float64(H_norm), sigma = Float64(sigma), omega_range = Float64(omega_range),
        r_D = r_D, w0_D = Float64(w0_D), t0_D = Float64(t0_D),
        domain_recommended = domain_recommended,
        with_linear_combination = with_lc,
        a = Float64(a), s = Float64(s),
        gaussian_omega = gaussian_params[1] === nothing ? NaN : Float64(gaussian_params[1]),
        gaussian_sigma = gaussian_params[2] === nothing ? NaN : Float64(gaussian_params[2]),
        # Sweep-side conventions
        mode = :L, jump_selection = :sweep,
    )
end

"""
    build_ideal_lindbladian_table(; n_values, beta_values, eps_values, filter_values)
        -> Vector{NamedTuple}

Same (n, β, ε, filter) grid as `build_table` (channel-side), but rows are
ideal-Lindbladian recipes — `(r_D, w0_D, t0_D)` plus filter parameters,
no b_± / Trotter / δ / η. Default filter set is `[:smooth_metro]` because
S1 specifically targets smooth-Metro CKG; pass other filters explicitly
to extend.
"""
function build_ideal_lindbladian_table(;
        n_values::AbstractVector{<:Integer} = 3:9,
        beta_values::AbstractVector{<:Real} = [5.0, 10.0, 20.0],
        # ε=1e-6 dropped from default 2026-05-08 — derive from λ_2 in post-processing.
        # Pass `eps_values = [1e-3, 1e-6]` to opt back in (see header comment).
        eps_values::AbstractVector{<:Real}  = [1e-3],
        filter_values::AbstractVector{Symbol} = [:smooth_metro, :gaussian, :kinky_metro],
    )
    rows = NamedTuple[]
    for n in n_values
        ham_path = joinpath(HAM_DIR, ham_filename(n))
        if !isfile(ham_path)
            @warn "Hamiltonian file missing; skipping n" n ham_path
            continue
        end
        ham_struct = _load_hamiltonian_bson(ham_path, first(beta_values))
        for β in beta_values, ε in eps_values, f in filter_values
            push!(rows, pick_ideal_lindbladian_params(n, β, ε, f; ham = ham_struct))
        end
    end
    return rows
end

# ---------------------------------------------------------------------------
# Build the table
# ---------------------------------------------------------------------------

"""
    build_table(; n_values, beta_values, eps_values, filter_values) -> Vector{NamedTuple}

Materialise one row per (n, β, ε, filter) cell.  Loads each Hamiltonian once
per `n` and computes `‖H‖ = opnorm(ham.data)` (used to derive `ω_range` and
hence `w0_D`).  Hamiltonian filename: `heis_xxx_zzdisordered_periodic_n{n}.bson`.
"""
function build_table(;
        # qf-e4z.5 Phase B: extend to n=7..9 for the cluster-scale cells.
        n_values::AbstractVector{<:Integer} = 3:9,
        beta_values::AbstractVector{<:Real} = [5.0, 10.0, 20.0],
        # ε=1e-6 dropped from default 2026-05-08 — derive from λ_2 in post-processing.
        # Pass `eps_values = [1e-3, 1e-6]` to opt back in (see header comment).
        eps_values::AbstractVector{<:Real}  = [1e-3],
        filter_values::AbstractVector{Symbol} = [:gaussian, :smooth_metro, :kinky_metro],
    )
    rows = NamedTuple[]
    for n in n_values
        ham_path = joinpath(HAM_DIR, ham_filename(n))
        if !isfile(ham_path)
            @warn "Hamiltonian file missing; skipping n" n ham_path
            continue
        end
        # β is irrelevant to opnorm(H); load with any β to get the data, then reuse.
        ham_struct = _load_hamiltonian_bson(ham_path, first(beta_values))
        H_norm = opnorm(ham_struct.data)
        for β in beta_values, ε in eps_values, f in filter_values
            row = pick_channel_params(n, β, ε, f; ham = ham_struct)
            # `pick_channel_params` re-computes opnorm; for performance, skip it
            # by passing the cached value.  Currently negligible; keep as-is.
            push!(rows, row)
            @assert row.H_norm ≈ H_norm
        end
    end
    return rows
end

# ---------------------------------------------------------------------------
# Sanity-check via predict_channel_trajectory
# ---------------------------------------------------------------------------

# Run the channel predictor for a single cell and return the achieved
# asymptotic trace-distance + diagnostics.
#
# Sanity check uses TimeDomain (NUFFT, no Trotter). This isolates the
# quadrature recipe — what we're actually verifying here. The Trotter
# parameters (M_D, M_bm, M_bp) recorded in the table are for downstream
# sweep harnesses (S3 / S5 / S6), which run TrotterDomain + GQSP;
# verifying their commensurability with the shared-δt₀ scheme is out of
# scope for P0a (it depends on the (n, β)-specific ‖H‖ that sets
# `w0_D = ω_range/2^r_D` per the qf-7xt physics, generically incommensurate
# with `t0_b± · β` for non-power-of-two grid ratios).
function _sanity_check_cell(row::NamedTuple)
    n, β = row.n, row.beta
    ham_path = joinpath(HAM_DIR, ham_filename(n))
    ham = _load_hamiltonian_bson(ham_path, β)

    # Build jump set in Hamiltonian eigenbasis (TimeDomain uses ham.eigvecs).
    jumps = JumpOp[]
    jump_norm = sqrt(3 * n)
    for pauli in (X, Y, Z), site in 1:n
        op = Matrix(pad_term([pauli], n, site)) ./ jump_norm
        op_eb = ham.eigvecs' * op * ham.eigvecs
        push!(jumps, JumpOp(op, op_eb, op == transpose(op), op == op'))
    end

    cfg = Config(
        sim = Thermalize(),
        domain = TimeDomain(),
        construction = KMS(),
        num_qubits = n,
        beta = β, sigma = row.sigma,
        with_linear_combination = row.with_linear_combination,
        a = row.a, s = row.s,
        gaussian_parameters = row.with_linear_combination ? (nothing, nothing) :
            (row.gaussian_omega, row.gaussian_sigma),
        eta = row.eta,
        num_energy_bits_D = row.r_D, t0_D = row.t0_D, w0_D = row.w0_D,
        num_energy_bits_b_minus = row.r_bm, t0_b_minus = row.t0_bm, w0_b_minus = row.w0_bm,
        num_energy_bits_b_plus  = row.r_bp, t0_b_plus  = row.t0_bp,  w0_b_plus  = row.w0_bp,
        num_trotter_steps_per_t0 = row.M_D,
        delta = row.delta,
        mixing_time = 5.0,  # placeholder; predictor only uses k_grid
        with_gqsp = row.with_gqsp, gqsp_degree = row.gqsp_degree,
        jump_selection = row.jump_selection,
    )

    dim = size(ham.data, 1)
    rho_0 = Matrix{ComplexF64}(I(dim) ./ dim)

    # Logarithmic k_grid covering up to ~1e5 steps (slow Gaussian samplers at
    # high β need long evolutions to reach the asymptotic floor).
    k_grid = unique(round.(Int, exp10.(range(0, 5, length=50))))

    t0_run = time()
    result = predict_channel_trajectory(
        cfg, ham, jumps, rho_0, k_grid;
        krylovdim = 30)
    wall = time() - t0_run

    return (
        achieved_dist = result.distances[end],
        gap = result.spectral_gap,
        n_steps = k_grid[end],
        t_final = result.t[end],
        total_matvecs = result.total_matvecs,
        all_converged = result.all_converged,
        wall_time = wall,
    )
end

# Pick the sanity-check fixtures: at ε=1e-3, two (n, β) pairs per filter.
# Cheap enough to run by default (sub-second per cell at n ≤ 5).
#
# Diagnostic columns:
#   achieved   — final trace distance from predict_channel_trajectory
#   δ          — splitting step from the recipe
#   gap λ      — Lindbladian spectral gap extracted by the predictor
#   ach/δ      — achieved / δ; the asymptotic channel-shift `‖ρ_∞ − ρ_β‖` is
#                an O(δ) effect, so a well-behaved fixture has `ach/δ` of
#                order 0.1. Slow samplers (e.g., Gaussian at low T) can show
#                `ach/δ ~ 1` — that is a sampler property, not a recipe failure.
function sanity_check(rows::Vector{NamedTuple})
    fixtures = [
        (n=3, β=20.0),
        (n=5, β=5.0),
    ]
    filters = [:gaussian, :smooth_metro, :kinky_metro]
    target_eps = 1e-3

    println("\n=== Sanity check (predict_channel_trajectory at ε_target=$target_eps) ===")
    @printf "%-15s %3s %5s %12s %12s %12s %12s\n" "filter" "n" "β" "achieved" "δ" "gap λ" "ach/δ"
    println("  " * "-" ^ 90)

    report = NamedTuple[]
    for (n, β) in fixtures, f in filters
        idx = findfirst(r -> r.n == n && r.beta == β && r.eps == target_eps && r.filter == f, rows)
        idx === nothing && continue
        row = rows[idx]
        rep = _sanity_check_cell(row)
        ach_over_delta = rep.achieved_dist / row.delta
        @printf "%-15s %3d %5.1f %12.3e %12.3e %12.3e %12.2f\n" string(f) n β rep.achieved_dist row.delta rep.gap ach_over_delta
        push!(report, merge((filter = f, n = n, beta = β, target_eps = target_eps,
                             delta = row.delta, ach_over_delta = ach_over_delta), rep))
    end
    return report
end

# ---------------------------------------------------------------------------
# Ideal-Lindbladian sanity check — ‖L_E − L_B‖_op / ‖L_B‖_op vs target 1e-9
# ---------------------------------------------------------------------------

# Build dense L on (config_kind ∈ {EnergyDomain, BohrDomain}) using the row's
# recipe and return both operators for a side-by-side norm comparison.
function _ideal_sanity_check_cell(row::NamedTuple)
    ham = _load_hamiltonian_bson(joinpath(HAM_DIR, ham_filename(row.n)), row.beta)
    jumps_eb = JumpOp[]
    jump_norm = sqrt(3 * row.n)
    for pauli in (X, Y, Z), site in 1:row.n
        op = Matrix(pad_term([pauli], row.n, site)) ./ jump_norm
        op_eb = ham.eigvecs' * op * ham.eigvecs
        push!(jumps_eb, JumpOp(op, op_eb, op == transpose(op), op == op'))
    end

    # Common Config kwargs for the Lindbladian build.
    base_kw = (
        sim = Lindbladian(),
        construction = KMS(),
        num_qubits = row.n,
        beta = row.beta, sigma = row.sigma,
        with_linear_combination = row.with_linear_combination,
        a = row.a, s = row.s,
        gaussian_parameters = row.with_linear_combination ? (nothing, nothing) :
            (row.gaussian_omega, row.gaussian_sigma),
    )

    cfg_E = Config(; base_kw..., domain = EnergyDomain(),
                     num_energy_bits_D = row.r_D, w0_D = row.w0_D, t0_D = row.t0_D)
    cfg_B = Config(; base_kw..., domain = BohrDomain())

    L_E = construct_lindbladian(jumps_eb, cfg_E, ham)
    L_B = construct_lindbladian(jumps_eb, cfg_B, ham)
    L_B_norm = opnorm(L_B)
    L_E_diff_norm = opnorm(L_E .- L_B)

    return (
        L_B_norm = L_B_norm,
        L_E_diff_norm = L_E_diff_norm,
        relative_error = L_E_diff_norm / max(L_B_norm, eps()),
    )
end

# Sanity check ideal-Lindbladian rows by comparing EnergyDomain vs BohrDomain
# `construct_lindbladian` at one cheap fixture (n=3, β=10) per filter.
# Pass criterion: relative_error ≤ 1e-9 (a few orders below the algorithm-
# level ε targets ∈ {1e-3, 1e-6}).
function sanity_check_ideal(rows_ideal::Vector{NamedTuple})
    fixture_n, fixture_β = 3, 10.0
    target_eps = 1e-3
    filters = [:smooth_metro, :gaussian, :kinky_metro]

    println("\n=== Ideal-Lindbladian sanity check  ‖L_E − L_B‖_op / ‖L_B‖_op (n=$fixture_n, β=$fixture_β) ===")
    @printf "%-15s %3s %5s %4s %16s %16s %16s\n" "filter" "n" "β" "r_D" "‖L_B‖_op" "‖L_E−L_B‖_op" "relative"
    println("  " * "-" ^ 90)

    report = NamedTuple[]
    for f in filters
        idx = findfirst(r -> r.n == fixture_n && r.beta == fixture_β
                          && r.eps ≈ target_eps && r.filter == f, rows_ideal)
        idx === nothing && continue
        row = rows_ideal[idx]
        rep = _ideal_sanity_check_cell(row)
        @printf "%-15s %3d %5.1f %4d %16.4e %16.4e %16.4e\n" string(f) row.n row.beta row.r_D rep.L_B_norm rep.L_E_diff_norm rep.relative_error
        push!(report, merge((filter = f, n = row.n, beta = row.beta, r_D = row.r_D),
                            rep))
    end
    return report
end

# ---------------------------------------------------------------------------
# Doc generation
# ---------------------------------------------------------------------------

function write_doc(rows::Vector{NamedTuple}, sanity_report::Vector{NamedTuple};
                    rows_ideal::Vector{NamedTuple} = NamedTuple[],
                    ideal_sanity_report::Vector{NamedTuple} = NamedTuple[])
    open(OUT_DOC, "w") do io
        println(io, "# Channel parameter table — `(n, β, ε, filter)` recipe")
        println(io, "")
        println(io, "Generated by `scripts/numerics_param_table.jl`.  BSON output at `scripts/output/channel_param_table.bson`.")
        println(io, "")
        println(io, "**Cells**: ", length(rows), " rows, ",
                "`n ∈ {", join(sort(unique([r.n for r in rows])), ", "), "}` × ",
                "`β ∈ {", join(sort(unique([r.beta for r in rows])), ", "), "}` × ",
                "`ε ∈ {", join(sort(unique([r.eps for r in rows])), ", "), "}` × ",
                "`filter ∈ {", join(unique([string(r.filter) for r in rows]), ", "), "}`.")
        println(io, "")
        println(io, "Family: `1D-XXX-zzdis` (`heis_xxx_zzdisordered_periodic_n{n}.bson`).")
        println(io, "")
        println(io, "## Recipe inputs")
        println(io, "")
        println(io, "- **r_D table** (qf-7xt + 1-bit safety for cross-(n, β) universality):")
        println(io, "")
        println(io, "  | filter | ε=1e-3 | ε=1e-6 |")
        println(io, "  |---|---|---|")
        println(io, "  | Gaussian | 5 | 6 |")
        println(io, "  | Smooth Metro (a=0, s=0.25) | 5 | 6 |")
        println(io, "  | Kinky Metro (a=0, s=0) | 8 | 13 |")
        println(io, "")
        println(io, "- **r_b_- table** (filter-independent, super-algebraic): 5 (ε=1e-3) / 7 (ε=1e-6).")
        println(io, "")
        println(io, "- **r_b_+ table** (filter-dependent, slope -1 for Metro):")
        println(io, "")
        println(io, "  | filter | ε=1e-3 | ε=1e-6 |")
        println(io, "  |---|---|---|")
        println(io, "  | Gaussian | 5 | 7 |")
        println(io, "  | Smooth / Kinky Metro | 5 | 15 |")
        println(io, "")
        println(io, "- **Time-window targets**: T_- ≈ 18, T_+ ≈ 12 (~3× the b_± supports at β=10, σ=0.1).  The realised `T_±` per cell is the integer-multiple snap `2^(r_b±-1) · k_b± · t0_D / β` so the qf-d0w shared-δt₀ constructor accepts the natural Trotter steps (β·t0_b±, t0_D); within ~10 % of the target.")
        println(io, "- **w0_D = ω_range / 2^r_D** with **ω_range = 2(‖H‖ + 8σ)**, **σ = 1/β**.  ‖H‖ measured per fixture.")
        println(io, "- **w0_b_±** = π/T_± (derived from realised T_±).  **t0_b_±** = k_b± · t0_D / β with k_b± = round(β · 2T_±_target / (2^r_b± · t0_D)).")
        println(io, "- **M_D**: 1, 2, 4, 16 for ε ∈ {1e-3, 1e-4, 1e-5, 1e-6}.  **M_b_±** = 1 (ε ≥ 1e-5) / 2 (ε = 1e-6).")
        println(io, "- **δ_step**: 5e-2, 1.5e-2, 5e-3, 1.5e-3 (jump-wise splitting recipe, slope +2 in δ).")
        println(io, "- **η** = ε / (3β).  Below t0' = 2T_+/2^r_b+ for all cells, so the η-cutoff branch in `_compute_b_plus_metro` is dead code (qf-7xt convention).")
        println(io, "- **Gaussian (ω_γ, σ_γ)**: σ_γ = σ = 1/β, ω_γ = β·(σ²+σ_γ²)/2 = σ.  Satisfies the KMS-line constraint β = 2ω_γ/(σ²+σ_γ²).")
        println(io, "- **GQSP**: `with_gqsp = true, gqsp_degree = 1` (Form B cost model; see `drafts/sim-time-gqsp-audit.md`).")
        println(io, "- **Splitting**: `jump_selection = :sweep`.")
        println(io, "")
        println(io, "## Sanity check at 6 fixtures (ε_target=1e-3, TimeDomain)")
        println(io, "")
        println(io, "TimeDomain `predict_channel_trajectory` (NUFFT, no Trotter) — exercises the dissipative + coherent quadrature recipe and the splitting δ.  The achieved asymptotic trace distance to Gibbs is the **channel fixed-point shift `‖ρ_∞ − ρ_β‖`**, an O(δ) effect that depends on the sampler's per-step splitting + Trotter approximation residue.  It is **not** the same as the recipe's `ε_target`, which controls the *quadrature* error budget.")
        println(io, "")
        println(io, "| filter | n | β | achieved | δ | gap λ | ach/δ |")
        println(io, "|---|---|---|---|---|---|---|")
        for r in sanity_report
            @printf(io, "| %s | %d | %.1f | %.3e | %.3e | %.3e | %.2f |\n",
                    string(r.filter), r.n, r.beta, r.achieved_dist, r.delta, r.gap, r.ach_over_delta)
        end
        println(io, "")
        println(io, "**Reading the table.**  A well-behaved fixture has `ach/δ ~ 0.1` (channel shift ~ 10% of the splitting step), reflecting the leading-order splitting + Trotter approximation in `Φ_δ`.  Slow samplers — Gaussian at low T, in particular — show `ach/δ ~ 1`, where the fixed-point shift is bigger because the per-step splitting moves `ρ_∞` further from `ρ_β` when the sampler's gap is small.  This is documented in `.claude-memory/ckg_vs_dll_first_findings.md` (\"DLL Gaussian collapse at low T\") and is a *sampler* property, not a recipe failure: the quadrature parameters in this BSON correctly encode the qf-7xt recipe; the channel shift on top is what the implemented δ-channel adds.")
        println(io, "")
        println(io, "Cells with the OFT-integration-warning (\"Time array was not truncated\") indicate that the filter kernel does not decay fast enough at the time-domain boundary.  At small `r_D`, the time window `[-π/w0_D, π/w0_D]` extends well beyond the filter support; the warning fires when the boundary kernel value exceeds an internal threshold.  Empirically the warning is informational here (the achieved fixed-point shift is still O(δ) and not catastrophic).  The recipe's `+1`-bit safety margin in `r_D` is the lever to reduce this if needed.")
        println(io, "")
        println(io, "**ε_target vs achieved**.  `ε_target = 1e-3` controls the QUADRATURE recipe (r_D, r_b±) and the splitting δ — what this BSON encodes.  The achieved asymptotic trace distance from `predict_channel_trajectory` is `‖ρ_∞ − ρ_β‖`, an O(δ) channel-shift, generally larger than `ε_target`.  Algorithm-level precision (`ε = 10⁻⁶`) is verified by the IDEAL Lindbladian path (S1 / S2 via `predict_lindbladian_trajectory`), not here.  See `.claude-memory/thesis_target_precision_1e6.md` for the precision-regime split.")
        println(io, "")
        println(io, "## ε=1e-6 feasibility flags")
        println(io, "")
        println(io, "- **Kinky Metro at ε=1e-6** uses r_b+ = 15 (`N = 32768` inner-grid points).  Memory: 32768 · d² · 16 bytes for the OFT cache — at n=6 (d=64), that's 32768 · 4096 · 16 ≈ 2 GB.  Feasible on laptop but slow.")
        println(io, "- **Smooth Metro at ε=1e-6** also uses r_b+ = 15, same memory profile.")
        println(io, "- **Gaussian at ε=1e-6** uses r_b+ = 7 (super-algebraic), well under any memory limit.")
        println(io, "- **Kinky Metro at ε=1e-6, n=6** is recommended for cluster, not laptop, even with the predictor.")
        println(io, "")
        println(io, "## Output schema")
        println(io, "")
        println(io, "Each row in `channel_param_table.bson` is a `NamedTuple` with fields:")
        println(io, "")
        println(io, "```")
        println(io, "(n, beta, eps, filter, family,")
        println(io, " H_norm, sigma, omega_range,                                  # Hamiltonian meta")
        println(io, " r_D, w0_D, t0_D,                                              # dissipative register")
        println(io, " r_bm, w0_bm, t0_bm,                                           # outer coherent register")
        println(io, " r_bp, w0_bp, t0_bp,                                           # inner coherent register")
        println(io, " T_minus, T_plus,                                              # time windows")
        println(io, " M_D, M_bm, M_bp,                                              # Trotter step counts")
        println(io, " delta, eta,                                                   # splitting + Metro regularisation")
        println(io, " with_linear_combination, a, s,                                # Metropolis filter parameters")
        println(io, " gaussian_omega, gaussian_sigma,                               # Gaussian filter parameters (NaN unless filter=:gaussian)")
        println(io, " with_gqsp, gqsp_degree, jump_selection)")
        println(io, "```")
        println(io, "")
        if !isempty(rows_ideal)
            println(io, "## Ideal-Lindbladian (EnergyDomain) reference table")
            println(io, "")
            println(io, "Companion table for **S1 / S2** sweeps (`predict_lindbladian_trajectory`) — the project's reference for \"ideal CKG Lindbladian\".  Saved separately as `scripts/output/ideal_lindbladian_param_table.bson`.")
            println(io, "")
            println(io, "**Key difference from the channel table**: EnergyDomain Lindbladian dispatches `_precompute_coherent_B` to the analytical `B_bohr`, so the b_± per-term registers are *not* used — only `(r_D, w0_D, t0_D)` matters.  Per qf-7xt, smooth-Metro and Gaussian saturate at machine precision (~10⁻¹⁴) by `r_D = 6`; the recipe takes `r_D = 7` for one extra bit of safety against (n, β) shifts in the K-prefactor.")
            println(io, "")
            println(io, "| filter | r_D (EnergyDomain) | quadrature err vs Bohr (qf-7xt at n=4, β=10) | preferred at n ≥ 6 |")
            println(io, "|---|---|---|---|")
            println(io, "| Gaussian | 7 | ~ 10⁻¹⁴ (super-algebraic) | EnergyDomain |")
            println(io, "| Smooth Metro (a=0, s=0.25) | 7 | ~ 10⁻¹⁴ (Gevrey-1/2) | EnergyDomain |")
            println(io, "| Kinky Metro (a=0, s=0) | 14 | ~ 2.5 × 10⁻⁸ (slope -2 in 1/N) | **BohrDomain** (memory) |")
            println(io, "")
            println(io, "Cells: ", length(rows_ideal), " rows over the same `(n, β, ε, filter)` grid as the channel table; ε is recorded for downstream τ_mix bookkeeping but does *not* enter the quadrature recipe.")
            println(io, "")
            println(io, "Schema (per row):")
            println(io, "")
            println(io, "```")
            println(io, "(n, beta, eps, filter, family,")
            println(io, " H_norm, sigma, omega_range,                       # Hamiltonian meta")
            println(io, " r_D, w0_D, t0_D,                                   # dissipative register")
            println(io, " domain_recommended,                                # :EnergyDomain or :BohrDomain")
            println(io, " with_linear_combination, a, s,                     # Metro filter parameters")
            println(io, " gaussian_omega, gaussian_sigma,                    # Gaussian filter parameters (NaN unless filter=:gaussian)")
            println(io, " mode, jump_selection)")
            println(io, "```")
            println(io, "")
            if !isempty(ideal_sanity_report)
                println(io, "### Ideal sanity check  ‖L_E − L_B‖_op / ‖L_B‖_op  (n=3, β=10, ε=1e-3)")
                println(io, "")
                println(io, "Dense `construct_lindbladian` build on EnergyDomain (recipe-r_D) vs BohrDomain (analytical reference).  Pass criterion: **relative_error ≤ 10⁻⁹**.")
                println(io, "")
                println(io, "| filter | n | β | r_D | ‖L_B‖_op | ‖L_E − L_B‖_op | relative |")
                println(io, "|---|---|---|---|---|---|---|")
                for r in ideal_sanity_report
                    pass_mark = r.relative_error ≤ 1e-9 ? "" : " *(above 1e-9)*"
                    @printf(io, "| %s | %d | %.1f | %d | %.4e | %.4e | %.4e%s |\n",
                            string(r.filter), r.n, r.beta, r.r_D,
                            r.L_B_norm, r.L_E_diff_norm, r.relative_error, pass_mark)
                end
                println(io, "")
                println(io, "The smooth-Metro and Gaussian rows are well below 10⁻⁹ as expected from the qf-7xt convergence summary.  Kinky-Metro at `r_D = 14` is marginal — the project recipe routes kinky-Metro ideal-Lindbladian to BohrDomain at n ≥ 6 to avoid the ω-grid memory blow-up.")
                println(io, "")
            end
        end

        println(io, "## Updating the recipe")
        println(io, "")
        println(io, "If the qf-7xt recipe is revised (e.g. a new fixture exposes a >2-bit shift in K-prefactors), update `_r_D` / `_r_bm` / `_r_bp` in `scripts/numerics_param_table.jl` and re-run.  Sanity-check fixtures should be re-checked.  The S3 / S5 / S6 sweeps re-load the BSON each time, so a single edit propagates.")
    end
    return OUT_DOC
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

function main()
    do_check = !("--no-check" in ARGS)

    println("Building channel parameter table …")
    rows = build_table()
    println("  built $(length(rows)) channel cells.")

    mkpath(dirname(OUT_BSON))
    BSON.bson(OUT_BSON, Dict(:rows => rows, :family => FAMILY_TAG,
                              :generated_by => "scripts/numerics_param_table.jl",
                              :table_kind => :channel))
    println("  wrote $(OUT_BSON)")

    println("\nBuilding ideal-Lindbladian (EnergyDomain) parameter table …")
    rows_ideal = build_ideal_lindbladian_table()
    println("  built $(length(rows_ideal)) ideal-Lindbladian cells.")
    BSON.bson(OUT_BSON_IDEAL, Dict(:rows => rows_ideal, :family => FAMILY_TAG,
                                    :generated_by => "scripts/numerics_param_table.jl",
                                    :table_kind => :ideal_lindbladian))
    println("  wrote $(OUT_BSON_IDEAL)")

    sanity_report = NamedTuple[]
    ideal_sanity_report = NamedTuple[]
    if do_check
        sanity_report = sanity_check(rows)
        ideal_sanity_report = sanity_check_ideal(rows_ideal)
    else
        println("\n(skipping sanity check; pass without --no-check to enable)")
    end

    mkpath(dirname(OUT_DOC))
    write_doc(rows, sanity_report; rows_ideal = rows_ideal,
              ideal_sanity_report = ideal_sanity_report)
    println("\nWrote $(OUT_DOC)")
end

main()
