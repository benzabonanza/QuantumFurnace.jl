# OFT filter abstraction.
#
# A filter is the function `f(t) ↔ f̂(ν)` plugged into the Operator Fourier
# Transform (OFT). The CKG construction (Chen–Kastoryano–Gilyén) uses a real
# Gaussian filter `f(t) = exp(-σ² t²)`. The DLL construction (Ding–Li–Lin 2024,
# "Efficient quantum Gibbs samplers with Kubo–Martin–Schwinger detailed balance")
# uses a Gaussian-type Gevrey filter that already absorbs the KMS factor
# `e^{-βν/4}` into `f̂(ν)`; its time-domain kernel is complex-valued.
#
# Phase 50 / DLL-1 integration scope: Time/Trotter NUFFT path only. The
# EnergyDomain `oft!` path stays Gaussian (DLL EnergyDomain integration is
# DLL-2/3). See `.planning/phases/50-dll-1-filter-integration/PLAN.md` and the
# verified prototype at `scripts/scratch_dll_filter.jl`.
#
# Fourier convention (matches DLL paper Eq. 3.3 and the NUFFT prefactor pipeline
# in `_prepare_oft_nufft_prefactors`):
#   inverse:  f(t)  = (1/2π) ∫ f̂(ν) e^{-i t ν} dν     (DLL Eq. 3.3)
#   forward:  f̂(ν) = ∫ f(t) e^{+i ν t} dt

"""
    AbstractFilter

Supertype for OFT filtering functions `f(t) ↔ f̂(ν)`. Concrete subtypes
implement `time_kernel(filter, t)` and `freq_kernel(filter, ν)`, plus a
`filter_time_cutoff(filter, tol)` helper used by the OFT time-label
truncation in `_truncate_time_labels_for_oft`.
"""
abstract type AbstractFilter end

"""
    GaussianFilter{T<:AbstractFloat}(sigma::T)

CKG (Chen–Kastoryano–Gilyén) Gaussian filter: `f(t) = exp(-σ² t²)`,
`f̂(ν) ∝ exp(-ν²/(4σ²))`. Mirrors the hard-coded form in `src/ofts.jl`
and `src/nufft.jl`. `time_kernel` returns a real `T`-typed value.
"""
struct GaussianFilter{T<:AbstractFloat} <: AbstractFilter
    sigma::T
end

"""
    DLLGaussianFilter{T<:AbstractFloat}(beta::T)

DLL Gaussian-type Gevrey filter (Ding–Li–Lin 2024, Eq. 3.21–3.22):

    q(ν)   = exp(-(βν)²/8)                              (Eq. 3.21, w ≡ 1)
    f̂(ν)  = q(ν) e^{-βν/4} = e^{1/8} exp(-(βν+1)²/8)   (Eq. 3.22)
    f(t)   = (e^{1/8} √(2/π)/β) exp(-2t²/β² + i t/β)    (Eq. 3.3 inverse FT)

`f̂(ν)` already includes the KMS factor `e^{-βν/4}`, so `freq_kernel` returns
the full DLL filter in frequency space (no extra normalisation). The time
kernel is complex-valued (centred frequency `ν_* = -1/β`, width `2/β`).

# PHYSICS CHECK: `w ≡ 1` (no compact-support bump) violates Assumption 15
of the paper but matches Fig. 1 and Eq. 3.22 verbatim. The compact-support
bump is needed only for the rigorous Paley–Wiener decay argument
(Lemma 30); we omit it here as in the paper's numerical illustration.
"""
struct DLLGaussianFilter{T<:AbstractFloat} <: AbstractFilter
    beta::T
end

# Filter element type: filters always produce complex-valued time kernels at
# the OFT call site (the Gaussian is real but stored as Complex{T} in the
# NUFFT prefactor stack; the DLL kernel is genuinely complex).
Base.eltype(::GaussianFilter{T}) where {T} = Complex{T}
Base.eltype(::DLLGaussianFilter{T}) where {T} = Complex{T}

# ---------------------------------------------------------------------------
# Generic stubs (concrete methods below)
# ---------------------------------------------------------------------------

function time_kernel end
function freq_kernel end
function filter_time_cutoff end

# ---------------------------------------------------------------------------
# CKG Gaussian filter
# ---------------------------------------------------------------------------

"""
    time_kernel(filter::GaussianFilter, t)

Returns `exp(-σ² t²)`. The form `exp(-(σ²)·t²)` (NOT `exp(-(σ·t)²)`) is
chosen to match the hardcoded form in `src/nufft.jl::_prepare_oft_nufft_prefactors`
(`exp.(-(sigma_f64^2) .* (time_labels_f64 .^ 2))`) byte-for-byte; the two
expressions are mathematically equal but can differ at the ULP level.
"""
@inline time_kernel(f::GaussianFilter{T}, t::Real) where {T} =
    exp(-(f.sigma^2) * t^2)

"""
    freq_kernel(filter::GaussianFilter, ν)

Returns the unnormalised Gaussian `exp(-ν²/(4σ²))`. The full Fourier
transform of `time_kernel` is `√π/σ · freq_kernel`, with the `√π/σ` factor
absorbed into `oft_domain_prefactor` (see `src/furnace_utensils.jl`).
"""
@inline freq_kernel(f::GaussianFilter{T}, nu::Real) where {T} =
    exp(-nu^2 / (4 * f.sigma^2))

"""
    _time_oft_prefactor_gaussian(filter::GaussianFilter)

The Gaussian-specific scalar `√(σ √(2/π)/(2π))` used inside the OFT
time-label truncation (`_truncate_time_labels_for_oft`). Matches the
existing `time_oft_prefactor` formula in `src/time_domain.jl`.
"""
@inline _time_oft_prefactor_gaussian(f::GaussianFilter{T}) where {T} =
    sqrt(f.sigma * sqrt(T(2) / T(pi)) / (2 * T(pi)))

"""
    filter_time_cutoff(filter::GaussianFilter, tol)

Solves `exp(-σ² t²) ≤ tol / prefactor` for `|t|`, giving
`cutoff = √(log(prefactor / tol)) / σ`. Matches the existing CKG cutoff
formula in `src/time_domain.jl`.
"""
@inline filter_time_cutoff(f::GaussianFilter{T}, tol::Real) where {T} =
    sqrt(log(_time_oft_prefactor_gaussian(f) / tol)) / f.sigma

# ---------------------------------------------------------------------------
# DLL Gaussian-type Gevrey filter (Ding–Li–Lin 2024, Sec. 3.2)
# ---------------------------------------------------------------------------

"""
    q_weight(filter::DLLGaussianFilter, ν)

The DLL paper's `q(ν) = exp(-(βν)²/8)` (Eq. 3.21, with the compact-support
bump `w ≡ 1`). This is the DB weight before the KMS factor `e^{-βν/4}` is
applied; `freq_kernel(filter, ν) == q_weight(filter, ν) * exp(-β ν / 4)`.
"""
@inline q_weight(f::DLLGaussianFilter{T}, nu::Real) where {T} =
    exp(-(f.beta * nu)^2 / 8)

"""
    freq_kernel(filter::DLLGaussianFilter, ν)

Returns `f̂(ν) = q(ν) e^{-βν/4} = e^{1/8} exp(-(βν+1)²/8)` from
Eq. 3.22 of Ding–Li–Lin 2024 (rearranged by completing the square). This
is the FULL DLL filter in frequency space — the KMS factor is already
included and no further normalisation is needed at the OFT stage.
"""
@inline freq_kernel(f::DLLGaussianFilter{T}, nu::Real) where {T} =
    exp(T(1) / 8) * exp(-(f.beta * nu + 1)^2 / 8)

"""
    time_kernel(filter::DLLGaussianFilter, t)

Closed-form inverse FT of `freq_kernel` (Eq. 3.3):

    f(t) = (e^{1/8} √(2/π) / β) · exp(-2 t² / β² + i t / β)

Returns `Complex{T}`. Frequency-domain centre is `ν_* = -1/β` and width
`σ_ν = 2/β`; in the time domain the modulus decays as `exp(-2t²/β²)`,
i.e. effective width `β/2`.
"""
@inline function time_kernel(f::DLLGaussianFilter{T}, t::Real) where {T}
    pref = exp(T(1) / 8) * sqrt(T(2) / T(pi)) / f.beta
    decay = exp(-2 * t^2 / f.beta^2)
    phase = cis(t / f.beta)  # cis(x) = exp(ix)
    return pref * decay * phase
end

"""
    _time_oft_prefactor_dll(filter::DLLGaussianFilter)

The DLL time-domain prefactor `e^{1/8} √(2/π)/β`, equal to `|f(0)|` from
the closed form above. Used by `filter_time_cutoff` to bound the residual
in `_truncate_time_labels_for_oft`.
"""
@inline _time_oft_prefactor_dll(f::DLLGaussianFilter{T}) where {T} =
    exp(T(1) / 8) * sqrt(T(2) / T(pi)) / f.beta

"""
    filter_time_cutoff(filter::DLLGaussianFilter, tol)

Solves `|f(t)| = prefactor · exp(-2 t² / β²) ≤ tol` for `|t|`:
`cutoff = (β/√2) · √(log(prefactor / tol))`. Strictly larger than the
Gaussian cutoff at the same `tol` because the DLL kernel is wider in time.
"""
@inline filter_time_cutoff(f::DLLGaussianFilter{T}, tol::Real) where {T} =
    (f.beta / sqrt(T(2))) * sqrt(log(_time_oft_prefactor_dll(f) / tol))

# ---------------------------------------------------------------------------
# Hörmander mollifier — Gevrey bump for the DLL Metropolis-type filter
# ---------------------------------------------------------------------------
#
# Ding–Li–Lin 2024 Assumption 15 (p. 14) and Eq. 3.19 require a compactly
# supported Gevrey bump `w` with `supp(w) ⊂ [-1, 1]` and `w ≡ 1` on
# `[-1/2, 1/2]`. The canonical construction (Hörmander, also [AHR17, Cor 2.8]):
#
#     η(t) = exp(-1/t)                    for t > 0,  0 otherwise
#     φ(t) = η(t) / (η(t) + η(1-t))       smooth from 0 (at t≤0) to 1 (at t≥1)
#     w(x) = φ(2(1-|x|))                  even, =1 on |x|≤1/2, =0 on |x|≥1
#
# The result is C^∞, even, satisfies `w ≡ 1` on `[-1/2, 1/2]` (so the bump is
# *invisible* on the flat-top region), and decays to 0 at the boundary
# `|x| = 1`. This is the Gevrey witness referenced in Assumption 15 — the
# bump's order `s_w = 2` (because `η(t) = exp(-1/t)` has derivative growth
# `(α!)^2`) gives the Paley-Wiener decay used in Lemma 30.
#
# Floating-point details: for `t ≲ 0`, `exp(-1/t)` underflows to 0; the
# branches below short-circuit those cases to avoid spurious NaN. For
# `t ≳ 1`, `η(1-t)` underflows similarly and `φ` saturates at 1.

"""
    _hormander_eta(t)

Helper bump factor `η(t) = exp(-1/t)` for `t > 0`, `0` for `t ≤ 0`. Pure
floating-point — no special-function calls.
"""
@inline function _hormander_eta(t::T) where {T<:AbstractFloat}
    return t > zero(T) ? exp(-one(T) / t) : zero(T)
end

"""
    _hormander_phi(t)

Smooth step `φ(t) = η(t) / (η(t) + η(1-t))`, `φ(t) = 0` for `t ≤ 0`,
`φ(t) = 1` for `t ≥ 1`, smooth on `[0, 1]`. Used to build the
unit-interval Hörmander bump.
"""
@inline function _hormander_phi(t::T) where {T<:AbstractFloat}
    if t <= zero(T)
        return zero(T)
    elseif t >= one(T)
        return one(T)
    end
    a = _hormander_eta(t)
    b = _hormander_eta(one(T) - t)
    return a / (a + b)
end

"""
    _hormander_bump(x)

Hörmander mollifier on `[-1, 1]`: even, `w(x) = 1` on `|x| ≤ 1/2`, smooth
decay to 0 at `|x| = 1`, `w(x) = 0` for `|x| ≥ 1`. Matches Ding–Li–Lin
Assumption 15 exactly. Type-preserving in `T <: AbstractFloat`.
"""
@inline function _hormander_bump(x::T) where {T<:AbstractFloat}
    ax = abs(x)
    if ax >= one(T)
        return zero(T)
    elseif ax <= one(T) / 2
        return one(T)
    end
    return _hormander_phi(T(2) * (one(T) - ax))
end

# Promote-to-Float64 fallback for non-AbstractFloat inputs (e.g. `Int`).
@inline _hormander_bump(x::Real) = _hormander_bump(float(x))

# ---------------------------------------------------------------------------
# DLL Metropolis-type Gevrey filter (Ding–Li–Lin 2024, Eq. 3.19–3.20)
# ---------------------------------------------------------------------------
#
# The Metropolis-type weighting (Eq. 3.19) realises a smoothed Metropolis
# acceptance: on the flat-top region |ν| ≤ S/2 of the bump, Eq. 3.20 gives
#
#   f̂(ν) = q(ν) · exp(-βν/4) ≈ min{1, exp(-βν/2)},
#
# matching the classical Metropolis weight. Compared to the Gaussian-type
# filter (DLLGaussianFilter, Eq. 3.21), the Metropolis support stays O(1) in
# β; in particular, |α_{ν,ν'}| = |f̂(ν) f̂(ν')*| stays O(1) for ν, ν' ≤ 0
# even at large β, whereas the Gaussian shrinks as O(β^{-1}).
#
# Fields:
#   beta — inverse temperature β > 0
#   S    — bump support radius. The bump w(ν/S) has flat top on |ν| ≤ S/2
#          and zero outside |ν| ≥ S, so f̂ is genuinely compactly supported.
#          Default S=2 matches the rescaled spectra used in the test fixtures
#          (eigvals ∈ [-0.45, 0.45], so Bohr freqs |ν| ≤ 0.9 < S/2 = 1, and
#          the bump is *invisible* on the relevant grid).

"""
    DLLMetropolisFilter{T<:AbstractFloat}(beta::T; S::Real = T(2))

DLL Metropolis-type Gevrey filter (Ding–Li–Lin 2024, Eq. 3.19–3.20):

    u(x)   = exp(-√(1+x²)/4)                                 (Eq. 3.19, smooth |x|/4)
    q(ν)   = u(βν) · w(ν/S) = exp(-√(1+(βν)²)/4) · w(ν/S)    (Eq. 3.19)
    f̂(ν)  = q(ν) e^{-βν/4} ≈ min{1, e^{-βν/2}}              (Eq. 3.20, on |ν| ≤ S/2)

`f̂(ν)` already includes the KMS factor `e^{-βν/4}`, so `freq_kernel` returns
the full DLL filter in frequency space. The bump `w` is the Hörmander
mollifier (`_hormander_bump`), which is `≡ 1` on `|ν/S| ≤ 1/2` (i.e. on
`|ν| ≤ S/2`) and vanishes for `|ν| ≥ S`. Default `S = 2` puts the flat top
at `[-1, 1]`, which contains the Bohr frequencies of the test fixtures
(rescaled `‖H‖ ≤ 0.45`); the bump is then mathematically present but
numerically invisible.

# PHYSICS CHECK: `S = 2` keyword default. The caller must ensure
`S/2 ≥ max|ν_BH|` for the chosen Hamiltonian — otherwise the bump bites
the Lindbladian (only the central region is Metropolis; outside `[-S/2, S/2]`
the asymptote `min{1, exp(-βν/2)}` does not hold). `validate_config!` only
checks `S > 0` (Hamiltonian-aware bound check is deferred to the caller).

# Examples
```julia
filter = DLLMetropolisFilter(5.0)            # β=5, S=2
filter = DLLMetropolisFilter(10.0; S = 5.0)  # β=10, S=5 (wider support)
```
"""
struct DLLMetropolisFilter{T<:AbstractFloat} <: AbstractFilter
    beta::T
    S::T
end

DLLMetropolisFilter(beta::T; S::Real = T(2)) where {T<:AbstractFloat} =
    DLLMetropolisFilter{T}(beta, T(S))

# Like DLLGaussianFilter, the time kernel is complex-valued (asymmetric
# `f̂(ν)` ⇒ complex `f(t)` via inverse Fourier transform).
Base.eltype(::DLLMetropolisFilter{T}) where {T} = Complex{T}

"""
    q_weight(filter::DLLMetropolisFilter, ν)

The DLL Metropolis weight `q(ν) = exp(-√(1+(βν)²)/4) · w(ν/S)` from
Eq. 3.19 of Ding–Li–Lin 2024. This is the DB weight before the KMS factor
`e^{-βν/4}` is applied; `freq_kernel(filter, ν) == q_weight(filter, ν) * exp(-β ν / 4)`.

The bump `w(ν/S)` strictly enforces `supp(q) ⊂ [-S, S]` (Assumption 15);
on the flat-top region `|ν| ≤ S/2` it equals 1 and the smooth Metropolis
shape `exp(-√(1+(βν)²)/4)` is uncovered.
"""
@inline q_weight(f::DLLMetropolisFilter{T}, nu::Real) where {T} =
    exp(-sqrt(one(T) + (f.beta * nu)^2) / 4) * _hormander_bump(nu / f.S)

"""
    freq_kernel(filter::DLLMetropolisFilter, ν)

Returns the full Metropolis filter in frequency space:

    f̂(ν) = q(ν) · e^{-βν/4} = exp(-(√(1+(βν)²) + βν)/4) · w(ν/S)        (Eq. 3.20)

On the flat-top region `|ν| ≤ S/2` (where `w = 1`):
- ν ≪ 0: `√(1+(βν)²) + βν → -βν + βν = 0`, so `f̂(ν) → 1`     (accept down-jump)
- ν ≫ 0: `√(1+(βν)²) + βν → 2βν`, so `f̂(ν) → e^{-βν/2}`      (Metropolis up-jump)
- ν = 0: `f̂(0) = e^{-1/4} ≈ 0.778`                          (smoothed corner)

Outside `|ν| ≥ S`, `f̂(ν) = 0` exactly (compact support).
"""
@inline freq_kernel(f::DLLMetropolisFilter{T}, nu::Real) where {T} =
    q_weight(f, nu) * exp(-f.beta * nu / 4)

"""
    time_kernel(filter::DLLMetropolisFilter, t)

Numerical inverse Fourier transform of `freq_kernel` (Eq. 3.3 of DLL):

    f(t) = (1/(2π)) ∫_{-S}^{S} f̂(ν) e^{-i t ν} dν

Computed via `QuadGK.quadgk` over the compact support `[-S, S]` (since
`freq_kernel ≡ 0` outside). Returns `Complex{T}`.

The integrand `f̂(ν) e^{-i t ν}` is `C^∞` thanks to the Hörmander bump
(all derivatives of `f̂` vanish at `±S`), so QuadGK converges fast — at
high `|t|`, oscillation slows convergence but the absolute tolerance
`atol = eps(T) · 64` is reached in tens of nodes per panel. Allocation
cost ~10 μs at typical settings; the OFT call site (`_prepare_oft_nufft_prefactors`)
calls this once per time-grid point in a single broadcast, so the per-build
cost scales like `O(N_t · QuadGK_cost)` ≈ tens of milliseconds for
`N_t ~ 4000`.

# Implementation note
We use a relative tolerance of `1e-12` plus an absolute floor of `64 ε`
to guard against catastrophic cancellation at large `|t|` (where the
integral is small). The Float32 path uses `Float64` arithmetic internally
and converts back at the end (QuadGK has no Float32 specialisation).
"""
function time_kernel(f::DLLMetropolisFilter{T}, t::Real) where {T}
    integrand = nu -> freq_kernel(f, nu) * cis(-t * nu)
    val, _ = quadgk(integrand, -f.S, f.S;
                    rtol = T(1e-12), atol = T(64) * eps(T))
    return Complex{T}(val / T(2 * π))
end

"""
    filter_time_cutoff(filter::DLLMetropolisFilter, tol)

Returns a time cutoff `t_c` such that `|f(t)| ≤ tol` for all `t ≥ t_c`.

No closed form: we use a doubling search starting from `8β`, sampling
`|time_kernel(filter, t)|` at four equispaced points within an
oscillation window `4π/S` (the natural period of `e^{-i t ν}` over the
support). The cutoff is the smallest doubled `t` for which all four
samples fall below `tol`. This handles the modulus oscillation pattern
of the IFT of an asymmetric compactly-supported `f̂`.

The search converges in `≲ 30` doublings (the IFT envelope decays
super-polynomially because `f̂ ∈ G^1_c`, with empirical rate
`≈ 0.07/√β` to `0.25/√β` per unit time). Each evaluation is one
`time_kernel` call (one `QuadGK` integration).

The returned cutoff is conservative — it is the upper bracket of the
bisection, which can be up to `2×` the true crossing point. This is
intentional: the OFT time-grid construction uses this to truncate
labels, and a moderate over-estimate is cheaper than risking a missed
boundary contribution.
"""
function filter_time_cutoff(f::DLLMetropolisFilter{T}, tol::Real) where {T}
    osc_period = T(4 * π) / f.S  # natural oscillation scale of e^{-i t ν}
    t = max(T(8) * f.beta, T(20))
    iter = 0
    while true
        # Sample the modulus at four offsets within one oscillation period.
        passed = true
        for k in 0:3
            if abs(time_kernel(f, t + k * osc_period / 4)) > tol
                passed = false
                break
            end
        end
        passed && return t
        t *= 2
        iter += 1
        if iter > 30
            error("filter_time_cutoff: doubling search exceeded 30 iterations " *
                  "(β=$(f.beta), S=$(f.S), tol=$tol).")
        end
    end
end

# ---------------------------------------------------------------------------
# Multi-channel DLL filter (qf-7go.1) — parked
# ---------------------------------------------------------------------------
#
# `DLLMultiChannelFilter`, `ShiftedSymmetricFilter`, and the
# `dll_multichannel_translates` factory live in `src/dll_multichannel.jl`
# alongside their operator overloads and simulator-coupling helpers. The
# multi-channel construction is preserved as a diagnostic that backs a few
# thesis paragraphs (qf-7go epic / `drafts/dll-multirank-taumix-findings.md`)
# but is not part of mainline simulations. New mainline DLL work should use
# the single-channel filters defined above.
