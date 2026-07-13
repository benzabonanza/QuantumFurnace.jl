# ============================================================================
# Multi-channel DLL filter (qf-7go epic) — diagnostic / parked
# ============================================================================
#
# Parametric multi-rank extension of the single-channel DLL construction:
# `q^a(ν) = Σ_ℓ q_ℓ(ν)` with one Lindblad operator per channel,
#
#     B_a = Σ_ℓ B_a^{(ℓ)} = Σ_ℓ Σ_ν f̂_ℓ(ν) A^a_ν.
#
# The corresponding Kossakowski decomposes linearly:
#
#     α^{multi}_{ν, ν'} = Σ_ℓ e^{-β(ν+ν')/4} q_ℓ(ν) conj(q_ℓ(ν'))
#                       = Σ_ℓ α^{(ℓ)}_{ν, ν'},
#
# so each per-channel α^{(ℓ)} is a rank-1 outer product and the sum has rank
# ≤ k. The dissipator and coherent-G operators sum at the operator level — the
# multi-channel `q_weight` / `freq_kernel` / `time_kernel` methods are
# diagnostic only and **never** called from the dissipator path.
#
# Status (2026-05-04): this module is no longer part of mainline simulations.
# It is preserved because the multi-rank τ_mix sweep (`drafts/dll-multirank-
# taumix-findings.md`, beads epic qf-7go) backs a few thesis paragraphs and the
# diagnostic should remain reproducible. New mainline DLL work should use the
# single-channel `DLLGaussianFilter` / `DLLMetropolisFilter` paths in
# `src/filters.jl` and `src/dll.jl` directly.
#
# All method definitions that dispatch on `DLLMultiChannelFilter` or
# `ShiftedSymmetricFilter` live in this file. Single-channel base methods
# remain in `src/filters.jl`, `src/dll.jl`, `src/jump_workers.jl`, and
# `src/furnace_utensils.jl`.
#
# Design choice: sub-filters share a single inverse-temperature `beta`. The
# construction's KMS factor `e^{-β(ν+ν')/4}` is locked across channels, so
# allowing per-channel β would break KMS-DBC at the multi-channel level. The
# constructor enforces β-consistency.

# ---------------------------------------------------------------------------
# Section 1 — DLLMultiChannelFilter struct
# ---------------------------------------------------------------------------

"""
    DLLMultiChannelFilter{T<:AbstractFloat, F<:AbstractFilter}(channels, beta)

Parametric multi-rank DLL filter (qf-7go epic). Wraps `k = length(channels)`
single-channel DLL filters `{q_ℓ}` (typically `DLLGaussianFilter`,
`DLLMetropolisFilter`, or shifted/symmetrised variants from
`dll_multichannel_translates`) into one composite filter whose multi-channel
Kossakowski is the sum of the per-channel rank-1 Kossakowskis.

# Fields
- `channels::Vector{F}`: `k ≥ 1` sub-filters; each must satisfy
  `channels[ℓ].beta ≈ beta`.
- `beta::T`: global inverse temperature (mirrored on every sub-filter).

# Construction

The dissipator (`dll_lindblad_op_*`) and coherent-G (`dll_coherent_op_*`)
operators iterate over `channels` and sum per-channel operators — the
multi-channel `q_weight`/`freq_kernel`/`time_kernel` are **diagnostic
helpers** for tests and plotting. They return the literal sum over channels
and are *not* used to assemble the dissipator (which would be wrong: the
dissipator sums per-channel α^{(ℓ)}, not the per-channel q's).

The KMS skew-symmetry α(ν, ν') = α(-ν', -ν) e^{-β(ν+ν')/2} (Eq. 4.7) holds
on the multi-channel α as long as it holds on each per-channel α^{(ℓ)} (the
sum of KMS-skew-symmetric matrices is again skew-symmetric).

# Examples

```julia
# Direct construction from a vector of channels (β must agree across channels).
ch = [DLLMetropolisFilter(5.0; S = 2.0), DLLMetropolisFilter(5.0; S = 2.0)]
multi = DLLMultiChannelFilter(ch, 5.0)

# k = 1 case: equivalent to the wrapped single channel for all per-channel
# methods (the sum has a single term).
single = DLLGaussianFilter(2.0)
DLLMultiChannelFilter([single], 2.0)
```

For the symmetrised-translates parametrisation
`q_ℓ(ν) = (q_base(ν - ν_ℓ) + q_base(ν + ν_ℓ))/√2`, see
`dll_multichannel_translates` (qf-7go.5).
"""
struct DLLMultiChannelFilter{T<:AbstractFloat, F<:AbstractFilter} <: AbstractFilter
    channels::Vector{F}
    beta::T

    function DLLMultiChannelFilter{T, F}(channels::Vector{F}, beta::T) where
            {T<:AbstractFloat, F<:AbstractFilter}
        if isempty(channels)
            throw(ArgumentError("DLLMultiChannelFilter requires at least one channel."))
        end
        beta_tol = T(10) * eps(T)
        for (ℓ, c) in enumerate(channels)
            if !hasproperty(c, :beta)
                throw(ArgumentError("DLLMultiChannelFilter channel $ℓ ($(typeof(c))) " *
                                    "lacks a `beta` field — only DLL-style filters supported."))
            end
            if !isapprox(c.beta, beta; atol = beta_tol)
                throw(ArgumentError("DLLMultiChannelFilter channel $ℓ has beta $(c.beta), " *
                                    "expected $beta (atol=$beta_tol)."))
            end
        end
        return new{T, F}(channels, beta)
    end
end

DLLMultiChannelFilter(channels::Vector{F}, beta::T) where
        {T<:AbstractFloat, F<:AbstractFilter} =
    DLLMultiChannelFilter{T, F}(channels, beta)

# Multi-channel time kernels are sums of per-channel `Complex{T}` kernels;
# `Complex{T}` is the right element type regardless of which sub-filters appear.
Base.eltype(::DLLMultiChannelFilter{T}) where {T} = Complex{T}

"""
    q_weight(filter::DLLMultiChannelFilter, ν) -> sum of per-channel q_weight

Diagnostic-only multi-channel q. Equals `Σ_ℓ q_ℓ(ν)` (the rank-1-equivalent
form of the multi-channel weighting). Each sub-filter must support
`q_weight`. **Not** used by the dissipator path — `dll_lindblad_op_*`
iterates over channels at the operator level and sums per-channel `B_a^{(ℓ)}`.
"""
@inline function q_weight(f::DLLMultiChannelFilter{T}, nu::Real) where {T}
    s = zero(T)
    @inbounds for c in f.channels
        s += T(q_weight(c, nu))
    end
    return s
end

"""
    freq_kernel(filter::DLLMultiChannelFilter, ν) -> sum of per-channel freq_kernel

Diagnostic-only multi-channel `f̂(ν) = Σ_ℓ f̂_ℓ(ν)`. **Not** used to assemble
the dissipator (the dissipator sums per-channel α^{(ℓ)}, which is *not* the
same as building `α` from the multi-channel `f̂`).
"""
@inline function freq_kernel(f::DLLMultiChannelFilter{T}, nu::Real) where {T}
    s = zero(T)
    @inbounds for c in f.channels
        s += T(freq_kernel(c, nu))
    end
    return s
end

"""
    time_kernel(filter::DLLMultiChannelFilter, t) -> sum of per-channel time_kernel

Diagnostic-only multi-channel `f(t) = Σ_ℓ f_ℓ(t)`. The actual time-domain
operators in `dll_lindblad_op_time` / `dll_coherent_op_time` are built
per-channel and summed at the operator level.
"""
@inline function time_kernel(f::DLLMultiChannelFilter{T}, t::Real) where {T}
    CT = Complex{T}
    s = zero(CT)
    @inbounds for c in f.channels
        s += CT(time_kernel(c, t))
    end
    return s
end

"""
    filter_time_cutoff(filter::DLLMultiChannelFilter, tol) -> conservative cutoff

Returns `max_ℓ filter_time_cutoff(channels[ℓ], tol/k)`: distribute the
tolerance budget across channels (so the *sum* of per-channel residuals stays
below `tol`), then take the largest per-channel cutoff. Conservative because
each per-channel cutoff is itself an over-estimate for an asymmetric `f̂`,
and because individual channels could cancel in the sum.
"""
@inline function filter_time_cutoff(f::DLLMultiChannelFilter{T}, tol::Real) where {T}
    k = length(f.channels)
    per_tol = T(tol) / T(k)
    tc = zero(T)
    @inbounds for c in f.channels
        tcc = T(filter_time_cutoff(c, per_tol))
        tcc > tc && (tc = tcc)
    end
    return tc
end

# ---------------------------------------------------------------------------
# Section 2 — ShiftedSymmetricFilter (symmetrised-translate channel, qf-7go.5)
# ---------------------------------------------------------------------------
#
# Each translated channel ℓ wraps a base DLL filter (`DLLGaussianFilter` or
# `DLLMetropolisFilter`) with a symmetric ν-shift `ν_ℓ` and a scalar weight
# `w_ℓ`:
#
#     q_ℓ(ν) = √(w_ℓ / 2) · [q_base(ν − ν_ℓ) + q_base(ν + ν_ℓ)]   (ν_ℓ ≠ 0)
#     q_ℓ(ν) = √(w_ℓ)     · q_base(ν)                              (ν_ℓ = 0)
#
# The two-sided shift makes `q_ℓ` real-even (since `q_base` is real-even),
# so the per-channel Kossakowski α^{(ℓ)} = q_ℓ ⊗ q_ℓ · e^{-β(ν+ν')/4}
# satisfies KMS skew-symmetry (Eq. 4.7) automatically.
#
# Pulling the KMS factor through the shifts gives the exact f̂_ℓ and f_ℓ
# expressions used below:
#
#     q_base(ν − ν_ℓ) · e^{-β ν / 4}
#         = e^{-β ν_ℓ / 4} · q_base(u) · e^{-β u / 4}      (u = ν − ν_ℓ)
#         = e^{-β ν_ℓ / 4} · f̂_base(ν − ν_ℓ),
#
# and similarly `q_base(ν + ν_ℓ) · e^{-β ν / 4} = e^{+β ν_ℓ / 4} · f̂_base(ν + ν_ℓ)`.
# Inverse-FT'ing each shifted f̂_base gives `e^{∓i ν_ℓ t} · f_base(t)`, so
#
#     f_ℓ(t) = √(w_ℓ / 2) · f_base(t) · 2 cosh(β ν_ℓ / 4 + i ν_ℓ t).

"""
    ShiftedSymmetricFilter{T<:AbstractFloat, F<:AbstractFilter}(base, shift, weight, beta)

Symmetrised-translate channel for the multi-rank DLL construction
(qf-7go.5). Wraps a base DLL filter with a symmetric ν-shift `shift`
and a scalar weight `weight` so that
`q_ℓ(ν) = √(weight/2) · [q_base(ν − shift) + q_base(ν + shift)]`
(or `√(weight) · q_base(ν)` when `shift = 0`). `q_ℓ` is real-even, so
the per-channel KMS skew-symmetry follows from `q_base`'s evenness.

Construct via `dll_multichannel_translates(base; centers, weights)`.
"""
struct ShiftedSymmetricFilter{T<:AbstractFloat, F<:AbstractFilter} <: AbstractFilter
    base::F
    shift::T
    weight::T
    beta::T
end

ShiftedSymmetricFilter(base::F, shift::T, weight::T) where
        {T<:AbstractFloat, F<:AbstractFilter} =
    ShiftedSymmetricFilter{T, F}(base, shift, weight, T(base.beta))

Base.eltype(::ShiftedSymmetricFilter{T}) where {T} = Complex{T}

"""
    q_weight(filter::ShiftedSymmetricFilter, ν) -> Real

`q_ℓ(ν) = √(w/2) · [q_base(ν − shift) + q_base(ν + shift)]` for `shift ≠ 0`,
or `√w · q_base(ν)` when `shift = 0`. Real-valued (because `q_base` is
real-valued and the symmetrisation preserves that).
"""
@inline function q_weight(f::ShiftedSymmetricFilter{T}, nu::Real) where {T}
    if iszero(f.shift)
        return T(sqrt(f.weight)) * T(q_weight(f.base, nu))
    end
    qm = T(q_weight(f.base, T(nu) - f.shift))
    qp = T(q_weight(f.base, T(nu) + f.shift))
    return T(sqrt(f.weight / T(2))) * (qm + qp)
end

"""
    freq_kernel(filter::ShiftedSymmetricFilter, ν) -> Real

`f̂_ℓ(ν) = q_ℓ(ν) · e^{-β ν / 4}`. Equivalently,
`f̂_ℓ(ν) = √(w/2) · [e^{-β shift / 4} f̂_base(ν − shift) +
                     e^{+β shift / 4} f̂_base(ν + shift)]`
(real, by `f̂_base` being real-valued for both supported DLL filters).
"""
@inline function freq_kernel(f::ShiftedSymmetricFilter{T}, nu::Real) where {T}
    return q_weight(f, nu) * exp(-f.beta * T(nu) / T(4))
end

"""
    time_kernel(filter::ShiftedSymmetricFilter, t) -> Complex

`f_ℓ(t) = √(w/2) · f_base(t) · 2 cosh(β shift / 4 + i shift t)` for
`shift ≠ 0`, or `√w · f_base(t)` for `shift = 0`. Derived by pulling the
KMS factor through each ν-shifted f̂_base and inverse-FT'ing.
"""
@inline function time_kernel(f::ShiftedSymmetricFilter{T}, t::Real) where {T}
    CT = Complex{T}
    fb = CT(time_kernel(f.base, t))
    if iszero(f.shift)
        return T(sqrt(f.weight)) * fb
    end
    z = Complex{T}(f.beta * f.shift / T(4), f.shift * T(t))
    return T(sqrt(f.weight / T(2))) * fb * (T(2) * cosh(z))
end

"""
    filter_time_cutoff(filter::ShiftedSymmetricFilter, tol) -> Real

The `cosh(β shift / 4 + i shift t)` envelope grows the modulus by at
most `2 · cosh(β shift / 4)` over the bare `|f_base(t)|`. Tighten the
per-base tolerance by that factor (and the `√(w/2) · 2` prefactor) so
that `|time_kernel(filter, t)| ≤ tol` at the returned cutoff.
"""
@inline function filter_time_cutoff(f::ShiftedSymmetricFilter{T}, tol::Real) where {T}
    if iszero(f.shift)
        return T(filter_time_cutoff(f.base, T(tol) / T(sqrt(f.weight))))
    end
    envelope = T(sqrt(f.weight / T(2))) * T(2) * cosh(f.beta * f.shift / T(4))
    return T(filter_time_cutoff(f.base, T(tol) / envelope))
end

# ---------------------------------------------------------------------------
# Section 3 — dll_multichannel_translates factory
# ---------------------------------------------------------------------------

"""
    dll_multichannel_translates(base::AbstractFilter;
                                 centers::AbstractVector{<:Real} = [0.0],
                                 weights::Union{Nothing, AbstractVector{<:Real}} = nothing)
        -> DLLMultiChannelFilter

Build a `DLLMultiChannelFilter` from symmetrised translates of `base`
(qf-7go.5):

    q_ℓ(ν) = √(weights[ℓ]/2) · [q_base(ν − centers[ℓ]) + q_base(ν + centers[ℓ])]
    q_ℓ(ν) = √(weights[ℓ])   · q_base(ν)                                          (centers[ℓ] = 0)

Each `q_ℓ` is real-even by construction, so KMS skew-symmetry of the
multi-channel α follows automatically.

# Arguments

- `base`: base DLL filter (`DLLGaussianFilter` or `DLLMetropolisFilter`).
  Must have a `beta` field.
- `centers`: ν-shifts. `centers[1] = 0` recovers the standard
  single-channel filter as the first channel. Default `[0.0]` (k = 1
  reduction; equivalent to `DLLMultiChannelFilter([base], base.beta)`).
- `weights`: per-channel scalar weights. Default `[1, 1, …, 1]` (uniform,
  no rescaling). The per-channel Kossakowski scales as `α^(ℓ) ∝ w_ℓ`.

# Constraints

For Metropolis-style bases with bump support `[-S, S]`, all centers
must lie inside the flat-top region: `|centers[ℓ]| ≤ base.S / 2`. This
keeps both shifted copies of the bump entirely inside the support
(`q_base(ν ± center) ≠ 0` only on `|ν ± center| ≤ base.S`, so the sum
is a valid Gevrey filter on `[-S - |center|, S + |center|]`). For
Gaussian bases the constraint is vacuous (no compact support).

# Examples

```julia
base = DLLMetropolisFilter(5.0; S = 2.0)
multi = dll_multichannel_translates(base; centers = [0.0, 0.5])
# k = 2 channels; first is the standard DLL Metropolis, second is a
# symmetrised pair of ν-translated copies at ±0.5.
```
"""
function dll_multichannel_translates(
    base::AbstractFilter;
    centers::AbstractVector{<:Real} = [0.0],
    weights::Union{Nothing, AbstractVector{<:Real}} = nothing,
)
    if !hasproperty(base, :beta)
        throw(ArgumentError("base filter $(typeof(base)) lacks a `beta` field — " *
                            "only DLL-style filters supported."))
    end
    if isempty(centers)
        throw(ArgumentError("centers must be non-empty."))
    end
    T = typeof(float(base.beta))
    k = length(centers)
    ws = if weights === nothing
        ones(T, k)
    else
        if length(weights) != k
            throw(ArgumentError("length(weights)=$(length(weights)) must equal " *
                                "length(centers)=$k."))
        end
        T.(weights)
    end
    if any(w -> w <= 0, ws)
        throw(ArgumentError("weights must be strictly positive."))
    end
    if base isa DLLMetropolisFilter
        S = base.S
        for (ℓ, c) in enumerate(centers)
            if abs(T(c)) > S / 2
                throw(ArgumentError("centers[$ℓ]=$c lies outside the bump flat-top " *
                                    "[-S/2, S/2] = [-$(S/2), $(S/2)]."))
            end
        end
    end

    channels = ShiftedSymmetricFilter{T, typeof(base)}[]
    sizehint!(channels, k)
    for ℓ in 1:k
        push!(channels, ShiftedSymmetricFilter{T, typeof(base)}(
            base, T(centers[ℓ]), ws[ℓ], T(base.beta)))
    end
    return DLLMultiChannelFilter{T, ShiftedSymmetricFilter{T, typeof(base)}}(
        channels, T(base.beta))
end

# ---------------------------------------------------------------------------
# Section 4 — Operator-level overloads (Vector{Matrix} dissipator, summed G)
# ---------------------------------------------------------------------------

"""
    dll_lindblad_op_bohr(jump, hamiltonian, filter::DLLMultiChannelFilter)
        -> Vector{Matrix}

Multi-channel Bohr-domain Lindblad operators (qf-7go.2). Returns a length-`k`
vector `[L^(1), …, L^(k)]` where each `L^(ℓ)` is the standard single-channel
Lindblad operator built from `filter.channels[ℓ]`.

# Why a vector and not a sum

The dissipator
`Σ_ℓ L^(ℓ) ρ (L^(ℓ))† − ½ {(L^(ℓ))† L^(ℓ), ρ}` is *not* equal to the
naive single-Lindblad sandwich `(Σ_ℓ L^(ℓ)) ρ (Σ_ℓ L^(ℓ))† − …` — the
cross-terms `L^(i) ρ (L^(j))†` (i ≠ j) appear in the latter but should be
absent (the multi-channel α^multi_{ν,ν'} = Σ_ℓ α^(ℓ)_{ν,ν'} has no cross
contribution between channels). Callers must therefore iterate over the
returned vector and accumulate per-channel dissipator contributions.

`_jump_contribution!` handles this automatically (see DLL-MR.4 wiring).
"""
function dll_lindblad_op_bohr(
    jump::JumpOp,
    hamiltonian::HamHam{T},
    filter::DLLMultiChannelFilter{T},
) where {T<:AbstractFloat}
    return [dll_lindblad_op_bohr(jump, hamiltonian, c) for c in filter.channels]
end

"""
    dll_lindblad_op_time(jump, hamiltonian, time_labels, filter::DLLMultiChannelFilter, t0)
        -> Vector{Matrix}

Multi-channel time-domain Lindblad operators (qf-7go.3). Returns a length-`k`
vector `[L^(1), …, L^(k)]` of single-channel time-domain operators (same
shape as `dll_lindblad_op_bohr`'s multi-channel overload). Caller must
iterate per-channel for the dissipator (cross-terms `L^(i) ρ (L^(j))†` are
absent in the multi-channel α).

Performance: k separate single-channel calls; wall time scales linearly
in `k` (the per-channel quadrature dominates). Channels are independent —
no NUFFT fusion attempted (k = O(1) in the τ_mix sweep).
"""
function dll_lindblad_op_time(
    jump::JumpOp,
    hamiltonian::HamHam{T},
    time_labels::AbstractVector{<:Real},
    filter::DLLMultiChannelFilter{T},
    t0::Real,
) where {T<:AbstractFloat}
    return [dll_lindblad_op_time(jump, hamiltonian, time_labels, c, t0)
            for c in filter.channels]
end

"""
    dll_coherent_op_bohr(jumps, hamiltonian, filter::DLLMultiChannelFilter, beta) -> Matrix

Multi-channel Bohr-domain coherent operator (qf-7go.2). The DLL coherent
kernel ĝ^a(ν, ν') is *linear* in the Kossakowski outer product
`f̂(ν) conj(f̂(ν'))`, so `G^multi = Σ_ℓ G^(ℓ)` with each `G^(ℓ)` built from
the single-channel `dll_coherent_op_bohr`. Returns a single accumulator
matrix (not a vector — unlike the dissipator path, the coherent term has
no cross-channel structure and adds linearly).
"""
function dll_coherent_op_bohr(
    jumps::AbstractVector{<:JumpOp},
    hamiltonian::HamHam{T},
    filter::DLLMultiChannelFilter{T},
    beta::Real,
) where {T<:AbstractFloat}
    G = dll_coherent_op_bohr(jumps, hamiltonian, filter.channels[1], beta)
    @inbounds for ℓ in 2:length(filter.channels)
        G .+= dll_coherent_op_bohr(jumps, hamiltonian, filter.channels[ℓ], beta)
    end
    return G
end

"""
    dll_coherent_op_time(jumps, hamiltonian, time_labels,
                         filter::DLLMultiChannelFilter, beta, τ) -> Matrix

Multi-channel time-domain coherent operator (qf-7go.3). Sums per-channel
coherent operators built via the existing `dll_coherent_op_time` overloads
(closed-form for Gaussian channels, 2D-NUFFT-from-grid for Metropolis
channels). Returns a single accumulator matrix — the coherent kernel
ĝ(ν, ν') is linear in `f̂(ν) f̂(ν')*`, so `G^multi = Σ_ℓ G^(ℓ)`.

Performance: one single-channel call per channel; wall time ≈ k× the
single-channel cost. The per-channel NUFFT plans are not fused (k = O(1)
in the τ_mix sweep).
"""
function dll_coherent_op_time(
    jumps::AbstractVector{<:JumpOp},
    hamiltonian::HamHam{T},
    time_labels::AbstractVector{<:Real},
    filter::DLLMultiChannelFilter{T},
    beta::Real,
    τ::Real;
    kwargs...,
) where {T<:AbstractFloat}
    G = dll_coherent_op_time(jumps, hamiltonian, time_labels,
                              filter.channels[1], beta, τ; kwargs...)
    @inbounds for ℓ in 2:length(filter.channels)
        G .+= dll_coherent_op_time(jumps, hamiltonian, time_labels,
                                    filter.channels[ℓ], beta, τ; kwargs...)
    end
    return G
end

# ---------------------------------------------------------------------------
# Section 5 — Simulator-coupling overloads
# ---------------------------------------------------------------------------

# Multi-channel DLL filter (qf-7go.1): k Lindblad operators per coupling.
# Sums per-channel dissipators — no cross terms in the multi-channel α.
@inline function _accumulate_dll_bohr_dissipator!(
    L_target::AbstractMatrix{<:Complex},
    jump::JumpOp,
    hamiltonian::HamHam,
    filter::DLLMultiChannelFilter,
    ws::Workspace{Lindbladian},
)
    Ls = dll_lindblad_op_bohr(jump, hamiltonian, filter)
    @inbounds for L_a in Ls
        _vectorize_liouv_diss_and_add!(L_target, L_a, 1.0, ws)
    end
    return L_target
end

# Multi-channel DLL filter: TimeDomain OFT-prefactor enumeration. The
# dissipator path then sums `L^(ℓ) ρ (L^(ℓ))† − …` per channel (no cross
# terms in the multi-channel α).
@inline _filter_channels_for_dll_oft(filter::DLLMultiChannelFilter) = filter.channels
