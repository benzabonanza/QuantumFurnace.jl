# Multi-channel DLL diagnostics. Production simulations normally use a
# single `DLLGaussianFilter` or `DLLMetropolisFilter`.
# Math: $alpha_(nu,nu') = sum_l f_l(nu) conj(f_l(nu'))$, so the dissipator
# must sum channel contributions without cross terms. All channels share beta.

"""
    DLLMultiChannelFilter{T<:AbstractFloat, F<:AbstractFilter}(channels, beta)

Group single-channel DLL filters with a common inverse temperature.

# Fields
- `channels`: Non-empty vector of filters with inverse temperature `beta`.
- `beta`: Shared algorithm-frame inverse temperature.

The scalar kernel methods are diagnostic sums. Operator constructors instead
sum per-channel dissipators to avoid spurious cross terms.
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
    q_weight(filter::DLLMultiChannelFilter, nu) -> Real

Return the diagnostic sum of the channel weights.
"""
@inline function q_weight(f::DLLMultiChannelFilter{T}, nu::Real) where {T}
    s = zero(T)
    @inbounds for c in f.channels
        s += T(q_weight(c, nu))
    end
    return s
end

"""
    freq_kernel(filter::DLLMultiChannelFilter, nu) -> Real

Return the diagnostic sum of the channel frequency kernels.
"""
@inline function freq_kernel(f::DLLMultiChannelFilter{T}, nu::Real) where {T}
    s = zero(T)
    @inbounds for c in f.channels
        s += T(freq_kernel(c, nu))
    end
    return s
end

"""
    time_kernel(filter::DLLMultiChannelFilter, t) -> Complex

Return the diagnostic sum of the channel time kernels.
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

Return the largest channel cutoff after dividing `tol` across channels.
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

# Symmetric frequency translates preserve the real-even DLL weight.
# Math: $q_l(nu) = sqrt(w_l/2) [q(nu-nu_l) + q(nu+nu_l)]$.

"""
    ShiftedSymmetricFilter{T<:AbstractFloat, F<:AbstractFilter}(base, shift, weight, beta)

Represent one weighted, symmetric frequency translate of a DLL filter.

# Fields
- `base`: Unshifted DLL filter.
- `shift`: Non-negative centre frequency.
- `weight`: Positive channel weight.
- `beta`: Algorithm-frame inverse temperature inherited from `base`.
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

Evaluate the weighted symmetric DLL channel weight.
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

Evaluate the shifted frequency kernel with its KMS factor.
"""
@inline function freq_kernel(f::ShiftedSymmetricFilter{T}, nu::Real) where {T}
    return q_weight(f, nu) * exp(-f.beta * T(nu) / T(4))
end

"""
    time_kernel(filter::ShiftedSymmetricFilter, t) -> Complex

Evaluate the shifted time kernel.
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

Return a cutoff that includes the shifted kernel's `cosh` envelope.
"""
@inline function filter_time_cutoff(f::ShiftedSymmetricFilter{T}, tol::Real) where {T}
    if iszero(f.shift)
        return T(filter_time_cutoff(f.base, T(tol) / T(sqrt(f.weight))))
    end
    envelope = T(sqrt(f.weight / T(2))) * T(2) * cosh(f.beta * f.shift / T(4))
    return T(filter_time_cutoff(f.base, T(tol) / envelope))
end

"""
    dll_multichannel_translates(base::AbstractFilter;
                                 centers::AbstractVector{<:Real} = [0.0],
                                 weights::Union{Nothing, AbstractVector{<:Real}} = nothing)
        -> DLLMultiChannelFilter

Build a multi-channel DLL filter from symmetric translates of `base`.

# Arguments

- `base`: DLL filter with a `beta` field.
- `centers`: Frequency shifts; zero reproduces the base channel.
- `weights`: Positive channel weights; defaults to one per centre.

# Constraints

For a Metropolis base, each centre must satisfy `abs(center) <= S/2`.
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

"""
    dll_lindblad_op_bohr(jump, hamiltonian, filter::DLLMultiChannelFilter)
        -> Vector{Matrix}

Construct one Bohr-domain Lindblad operator per DLL channel.

# Returns
A vector of operators. Callers must accumulate their dissipators separately;
summing the operators first would introduce cross terms.
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

Construct one time-domain Lindblad operator per DLL channel.

# Returns
A vector of operators whose dissipators must be accumulated separately.
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

Construct the sum of the channels' Bohr-domain coherent operators.
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

Construct the sum of the channels' time-domain coherent operators.
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

# Sum per-channel dissipators; the multi-channel Kossakowski matrix has no
# cross-channel terms.
@inline function _accumulate_dll_bohr_dissipator!(
    L_target::AbstractMatrix{<:Complex},
    jump::JumpOp,
    hamiltonian::HamHam,
    filter::DLLMultiChannelFilter,
    ws::DenseLindbladianWorkspace,
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
