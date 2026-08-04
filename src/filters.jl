# Operator Fourier-transform filters.
# Math: $f(t) = 1/(2 pi) integral hat(f)(nu) exp(-i t nu) dif nu$.

"""
    AbstractFilter

Supertype for time/frequency kernels used by the operator Fourier transform.
"""
abstract type AbstractFilter end

"""
    GaussianFilter{T<:AbstractFloat}(sigma::T)

CKG Gaussian filter with energy width `sigma`.

Math: `\$f(t) = exp(-sigma^2 t^2)\$` and
`\$hat(f)(nu) prop exp(-nu^2/(4 sigma^2))\$`.
"""
struct GaussianFilter{T<:AbstractFloat} <: AbstractFilter
    sigma::T
end

"""
    DLLGaussianFilter{T<:AbstractFloat}(beta::T)

DLL Gaussian-type filter at inverse temperature `beta`.

The frequency kernel includes the KMS factor:
`\$hat(f)(nu) = exp(1/8) exp(-(beta nu + 1)^2/8)\$`. This numerical variant uses
`w = 1`, so it does not have the compact support required by the rigorous
Paley–Wiener bound.
"""
struct DLLGaussianFilter{T<:AbstractFloat} <: AbstractFilter
    beta::T
end

# NUFFT prefactors store every kernel as complex data.
Base.eltype(::GaussianFilter{T}) where {T} = Complex{T}
Base.eltype(::DLLGaussianFilter{T}) where {T} = Complex{T}

function time_kernel end
function freq_kernel end
function filter_time_cutoff end

"""
    time_kernel(filter::GaussianFilter, t) -> Real

Evaluate `\$exp(-sigma^2 t^2)\$` with deterministic operation ordering.
"""
@inline time_kernel(f::GaussianFilter{T}, t::Real) where {T} =
    exp(-(f.sigma^2) * t^2)

"""
    freq_kernel(filter::GaussianFilter, ν)

Evaluate the unnormalised Gaussian `\$exp(-nu^2/(4 sigma^2))\$`.
"""
@inline freq_kernel(f::GaussianFilter{T}, nu::Real) where {T} =
    exp(-nu^2 / (4 * f.sigma^2))

"""
    _time_oft_prefactor_gaussian(filter::GaussianFilter)

Return the Gaussian OFT normalisation used for time-grid truncation.
"""
@inline _time_oft_prefactor_gaussian(f::GaussianFilter{T}) where {T} =
    sqrt(f.sigma * sqrt(T(2) / T(pi)) / (2 * T(pi)))

"""
    filter_time_cutoff(filter::GaussianFilter, tol)

Return a time cutoff whose Gaussian tail is below `tol`.
"""
@inline filter_time_cutoff(f::GaussianFilter{T}, tol::Real) where {T} =
    sqrt(log(_time_oft_prefactor_gaussian(f) / tol)) / f.sigma

"""
    q_weight(filter::DLLGaussianFilter, nu) -> Real

Evaluate the DLL balance weight `\$q(nu) = exp(-(beta nu)^2/8)\$`.
"""
@inline q_weight(f::DLLGaussianFilter{T}, nu::Real) where {T} =
    exp(-(f.beta * nu)^2 / 8)

"""
    freq_kernel(filter::DLLGaussianFilter, ν)

Evaluate the full DLL frequency kernel, including the KMS factor.
"""
@inline freq_kernel(f::DLLGaussianFilter{T}, nu::Real) where {T} =
    exp(T(1) / 8) * exp(-(f.beta * nu + 1)^2 / 8)

"""
    time_kernel(filter::DLLGaussianFilter, t)

Evaluate the complex closed-form inverse Fourier transform.

Math: `\$f(t) = exp(1/8) sqrt(2/pi)/beta exp(-2t^2/beta^2 + i t/beta)\$`.
"""
@inline function time_kernel(f::DLLGaussianFilter{T}, t::Real) where {T}
    pref = exp(T(1) / 8) * sqrt(T(2) / T(pi)) / f.beta
    decay = exp(-2 * t^2 / f.beta^2)
    phase = cis(t / f.beta)
    return pref * decay * phase
end

"""
    _time_oft_prefactor_dll(filter::DLLGaussianFilter)

Return `|f(0)|` for DLL time-grid truncation.
"""
@inline _time_oft_prefactor_dll(f::DLLGaussianFilter{T}) where {T} =
    exp(T(1) / 8) * sqrt(T(2) / T(pi)) / f.beta

"""
    filter_time_cutoff(filter::DLLGaussianFilter, tol)

Return a time cutoff whose DLL Gaussian tail is below `tol`.
"""
@inline filter_time_cutoff(f::DLLGaussianFilter{T}, tol::Real) where {T} =
    (f.beta / sqrt(T(2))) * sqrt(log(_time_oft_prefactor_dll(f) / tol))

# Hörmander bump used to compactly support the DLL Metropolis filter.
# Math: $eta(t) = exp(-1/t)$ for $t > 0$ and
# $w(x) = phi(2(1-abs(x)))$, with a flat top on $abs(x) <= 1/2$.

"""
    _hormander_eta(t)

Evaluate the one-sided smooth bump factor `eta`.
"""
@inline function _hormander_eta(t::T) where {T<:AbstractFloat}
    return t > zero(T) ? exp(-one(T) / t) : zero(T)
end

"""
    _hormander_phi(t)

Evaluate the smooth step from zero to one on `[0, 1]`.
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

Evaluate the even compact bump with unit plateau on `[-1/2, 1/2]`.
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

# Promote non-floating inputs before evaluating the bump.
@inline _hormander_bump(x::Real) = _hormander_bump(float(x))

"""
    DLLMetropolisFilter{T<:AbstractFloat}(beta::T; S::Real = T(2))

Compactly supported DLL Metropolis-type filter.

# Arguments
- `beta`: Inverse temperature.

# Keywords
- `S`: Support radius; the bump is flat on `abs(nu) <= S/2` and zero for
  `abs(nu) >= S`.

The caller must ensure `\$S/2 >= max abs(nu_BH)\$` when the whole Bohr spectrum
should see the Metropolis plateau.
"""
struct DLLMetropolisFilter{T<:AbstractFloat} <: AbstractFilter
    beta::T
    S::T
end

DLLMetropolisFilter(beta::T; S::Real = T(2)) where {T<:AbstractFloat} =
    DLLMetropolisFilter{T}(beta, T(S))

Base.eltype(::DLLMetropolisFilter{T}) where {T} = Complex{T}

"""
    q_weight(filter::DLLMetropolisFilter, ν)

Evaluate the compactly supported DLL balance weight `q(nu)`.
"""
@inline q_weight(f::DLLMetropolisFilter{T}, nu::Real) where {T} =
    exp(-sqrt(one(T) + (f.beta * nu)^2) / 4) * _hormander_bump(nu / f.S)

"""
    freq_kernel(filter::DLLMetropolisFilter, ν)

Evaluate `\$hat(f)(nu) = q(nu) exp(-beta nu/4)\$`.

On the flat top this approaches one for downward transitions and
`\$exp(-beta nu/2)\$` for upward transitions; it is zero outside `[-S, S]`.
"""
@inline freq_kernel(f::DLLMetropolisFilter{T}, nu::Real) where {T} =
    q_weight(f, nu) * exp(-f.beta * nu / 4)

"""
    time_kernel(filter::DLLMetropolisFilter, t)

Numerically evaluate the compact-support inverse Fourier transform.

# Returns
A `Complex{T}` value. Quadrature uses relative tolerance `1e-12` and absolute
floor `64eps(T)` to control cancellation at large `abs(t)`.
"""
function time_kernel(f::DLLMetropolisFilter{T}, t::Real) where {T}
    # Math: $f(t) = 1/(2 pi) integral_(-S)^S hat(f)(nu) exp(-i t nu) dif nu$.
    integrand = nu -> freq_kernel(f, nu) * cis(-t * nu)
    val, _ = quadgk(integrand, -f.S, f.S;
                    rtol = T(1e-12), atol = T(64) * eps(T))
    return Complex{T}(val / T(2 * π))
end

"""
    filter_time_cutoff(filter::DLLMetropolisFilter, tol)

Return a conservative time cutoff for the compact-support kernel.

A doubling search samples four points per oscillation window and returns the
first window below `tol`, failing after 30 doublings.
"""
function filter_time_cutoff(f::DLLMetropolisFilter{T}, tol::Real) where {T}
    osc_period = T(4 * π) / f.S
    t = max(T(8) * f.beta, T(20))
    iter = 0
    while true
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

# Multi-channel filter types live in `dll_multichannel.jl`.
