# Mixing-time post-processing for sampled curves and Krylov spectral data.

"""
    MixingTimeEstimate

Result of fitting a sampled distance-to-equilibrium curve.

# Fields
- `fitted_gap`, `amplitude`, `offset`: Slow fitted mode.
- `gap_ci`, `gap_se`, `r_squared`, `converged`: Fit diagnostics.
- `mixing_time`: Selected actual or extrapolated answer.
- `mixing_time_extrapolated`, `mixing_time_actual`, `target_epsilon`: Crossing details.
- `fit_result`, `biexp_fit_result`, `model_used`: Underlying model results.
"""
struct MixingTimeEstimate
    # Fit parameters (from FitResult)
    fitted_gap::Float64
    amplitude::Float64
    offset::Float64
    gap_ci::Tuple{Float64, Float64}
    gap_se::Float64
    r_squared::Float64
    converged::Bool
    # Mixing time results
    mixing_time::Float64
    mixing_time_extrapolated::Union{Nothing, Float64}
    mixing_time_actual::Union{Nothing, Float64}
    target_epsilon::Union{Nothing, Float64}
    # Full fit for advanced users
    fit_result::FitResult
    # Fitted model and optional two-rate details.
    model_used::Symbol
    biexp_fit_result::Union{Nothing, BiexpFitResult}
end

"""
    _check_fit_quality(fit::FitResult, target_epsilon)

Issue `@warn` messages for quality gate violations. Does not throw.
"""
function _check_fit_quality(fit::FitResult, target_epsilon::Union{Nothing, Float64})
    if fit.r_squared < 0.95
        @warn "Fit R-squared = $(fit.r_squared) < 0.95. Single-exponential model may not describe the data well."
    end
    if target_epsilon !== nothing && fit.offset > 0.1 * target_epsilon
        @warn "Fit offset C = $(fit.offset) is large relative to target epsilon = $(target_epsilon). Extrapolation may be unreliable."
    end
    if !fit.converged
        @warn "Levenberg-Marquardt optimization did not converge. Fit results are unreliable."
    end
    if fit.gap_se > 0.5 * fit.gap && isfinite(fit.gap_se)
        @warn "Gap standard error ($(fit.gap_se)) exceeds 50% of fitted gap ($(fit.gap)). Confidence interval is very wide."
    end
    return nothing
end

"""
    _find_actual_mixing_time(times, dists, target_epsilon)

Return the first sampled threshold crossing, or `nothing` if none exists.
"""
function _find_actual_mixing_time(
    times::AbstractVector{<:Real},
    dists::AbstractVector{<:Real},
    target_epsilon::Union{Nothing, Float64},
)
    target_epsilon === nothing && return nothing
    for i in eachindex(dists)
        if dists[i] <= target_epsilon
            return Float64(times[i])
        end
    end
    return nothing  # target not reached in data
end

"""
    _extrapolate_mixing_time(fit::FitResult, target_epsilon)

Solve the single-exponential fit for a target crossing.

Return `nothing` when the parameters do not define a future crossing.
"""
function _extrapolate_mixing_time(fit::FitResult, target_epsilon::Union{Nothing, Float64})
    target_epsilon === nothing && return nothing
    fit.gap <= 0.0 && return nothing
    fit.amplitude <= 0.0 && return nothing

    # Math: $t = -log((epsilon-C)/A) / Delta$.
    effective_target = target_epsilon - fit.offset
    effective_target <= 0.0 && return nothing   # offset exceeds target
    effective_target >= fit.amplitude && return nothing  # already below target at t=0

    return -log(effective_target / fit.amplitude) / fit.gap
end

"""
    _extrapolate_mixing_time_biexp(bifit::BiexpFitResult, target_epsilon)

Solve the bi-exponential fit for a target crossing by bisection.
"""
function _extrapolate_mixing_time_biexp(bifit::BiexpFitResult, target_epsilon::Union{Nothing, Float64})
    target_epsilon === nothing && return nothing
    bifit.gap <= 0.0 && return nothing

    # Check: f(0) = A1 + A2 + C; need f(0) > target for a crossing to exist
    f0 = bifit.amplitude_fast + bifit.amplitude + bifit.offset
    f0 <= target_epsilon && return nothing  # already below target at t=0

    # Check: asymptotic value C must be below target
    bifit.offset >= target_epsilon && return nothing

    # Define the function to find root of: f(t) - epsilon = 0
    function biexp_residual(t)
        return bifit.amplitude_fast * exp(-bifit.gap_fast * t) +
               bifit.amplitude * exp(-bifit.gap * t) +
               bifit.offset - target_epsilon
    end

    # Upper bracket: use slow-mode estimate with 3x safety margin
    # From slow mode alone: t_slow = -ln((eps - C) / A_slow) / gap_slow
    eff = target_epsilon - bifit.offset
    if eff <= 0.0
        return nothing
    end
    t_slow_est = if bifit.amplitude > 0.0 && eff < bifit.amplitude
        -log(eff / bifit.amplitude) / bifit.gap
    else
        # Fallback: use a large upper bracket
        100.0 / bifit.gap
    end
    t_upper = max(t_slow_est * 3.0, 10.0 / bifit.gap)

    # Ensure f(t_upper) < target (bracket is valid)
    if biexp_residual(t_upper) > 0.0
        # Expand bracket
        t_upper *= 3.0
        if biexp_residual(t_upper) > 0.0
            return nothing  # cannot bracket
        end
    end

    try
        t_mix = Roots.find_zero(biexp_residual, (0.0, t_upper), Roots.Bisection())
        return t_mix
    catch
        return nothing
    end
end

"""
    _biexp_to_single_fit_result(bifit::BiexpFitResult) -> FitResult

Project a `BiexpFitResult` onto its slow-mode `FitResult` fields.
"""
function _biexp_to_single_fit_result(bifit::BiexpFitResult)
    return FitResult(
        bifit.gap,          # gap (slow mode)
        bifit.amplitude,    # amplitude (slow mode)
        bifit.offset,       # offset
        bifit.gap_ci,       # gap_ci (slow mode)
        bifit.gap_se,       # gap_se (slow mode)
        bifit.r_squared,    # r_squared (from full bi-exp fit)
        bifit.converged,    # converged
        bifit.residuals,    # residuals (from full bi-exp fit)
        bifit.times_used,   # times_used
        bifit.values_used,  # values_used
    )
end

"""
    estimate_mixing_time(times, distances; kwargs...) -> MixingTimeEstimate

Fit a decay curve and optionally extrapolate its target crossing.

# Arguments
- `times`: Sample times.
- `distances`: Non-negative distance or observable values.

# Keywords
- `skip_initial`: Fraction of leading samples to discard.
- `target_epsilon`: Target value; required for extrapolation.
- `extrapolate`: Solve the fitted model rather than use a sampled crossing.
- `level`: Confidence level for the gap interval.
- `model`: `:single` or `:biexp`.

# Returns
A [`MixingTimeEstimate`](@ref) with crossing and fit diagnostics.
"""
function estimate_mixing_time(
    times::AbstractVector{<:Real},
    distances::AbstractVector{<:Real};
    skip_initial::Real = 0.2,
    target_epsilon::Union{Nothing, Real} = nothing,
    extrapolate::Bool = false,
    level::Real = 0.95,
    model::Symbol = :biexp,
)::MixingTimeEstimate
    length(times) == length(distances) || throw(ArgumentError(
        "times and distances must have the same length (got $(length(times)) and $(length(distances)))"))
    length(times) >= 10 || throw(ArgumentError(
        "Need at least 10 data points for mixing time estimation (got $(length(times)))"))

    if extrapolate && target_epsilon === nothing
        throw(ArgumentError("target_epsilon required when extrapolate=true"))
    end

    model in (:single, :biexp) || throw(ArgumentError(
        "model must be :single or :biexp (got :$model)"))

    skip_initial_f = Float64(skip_initial)
    level_f        = Float64(level)
    target_eps_f   = target_epsilon === nothing ? nothing : Float64(target_epsilon)

    t_mix_actual = _find_actual_mixing_time(times, distances, target_eps_f)

    if model == :single
        fit = fit_exponential_decay(Float64.(times), Float64.(distances);
            skip_initial=skip_initial_f, level=level_f)

        _check_fit_quality(fit, target_eps_f)

        t_mix_extrap = extrapolate ? _extrapolate_mixing_time(fit, target_eps_f) : nothing

        mixing_time = if extrapolate
            t_mix_extrap !== nothing ? t_mix_extrap : NaN
        elseif target_eps_f !== nothing
            t_mix_actual !== nothing ? t_mix_actual : NaN
        else
            Float64(last(times))
        end

        return MixingTimeEstimate(
            fit.gap, fit.amplitude, fit.offset,
            fit.gap_ci, fit.gap_se, fit.r_squared, fit.converged,
            mixing_time, t_mix_extrap, t_mix_actual, target_eps_f,
            fit, :single, nothing,
        )

    else  # model == :biexp
        bifit = fit_biexponential_decay(Float64.(times), Float64.(distances);
            skip_initial=skip_initial_f, level=level_f)

        t_mix_extrap = extrapolate ? _extrapolate_mixing_time_biexp(bifit, target_eps_f) : nothing

        mixing_time = if extrapolate
            t_mix_extrap !== nothing ? t_mix_extrap : NaN
        elseif target_eps_f !== nothing
            t_mix_actual !== nothing ? t_mix_actual : NaN
        else
            Float64(last(times))
        end

        synthetic_fit = _biexp_to_single_fit_result(bifit)

        return MixingTimeEstimate(
            bifit.gap, bifit.amplitude, bifit.offset,
            bifit.gap_ci, bifit.gap_se, bifit.r_squared, bifit.converged,
            mixing_time, t_mix_extrap, t_mix_actual, target_eps_f,
            synthetic_fit, :biexp, bifit,
        )
    end
end

"""
    estimate_mixing_time(result::ThermalizeResults; kwargs...) -> MixingTimeEstimate

Estimate mixing from a completed full-density-matrix simulation.

# Arguments
- `result`: Simulation result containing sample times and trace distances.

# Keywords
- `skip_initial`, `target_epsilon`, `extrapolate`, `level`: As in the vector method.
- `model`: Fit model; defaults to `:single` for this overload.

# Returns
A [`MixingTimeEstimate`](@ref). This is post-processing and does not rerun the simulation.
"""
function estimate_mixing_time(
    result::ThermalizeResults;
    skip_initial::Real = 0.2,
    target_epsilon::Union{Nothing, Real} = nothing,
    extrapolate::Bool = false,
    level::Real = 0.95,
    model::Symbol = :single,
)::MixingTimeEstimate
    return estimate_mixing_time(
        result.time_steps, result.trace_distances;
        skip_initial    = skip_initial,
        target_epsilon  = target_epsilon,
        extrapolate     = extrapolate,
        level           = level,
        model           = model,
    )
end

"""
    estimate_mixing_time(integrator_result::NamedTuple; kwargs...) -> MixingTimeEstimate

Estimate mixing from an integrator result containing `t` and `distances`.

Krylov predictor results are rejected because their spectral data should be
passed to [`eigenmode_mixing_time`](@ref) instead of curve fitting.
"""
function estimate_mixing_time(integrator_result::NamedTuple; kwargs...)::MixingTimeEstimate
    # Predictor spectral data require exact subspace bisection, not a curve fit.
    if haskey(integrator_result, :R_modes) || haskey(integrator_result, :eigenvalues)
        throw(ArgumentError(
            "estimate_mixing_time received a Krylov trajectory-predictor result " *
            "(it carries the spectral fields :R_modes / :eigenvalues). The " *
            "(bi)exponential curve fit must not be used on " *
            "predict_lindbladian_trajectory / predict_channel_trajectory output — " *
            "call `eigenmode_mixing_time(traj, target_epsilon)` (exact closed-form " *
            "bisection on the Krylov spectrum) instead. The curve-fit estimator " *
            "remains the right tool for the matrix-free integrator outputs " *
            "(lindblad_action_integrate / discriminant_action_integrate) and " *
            "ThermalizeResults."))
    end
    return estimate_mixing_time(
        integrator_result.t, integrator_result.distances;
        kwargs...,
    )
end

# Exact mixing-time bisection on a captured Krylov spectral expansion.

const _EIGENMODE_ZERO_TOL = 1e-10  # |λ| below this counts as the steady mode

"""
    eigenmode_mixing_time(eigenvalues, c, R_modes, rho_inf, sigma_beta,
                          target_epsilon; t_upper, atol, max_iters,
                          eigenvalue_zero_tol)
        -> NamedTuple

Bisect the trace-distance crossing of a biorthogonal Krylov expansion.

# Arguments
- `eigenvalues`: Lindbladian rates, including the steady mode.
- `c`: Biorthogonal expansion coefficients.
- `R_modes`: Right eigenmodes as matrices.
- `rho_inf`: Captured stationary state.
- `sigma_beta`: Reference Gibbs state in the same basis.
- `target_epsilon`: Trace-distance threshold.

# Keywords
- `t_upper`: Upper bracket; zero selects `max(10/gap, 50)`.
- `atol`: Absolute time-axis tolerance.
- `max_iters`: Maximum bisection iterations.
- `eigenvalue_zero_tol`: Magnitude below which a rate is stationary.

# Returns
A named tuple with `mixing_time`, `gap`, `floor_distance`, `source`, and
`n_evals`. `mixing_time` is `Inf` when the target is below the asymptotic
floor and `NaN` for degenerate input or a failed bracket.
"""
function eigenmode_mixing_time(
    eigenvalues::AbstractVector{<:Complex},
    c::AbstractVector{<:Complex},
    R_modes::AbstractVector{<:AbstractMatrix},
    rho_inf::AbstractMatrix,
    sigma_beta::AbstractMatrix,
    target_epsilon::Real;
    t_upper::Real = 0.0,
    atol::Real = 1e-3,
    max_iters::Int = 64,
    eigenvalue_zero_tol::Real = _EIGENMODE_ZERO_TOL,
)::NamedTuple
    h = length(eigenvalues)
    h == length(c) == length(R_modes) || throw(ArgumentError(
        "eigenvalues, c, R_modes must have the same length (got $h, $(length(c)), $(length(R_modes)))"))
    target_epsilon > 0.0 || throw(ArgumentError(
        "target_epsilon must be positive (got $target_epsilon)"))

    # Math: $d(infinity) = norm(rho_infinity-sigma_beta)_1 / 2$.
    floor_distance = sum(svdvals(rho_inf .- sigma_beta)) / 2

    # Math: $Delta = min_(lambda_i != 0) abs(Re(lambda_i))$.
    gap = Inf
    for i in 1:h
        abs(eigenvalues[i]) < eigenvalue_zero_tol && continue
        gi = abs(real(eigenvalues[i]))
        gi < gap && (gap = gi)
    end

    # Degenerate input: no non-steady modes captured.
    if !isfinite(gap) || gap <= 0.0
        return (
            mixing_time     = NaN,
            gap             = isfinite(gap) ? gap : NaN,
            floor_distance  = floor_distance,
            source          = :nan,
            n_evals         = 0,
        )
    end

    # Floor branch: target below asymptotic floor → no crossing.
    if floor_distance >= target_epsilon
        return (
            mixing_time     = Inf,
            gap             = gap,
            floor_distance  = floor_distance,
            source          = :floor,
            n_evals         = 0,
        )
    end

    d = size(rho_inf, 1)
    T = promote_type(eltype(rho_inf), eltype(sigma_beta), ComplexF64)
    rho_t = Matrix{T}(undef, d, d)
    floor_residual = Matrix{T}(rho_inf .- sigma_beta)
    eval_count = Ref(0)

    function d_at(t::Real)::Float64
        copyto!(rho_t, floor_residual)
        @inbounds for i in 1:h
            abs(eigenvalues[i]) < eigenvalue_zero_tol && continue
            phase = exp(eigenvalues[i] * t)
            rho_t .+= (c[i] * phase) .* R_modes[i]
        end
        # Defensive Hermitisation (mirrors predict_*_trajectory loop).
        @inbounds for j in 1:d, k in 1:d
            rho_t[k, j] = (rho_t[k, j] + conj(rho_t[j, k])) / 2
        end
        eval_count[] += 1
        return sum(svdvals(rho_t)) / 2
    end

    # Bisection bracket. The eigenmode formula is monotonically decreasing in
    # t on the slow tail, so [0, t_upper] with d(t_upper) < ε brackets a root.
    t_hi = float(t_upper)
    if t_hi <= 0.0
        t_hi = max(10.0 / gap, 50.0)
    end

    d_hi = d_at(t_hi)
    if d_hi > target_epsilon
        # One 3× expansion. The slow mode decays like e^{-gap * t}, so
        # log(d(t)/d(t_hi)) ≈ -gap (t - t_hi); 3× is enough margin in
        # practice for any well-behaved spectrum.
        t_hi *= 3.0
        d_hi = d_at(t_hi)
        if d_hi > target_epsilon
            return (
                mixing_time     = NaN,
                gap             = gap,
                floor_distance  = floor_distance,
                source          = :nan,
                n_evals         = eval_count[],
            )
        end
    end

    # d(0) ≥ target_epsilon must hold for a crossing on [0, t_hi]; if not,
    # the trajectory is below target at t=0 — degenerate (or trivially mixed).
    d_zero = d_at(0.0)
    if d_zero <= target_epsilon
        return (
            mixing_time     = 0.0,
            gap             = gap,
            floor_distance  = floor_distance,
            source          = :extrapolated,
            n_evals         = eval_count[],
        )
    end

    residual(t::Float64) = d_at(t) - float(target_epsilon)

    t_mix = try
        # `atol` is the t-axis tolerance ⇒ pass as `xatol` to Roots (their
        # `atol` kw governs f(x)=d(t)-ε, which would couple to the slow
        # mode's amplitude and hide stride O(ε / |d'|) ≫ atol).
        Roots.find_zero(residual, (0.0, t_hi), Roots.Bisection();
                         xatol=float(atol), maxiters=max_iters)
    catch
        return (
            mixing_time     = NaN,
            gap             = gap,
            floor_distance  = floor_distance,
            source          = :nan,
            n_evals         = eval_count[],
        )
    end

    return (
        mixing_time     = float(t_mix),
        gap             = gap,
        floor_distance  = floor_distance,
        source          = :extrapolated,
        n_evals         = eval_count[],
    )
end

"""
    eigenmode_mixing_time(traj::NamedTuple, target_epsilon; kwargs...) -> NamedTuple

Compute mixing time directly from a Krylov predictor result.

Channel eigenvalues are converted to rates as `log(mu)/delta`; Lindbladian
eigenvalues are already rates. Remaining keywords forward to the spectral-data
method.
"""
function eigenmode_mixing_time(traj::NamedTuple, target_epsilon::Real; kwargs...)::NamedTuple
    for f in (:eigenvalues, :c, :R_modes, :rho_inf, :sigma_beta)
        haskey(traj, f) || throw(ArgumentError(
            "eigenmode_mixing_time(traj::NamedTuple, …) expects a Krylov " *
            "trajectory-predictor result with fields :eigenvalues, :c, :R_modes, " *
            ":rho_inf, :sigma_beta (missing :$f). Pass the output of " *
            "predict_lindbladian_trajectory / predict_channel_trajectory."))
    end
    eigenvalues = if haskey(traj, :delta_used)
        # Math: $lambda_eff = log(mu)/delta$, so $exp(lambda_eff k delta) = mu^k$.
        δ = traj.delta_used
        ComplexF64[log(complex(μ)) / δ for μ in traj.eigenvalues]
    else
        ComplexF64.(traj.eigenvalues)
    end
    return eigenmode_mixing_time(
        eigenvalues, traj.c, traj.R_modes, traj.rho_inf, traj.sigma_beta,
        target_epsilon; kwargs...)
end
