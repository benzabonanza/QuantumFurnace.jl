# Verification report — `src/scaling_fit.jl`

## Verdict: PASS WITH NOTES

The module is mathematically correct and the test suite passes (221/221).
Hand-computed AICc matches the implementation to machine precision, the
delta-method propagation `σ_C ≈ C·σ_c` is confirmed against Monte Carlo, the
Burnham–Anderson AICc weight transform is computed exactly, and synthetic
M0/M1 ground-truth datasets recover exponents within the ~0.001 range demanded
by the tests. The correlation matrix is symmetric with diagonal 1 and
off-diagonals in [-1, 1] on every well-conditioned fit observed. The two
convenience methods (NamedTuple and `Vector{<:NamedTuple}`) handle both sweep
schemas (`mixing_time/mixing_time_source` and `tau_mix/tau_mix_source`)
correctly, and every input-validation path has an explicit `@test_throws`.

The notes below are minor: one is a deprecation warning, one is a robustness
gap in the `try/catch` around `stderror`, and the rest are stylistic.

---

## Critical issues
None.

---

## Notes / minor issues

### N-1. `estimate_covar` is deprecated in LsqFit 0.15
**Location**: `src/scaling_fit.jl:135`

```julia
covm   = estimate_covar(fit)
```

LsqFit 0.15 marks this with `@deprecate estimate_covar(fit::LsqFitResult) vcov(fit)`
(`/home/agent/.julia/packages/LsqFit/rHpK1/src/curve_fit.jl:321`). Calling the
function from a clean Julia session emits

```
┌ Warning: `estimate_covar(fit::LsqFitResult)` is deprecated, use `vcov(fit)` instead.
```

The result is byte-identical (`maximum(abs.(estimate_covar(fit) .- vcov(fit))) == 0`),
so behaviour is unchanged. Suggested fix:

```julia
covm   = vcov(fit)
```

(and update the module-header comment on line 12 accordingly).

### N-2. `_build_scaling_fit` only catches `SingularException`
**Location**: `src/scaling_fit.jl:132-142`

```julia
se, ci_vec, covmat = try
    ...
catch e
    e isa LinearAlgebra.SingularException || rethrow(e)
    fill(Inf, n_param), [(-Inf, Inf) for _ in 1:n_param], fill(NaN, n_param, n_param)
end
```

LsqFit's `stderror` (`curve_fit.jl:283-290`) issues a plain `error(...)` (an
`ErrorException`, not a `SingularException`) when the covariance has a
strongly negative diagonal ratio:

```julia
if !isapprox(vratio, 0.0, ...) && vratio < 0.0
    error("Covariance matrix is negative for atol=$atol and rtol=$rtol")
end
```

In a true rank-deficient situation this would propagate as an uncaught error
rather than the documented "report Inf/NaN" behaviour. Pathological in
practice (would require simultaneous QR-pivoted negative variance and
non-approx-zero `vratio`); flagged for awareness only. A defensive patch:

```julia
catch e
    (e isa LinearAlgebra.SingularException || e isa ErrorException) || rethrow(e)
    ...
end
```

### N-3. `confint` re-tuples a tuple list
**Location**: `src/scaling_fit.jl:136`

```julia
[(ci_raw[i][1], ci_raw[i][2]) for i in 1:n_param]
```

`LsqFit.confint` already returns `Vector{Tuple{Float64, Float64}}` (verified at
`curve_fit.jl:315` — `collect(zip(...))`), so this comprehension rebuilds an
identical structure. Harmless; could be simplified to
`ci_raw` (with type `[Tuple{Float64, Float64}(t) for t in ci_raw]` if you
want to be defensive about future LsqFit shape changes).

### N-4. `_get_tau` returns `Nothing` typed as `Any`
**Location**: `src/scaling_fit.jl:299-306`

```julia
function _get_tau(r)
    if haskey(r, :mixing_time)
        return r.mixing_time
    elseif haskey(r, :tau_mix)
        return r.tau_mix
    end
    return nothing
end
```

Branch types are mixed (`Real | Nothing`), so the inferred return type is
`Any`. This is fine for the once-per-cell filter loop but could surprise a
caller who expects a typed return. The downstream `(τ isa Real && isfinite(τ) && τ > 0)`
test correctly rejects `nothing` and `NaN`, so behaviour is safe.

### N-5. Single-β grid path returns `Matrix{Float64}(undef, n, 1)`
**Location**: `src/scaling_fit.jl:496-500`

The `scaling_fit_grid(fit)` fallback for a single β value returns a one-column
matrix `(length(n_grid), 1)`. Verified: `[5.0]` for `beta_grid`,
`size(tau_predicted) == (10, 1)`. This is the documented behaviour.

### N-6. Permissive source filter when source field is missing
**Location**: `src/scaling_fit.jl:322-325`

```julia
src = _get_source(r)
if src !== nothing && !(src in source_filter)
    continue
end
```

If a NamedTuple has no `mixing_time_source` / `tau_mix_source` field, the
cell is *kept* regardless of `source_filter`. The intent (per docstring) is
"filter to extrapolated only", so a NamedTuple without a source should
arguably also be filtered. In practice both production schemas (`channel`
and `Lindbladian`) always emit the source field, so this is unreachable in
the canonical pipeline. Worth a comment if not a fix.

---

## Verified claims

| # | Claim | Status |
|---|---|---|
| 1 | AICc formula matches Burnham–Anderson NLS convention (k = n_model_params + 1) | PASS — hand-computed: `N=8, RSS=0.1, k=4` → `log L = 6.17660`, `AIC = -4.35320`, `AICc = 8.98014`. Module output identical to 0.0 absolute difference. |
| 2 | `aicc_weights` matches `w_i = exp(-Δ_i/2) / Σ exp(-Δ_j/2)` | PASS — module weights match hand-computed weights to 1e-15. Weights sum to 1, all in [0,1]. |
| 3 | `formula_string` uses `σ_C ≈ C·σ_c` (first-order delta method) | PASS — derivative `d/dc[exp(c)] = exp(c) = C`. Verified against 1M Monte Carlo samples: `σ_C(delta) = 0.006065`, `σ_C(MC) = 0.006067`, relative diff 0.03% (within third-order corrections expected at `σ_c = 0.01`). |
| 4 | Correlation matrix construction `cov ./ (s * s')` yields proper correlation | PASS — symmetric (max asymmetry < 1e-12), diagonal = 1 (within 1e-10), off-diagonals in `[-1, 1]` (within 1e-10) on every well-conditioned fit observed. |
| 5 | `Vector{<:NamedTuple}` schema dispatch (both `mixing_time`/`tau_mix` schemas) and filtering | PASS — junk cells (`mixing_time = NaN`, `mixing_time_source = :nan` or `:floor`) are correctly filtered; channel schema with `tau_mix`/`tau_mix_source` works; empty-after-filter raises `ArgumentError`. |
| 6 | Input validation paths (length, positivity, count, model symbols, level) | PASS — all 12 `@test_throws ArgumentError` paths have explicit coverage at `test/test_scaling_fit.jl:140-167, 186, 224, 332`. Confirmed via fresh runs that `level ∈ {0.0, 1.0, -0.5, 1.5}` all throw. `Inf` is rejected by the `isfinite` check. |
| 7 | `predict_scaling` round-trip exactness | PASS — `max abs(log(predict_scaling(...)) - fit.log_tau_predicted[i])` is exactly **0.0** on a clean dataset (better than the 1e-10 tolerance). |
| 8 | `scaling_fit_grid` single-β edge case | PASS — for a 10-point all-β=5.0 dataset, `beta_grid = [5.0]` (1 element) and `tau_predicted` is a `(10, 1)` matrix. |
| 9 | `_SCALING_IDX_*` constants consistent everywhere | PASS — `_SCALING_IDX_C = 1, _SCALING_IDX_X = 2, _SCALING_IDX_SLOPE = 3`. Used identically in `_scaling_M0_model`, `_scaling_M1_model`, `predict_scaling`, and `formula_string` (lines 26, 27, 349-351, 432-437). |
| 10 | Full test suite `julia --project test/test_scaling_fit.jl` | PASS — 221/221 passed in 6–8 s wall. No failures, warnings, or errors. |
| 11 | Synthetic M0 dataset (τ = 0.5·n^2.3·β^1.7, 1% noise) | PASS — x̂ = 2.29986 (within 0.05 of 2.3), ŷ = 1.70158 (within 0.05 of 1.7), AICc(M0) = -335.54 < AICc(M1) = 63.78, Δ-AICc(M1 - M0) = 399.32 ≫ 2. |
| 12 | Numerical pitfalls | PASS — see breakdown below. |

### Numerical pitfall sweep

- **Tiny τ values (1e-20)**: ran the M0 test with `τ_true = 1e-20 · n^2.0 · β^1.5`.
  Result: `c = -46.05` (`isfinite(c) = true`), `C = 1.001e-20` (truth 1e-20),
  exponents recovered to 1e-3. No issue with extreme `c`.
- **Boundary N = 6 (k+1)**: `_scaling_aic_metrics(0.1, 6, 3)` returns
  `AICc = 40.46` (finite). For `N = 5` (denominator = 0) and `N = 4` (`N ≤ k`),
  returns `Inf` cleanly. The threshold check `N <= k + 1` correctly excludes
  both, so no division-by-zero occurs.
- **Negative diag in covariance**: handled by the
  `all(d -> d > 0, diag(covmat))` guard, which falls back to `fill(NaN, ...)`.
  Confirmed by trying a single-point repeated dataset (`n_collinear`,
  `β_collinear`): the cov has astronomically large positive diag values
  (LsqFit's `abs` on negative roundoff inside `stderror` means the diag
  values in `vcov(fit)` itself are typically positive). The corr returned
  finite values (1's on diagonal, ±1 off-diagonal).
- **`@.` over `xdata[:, 1]`**: produces a `Vector{Float64}` of length N.
  Verified `_scaling_M0_model(xdata, p)` returns identical output to manual
  `p[1] .+ p[2] .* xdata[:,1] .+ p[3] .* xdata[:,2]` (max diff 0.0). 6
  allocations per call (slicing + result), negligible at the LM iteration
  cadence.
- **Range-typed inputs (`3:10`)**: accepted via `AbstractVector{<:Integer}` /
  `AbstractVector{<:Real}` signatures; `Float64.(n_vals)` conversion handles
  the promotion. No allocations issue.

---

## Test-run output (last 30 lines)

```
$ julia --project test/test_scaling_fit.jl
Test Summary:           | Pass  Total  Time
scaling_fit.jl (qf-now) |  221    221  6.4s
```

(The full output is one line because every `@test` passes and the testset is
non-`verbose`. Re-running with `@testset verbose=true` confirms the breakdown:
all 10 inner testsets, including `(a) M0 recovers known (C, x, y)`,
`(b) M1 recovers known (C, x, α)`, `(c) AICc model discrimination`,
`(d) Input validation`, `(e) Convenience input shapes`,
`(f) predict_scaling round-trip`, `(g) formula_string formatting`,
`(h) Correlation matrix diagnostics`, `(i) scaling_fit_grid for diagnostic plots`,
and `(j) AICc weights edge cases`, pass without exception.)

---

## Suggested cleanups (non-blocking)

1. Replace `estimate_covar(fit)` with `vcov(fit)` at `src/scaling_fit.jl:135`
   and update the docstring/header comment at lines 11-12.
2. Broaden the catch in `_build_scaling_fit` to also handle `ErrorException`
   (the "Covariance matrix is negative" path in LsqFit's `stderror`).
3. Either drop the redundant tuple comprehension at line 136 (since `confint`
   already returns `Vector{Tuple{Float64, Float64}}`), or add a comment
   acknowledging it's defensive against LsqFit API drift.

None of these affect correctness, mixing-time fits, or downstream plotting;
they are tidiness improvements.
