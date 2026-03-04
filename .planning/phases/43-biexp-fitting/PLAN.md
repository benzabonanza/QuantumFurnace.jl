# Phase 43: Bi-Exponential Fitting for Improved Mixing Time Extrapolation

**Status:** Planned, not yet executed
**Date:** 2026-03-03

## Context

The current single-exponential fit `d(t) = A*exp(-gap*t) + C` overestimates the offset C by ~47% (9.9e-5 vs true 6.8e-5). When the target epsilon is close to C, the extrapolation formula `t = -ln((ε-C)/A)/gap` is hypersensitive to C: a 47% error in C causes a 32x error in `ε-C`, which translates to ~26% error in predicted mixing time. A bi-exponential model captures the multi-timescale Liouvillian decay, yielding a more accurate C and better extrapolation.

## Design Decision

User chose `:biexp` as an explicit option only (no `:auto` mode with AICc). Keep it simple.

## Changes

### 1. `src/fitting.jl` — Add bi-exponential fitting (append after existing code)

**New struct `BiexpFitResult`:**
- Fields: `gap` (slow), `gap_fast`, `amplitude` (slow), `amplitude_fast`, `offset`, `gap_ci`, `gap_se`, `r_squared`, `converged`, `residuals`, `times_used`, `values_used`
- `gap` is the slow mode (spectral gap estimate, used for extrapolation)

**New model:** `_biexp_decay_model(t, p) = A1*exp(-g1*t) + A2*exp(-g2*t) + C` (5 params)

**New function `fit_biexponential_decay(times, values; skip_initial, p0, level)`:**
- Requires >= 8 data points (5 params need margin)
- Initial guess strategy: fit single-exp first, analyze residuals to seed the fast mode
- Bounds: both gaps >= 0
- After fitting, sort so g1 >= g2 (fast >= slow); extract SE/CI for the slow gap using the pre-swap index
- Helper: `_biexp_initial_guess(times, values, single_fit::FitResult)`

### 2. `src/mixing.jl` — Extend `estimate_mixing_time` with `model` keyword

**Modified `MixingTimeEstimate` struct:** Add two fields at end:
- `model_used::Symbol` — `:single` or `:biexp`
- `biexp_fit_result::Union{Nothing, BiexpFitResult}` — `nothing` when `model=:single`

**Modified `estimate_mixing_time` signature:** Add `model::Symbol = :single`
- `:single` — current behavior (backward compatible default)
- `:biexp` — use bi-exponential fit; populate `fitted_gap`/`amplitude`/`offset` from slow mode; store full `BiexpFitResult` in `biexp_fit_result`; construct a synthetic `FitResult` for `fit_result` field (backward compat)

**New helper `_extrapolate_mixing_time_biexp(fit::BiexpFitResult, target_epsilon)`:**
- Solves `A1*exp(-g1*t) + A2*exp(-g2*t) + C = ε` via `Roots.find_zero` with `Bisection()` (Roots.jl already a dependency)
- Bracket: `[0, t_upper]` where `t_upper` from slow-mode estimate with 2x safety margin
- Guards: offset >= target, f(0) <= target, etc.

**Helper `_biexp_to_single_fit_result(bifit)`:** Constructs a `FitResult` from `BiexpFitResult` slow-mode params for the `fit_result` field.

### 3. `src/QuantumFurnace.jl` — Add exports
- `export fit_biexponential_decay, BiexpFitResult` (line 78)

### 4. `test/test_fitting.jl` — Add bi-exponential fitting tests
- **BIEXP-01:** Clean bi-exp data recovery (gap, gap_fast, offset, R² > 0.999)
- **BIEXP-02:** Offset accuracy comparison — bi-exp offset closer to true C than single-exp offset (the key validation)
- **BIEXP-03:** skip_initial works with bi-exp

### 5. `test/test_mixing.jl` — Add bi-exp mixing time tests
- **BIEXP-MIX-01:** Extrapolation accuracy <5% on synthetic bi-exp data (the acceptance criterion)
- **BIEXP-MIX-02:** Backward compat — default `model=:single` produces identical results, `biexp_fit_result === nothing`
- Update MIX-05 field check to include new `model_used` and `biexp_fit_result` fields

### 6. `scripts/mixing_time_extrapolate_verify.jl` — Update demo
- Add a second extrapolation call with `model=:biexp` alongside the existing single-exp call
- Print both predictions and compare errors

## Key Implementation Notes

- `Roots.jl` v2.x is already a dependency (`using Roots` in QuantumFurnace.jl)
- `LsqFit.jl` v0.15.x is already imported
- After fitting bi-exp, sort params so g1 >= g2; track pre-swap indices for SE/CI extraction
- The Phase 42 research labeled bi-exp as anti-feature — this is now overridden by empirical evidence of 26% extrapolation error

## Acceptance Criteria
1. All existing tests pass unchanged
2. Bi-exp extrapolation achieves <5% error on synthetic bi-exp data
3. `model=:single` (default) preserves exact current behavior

## Verification
1. Run `julia --project -e 'using Pkg; Pkg.test()'` — all existing + new tests pass
2. Run `julia --project scripts/mixing_time_extrapolate_verify.jl` — bi-exp shows <5% error vs 26% from single-exp
