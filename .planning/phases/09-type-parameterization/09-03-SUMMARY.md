---
phase: 09-type-parameterization
plan: 03
subsystem: simulation-functions
tags: [type-parameterization, AbstractFloat, generics, simulation-pipeline, function-signatures, cross-struct-validation]

# Dependency graph
requires:
  - phase: 09-type-parameterization
    plan: 01
    provides: HamHam{T} and TrottTrott{T} parameterized structs
  - phase: 09-type-parameterization
    plan: 02
    provides: Config{D,T}, LindbladianWorkspace{T}, KrausScratch{<:Complex}, NUFFTPrefactors{T}
provides:
  - All simulation functions (run_lindbladian, run_thermalization, run_trajectories) accept generic T
  - Cross-struct T mismatch detection between HamHam{T} and Config{D,T}
  - Generic function signatures across jump_workers, qi_tools, coherent, ofts, domain files
  - Complete type parameterization pipeline: HamHam{T} -> Config{D,T} -> computation -> results{D,T}
affects: [downstream-Float32-support, future-precision-benchmarks]

# Tech tracking
tech-stack:
  added: []
  patterns: [generic-dispatch-via-<:Complex, eltype-inference-for-allocations, cross-struct-T-validation]

key-files:
  created: []
  modified:
    - src/jump_workers.jl
    - src/qi_tools.jl
    - src/coherent.jl
    - src/ofts.jl
    - src/energy_domain.jl
    - src/bohr_domain.jl
    - src/furnace.jl
    - src/trajectories.jl
    - src/time_domain.jl

key-decisions:
  - "Function signatures use <:Complex for dispatch, eltype() for internal allocations"
  - "Cross-struct T mismatch check added at run_lindbladian and run_thermalization entry points"
  - "TrajectoryFramework step parameters (delta, delta_eff, alpha) stay Float64 since they control numerical stepping"
  - "Domain helper functions (create_alpha, create_f, etc.) widened from Float64 to Real for generic acceptance"

patterns-established:
  - "Generic dispatch pattern: function f(x::Matrix{<:Complex}) rather than f(x::Matrix{ComplexF64})"
  - "Allocation inference pattern: CT = eltype(existing_matrix); zeros(CT, dim, dim)"
  - "Entry-point validation: top-level simulation functions check type parameter consistency"

# Metrics
duration: 9min
completed: 2026-02-15
---

# Phase 9 Plan 03: Simulation Function Generic Type Propagation Summary

**All simulation function signatures generalized from hardcoded ComplexF64/Float64 to generic T with cross-struct type mismatch validation**

## Performance

- **Duration:** 9 min
- **Started:** 2026-02-15T12:06:44Z
- **Completed:** 2026-02-15T12:15:49Z
- **Tasks:** 2
- **Files modified:** 9

## Accomplishments
- All function signatures across jump_workers, qi_tools, coherent, ofts, and domain files generalized from hardcoded ComplexF64/Float64 to generic <:Complex/Real types
- Simulation pipeline (run_lindbladian, run_thermalization, run_trajectories) infers element types from inputs and allocates matching buffers
- Cross-struct T mismatch detection added: HamHam{Float32} + Config{D,Float64} produces clear error message
- Full backward compatibility: all 224 existing tests pass unchanged with default Float64 types
- Phase 9 type parameterization complete: HamHam{T}, Config{D,T}, LindbladianWorkspace{T}, and all functions work generically

## Task Commits

Each task was committed atomically:

1. **Task 1: Update jump_workers, qi_tools, coherent, ofts, and domain files for generic T** - `45dc633` (feat)
2. **Task 2: Update furnace, trajectories, and verify full pipeline** - `8b87dc9` (feat)

**Plan metadata:** [pending] (docs: complete plan)

## Files Created/Modified
- `src/jump_workers.jl` - All jump_contribution!, apply_coherent_unitary!, finalize_kraus_step! accept <:Complex instead of ComplexF64
- `src/qi_tools.jl` - vectorize_liouv_diss_and_add!, vectorize_liouvillian_coherent! accept <:Complex; is_density_matrix generalized
- `src/coherent.jl` - compute_b_minus/plus/metro/smooth accept Real; get_truncated_indices/compute_truncated_func generalized
- `src/ofts.jl` - oft!, time_oft!, trotter_oft! accept <:Complex matrices and Real scalars
- `src/energy_domain.jl` - create_energy_labels and truncate_energy_labels accept generic Real/Integer
- `src/bohr_domain.jl` - coherent_bohr infers CT from HamHam{T}; create_f/alpha/alpha_gns/alpha_gauss accept Real
- `src/furnace.jl` - run_lindbladian/run_thermalization infer CT from HamHam{T}; cross-struct T mismatch check added
- `src/trajectories.jl` - build_trajectoryframework, precompute_R, run_trajectories, step_along_trajectory! all accept generic types
- `src/time_domain.jl` - truncate_time_labels_for_oft accepts AbstractVector{<:Real}

## Decisions Made
- **Generic dispatch via <:Complex**: Function signatures use `Matrix{<:Complex}` rather than `Matrix{Complex{T}}` for dispatch flexibility. Internal allocations use `eltype()` to infer the concrete type. This follows Julia best practices for generic programming.
- **Cross-struct T validation at entry points**: Type mismatch check added only at `run_lindbladian` and `run_thermalization` (the user-facing entry points), not at internal functions. This avoids redundant checks in hot paths.
- **TrajectoryFramework step parameters stay Float64**: The `delta`, `delta_eff`, and `alpha` fields in TrajectoryFramework remain Float64 since they control numerical stepping precision and are always double precision regardless of the simulation's element type.
- **Domain helper functions widened to Real**: Functions like `create_alpha`, `create_f`, `compute_b_minus` etc. changed from `Float64` to `Real` annotations so they accept any AbstractFloat subtype.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Added time_domain.jl generalization**
- **Found during:** Task 2 (furnace_utensils/trajectories pipeline)
- **Issue:** `truncate_time_labels_for_oft` in time_domain.jl had hardcoded `Vector{Float64}` and `Float64` parameters, blocking the generic pipeline
- **Fix:** Changed to `AbstractVector{<:Real}` and `Real` annotations
- **Files modified:** src/time_domain.jl
- **Verification:** All 224 tests pass
- **Committed in:** 8b87dc9 (Task 2 commit)

---

**Total deviations:** 1 auto-fixed (1 blocking)
**Impact on plan:** Auto-fix necessary to complete the generic pipeline. time_domain.jl was not listed in the plan's files but contained a hardcoded Float64 function on the critical path. No scope creep.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Phase 9 type parameterization is complete: the entire pipeline from HamHam{T} through computation to results is type-consistent
- Float32 end-to-end paths are now structurally supported (all types parameterized, all functions generic)
- Actual Float32 testing and optimization would be a future phase concern
- All existing Float64 code continues to work without modification
- The `precision=Float32` kwarg on HamHam and Config constructors is ready for use

## Self-Check: PASSED

- All 9 modified files exist on disk
- SUMMARY.md created at expected path
- Commit 45dc633 (Task 1) found in git log
- Commit 8b87dc9 (Task 2) found in git log
- All 224 tests pass

---
*Phase: 09-type-parameterization*
*Completed: 2026-02-15*
