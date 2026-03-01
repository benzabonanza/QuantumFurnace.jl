---
phase: 39-per-jump-precomputation
plan: 02
subsystem: thermalize
tags: [precomputation, cptp-channel, hot-loop, eigendecomposition-elimination, bohr-domain]

# Dependency graph
requires:
  - phase: 39-per-jump-precomputation
    plan: 01
    provides: "_precompute_per_jump_channels, _accumulate_rho_jump!, _apply_precomputed_channel!"
provides:
  - "run_thermalize with precomputed channel application (no eigen() in hot loop)"
  - "CPTP completeness verified for BohrDomain via _precompute_per_jump_channels"
  - "DM precomputation path matches trajectory workspace K0s/U_residuals within 1e-15"
affects: [run_thermalize, phase-40-save-every, phase-41-blas-threading]

# Tech tracking
tech-stack:
  added: []
  patterns: ["precomputed CPTP channel application in DM hot loop"]

key-files:
  created: []
  modified:
    - src/furnace.jl
    - src/furnace_utensils.jl
    - src/bohr_domain.jl
    - src/coherent.jl
    - test/test_cptp.jl

key-decisions:
  - "jump_weight_scaling precomputed before hot loop: gamma_norm_factor / p_jump when rescale_by_inv_prob=true"
  - "B_bohr/B_time/B_trotter signatures relaxed from Vector{JumpOp} to AbstractVector{<:JumpOp} for type invariance correctness"

patterns-established:
  - "run_thermalize hot loop: coherent -> accumulate_rho_jump -> apply_precomputed_channel (3-step pattern)"

# Metrics
duration: 9min
completed: 2026-03-01
---

# Phase 39 Plan 02: Per-Jump Precomputation Integration Summary

**Integrated per-jump precomputed CPTP channels into run_thermalize, eliminating eigendecomposition from the hot loop across all 4 domains**

## Performance

- **Duration:** 9 min
- **Started:** 2026-03-01T09:34:59Z
- **Completed:** 2026-03-01T09:44:51Z
- **Tasks:** 2
- **Files modified:** 5

## Accomplishments
- Refactored `run_thermalize` to call `_precompute_per_jump_channels` before the hot loop, then use `_apply_coherent_unitary!` + `_accumulate_rho_jump!` + `_apply_precomputed_channel!` inside the loop
- Verified numerical equivalence across all 4 domains (Energy, Time, Trotter, Bohr) -- trace distances decrease monotonically from maximally mixed state
- Extended CPTP completeness tests to cover BohrDomain via `_precompute_per_jump_channels`
- Added DM precomputation vs trajectory workspace cross-validation (K0s/U_residuals match within 1e-15)
- Full test suite passes: 1141 tests green

## Task Commits

Each task was committed atomically:

1. **Task 1: Refactor run_thermalize to use precomputed channels** - `c720bfe` (feat)
2. **Task 2: CPTP test extension and full test suite validation** - `84134fc` (test)

## Files Created/Modified
- `src/furnace.jl` - Refactored `run_thermalize` hot loop to use precomputed channels instead of `_jump_contribution!`
- `src/furnace_utensils.jl` - Fixed `_precompute_per_jump_channels` type handling for TrottTrott and JumpOp invariance
- `src/bohr_domain.jl` - Fixed `B_bohr` signature from `Vector{JumpOp}` to `AbstractVector{<:JumpOp}`
- `src/coherent.jl` - Fixed `B_time`/`B_trotter` signatures from `Vector{JumpOp}` to `AbstractVector{<:JumpOp}`
- `test/test_cptp.jl` - Added BohrDomain CPTP completeness test and DM-vs-trajectory cross-validation

## Decisions Made
- `jump_weight_scaling` is precomputed once before the hot loop as `gamma_norm_factor / p_jump` (when `rescale_by_inv_prob=true`), matching the existing `_jump_contribution!` convention
- DM-vs-trajectory cross-validation uses `atol=1e-15` (near machine precision) rather than bit-exact `atol=0` to accommodate potential copy/construction order differences

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed B_bohr/B_time/B_trotter type invariance for parametric JumpOp**
- **Found during:** Task 1 (run_thermalize verification)
- **Issue:** `B_bohr(hamiltonian, [jump], config)` fails because `[jump]` creates `Vector{JumpOp{Matrix{ComplexF64}}}` which is not a subtype of `Vector{JumpOp}` due to Julia's type invariance. Same issue in `B_time` and `B_trotter`. Pre-existing bug that was never triggered because `run_thermalize` had no tests.
- **Fix:** Changed signatures from `Vector{JumpOp}` to `AbstractVector{<:JumpOp}` in `B_bohr`, `B_time`, `B_trotter`
- **Files modified:** `src/bohr_domain.jl`, `src/coherent.jl`
- **Verification:** All 4 domains now execute `run_thermalize` without error
- **Committed in:** `c720bfe` (Task 1 commit)

**2. [Rule 1 - Bug] Fixed _precompute_per_jump_channels TrottTrott field access**
- **Found during:** Task 1 (TrotterDomain verification)
- **Issue:** `_precompute_per_jump_channels` accessed `ham_or_trott.eigvals` and `ham_or_trott.data` for TrottTrott, but TrottTrott has `eigvals_t0` (not `eigvals`) and no `data` field (uses `eigvecs` for dimension). Bug from Plan 39-01.
- **Fix:** Changed to `ham_or_trott.eigvals_t0` for type extraction and `ham_or_trott.eigvecs` for dimension; also fixed `Complex{Complex{T}}` issue by using `eltype(eigvecs)` directly
- **Files modified:** `src/furnace_utensils.jl`
- **Verification:** TrotterDomain `run_thermalize` completes successfully
- **Committed in:** `c720bfe` (Task 1 commit)

**3. [Rule 1 - Bug] Fixed _precompute_per_jump_channels JumpOp type invariance**
- **Found during:** Task 1 (code review during fix #1)
- **Issue:** Same `Vector{JumpOp}` invariance issue as fix #1
- **Fix:** Changed to `AbstractVector{<:JumpOp}`
- **Files modified:** `src/furnace_utensils.jl`
- **Committed in:** `c720bfe` (Task 1 commit)

---

**Total deviations:** 3 auto-fixed (3 bugs: 2 pre-existing, 1 from Plan 39-01)
**Impact on plan:** All auto-fixes necessary for correctness. No scope creep.

## Issues Encountered

None beyond the auto-fixed bugs documented above.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- Phase 39 (Per-Jump Precomputation) is complete
- `run_thermalize` no longer calls `eigen()` or `_build_cptp_channel` in the per-step hot loop
- All 4 domains produce correct, monotonically-converging trace distances
- Full test suite (1141 tests) passes
- Ready for Phase 40 (save_every optimization)

## Self-Check: PASSED

All files and commits verified.

---
*Phase: 39-per-jump-precomputation*
*Completed: 2026-03-01*
