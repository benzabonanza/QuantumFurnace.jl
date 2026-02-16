---
phase: 19-logic-simplification
verified: 2026-02-16T14:49:25Z
status: passed
score: 17/17 must-haves verified
re_verification: false
---

# Phase 19: Logic Simplification Verification Report

**Phase Goal:** Simplify overly complex logic accumulated during v1.2 development: flatten trajectory call chain from 5 to 3 levels, eliminate redundant jump basis transforms, and simplify result struct hierarchy

**Verified:** 2026-02-16T14:49:25Z
**Status:** passed
**Re-verification:** No - initial verification

## Goal Achievement

### Observable Truths

| #   | Truth | Status | Evidence |
|-----|-------|--------|----------|
| 1   | JumpOp construction uses domain-aware basis (trotter.eigvecs for TrotterDomain, hamiltonian.eigvecs otherwise) | ✓ VERIFIED | test/test_helpers.jl:126 selects basis via trotter kwarg, experiments/run_sweep.jl:69 uses trotter.eigvecs |
| 2   | No internal code path calls transform_jumps_to_basis for TrotterDomain basis conversion | ✓ VERIFIED | grep shows zero hits in trajectories.jl, furnace.jl, coherent.jl |
| 3   | All 539 existing tests pass with updated test call sites | ✓ VERIFIED | Full test suite passes: 539/539 |
| 4   | Test helpers produce both hamiltonian-basis and trotter-basis JumpOps via optional trotter parameter | ✓ VERIFIED | TEST_TROTTER_JUMPS and SMALL_TROTTER_JUMPS constants exist at test/test_helpers.jl:217,280 |
| 5   | TrajectoryFramework is built once and reused across batches in adaptive/convergence runners | ✓ VERIFIED | src/convergence.jl:169 builds fw once, line 196 reuses in loop |
| 6   | The _evolve_along_trajectory! wrapper is eliminated -- step loop is inlined into chunk runners | ✓ VERIFIED | Function deleted (only found in trajectories.jl comments), step loop at src/trajectories.jl:376-378 |
| 7   | The call chain from public API to physics step is 3 levels deep, not 5 | ✓ VERIFIED | run_trajectories_adaptive -> _run_batch_no_obs! -> _run_chunk_no_obs! -> step_along_trajectory! (3 levels) |
| 8   | All 539 existing tests pass unchanged | ✓ VERIFIED | Full test suite passes: 539/539 |
| 9   | ConvergenceData is embedded inside TrajectoryResult (not returned as a separate tuple element) | ✓ VERIFIED | src/trajectories.jl:29 convergence field exists |
| 10  | HotSpectralResults is renamed to LindbladianResult with reduced fields (no hamiltonian/config baggage) | ✓ VERIFIED | src/structs.jl:274 LindbladianResult with 4 fields, zero HotSpectralResults in src/ |
| 11  | HotAlgorithmResults is renamed to DMSimulationResult with reduced fields (no hamiltonian/config baggage) | ✓ VERIFIED | src/structs.jl:254 DMSimulationResult with 3 fields, zero HotAlgorithmResults in src/ |
| 12  | Convergence/adaptive runners return TrajectoryResult (single value, not tuple) | ✓ VERIFIED | src/convergence.jl:231,385 return TrajectoryResult, zero tuple destructuring in experiments/ |
| 13  | All existing tests pass (updated to match new struct shapes) | ✓ VERIFIED | Full test suite passes: 539/539 |
| 14  | BSON serialization round-trips correctly with new struct layouts | ✓ VERIFIED | src/results.jl:129-133 backward-compatible loading, test_results.jl passes |

**Score:** 14/14 truths verified (17/17 total must-haves when counting artifacts and key links)

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| test/test_helpers.jl | make_test_system and make_small_test_system with optional trotter kwarg; SMALL_TROTTER_JUMPS / TEST_TROTTER_JUMPS constants | ✓ VERIFIED | Lines 112-136, 217, 280 - functions accept trotter kwarg, constants defined |
| experiments/run_sweep.jl | build_trotter_system uses trotter.eigvecs for JumpOp construction | ✓ VERIFIED | Line 69 uses trotter.eigvecs for basis selection |
| src/trajectories.jl (plan 01) | build_trajectoryframework without transform_jumps_to_basis for TrotterDomain | ✓ VERIFIED | Lines 91-93 - comments confirm jumps arrive in correct basis, no transform call |
| src/furnace.jl | construct_lindbladian and run_thermalization without transform_jumps_to_basis | ✓ VERIFIED | Zero hits for transform_jumps_to_basis in file |
| src/coherent.jl | Coherent term functions using jump.in_eigenbasis directly for TrotterDomain | ✓ VERIFIED | Zero hits for transform_jumps_to_basis in file |
| src/trajectories.jl (plan 02) | Flattened trajectory execution with _run_batch_no_obs! accepting pre-built framework, step loop inlined into chunk runners | ✓ VERIFIED | Lines 376-378 inline step loop, _run_batch_no_obs! at line 441 |
| src/convergence.jl | Convergence runners that build framework once and pass to batch execution | ✓ VERIFIED | Line 169 builds fw once, line 196 calls _run_batch_no_obs! with pre-built fw |
| src/structs.jl | LindbladianResult and DMSimulationResult replacing HotSpectralResults/HotAlgorithmResults | ✓ VERIFIED | Lines 254, 274 - new structs defined with slim fields |
| src/trajectories.jl (plan 03) | TrajectoryResult with convergence::Union{Nothing, ConvergenceData} field | ✓ VERIFIED | Line 29 - convergence field exists |
| src/results.jl | Updated serialization for new struct layouts, backward-compatible loading | ✓ VERIFIED | Lines 129-133 - backward-compatible convergence field loading |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|----|--------|---------|
| test/test_helpers.jl | JumpOp constructor | optional trotter kwarg selects basis_unitary | ✓ WIRED | Line 126: basis_unitary = trotter !== nothing ? trotter.eigvecs : hamiltonian.eigvecs |
| test/test_gns_trajectory.jl | test/test_helpers.jl | TrotterDomain tests use SMALL_TROTTER_JUMPS instead of SMALL_JUMPS | ✓ WIRED | Lines 18, 52, 68, 81 use SMALL_TROTTER_JUMPS (verified by passing tests) |
| src/trajectories.jl | build_trajectoryframework | jumps already in correct basis, no transform needed | ✓ WIRED | Line 93: direct conversion without transform |
| src/convergence.jl | src/trajectories.jl | Convergence runners call a lower-level batch function with pre-built framework | ✓ WIRED | Line 196: _run_batch_no_obs!(fw, psi0, batch_size, batch_seed, total_time) |
| src/trajectories.jl | step_along_trajectory! | Step loop directly in chunk runners, no _evolve_along_trajectory! wrapper | ✓ WIRED | Line 377: direct call to step_along_trajectory! in loop |
| src/convergence.jl | src/trajectories.jl | TrajectoryResult constructor with convergence field | ✓ WIRED | Lines 231, 385: TrajectoryResult(..., conv_data) |
| src/results.jl | src/trajectories.jl | _trajectory_to_dict serializes convergence field | ✓ WIRED | Lines 129-133: convergence field handling in deserialization |

### Requirements Coverage

No phase-specific requirements mapped in REQUIREMENTS.md - this is a refactoring/simplification phase.

### Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| - | - | - | - | - |

**No anti-patterns detected.** All code is substantive, properly wired, and free of stubs/placeholders.

### Human Verification Required

No human verification needed. All changes are internal refactoring with automated test coverage:
- Test suite comprehensively validates correctness (539 tests)
- Determinism tests verify seed behavior unchanged
- Convergence tests verify adaptive/convergence runners work correctly
- Serialization tests verify backward compatibility

### Summary

**All phase goals achieved:**

1. **Jump basis transforms eliminated** - All 6 internal `transform_jumps_to_basis` calls removed from src/trajectories.jl, src/furnace.jl, src/coherent.jl. JumpOps constructed with correct basis at source (trotter.eigvecs for TrotterDomain, hamiltonian.eigvecs otherwise).

2. **Call chain flattened from 5 to 3 levels** - Architecture now: (1) public API (run_trajectories/run_trajectories_adaptive/run_trajectories_convergence) → (2) batch execution (_run_batch_no_obs! / _run_chunk_no_obs!) → (3) step_along_trajectory!. The _evolve_along_trajectory! wrapper eliminated.

3. **Result struct hierarchy simplified** - HotSpectralResults → LindbladianResult (4 fields, no config/hamiltonian baggage). HotAlgorithmResults → DMSimulationResult (3 fields, no config/hamiltonian baggage). ConvergenceData embedded in TrajectoryResult.convergence instead of tuple return.

**Impact:**
- Zero functional regressions (539/539 tests pass)
- Framework building eliminated in convergence/adaptive batch loops (performance improvement)
- Cleaner API (single return value instead of tuples)
- Simpler struct definitions (only output data, not input references)

---

_Verified: 2026-02-16T14:49:25Z_
_Verifier: Claude (gsd-verifier)_
