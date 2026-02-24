---
phase: 29-eigensolver-integration
plan: 01
subsystem: krylov-eigensolver
tags: [KrylovKit, Arnoldi, spectral-gap, eigsolve, matrix-free, CPTP-channel]

# Dependency graph
requires:
  - phase: 27-core-matvec-infrastructure
    provides: "KrylovWorkspace, apply_lindbladian!, apply_adjoint_lindbladian! for all 4 domains"
  - phase: 28-domain-matvec-validation
    provides: "Validated matvec correctness for Energy, Time, Trotter, Bohr domains"
provides:
  - "KrylovGapResult struct for matrix-free spectral gap results"
  - "krylov_spectral_gap() with config-type dispatch (Lindbladian :LR and channel :LM)"
  - "apply_delta_channel! for CPTP channel E(rho) = rho + delta*L(rho)"
  - "_eigsolve_with_retry convergence retry with 1.5x krylovdim increase"
  - "_check_krylov_memory pre-flight advisory memory guard"
  - "_thermalize_to_liouv_config for ThermalizeConfig -> LiouvConfig conversion (2 methods)"
affects: [29-02-testing, 30-cross-validation, 31-production-runs]

# Tech tracking
tech-stack:
  added: [KrylovKit 0.8/0.9/0.10]
  patterns: [config-type dispatch for eigsolve path, KrylovKit function-based linear map, convergence retry wrapper]

key-files:
  created: [src/krylov_eigsolve.jl]
  modified: [Project.toml, src/QuantumFurnace.jl]

key-decisions:
  - "Arnoldi algorithm (not Lanczos) for non-Hermitian Lindbladian"
  - "copy(vec(ws.rho_out)) in KrylovKit closure to prevent aliasing"
  - "Exact linear formula lambda_L = (mu-1)/delta for channel eigenvalue conversion"
  - "50% krylovdim increase per retry (30->45->68->102), 3 retries max"
  - "Maximally mixed initial vector I/dim for stable overlap with steady state"

patterns-established:
  - "Config-type dispatch: AbstractLiouvConfig -> :LR targeting, AbstractThermalizeConfig -> :LM targeting"
  - "KrylovKit closure pattern: reshape input view, copy output to avoid aliasing"
  - "ThermalizeConfig -> LiouvConfig conversion for reusing Lindbladian matvec in channel path"

# Metrics
duration: 3min
completed: 2026-02-24
---

# Phase 29 Plan 01: Eigensolver Integration Summary

**Matrix-free spectral gap API wrapping KrylovKit Arnoldi with config-type dispatch for Lindbladian (:LR) and CPTP channel (:LM) paths, convergence retry, and memory guard**

## Performance

- **Duration:** 3 min
- **Started:** 2026-02-24T12:40:34Z
- **Completed:** 2026-02-24T12:43:52Z
- **Tasks:** 2
- **Files modified:** 3

## Accomplishments
- KrylovGapResult struct with 10 fields covering eigenvalues, spectral gap, fixed point, gap mode, convergence metadata, and channel-specific data
- Two-method krylov_spectral_gap: Lindbladian path (AbstractLiouvConfig, :LR) and channel path (AbstractThermalizeConfig, :LM with delta conversion)
- apply_delta_channel! implementing E(rho) = rho + delta*L(rho) reusing existing Lindbladian matvec
- Convergence retry (_eigsolve_with_retry) with 50% krylovdim increase per attempt up to 3 retries
- Pre-flight memory guard warning when estimated bytes exceed 80% of Sys.free_memory()
- ThermalizeConfig -> LiouvConfig conversion for both standard and GNS config types

## Task Commits

Each task was committed atomically:

1. **Task 1: Add KrylovKit dependency and create krylov_eigsolve.jl** - `a3b8adf` (feat)
2. **Task 2: Integrate krylov_eigsolve.jl into QuantumFurnace module** - `c11d160` (feat)

**Plan metadata:** TBD (docs: complete plan)

## Files Created/Modified
- `src/krylov_eigsolve.jl` - KrylovGapResult, krylov_spectral_gap (2 dispatch methods), apply_delta_channel!, _eigsolve_with_retry, _check_krylov_memory, _thermalize_to_liouv_config (2 methods)
- `Project.toml` - KrylovKit added to [deps] and [compat] (0.8, 0.9, 0.10)
- `src/QuantumFurnace.jl` - using KrylovKit, include krylov_eigsolve.jl, export KrylovGapResult/krylov_spectral_gap/apply_delta_channel!

## Decisions Made
- **Arnoldi over Lanczos:** Lindbladian is non-Hermitian, Arnoldi is the only valid Krylov algorithm. KrylovKit's Krylov-Schur variant with thick restarts provides numerical stability.
- **copy(vec(ws.rho_out)) aliasing prevention:** vec() returns a view into ws.rho_out which gets overwritten each call. Without copy, all Krylov basis vectors would alias the same memory, producing garbage.
- **Exact linear conversion lambda_L = (mu-1)/delta:** E = I + delta*L is linear (not exponential), so the eigenvalue relationship is exact. No log approximation needed.
- **50% krylovdim increase per retry:** Balanced between too-small (insufficient improvement) and doubling (memory waste). Sequence 30->45->68->102 gives 3.4x range.
- **Maximally mixed initial vector:** I/dim has guaranteed overlap with the steady state (trace-preserving property), unlike random vectors that might have poor overlap.
- **Removed redundant using KrylovKit from included file:** Module-level using already brings symbols into scope for all included files, consistent with codebase convention.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Removed redundant `using KrylovKit: eigsolve, Arnoldi` from krylov_eigsolve.jl**
- **Found during:** Task 2 (Module integration)
- **Issue:** The included file had its own `using KrylovKit` which is redundant since the module-level `using KrylovKit` already brings all symbols into scope. Other krylov files (krylov_matvec.jl, krylov_workspace.jl) do not have their own using statements.
- **Fix:** Removed the `using KrylovKit: eigsolve, Arnoldi` line from krylov_eigsolve.jl
- **Files modified:** src/krylov_eigsolve.jl
- **Verification:** Consistent with codebase convention
- **Committed in:** c11d160 (Task 2 commit)

---

**Total deviations:** 1 auto-fixed (1 blocking/consistency)
**Impact on plan:** Trivial cleanup for codebase consistency. No scope creep.

## Issues Encountered
- Julia is not installed in the sandbox environment, so the verification commands (`julia -e 'using QuantumFurnace; ...'`) could not be executed. The code was verified structurally: all required functions exist, type signatures are correct, include order respects dependencies, and exports are properly listed. Functional verification will occur when Julia is available or during Phase 29 Plan 02 (testing).

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- krylov_eigsolve.jl is ready for testing in Plan 02 (cross-validation against dense eigen)
- All four domain paths (Energy, Time, Trotter, Bohr) are supported through the existing apply_lindbladian! dispatch
- Channel path (ThermalizeConfig) is ready for validation against existing thermalization results
- KrylovKit needs to be installed (`Pkg.instantiate()`) before first use

## Self-Check: PASSED

- FOUND: src/krylov_eigsolve.jl
- FOUND: Project.toml (KrylovKit in deps and compat)
- FOUND: src/QuantumFurnace.jl (using, include, exports)
- FOUND: a3b8adf (Task 1 commit)
- FOUND: c11d160 (Task 2 commit)

---
*Phase: 29-eigensolver-integration*
*Completed: 2026-02-24*
