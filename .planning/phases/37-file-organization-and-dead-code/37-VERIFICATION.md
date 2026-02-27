---
phase: 37-file-organization-and-dead-code
verified: 2026-02-27T16:05:15Z
status: passed
score: 7/7 must-haves verified
re_verification: false
---

# Phase 37: File Organization and Dead Code Verification Report

**Phase Goal:** Source files are renamed for clarity with PRE/MID/POST logical grouping, dead code is removed, staging code is separated, and the module export list matches the new structure
**Verified:** 2026-02-27T16:05:15Z
**Status:** passed
**Re-verification:** No -- initial verification

## Scope Note

ROADMAP success criterion 1 ("All src/ files are renamed with clear names reflecting their role in the PRE/MID/POST architecture") was **explicitly descoped by the user** during the context-gathering phase. From `37-CONTEXT.md` line 9:

> "File renaming/PRE-MID-POST reorganization of src/ is **excluded** (user will handle manually)."

The phase plans correctly excluded file renaming. Verification is against the actual phase scope: dead code removal, staging separation, Distributed/dead-import removal, diagnostics independence, and export list reorganization.

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | Old entry points run_lindbladian() and run_thermalization() no longer exist in active source | VERIFIED | `grep` for `run_lindbladian\|run_thermalization` in `src/` returns zero matches. furnace.jl contains only construct_lindbladian, run_lindblad, run_thermalize. |
| 2 | Dead structs DMSimulationResult, LindbladianResult, LSIFramework no longer exist in active source | VERIFIED | `grep` for these names in `src/` returns only: (a) STAGING comment in QuantumFurnace.jl:78 referencing LSIFramework, (b) commented-out code in src/staging/log_sobolev.jl:18. No active struct definitions remain. |
| 3 | gap_estimation.jl, fitting.jl, log_sobolev.jl are in src/staging/ | VERIFIED | `ls src/staging/` confirms all three files present. File sizes: fitting.jl=217 lines, gap_estimation.jl=312 lines, log_sobolev.jl=206 lines. Not empty stubs. |
| 4 | Staging tests are in test/staging/ and excluded from runtests.jl | VERIFIED | `ls test/staging/` confirms test_gap_estimation.jl (279 lines), test_fitting.jl (155 lines). runtests.jl contains no include for either file. |
| 5 | Distributed, LsqFit, Optim, SharedArrays removed from module and Project.toml | VERIFIED | `grep` for `using Distributed\|using Optim\|using LsqFit\|using SharedArrays` in QuantumFurnace.jl returns no matches. `grep` for `Distributed\|LsqFit\|Optim\|SharedArrays` in Project.toml returns no matches. No `nprocs` or `@distributed` anywhere in src/. |
| 6 | Diagnostics remains as a separate analysis module | VERIFIED | diagnostics.jl is 573 lines with its own result structs, 6 DIAG functions, and `run_exact_diagnostics()` bundle. Exported in dedicated `# --- Diagnostics ---` section. Not folded into furnace.jl or any simulation path. |
| 7 | Module export list is organized by simulation type with clean groupings | VERIFIED | QuantumFurnace.jl exports organized into 6 sections: `# --- Lindbladian ---`, `# --- Thermalize ---`, `# --- Krylov ---`, `# --- Trajectory ---`, `# --- Diagnostics ---`, `# --- Common ---`. Dormant exports preserved as `# STAGING:` commented block (lines 76-78). |

**Score:** 7/7 truths verified

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `src/staging/gap_estimation.jl` | Staged gap estimation code | VERIFIED | 312 lines, moved from src/ via git mv |
| `src/staging/fitting.jl` | Staged exponential fitting code | VERIFIED | 217 lines, moved from src/ via git mv |
| `src/staging/log_sobolev.jl` | Staged LSI code (commented out) | VERIFIED | 206 lines, moved from src/ via git mv |
| `test/staging/test_gap_estimation.jl` | Staged gap estimation tests | VERIFIED | 279 lines, moved from test/ via git mv |
| `test/staging/test_fitting.jl` | Staged fitting tests | VERIFIED | 155 lines, moved from test/ via git mv |
| `src/QuantumFurnace.jl` | Clean module definition with organized exports | VERIFIED | 104 lines, 6 export sections, STAGING block, no dead includes/imports |
| `test/runtests.jl` | Updated test runner without staged test includes | VERIFIED | No references to test_fitting or test_gap_estimation |
| `simulations/main_krylov.jl` | Krylov spectrum simulation script | VERIFIED | 121 lines, uses run_krylov_spectrum, substantive script with main() function |
| `src/errors.jl` | DELETED (empty placeholder) | VERIFIED | File does not exist |
| `src/kraus.jl` | DELETED (empty placeholder) | VERIFIED | File does not exist |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `src/furnace.jl` | `src/structs.jl` | LindbladResults, ThermalizeResults usage | WIRED | run_lindblad returns LindbladResults{T} (line 110), run_thermalize returns ThermalizeResults{T} (line 216) |
| `src/convergence.jl` | `src/trajectories.jl` | _run_trajectory_convergence, _run_trajectory_adaptive | WIRED | trajectories.jl calls _run_trajectory_convergence (line 1223) and _run_trajectory_adaptive (line 1205); convergence.jl defines both functions (lines 408, 448) |
| `src/QuantumFurnace.jl` | All active src/*.jl | include() statements | WIRED | 22 include() statements for 22 active source files. No includes for deleted/staged files. |
| `test/runtests.jl` | Active test files | include() statements | WIRED | 18 test includes. No includes for staged tests. |
| `simulations/main_krylov.jl` | `src/QuantumFurnace.jl` | using QuantumFurnace + run_krylov_spectrum | WIRED | Lines 20-21 load module; line 97 calls run_krylov_spectrum |

### Requirements Coverage

| Requirement | Status | Notes |
|-------------|--------|-------|
| ORG-01 (Rename src/ files PRE/MID/POST) | DESCOPED | User explicitly excluded from this phase: "user will handle manually" |
| ORG-02 (Move gap/fitting to staging) | SATISFIED | src/staging/ and test/staging/ created with all 5 files |
| ORG-03 (Remove @distributed and using Distributed) | SATISFIED | No @distributed, nprocs, Distributed, or SharedArrays anywhere in active source |
| ORG-07 (4 simulation scripts matching 4 run_* entry points) | SATISFIED | main_liouv.jl, main_thermalize.jl, main_trajectory.jl, main_krylov.jl all exist |
| ORG-08 (Clean module export list by simulation type) | SATISFIED | 6 sections: Lindbladian, Thermalize, Krylov, Trajectory, Diagnostics, Common |
| ORG-09 (Diagnostics as separate analysis module) | SATISFIED | diagnostics.jl is 573-line standalone module with own export section |

### Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| `src/structs.jl` | 25 | "placeholder for Ding et al." comment | Info | Design comment about DLL construction type, not a code stub |

No blocker or warning-level anti-patterns found.

### Human Verification Required

### 1. Module loads without error

**Test:** Run `julia --project=. -e 'using Pkg; Pkg.activate("."); using QuantumFurnace'`
**Expected:** Module loads successfully with no errors
**Why human:** Requires Julia runtime to verify all includes resolve and exports are valid

### 2. Full test suite passes

**Test:** Run `julia --project=. -e 'using Pkg; Pkg.test()'`
**Expected:** All tests pass (SUMMARY reports 4 pre-existing flaky allocation tests in test_krylov_matvec.jl)
**Why human:** Requires Julia runtime execution

### 3. SharedArrays removal is safe

**Test:** Verify no runtime code path needs SharedArrays
**Expected:** After removing the `nprocs() > 1` dead branches, no code needs SharedArrays
**Why human:** The ROADMAP said "SharedArrays import stays" but Plan 02 correctly identified and removed it since its only usage was behind the dead Distributed branch. User should confirm this deviation is acceptable.

### Gaps Summary

No gaps found. All 7 observable truths verified. All artifacts exist, are substantive, and are properly wired. All key links confirmed. The one ROADMAP success criterion that was not achieved (file renaming, criterion 1) was explicitly descoped by the user during context gathering.

The deviation from ROADMAP criterion 3 ("SharedArrays import stays") is an improvement -- SharedArrays was correctly removed because its only usage was within the dead `nprocs() > 1` branches that were cleaned up after Distributed removal. This is documented in the 37-02-SUMMARY.md.

---

_Verified: 2026-02-27T16:05:15Z_
_Verifier: Claude (gsd-verifier)_
