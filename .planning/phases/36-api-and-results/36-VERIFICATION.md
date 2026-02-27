---
phase: 36-api-and-results
verified: 2026-02-27T12:30:43Z
status: passed
score: 7/7 must-haves verified
---

# Phase 36: API and Results Verification Report

**Phase Goal:** Four clean public entry points (`run_lindblad`, `run_thermalize`, `run_krylov_spectrum`, `run_trajectory`) each return a typed Result struct with optional BSON save capability
**Verified:** 2026-02-27T12:30:43Z
**Status:** passed
**Re-verification:** No -- initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | 4 public entry points exist, each dispatching on Config{S,D,C,T} | VERIFIED | `run_lindblad` (furnace.jl:192), `run_thermalize` (furnace.jl:264), `run_krylov_spectrum` (krylov_eigsolve.jl:603), `run_trajectory` (trajectories.jl:1173). All dispatch on Config type parameters. |
| 2 | Each entry point returns a typed Result struct containing config and metadata | VERIFIED | `run_lindblad` returns `LindbladResults{Tc}` (furnace.jl:231), `run_thermalize` returns `ThermalizeResults{Tc}` (furnace.jl:341), `run_krylov_spectrum` returns `KrylovSpectrumResults{Float64}` (krylov_eigsolve.jl:628), `run_trajectory` returns `TrajectoryResults{Float64}` in all 4 code paths (trajectories.jl:1215,1231,1246,1262). All include config and metadata fields populated with `_capture_metadata(wall_time_seconds=...)`. |
| 3 | save_result serializes any Result to BSON with companion .txt | VERIFIED | `save_result(::AbstractResults, ::String)` at results.jl:620 calls `_result_to_dict` (dispatch to 4 type-specific converters) then `BSON.bson` then `_write_result_companion_txt`. Companion .txt writer at results.jl:660 has type-specific sections for all 4 Result types. |
| 4 | load_result round-trips correctly via result_type tag auto-detection | VERIFIED | `load_result(::String)` at results.jl:634 reads `:result_type` tag and dispatches to `_dict_to_*_results` for all 4 types. Round-trip tests in test/test_results.jl:364-575 cover all 4 types (7 testsets). |
| 5 | All 4 Result structs are subtypes of AbstractResults with correct fields | VERIFIED | `abstract type AbstractResults end` at structs.jl:217. `LindbladResults{T} <: AbstractResults` (structs.jl:226), `ThermalizeResults{T} <: AbstractResults` (structs.jl:241), `KrylovSpectrumResults{T} <: AbstractResults` (structs.jl:255), `TrajectoryResults{T} <: AbstractResults` (structs.jl:277). All have metadata::Dict{Symbol,Any} field. TrajectoryResults has Union{Nothing,...} for optional fields. |
| 6 | All new types and functions are exported | VERIFIED | QuantumFurnace.jl:53 exports `AbstractResults, LindbladResults, ThermalizeResults, KrylovSpectrumResults, TrajectoryResults`. Line 54 exports `save_result, load_result`. Line 77 exports `run_lindblad, run_thermalize, run_krylov_spectrum, run_trajectory`. |
| 7 | Simulation scripts demonstrate the new API entry points | VERIFIED | main_liouv.jl:130 uses `run_lindblad`, main_thermalize.jl:130 uses `run_thermalize`, main_krylov_benchmark.jl:614 has commented `run_krylov_spectrum` example. `run_trajectory` has no dedicated script (descoped; no pre-existing main_trajectory.jl). |

**Score:** 7/7 truths verified

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `src/structs.jl` | AbstractResults + 4 Result structs | VERIFIED | Lines 217-286: abstract type + 4 concrete structs with all planned fields |
| `src/results.jl` | save_result, load_result, _result_to_dict, _dict_to_*_results | VERIFIED | Lines 448-736: type tags, to_dict (4 types), from_dict (4 types), save/load, companion .txt, filename gen |
| `src/furnace.jl` | run_lindblad, run_thermalize | VERIFIED | run_lindblad (line 192-239, ~47 lines), run_thermalize (line 264-348, ~84 lines). Both substantive with full simulation logic, timing, and Result construction. |
| `src/krylov_eigsolve.jl` | run_krylov_spectrum | VERIFIED | Line 603-642, ~39 lines. Delegates to krylov_spectral_gap, wraps result into KrylovSpectrumResults with metadata. |
| `src/trajectories.jl` | run_trajectory | VERIFIED | Line 1173-1266, ~93 lines. 4-path keyword dispatch (adaptive/convergence/observable/default), all returning TrajectoryResults. |
| `src/convergence.jl` | _run_trajectory_convergence, _run_trajectory_adaptive | VERIFIED | Lines 408-481. Both delegate to existing convergence functions after computing Gibbs and auto-building observables. |
| `src/QuantumFurnace.jl` | Updated exports | VERIFIED | Lines 52-54 (types + save/load), line 77 (run_* functions) |
| `test/test_results.jl` | Round-trip tests for all 4 types | VERIFIED | Lines 364-575: 7 testsets covering LindbladResults, ThermalizeResults, KrylovSpectrumResults (plain + channel), TrajectoryResults (plain + convergence), metadata exclusion |
| `simulations/main_liouv.jl` | Uses run_lindblad | VERIFIED | Line 130: `run_lindblad(jumps, config, hamiltonian, trotter)` |
| `simulations/main_thermalize.jl` | Uses run_thermalize | VERIFIED | Line 130: `run_thermalize(jumps, config, hamiltonian, trotter; initial_dm=initial_dm)` |
| `simulations/main_krylov_benchmark.jl` | run_krylov_spectrum example | VERIFIED | Line 614: commented example showing usage |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| src/furnace.jl | src/structs.jl | LindbladResults constructor | WIRED | `LindbladResults{Tc}(...)` at furnace.jl:231 |
| src/furnace.jl | src/structs.jl | ThermalizeResults constructor | WIRED | `ThermalizeResults{Tc}(...)` at furnace.jl:341 |
| src/krylov_eigsolve.jl | src/structs.jl | KrylovSpectrumResults constructor | WIRED | `KrylovSpectrumResults{Float64}(...)` at krylov_eigsolve.jl:628 |
| src/trajectories.jl | src/structs.jl | TrajectoryResults constructor | WIRED | `TrajectoryResults{Float64}(...)` at trajectories.jl:1215, 1231, 1246, 1262 |
| src/results.jl | src/structs.jl | AbstractResults used by save/load | WIRED | `save_result(result::AbstractResults, ...)` dispatches on all 4 subtypes via `_result_to_dict` |
| src/trajectories.jl | src/convergence.jl | run_trajectory delegates convergence | WIRED | `_run_trajectory_convergence(...)` at trajectories.jl:1223, `_run_trajectory_adaptive(...)` at trajectories.jl:1205. Helpers defined in convergence.jl:408,448 and resolved via Julia late binding. |
| src/QuantumFurnace.jl | src/results.jl | exports save_result/load_result | WIRED | `export save_result, load_result` at QuantumFurnace.jl:54 |
| test/test_results.jl | src/results.jl | Tests call save_result/load_result | WIRED | 12 calls to save_result/load_result in test_results.jl:364-575 |

### Requirements Coverage

| Requirement | Status | Blocking Issue |
|-------------|--------|----------------|
| WORK-03: Define 4 clean run_* entry points | SATISFIED | All 4 exist with uniform positional signature |
| WORK-04: Define 4 Result structs | SATISFIED | All 4 defined as subtypes of AbstractResults |
| WORK-05: Add save capability with metadata | SATISFIED | save_result/load_result with BSON + companion .txt; metadata includes git hash, timestamp, n_threads, wall_time_seconds |

### Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| src/results.jl | 396 | `julia_version` reference in `_write_companion_txt` for ExperimentResult | Info | Only applies to OLD ExperimentResult companion .txt (not new Result types). No impact on Phase 36 goal -- old code preserved for backward compat. |

No TODO/FIXME/placeholder/stub anti-patterns found in any of the Phase 36 modified files.

### Human Verification Required

### 1. Full Test Suite Execution

**Test:** Run `julia --project=@. -e "using Pkg; Pkg.test()"` in the project directory.
**Expected:** All tests pass (summaries claim 1254 tests). The round-trip serialization tests for all 4 Result types pass.
**Why human:** Cannot run Julia in this verification environment; test execution requires the Julia runtime with all dependencies installed.

### 2. Simulation Script Syntax Validity

**Test:** Run `julia --project=@. -e 'include("simulations/main_liouv.jl")'` (and similarly for main_thermalize.jl).
**Expected:** Scripts parse without errors (they will fail at runtime without Hamiltonian data, but should not have syntax errors).
**Why human:** Requires Julia runtime to parse and check syntax.

### 3. run_trajectory Keyword Dispatch Modes

**Test:** Exercise `run_trajectory` with (a) default mode, (b) observables provided, (c) convergence=true, (d) convergence=true + adaptive=true, using a small system.
**Expected:** Each mode returns a TrajectoryResults with the appropriate fields populated/nothing.
**Why human:** Requires full simulation setup (Hamiltonian, jump operators, initial state) and Julia runtime.

### Gaps Summary

No gaps found. All 7 observable truths are verified. All artifacts exist, are substantive (not stubs), and are wired together. All key links are connected. Requirements WORK-03, WORK-04, WORK-05 are satisfied.

Minor note: Success criterion 4 ("Simulation scripts demonstrate all 4 entry points") is 3/4 -- `run_trajectory` has no dedicated simulation script. This was explicitly descoped in Plan 04 (no pre-existing main_trajectory.jl; creating one would exceed minimal demonstration scope). ORG-07 (create/update 4 simulation scripts) is assigned to Phase 37, not Phase 36, per REQUIREMENTS.md traceability. The 3 updated scripts adequately demonstrate the new API pattern.

---

_Verified: 2026-02-27T12:30:43Z_
_Verifier: Claude (gsd-verifier)_
