---
phase: 10-api-surface-cleanup
verified: 2026-02-15T14:30:00Z
status: passed
score: 5/5 must-haves verified
---

# Phase 10: API Surface Cleanup Verification Report

**Phase Goal:** Public API exposes exactly what users and researchers need -- building blocks for pedagogy are exported, implementation details are internal

**Verified:** 2026-02-15T14:30:00Z
**Status:** passed
**Re-verification:** No -- initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | Non-mutating `oft` wrapper and dead struct names are no longer exported | ✓ VERIFIED | LindbladianJumpCaches, LiouvLiouv, LindbladWorkspace not in export list; oft! kept for pedagogy (explicit comment in export block) |
| 2 | ~18 implementation-detail exports removed from public API | ✓ VERIFIED | OFTCaches, NUFFTPrefactors, KrausScratch, TrajectoryWorkspace, prepare_oft_nufft_prefactors, prefactor_view, precompute_* functions, jump_contribution!, generate_filename all absent from exports |
| 3 | `trace_distance_h` is exported and accessible for convergence analysis workflows | ✓ VERIFIED | `trace_distance_h` in export list (line 46 of QuantumFurnace.jl); accessible via `using QuantumFurnace` without qualification |
| 4 | All 224 existing tests pass with no regressions | ✓ VERIFIED | Test suite passes: "Test Summary: QuantumFurnace.jl \| Pass  224    224  41.6s" |
| 5 | All cross-file call sites use _-prefixed names, no old unprefixed internal function calls remain | ✓ VERIFIED | Verified via grep: precompute_data, vectorize_liouv_diss_and_add!, print_press, generate_filename all return "No old unprefixed calls found" |

**Score:** 5/5 truths verified

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `src/QuantumFurnace.jl` | Reorganized export block with labeled groups, internal exports removed, physics exports added | ✓ VERIFIED | Contains "# --- Public API ---" header with labeled groups (Types: Simulation, QI Tools, Gibbs & Hamiltonian, etc.) |
| `src/furnace.jl` | All cross-file calls updated to _-prefixed names | ✓ VERIFIED | Contains _precompute_data calls at lines 58, 123 |
| `src/furnace_utensils.jl` | All cross-file calls updated to _-prefixed names | ✓ VERIFIED | Contains _truncate_time_labels_for_oft at line 87 |
| `src/jump_workers.jl` | All cross-file calls updated to _-prefixed names | ✓ VERIFIED | Contains _vectorize_liouv_diss_and_add! at lines 42, 74, 77, 84, 119, 122, 130 |
| `test/test_dm_scaling.jl` | Qualified access updated to _-prefixed names | ✓ VERIFIED | Contains QuantumFurnace._create_energy_labels at lines 152, 203 |
| `src/qi_tools.jl` | Internal function definitions prefixed with _ | ✓ VERIFIED | Contains function _kron!, _vectorize_liouv_diss_and_add!, _vectorize_liouvillian_coherent! |
| `src/coherent.jl` | Internal function definitions prefixed with _ | ✓ VERIFIED | Contains function _precompute_coherent_total_B, _compute_b_minus, _compute_b_plus, etc. |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|----|--------|---------|
| src/QuantumFurnace.jl | src/qi_tools.jl | Export trace_distance_h, fidelity, etc. | ✓ WIRED | trace_distance_h found in export list line 46 |
| src/QuantumFurnace.jl | src/bohr_domain.jl | Export create_f, create_alpha_gauss, etc. | ✓ WIRED | create_f found in export list line 57 |
| src/furnace.jl | src/furnace_utensils.jl | Calls _precompute_data instead of precompute_data | ✓ WIRED | _precompute_data called at lines 58, 123 in furnace.jl |
| src/jump_workers.jl | src/qi_tools.jl | Calls _vectorize_liouv_diss_and_add! | ✓ WIRED | _vectorize_liouv_diss_and_add! called at 7 sites in jump_workers.jl |
| src/furnace_utensils.jl | src/coherent.jl | Calls _compute_b_minus, _compute_b_plus, etc. | ✓ WIRED | _compute_b_minus called via _compute_truncated_func at line 94 |
| test/test_dm_scaling.jl | src/energy_domain.jl | Qualified access QuantumFurnace._create_energy_labels | ✓ WIRED | Pattern found at lines 152, 203 |

### Requirements Coverage

Phase 10 addresses requirements API-01, API-02, and API-03 from REQUIREMENTS.md.

| Requirement | Status | Blocking Issue |
|-------------|--------|----------------|
| API-01: Remove dead/deprecated exports (non-mutating oft, dead structs) | ✓ SATISFIED | LindbladianJumpCaches, LiouvLiouv, LindbladWorkspace removed in Phase 6, confirmed absent from export list. oft! kept for pedagogy (documented). |
| API-02: Internalize implementation-detail exports (~18 items) | ✓ SATISFIED | All workspace types (OFTCaches, NUFFTPrefactors, KrausScratch, TrajectoryWorkspace), precompute helpers (prepare_oft_nufft_prefactors, prefactor_view, precompute_data, etc.), and internal dispatch functions (jump_contribution!, generate_filename) removed from exports. All ~45 internal functions renamed with _ prefix. |
| API-03: Export trace_distance_h for convergence analysis | ✓ SATISFIED | trace_distance_h exported in QI Tools group; accessible via `using QuantumFurnace` without qualification. |

### Anti-Patterns Found

No anti-patterns found in Phase 10 modified files.

**Pre-existing TODOs** (not blockers, outside Phase 10 scope):
- src/jump_workers.jl:454 - `#TODO: test it; set BLAS threads to 1, let julia threads be more.` (threading optimization note)
- src/coherent.jl:280 - `#TODO: Reintroduce sigmas here` (future enhancement)

These are pre-existing technical debt comments not introduced by Phase 10 work.

### Human Verification Required

None. All verifiable items passed automated checks:
- Export list reorganization: programmatically verified
- Function accessibility: tested via Julia REPL
- Internal function renaming: verified via grep and module loading
- Test passage: verified via `Pkg.test()`
- Aqua.jl quality checks: included in test suite, passed

---

## Phase 10 Execution Summary

Phase 10 was executed across 3 plans:

**10-01: Export Block Reorganization** (Plan 01)
- Reorganized export block into labeled groups (Types: Simulation, QI Tools, etc.)
- Added 10+ new physics function exports (trace_distance_h, create_f, trotterize, etc.)
- Removed 12+ implementation-detail exports (workspaces, precompute helpers, internal dispatch)
- Updated tests to use unqualified trace_distance_h (now exported)
- Commits: 9eecd80 (export reorganization), 6c2243e (test updates)

**10-02: Internal Function Prefix Rename** (Plan 02)
- Renamed ~45 internal function definitions to _-prefixed names across 14 source files
- Updated all intra-file call sites to use new _-prefixed names
- Preserved all exported function names (no changes)
- Preserved all type names (no changes)
- Commits: f0739b7 (10 self-contained files), 05734d6 (4 worker/integration files)

**10-03: Cross-File Call Site Update** (Plan 03)
- Updated all cross-file call sites in 5 source files to _-prefixed names (~30 call sites)
- Updated test qualified access patterns in 7 test files
- Fixed trajectory_validation standalone scripts (obsolete 3-arg precompute_data -> 2-arg)
- Updated coherent.jl docstrings to match _-prefixed function names
- Commits: 6b8ee6e (source file updates), f5b29a4 (test file updates)

**All 224 tests pass** including Aqua.jl quality checks (no undefined exports, no ambiguities flagged).

---

## Detailed Verification Results

### 1. Export List Verification (API-01, API-02, API-03)

**Dead exports (API-01):** VERIFIED ABSENT
```bash
$ grep -E "export.*LindbladianJumpCaches|export.*LiouvLiouv|export.*LindbladWorkspace" src/QuantumFurnace.jl
# No matches - confirmed removed
```

**Implementation-detail exports (API-02):** VERIFIED REMOVED
```bash
$ grep -E "export.*OFTCaches|export.*NUFFTPrefactors|export.*precompute_" src/QuantumFurnace.jl
# No matches - confirmed removed
```

**Physics exports (API-03):** VERIFIED PRESENT
```bash
$ grep "export.*trace_distance_h" src/QuantumFurnace.jl
export trace_distance_h, trace_distance_nh, trace_norm_h, trace_norm_nh,
```

### 2. Internal Function Naming Convention

**Sample verification of _-prefixed definitions:**
- src/qi_tools.jl: `function _kron!`, `function _vectorize_liouv_diss_and_add!`, `function _vectorize_liouvillian_coherent!`
- src/hamiltonian.jl: `function _construct_base_ham`, `function _construct_disordering_terms`, `function _rescaling_and_shift_factors`
- src/coherent.jl: `function _precompute_coherent_total_B`, `function _compute_b_minus`, `function _compute_b_plus`
- src/furnace_utensils.jl: `function _precompute_data`, `function _precompute_labels`, `function _select_b_plus_calculator`

**Verified no old unprefixed calls remain:**
```bash
$ grep -rn "\bprecompute_data\b" src/ --include="*.jl" | grep -v "function\|#\|_precompute_data"
# No matches

$ grep -rn "\bvectorize_liouv_diss_and_add!\b" src/ --include="*.jl" | grep -v "function\|#\|_vectorize_liouv_diss_and_add!"
# No matches

$ grep -rn "\bprint_press\b" src/ --include="*.jl" | grep -v "function\|#\|_print_press"
# No matches

$ grep -rn "\bgenerate_filename\b" src/ --include="*.jl" | grep -v "function\|#\|_generate_filename"
# No matches
```

### 3. Cross-File Wiring Verification

**furnace.jl → furnace_utensils.jl (_precompute_data):**
```bash
$ grep -n "_precompute_data" src/furnace.jl
58:    precomputed_data = _precompute_data(config, ham_or_trott)
123:    precomputed_data = _precompute_data(config, ham_or_trott)
```

**jump_workers.jl → qi_tools.jl (_vectorize_liouv_diss_and_add!):**
```bash
$ grep -n "_vectorize_liouv_diss_and_add!" src/jump_workers.jl
42:        _vectorize_liouv_diss_and_add!(L_target, alpha_A_nu1, A_nu_2_dag, gamma_norm_factor, ws)
74:            _vectorize_liouv_diss_and_add!(L_target, jump_oft, scalar_w, ws)
77:                _vectorize_liouv_diss_and_add!(L_target, jump_oft', scalar_negative_w, ws)
84:            _vectorize_liouv_diss_and_add!(L_target, jump_oft, scalar_w, ws)
119:            _vectorize_liouv_diss_and_add!(L_target, jump_oft, scalar_w, ws)
122:                _vectorize_liouv_diss_and_add!(L_target, jump_oft', scalar_negative_w, ws)
130:            _vectorize_liouv_diss_and_add!(L_target, jump_oft, scalar_w, ws)
```

**furnace_utensils.jl → coherent.jl (_compute_b_minus):**
```bash
$ grep -n "_compute_b_minus" src/furnace_utensils.jl
94:        _b_minus = _compute_truncated_func(_compute_b_minus, time_labels, config.beta, config.sigma)
```

**test/test_dm_scaling.jl → energy_domain.jl (QuantumFurnace._create_energy_labels):**
```bash
$ grep -n "QuantumFurnace\._create_energy_labels" test/test_dm_scaling.jl
152:    energy_labels = QuantumFurnace._create_energy_labels(NUM_ENERGY_BITS, W0)
203:    energy_labels = QuantumFurnace._create_energy_labels(NUM_ENERGY_BITS, W0)
```

### 4. Module Loading and Test Execution

**Module loads successfully:**
```bash
$ julia --project -e 'using QuantumFurnace; println("Module loaded successfully")'
Module loaded successfully
```

**Public API functions accessible without qualification:**
```bash
$ julia --project -e 'using QuantumFurnace; println("trace_distance_h available: ", isdefined(Main, :trace_distance_h)); println("create_f available: ", isdefined(Main, :create_f)); println("trotterize available: ", isdefined(Main, :trotterize))'
trace_distance_h available: true
create_f available: true
trotterize available: true
```

**Internal functions require qualification:**
```bash
$ julia --project -e 'using QuantumFurnace; try; precompute_data; println("FAIL: precompute_data accessible"); catch; println("PASS: precompute_data not accessible"); end; println("_precompute_data via qualified: ", isdefined(QuantumFurnace, :_precompute_data))'
PASS: precompute_data not accessible
_precompute_data via qualified: true
```

**All 224 tests pass:**
```bash
$ julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -5
Test Summary:     | Pass  Total   Time
QuantumFurnace.jl |  224    224  41.6s
     Testing QuantumFurnace tests passed
```

### 5. Commit Verification

All commits documented in SUMMARYs exist and match stated changes:

**Phase 10-01:**
- 9eecd80: Export block reorganization
- 6c2243e: Test updates for unqualified trace_distance_h

**Phase 10-02:**
- f0739b7: Rename internal functions in 10 self-contained files
- 05734d6: Rename internal functions in 4 worker/integration files

**Phase 10-03:**
- 6b8ee6e: Update cross-file call sites to _-prefixed names
- f5b29a4: Update test qualified access to _-prefixed names

All commits verified via `git show --stat`.

---

## Conclusion

**Phase 10 goal ACHIEVED:**

✓ Public API exposes exactly what users and researchers need
✓ Building blocks for pedagogy are exported (trace_distance_h, create_f, trotterize, hermitianize!, etc.)
✓ Implementation details are internal (_-prefixed functions, non-exported workspaces)
✓ All 224 tests pass with no regressions
✓ Aqua.jl quality checks pass

The API surface is now clean and well-organized:
- 60+ physics functions exported (simulation, QI tools, Hamiltonian building blocks, Trotter/Pauli helpers)
- 45+ internal functions hidden behind _ prefix convention
- 12+ workspace/precompute helpers removed from public API
- Export list organized into labeled groups for discoverability

Phase 10 is complete and ready to proceed to Phase 11 (Allocation Optimization).

---

_Verified: 2026-02-15T14:30:00Z_
_Verifier: Claude (gsd-verifier)_
