---
phase: 41-threading
verified: 2026-03-01T12:18:05Z
status: passed
score: 5/5 must-haves verified
---

# Phase 41: Threading Verification Report

**Phase Goal:** Users get multi-threaded BLAS speedup on DM thermalization automatically, with optional omega-loop parallelism for large frequency grids, without data races or BLAS thread state leaks
**Verified:** 2026-03-01T12:18:05Z
**Status:** passed
**Re-verification:** No -- initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | run_thermalize with multi-threaded BLAS produces results matching single-threaded execution within expected floating-point tolerance (atol < 1e-10 for accumulated multi-step) | VERIFIED | furnace.jl:204-241 wraps hot loop with `old_blas = BLAS.get_num_threads()` / `BLAS.set_num_threads(Threads.nthreads())` / `try ... finally BLAS.set_num_threads(old_blas)`. Test "DM serial-threaded BLAS agreement" (test_threading.jl:222-248) runs serial-BLAS vs multi-BLAS with same RNG seed and asserts `isapprox(atol=1e-10)` on both trace_distances and final_dm |
| 2 | BLAS thread count is always restored after run_thermalize returns, even if an error occurs (try/finally pattern verified) | VERIFIED | furnace.jl:204 saves `old_blas`, line 205 opens `try`, line 239-241 has `finally BLAS.set_num_threads(old_blas) end`. Three omega-loop threaded variants (jump_workers.jl:699-709, 783-793, 862-872) each have their own BLAS=1 save/restore try/finally. Three _precompute_R threaded variants (trajectories.jl:414-424, 497-507, 574-584) also have BLAS=1 save/restore. Test "DM BLAS thread restoration" (test_threading.jl:200-219) verifies across Energy, Time, and Trotter domains |
| 3 | Threaded omega-loop rho_jump accumulation (when enabled) produces results matching serial accumulation within floating-point tolerance | VERIFIED | jump_workers.jl:492-494 (EnergyDomain), 555-557 (TimeDomain/TrotterDomain), 623-624 (BohrDomain) each dispatch to threaded variant when nthreads > 1 && count >= OMEGA_THREAD_THRESHOLD(50). Per-task ThermalizeScratch provides data race isolation (each spawned task gets private rho_jump/jump_oft/sandwich_tmp). Test "Omega-loop threading determinism" (test_threading.jl:274-300) verifies all 4 domains. Test "Serial vs threaded omega-loop agreement" (test_threading.jl:302-336) verifies atol=1e-10 |
| 4 | BohrDomain benefits from threading where applicable (Bohr bucket iteration, sandwich accumulation) | VERIFIED | jump_workers.jl:840-882 implements `_accumulate_rho_jump_threaded_bohr!` with per-task scratch, BLAS=1, partition of bohr_keys, and chunk worker `_accumulate_rho_jump_chunk_bohr!` (lines 884-936). trajectories.jl:552-640 implements `_precompute_R_threaded_bohr!` with identical pattern. Both gate on `nthreads > 1 && n_keys >= OMEGA_THREAD_THRESHOLD`. Test "Omega-loop threading determinism" includes BohrDomain explicitly (test_threading.jl:280) |
| 5 | No BLAS thread oversubscription occurs when run_thermalize is called independently or in sequence with trajectory functions | VERIFIED | run_thermalize sets BLAS=nthreads for BLAS-level parallelism (furnace.jl:206), while omega-loop threaded variants set BLAS=1 during Julia-level parallelism (jump_workers.jl:700, 784, 863; trajectories.jl:415, 498, 575). These are mutually exclusive: BLAS-level parallelism in the outer hot loop, BLAS=1 during inner omega-loop spawns. Test "DM-trajectory BLAS isolation" (test_threading.jl:250-268) verifies sequential DM + trajectory calls restore BLAS state. Test "Omega-loop BLAS restoration" (test_threading.jl:338-345) verifies omega-loop path restores BLAS |

**Score:** 5/5 truths verified

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `src/furnace.jl` | BLAS try/finally save/restore wrapping run_thermalize hot loop | VERIFIED | Lines 204-241: `old_blas = BLAS.get_num_threads()`, `try`, `BLAS.set_num_threads(Threads.nthreads())`, hot loop, `finally BLAS.set_num_threads(old_blas) end`. Pattern matches plan spec exactly |
| `src/jump_workers.jl` | Threaded _accumulate_rho_jump! for all 4 domains, OMEGA_THREAD_THRESHOLD, _partition_range | VERIFIED | OMEGA_THREAD_THRESHOLD=50 (line 446), _partition_range (lines 454-467), threaded variants for Energy (671-753), TimeTrot (757-836), Bohr (840-936). Each has per-task ThermalizeScratch, BLAS=1, try/finally, reduce. Serial fallback via gate check |
| `src/trajectories.jl` | Threaded _precompute_R for Energy, TimeTrot, Bohr domains | VERIFIED | Threaded variants for Energy (386-467), TimeTrot (471-548), Bohr (552-640). Each called from corresponding _precompute_R serial function via gate (lines 204, 265, 341). Same pattern: per-task scratch, BLAS=1, try/finally, reduce |
| `test/test_threading.jl` | DM BLAS restoration, serial-threaded agreement, DM-trajectory isolation, omega-loop determinism, omega-loop agreement, omega-loop BLAS restoration | VERIFIED | 6 new test blocks (lines 200-345): DM BLAS thread restoration, DM serial-threaded BLAS agreement, DM-trajectory BLAS isolation, Omega-loop threading determinism, Serial vs threaded omega-loop agreement, Omega-loop BLAS restoration. All test substantive assertions (not stubs) |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `src/furnace.jl` | `LinearAlgebra.BLAS` | `BLAS.get_num_threads/set_num_threads with try/finally` | WIRED | furnace.jl:204 `BLAS.get_num_threads()`, :206 `BLAS.set_num_threads(Threads.nthreads())`, :240 `BLAS.set_num_threads(old_blas)` inside try/finally block |
| `src/jump_workers.jl` | `src/furnace.jl` | `_accumulate_rho_jump! called from run_thermalize hot loop` | WIRED | furnace.jl:219-222 calls `_accumulate_rho_jump!(scratch, evolving_dm, jump, ham_or_trott, config, precomputed_data; ...)` which dispatches to domain-specific variants in jump_workers.jl |
| `src/trajectories.jl` | `src/trajectories.jl` | `_precompute_R called from _precompute_per_jump_channels` | WIRED | trajectories.jl:113 calls `_precompute_R([jumps_for_diss[a]], ham_or_trott, config, precomputed_data, builder_scratch)` which dispatches to domain-specific _precompute_R functions with threading gates at :204, :265, :341 |
| `src/jump_workers.jl` (OMEGA_THREAD_THRESHOLD, _partition_range) | `src/trajectories.jl` | Module include order in QuantumFurnace.jl | WIRED | QuantumFurnace.jl:93 includes jump_workers.jl before :94 includes trajectories.jl, so constants and functions from jump_workers.jl are available in trajectories.jl |

### Requirements Coverage

| Requirement | Status | Blocking Issue |
|-------------|--------|----------------|
| THREAD-01: Dissipative sandwich omega-loop runs multi-threaded with per-task accumulators for Energy, Time, Trotter | SATISFIED | All three domain variants implemented in jump_workers.jl with per-task ThermalizeScratch |
| THREAD-02: R^a precomputation omega-loop runs multi-threaded with per-task accumulators | SATISFIED | All three domain families (Energy, TimeTrot, Bohr) have threaded _precompute_R in trajectories.jl |
| THREAD-03: Multi-threaded BLAS enabled during DM thermalization | SATISFIED | furnace.jl:206 `BLAS.set_num_threads(Threads.nthreads())` inside run_thermalize try block |
| THREAD-04: BohrDomain multi-threaded where beneficial | SATISFIED | _accumulate_rho_jump_threaded_bohr! (jump_workers.jl:840-882) and _precompute_R_threaded_bohr! (trajectories.jl:552-594) both implemented |
| THREAD-05: BLAS thread control follows try/finally save/restore pattern | SATISFIED | run_thermalize (furnace.jl:204-241), all 3 omega-loop threaded variants (jump_workers.jl), all 3 _precompute_R threaded variants (trajectories.jl) each have BLAS save/try/finally/restore |

### Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| test/test_threading.jl | 31 | `@test true  # placeholder so testset is not empty` | Info | Standard Julia testing pattern for skipped testsets (when nthreads==1). Not a code stub -- the tests are fully implemented and guarded by nthreads > 1 check |

No TODO/FIXME/HACK/PLACEHOLDER found in any source files (furnace.jl, jump_workers.jl, trajectories.jl). No empty implementations, no console-log-only handlers, no stub returns.

### Human Verification Required

### 1. Multi-threaded BLAS speedup for large systems

**Test:** Run run_thermalize with a 6+ qubit system (dim >= 64) and compare wall-clock time between `-t 1` and `-t 4` Julia invocations.
**Expected:** Multi-threaded BLAS should provide observable speedup (2-4x) for dim >= 64 where BLAS gemm threading becomes effective.
**Why human:** Speedup magnitude depends on hardware, system size, and BLAS implementation. Cannot be verified via static analysis.

### 2. Omega-loop threading speedup for large frequency grids

**Test:** Run run_thermalize with a system that generates > 50 energy labels and compare timing with OMEGA_THREAD_THRESHOLD set high (no threading) vs default (threading active).
**Expected:** Measurable speedup for large frequency grids when omega-loop parallelism activates.
**Why human:** Whether the threshold triggers depends on system parameters, and speedup depends on hardware.

### Gaps Summary

No gaps found. All 5 observable truths verified against the actual codebase. All artifacts exist, are substantive (not stubs), and are properly wired. All 5 THREAD requirements are satisfied. All commits exist and correspond to the claimed changes.

The implementation follows a clean pattern:
- **BLAS threading (Plan 01):** run_thermalize hot loop sets BLAS threads high for BLAS-level parallelism, restores on exit via try/finally.
- **Omega-loop threading (Plan 02):** _accumulate_rho_jump! and _precompute_R for all 4 domains dispatch to threaded variants when nthreads > 1 and work count >= 50. Per-task ThermalizeScratch provides data race isolation. BLAS set to 1 during Julia-level parallelism (mutual exclusion with BLAS-level threading).
- **Tests:** 6 new test blocks covering BLAS restoration, serial-threaded agreement, DM-trajectory isolation, omega-loop determinism, and omega-loop BLAS restoration.

---

_Verified: 2026-03-01T12:18:05Z_
_Verifier: Claude (gsd-verifier)_
