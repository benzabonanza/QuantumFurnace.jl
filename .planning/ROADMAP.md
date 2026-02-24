# Roadmap: QuantumFurnace.jl

## Milestones

- ✅ **v1.0 Trajectories** -- Phases 1-5 (shipped 2026-02-14)
- ✅ **v1.1 Reduce** -- Phases 6-11 (shipped 2026-02-15)
- ✅ **v1.2 Multi-threading** -- Phases 12-19 (shipped 2026-02-16)
- ✅ **v1.3 Mixing Time Estimation** -- Phases 20-25 (shipped 2026-02-19)
- ✅ **v1.4 Spectral Gap Refinement** -- Phase 26 (partial, shipped 2026-02-20)
- 🚧 **v1.5 Krylov Gap Estimation** -- Phases 27-31 (in progress)

## Phases

<details>
<summary>✅ v1.0 Trajectories (Phases 1-5) -- SHIPPED 2026-02-14</summary>

- [x] Phase 1: Foundation and Compilation (1/1 plans) -- completed 2026-02-13
- [x] Phase 2: Trajectory Bug Fixes (2/2 plans) -- completed 2026-02-14
- [x] Phase 3: DM Reference Test Suite (3/3 plans) -- completed 2026-02-14
- [x] Phase 4: Trajectory Cross-Validation (2/2 plans) -- completed 2026-02-14
- [x] Phase 5: Statistical Validation and Regression (2/2 plans) -- completed 2026-02-14

Full details: [milestones/v1.0-ROADMAP.md](milestones/v1.0-ROADMAP.md)

</details>

<details>
<summary>✅ v1.1 Reduce (Phases 6-11) -- SHIPPED 2026-02-15</summary>

- [x] Phase 6: Dead Code Pruning (2/2 plans) -- completed 2026-02-15
- [x] Phase 7: DRY Refactoring (2/2 plans) -- completed 2026-02-15
- [x] Phase 8: Struct Simplification (3/3 plans) -- completed 2026-02-15
- [x] Phase 9: Type Parameterization (3/3 plans) -- completed 2026-02-15
- [x] Phase 10: API Surface Cleanup (3/3 plans) -- completed 2026-02-15
- [x] Phase 11: Allocation Optimization (3/3 plans) -- completed 2026-02-15

Full details: [milestones/v1.1-ROADMAP.md](milestones/v1.1-ROADMAP.md)

</details>

<details>
<summary>✅ v1.2 Multi-threading (Phases 12-19) -- SHIPPED 2026-02-16</summary>

- [x] Phase 12: Workspace Refactor (2/2 plans) -- completed 2026-02-15
- [x] Phase 13: Multi-Threaded Trajectory Engine (2/2 plans) -- completed 2026-02-16
- [x] Phase 14: GNS Trajectory Path (1/1 plans) -- completed 2026-02-16
- [x] Phase 15: Data Architecture (2/2 plans) -- completed 2026-02-16
- [x] Phase 16: Convergence Tracking (2/2 plans) -- completed 2026-02-16
- [x] Phase 17: Adaptive Sampling (2/2 plans) -- completed 2026-02-16
- [x] Phase 18: KMS-vs-GNS Experiments (1/1 plans) -- completed 2026-02-16
- [x] Phase 19: Logic Simplification (3/3 plans) -- completed 2026-02-16

Full details: [milestones/v1.2-ROADMAP.md](milestones/v1.2-ROADMAP.md)

</details>

<details>
<summary>✅ v1.3 Mixing Time Estimation (Phases 20-25) -- SHIPPED 2026-02-19</summary>

- [x] Phase 20: Observable Infrastructure (1/1 plans) -- completed 2026-02-17
- [x] Phase 21: Exponential Fitting (1/1 plans) -- completed 2026-02-17
- [x] Phase 22: Observable-Only Trajectory Runner (1/1 plans) -- completed 2026-02-17
- [x] Phase 23: Gap Estimation API (1/1 plans) -- completed 2026-02-17
- [x] Phase 24: Cross-Validation (3/3 plans) -- completed 2026-02-17
- [x] Phase 25: Spectral Gap Validation Overhaul (3/3 plans) -- completed 2026-02-18

Full details: [milestones/v1.3-ROADMAP.md](milestones/v1.3-ROADMAP.md)

</details>

<details>
<summary>✅ v1.4 Spectral Gap Refinement (Phase 26) -- PARTIAL, SHIPPED 2026-02-20</summary>

- [x] Phase 26: Exact Reference and Structural Diagnostics (2/2 plans) -- completed 2026-02-19

Phases 27-30 deferred (two-exponential fitting, effective rates, bootstrap, dashboard).

Full details: [milestones/v1.4-ROADMAP.md](milestones/v1.4-ROADMAP.md)

</details>

### 🚧 v1.5 Krylov Gap Estimation (In Progress)

**Milestone Goal:** Matrix-free spectral gap estimation via KrylovKit.jl, enabling gap computation for system sizes (up to ~12 qubits) where the full Lindbladian matrix cannot be stored.

- [x] **Phase 27: Core Matvec Infrastructure** - KrylovWorkspace, EnergyDomain matvec, coherent term, adjoint, round-trip correctness gate
- [x] **Phase 28: Domain Extension** - Matrix-free Lindbladian action for TimeDomain, TrotterDomain, and BohrDomain
- [x] **Phase 29: Eigensolver Integration** - KrylovKit eigsolve wrapper, CPTP channel path, result struct, memory guard, convergence handling
- [x] **Phase 30: Cross-Validation** - Krylov gap vs dense eigen() at n=4,6; L-vs-E consistency; KMS-vs-GNS comparison
- [ ] **Phase 31: Scaling Benchmarks** - Timing and memory at n=3-7, 4^n scaling fit, per-matvec breakdown, extrapolation to n=10,12

## Phase Details

### Phase 27: Core Matvec Infrastructure
**Goal**: Users can apply the Lindbladian superoperator to a density matrix without forming the full dim^2 x dim^2 matrix, validated against dense construction at n=4
**Depends on**: Phase 26 (v1.4 diagnostics provide the dense reference)
**Requirements**: MATVEC-01, MATVEC-05, MATVEC-07, MATVEC-08, MATVEC-09
**Success Criteria** (what must be TRUE):
  1. `apply_lindbladian!(ws, rho, config, ham)` returns L(rho) for EnergyDomain with both KMS and GNS configs, matching dense `construct_lindbladian()` to < 1e-12 for 10 random density matrices at n=4
  2. Coherent term -i[B, rho] is included when `with_coherent=true` for KMS configs, and the round-trip test passes with coherent term enabled
  3. Adjoint Lindbladian L'(rho) (dissipator with A <-> A' swap) passes its own round-trip test against dense adjoint at n=4
  4. KrylovWorkspace pre-allocates all scratch matrices at construction time and the matvec hot path produces zero heap allocations (verified by `@allocated`)
**Plans:** 2 plans -- completed 2026-02-20

Plans:
- [x] 27-01-PLAN.md -- KrylovWorkspace struct + constructor, apply_lindbladian!, apply_adjoint_lindbladian!, _accumulate_dissipator!
- [x] 27-02-PLAN.md -- Module integration, test helpers, round-trip correctness tests, allocation regression tests

### Phase 28: Domain Extension
**Goal**: Matrix-free Lindbladian action works for all four domains, with proper NUFFT prefactors (Time), Trotter eigenbasis (Trotter), and Bohr bucket iteration (Bohr)
**Depends on**: Phase 27 (shared dissipator helpers and workspace pattern)
**Requirements**: MATVEC-02, MATVEC-03, MATVEC-04
**Success Criteria** (what must be TRUE):
  1. `apply_lindbladian!` for TimeDomain passes round-trip test against dense at n=4 (KMS + GNS), using NUFFT prefactors from existing `_precompute_data()`
  2. `apply_lindbladian!` for TrotterDomain passes round-trip test at n=4 (KMS + GNS), operating in Trotter eigenbasis (trotter.eigvecs) and producing results consistent with dense construction
  3. `apply_lindbladian!` for BohrDomain passes round-trip test at n=4 (KMS + GNS), using Bohr frequency bucket iteration with the generalized two-operator dissipator
**Plans:** 2 plans -- completed 2026-02-20

Plans:
- [x] 28-01-PLAN.md -- Time/Trotter forward+adjoint matvec with NUFFT prefactors, round-trip and allocation tests
- [x] 28-02-PLAN.md -- BohrDomain forward+adjoint matvec with 2op dissipator, bucket iteration, round-trip and duality tests

### Phase 29: Eigensolver Integration
**Goal**: Users can compute spectral gaps via a single `krylov_spectral_gap()` call that wraps KrylovKit eigsolve, with both Lindbladian (:LR) and CPTP channel (:LM) targeting
**Depends on**: Phase 28 (all domains must have correct matvec before eigsolve)
**Requirements**: MATVEC-06, KRYLOV-01, KRYLOV-02, KRYLOV-03, KRYLOV-04, KRYLOV-05
**Success Criteria** (what must be TRUE):
  1. `krylov_spectral_gap(config, ham)` returns a `KrylovGapResult` containing eigenvalues sorted by real part, the spectral gap, convergence info, matvec count, and the fixed-point density matrix
  2. CPTP channel path `apply_delta_channel!(ws, rho, delta, config, ham)` computes E(rho) = (I + delta*L)(rho) correctly, and `krylov_spectral_gap` with `:LM` targeting finds channel eigenvalues
  3. Pre-flight memory estimation warns when `krylovdim * 4^n * 16 * 1.5` exceeds a configurable threshold, before any computation begins
  4. When KrylovKit reports `info.converged < howmany`, the solver retries with increased krylovdim and issues a warning; if still unconverged, an error is raised
**Plans:** 2 plans -- completed 2026-02-24

Plans:
- [x] 29-01-PLAN.md -- KrylovKit dependency, KrylovGapResult struct, apply_delta_channel!, krylov_spectral_gap (Lindbladian + channel), retry logic, memory guard
- [x] 29-02-PLAN.md -- Tests: apply_delta_channel! round-trip, eigsolve accuracy, channel path, all-domain coverage, guard rails

### Phase 30: Cross-Validation
**Goal**: Krylov spectral gap results are validated against dense eigen() reference values across all domains and balance types, establishing trust for n>6 production use
**Depends on**: Phase 29 (eigensolver must be functional)
**Requirements**: XVAL-01, XVAL-02, XVAL-03, XVAL-04
**Success Criteria** (what must be TRUE):
  1. Krylov gap matches dense eigen() gap to < 1e-8 at n=4 for all 4 domains with KMS balance
  2. Krylov gap matches dense eigen() gap to < 1e-6 at n=6 for all 4 domains with KMS balance
  3. Krylov Lindbladian gap and Krylov channel gap are consistent: gap_L approximately equals -log(|lambda_2(E)|)/delta, within O(delta^2) tolerance
  4. Krylov KMS vs GNS gap comparison at n=4 produces results consistent with existing dense-method gap values
**Plans:** 2 plans -- completed 2026-02-24

Plans:
- [x] 30-01-PLAN.md -- n=4 cross-validation (KMS + GNS all domains), L-vs-E convergence, helper functions, runtests.jl include
- [x] 30-02-PLAN.md -- n=6 env-gated cross-validation (KMS all domains), n=6 test system factories

### Phase 31: Scaling Benchmarks
**Goal**: Empirical timing and memory data at n=3-7 establishes the 4^n scaling law and produces resource estimates for n=10,12 production runs
**Depends on**: Phase 30 (cross-validation must pass before trusting large-n results)
**Requirements**: BENCH-01, BENCH-02, BENCH-03, BENCH-04
**Success Criteria** (what must be TRUE):
  1. Wall-clock timing benchmarks exist for n=3,4,5,6 (n=7 if feasible within 10-minute timeout) with 4 BLAS threads, for at least EnergyDomain KMS
  2. Peak memory usage is measured at each system size, confirming it stays within the pre-flight estimate
  3. A power-law fit to the timing data confirms approximately 4^n scaling, with extrapolated wall-clock estimates for n=10 and n=12
  4. Per-matvec timing breakdown separates BLAS matrix multiplication cost from precomputation lookup and Krylov iteration overhead
**Plans:** 2 plans

Plans:
- [ ] 31-01-PLAN.md -- Complete benchmark script: system factory, timing/memory measurement, krylovdim probe at n=6, scaling fit, report generation
- [ ] 31-02-PLAN.md -- Execute benchmark, verify report output, validate scaling fits and extrapolation

## Progress

**Execution Order:**
Phases execute in numeric order: 27 -> 28 -> 29 -> 30 -> 31

| Phase | Milestone | Plans Complete | Status | Completed |
|-------|-----------|----------------|--------|-----------|
| 1. Foundation and Compilation | v1.0 | 1/1 | Complete | 2026-02-13 |
| 2. Trajectory Bug Fixes | v1.0 | 2/2 | Complete | 2026-02-14 |
| 3. DM Reference Test Suite | v1.0 | 3/3 | Complete | 2026-02-14 |
| 4. Trajectory Cross-Validation | v1.0 | 2/2 | Complete | 2026-02-14 |
| 5. Statistical Validation and Regression | v1.0 | 2/2 | Complete | 2026-02-14 |
| 6. Dead Code Pruning | v1.1 | 2/2 | Complete | 2026-02-15 |
| 7. DRY Refactoring | v1.1 | 2/2 | Complete | 2026-02-15 |
| 8. Struct Simplification | v1.1 | 3/3 | Complete | 2026-02-15 |
| 9. Type Parameterization | v1.1 | 3/3 | Complete | 2026-02-15 |
| 10. API Surface Cleanup | v1.1 | 3/3 | Complete | 2026-02-15 |
| 11. Allocation Optimization | v1.1 | 3/3 | Complete | 2026-02-15 |
| 12. Workspace Refactor | v1.2 | 2/2 | Complete | 2026-02-15 |
| 13. Multi-Threaded Trajectory Engine | v1.2 | 2/2 | Complete | 2026-02-16 |
| 14. GNS Trajectory Path | v1.2 | 1/1 | Complete | 2026-02-16 |
| 15. Data Architecture | v1.2 | 2/2 | Complete | 2026-02-16 |
| 16. Convergence Tracking | v1.2 | 2/2 | Complete | 2026-02-16 |
| 17. Adaptive Sampling | v1.2 | 2/2 | Complete | 2026-02-16 |
| 18. KMS-vs-GNS Experiments | v1.2 | 1/1 | Complete | 2026-02-16 |
| 19. Logic Simplification | v1.2 | 3/3 | Complete | 2026-02-16 |
| 20. Observable Infrastructure | v1.3 | 1/1 | Complete | 2026-02-17 |
| 21. Exponential Fitting | v1.3 | 1/1 | Complete | 2026-02-17 |
| 22. Observable-Only Trajectory Runner | v1.3 | 1/1 | Complete | 2026-02-17 |
| 23. Gap Estimation API | v1.3 | 1/1 | Complete | 2026-02-17 |
| 24. Cross-Validation | v1.3 | 3/3 | Complete | 2026-02-17 |
| 25. Spectral Gap Validation Overhaul | v1.3 | 3/3 | Complete | 2026-02-18 |
| 26. Exact Reference and Structural Diagnostics | v1.4 | 2/2 | Complete | 2026-02-19 |
| 27. Core Matvec Infrastructure | v1.5 | 2/2 | Complete | 2026-02-20 |
| 28. Domain Extension | v1.5 | 2/2 | Complete | 2026-02-20 |
| 29. Eigensolver Integration | v1.5 | 2/2 | Complete | 2026-02-24 |
| 30. Cross-Validation | v1.5 | 2/2 | Complete | 2026-02-24 |
| 31. Scaling Benchmarks | v1.5 | 0/TBD | Not started | - |
