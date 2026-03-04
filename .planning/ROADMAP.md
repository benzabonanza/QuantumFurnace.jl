# Roadmap: QuantumFurnace.jl

## Milestones

- ✅ **v1.0 Trajectories** -- Phases 1-5 (shipped 2026-02-14)
- ✅ **v1.1 Reduce** -- Phases 6-11 (shipped 2026-02-15)
- ✅ **v1.2 Multi-threading** -- Phases 12-19 (shipped 2026-02-16)
- ✅ **v1.3 Mixing Time Estimation** -- Phases 20-25 (shipped 2026-02-19)
- ✅ **v1.4 Spectral Gap Refinement** -- Phase 26 (partial, shipped 2026-02-20)
- ✅ **v1.5 Krylov Gap Estimation** -- Phases 27-32 (shipped 2026-02-25)
- ✅ **v2.0 Restructure** -- Phases 33-38 (shipped 2026-02-28)
- ✅ **v2.1 Speedup & Mixing Time** -- Phases 39-43 (shipped 2026-03-04)
- **v2.2 Ham Sim Time Counting** -- Phases 44-47 (in progress)

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

<details>
<summary>✅ v1.5 Krylov Gap Estimation (Phases 27-32) -- SHIPPED 2026-02-25</summary>

- [x] Phase 27: Core Matvec Infrastructure (2/2 plans) -- completed 2026-02-20
- [x] Phase 28: Domain Extension (2/2 plans) -- completed 2026-02-20
- [x] Phase 29: Eigensolver Integration (2/2 plans) -- completed 2026-02-24
- [x] Phase 30: Cross-Validation (2/2 plans) -- completed 2026-02-24
- [x] Phase 31: Scaling Benchmarks (2/2 plans) -- completed 2026-02-24
- [x] Phase 32: Krylov Simulator Speedup (2/2 plans) -- completed 2026-02-25

Full details: [milestones/v1.5-ROADMAP.md](milestones/v1.5-ROADMAP.md)

</details>

<details>
<summary>✅ v2.0 Restructure (Phases 33-38) -- SHIPPED 2026-02-28</summary>

- [x] Phase 33: Type Foundation (4/4 plans) -- completed 2026-02-25
- [x] Phase 34: Code Deduplication (2/2 plans) -- completed 2026-02-26
- [x] Phase 35: Workspace and Channel Consolidation (2/2 plans) -- completed 2026-02-27
- [x] Phase 36: API and Results (4/4 plans) -- completed 2026-02-27
- [x] Phase 37: File Organization and Dead Code (2/2 plans) -- completed 2026-02-27
- [x] Phase 38: Test Cleanup (5/5 plans) -- completed 2026-02-28

Full details: [milestones/v2.0-ROADMAP.md](milestones/v2.0-ROADMAP.md)

</details>

<details>
<summary>✅ v2.1 Speedup & Mixing Time (Phases 39-43) -- SHIPPED 2026-03-04</summary>

- [x] Phase 39: Per-Jump Precomputation (2/2 plans) -- completed 2026-03-01
- [x] Phase 40: Save Every (1/1 plans) -- completed 2026-03-01
- [x] Phase 41: Threading (2/2 plans) -- completed 2026-03-01
- [x] Phase 42: Mixing Time Estimation (2/2 plans) -- completed 2026-03-01
- [x] Phase 43: Bi-Exponential Fitting (1/1 plan) -- completed 2026-03-04

Full details: [milestones/v2.1-ROADMAP.md](milestones/v2.1-ROADMAP.md)

</details>

## v2.2 Ham Sim Time Counting

**Milestone Goal:** Build functions that count the total Hamiltonian simulation time needed for Chen's quantum Gibbs sampling algorithm to reach epsilon-distance from the Gibbs state, using QPE-parameter-based time grids.

- [ ] **Phase 44: Struct and Grid Infrastructure** - Result container and QPE grid utilities
- [ ] **Phase 45: OFT Time Counting** - Operator Fourier transform Hamiltonian simulation time from QPE grid
- [ ] **Phase 46: B Coherent Term Counting** - B term double-sum Hamiltonian simulation time with KMS/GNS dispatch
- [ ] **Phase 47: Integration API and Validation** - Public API, convenience overloads, cost assembly, and end-to-end validation

## Phase Details

### Phase 44: Struct and Grid Infrastructure
**Goal**: Users have a well-defined output container and QPE grid utilities that all subsequent cost computations build upon
**Depends on**: Nothing (first phase of v2.2)
**Requirements**: GRID-01, GRID-02, COST-04
**Success Criteria** (what must be TRUE):
  1. `SimulationTimeBudget` immutable struct exists with per-component fields (oft_time, b_time, per_step_time, n_steps, total_time) and grid info (r, N, w0, t0, energy_range)
  2. `_qpe_grid_info(r, w0)` returns correct N=2^r grid size, t0 satisfying the Fourier relation t0*w0=2*pi/N, and energy range
  3. Unit tests confirm grid arithmetic: N=2^r for various r values, Fourier relation holds to machine precision
  4. Struct constructor enforces immutability and stores all parameters needed for paper resource estimation tables
**Plans**: TBD

### Phase 45: OFT Time Counting
**Goal**: Users can compute the OFT Hamiltonian simulation time -- the dominant cost component -- with closed-form validation and all three transition weight functions
**Depends on**: Phase 44
**Requirements**: OFT-01, OFT-02, OFT-03, OFT-04, VAL-01
**Success Criteria** (what must be TRUE):
  1. `_oft_hamiltonian_time` computes the sum of |t_k| over the full 2^r QPE grid (not the truncated simulation grid) and the result is independent of sigma
  2. OFT time counting works with Gaussian, Metropolis (a=0), and smooth Metropolis (a>0, b>0) transition weight functions
  3. Closed-form sanity check passes: unweighted QPE grid sum equals t0*N^2/4 to machine precision
  4. The inner QPE time sum does not apply filter function weights (changing sigma does not change the time sum)
**Plans**: TBD

### Phase 46: B Coherent Term Counting
**Goal**: Users can compute the B coherent term Hamiltonian simulation time as a double sum over truncated b-dicts, with correct KMS/GNS dispatch
**Depends on**: Phase 45
**Requirements**: BTERM-01, BTERM-02, BTERM-03, VAL-03
**Success Criteria** (what must be TRUE):
  1. `_b_hamiltonian_time` computes the double sum over b_minus (outer) and b_plus (inner) truncated dictionaries, with separate inner and outer contributions
  2. B term returns exactly 0.0 for GNS construction (no coherent term)
  3. B term works with all three transition weight functions via correct b_plus variant dispatch
  4. B term result increases with beta (larger beta means more Hamiltonian simulation time for the coherent correction)
**Plans**: TBD

### Phase 47: Integration API and Validation
**Goal**: Users can call a single function to get a complete Hamiltonian simulation time budget from HamHam + scalar parameters, with convenience overloads for MixingTimeEstimate chaining
**Depends on**: Phase 46
**Requirements**: COST-01, COST-02, COST-03, API-01, API-02, API-03, API-04, VAL-02, VAL-04
**Success Criteria** (what must be TRUE):
  1. `compute_simulation_time(ham, r, delta, mixing_time; ...)` returns a `SimulationTimeBudget` with per_step_time = 2*oft_time + b_time (not 2*(oft_time + b_time))
  2. Total step count is ceil(mixing_time / delta) and total_time = n_steps * per_step_time
  3. HamHam convenience overload extracts rescaling_factor and n_qubits automatically; MixingTimeEstimate overload extracts mixing_time from the estimate
  4. Validation test cases with r=12 on n=3,4 systems produce finite, positive, physically plausible results
  5. Module exports `SimulationTimeBudget` and `compute_simulation_time`; `src/simulation_time.jl` included in main module
**Plans**: TBD

## Progress

**Execution Order:**
Phases execute in numeric order: 44 -> 45 -> 46 -> 47

| Phase | Milestone | Plans Complete | Status | Completed |
|-------|-----------|----------------|--------|-----------|
| 1-5 | v1.0 | 10/10 | Complete | 2026-02-14 |
| 6-11 | v1.1 | 16/16 | Complete | 2026-02-15 |
| 12-19 | v1.2 | 15/15 | Complete | 2026-02-16 |
| 20-25 | v1.3 | 10/10 | Complete | 2026-02-18 |
| 26 | v1.4 | 2/2 | Complete | 2026-02-19 |
| 27-32 | v1.5 | 12/12 | Complete | 2026-02-25 |
| 33-38 | v2.0 | 19/19 | Complete | 2026-02-28 |
| 39-43 | v2.1 | 8/8 | Complete | 2026-03-04 |
| 44. Struct and Grid Infrastructure | v2.2 | 0/TBD | Not started | - |
| 45. OFT Time Counting | v2.2 | 0/TBD | Not started | - |
| 46. B Coherent Term Counting | v2.2 | 0/TBD | Not started | - |
| 47. Integration API and Validation | v2.2 | 0/TBD | Not started | - |
