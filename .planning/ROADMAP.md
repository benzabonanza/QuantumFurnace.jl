# Roadmap: QuantumFurnace.jl v1.0 Trajectories

## Overview

This milestone takes trajectory simulation from broken (compilation errors, ordering bugs) to validated (matches density matrix evolution within statistical noise). The work flows through five phases: set up test infrastructure and fix the blocking compilation bug, fix the remaining trajectory code bugs with CPTP verification, establish DM reference tests as ground truth, cross-validate trajectory averages against DM for three approximation domains, and lock down results with statistical convergence tests and regression data.

## Phases

**Phase Numbering:**
- Integer phases (1, 2, 3): Planned milestone work
- Decimal phases (2.1, 2.2): Urgent insertions (marked with INSERTED)

Decimal phases appear between their surrounding integers in numeric order.

- [x] **Phase 1: Foundation and Compilation** - Test infrastructure setup and trajectory compilation fix
- [x] **Phase 2: Trajectory Bug Fixes** - Fix trajectory code bugs and verify CPTP channel correctness
- [ ] **Phase 3: DM Reference Test Suite** - Establish density matrix ground truth across all domains
- [ ] **Phase 4: Trajectory Cross-Validation** - Validate trajectory averages match DM evolution per domain
- [ ] **Phase 5: Statistical Validation and Regression** - Convergence properties and frozen reference data

## Phase Details

### Phase 1: Foundation and Compilation
**Goal**: Trajectory code compiles and test infrastructure exists for all subsequent phases
**Depends on**: Nothing (first phase)
**Requirements**: TINF-04, TINF-01, TFIX-01
**Success Criteria** (what must be TRUE):
  1. `using QuantumFurnace` loads without errors and `build_trajectoryframework` can be called without compilation failures
  2. `Pkg.test()` runs a test suite that includes shared fixtures (`make_test_system`, tiered tolerance constants matching error hierarchy)
  3. Project.toml has test-only dependencies (StableRNGs, HypothesisTests, StatsBase, Aqua) in `[extras]` section, not polluting production `[deps]`
**Plans:** 1 plan

Plans:
- [x] 01-01-PLAN.md — Fix compilation bugs, clean Project.toml, create test infrastructure with fixtures and smoke tests

### Phase 2: Trajectory Bug Fixes
**Goal**: Trajectory simulation runs correctly with proper jump sampling, normalization guards, and CPTP channel verification
**Depends on**: Phase 1
**Requirements**: TFIX-02, TFIX-03, TFIX-04, TFIX-05, TVAL-01
**Success Criteria** (what must be TRUE):
  1. Coherent unitary U_B is applied after branch selection in `step_along_trajectory!`, matching DM code ordering
  2. Running a single trajectory step triggers a normalization assertion that verifies p_nojump + p_res + p_jump_total is approximately 1.0
  3. `build_trajectoryframework` with a non-PSD S matrix does not crash on Cholesky but falls back gracefully (eigenvalue guard)
  4. Trajectory jump sampling faithfully implements Chen's weak measurement scheme, verified by comparing channel structure against DM code (`jump_contribution!`) and the paper's algorithm
  5. CPTP verification test confirms K0*K0 + delta*R + U_res*U_res = I to machine precision for Energy, Time, and Trotter domains
**Plans:** 2 plans

Plans:
- [x] 02-01-PLAN.md — Fix trajectory bugs (U_B ordering, normalization warning, PSD guard, jump sampling verification) and create single-step tests
- [x] 02-02-PLAN.md — CPTP completeness verification test across all three domains

### Phase 3: DM Reference Test Suite
**Goal**: Density matrix simulation has a comprehensive correctness test suite establishing ground truth for all approximation domains
**Depends on**: Phase 1
**Requirements**: DMTST-01, DMTST-02, DMTST-03, DMTST-04, DMTST-05, DMTST-06, TINF-03
**Success Criteria** (what must be TRUE):
  1. BohrDomain with coherent term B produces the exact Gibbs state as its fixed point (trace distance < 1e-10) for a 3-qubit Heisenberg system
  2. Domain error hierarchy is verified: dist(Gibbs, Bohr) <= dist(Gibbs, Energy) <= dist(Gibbs, Time) <= dist(Gibbs, Trotter) for matched parameters
  3. Single DM step error scales as O(delta^2) and multi-step error accumulates as O(delta) over full evolution, verified empirically with a delta-sweep
  4. Coherent term B is consistent across domains: B_bohr vs B_time match up to time quadrature errors, B_trotter matches with additional Trotter errors
  5. OFT consistency verified: `oft!()` and `time_oft!()` produce matching results up to time quadrature errors; `trotter_oft!()` matches with additional Trotter errors
**Plans:** 3 plans

Plans:
- [ ] 03-01-PLAN.md — Detailed balance (DMTST-01) and domain error hierarchy (DMTST-02) tests with 3-qubit fixture
- [ ] 03-02-PLAN.md — DM step error scaling (DMTST-03/04), coherent term B consistency (DMTST-05), and OFT consistency (DMTST-06)
- [ ] 03-03-PLAN.md — Aqua.jl package quality checks (TINF-03)

### Phase 4: Trajectory Cross-Validation
**Goal**: Trajectory-averaged density matrix matches DM evolution for Energy, Time, and Trotter domains, and coherent term produces correct Gibbs convergence
**Depends on**: Phase 2, Phase 3
**Requirements**: TVAL-02, TVAL-03, TVAL-04, TVAL-06
**Success Criteria** (what must be TRUE):
  1. For EnergyDomain on 3-4 qubit Heisenberg: trajectory-averaged rho matches DM rho within statistical tolerance (trace distance bounded by C/sqrt(N_traj))
  2. For TimeDomain on 3-4 qubit Heisenberg: trajectory-averaged rho matches DM rho within statistical tolerance
  3. For TrotterDomain with_coherent=true on 3-4 qubit Heisenberg: trajectory-averaged rho matches DM rho within statistical tolerance
  4. TrotterDomain with_coherent=true reaches trace distance to Gibbs state of 1e-6 or less, confirming coherent term correctness in trajectory mode
**Plans**: TBD

Plans:
- [ ] 04-01: Trajectory vs DM cross-validation for Energy, Time, Trotter domains
- [ ] 04-02: Coherent term Gibbs convergence test

### Phase 5: Statistical Validation and Regression
**Goal**: Trajectory convergence properties are verified and known-good numerical results are frozen for regression testing
**Depends on**: Phase 4
**Requirements**: TVAL-05, TINF-02
**Success Criteria** (what must be TRUE):
  1. Trajectory error (trace distance to DM) decreases as 1/sqrt(N_traj) when doubling trajectory count, verified with a geometric progression of N_traj values
  2. Frozen reference data files exist for known-good DM and trajectory results, and regression tests compare current output against these references within tight tolerances
**Plans**: TBD

Plans:
- [ ] 05-01: Statistical convergence test (1/sqrt(N) scaling)
- [ ] 05-02: Regression test framework with frozen reference data

## Progress

**Execution Order:**
Phases execute in numeric order: 1 -> 2 -> 3 -> 4 -> 5
(Note: Phase 2 and Phase 3 share only a Phase 1 dependency and could overlap, but Phase 4 requires both.)

| Phase | Plans Complete | Status | Completed |
|-------|----------------|--------|-----------|
| 1. Foundation and Compilation | 1/1 | ✓ Complete | 2026-02-13 |
| 2. Trajectory Bug Fixes | 2/2 | ✓ Complete | 2026-02-14 |
| 3. DM Reference Test Suite | 0/3 | Not started | - |
| 4. Trajectory Cross-Validation | 0/2 | Not started | - |
| 5. Statistical Validation and Regression | 0/2 | Not started | - |
