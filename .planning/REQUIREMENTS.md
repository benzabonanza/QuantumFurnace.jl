# Requirements: QuantumFurnace.jl v1.0 Trajectories

**Defined:** 2026-02-13
**Core Value:** Correct and efficient classical simulation of Lindbladian-based quantum Gibbs samplers, validated against Chen et al. (2023, 2025)

## v1 Requirements

Requirements for milestone v1.0 Trajectories. Each maps to roadmap phases.

### Trajectory Bug Fixes

- [ ] **TFIX-01**: Fix `build_trajectoryframework` compilation bugs (undefined `trotter` variable, uninitialized `B_total`)
- [ ] **TFIX-02**: Fix coherent unitary U_B ordering in `step_along_trajectory!` to match DM code
- [ ] **TFIX-03**: Add probability normalization assertion (p_nojump + p_res + p_jump_total ≈ 1.0)
- [ ] **TFIX-04**: Add S matrix PSD guard before Cholesky factorization
- [ ] **TFIX-05**: Ensure trajectory jump sampling faithfully implements Chen's weak measurement scheme, matching how the DM simulator applies the CPTP channel — verify by comparing trajectory channel structure against `jump_contribution!` in DM code and the paper's algorithm

### Trajectory Validation

- [ ] **TVAL-01**: CPTP verification test confirms K0†K0 + delta*R + U_res†U_res = I to machine precision
- [ ] **TVAL-02**: Trajectory-averaged rho matches DM evolution for EnergyDomain (3-4 qubit Heisenberg)
- [ ] **TVAL-03**: Trajectory-averaged rho matches DM evolution for TimeDomain (3-4 qubit Heisenberg)
- [ ] **TVAL-04**: Trajectory-averaged rho matches DM evolution for TrotterDomain with_coherent=true (3-4 qubit Heisenberg)
- [ ] **TVAL-05**: Statistical 1/sqrt(N) convergence test shows trajectory error decreases as expected with N_traj
- [ ] **TVAL-06**: Coherent term correctness: with_coherent=true TrotterDomain reaches ≤ 1e-6 distance to Gibbs state

### DM Correctness Tests

- [ ] **DMTST-01**: BohrDomain with coherent term B has exact detailed balance (Gibbs state is fixed point to machine precision)
- [ ] **DMTST-02**: Domain error hierarchy verified: dist_bohr ≤ dist_energy ≤ dist_time ≤ dist_trotter
- [ ] **DMTST-03**: Single-step DM error scales as delta^2 (empirical verification of Chen Theorem III.1)
- [ ] **DMTST-04**: Multi-step DM error accumulates as O(delta) over full evolution
- [ ] **DMTST-05**: Coherent term B consistency across domains: B_bohr vs B_time match up to time quadrature errors; B_trotter matches with additional Trotter errors
- [ ] **DMTST-06**: OFT consistency: oft!() and time_oft!() match up to time quadrature errors; trotter_oft!() matches with additional Trotter errors

### Test Infrastructure

- [ ] **TINF-01**: Test helpers module with shared fixtures (make_test_system, configs, tiered tolerance constants matching error hierarchy)
- [ ] **TINF-02**: Regression tests with frozen reference data for known-good numerical results
- [ ] **TINF-03**: Aqua.jl package quality checks (ambiguities, stale deps, method piracy)
- [ ] **TINF-04**: Project.toml cleanup (move dev deps to [extras], add test-only deps: StableRNGs, HypothesisTests, StatsBase, Aqua)

## v2 Requirements

Deferred to future milestones. Tracked but not in current roadmap.

### Trajectory Enhancements

- **TENH-01**: Multi-threaded trajectory sampling with shared precomputed data
- **TENH-02**: Per-observable trajectory convergence (<Z_i> from trajectories matches DM)
- **TENH-03**: Jump statistics histogram (empirical rates match theoretical predictions)
- **TENH-04**: Confidence interval reporting (bootstrap CIs on trajectory-vs-DM trace distance)

### Additional Validation

- **AVAL-01**: Diamond norm channel comparison between DM and trajectory channels
- **AVAL-02**: Non-Hermitian jump operator support in trajectory tests
- **AVAL-03**: R matrix cross-validation (trajectory R matches Liouvillian R)
- **AVAL-04**: Normalization drift monitoring (pre-normalization ||psi||^2 per step)

## Out of Scope

Explicitly excluded. Documented to prevent scope creep.

| Feature | Reason |
|---------|--------|
| Continuous-time MCWF (adaptive timestep) | QuantumFurnace uses discrete-time CPTP maps per Chen's weak measurement scheme; adaptive timestep requires fundamentally different formulation |
| Bohr domain trajectories | Diagonalizing Kossakowski matrix prohibitive for >3 qubits; validate Bohr via DM only |
| Large system benchmarks (>6 qubits) | Validation requires DM comparison which is the bottleneck; correctness before scale |
| Quantum state tomography from trajectories | Direct psi*psi^dag averaging is the correct and cheaper approach |
| Multi-threaded trajectory validation | Correctness first, performance deferred to v2 |

## Traceability

Which phases cover which requirements. Updated during roadmap creation.

| Requirement | Phase | Status |
|-------------|-------|--------|
| TFIX-01 | Phase 1 | Pending |
| TFIX-02 | Phase 2 | Pending |
| TFIX-03 | Phase 2 | Pending |
| TFIX-04 | Phase 2 | Pending |
| TFIX-05 | Phase 2 | Pending |
| TVAL-01 | Phase 2 | Pending |
| TVAL-02 | Phase 4 | Pending |
| TVAL-03 | Phase 4 | Pending |
| TVAL-04 | Phase 4 | Pending |
| TVAL-05 | Phase 5 | Pending |
| TVAL-06 | Phase 4 | Pending |
| DMTST-01 | Phase 3 | Pending |
| DMTST-02 | Phase 3 | Pending |
| DMTST-03 | Phase 3 | Pending |
| DMTST-04 | Phase 3 | Pending |
| DMTST-05 | Phase 3 | Pending |
| DMTST-06 | Phase 3 | Pending |
| TINF-01 | Phase 1 | Pending |
| TINF-02 | Phase 5 | Pending |
| TINF-03 | Phase 3 | Pending |
| TINF-04 | Phase 1 | Pending |

**Coverage:**
- v1 requirements: 21 total
- Mapped to phases: 21
- Unmapped: 0

---
*Requirements defined: 2026-02-13*
*Last updated: 2026-02-13 after roadmap creation*
