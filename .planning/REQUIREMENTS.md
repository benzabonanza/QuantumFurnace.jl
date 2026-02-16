# Requirements: QuantumFurnace.jl

**Defined:** 2026-02-16
**Core Value:** Correct and efficient classical simulation of Lindbladian-based quantum Gibbs samplers

## v1.3 Requirements

Requirements for milestone v1.3 Mixing Time Estimation. Each maps to roadmap phases.

### Observables

- [ ] **OBS-01**: User can build total magnetization observable (M_z = sum Z_i) in Hamiltonian eigenbasis
- [ ] **OBS-02**: User can build total magnetization observable in Trotter eigenbasis (for TrotterDomain)
- [ ] **OBS-03**: User can build combined gap estimation observables (H, M_z, all ZZ correlations) via single function call

### Trajectory Runner

- [ ] **TRAJ-01**: User can run trajectories measuring time-resolved observables without mid-simulation DM reconstruction
- [ ] **TRAJ-02**: User can optionally reconstruct the averaged DM once at the end of the observable-only run
- [ ] **TRAJ-03**: Observable-only runner supports multi-threaded execution (same pattern as existing run_trajectories)

### Fitting

- [ ] **FIT-01**: User can fit single-exponential decay model (A * exp(-gap * t) + C) to observable time series using LsqFit.jl
- [ ] **FIT-02**: Fitting auto-generates initial guess via log-linear estimate (no manual tuning required for typical data)
- [ ] **FIT-03**: User can specify fitting window (skip_initial fraction) to exclude early multi-exponential transient
- [ ] **FIT-04**: Fit returns quality metrics: R-squared, confidence interval on gap, standard error
- [ ] **FIT-05**: Fit enforces gap > 0 via parameter bounds

### Gap Estimation API

- [ ] **GAP-01**: User can call `estimate_spectral_gap` to run trajectories + fit all observables + select best estimate
- [ ] **GAP-02**: Result struct `SpectralGapResult` contains gap estimate, CI, per-observable gaps, R-squared values, and fit metadata
- [ ] **GAP-03**: Best observable selected by fit quality (R-squared, convergence, gap > 0)

### Cross-Validation

- [ ] **VAL-01**: User can cross-validate fitted gap against exact Liouvillian spectral gap using `abs(real(spectral_gap))`
- [ ] **VAL-02**: Cross-validation warns when imaginary part of exact eigenvalue is significant (|Im/Re| > 0.1)
- [ ] **VAL-03**: Validation script demonstrates gap estimation agreement for n=4 and n=6

## Future Requirements

Deferred to future milestones. Tracked but not in current roadmap.

### Fitting Extensions

- **FIT-EXT-01**: Variance-weighted fitting using per-time-point trajectory variance
- **FIT-EXT-02**: Damped oscillation model (A * exp(-gamma * t) * cos(omega * t + phi) + C) for complex eigenvalues
- **FIT-EXT-03**: Automatic model selection (BIC/AIC) between exponential and damped oscillation
- **FIT-EXT-04**: Multi-exponential fit (A1 * exp(-g1 * t) + A2 * exp(-g2 * t) + C) for resolving multiple modes
- **FIT-EXT-05**: Bootstrap confidence intervals on gap from trajectory batch resampling

### Diagnostics

- **DIAG-01**: Observable overlap diagnostic (compute |tr(O * gap_mode)| for each observable)
- **DIAG-02**: Gap vs beta scaling plots (paper figures)
- **DIAG-03**: Gap vs n scaling plots (test system-size independence prediction)

## Out of Scope

Explicitly excluded. Documented to prevent scope creep.

| Feature | Reason |
|---------|--------|
| Prony's method for multi-exponential extraction | Severely ill-conditioned for noisy trajectory data; single-exponential with window selection is more robust |
| Full Liouvillian construction for n>6 | 65536x65536 at n=8; trajectory-based estimation is the whole point |
| Autocorrelation-based gap estimation | Requires stationarity; our trajectories start far from equilibrium |
| Stochastic trace estimation for Liouvillian gap | Requires Liouvillian matvec primitive we don't have |
| GPU-accelerated fitting | Fitting is tiny computation; trajectory generation is the bottleneck |
| Real-time plotting during fitting | Post-process from saved data; paper plots are separate |

## Traceability

Which phases cover which requirements. Updated during roadmap creation.

| Requirement | Phase | Status |
|-------------|-------|--------|
| OBS-01 | — | Pending |
| OBS-02 | — | Pending |
| OBS-03 | — | Pending |
| TRAJ-01 | — | Pending |
| TRAJ-02 | — | Pending |
| TRAJ-03 | — | Pending |
| FIT-01 | — | Pending |
| FIT-02 | — | Pending |
| FIT-03 | — | Pending |
| FIT-04 | — | Pending |
| FIT-05 | — | Pending |
| GAP-01 | — | Pending |
| GAP-02 | — | Pending |
| GAP-03 | — | Pending |
| VAL-01 | — | Pending |
| VAL-02 | — | Pending |
| VAL-03 | — | Pending |

**Coverage:**
- v1.3 requirements: 17 total
- Mapped to phases: 0
- Unmapped: 17 ⚠️

---
*Requirements defined: 2026-02-16*
*Last updated: 2026-02-16 after initial definition*
