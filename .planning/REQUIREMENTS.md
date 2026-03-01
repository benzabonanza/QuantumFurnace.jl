# Requirements: QuantumFurnace.jl

**Defined:** 2026-03-01
**Core Value:** Correct and efficient classical simulation of Lindbladian-based quantum Gibbs samplers

## v2.1 Requirements

Requirements for v2.1 Speedup & Mixing Time. Each maps to roadmap phases.

### Per-Jump Precomputation

- [ ] **PRECOMP-01**: R^a (rate-weighted sum over frequencies of L_w'*L_w) precomputed once per jump at simulation start for Energy, Time, and Trotter domains
- [ ] **PRECOMP-02**: K0^a and U_residual^a derived from precomputed R^a via _build_cptp_channel, stored per jump
- [ ] **PRECOMP-03**: Hot loop uses precomputed K0^a, U_residual^a to apply CPTP channel — no per-step eigendecomposition
- [ ] **PRECOMP-04**: BohrDomain receives general speedups and threading but no per-Bohr-frequency vector precomputation (Bohr frequency count grows too fast for large systems)

### Multi-Threading

- [ ] **THREAD-01**: Dissipative sandwich ω-loop (rho_jump accumulation) runs multi-threaded with per-task accumulators for Energy, Time, and Trotter domains
- [ ] **THREAD-02**: R^a precomputation ω-loop runs multi-threaded with per-task accumulators
- [ ] **THREAD-03**: Multi-threaded BLAS enabled during DM thermalization channel application (no trajectory-level parallelism competing)
- [ ] **THREAD-04**: BohrDomain multi-threaded where beneficial (Bohr bucket iteration, sandwich accumulation)
- [ ] **THREAD-05**: BLAS thread control follows existing try/finally save/restore pattern from trajectory engine

### Save Every

- [ ] **SAVE-01**: save_every keyword in run_thermalize controls how often trace distance to Gibbs state is computed and stored
- [ ] **SAVE-02**: Default save_every=1 preserves backward compatibility (trace distance computed every step)

### Mixing Time Estimation

- [ ] **MIX-01**: Exponential fit A*exp(-gap*t)+C applied to trace distance convergence curve after skip_initial burn-in fraction
- [ ] **MIX-02**: extrapolate=true keyword causes simulation to stop when fit is reliable, returns estimated mixing time to target epsilon without full convergence
- [ ] **MIX-03**: extrapolate=false (default) runs full convergence simulation, reports actual number of steps/time to reach target trace distance
- [ ] **MIX-04**: skip_initial keyword (default 0.2) controls fraction of trace distance data excluded from exponential fit as burn-in
- [ ] **MIX-05**: ThermalizeResults extended with mixing_time, fitted_gap, and fit quality metrics (R-squared, confidence interval)
- [ ] **MIX-06**: Final density matrix for extrapolation case is the last actual simulated DM before extrapolation stopped
- [ ] **MIX-07**: Quality gates on extrapolation: warn when R^2 < 0.95, warn when offset C is large relative to target epsilon
- [ ] **MIX-08**: LsqFit.jl re-added as dependency and fitting code promoted from src/staging/ to active source

## Future Requirements

### Mixing Time Refinements

- **MIX-F01**: Automatic burn-in detection from trace distance curve shape (replace fixed skip_initial)
- **MIX-F02**: Effective rate diagnostic λ_eff(t) for model-free validation of fitting window
- **MIX-F03**: Two-exponential fitting for transient contamination absorption
- **MIX-F04**: Bootstrap uncertainty quantification on fitted gap

### Performance Refinements

- **PERF-F01**: BohrDomain per-Bohr-frequency precomputation (if needed for specific use cases)
- **PERF-F02**: Adaptive threading strategy (auto-select ω-parallel vs BLAS-parallel based on system size)

## Out of Scope

| Feature | Reason |
|---------|--------|
| Richardson extrapolation in delta | Deferred from v1.4; not needed for DM-based mixing time (no Trotter bias in DM trace distance) |
| Dashboard/plotting for diagnostics | Paper-ready plotting is a separate future milestone |
| Trajectory-based mixing time estimation | DM trace distance is noiseless and sufficient; trajectory gap estimation remains in staging |
| BohrDomain per-frequency precomputation | Bohr frequency count grows too fast for large systems; would be anti-feature for target sizes |

## Traceability

| Requirement | Phase | Status |
|-------------|-------|--------|
| PRECOMP-01 | Phase 39 | Pending |
| PRECOMP-02 | Phase 39 | Pending |
| PRECOMP-03 | Phase 39 | Pending |
| PRECOMP-04 | Phase 39 | Pending |
| THREAD-01 | Phase 41 | Pending |
| THREAD-02 | Phase 41 | Pending |
| THREAD-03 | Phase 41 | Pending |
| THREAD-04 | Phase 41 | Pending |
| THREAD-05 | Phase 41 | Pending |
| SAVE-01 | Phase 40 | Pending |
| SAVE-02 | Phase 40 | Pending |
| MIX-01 | Phase 42 | Pending |
| MIX-02 | Phase 42 | Pending |
| MIX-03 | Phase 42 | Pending |
| MIX-04 | Phase 42 | Pending |
| MIX-05 | Phase 42 | Pending |
| MIX-06 | Phase 42 | Pending |
| MIX-07 | Phase 42 | Pending |
| MIX-08 | Phase 42 | Pending |

**Coverage:**
- v2.1 requirements: 19 total
- Mapped to phases: 19
- Unmapped: 0

---
*Requirements defined: 2026-03-01*
*Last updated: 2026-03-01 after roadmap creation*
