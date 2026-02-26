# Requirements: QuantumFurnace.jl

**Defined:** 2026-02-25
**Core Value:** Correct and efficient classical simulation of Lindbladian-based quantum Gibbs samplers

## v2.0 Requirements

Requirements for the v2.0 Restructure milestone. Each maps to roadmap phases.

### Type Foundation

- [ ] **TYPE-01**: Define `AbstractSimulation` hierarchy with `Lindbladian`, `Thermalize`, `KrylovSpectrum`, `Trajectory` singleton types
- [ ] **TYPE-02**: Define `AbstractConstruction` hierarchy with `KMS`, `GNS` singleton types (+ `DLL` placeholder for future Ding et al.)
- [ ] **TYPE-03**: Define unified `Config{S,D,C,T}` struct replacing `LiouvConfig`, `LiouvConfigGNS`, `ThermalizeConfig`, `ThermalizeConfigGNS`
- [ ] **TYPE-04**: Derive `with_coherent` from construction type (`KMS` -> true, `GNS` -> false) instead of storing as field
- [ ] **TYPE-05**: Provide backward-compatible type aliases during migration (`const LiouvConfig{D,T} = Config{Lindbladian,D,KMS,T}` etc.)
- [ ] **TYPE-06**: Migrate all dispatch sites (run_*, _precompute_data, _jump_contribution!, krylov, trajectories) to use new `Config{S,D,C,T}`

### Code Deduplication

- [ ] **DEDUP-01**: Extract `domain_prefactor()` function replacing 16 identical formula copies across jump_workers.jl, trajectories.jl, krylov_workspace.jl, krylov_matvec.jl, krylov_eigsolve.jl
- [ ] **DEDUP-02**: ~~Extract `foreach_frequency()` iterator replacing 16 hermitian half-grid branching patterns across all simulation paths~~ — **DEFERRED** (user decision: keep explicit for-loop patterns as-is; see Phase 34 CONTEXT.md)
- [ ] **DEDUP-03**: Unify `oft!()` with best logic from `_krylov_oft!` (keep `time_oft!`/`trotter_oft!` as test/debug utilities with clear marking)

### Workspace & API

- [ ] **WORK-01**: Consolidate `KrylovWorkspace` + `KrausScratch` + `LindbladianWorkspace` into unified workspace struct (`TrajectoryWorkspace` stays separate)
- [ ] **WORK-02**: Unify R/K0/U_residual computation: extract shared CPTP channel formulas into helper functions; per-jump (R^a, K0^a, U_residual^a) for DM/Trajectory path, summed (R_total, K0_total, U_residual_total) for Krylov path
- [ ] **WORK-03**: Define 4 clean `run_*` entry points: `run_lindblad()`, `run_thermalize()`, `run_krylov_spectrum()`, `run_trajectory()`
- [ ] **WORK-04**: Define 4 Result structs: `LindbladResults`, `ThermalizeResults`, `KrylovSpectrumResults`, `TrajectoryResults`
- [ ] **WORK-05**: Add save capability to all Result structs (optional BSON serialization with metadata: git hash, timestamp, config)

### Organization & Cleanup

- [ ] **ORG-01**: Rename src/ files to match PRE/MID/POST architecture grouping (stay flat, no subdirectories)
- [ ] **ORG-02**: Move gap estimation/fitting/convergence code to staging area separated from active source (src/staging/ or similar)
- [ ] **ORG-03**: Remove `@distributed` dead code from furnace.jl and remove `using Distributed` import from module
- [ ] **ORG-04**: Consolidate test helpers: parametrize config factories by system size, eliminate duplicate setup patterns across test files
- [ ] **ORG-05**: Add `@info` printout to all test blocks showing what is being tested and key numerical results
- [ ] **ORG-06**: Review and fix dubious test thresholds/parameters identified during restructure
- [ ] **ORG-07**: Create/update 4 simulation scripts in `simulations/` matching the 4 `run_*` entry points
- [ ] **ORG-08**: Clean module export list matching new structure (organized by simulation type)
- [ ] **ORG-09**: Keep diagnostics module as a separate analysis capability (not folded into Lindblad path)

## Future Requirements

### DLL Construction (Ding et al. 2024)

- **DLL-01**: `DLL` construction type with different filter functions, no energy labels
- **DLL-02**: BohrDomain, TimeDomain, TrotterDomain variants for DLL
- **DLL-03**: KMS-vs-DLL comparison experiments

### Error Estimation

- **ERR-01**: Quadrature error estimation functions in errors.jl
- **ERR-02**: Trotter error estimation functions in errors.jl

### Analysis Tools

- **ANAL-01**: Quantum discriminant / anti-Hermitian norm analysis
- **ANAL-02**: Kossakowski matrix analysis and visualization
- **ANAL-03**: Lieb-Robinson bound support computation

### Infrastructure

- **INFRA-01**: Log-Sobolev constant computation (rewrite with apply_lindbladian!)
- **INFRA-02**: Hamiltonian simulation time accumulator
- **INFRA-03**: Qiskit circuit generation for gate estimation

## Out of Scope

| Feature | Reason |
|---------|--------|
| DLL construction implementation | Just the type placeholder; full implementation is a separate milestone |
| Sandwich function consolidation | 4 functions, only 2 unique, but not selected for this milestone |
| GPU acceleration | Not needed for current system sizes |
| Float32 path testing | Type params enable it but testing deferred |
| log_sobolev.jl rewrite | Keep as-is, future project |
| errors.jl population | Keep stub, populate in future error estimation milestone |
| Subdirectories in src/ | Decided to stay flat with renamed files |
| SharedArrays removal | Keep for potential future multi-process NUFFT |
| OFTCaches removal | Keep for testing reference |

## Traceability

| Requirement | Phase | Status |
|-------------|-------|--------|
| TYPE-01 | Phase 33 | Pending |
| TYPE-02 | Phase 33 | Pending |
| TYPE-03 | Phase 33 | Pending |
| TYPE-04 | Phase 33 | Pending |
| TYPE-05 | Phase 33 | Pending |
| TYPE-06 | Phase 33 | Pending |
| DEDUP-01 | Phase 34 | Pending |
| DEDUP-02 | ~~Phase 34~~ | Deferred |
| DEDUP-03 | Phase 34 | Pending |
| WORK-01 | Phase 35 | Pending |
| WORK-02 | Phase 35 | Pending |
| WORK-03 | Phase 36 | Pending |
| WORK-04 | Phase 36 | Pending |
| WORK-05 | Phase 36 | Pending |
| ORG-01 | Phase 37 | Pending |
| ORG-02 | Phase 37 | Pending |
| ORG-03 | Phase 37 | Pending |
| ORG-04 | Phase 38 | Pending |
| ORG-05 | Phase 38 | Pending |
| ORG-06 | Phase 38 | Pending |
| ORG-07 | Phase 37 | Pending |
| ORG-08 | Phase 37 | Pending |
| ORG-09 | Phase 37 | Pending |

**Coverage:**
- v2.0 requirements: 23 total
- Mapped to phases: 23
- Unmapped: 0

---
*Requirements defined: 2026-02-25*
*Last updated: 2026-02-25 after roadmap creation*
