# Requirements: QuantumFurnace.jl

**Defined:** 2026-02-15
**Core Value:** Correct and efficient classical simulation of Lindbladian-based quantum Gibbs samplers

## v1.2 Requirements

Requirements for multi-threaded trajectory engine and KMS-vs-GNS comparison experiments. Each maps to roadmap phases.

### Threading

- [ ] **THRD-01**: TrajectoryWorkspace is separated from TrajectoryFramework and passed explicitly to step_along_trajectory!
- [ ] **THRD-02**: Multi-threaded trajectory sampling runs N trajectories across threads with per-thread workspace, seeding task-local RNG per trajectory for reproducibility
- [ ] **THRD-03**: BLAS thread count is set to 1 during multi-threaded trajectory execution to avoid oversubscription
- [ ] **THRD-04**: Trajectory results are deterministic given a master seed and same thread count

### GNS Path

- [ ] **GNS-01**: GNS trajectory simulation works via ThermalizeConfigGNS dispatch (no B term, approximate detailed balance)
- [ ] **GNS-02**: GNS trajectories produce valid density matrices that converge toward the GNS fixed point

### Convergence

- [ ] **CONV-01**: Trace distance to Gibbs state tracked at batch checkpoints during trajectory sampling
- [ ] **CONV-02**: Per-observable convergence tracked for nearest-neighbor correlations <Z_iZ_{i+1}>
- [ ] **CONV-03**: Per-observable convergence tracked for energy <H>
- [ ] **CONV-04**: Adaptive sampling runs trajectory batches until convergence criterion met (relative change <1% for 3 consecutive batches)
- [ ] **CONV-05**: Hard cap on maximum trajectories prevents infinite adaptive loops

### Data

- [ ] **DATA-01**: Experiment results saved to BSON with config, convergence curves, observables, and final density matrix
- [ ] **DATA-02**: Experiment results loadable and reproducible from saved BSON files

### Experiments

- [ ] **EXPT-01**: KMS-vs-GNS experiment driver runs matched experiments (same Hamiltonian, beta, delta) with KMS and GNS(sigma=1/beta)
- [ ] **EXPT-02**: GNS experiments also run at sigma=0.5/beta for cost-accuracy comparison
- [ ] **EXPT-03**: Experiments sweep across n=4,6,8 and beta=5,10,20 using TrotterDomain
- [ ] **EXPT-04**: Experiment results show KMS converges closer to Gibbs state than GNS at sigma=1/beta

## Future Requirements

### Resource Estimation

- **RSRC-01**: Hamiltonian simulation time counter tracks total simulation cost per Gibbs sample
- **RSRC-02**: Gate complexity counter for Trotter-based circuit implementations
- **RSRC-03**: Mixing time estimation from trajectory convergence (for n>8 where full Lindbladian infeasible)

### Additional Hamiltonians

- **HAM-01**: 1D Ising model Hamiltonian generation
- **HAM-02**: 2D Heisenberg Hamiltonian generation (lattice graph support)
- **HAM-03**: General k-local Hamiltonian construction on arbitrary graphs

### Documentation & Plotting

- **DOC-01**: API docs via Documenter.jl, theory tutorials via Literate.jl
- **PLOT-01**: Paper-ready convergence plots, sim time plots, gate complexity, mixing time

## Out of Scope

| Feature | Reason |
|---------|--------|
| GPU acceleration | Not needed for current system sizes (n<=12); dim=256 too small for GPU advantage |
| MPI distributed trajectories | Single-node multi-core sufficient; cluster nodes have enough RAM |
| Float32 trajectory testing | Type params enable it but precision insufficient for KMS convergence targets (~1e-6) |
| Continuous-time adaptive timestep | QuantumFurnace uses discrete-time CPTP maps per Chen's formulation |
| Bootstrap confidence intervals | Deferred to later milestone for paper polish |
| GNS sigma sweep beyond two points | Two comparison points (1/beta, 0.5/beta) sufficient for this milestone |
| Local magnetization <Z_i> tracking | Nearest-neighbor correlations and energy sufficient for paper |
| Convergence plotting | Data architecture only; plotting deferred to later milestone |
| Qiskit circuit generation | Future milestone for resource estimation |

## Traceability

| Requirement | Phase | Status |
|-------------|-------|--------|
| THRD-01 | — | Pending |
| THRD-02 | — | Pending |
| THRD-03 | — | Pending |
| THRD-04 | — | Pending |
| GNS-01 | — | Pending |
| GNS-02 | — | Pending |
| CONV-01 | — | Pending |
| CONV-02 | — | Pending |
| CONV-03 | — | Pending |
| CONV-04 | — | Pending |
| CONV-05 | — | Pending |
| DATA-01 | — | Pending |
| DATA-02 | — | Pending |
| EXPT-01 | — | Pending |
| EXPT-02 | — | Pending |
| EXPT-03 | — | Pending |
| EXPT-04 | — | Pending |

**Coverage:**
- v1.2 requirements: 17 total
- Mapped to phases: 0
- Unmapped: 17 (awaiting roadmap)

---
*Requirements defined: 2026-02-15*
*Last updated: 2026-02-15 after initial definition*
