# Requirements: QuantumFurnace.jl v2.2

**Defined:** 2026-03-04
**Core Value:** Correct and efficient classical simulation of Lindbladian-based quantum Gibbs samplers

## v2.2 Requirements

Requirements for Hamiltonian simulation time counting. Each maps to roadmap phases.

### Grid Infrastructure

- [ ] **GRID-01**: QPE grid generation from r estimating qubits producing full 2^r time/energy labels (not truncated)
- [ ] **GRID-02**: Grid parameter reporting: N, w0, t0, t_max, energy range stored in result struct

### OFT Time Counting

- [ ] **OFT-01**: OFT Hamiltonian simulation time computed as sum of |t_k| over full QPE grid (unweighted by filter function)
- [ ] **OFT-02**: OFT time counting supports Gaussian transition weight function
- [ ] **OFT-03**: OFT time counting supports Metropolis transition weight function (a=0)
- [ ] **OFT-04**: OFT time counting supports smooth Metropolis transition weight function (a>0, b>0)

### B Coherent Term

- [ ] **BTERM-01**: B coherent term Ham sim time computed as double sum over b_minus (outer) and b_plus (inner) truncated dictionaries
- [ ] **BTERM-02**: B term returns 0.0 for GNS construction (no coherent term)
- [ ] **BTERM-03**: B term supports all 3 transition weight functions via b_plus variant dispatch

### Cost Assembly

- [ ] **COST-01**: Per-step cost computed as 2 × OFT_time + B_time (factor of 2 for dissipative U + controlled-U† only)
- [ ] **COST-02**: Total step count computed as ceil(mixing_time / delta)
- [ ] **COST-03**: Total Ham sim time computed as n_steps × per_step_cost
- [ ] **COST-04**: Immutable result struct with per-component breakdown (oft_time, b_time, per_step_time, n_steps, total_time, grid info, parameters)

### API

- [ ] **API-01**: Primary function accepts HamHam + scalar parameters (r, sigma, beta, delta, mixing_time)
- [ ] **API-02**: HamHam convenience overload extracts rescaling_factor and n_qubits automatically
- [ ] **API-03**: MixingTimeEstimate convenience overload extracts mixing_time from estimate
- [ ] **API-04**: Mixing time accepted as explicit user input (not derived from simulation state)

### Validation

- [ ] **VAL-01**: Closed-form sanity check: unweighted QPE grid sum equals t0 × N²/4
- [ ] **VAL-02**: Validation test cases with r=12 on small systems (n=3,4)
- [ ] **VAL-03**: B_time = 0.0 verified for GNS construction
- [ ] **VAL-04**: per_step_time = 2 × oft_time + b_time verified (not 2 × (oft_time + b_time))

## Future Requirements

### Gate Complexity (v2.3+)

- **GATE-01**: Trotter step counting for quantum algorithm's internal Hamiltonian simulation
- **GATE-02**: Gate-level circuit compilation cost estimation
- **GATE-03**: QSP cost for B coherent term application

### Scaling Studies (separate milestone)

- **SCALE-01**: Mixing time scaling with (n, beta) for Heisenberg chains via bi-exp DM thermalization
- **SCALE-02**: Ham sim time scaling with (n, beta) for given Hamiltonians
- **SCALE-03**: Parameter sweep scripts and plotting infrastructure

## Out of Scope

| Feature | Reason |
|---------|--------|
| Domain dispatch (Energy/Trotter) in cost counting | Quantum algorithm cost is domain-independent; domain is a classical simulation distinction |
| Config-based API | Config carries simulation infrastructure irrelevant to quantum cost counting |
| Running actual simulations for cost estimation | Cost counting is purely analytical arithmetic on parameters |
| Modifying HamHam or Config structs | Breaks BSON serialization; new standalone struct instead |
| Plotting or visualization | Presentation is a separate concern; return data struct |
| Asymptotic complexity analysis | This counter gives exact numerical values, not big-O bounds |
| Gate-level compilation | Depends on Hamiltonian simulation method choice (Trotter order, QSP, etc.) |

## Traceability

| Requirement | Phase | Status |
|-------------|-------|--------|
| GRID-01 | Phase 44 | Pending |
| GRID-02 | Phase 44 | Pending |
| OFT-01 | Phase 45 | Pending |
| OFT-02 | Phase 45 | Pending |
| OFT-03 | Phase 45 | Pending |
| OFT-04 | Phase 45 | Pending |
| BTERM-01 | Phase 46 | Pending |
| BTERM-02 | Phase 46 | Pending |
| BTERM-03 | Phase 46 | Pending |
| COST-01 | Phase 47 | Pending |
| COST-02 | Phase 47 | Pending |
| COST-03 | Phase 47 | Pending |
| COST-04 | Phase 44 | Pending |
| API-01 | Phase 47 | Pending |
| API-02 | Phase 47 | Pending |
| API-03 | Phase 47 | Pending |
| API-04 | Phase 47 | Pending |
| VAL-01 | Phase 45 | Pending |
| VAL-02 | Phase 47 | Pending |
| VAL-03 | Phase 46 | Pending |
| VAL-04 | Phase 47 | Pending |

**Coverage:**
- v2.2 requirements: 21 total
- Mapped to phases: 21
- Unmapped: 0

---
*Requirements defined: 2026-03-04*
*Last updated: 2026-03-04 after roadmap creation*
