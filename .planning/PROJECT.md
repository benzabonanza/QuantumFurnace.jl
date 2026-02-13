# QuantumFurnace.jl

## What This Is

A Julia package for simulating quantum Gibbs sampling via Lindbladian evolution. It implements multiple Lindbladian constructions (GNS and KMS detailed balance) across a hierarchy of approximation domains (Bohr, Energy, Time, Trotter), enabling both density matrix and trajectory-based simulation of quantum thermalization. The package targets researchers and students working on quantum Gibbs sampling, providing both fast numerical simulations and pedagogic documentation grounded in the theoretical literature.

## Core Value

Correct and efficient classical simulation of Lindbladian-based quantum Gibbs samplers, validated against the mathematical constructions in Chen et al. (2023, 2025) and Ding et al. (2024), so that researchers can trust the numerical results for publication and learning.

## Requirements

### Validated

- Density matrix Lindbladian construction for GNS detailed balance (Chen 2023 approximate construction) -- existing
- Density matrix Lindbladian construction for KMS detailed balance (Chen 2025 exact construction) -- existing
- Four approximation domains: Bohr, Energy, Time, Trotter -- existing
- Gaussian and Metropolis transition weight functions -- existing
- Coherent term (B operator) computation for exact KMS detailed balance -- existing
- `run_lindbladian()` computing full Liouvillian spectrum and steady state -- existing
- `run_thermalize()` step-by-step thermalization with trace distance tracking -- existing
- 1D Heisenberg Hamiltonian with optional external field -- existing
- Hamiltonian eigendecomposition, Bohr frequency computation, Gibbs state construction -- existing
- NUFFT-based precomputation for Time and Trotter domains -- existing
- Workspace/precomputation caching pattern for performance -- existing
- BSON serialization for Hamiltonian objects -- existing

### Active

- [ ] Trajectory simulation: debug and validate `TrajectoryFramework` against density matrix results
- [ ] 1D Ising model Hamiltonian generation
- [ ] 2D Heisenberg Hamiltonian generation (lattice graph support)
- [ ] General k-local Hamiltonian construction on arbitrary graphs
- [ ] Comprehensive test suite: unit tests for each domain, cross-validation between DM and trajectory methods
- [ ] Error/convergence analysis functions: trace distance to Gibbs, fidelity, relative entropy vs. time/steps
- [ ] Qiskit circuit generation for resource estimation (gate count, circuit depth)
- [ ] Hamiltonian simulation time counter (total Hamiltonian simulation cost per Gibbs sample)
- [ ] Gate complexity counter for Trotter-based circuit implementations
- [ ] Multi-threaded trajectory sampling with shared precomputed data
- [ ] Documentation: API docs via Documenter.jl, theory tutorials via Literate.jl
- [ ] Paper-ready numerical results: convergence plots, mixing time scaling, domain comparison
- [ ] Ding et al. (2024) KMS Lindbladian construction with discrete jump operators (future addition)

### Out of Scope

- GPU acceleration -- not needed for current system sizes (~12 qubits); revisit if scaling beyond 14+ qubits
- Real quantum hardware execution -- circuits are for resource estimation only, not viable for current hardware
- Quantum circuit simulation via Qiskit -- simulation done by QuantumFurnace's own simulators
- MPI distributed computing for single simulations -- single-node multi-core sufficient; cluster used for parameter sweeps via independent Julia processes
- Classical MCMC / classical Gibbs sampling algorithms -- focus is quantum Lindbladian methods only
- Non-Markovian dynamics -- all constructions assume Markovian (Lindblad) framework

## Context

**Theoretical Foundation:**
The package implements quantum Gibbs samplers from three key papers:
1. **Chen, Kastoryano, Brandao, Gilyen (2023)** "Quantum Thermal State Preparation" -- introduces smoothed Davies generators with approximate GNS detailed balance via Gaussian-weighted operator Fourier transforms. The Lindbladian's fixed point is approximately the Gibbs state. Energy uncertainty ~O(sigma_E) from finite-time Fourier transform.
2. **Chen, Kastoryano, Gilyen (2025)** "An efficient and exact noncommutative quantum Gibbs sampler" -- constructs the first exactly KMS-detailed-balanced Lindbladian with Gaussian transition weight gamma(omega) and a coherent term B that exactly cancels errors. Fixed point is exactly the Gibbs state. Cost: O_tilde(beta) Hamiltonian simulation time per unit Lindbladian evolution.
3. **Ding, Li, Lin (2024)** "Efficient quantum Gibbs samplers with KMS detailed balance condition" -- generalizes KMS construction to discrete jump operators (as few as one), simplifying implementation. Uses Fagnola-Umanita structural characterization. Planned for future integration.

**Existing Codebase State:**
- Core density matrix simulation (`run_lindbladian`, `run_thermalize`) works and produces results
- Trajectory simulation recently refactored (`TrajectoryFramework`) but not yet validated -- commits indicate "updates without test runs" and "definitely not yet done"
- Coherent term computation has a TODO about sigma parameters that may affect correctness
- ~200 lines of commented-out old code in jump_workers.jl and trajectories.jl
- Hardcoded tolerances (1e-12) scattered without centralized configuration
- Test infrastructure exists but focuses on timing/benchmarks rather than correctness validation

**Target System Sizes:**
- Primary: up to ~12 qubits (density matrix 4096x4096, Liouvillian 16M entries)
- Memory driver: precomputed Kraus/NUFFT prefactor matrices, not the density matrix itself
- Trajectories need thousands of samples for statistics -- parallelism over samples is key
- Single-node execution on cluster (up to 512 GB RAM), multi-core parallelism

**Paper Goals:**
Results needed for publication: convergence curves (trace distance vs. steps), mixing time scaling with system size and temperature, comparison across approximation domains (Bohr vs Energy vs Time vs Trotter), trajectory vs density matrix agreement.

## Constraints

- **Language**: Julia -- non-negotiable, leverages Julia's multiple dispatch for domain hierarchy and performance
- **Correctness**: All Lindbladian constructions must be mathematically faithful to the source papers; detailed balance properties must hold to machine precision (for exact KMS) or documented approximation bounds (for GNS)
- **System size**: Practical limit ~12-14 qubits due to exponential Hilbert space scaling (density matrix 2^n x 2^n)
- **Dependencies**: Minimize external dependencies; rely on Julia stdlib + established numerical packages (Arpack, FINUFFT)
- **Dual purpose**: Code must be both performant (for research results) and readable (for community pedagogy) -- well-documented functions with mathematical cross-references to papers

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| Four-domain approximation hierarchy (Bohr/Energy/Time/Trotter) | Maps directly to the theoretical approximation levels in Chen 2023/2025; enables comparing accuracy vs. cost | -- Pending evaluation |
| KMS detailed balance as primary construction | Exact stationarity (Chen 2025) eliminates accumulated trajectory error vs approximate GNS | -- Pending validation |
| Density matrix + trajectory dual simulation | DM gives exact dynamics for small systems; trajectories scale better and model the actual quantum algorithm | -- Pending cross-validation |
| Qiskit for circuit generation (Python interop) | Qiskit is the standard for quantum circuit representation; Julia quantum circuit ecosystem less mature | -- Pending |
| Single-node multi-core for trajectories | Shared memory for precomputed data avoids serialization overhead; cluster nodes have enough RAM | -- Pending benchmarks |

---
*Last updated: 2026-02-13 after initialization*
