# Feature Research

**Domain:** Quantum Gibbs sampling / open quantum systems simulation (Julia package, paper companion + community tool)
**Researched:** 2026-02-13
**Confidence:** MEDIUM-HIGH (direct codebase analysis + verified external landscape; some circuit generation details LOW)

## Feature Landscape

### Table Stakes (Users Expect These)

Features users assume exist. Missing these = product feels incomplete or unpublishable.

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| **Correctness tests against known Gibbs states** | Every quantum simulation paper must demonstrate convergence to analytically known thermal states. Reviewers will reject without this. | MEDIUM | QuantumFurnace has `distances_to_gibbs` tracking but no automated test suite asserting convergence. Need tests for 1-4 qubit systems with known analytics (single qubit thermal state, 2-qubit Heisenberg, etc.). |
| **Multiple error/distance metrics** | Trace distance alone is insufficient. Community expects fidelity, relative entropy, and operator norm comparisons. QuTiP provides trace_distance, fidelity, entropy, Bures distance. | LOW | `qi_tools.jl` has `trace_distance_h`. Add fidelity (F = tr(sqrt(sqrt(rho)*sigma*sqrt(rho)))^2), relative entropy (S(rho||sigma)), and diamond norm for channel comparison. |
| **Documented API with docstrings** | QuantumOptics.jl and QuantumToolbox.jl both ship with full Documenter.jl sites. Papers link to docs. Without them the package is unusable by anyone but the author. | MEDIUM | Documenter.jl infrastructure exists (`docs/make.jl`, literate tutorials). Many functions lack docstrings. Need systematic pass to document every exported function. |
| **Reproducible examples / tutorials** | Every competing package (QuantumOptics.jl, QuantumToolbox.jl, QuTiP) has step-by-step tutorials. A paper companion must have notebooks or literate scripts that reproduce paper figures. | MEDIUM | Literate tutorials exist in `docs/src/literate/` but several are stubs or reference outdated APIs. Need working end-to-end examples for each domain (Bohr, Energy, Time, Trotter). |
| **Proper Julia package structure** | `] add QuantumFurnace` must work. Tests must run via `] test`. CI must pass. This is the bare minimum for a registered Julia package. | LOW | `Project.toml` exists. No CI configuration found (no `.github/workflows/`). Test files exist but use `includet` (Revise-style) rather than standard `Test` framework. |
| **Hamiltonian library beyond Heisenberg** | Users expect to test with standard models: Ising, XY, XXZ, random local Hamiltonians. A single Hamiltonian makes the package feel like a homework assignment, not a research tool. | LOW-MEDIUM | Currently only Heisenberg with Z-disorder. The `HamHam` constructor is general enough -- need convenience constructors for Ising (transverse field), XY, XXZ, and a random-local-Hamiltonian generator. |
| **Density matrix evolution (Lindbladian master equation)** | The core simulation modality for open quantum systems. This is what `run_lindbladian` and `run_thermalization` provide. | ALREADY DONE | Working for all four domains (Bohr, Energy, Time, Trotter). Both spectral analysis and step-by-step evolution implemented. |
| **Quantum trajectory / Monte Carlo wave function** | State-vector unraveling of the Lindbladian. Essential for scaling beyond small systems and for physical interpretation. Both QuantumOptics.jl and QuantumToolbox.jl provide this. | MEDIUM | `trajectories.jl` has `step_along_trajectory!` and `run_trajectories` but the code has unresolved issues (commented-out older framework, `trotter` variable not in scope in `build_trajectoryframework`). Needs validation against density matrix results. |
| **Spectral gap / mixing time analysis** | Users studying thermalization need to analyze the Liouvillian spectrum. This is core to the Chen/Kastoryano/Gilyen theory. | ALREADY DONE | `run_lindbladian` computes eigenvalues via Arpack, extracts spectral gap and steady state. `HotSpectralResults` stores these. |
| **Save/load results** | Users need to persist simulation results (density matrices, convergence curves) for later analysis and paper figures. | LOW | BSON is a dependency. `generate_filename` exists. But no standardized save/load for `HotAlgorithmResults` or `HotSpectralResults`. |

### Differentiators (Competitive Advantage)

Features that set QuantumFurnace apart. Not expected in generic packages, but high value for the target audience.

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| **Four approximation domains (Bohr/Energy/Time/Trotter)** | No other package implements the full Chen/Kastoryano/Gilyen hierarchy of approximation domains. This is the core intellectual contribution -- users can study how errors propagate from theory (Bohr) down to circuit-implementable (Trotter). | ALREADY DONE | Unique to QuantumFurnace. The domain dispatch pattern (`BohrDomain`, `EnergyDomain`, `TimeDomain`, `TrotterDomain`) is well-structured. Needs cross-domain comparison tooling (run same system in all domains, compare). |
| **Exact vs approximate detailed balance** | Supporting both KMS-DB (exact, with coherent term B) and GNS-DB (approximate, without coherent term) Lindbladians is unique. The theory paper distinguishes these and users will want to compare. | ALREADY DONE | `LiouvConfig` vs `LiouvConfigGNS` and `ThermalizeConfig` vs `ThermalizeConfigGNS`. Well-factored. |
| **Convex combination / transition function library** | The Metropolis-like and Glauber-like transition functions (parametrized by a, b) are a theoretical contribution. No other package provides this. | MOSTLY DONE | `pick_transition` and the (a, b) parametrization exist. Need documentation explaining the physics and when to use each. |
| **Log-Sobolev inequality (LSI) constant computation** | Computing the LSI constant alpha_2 via optimization gives rigorous mixing time bounds. This is cutting-edge -- most packages only provide heuristic convergence analysis. | ALREADY DONE | `compute_LSI_alpha2` in `log_sobolev.jl` with LBFGS optimization. Needs testing and documentation. |
| **Qiskit circuit generation from Trotter domain** | Translating the TrotterDomain simulation into actual quantum circuits bridges theory and experiment. No existing quantum Gibbs sampling package does this. Enables resource estimation (gate counts, circuit depth, qubit counts). | HIGH | Listed as "upcoming" in README. Not yet implemented. Requires mapping Trotter steps to gate sequences, handling controlled-time-evolution subroutines, and QPE-like energy estimation circuits. Depends on a Qiskit or Yao.jl integration. |
| **Cross-domain error analysis** | Quantifying the error introduced at each approximation level (Bohr->Energy, Energy->Time, Time->Trotter) is a key theoretical contribution. | MEDIUM | `errors.jl` has skeleton code. `compute_energy_quadrature_error` partially implemented. Time and Trotter error functions are stubs. Completing this enables the core paper narrative. |
| **NUFFT-accelerated operator filtering** | Using non-uniform FFT (FINUFFT) for the operator Fourier transform (OFT) in the Time/Trotter domains is a performance differentiator. | ALREADY DONE | `nufft.jl` and `ofts.jl` implement this. Unique to QuantumFurnace. |
| **Paper figure reproduction pipeline** | A companion package that exactly reproduces all figures in the associated paper is extremely valuable for the community and for peer review. | MEDIUM | No pipeline yet. Need scripts/notebooks that produce each paper figure from simulation data. |
| **Coherent term (B operator) computation** | The coherent correction term that makes the Lindbladian exactly detailed balanced. Computing B efficiently is non-trivial and unique to this algorithm. | ALREADY DONE | `coherent.jl` with `coherent_bohr`, `B_time`, `B_trotter`, and NUFFT-accelerated variants. |

### Anti-Features (Commonly Requested, Often Problematic)

Features that seem good but create problems. Deliberately NOT building these.

| Feature | Why Requested | Why Problematic | Alternative |
|---------|---------------|-----------------|-------------|
| **General-purpose Lindblad master equation solver** | "Why not make it work for any Lindbladian, not just Gibbs sampling ones?" | QuantumOptics.jl and QuantumToolbox.jl already do this excellently. Competing on general ODE integration is a losing battle against DifferentialEquations.jl backends. It would dilute the package's unique value. | Stay focused on the Chen/Kastoryano/Gilyen algorithm. For general Lindbladians, recommend users use QuantumToolbox.jl and feed the steady state back for comparison. |
| **GPU acceleration** | QuantumToolbox.jl supports GPU. Seems like a performance win. | At 12 qubits (4096x4096 matrices), GPU overhead dominates. The NUFFT operations use FINUFFT which is CPU-optimized. GPU would require rewriting the core kernels. Premature optimization for the current scale. | Optimize CPU path first (threading, BLAS tuning). Consider GPU only if targeting >14 qubits later. |
| **Symbolic Hamiltonian construction** | "Let me write H = J*XX + h*Z symbolically." | Adds massive dependency (SymbolicUtils.jl, Symbolics.jl). The existing `HamHam` constructor with term/coefficient vectors is clear enough. Symbolic manipulation is a separate problem domain. | Provide convenience constructors (`heisenberg_chain`, `transverse_ising`, `xxz_chain`) that hide the term/coefficient boilerplate. |
| **Real-time visualization / live plotting** | "Show me the state evolving in real time." | Adds Plots.jl or Makie.jl as hard dependencies, massively increasing load time. Most users will post-process in notebooks. | Return convergence data as arrays. Provide plotting recipes (optional extension via `QuantumFurnacePlots.jl` or Plots recipe if needed). Keep core package plot-free. |
| **Interoperability with QuantumOptics.jl types** | "Let me pass a QuantumOptics.jl Operator into QuantumFurnace." | Coupling to another package's type system creates maintenance burden and version conflicts. QuantumOptics.jl and QuantumToolbox.jl have different type hierarchies. | Use standard Julia types (Matrix{ComplexF64}, Vector{ComplexF64}). Provide conversion examples in docs. |
| **Distributed computing / cluster support** | ClusterManagers.jl and Distributed.jl are already dependencies. | Distributed trajectory averaging is theoretically easy but practically complex (data serialization, fault tolerance, load balancing). The current scale (12 qubits) does not require it. | Use `Base.Threads` for multi-threaded trajectory averaging on a single machine. Remove ClusterManagers/Distributed dependencies unless specifically needed. |
| **Automatic parameter optimization** | "Find me the best w0, t0, num_energy_bits for my system." | This is a research question, not a software feature. The optimal parameters depend on the Hamiltonian spectrum and desired accuracy in ways that are not yet fully understood theoretically. | Provide parameter sensitivity analysis tools and guidelines in documentation. Let the user make informed choices. |

## Feature Dependencies

```
[Hamiltonian Library]
    |
    v
[Correctness Tests] ──requires──> [Multiple Error Metrics]
    |
    v
[Cross-Domain Error Analysis] ──requires──> [All Four Domains Working]
    |                                             |
    v                                             v
[Paper Figure Pipeline] ──requires──> [Save/Load Results]
                                           |
                                           v
                              [Documented API + Tutorials]

[Trajectory Validation] ──requires──> [Correctness Tests]
    |
    v
[Qiskit Circuit Generation] ──requires──> [TrotterDomain Validated]
                                               |
                                               v
                                    [Resource Estimation]

[LSI Constant Computation] ──enhances──> [Spectral Gap Analysis]

[Exact DB] ──conflicts──> [GNS DB] (same simulation; user picks one)
```

### Dependency Notes

- **Correctness Tests require Multiple Error Metrics:** Cannot validate convergence without fidelity/trace distance/relative entropy to compare against known states.
- **Cross-Domain Error Analysis requires All Four Domains Working:** The whole point is comparing across domains; if any domain is broken, the comparison is meaningless.
- **Qiskit Circuit Generation requires TrotterDomain Validated:** Generating circuits from incorrect Trotter simulations is worse than useless. Must validate Trotter convergence first.
- **Paper Figure Pipeline requires Save/Load + Documented API:** Reproducible figures need persistent data and documented parameters.
- **Trajectory Validation requires Correctness Tests:** Must prove trajectory-averaged density matrices converge to the same state as the master equation approach.
- **LSI Constant Computation enhances Spectral Gap Analysis:** LSI gives tighter mixing time bounds than spectral gap alone, but spectral gap is the prerequisite.

## MVP Definition

### Launch With (v1 -- Paper Companion)

Minimum viable product -- what's needed for publication and first community users.

- [ ] **Correctness tests for all four domains** -- Reviewers will ask "how do you know this is correct?" Must demonstrate convergence to exact Gibbs state for 2-4 qubit systems.
- [ ] **Multiple error metrics (fidelity, trace distance, relative entropy)** -- Standard quantum information metrics that appear in every paper in this field.
- [ ] **Cross-domain error analysis (complete `errors.jl`)** -- Core narrative of the paper: show how errors accumulate Bohr -> Energy -> Time -> Trotter.
- [ ] **Trajectory validation** -- Fix `build_trajectoryframework` bugs, prove agreement with density matrix evolution.
- [ ] **Documented API (docstrings for all exports)** -- Every exported function needs a docstring. Minimum bar for a publishable package.
- [ ] **Working tutorials for each domain** -- One tutorial per domain showing thermalization of 4-qubit Heisenberg.
- [ ] **CI pipeline (GitHub Actions)** -- Tests run automatically on push. Non-negotiable for a Julia package.
- [ ] **2-3 additional Hamiltonians (Ising, XY, random)** -- Demonstrates generality beyond Heisenberg.

### Add After Validation (v1.x)

Features to add once core is working and paper is submitted.

- [ ] **Paper figure reproduction scripts** -- Once the paper draft is stable, lock down the exact scripts that produce each figure.
- [ ] **Qiskit circuit generation (TrotterDomain -> circuits)** -- Add after Trotter domain is validated. High impact but high complexity.
- [ ] **Resource estimation (gate counts, depth, qubit overhead)** -- Builds on circuit generation. Quantifies the "cost" of the quantum algorithm.
- [ ] **Save/load with metadata (parameters, git hash, timestamp)** -- For long-running simulations and reproducibility.
- [ ] **Parameter sensitivity analysis tools** -- How do results change with w0, t0, num_energy_bits? Sweep tools.

### Future Consideration (v2+)

Features to defer until paper is published and community feedback arrives.

- [ ] **Large-scale trajectory simulation (MPS-based)** -- For >12 qubits, dense matrices are infeasible. Would need ITensors.jl integration. Fundamentally different algorithm.
- [ ] **Automatic differentiation through the Lindbladian** -- For gradient-based optimization of Hamiltonian parameters. Requires Zygote/Enzyme compatibility.
- [ ] **Plotting recipes / visualization extension** -- Optional `QuantumFurnacePlots.jl` with convergence plots, Bloch sphere, spectrum visualizations.
- [ ] **Non-Markovian extensions** -- HEOM or other memory-kernel approaches. Different physics entirely.
- [ ] **Fermionic Hamiltonians (Jordan-Wigner)** -- Recent work on Gibbs sampling for Fermi-Hubbard. Would broaden audience significantly but is a major extension.

## Feature Prioritization Matrix

| Feature | User Value | Implementation Cost | Priority |
|---------|------------|---------------------|----------|
| Correctness tests | HIGH | MEDIUM | P1 |
| Multiple error metrics | HIGH | LOW | P1 |
| Documented API | HIGH | MEDIUM | P1 |
| Working tutorials | HIGH | MEDIUM | P1 |
| CI pipeline | HIGH | LOW | P1 |
| Cross-domain error analysis | HIGH | MEDIUM | P1 |
| Trajectory validation | HIGH | MEDIUM | P1 |
| Additional Hamiltonians | MEDIUM | LOW | P1 |
| Paper figure pipeline | HIGH | MEDIUM | P2 |
| Save/load with metadata | MEDIUM | LOW | P2 |
| Qiskit circuit generation | HIGH | HIGH | P2 |
| Resource estimation | HIGH | HIGH | P2 |
| Parameter sensitivity tools | MEDIUM | MEDIUM | P2 |
| LSI documentation + tests | MEDIUM | LOW | P2 |
| MPS-based trajectories | MEDIUM | HIGH | P3 |
| AD through Lindbladian | LOW | HIGH | P3 |
| Plotting extension | LOW | MEDIUM | P3 |
| Non-Markovian extensions | LOW | HIGH | P3 |
| Fermionic Hamiltonians | MEDIUM | HIGH | P3 |

**Priority key:**
- P1: Must have for paper submission
- P2: Should have, add when core is solid
- P3: Nice to have, future consideration

## Competitor Feature Analysis

| Feature | QuantumOptics.jl | QuantumToolbox.jl | QuTiP | QuantumFurnace.jl |
|---------|-----------------|-------------------|-------|-------------------|
| General Lindblad ME solver | Full (any Lindbladian) | Full (any Lindbladian) | Full (any Lindbladian) | Specialized (Gibbs sampling Lindbladians only) |
| Quantum trajectory / MCWF | Standard MCWF | Standard MCWF + distributed | Standard MCWF | Custom unraveling for Chen algorithm (unique) |
| Hamiltonian library | Spin, Fock, NLevel, ManyBody | Same as QuTiP | Extensive predefined | Heisenberg only (needs expansion) |
| Approximation domain hierarchy | N/A | N/A | N/A | **Unique: Bohr/Energy/Time/Trotter** |
| Exact detailed balance | N/A | N/A | N/A | **Unique: KMS-DB with coherent term** |
| Circuit generation | N/A | N/A | Via Qiskit (separate) | Planned (TrotterDomain -> Qiskit) |
| Error metrics | trace_distance, fidelity, entropy, negativity | Same as QuTiP | Comprehensive suite | trace_distance only (needs expansion) |
| Steady state solvers | Eigenvalue, direct, iterative | Eigenvalue, direct, Fourier | Eigenvalue, direct, iterative | Eigenvalue (Arpack eigs) |
| Spectral analysis | Eigenspectrum | Liouvillian eigenspectrum | Eigenspectrum | Spectral gap + LSI constant (unique) |
| GPU support | No | Yes (CUDA.jl) | Partial | No (not needed at current scale) |
| Documentation quality | Excellent (Documenter.jl + tutorials) | Excellent (QuTiP-style) | Excellent (Sphinx) | Partial (structure exists, content sparse) |
| CI/testing | Full (GitHub Actions) | Full (GitHub Actions) | Full | None (critical gap) |
| NUFFT acceleration | N/A | N/A | N/A | **Unique: FINUFFT for OFT** |
| Transition function library | N/A | N/A | N/A | **Unique: Gaussian, Metropolis, Glauber** |
| LSI mixing time bounds | N/A | N/A | N/A | **Unique** |

## Sources

### Verified (HIGH confidence)
- [QuantumToolbox.jl GitHub](https://github.com/qutip/QuantumToolbox.jl) -- feature list, architecture
- [QuantumToolbox.jl paper (Quantum, 2025)](https://quantum-journal.org/papers/q-2025-09-29-1866/) -- published capabilities
- [QuantumOptics.jl documentation](https://docs.qojulia.org/api/) -- full API reference
- [QuantumOptics.jl paper](https://arxiv.org/pdf/1707.01060) -- design philosophy
- [QuTiP 5.3 documentation](https://qutip.readthedocs.io/en/stable/) -- feature comparison baseline
- [Chen, Kastoryano, Gilyen 2023 (arXiv:2311.09207)](https://arxiv.org/abs/2311.09207) -- the algorithm QuantumFurnace implements
- [Chen, Kastoryano, Brandao, Gilyen 2023 (arXiv:2303.18224)](https://arxiv.org/abs/2303.18224) -- quantum thermal state preparation theory
- QuantumFurnace.jl codebase -- direct analysis of all source files

### Partially verified (MEDIUM confidence)
- [QuantumToolbox.jl introduction page](https://qutip.org/QuantumToolbox.jl/) -- GPU and AD claims verified by docs but not tested
- [Efficient Quantum Gibbs Samplers with KMS (arXiv:2404.05998)](https://arxiv.org/abs/2404.05998) -- Ding 2024 extensions referenced in project context
- [Lindblad engineering for Gibbs state preparation (Quantum, 2025)](https://quantum-journal.org/papers/q-2025-08-29-1843/) -- competitor algorithmic approaches

### Unverified (LOW confidence)
- Qiskit circuit generation feasibility for the Chen algorithm -- no existing implementation found. Complexity estimate based on general Trotter-to-circuit mapping knowledge.
- MPS-based trajectory scaling claims -- based on the tensor jump method paper but not validated for this specific algorithm.

---
*Feature research for: Quantum Gibbs sampling / open quantum systems simulation*
*Researched: 2026-02-13*
