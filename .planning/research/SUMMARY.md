# Project Research Summary

**Project:** QuantumFurnace.jl — Quantum Gibbs Sampling Simulation Package
**Domain:** Computational quantum physics (Julia scientific computing package)
**Researched:** 2026-02-13
**Confidence:** HIGH

## Executive Summary

QuantumFurnace.jl is a Julia package implementing the Chen/Kastoryano/Gilyen quantum Gibbs sampling algorithm for open quantum systems. The package simulates thermalization through Lindbladian dynamics across four approximation domains (Bohr, Energy, Time, Trotter), providing a unique research tool that bridges theoretical quantum algorithms with practical implementation. The current codebase is feature-complete for the core algorithm but requires systematic testing infrastructure, trajectory validation, and stabilization before publication.

The recommended approach is to focus on validation and stabilization first, building a comprehensive test suite around the existing simulation engine, then expanding documentation and adding comparative features (additional Hamiltonians, error metrics). The trajectory simulation path has critical bugs that mask normalization errors and must be fixed before any trajectory-based results can be trusted. The package should target the Julia General Registry after validation, positioning itself as a specialized tool complementing existing general-purpose quantum packages (QuantumOptics.jl, QuantumToolbox.jl) rather than competing with them.

Key risks include: (1) trajectory normalization bugs hidden by fallback logic, (2) numerical instability in the Cholesky-based CPTP channel construction, (3) hardcoded tolerances that break at larger qubit counts, and (4) missing test infrastructure preventing regression detection. These are addressable through systematic testing, careful numerical validation against known Gibbs states, and refactoring tolerance values to be parameter-dependent rather than hardcoded.

## Key Findings

### Recommended Stack

QuantumFurnace.jl has a solid foundation of standard Julia scientific computing packages. The critical recommendation is to replace Arpack.jl with KrylovKit.jl for sparse eigensolve (thread-safety and Hermitian specialization), add ExponentialUtilities.jl for Krylov-based matrix exponentials (scaling to larger systems), and modernize Python interop by using PythonCall.jl instead of PyCall.jl for future Qiskit circuit generation. The existing dependencies (LinearAlgebra, SparseArrays, FINUFFT, LinearMaps, QuadGK, Optim) are well-chosen and should be retained.

**Core technologies (keep/already present):**
- Julia 1.12+ — already in use, provides redefinable types and improved compilation
- FINUFFT.jl >= 3.4.2 — non-uniform FFT, critical for time/energy domain OFT acceleration
- LinearMaps.jl — matrix-free Liouvillian for scaling beyond 10 qubits
- SparseArrays + LinearAlgebra (stdlib) — standard for quantum Hamiltonians
- Optim.jl — log-Sobolev bound computation
- BSON.jl — result serialization (consider migrating to JLD2.jl long-term)

**Add for stabilization and scaling:**
- KrylovKit.jl >= 0.8 — replace Arpack (thread-safe, Hermitian-specialized eigensolve)
- ExponentialUtilities.jl >= 1.27 — Krylov expmv for trajectory simulation without dense matrix exponentials
- PythonCall.jl + CondaPkg.jl — modern Julia-Python bridge for Qiskit circuit generation
- OrdinaryDiffEqTsit5.jl — lightweight ODE solver for master equation validation
- Aqua.jl >= 0.8 — automated package quality tests (method ambiguities, stale deps)
- SafeTestsets.jl — isolated test sets preventing state leakage

**Remove or move:**
- Arpack.jl — replace with KrylovKit (not thread-safe, lacks Hermitian specialization)
- Revise/Debugger/BenchmarkTools — move from [deps] to [extras] (dev-only tools)
- Plots.jl — consider moving to extras or separate extension (heavy dependency)

### Expected Features

The feature landscape reveals a package that implements cutting-edge theory (four-domain hierarchy, exact KMS detailed balance, NUFFT-accelerated OFT) but lacks the packaging and validation infrastructure to be usable by the community.

**Must have (table stakes for publication):**
- Correctness tests against known Gibbs states (2-4 qubit systems with analytical solutions)
- Multiple error metrics (fidelity, trace distance, relative entropy) — currently only trace distance exists
- Documented API with docstrings on all exported functions — many functions lack docstrings
- CI pipeline (GitHub Actions) — infrastructure exists but no `test/runtests.jl` to run
- Working tutorials for each domain — Literate.jl infrastructure exists, content is sparse
- Trajectory validation — density matrix and trajectory paths must converge to same steady state
- Additional Hamiltonians (Ising, XY, XXZ, random local) — currently only Heisenberg

**Should have (differentiators):**
- Cross-domain error analysis (complete `errors.jl` stubs) — quantify Bohr→Energy→Time→Trotter errors
- Paper figure reproduction pipeline — scripts that generate all paper figures
- Qiskit circuit generation from TrotterDomain — bridges theory and experiment, enables resource estimation
- LSI constant computation documentation + tests — unique feature, needs validation

**Defer (v2+, post-publication):**
- MPS-based trajectory simulation for >12 qubits — different algorithm class, needs ITensors.jl
- GPU acceleration — overhead dominates at <=12 qubits, CPU path is appropriate scale
- Automatic differentiation through Lindbladian — specialized use case, Zygote/Enzyme compatibility uncertain
- Plotting recipes / visualization extension — keep core package plot-free, provide examples

### Architecture Approach

The package follows a well-factored domain-dispatch pattern using Julia's multiple dispatch for extensibility. Four singleton domain types (BohrDomain, EnergyDomain, TimeDomain, TrotterDomain) drive method selection, with precomputation caching to eliminate inner-loop allocation. The architecture supports dual simulation modes: full Liouvillian construction for spectral analysis, and stochastic trajectory simulation for sampling statistics. This dual-mode pattern matches established quantum packages (QuantumOptics.jl, QuantumToolbox.jl).

**Major components:**
1. **Domain Dispatch Layer** (`jump_workers.jl`) — multiple `jump_contribution!` methods specialized per domain, central hub for Liouvillian assembly and Kraus evolution
2. **Precomputation Layer** (`furnace_utensils.jl`, `nufft.jl`, `ofts.jl`) — caches NUFFT prefactors, transition functions, energy labels to avoid allocation in tight loops
3. **Data Model** (`structs.jl`, `hamiltonian.jl`, `trotter_domain.jl`) — `HamHam` (spectral decomposition + Bohr frequencies), `TrottTrott` (Trotterized evolution), domain/config types, workspaces
4. **Coherent Term Engine** (`coherent.jl`) — B operator computation for exact KMS detailed balance, per-domain methods
5. **Entry Points** (`furnace.jl`, `trajectories.jl`) — `run_lindbladian` (spectral analysis), `run_thermalization` (density matrix), `run_trajectories` (state vector)

**Key architectural patterns:**
- Domain dispatch via singleton types — open/closed principle, compiler specialization
- Workspace/precomputation cache — zero-allocation inner loops, named tuples for self-documentation
- Dual simulation modes — density matrix master equation + trajectory unraveling with cross-validation

**Scaling implications:**
- 3-10 qubits: current dense architecture works well
- 11-12 qubits: Liouvillian is 4M-16M entries (128MB-1GB), trajectory mode becomes preferred
- 13+ qubits: must use LinearMaps (matrix-free) or trajectory-only, dense Liouvillian infeasible (>64GB)

### Critical Pitfalls

Six critical pitfalls identified from codebase analysis, prioritized by severity:

1. **Trajectory fallback logic silently masks probability normalization bugs** — The `!chosen` fallback in `step_along_trajectory!` applies the last candidate state when cumulative probability scan fails, hiding inconsistencies between `precompute_R()` and per-jump rate formulas. Add assertion: `|p_nojump + p_res + p_jump_total - 1.0| < epsilon`.

2. **Cholesky residual S matrix can silently go non-PSD, breaking CPTP guarantee** — `S = (2*alpha - delta)*R - alpha^2 * R^2` with `check=false` Cholesky produces garbage when S has negative eigenvalues. The `10*eps(Float64)` shift is too small for strong coupling. Replace with eigendecomposition-based square root or adaptive epsilon shift.

3. **Hardcoded 1e-12 tolerance used for incompatible purposes** — Same value serves as frequency cutoff, integration tolerance, NUFFT accuracy, and LSI regularization. Breaks when `nu_min < 1e-12` at larger qubit counts. Define named constants: `FREQ_ZERO_TOL = nu_min/10`, `INTEGRATION_TOL = 1e-12`, `NUFFT_EPS = 1e-12`.

4. **Density matrix and trajectory paths use different delta scaling conventions** — DM thermalization applies `delta_scale = 1.0/p_jump` per-jump rescaling, trajectory framework does not. Mathematically equivalent channels require matched conventions. Add channel equivalence test before cross-validation.

5. **`unsafe_wrap` pointer arithmetic in log-Sobolev optimization can corrupt memory** — `pointer(A_flat, n_params * sizeof(Float64) + 1)` confuses element indexing with byte offsets, overshooting by 8x. Replace with `reshape(view(...))` or fix to `pointer(A_flat, n_params + 1)`.

6. **Coherent term B has undefined `trotter` variable in trajectory framework** — `build_trajectoryframework` passes `trotter=trotter` keyword but variable undefined. Crashes when `with_coherent=true`. Remove keyword, use `ham_or_trott` directly.

## Implications for Roadmap

Based on dependencies discovered in research, the optimal phase structure addresses foundational stability first, then systematic validation, then community-facing features. The current codebase is feature-complete for the core algorithm but unstable for research use.

### Phase 1: Core Stabilization and Testing Infrastructure
**Rationale:** Cannot validate or extend a package without a working test suite. Critical bugs in trajectories and numerical stability must be fixed before any results can be trusted. The architecture research shows test files use `includet` (interactive scripts) instead of proper `@testset` blocks — this prevents CI and regression detection.

**Delivers:**
- `test/runtests.jl` with `@testset` structure runnable via `Pkg.test()`
- Fix trajectory probability normalization (Pitfall 1)
- Fix Cholesky non-PSD handling (Pitfall 2)
- Fix coherent term `trotter` undefined variable (Pitfall 6)
- Fix `unsafe_wrap` memory corruption in LSI optimization (Pitfall 5)
- CI configuration running tests on Julia 1.12 and latest

**Addresses features:**
- CI pipeline (table stakes)
- Proper Julia package structure (table stakes)

**Avoids pitfalls:**
- All six critical pitfalls require test coverage to detect and prevent regression
- "Looks done but isn't" trap — trajectory framework appears functional but crashes with `with_coherent=true`

**Research flags:** None — standard Julia testing patterns, well-documented.

### Phase 2: Validation and Correctness
**Rationale:** The core algorithm's theoretical correctness must be validated against analytically known systems. Architecture research shows this is the validation burden for dual-mode simulation — density matrix and trajectory paths must produce identical steady states. Feature research identifies this as table stakes for publication.

**Delivers:**
- Correctness tests for all four domains against known Gibbs states (1-4 qubit systems)
- Multiple error metrics (fidelity, relative entropy, operator norm) in `qi_tools.jl`
- Cross-domain validation (BohrDomain result matches EnergyDomain for same parameters)
- Trajectory validation (trajectory-averaged rho matches DM-evolved rho)
- Channel equivalence tests (DM and trajectory implement same CPTP map)

**Addresses features:**
- Correctness tests (table stakes)
- Multiple error metrics (table stakes)
- Trajectory validation (table stakes)

**Uses stack:**
- Add OrdinaryDiffEqTsit5 for master equation reference integration
- Keep FINUFFT, LinearMaps, existing numerics stack

**Avoids pitfalls:**
- Trajectory fallback masking bugs (Pitfall 1) — probability sum assertion catches this
- Delta scaling convention mismatch (Pitfall 4) — channel equivalence test detects divergence
- Testing only final Gibbs state without approach dynamics — decay curves must match

**Research flags:** Minor — ODE solver setup is standard, but cross-validation between stochastic trajectories and master equation needs careful statistical error analysis.

### Phase 3: Refactoring and Scaling
**Rationale:** With validated correctness, address technical debt that limits scaling. Hardcoded tolerances (Pitfall 3) break at larger qubit counts, and the dense Liouvillian approach becomes infeasible beyond 10 qubits. Stack research recommends KrylovKit and ExponentialUtilities for scaling.

**Delivers:**
- Replace hardcoded `1e-12` with named constants parameterized by system properties (`nu_min`)
- Replace Arpack with KrylovKit for eigensolve
- Add ExponentialUtilities for Krylov expmv in trajectory evolution
- Optionally refactor config structs (4 variants with 90% overlap is technical debt)
- Move dev-only deps (Revise, Debugger) from [deps] to [extras]

**Addresses features:**
- Improved scaling (supporting feature for larger systems)
- Cleaner package structure (Julia registry preparation)

**Uses stack:**
- Add KrylovKit.jl >= 0.8
- Add ExponentialUtilities.jl >= 1.27
- Remove Arpack.jl

**Avoids pitfalls:**
- Hardcoded tolerance pitfall (Pitfall 3) — frequency cutoff becomes parameter-dependent
- Dense Liouvillian scaling trap — prepare for LinearMaps transition

**Research flags:** None — KrylovKit and ExponentialUtilities are well-documented with clear migration paths from Arpack.

### Phase 4: Documentation and API Stabilization
**Rationale:** With correctness validated and API stable, systematic documentation enables community use. Architecture research shows Documenter.jl + Literate.jl infrastructure exists but content is sparse. Feature research identifies this as table stakes — every competing package has full API docs and tutorials.

**Delivers:**
- Docstrings on all exported functions (API reference via `@autodocs`)
- Working Literate.jl tutorials for each domain (Bohr, Energy, Time, Trotter)
- Theory pages explaining KMS detailed balance, transition functions, LSI bounds
- CI auto-deployment of documentation site
- Add Aqua.jl tests for package quality (method ambiguities, missing compat entries)

**Addresses features:**
- Documented API (table stakes)
- Working tutorials (table stakes)

**Uses stack:**
- Documenter.jl >= 1.16 (already present)
- Literate.jl >= 2.20 (already present)
- Add Aqua.jl >= 0.8

**Avoids pitfalls:**
- Undocumented conventions (e.g., delta scaling differences between DM and trajectory modes)
- Incomplete examples hiding API surface bugs

**Research flags:** None — Documenter + Literate is standard Julia documentation stack.

### Phase 5: Feature Expansion
**Rationale:** With stable, documented core, add differentiators and comparative features. Feature research identifies cross-domain error analysis and additional Hamiltonians as high-value, moderate-complexity additions.

**Delivers:**
- Complete `errors.jl` stubs (cross-domain error quantification)
- Additional Hamiltonian constructors (Ising, XY, XXZ, random local Hamiltonians)
- Save/load with metadata (parameters, git hash, timestamp)
- Parameter sensitivity analysis tools (sweep `w0`, `t0`, `num_energy_bits`)
- LSI constant computation tests and documentation

**Addresses features:**
- Cross-domain error analysis (differentiator, key paper narrative)
- Hamiltonian library beyond Heisenberg (table stakes for generality)
- Save/load results (table stakes for reproducibility)
- LSI documentation + tests (differentiator, unique feature)

**Avoids pitfalls:**
- Testing only Heisenberg chain (limits generalization claims)
- No persistent result storage (prevents paper figure reproducibility)

**Research flags:** Low — cross-domain error formulas are derived in Chen/Kastoryano/Gilyen papers, implementation is standard numerical integration.

### Phase 6: Qiskit Interop and Circuit Generation (Optional, Post-Publication)
**Rationale:** Bridges theory and experiment, enables resource estimation. High complexity, depends on validated TrotterDomain. Feature research identifies this as a key differentiator but not required for paper submission. Stack research recommends PythonCall.jl over PyCall.jl.

**Delivers:**
- PythonCall.jl + CondaPkg.jl setup
- TrotterDomain gate sequence extraction
- Qiskit QuantumCircuit construction from gate sequences
- Resource estimation (gate counts, circuit depth, qubit overhead)
- Cross-validation test (Julia Trotter unitary == Python Trotter unitary)

**Addresses features:**
- Qiskit circuit generation (differentiator, high impact)
- Resource estimation (differentiator, quantifies algorithm cost)

**Uses stack:**
- Add PythonCall.jl >= 0.9.23
- Add CondaPkg.jl
- Qiskit >= 1.0 (via CondaPkg)

**Avoids pitfalls:**
- Python convention mismatch (rescaling factors must be exported from Julia)
- Using PyCall.jl (deprecated, type-unstable, single-threaded GIL)

**Research flags:** MEDIUM — PythonCall.jl is well-documented, but specific Qiskit integration patterns for Julia-generated circuits are not widely documented. May need research-phase for circuit generation API patterns.

### Phase Ordering Rationale

- **Phase 1 before all others:** Cannot validate correctness or extend features without a working test suite. Critical bugs make trajectory results untrustworthy.
- **Phase 2 before refactoring/features:** Must know the algorithm is correct before optimizing or extending it. Cross-validation between DM and trajectory modes is the fundamental correctness check.
- **Phase 3 before documentation:** API may change during refactoring (e.g., config struct consolidation). Documenting unstable APIs creates maintenance burden.
- **Phase 4 before expansion:** Stable, documented API prevents breaking changes when adding features. New Hamiltonians and error analysis build on existing primitives.
- **Phase 6 as optional post-publication:** Qiskit interop is high-value but high-complexity. Paper submission does not require it — simulation results stand alone. Defer to v1.x after community feedback.

**Dependency chain:** Testing → Validation → Refactoring → Documentation → Features → Qiskit (each phase depends on previous completing)

### Research Flags

**Phases needing deeper research during planning:**
- **Phase 6 (Qiskit interop):** Python-Julia bridge for circuit generation is not widely documented. May need `/gsd:research-phase` to investigate Qiskit API patterns, gate set mappings, and resource estimation tools.

**Phases with standard patterns (skip research-phase):**
- **Phase 1 (Testing):** Julia Test module is standard, documented extensively
- **Phase 2 (Validation):** Numerical validation against known solutions is standard scientific computing
- **Phase 3 (Refactoring):** KrylovKit migration from Arpack is well-documented, named constants are trivial
- **Phase 4 (Documentation):** Documenter + Literate is the Julia ecosystem standard
- **Phase 5 (Features):** Cross-domain error formulas are derived in theory papers, Hamiltonian construction is standard quantum mechanics

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Stack | HIGH | Core dependencies (FINUFFT, LinearMaps, Optim) verified from codebase. KrylovKit/ExponentialUtilities recommendations verified from official docs and ecosystem usage. PythonCall vs PyCall comparison verified from official comparison page. |
| Features | HIGH | Table stakes derived from direct comparison with QuantumOptics.jl, QuantumToolbox.jl, QuTiP feature sets. Differentiators (four-domain hierarchy, KMS detailed balance, NUFFT OFT) verified from codebase analysis. MVP definition grounded in publication requirements. |
| Architecture | HIGH | Domain dispatch pattern, dual-mode simulation, workspace caching all verified from source code. Comparison with QuantumOptics.jl and QuantumToolbox.jl architectures verified from GitHub repos and published papers. |
| Pitfalls | HIGH | All six critical pitfalls derived from direct codebase analysis with line numbers. Trajectory normalization fallback, Cholesky non-PSD, `unsafe_wrap` pointer arithmetic, undefined `trotter` variable all verified in source. |

**Overall confidence:** HIGH

Research is grounded in codebase analysis (direct file inspection) supplemented by verified external sources (official documentation, published papers, GitHub repositories of comparable packages). Recommendations are specific and actionable with clear rationale.

### Gaps to Address

**Gaps requiring validation during implementation:**

- **Trajectory statistical error analysis:** Research identifies that trajectory-averaged density matrices must match DM evolution within `O(1/sqrt(ntraj))`, but the exact prefactor and how many trajectories are needed for convergence is system-dependent. During Phase 2, empirically determine trajectory count requirements for 2-4 qubit test systems.

- **NUFFT accuracy vs. performance tradeoff:** Stack research notes `eps=1e-12` may be tighter than needed. During Phase 3, profile FINUFFT time vs. accuracy for 8-12 qubit systems to determine if `eps=1e-8` suffices (Gaussian filter provides natural smoothing in energy domain).

- **LinearMaps scaling transition point:** Architecture research suggests dense Liouvillian becomes infeasible at 13+ qubits, but the exact transition depends on available memory and iterative solver convergence. During Phase 3, benchmark to determine when to switch from dense to matrix-free.

- **Qiskit gate set mapping:** Feature and stack research recommend Qiskit interop via PythonCall.jl, but the specific mapping from Trotter unitary decomposition to Qiskit gate sequences is not documented. Phase 6 (if pursued) will need research-phase to investigate Qiskit's gate synthesis tools and best practices for controlled time evolution circuits.

**Technical debt to track:**

- ~200 lines of commented-out code in `trajectories.jl` and `jump_workers.jl` should be deleted and git history used for reference (noted in pitfalls research)
- Four config struct variants with 90% field overlap should be consolidated or justified (noted in architecture anti-patterns)
- `verify_completeness` is exported but entirely commented out (noted in "looks done but isn't" checklist)

## Sources

### Primary (HIGH confidence)
- QuantumFurnace.jl codebase — direct analysis of all source files (stack, features, architecture, pitfalls)
- [KrylovKit.jl official docs](https://jutho.github.io/KrylovKit.jl/stable/) — Hermitian Lanczos specialization, API
- [PythonCall.jl vs PyCall comparison](https://juliapy.github.io/PythonCall.jl/v0.2/pycall/) — feature comparison, type stability
- [QuantumToolbox.jl paper (Quantum, 2025)](https://quantum-journal.org/papers/q-2025-09-29-1866/) — architecture patterns, benchmarks, GPU support
- [QuantumOptics.jl GitHub](https://github.com/qojulia/QuantumOptics.jl) — src structure, master.jl + mcwf.jl dual mode
- [QuantumToolbox.jl GitHub](https://github.com/qutip/QuantumToolbox.jl) — feature list, Qobj design, multiple dispatch
- [Documenter.jl v1.16.1](https://documenter.juliadocs.org/stable/) — documentation generation
- [Julia Test module docs](https://docs.julialang.org/en/v1/stdlib/Test/) — testing best practices
- [Chen, Kastoryano, Gilyen 2023 (arXiv:2311.09207)](https://arxiv.org/abs/2311.09207) — the algorithm QuantumFurnace implements

### Secondary (MEDIUM confidence)
- [ExponentialUtilities.jl (SciML docs)](https://docs.sciml.ai/ExponentialUtilities/stable/) — Krylov expmv verified but version from training data
- [OrdinaryDiffEq.jl split packages](https://github.com/SciML/OrdinaryDiffEq.jl) — split approach verified, version from training data
- [Aqua.jl (JuliaTesting, v0.8.14)](https://github.com/JuliaTesting/Aqua.jl) — package quality checks
- [FINUFFT.jl v3.x docs](https://ludvigak.github.io/FINUFFT.jl/latest/) — NUFFT API and performance
- [Bridging Worlds: Julia-Python Interoperability (arXiv:2404.18170)](https://arxiv.org/html/2404.18170v1) — PythonCall vs PyCall analysis
- [Julia Package Testing Best Practices](https://blog.glcs.io/package-testing) — test organization patterns

### Tertiary (LOW confidence, needs validation)
- Qiskit circuit generation feasibility for Chen algorithm — no existing implementation found, complexity estimate based on general Trotter-to-circuit mapping knowledge
- MPS-based trajectory scaling claims — based on tensor jump method literature but not validated for this specific algorithm
- Exact trajectory count for statistical convergence — system-dependent, requires empirical determination

---
*Research completed: 2026-02-13*
*Ready for roadmap: yes*
