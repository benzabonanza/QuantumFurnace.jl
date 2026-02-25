# QuantumFurnace.jl

## What This Is

A Julia package for simulating quantum Gibbs sampling via Lindbladian evolution. It implements multiple Lindbladian constructions (GNS and KMS detailed balance) across a hierarchy of approximation domains (Bohr, Energy, Time, Trotter), enabling density matrix simulation, trajectory-based sampling, and matrix-free Krylov spectral gap estimation. The package features multi-threaded trajectory sampling with adaptive convergence-driven batching, spectral gap estimation from both trajectory-based observable decay and matrix-free KrylovKit eigsolve, exact Lindbladian diagnostics with biorthogonal eigenvector analysis, BSON experiment serialization, and validated KMS-vs-GNS comparison capabilities. Core structs are parameterized on element type `{T<:AbstractFloat}` for precision flexibility. The package targets researchers and students working on quantum Gibbs sampling, providing both fast numerical simulations and pedagogic documentation grounded in the theoretical literature.

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
- ✓ Trajectory compilation and correctness fixes (U_B ordering, PSD guard, normalization) -- v1.0
- ✓ CPTP channel completeness verification for Energy, Time, Trotter domains -- v1.0
- ✓ BohrDomain detailed balance (Gibbs fixed point to 1.6e-15 trace distance) -- v1.0
- ✓ Domain error hierarchy: dist_bohr <= dist_energy <= dist_time <= dist_trotter -- v1.0
- ✓ DM step error scaling: O(delta^2) single-step, O(delta) multi-step accumulation -- v1.0
- ✓ Coherent term B and OFT consistency across domains -- v1.0
- ✓ Trajectory-vs-DM cross-validation for Energy, Time, Trotter domains -- v1.0
- ✓ Per-operator Lie-Trotter splitting in trajectory simulation -- v1.0
- ✓ TrotterDomain Gibbs convergence with coherent term -- v1.0
- ✓ 1/sqrt(N) trajectory convergence rate verification -- v1.0
- ✓ DM-based regression tests (portable across Julia versions and platforms) -- v1.0
- ✓ Aqua.jl package quality gate with compat bounds -- v1.0
- ✓ Test infrastructure with shared fixtures and tiered tolerances -- v1.0
- ✓ Dead code pruning: ~987 lines commented code, ~35 unused functions, 3 dead structs removed -- v1.1
- ✓ DRY helpers: hermitianize!, apply_cptp_channel!, apply_coherent_unitary! -- v1.1
- ✓ Struct simplification: immutable TrottTrott, fully-initialized HamHam, deduplicated config constructors -- v1.1
- ✓ Type parameterization: HamHam{T}, TrottTrott{T}, Config{D,T}, LindbladianWorkspace{T} -- v1.1
- ✓ API surface cleanup: organized exports, ~18 internals _-prefixed, trace_distance_h exported -- v1.1
- ✓ Allocation optimization: index-based B_bohr, in-place phase rotations, precomputed basis transforms -- v1.1
- ✓ Multi-threaded trajectory sampling with per-thread workspace, BLAS control, and deterministic Xoshiro seeding -- v1.2
- ✓ GNS (approximate, no B term) trajectory path via ThermalizeConfigGNS, converging to GNS fixed point -- v1.2
- ✓ Adaptive trajectory sampling with convergence-driven batching (relative change <1% for 3 consecutive batches, hard cap) -- v1.2
- ✓ Convergence tracking: trace distance to Gibbs + per-observable (ZZ correlations, energy) at batch checkpoints -- v1.2
- ✓ BSON-based ExperimentResult serialization with metadata (git hash, timestamp, seed, thread count) -- v1.2
- ✓ KMS-vs-GNS parameter sweep experiments confirming KMS achieves lower trace distance than GNS -- v1.2
- ✓ Logic simplification: flattened call chain (5->3 levels), domain-aware JumpOp construction, simplified result structs -- v1.2
- ✓ Observable-only trajectory runner with optional DM reconstruction (`run_observable_trajectories`) -- v1.3
- ✓ 8-observable bundle (`build_preset_trajectory_observables`): H, Mz, XX_avg, YY_avg, ZZ_avg, Mz_stagg, Z1, XZ_stagg -- v1.3
- ✓ Single-exponential decay fitting with auto log-linear initial guess, CIs, R-squared (`fit_exponential_decay`) -- v1.3
- ✓ Single-call spectral gap estimation API (`estimate_spectral_gap` → `SpectralGapResult`) -- v1.3
- ✓ Eigenbasis overlap diagnostic (`eigenbasis_overlap_analysis` → `OverlapAnalysisResult`) -- v1.3
- ✓ n=4 spectral gap cross-validated to 0.72% against exact Liouvillian eigenvalue (20k trajectories) -- v1.3
- ✓ LsqFit.jl dependency for Levenberg-Marquardt fitting with bounded parameters -- v1.3
- ✓ Dense left+right eigenvector extraction with biorthonormal basis for exact Lindbladian diagnostics (DIAG-01) -- v1.4
- ✓ Lindbladian fixed point computation and trace distance to Gibbs state (DIAG-02) -- v1.4
- ✓ KMS similarity transform anti-Hermitian defect analysis with advisory warning (DIAG-03/04) -- v1.4
- ✓ Observable overlap coefficients via biorthogonal formula c_k = Tr[O R_k] * Tr[L_k^dagger(rho_0 - rho_beta)] (DIAG-05) -- v1.4
- ✓ Delta-Sz symmetry sector labeling with purity fractions and multiplet detection (DIAG-06) -- v1.4
- ✓ Canonical 6-observable set: Z1, X1, Z1_Zhalf, H, Rand_traceless, Mz_stagg (replacing v1.3 8-obs set) -- v1.4
- ✓ run_exact_diagnostics bundle with TrotterDomain support via basis_eigvecs keyword -- v1.4
- ✓ Matrix-free `apply_lindbladian!` for all 4 domains (Energy, Time, Trotter, Bohr) with zero-allocation BLAS.gemm! hot path -- v1.5
- ✓ KrylovWorkspace with precomputed G_left/G_right effective Hamiltonian matrices for optimized matvec -- v1.5
- ✓ KrylovKit-based `krylov_spectral_gap()` API with Lindbladian (:LR) and CPTP channel (:LM) targeting -- v1.5
- ✓ KrylovGapResult struct with eigenvalues, spectral gap, convergence info, matvec count, fixed point -- v1.5
- ✓ Convergence retry with 50% krylovdim increase and pre-flight memory guard -- v1.5
- ✓ Faithful Chen CPTP channel `apply_delta_channel!` with R_total, K0, U_residual precomputation -- v1.5
- ✓ Krylov gap cross-validated vs dense eigen() at n=4 (atol=1e-8) and n=6 (atol=1e-6) for all domains -- v1.5
- ✓ L-vs-E convergence order verification across 3 deltas for all domains -- v1.5
- ✓ Scaling benchmarks at n=3-7 confirming ~4^10 (Energy) and ~4^7.5 (Trotter) power-law with n=10/12 extrapolation -- v1.5
- ✓ G_left/G_right precomputation reducing per-matvec GEMM from 5N to 2+2N -- v1.5

### Active

**Future milestones:**
- [ ] 1D Ising model Hamiltonian generation
- [ ] 2D Heisenberg Hamiltonian generation (lattice graph support)
- [ ] General k-local Hamiltonian construction on arbitrary graphs
- [ ] Hamiltonian simulation time counter (total Hamiltonian simulation cost per Gibbs sample)
- [ ] Gate complexity counter for Trotter-based circuit implementations
- [ ] Qiskit circuit generation for resource estimation (gate count, circuit depth)
- [ ] Documentation: API docs via Documenter.jl, theory tutorials via Literate.jl
- [ ] Paper-ready plotting: convergence curves, sim time plots, gate complexity, mixing time
- [ ] Ding et al. (2024) KMS Lindbladian construction with discrete jump operators (future addition)
- [ ] Jump statistics histogram (empirical rates match theoretical predictions)

### Out of Scope

- GPU acceleration -- not needed for current system sizes (~12 qubits); revisit if scaling beyond 14+ qubits
- Real quantum hardware execution -- circuits are for resource estimation only, not viable for current hardware
- Quantum circuit simulation via Qiskit -- simulation done by QuantumFurnace's own simulators
- MPI distributed computing for single simulations -- single-node multi-core sufficient; cluster used for parameter sweeps via independent Julia processes
- Classical MCMC / classical Gibbs sampling algorithms -- focus is quantum Lindbladian methods only
- Non-Markovian dynamics -- all constructions assume Markovian (Lindblad) framework
- Continuous-time MCWF (adaptive timestep) -- discrete-time CPTP maps per Chen's weak measurement scheme
- Bohr domain trajectories -- diagonalizing Kossakowski matrix prohibitive for >3 qubits; validate Bohr via DM only
- Float32 simulation testing -- type params enable it but testing F32 paths is future work
- n>12 qubit support -- memory and time scale exponentially; 12 qubits is practical limit
- Sparse Lindbladian storage -- full Lindbladian is dense for these constructions; sparsity not exploitable

## Context

**Current State (v1.5 shipped, Phase 32 complete):**
- 8,312 LOC src + 5,071 LOC test (Julia)
- Tech stack: Julia, LinearAlgebra, KrylovKit, FINUFFT, Arpack, LsqFit, BSON, StableRNGs, LibGit2, Dates
- Three gap estimation methods: (1) trajectory-based observable decay fitting, (2) dense eigendecomposition, (3) matrix-free Krylov eigsolve
- Matrix-free Krylov eigsolve for all 4 domains with zero-allocation BLAS hot path and G_left/G_right precomputation
- Cross-validated at n=4 (<1e-8) and n=6 (<1e-6) against dense reference; scaling benchmarks at n=3-7
- EnergyDomain n=10 feasible on cluster (~111h, ~1.3 GB); n=12 infeasible (~12000h); TrotterDomain n=10+ needs on-the-fly NUFFT
- Both DM and trajectory simulation validated and cross-checked
- Multi-threaded trajectory engine with per-thread workspace/RNG, BLAS control, deterministic seeding
- Exact Lindbladian diagnostics: eigendata, fixed point, KMS defect, overlap coefficients, Sz sectors
- Core structs parameterized on `{T<:AbstractFloat}` (Float64 default, Float32-ready)
- API organized: physics building blocks exported, implementation details `_`-prefixed

**Known Limitations:**
- n=6 periodic Heisenberg chain: gap mode has translational + discrete symmetry protection (diagnosed via Sz sector labeling)
- Single-exponential fitting produces non-monotonic error at small delta due to multi-mode contamination
- BENCH-04 partial: total time and matvec count recorded but no isolated per-component timing breakdown
- Memory guard pre-flight estimate underestimates by 28-298x (calibration data available but formula not updated)
- Two-exponential fitting, bootstrap uncertainty, and Richardson extrapolation deferred from v1.4

**Theoretical Foundation:**
The package implements quantum Gibbs samplers from three key papers:
1. **Chen, Kastoryano, Brandao, Gilyen (2023)** "Quantum Thermal State Preparation" -- introduces smoothed Davies generators with approximate GNS detailed balance via Gaussian-weighted operator Fourier transforms. The Lindbladian's fixed point is approximately the Gibbs state. Energy uncertainty ~O(sigma_E) from finite-time Fourier transform.
2. **Chen, Kastoryano, Gilyen (2025)** "An efficient and exact noncommutative quantum Gibbs sampler" -- constructs the first exactly KMS-detailed-balanced Lindbladian with Gaussian transition weight gamma(omega) and a coherent term B that exactly cancels errors. Fixed point is exactly the Gibbs state. Cost: O_tilde(beta) Hamiltonian simulation time per unit Lindbladian evolution.
3. **Ding, Li, Lin (2024)** "Efficient quantum Gibbs samplers with KMS detailed balance condition" -- generalizes KMS construction to discrete jump operators (as few as one), simplifying implementation. Uses Fagnola-Umanita structural characterization. Planned for future integration.

**Target System Sizes:**
- Primary: up to ~12 qubits (density matrix 4096x4096, Liouvillian 16M entries)
- Krylov eigsolve extends practical range: n=10 feasible for EnergyDomain (~111h on cluster)
- Memory driver: precomputed Kraus/NUFFT prefactor matrices, not the density matrix itself
- Trajectories need thousands of samples for statistics -- parallelism over samples is key
- Single-node execution on cluster (up to 512 GB RAM), multi-core parallelism

**Paper Goals:**
Results needed for publication: convergence curves (trace distance vs. steps), mixing time scaling with system size and temperature, comparison across approximation domains (Bohr vs Energy vs Time vs Trotter), trajectory vs density matrix agreement, KMS-vs-GNS comparison data, spectral gap estimation validation, Krylov scaling benchmarks.

## Constraints

- **Language**: Julia -- non-negotiable, leverages Julia's multiple dispatch for domain hierarchy and performance
- **Correctness**: All Lindbladian constructions must be mathematically faithful to the source papers; detailed balance properties must hold to machine precision (for exact KMS) or documented approximation bounds (for GNS and for implemented KMS due to Trotterization or quadrature errors.)
- **System size**: Practical limit ~12-14 qubits due to exponential Hilbert space scaling (density matrix 2^n x 2^n)
- **Dependencies**: Minimize external dependencies; rely on Julia stdlib + established numerical packages (Arpack, FINUFFT, KrylovKit, LsqFit)
- **Dual purpose**: Code must be both performant (for research results) and readable (for community pedagogy) -- well-documented functions with mathematical cross-references to papers

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| Four-domain approximation hierarchy (Bohr/Energy/Time/Trotter) | Maps directly to the theoretical approximation levels in Chen 2023/2025; enables comparing accuracy vs. cost | ✓ Good -- error hierarchy verified empirically (DMTST-02) |
| KMS detailed balance as primary construction | Exact stationarity (Chen 2025) eliminates accumulated trajectory error vs approximate GNS | ✓ Good -- BohrDomain fixed point at 1.6e-15, TrotterDomain at ~9e-9 |
| Density matrix + trajectory dual simulation | DM gives exact dynamics for small systems; trajectories scale better and model the actual quantum algorithm | ✓ Good -- cross-validation confirms agreement; DM serves as ground truth |
| Per-operator Lie-Trotter splitting for trajectories | Matches paper's algorithm; per-operator branching gives correct delta_eff = delta * N_jumps | ✓ Good -- fixed TrotterDomain Gibbs distance from 0.004 to 9e-9 |
| Eigendecomposition instead of Cholesky for PSD guard | Silent clamp of negative eigenvalues; more robust than Cholesky fallback | ✓ Good -- handles non-PSD S matrices gracefully |
| DM-based regression instead of frozen BSON | Frozen trajectories not portable across Julia versions/platforms (BLAS differences) | ✓ Good -- DM comparison with delta=0.01 stays within atol=0.05 |
| Type parameterization on {T<:AbstractFloat} | Enables Float32 paths for future performance experiments without changing call sites | ✓ Good -- all Float64 tests pass unchanged, Float32 ready |
| Domain dispatch via config type parameter | AbstractConfig{D,T} replaces separate domain argument; cleaner multiple dispatch | ✓ Good -- eliminated redundant domain arguments across codebase |
| Immutable TrottTrott + fully-initialized HamHam | Eliminates two-step init patterns and mutable state bugs | ✓ Good -- cleaner construction, BSON backward compat preserved |
| Internal functions _-prefixed (not unexported) | Visible in qualified access for tests; clear public/internal boundary | ✓ Good -- 18 internals prefixed, tests use QuantumFurnace._name |
| Index-based sparse accumulation in B_bohr | Avoid per-iteration sparse matrix allocation in hot loop | ✓ Good -- allocation regression tests confirm zero sparse alloc |
| Read-only TrajectoryFramework + explicit workspace/RNG | Enables thread-safe concurrent trajectory execution from shared framework | ✓ Good -- workspace independence verified, threading works correctly |
| Per-trajectory Xoshiro(seed + traj_id) seeding | Deterministic reproducibility independent of thread count; same results serial vs threaded | ✓ Good -- bitwise identical results across runs |
| BLAS.set_num_threads(1) during threaded execution | Prevents thread oversubscription from BLAS competing with Julia threads | ✓ Good -- restored in try/finally, performance tests confirm speedup |
| Concrete-typed TrajectoryFramework fields (F,P params) | Zero-allocation hot path via function barrier; avoids abstract type access in inner loop | ✓ Good -- allocation regression tests confirm zero alloc in step |
| Separate run_trajectories_convergence (not modifying run_trajectories) | Avoids API bloat on core function; convergence is a distinct use case | ✓ Good -- clean separation, shared _run_batch_no_obs! backend |
| Dict-based BSON serialization for ExperimentResult | Avoids parametric struct pitfalls in BSON; Dict round-trips reliably | ✓ Good -- all round-trip tests pass |
| Domain-aware JumpOp basis at construction (trotter.eigvecs for TrotterDomain) | Eliminates redundant transform_jumps_to_basis calls downstream | ✓ Good -- one transform at source, flatter call chain |
| LindbladianResult/DMSimulationResult replacing HotSpectralResults/HotAlgorithmResults | Simpler 3-4 field structs without hamiltonian/config/trotter baggage | ✓ Good -- cleaner API, less coupling |
| Single-node multi-core for trajectories | Shared memory for precomputed data avoids serialization overhead; cluster nodes have enough RAM | ✓ Good -- multi-threaded engine operational |
| LsqFit.jl for exponential fitting (v1.3) | Levenberg-Marquardt with parameter bounds and covariance CIs; standard Julia choice | ✓ Good -- recovers synthetic decay rates within CIs; SingularException handled gracefully |
| Single-exponential model A*exp(-gap*t)+C (v1.3) | Simplest model capturing dominant decay; multi-exponential deferred | ⚠️ Revisit -- works for n=4 (<1% error) but non-monotonic at small delta due to multi-mode contamination |
| Smallest-gap selection for best observable (v1.3) | Gap mode is the slowest decaying; smallest fitted gap closest to true gap | ✓ Good -- reduces n=4 from ~1.6x to ~1.17x factor; Quick-23 validation |
| Consolidated single observable builder (v1.3 Phase 25) | build_preset_trajectory_observables replaces 4 separate builders | ✓ Good -- cleaner API, single entry point, 8-observable bundle |
| CrossValidationResult removed (v1.3 Phase 25) | Thin wrapper over manual comparison; eigenbasis_overlap_analysis more useful | ✓ Good -- simpler API surface, overlap analysis provides more insight |
| abs(real(spectral_gap)) for exact gap (v1.3) | Complex eigenvalues: real part gives decay rate; imaginary part gives oscillation frequency | ✓ Good -- locked decision, enforced throughout |
| Dense eigen() for exact diagnostics (v1.4) | Arpack cannot compute left eigenvectors needed for biorthogonal overlap formula | ✓ Good -- handles n<=6 (4096x4096) in seconds |
| Biorthogonal overlap formula c_k = Tr[O R_k] * Tr[L_k^dagger(rho_0 - rho_beta)] (v1.4) | Proper steady-state subtraction using exact left+right eigenvectors | ✓ Good -- explains zero-overlap mystery at n=6 |
| Canonical 6-observable set replacing 8-obs bundle (v1.4) | Z1, X1, Z1_Zhalf, H, Rand_traceless, Mz_stagg -- physically motivated from supplementary info | ✓ Good -- cleaner set, random traceless adds coverage |
| Advisory-only defect warning at 0.1 threshold (v1.4) | Defect ratio is informational, does not gate computation | ✓ Good -- non-blocking diagnostics |
| Left eigenvectors via transpose(inv(V_right)) (v1.4) | Transpose (not conjugate transpose) for proper biorthogonality Tr[L_k R_j] = delta_kj | ✓ Good -- validated against known systems |
| KrylovWorkspace{T,PD} dual type params (v1.5) | PD captures NamedTuple type for zero-overhead access to precomputed data | ✓ Good -- eliminates abstract field boxing |
| BLAS.gemm!/axpy! for zero-allocation matvec (v1.5) | mul! and broadcast allocate; raw BLAS calls do not | ✓ Good -- verified 0 bytes in all domains |
| kron(A,B)vec(X)=vec(B*X*A^T) convention (v1.5) | Consistent vectorization convention matching dense Lindbladian construction | ✓ Good -- resolved complex jump operator discrepancy (quick-35) |
| KrylovKit Arnoldi for non-Hermitian Lindbladian (v1.5) | Lanczos only for Hermitian; Arnoldi handles general matrices | ✓ Good -- converges reliably for all domains |
| Faithful Chen CPTP channel (v1.5) | Euler approximation E=I+delta*L had O(delta^2) error; faithful channel matches run_thermalization | ✓ Good -- channel eigsolve accuracy improved (quick-36) |
| G_left/G_right precomputation (v1.5) | Reduces per-matvec GEMM count from 5N to 2+2N by precomputing aggregate effective Hamiltonian | ✓ Good -- all correctness tests pass, 309 LOC dead code removed |
| 4 separate G matrices for BohrDomain adjoint (v1.5) | R_total is non-Hermitian for Bohr, so G_left_adj != G_right^T | ✓ Good -- correct adjoint verified via duality test |
| EnergyDomain scaling ~4^10, TrotterDomain ~4^7.5 (v1.5) | Empirical power-law fit at n=3-7; n=10 feasible for Energy, not for Trotter | ✓ Good -- calibrated extrapolation |
| Qiskit for circuit generation (Python interop) | Qiskit is the standard for quantum circuit representation; Julia quantum circuit ecosystem less mature | -- Pending |

---
*Last updated: 2026-02-25 after v1.5 milestone*
