# Milestones: QuantumFurnace.jl

## v1.0 Trajectories (Shipped: 2026-02-14)

**Started:** 2026-02-13 | **Shipped:** 2026-02-14
**Phases:** 1-5 (10 plans + 12 quick tasks) | **Tests:** 0 → 224 | **Commits:** 85
**Julia LOC:** 13,082 total (2,034 test) | **Files changed:** 18 (+1,793 / -243)
**Git range:** ecd17ae..5dc2f45

**Delivered:** Trajectory simulation fixed, validated against density matrix evolution, and locked down with a comprehensive 224-test correctness suite covering CPTP channels, detailed balance, error scaling, cross-validation, and regression.

**Key accomplishments:**
1. Fixed trajectory compilation (5 bugs) and the critical U_B ordering bug in EnergyDomain
2. CPTP channel verification: K0†K0 + δR + U_res†U_res = I to 1e-10 for all domains
3. DM ground truth: BohrDomain detailed balance at 1.6e-15, domain error hierarchy verified, O(δ²)/O(δ) scaling confirmed
4. Trajectory-vs-DM cross-validation: single-step O(δ²) scaling and multi-step Gibbs convergence for Energy, Time, Trotter
5. Per-operator Lie-Trotter splitting refactor: TrotterDomain Gibbs fixed point distance 0.004 → 9e-9
6. Portable regression framework: DM-based comparison replacing frozen BSON for cross-platform/version stability

**Archives:** [v1.0-ROADMAP.md](milestones/v1.0-ROADMAP.md) | [v1.0-REQUIREMENTS.md](milestones/v1.0-REQUIREMENTS.md)

---


## v1.1 Reduce (Shipped: 2026-02-15)

**Started:** 2026-02-15 | **Shipped:** 2026-02-15
**Phases:** 6-11 (16 plans + 5 quick tasks) | **Tests:** 224 → 231 | **Commits:** 93
**Julia LOC:** 4,422 src + 2,222 test | **Files changed:** 87 (+9,420 / -2,442)
**Git range:** chore(06-01)..docs(11-03)

**Delivered:** Codebase refactored and simplified — pruned ~987 lines of dead code and ~35 unused functions, extracted 4 DRY helpers, simplified core structs with type parameterization on `{T<:AbstractFloat}`, cleaned public API surface, and eliminated hot-path heap allocations. All 231 tests passing.

**Key accomplishments:**
1. Pruned ~987 lines of dead commented code and ~35 unused functions/3 dead structs
2. Extracted 4 DRY helpers: `hermitianize!`, `transform_jumps_to_basis`, `apply_cptp_channel!`, `apply_coherent_unitary!`
3. Simplified core structs: immutable `TrottTrott`, fully-initialized `HamHam`, deduplicated config constructors
4. Parameterized all core types on element type `{T<:AbstractFloat}` enabling future Float32 paths
5. Cleaned API surface: organized exports into labeled groups, internalized ~18 implementation details with `_` prefix
6. Eliminated hot-path allocations: index-based sparse accumulation, in-place phase rotations, precomputed basis transforms

**Archives:** [v1.1-ROADMAP.md](milestones/v1.1-ROADMAP.md) | [v1.1-REQUIREMENTS.md](milestones/v1.1-REQUIREMENTS.md)

---


## v1.2 Multi-threading (Shipped: 2026-02-16)

**Started:** 2026-02-15 | **Shipped:** 2026-02-16
**Phases:** 12-19 (15 plans + 3 quick tasks) | **Tests:** 231 → 539 | **Commits:** ~60
**Julia LOC:** 5,479 src + 3,796 test | **Timeline:** 2 days

**Delivered:** Multi-threaded trajectory sampling engine with GNS comparison path, adaptive convergence-driven batching, BSON experiment serialization, and paper-ready KMS-vs-GNS parameter sweep experiments. All 539 tests passing.

**Key accomplishments:**
1. Thread-safe workspace separation (read-only TrajectoryFramework, explicit workspace/RNG) enabling concurrent trajectory execution
2. Multi-threaded trajectory engine with per-thread Xoshiro seeding, BLAS thread control, and zero-allocation hot path via concrete-typed framework fields
3. GNS (approximate detailed balance) trajectory path verified end-to-end with convergence toward GNS fixed point
4. BSON-based ExperimentResult serialization with full metadata capture (git hash, timestamp, seed, thread count)
5. Convergence tracking: batch-level trace distance to Gibbs + per-observable (ZZ correlations, energy) monitoring
6. Adaptive sampling with automatic stopping (relative change <1% for 3 consecutive batches) and hard trajectory cap
7. KMS-vs-GNS parameter sweep (n=4 x beta=5,10,20 x {KMS, GNS@1/beta, GNS@0.5/beta}) confirming KMS achieves lower trace distance
8. Post-milestone cleanup: flattened 5-level call chain to 3, eliminated redundant basis transforms, simplified result structs (LindbladianResult, DMSimulationResult)

**Archives:** [v1.2-ROADMAP.md](milestones/v1.2-ROADMAP.md) | [v1.2-REQUIREMENTS.md](milestones/v1.2-REQUIREMENTS.md)

---


## v1.3 Mixing Time Estimation (Shipped: 2026-02-19)

**Started:** 2026-02-17 | **Shipped:** 2026-02-19
**Phases:** 20-25 (10 plans + 11 quick tasks) | **Tests:** 539 → 666 | **Commits:** 98 since v1.2
**Julia LOC:** 6,274 src + 4,366 test | **Files changed:** 93 (+16,804 / -2,329)
**Git range:** feat(20-01)..docs(quick-32)

**Delivered:** Spectral gap estimation from trajectory-based observable decay, with eigenbasis overlap diagnostics, cross-validated against exact Liouvillian eigenvalues. n=4 Heisenberg chain achieves 0.72% accuracy at 20k trajectories. Includes observable-only trajectory runner, exponential fitting, and unified validation pipeline.

**Key accomplishments:**
1. Single-call `estimate_spectral_gap` API orchestrating observable construction, trajectory simulation, exponential fitting, and best-observable selection
2. 8-observable eigenbasis overlap diagnostic (`eigenbasis_overlap_analysis`) decomposing observables into Lindbladian eigenmodes to identify gap-mode coupling
3. Observable-only trajectory runner (`run_observable_trajectories`) measuring time-resolved `<O>(t)` without per-trajectory DM reconstruction, bitwise cross-validated
4. Exponential decay fitting (`fit_exponential_decay`) via LsqFit.jl with auto log-linear initial guess, confidence intervals, and R-squared quality metrics
5. n=4 spectral gap validated to 0.72% accuracy (20k trajectories, beta=10); n=6 zero-overlap physics limitation diagnosed and documented
6. Systematic investigation (Quick-22 through Quick-32) of gap estimation accuracy: delta_eff fix, observable selection improvements, symmetry analysis, and confirmation that trajectory simulation is correct

**Archives:** [v1.3-ROADMAP.md](milestones/v1.3-ROADMAP.md) | [v1.3-REQUIREMENTS.md](milestones/v1.3-REQUIREMENTS.md)

---


## v1.4 Spectral Gap Refinement (Partial — Shipped: 2026-02-20)

**Started:** 2026-02-19 | **Shipped:** 2026-02-20 (partial — Phase 26 only, Phases 27-30 deferred)
**Phases:** 26 (1 of 5 planned) | **Plans:** 2 + 1 quick task | **Commits:** 21 since v1.3
**Julia LOC:** 6,869 src + 3,931 test | **Files changed:** 38 (+5,882 / -1,947)
**Git range:** feat(26-01)..035cfc4

**Delivered:** Exact Lindbladian diagnostics infrastructure — dense eigendecomposition with biorthonormal left+right eigenvectors, KMS similarity transform defect analysis, observable overlap coefficients via biorthogonal formula, Delta-Sz symmetry sector labeling, and canonical 6-observable set. TrotterDomain support added via basis_eigvecs keyword.

**Key accomplishments:**
1. Dense left+right eigenvector extraction enabling biorthogonal overlap formula c_k = Tr[O R_k] * Tr[L_k^dagger(rho_0 - rho_beta)]
2. Lindbladian fixed point computation with trace distance to Gibbs state verification
3. KMS diagonal similarity transform for anti-Hermitian defect diagnosis (advisory-only at 0.1 threshold)
4. Delta-Sz symmetry sector labeling explaining n=6 zero-overlap gap mode (translational + discrete symmetry protection)
5. Canonical 6-observable set (Z1, X1, Z1_Zhalf, H, Rand_traceless, Mz_stagg) replacing v1.3 8-obs bundle
6. run_exact_diagnostics bundle with TrotterDomain support via basis_eigvecs keyword

**Deferred to future milestone:** Two-exponential fitting (FIT), effective rate plots (RATE), bootstrap uncertainty (BOOT), Richardson extrapolation (RICH), validation dashboard (VAL) — 17 requirements

**Archives:** [v1.4-ROADMAP.md](milestones/v1.4-ROADMAP.md) | [v1.4-REQUIREMENTS.md](milestones/v1.4-REQUIREMENTS.md)

---


## v1.5 Krylov Gap Estimation (Shipped: 2026-02-25)

**Started:** 2026-02-20 | **Shipped:** 2026-02-25
**Phases:** 27-32 (12 plans + 3 quick tasks) | **Commits:** ~72
**Julia LOC:** 8,312 src + 5,071 test | **Files changed:** 13 (+2,920 / -423)
**Git range:** feat(27-01)..docs(quick-37)

**Delivered:** Matrix-free spectral gap estimation via KrylovKit.jl, enabling gap computation for system sizes (up to ~12 qubits) where the full Lindbladian matrix cannot be stored. Zero-allocation matvec for all 4 domains, cross-validated against dense eigen() at n=4 and n=6, with empirical scaling benchmarks at n=3-7.

**Key accomplishments:**
1. Matrix-free `apply_lindbladian!` for all 4 domains (Energy, Time, Trotter, Bohr) with zero-allocation BLAS.gemm! hot path, validated to <1e-12 against dense construction
2. KrylovKit-based `krylov_spectral_gap()` API with Lindbladian (:LR) and CPTP channel (:LM) targeting, convergence retry logic, and pre-flight memory guard
3. Cross-validated Krylov gap vs dense eigen() at n=4 (atol=1e-8) and n=6 (atol=1e-6) for all domains; L-vs-E consistency verified
4. Empirical scaling benchmarks at n=3-7 confirming ~4^10 (Energy) and ~4^7.5 (Trotter) power-law scaling with n=10/12 extrapolation
5. G_left/G_right precomputation reducing per-matvec GEMM count from 5N to 2+2N; 309 lines of dead code removed

**Tech debt (from audit):**
- BENCH-04 partial: per-matvec timing recorded but no isolated BLAS/precompute/Krylov overhead split
- Memory guard pre-flight estimate underestimates by 28-298x (calibration data available)
- Stale BohrDomain adjoint docstring (documentation only)
- Julia runtime tests not executed in sandbox (structural verification only)

**Archives:** [v1.5-ROADMAP.md](milestones/v1.5-ROADMAP.md) | [v1.5-REQUIREMENTS.md](milestones/v1.5-REQUIREMENTS.md)

---


## v2.0 Restructure (Shipped: 2026-02-28)

**Started:** 2026-02-25 | **Shipped:** 2026-02-28
**Phases:** 33-38 (19 plans + 2 quick tasks) | **Tests:** ~1199 | **Commits:** 91
**Julia LOC:** 8,299 src + 4,559 test | **Files changed:** 139 (+16,664 / -4,572)
**Git range:** `docs(33): capture phase context`..`docs(38-05): complete Krylov test @info`

**Delivered:** Major codebase restructure — unified Config{S,D,C,T} type hierarchy replacing 4 separate config types, eliminated code duplication (domain_prefactor, unified oft!, shared CPTP channel helper), consolidated 6 workspace types into 2, defined 4 clean run_* entry points with typed Result structs and BSON save/load, reorganized exports and removed dead code, and instrumented all tests with @info output and threshold rationale. Architecture now ready for DLL construction, error estimation, and gate complexity features.

**Key accomplishments:**
1. Unified Config{S,D,C,T} type hierarchy migrating 37+ files from 4 separate config types, enabling future DLL construction extensibility
2. Code deduplication: domain_prefactor() replacing 16 formula copies, unified oft!(), consolidated sandwich helpers, extracted _build_cptp_channel
3. Workspace consolidation: merged KrylovWorkspace + KrausScratch + LindbladianWorkspace into Workspace{S,D,C,T}; flattened trajectory per-operator Kraus data
4. Clean API: 4 run_* entry points (run_lindblad, run_thermalize, run_krylov_spectrum, run_trajectory) with typed Result structs and BSON save/load
5. File organization: dead code removed (@distributed, SharedArrays), staging directory created, exports reorganized by simulation type
6. Test infrastructure: consolidated helpers (make_config/make_test_system), 204+ @info outputs, 163 threshold rationale comments across all test files

**Tech debt (from audit):**
- DEDUP-02 deferred by user (keep explicit for-loops)
- ORG-01 descoped by user (file renaming handled manually)
- Phase 35 missing formal VERIFICATION.md (UAT conducted, all issues resolved)
- run_* entry points lack direct automated test coverage (exercised via simulation scripts and result round-trip tests)

**Archives:** [v2.0-ROADMAP.md](milestones/v2.0-ROADMAP.md) | [v2.0-REQUIREMENTS.md](milestones/v2.0-REQUIREMENTS.md)

---


## v2.1 Speedup & Mixing Time (Shipped: 2026-03-04)

**Started:** 2026-03-01 | **Shipped:** 2026-03-04
**Phases:** 39-43 (8 plans, 19 tasks) | **Tests:** 1141 → 1273 | **Commits:** ~35
**Julia LOC:** 10,242 src + 6,316 test | **Files changed:** 50 (+10,321 / -2,249)
**Git range:** b5ad26f..c415d6f

**Delivered:** Performance optimization and mixing time estimation -- per-jump CPTP channel precomputation eliminates eigendecomposition from the DM hot loop, multi-threaded BLAS and omega-loop parallelism for DM thermalization, save_every observation gating, single-exponential and bi-exponential fitting with quality gates and extrapolation support. Bi-exponential model reduces extrapolation error from ~26% to <0.001%.

**Key accomplishments:**
1. Per-jump CPTP channel precomputation eliminates eigendecomposition from run_thermalize hot loop across all 4 domains (Energy, Time, Trotter, Bohr)
2. save_every keyword for trace distance computation frequency with full backward compatibility (default save_every=1)
3. Multi-threaded BLAS + omega-loop parallelism for DM thermalization with per-task ThermalizeScratch isolation and BLAS try/finally save/restore
4. Mixing time estimation API via exponential fit on trace distance curve with quality gates (R^2, offset), extrapolation support, and MixingTimeEstimate struct
5. Bi-exponential decay fitting (BiexpFitResult, fit_biexponential_decay, model=:biexp) with Roots.Bisection extrapolation, reducing error from ~26% to <0.001%

**Tech debt (from audit):**
- src/staging/fitting.jl: dead code leftover from promotion to src/fitting.jl (no runtime impact)
- test/test_threading.jl:31: `@test true` placeholder for nthreads==1 guard (standard Julia pattern)

**Archives:** [v2.1-ROADMAP.md](milestones/v2.1-ROADMAP.md) | [v2.1-REQUIREMENTS.md](milestones/v2.1-REQUIREMENTS.md)

---

