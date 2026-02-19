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

