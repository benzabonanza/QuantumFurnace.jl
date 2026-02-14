# Project Research Summary

**Project:** QuantumFurnace.jl v1.0 Trajectories Milestone
**Domain:** Quantum trajectory simulation validation and correctness testing
**Researched:** 2026-02-13
**Confidence:** HIGH (direct codebase analysis + verified quantum simulation literature)

## Executive Summary

The v1.0 Trajectories milestone is about fixing trajectory sampling bugs and building a comprehensive test suite to validate quantum Monte Carlo wave function (MCWF) simulation against density matrix (DM) evolution. Research reveals that the current trajectory implementation has three specific compilation bugs in `build_trajectoryframework` and a coherent term ordering issue in `step_along_trajectory!`. The "flat vs two-stage sampling" question is a red herring for the uniform Heisenberg case — the real bugs are code errors that prevent execution.

The recommended validation strategy follows a strict hierarchy: CPTP verification must come first (before any trajectory runs), then DM-only tests establish reference values, and only then can trajectory-vs-DM cross-validation proceed. The error hierarchy (Bohr → Energy → Time → Trotter → delta-stepping → trajectory sampling) provides natural test tiers. The DM simulation (`run_thermalization`) is the ground truth oracle — no external validation libraries needed.

Critical insight: compare trajectory-averaged rho against DM rho, not against Gibbs. Both share the same discretization errors (quadrature, Trotter, delta-step); the only difference should be statistical sampling noise O(1/sqrt(N_traj)). Testing against Gibbs conflates multiple error sources and leads to false failures or false confidence. The stack needs only 4 test-only dependencies (StableRNGs, HypothesisTests, StatsBase, Aqua) — no production dependencies required.

## Key Findings

### Recommended Stack

**No new production dependencies needed.** The trajectory fix is code correction, not dependency addition. All numerical tools (LinearAlgebra, Random, Statistics, BLAS) are already present.

**Four test-only dependencies recommended:**
- **StableRNGs**: Reproducible RNG for regression tests — Julia's default RNG is not version-stable
- **HypothesisTests**: Statistical significance tests (OneSampleTTest) for trajectory vs DM comparison
- **StatsBase**: Standard error of mean, distribution distance metrics (L2dist, kldivergence)
- **Aqua**: Automated package quality checks (ambiguities, stale deps, method piracy)

These go in `[extras]` section of Project.toml, not `[deps]`. They do not increase dependency burden for users.

**What NOT to add:**
- OrdinaryDiffEq — overkill; DM simulation already provides reference
- ExponentialUtilities — trajectory code uses discrete Kraus maps, not matrix exponentials
- KrylovKit — eigendecomposition is a separate concern, defer to later
- QuantumOptics.jl / QuantumToolbox.jl — third-party frameworks with different type systems; DM simulation is the reference
- Distributions.jl — heavyweight (~15 deps); HypothesisTests already provides needed statistical tests

### Expected Features

**Must have (table stakes for milestone completion):**
- Fix `build_trajectoryframework` compilation bugs — `trotter` variable undefined, `B_total` used when unassigned
- CPTP verification test — verify K0†K0 + delta*R + U_res†U_res = I
- Trajectory-averaged rho vs DM evolution — primary correctness criterion for three domains (Energy, Time, Trotter)
- Domain error hierarchy test — Bohr < Energy < Time < Trotter
- Single-step delta^2 error scaling — empirical verification of Chen Theorem III.1
- Multi-step delta error scaling — O(delta) accumulation over full evolution
- Detailed balance test — BohrDomain fixed point = Gibbs to machine precision
- Coherent term correctness — with_coherent=true dramatically closer to Gibbs than false
- Statistical 1/sqrt(N) convergence test — trajectory averaging error decreases as expected

**Should have (strengthens validation):**
- OFT consistency across domains — Energy OFT matches Time/Trotter NUFFT OFT
- R matrix cross-validation — trajectory R matches Liouvillian R
- Per-observable trajectory convergence — individual expectation values match DM
- Jump statistics histogram — empirical jump rates match theoretical predictions
- Regression tests with frozen data — lock down reference results

**Defer (v2+):**
- Diamond norm channel comparison — expensive (SDP), only needed for publication rigor
- Multi-threaded trajectory validation — correctness first, performance later
- Non-Hermitian jump operator tests — all target systems use Hermitian jumps
- Large system benchmarks (>6 qubits) — validation is about correctness, not scale
- Bohr domain trajectories — diagonalizing Kossakowski matrix is prohibitive

### Architecture Approach

The validation architecture uses a **reference oracle pattern**: each approximation layer (Energy, Time, Trotter) is tested against the next-higher-fidelity layer. DM simulation is known-good and working; trajectory code is the new component being validated. Tests follow a six-tier structure that mirrors the error hierarchy.

**Major components and their roles:**

1. **Test helpers** (`test/test_helpers.jl`) — Shared fixtures: `make_test_system()` creates standard Heisenberg systems, `make_liouv_config()` / `make_therm_config()` generate configs, tiered tolerance constants match error levels (TOL_BOHR_EXACT, TOL_ENERGY_QUAD, etc.)

2. **DM reference tests** (Tiers 1-4) — Validate DM path first: detailed balance (Gibbs fixed point), coherent term consistency (B_bohr vs B_time vs B_trotter), OFT consistency (exact vs NUFFT), DM step error scaling (delta^2 per step, delta total). These establish the ground truth before any trajectory code runs.

3. **Trajectory core** (`trajectories.jl`) — Needs three bug fixes: (a) remove `where T` and `trotter=trotter` kwarg from `build_trajectoryframework`, (b) initialize `B_total = nothing` before conditional, (c) move U_B application to after branch selection (matches DM ordering). The two-stage sampling question is resolved: flat scan is mathematically correct for uniform jump operators, which is the current case.

4. **Trajectory validation tests** (Tiers 5-6) — After fixes: trajectory average matches DM (primary test), trajectory converges to Gibbs with correct error budget, error hierarchy preserved in trajectories, 1/sqrt(N) convergence verified.

**Critical architectural insight:** `run_thermalization` uses randomized jump selection (pick one jump per step, rescale by 1/p_jump). Trajectory code uses deterministic full map (all jumps contribute to R). These are different CPTP channels. For validation, compare both to the exact Lindbladian channel or create a `run_thermalization_deterministic` that matches the trajectory channel structure.

### Critical Pitfalls

1. **Comparing trajectory to Gibbs instead of to DM** — The error budget has 6+ sources (sampling, delta-step, quadrature, Trotter, coherent approx). DM shares all except sampling. Test `||rho_traj - rho_DM|| < C/sqrt(N_traj)` isolates trajectory correctness. Testing against Gibbs conflates errors and leads to false failures or false confidence. **Mitigation:** Primary test is always trajectory-vs-DM; secondary test compares both to Gibbs separately.

2. **Wrong delta-step error scaling expectation** — Chen Theorem III.1: O(delta^2) per step, O(delta) total over T/delta steps. Developers may expect either O(delta) per step (too pessimistic) or O(delta^2) total (too optimistic). **Mitigation:** Delta-scaling test verifies empirical O(delta) total error. Document "per step: O(delta^2), total: O(delta)" in all relevant tests.

3. **Stochastic test flakiness from fixed seeds hiding bugs** — Fixed seed may avoid code paths where bugs exist (e.g., fallback branch never triggered). **Mitigation:** Two-tier testing: (1) deterministic property tests (CPTP, normalization), (2) statistical convergence tests (many trajectories, no fixed seed, 3-sigma bounds).

4. **Probability normalization sum drift masked by fallback** — `p_nojump + p_res + p_jump_total` should equal 1.0. If `base_prefactor` formula differs between `precompute_R` and `step_along_trajectory!`, cumulative scan fails to reach target or overshoots, triggering silent fallback. **Mitigation:** Add debug assertion `|csum - p_jump_total| / p_jump_total < 1e-10`, extract `base_prefactor` to shared function.

5. **Coherent unitary applied to wrong state** — Current code applies U_B after computing branch probabilities but before selection, then branches use pre-U_B cached quantities (K0*psi, R*psi). Correct order: apply Kraus map first, then U_B. **Mitigation:** Move U_B application to after branch selection, matching DM code ordering (U_B applied to entire Kraus output).

## Implications for Roadmap

Based on research, the milestone requires **three sequential phases** with clear dependencies:

### Phase 1: Fix Trajectory Code and Add CPTP Guards
**Rationale:** Nothing trajectory-related works until compilation errors are fixed. CPTP verification must come before any trajectory runs — without it, non-unit total probability will bias all downstream results.

**Delivers:**
- `build_trajectoryframework` bugs fixed (undefined `trotter`, uninitialized `B_total`, remove `where T`)
- Coherent unitary ordering corrected (move U_B after branch selection)
- S matrix PSD check added (`issuccess(cholesky)` or eigenvalue guard)
- Probability normalization assertion added to `step_along_trajectory!`

**Addresses features:**
- build_trajectoryframework bug fix (must-have P0)
- Coherent term correctness (must-have)
- CPTP verification test (must-have P0)

**Avoids pitfalls:**
- Pitfall 5 (coherent unitary ordering)
- Pitfall 6 (normalization drift)
- Pitfall 9 (Cholesky non-PSD silent failure)

**Estimated effort:** LOW-MEDIUM (4 specific bugs, well-characterized)

### Phase 2: Build DM-Only Test Suite (Reference Establishment)
**Rationale:** The DM simulation path is known-good and working. Establishing reference values from DM tests first enables trajectory validation in Phase 3. If DM tests fail, the bug is in DM code (not trajectories), saving debugging time.

**Delivers:**
- `test/runtests.jl` — central test runner
- `test/test_helpers.jl` — shared fixtures (make_test_system, configs, tolerances)
- `test/test_detailed_balance.jl` — Gibbs fixed point for Bohr/Energy/Time/Trotter
- `test/test_coherent_term.jl` — B operator cross-domain consistency
- `test/test_oft_consistency.jl` — OFT consistency (exact vs NUFFT)
- `test/test_dm_step_errors.jl` — delta^2 per step, delta total scaling

**Addresses features:**
- Detailed balance test (must-have)
- Domain error hierarchy test (must-have)
- Single/multi-step DM error scaling (must-have)
- OFT consistency (should-have)
- R matrix cross-validation (should-have)

**Uses stack:**
- Test stdlib (already present)
- StableRNGs (for any deterministic regression tests)
- HypothesisTests (for statistical comparisons if needed)

**Avoids pitfalls:**
- Pitfall 2 (error budget confusion) — establishes DM-only error floor
- Pitfall 3 (wrong scaling expectation) — verifies empirical O(delta) total
- Pitfall 4 (seed flakiness) — uses two-tier pattern

**Estimated effort:** MEDIUM (6 test files, ~300-500 LOC total, well-structured)

### Phase 3: Trajectory Validation and Cross-Validation
**Rationale:** After DM reference is established (Phase 2) and trajectory bugs are fixed (Phase 1), this phase validates that trajectory-averaged rho matches DM rho within statistical noise. This is the milestone success criterion.

**Delivers:**
- `test/test_trajectory.jl` — trajectory average matches DM for Energy/Time/Trotter
- `test/test_traj_vs_dm.jl` — error hierarchy preserved in trajectories
- Statistical convergence verification (1/sqrt(N) scaling)
- Per-observable convergence tests (should-have)
- Jump statistics histogram (should-have)
- Regression tests with frozen data (should-have)

**Addresses features:**
- Trajectory-averaged rho vs DM evolution (must-have, CRITICAL)
- Statistical 1/sqrt(N) convergence test (must-have)
- Per-observable trajectory convergence (should-have)
- Jump statistics histogram (should-have)
- Regression tests (should-have)

**Uses stack:**
- StableRNGs (reproducible seeds for regression tests)
- HypothesisTests (OneSampleTTest for trajectory mean vs DM value)
- StatsBase (sem for standard error, distribution distances)

**Avoids pitfalls:**
- Pitfall 1 (flat vs two-stage) — research resolved: flat scan is correct for uniform jumps
- Pitfall 2 (error budget) — compares trajectory-vs-DM, not trajectory-vs-Gibbs
- Pitfall 8 (insufficient trajectories) — power analysis determines N_traj (500-2000 depending on system size)

**Estimated effort:** MEDIUM-HIGH (~400 LOC, requires multiple trajectory runs for convergence tests)

### Phase Ordering Rationale

**Sequential dependencies:**
- Phase 2 and Phase 3 both depend on Phase 1 (trajectory code must compile and run)
- Phase 3 depends on Phase 2 (DM reference values needed for comparison)
- Within Phase 2, tests are independent (can be written in any order)
- Within Phase 3, basic trajectory test must pass before statistical convergence tests

**Why this grouping:**
- Phase 1 is pure code fixes (no tests, but enables everything)
- Phase 2 validates the known-good path (DM) and establishes oracle
- Phase 3 validates the new path (trajectories) against oracle

**Parallelization opportunities:**
- Phase 2 test files can be written in parallel (they share only test_helpers.jl)
- Phase 3 cannot start until Phase 1 and Phase 2 complete

**Cleanup after validation:**
- Delete ~200 lines of commented code in jump_workers.jl, trajectories.jl
- Remove old KrausFramework references
- Project.toml cleanup (move Revise, Debugger, BenchmarkTools to proper sections) — can happen during Phase 2

### Research Flags

**No phases need deeper research during planning.** This is a bug-fix and validation milestone for existing code, not new feature development. All research is complete.

**Standard patterns (research done):**
- CPTP channel validation — well-understood quantum information pattern
- Monte Carlo convergence testing — standard statistical pattern (1/sqrt(N))
- Error hierarchy testing — reference oracle pattern, well-documented
- Julia Test.jl patterns — standard Julia testing practices

**Phase-specific notes:**
- Phase 1: All four bugs are characterized (file, line number, fix known)
- Phase 2: Test structure follows QuantumOptics.jl and Julia community standards
- Phase 3: Validation methodology follows MCWF literature (QDYN, QuantumOptics, QuantumToolbox)

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Stack | HIGH | 4 test-only deps verified via official docs and GitHub repos; all are JuliaStats/JuliaTesting org packages with stable APIs |
| Features | HIGH | Codebase analysis identified specific bugs; validation methodology verified across QuantumOptics.jl, QuantumToolbox.jl, QDYN docs |
| Architecture | HIGH | Direct analysis of all 23 source files, 7 test files, simulation scripts; error hierarchy from Chen papers (2023, 2025) |
| Pitfalls | HIGH | 9 critical/moderate pitfalls identified from codebase bugs + MCWF literature + Julia testing patterns |

**Overall confidence:** HIGH

### Gaps to Address

**No gaps requiring resolution during planning.** All technical questions are answered:

- Two-stage sampling question: RESOLVED — flat scan is mathematically correct for uniform jump operators (current case); implement normalization check, not restructure
- CPTP verification approach: RESOLVED — check sum of K†K, use S matrix eigenvalue guard
- DM vs trajectory channel difference: RESOLVED — both channels converge to same Lindbladian in delta→0 limit; test at matched delta or compare both to expm(L*delta)
- Error scaling expectations: RESOLVED — O(delta^2) per step, O(delta) total; test with delta-sweep
- Test infrastructure pattern: RESOLVED — two-tier (deterministic property + statistical convergence)

**Only validation needed:** Run the tests after Phase 1 fixes to confirm bugs are resolved.

## Test System Specifications

**Recommended test systems:**
- 2-qubit Heisenberg: dim=4, for large-N_traj statistical tests (fast)
- 3-qubit Heisenberg: dim=8, primary validation target (rich spectrum, still fast)
- 4-qubit Heisenberg: dim=16, secondary validation (approaching "real" complexity)

**Recommended parameters:**
- beta: 1.0, 5.0, 10.0 (low/medium/high temperature)
- sigma: 1/beta (standard Chen choice)
- delta: 0.01, 0.05, 0.1 (0.01 is safe, 0.1 is aggressive)
- num_energy_bits: 10-12 (1024-4096 energy grid points)
- N_traj: 100 (smoke test), 500 (validation), 1000-2000 (statistical convergence)
- total_time: 5-20 * (1/spectral_gap) (long enough for steady state)

**Error tolerance hierarchy:**
- TOL_MACHINE = 1e-12 (exact operations)
- TOL_BOHR_EXACT = 1e-10 (Bohr+B → Gibbs)
- TOL_ENERGY_QUAD = 1e-4 (Energy domain quadrature error)
- TOL_TIME_QUAD = 1e-3 (Time domain adds time quadrature)
- TOL_TROTTER = 1e-2 (Trotter discretization)
- TOL_TRAJECTORY_STAT = 0.05 (1/sqrt(N_traj) sampling noise)

## Sources

### Primary (HIGH confidence — direct verification)

**Codebase analysis:**
- `/Users/bence/code/QuantumFurnace.jl/src/trajectories.jl` — identified 3 bugs (lines 53, 96, 472-477)
- `/Users/bence/code/QuantumFurnace.jl/src/jump_workers.jl` — DM reference implementation
- `/Users/bence/code/QuantumFurnace.jl/src/furnace.jl` — run_lindbladian, run_thermalization APIs
- `/Users/bence/code/QuantumFurnace.jl/src/coherent.jl` — B operator computation
- `/Users/bence/code/QuantumFurnace.jl/src/qi_tools.jl` — trace_distance, fidelity metrics

**Papers (direct reading):**
- Chen, Kastoryano, Brandao, Gilyen (2023) "Quantum Thermal State Preparation" [arXiv:2303.18224] — Theorem III.1 (O(delta^2) per step, O(T^2/epsilon) cost), weak measurement scheme
- Chen, Kastoryano, Gilyen (2025) "An efficient and exact noncommutative quantum Gibbs sampler" [arXiv:2311.09207] — KMS detailed balance, coherent term B, Lindbladian stationarity

**Documentation:**
- [QDYN: Quantum Jump Method](https://ag-koch.gitpages.physik.fu-berlin.de/qdyn/main/concepts/mcwf.html) — ensemble averaging, statistical error, jump selection
- [QuantumOptics.jl: Quantum trajectories](https://docs.qojulia.org/timeevolution/mcwf/) — MCWF validation methodology
- [QuantumToolbox.jl: Monte Carlo Solver](https://qutip.org/QuantumToolbox.jl/stable/users_guide/time_evolution/mcsolve) — jump selection, convergence patterns

**Stack verification:**
- [StableRNGs.jl GitHub](https://github.com/JuliaRandom/StableRNGs.jl) — v1.0.1, LehmerRNG, Big Crush passed
- [HypothesisTests.jl docs](https://juliastats.org/HypothesisTests.jl/stable/parametric/) — OneSampleTTest, ChisqTest APIs
- [StatsBase.jl docs](https://juliastats.org/StatsBase.jl/stable/) — sem(), L2dist, kldivergence
- [Aqua.jl docs](https://juliatesting.github.io/Aqua.jl/dev/) — test_all() API
- [Julia 1.12 Test stdlib](https://docs.julialang.org/en/v1/stdlib/Test/) — @testset rng= keyword

### Secondary (MEDIUM confidence — community consensus)

- Abdelhafez et al. (2019) "The Monte Carlo wave-function method: A robust adaptive algorithm and a study in convergence" [arXiv:1803.08589] — convergence properties, discretization error analysis
- Rall et al. (2025) "A Randomized Method for Simulating Lindblad Equations and Thermal State Preparation" [arXiv:2407.06594] — per-step O(λ²τ²) error, multi-step O(T²/M) error
- Ding, Li, Lin (2024) "Efficient quantum Gibbs samplers with KMS detailed balance" [arXiv:2404.05998] — alternative KMS construction with discrete jumps

### Tertiary (Domain knowledge — standard quantum information)

- CPTP map completeness: sum K_i† K_i = I (fundamental quantum channels requirement)
- 1/sqrt(N) convergence of Monte Carlo ensemble averages (Central Limit Theorem)
- Trace distance and fidelity as state comparison metrics (Fuchs-van de Graaf inequalities)

---
*Research completed: 2026-02-13*
*Ready for roadmap: yes*
