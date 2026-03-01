# Project Research Summary

**Project:** QuantumFurnace.jl v2.1 — Speedup & Mixing Time
**Domain:** Julia performance optimization + quantum Lindbladian simulation
**Researched:** 2026-03-01
**Confidence:** HIGH

## Executive Summary

QuantumFurnace.jl v2.1 targets two tightly coupled goals: a significant performance improvement to the DM thermalization engine (`run_thermalize`) and a new scientific capability — mixing time estimation via exponential fitting of the trace distance convergence curve. The performance bottleneck is structural: the existing DM path calls `_build_cptp_channel` (which runs `eigen(Hermitian(S))`) at every single thermalization step, even though the result depends only on the static jump operators and fixed delta — not on the evolving density matrix. The trajectory engine already solved this exact problem via per-jump precomputation in `_build_trajectory_workspace`. Extending this pattern to the DM path is the single highest-impact change in the milestone.

The recommended approach follows the dependency order revealed by cross-research analysis: (1) per-jump precomputation to eliminate the per-step eigendecomposition, (2) `save_every` for trace distance to decouple observation cost from step cost, (3) multi-threaded BLAS management and optional omega-loop threading, and (4) mixing time estimation by promoting `staging/fitting.jl` to active code. This order is driven by pitfall severity — the precomputation must come first because it moves the allocation-heavy `eigen()` call out of the hot path, which is a prerequisite for safe omega-loop threading. All four research files converge on this ordering.

The key risk profile is dominated by Julia threading pitfalls, not by scientific modeling pitfalls. BLAS thread state is process-global, race conditions in per-task matrix accumulators are silent data corruptors, and `threadid()`-indexed buffers are unsafe due to task migration. All of these already have proven mitigations in the existing trajectory engine (`try/finally` BLAS scoping, per-chunk accumulators, `@spawn` over `@threads`). The mixing time estimation carries a separate risk: exponential fit model mismatch when the trace distance has not fully entered the single-exponential regime. The `skip_initial` parameter and R-squared quality gate (R^2 > 0.95 threshold) are the primary safeguards. The DM trace distance is a deterministic, noiseless signal — substantially cleaner than the noisy trajectory observables the fitting code was originally built for — making these safeguards effective.

## Key Findings

### Recommended Stack

No new fundamental dependencies are needed for the performance features. The only dependency addition is re-introducing **LsqFit.jl** (previously removed in v2.0 Phase 37 cleanup when `fitting.jl` was staged), plus moving `src/staging/fitting.jl` back into the active module. Julia 1.11+ is already required and fully supports the `@sync`/`Threads.@spawn` pattern and `BLAS.set_num_threads`. The BLAS threading strategy is adaptive: for dim <= 64, prefer Julia-level omega-loop parallelism with `BLAS.set_num_threads(1)`; for dim >= 128, prefer serial omega-loop with multi-threaded BLAS per GEMM call.

**Core technologies:**
- `LinearAlgebra.BLAS.set_num_threads` — adaptive BLAS thread control; mandatory try/finally scoping pattern already proven in trajectory engine
- `Threads.@spawn` + `@sync` — omega-loop parallelism; codebase idiom from trajectory path; avoids `@threads :static` composability issue
- `LsqFit.jl` (0.15+) — Levenberg-Marquardt exponential fitting; re-add to Project.toml; `staging/fitting.jl` is already feature-complete (217 LOC, comprehensive error handling)
- `OmegaLoopScratch` struct (new, internal) — per-task scratch buffers (3 dim×dim matrices each); independently allocated to avoid false sharing
- No GPU, no OhMyThreads.jl, no multi-exponential fitting — all evaluated and rejected as over-engineered for current system sizes

### Expected Features

**Must have (table stakes):**
- **TS-01: Per-jump CPTP channel precomputation** — eliminate `eigen(Hermitian(S))` from the hot loop; directly mirrors the existing trajectory workspace pattern; biggest single speedup
- **TS-02: Precomputed CPTP application** — replace `_finalize_kraus_step!` (allocates K0/U_residual every step) with `_apply_precomputed_channel!` (4 `mul!` calls, zero allocation); follows directly from TS-01
- **TS-03: `save_every` for `run_thermalize`** — control trace distance computation frequency; identical pattern to existing trajectory `save_every`; eliminates O(dim^3) eigendecomposition from every step
- **TS-04: Mixing time estimation via exponential fit** — activate `staging/fitting.jl`, write `estimate_mixing_time(result::ThermalizeResults)` as a post-processing function; returns `MixingTimeEstimate` with gap, R^2, CI, and optional extrapolated t_mix
- **TS-05: Multi-threaded BLAS for DM thermalization** — ensure BLAS thread count is at default (not 1) during `run_thermalize`; significant speedup for dim >= 64; requires no code changes if called independently of trajectory engine

**Should have (competitive):**
- **DIFF-02: Effective rate diagnostic** — compute and return `lambda_eff(t) = -log(d(t+tau)/d(t))/tau` alongside the fit; provides model-free validation of whether the single-exponential regime has been reached; low complexity, high scientific value
- **DIFF-03: Fit quality gates with actionable warnings** — structured `@warn` messages for R^2 < 0.95, negative offset C, gap uncertainty > 100%, and extrapolation ratio > 10x; prevents users from trusting bad mixing time estimates

**Defer (v2+):**
- **DIFF-01: Multi-threaded omega-loop for rho_jump accumulation** — requires per-task accumulators, BLAS single-thread scoping, and careful `hermitianize`-after-merge ordering; adds complexity without clear benefit until benchmarking confirms the omega-loop is the post-precomputation bottleneck
- Two-exponential fitting — ill-conditioned (Prony problem); `skip_initial` on single-exponential handles the same transient modes
- Richardson extrapolation in delta — diagnostic tool, not production feature; document instead
- GPU acceleration — not needed for n <= 12 (dim <= 4096); multi-threaded BLAS is sufficient

### Architecture Approach

The v2.1 architecture separates the DM thermalization engine into a one-time precomputation phase and a per-step hot loop, exactly mirroring the trajectory engine's established two-phase structure. The key insight from architecture research: the omega-loop in `_jump_contribution!` interleaves R accumulation (rho-independent) and rho_jump accumulation (rho-dependent). These must be surgically separated. New functions `_precompute_R` (shared with trajectories, moved to `furnace_utensils.jl`) and `_accumulate_rho_jump_only!` (domain-dispatched, 3 variants) replace the monolithic `_jump_contribution!` in the hot path. Mixing time estimation is explicitly a post-processing concern — kept separate from `run_thermalize` to avoid coupling simulation and analysis logic.

**Major components:**
1. **Precomputation phase** (`run_thermalize` setup) — `_precompute_R` + `_build_cptp_channel` per jump, stored as local `K0s`/`U_residuals` vectors; BohrDomain needs a new `_precompute_R` variant since the trajectory path never handles BohrDomain
2. **Hot loop** (per-step) — `_accumulate_rho_jump_only!` (omega-loop, rho-dependent only) + `_apply_precomputed_channel!` (4 BLAS GEMMs, zero allocation); plus conditional `trace_distance_h` at `save_every` intervals
3. **Post-processing** (`estimate_mixing_time`) — wraps `fit_exponential_decay` from promoted `fitting.jl`; separate call, not embedded in simulation; returns `MixingTimeEstimate` with quality gates and effective rate diagnostic

### Critical Pitfalls

1. **BLAS global thread state leaking across calls (CRIT-01)** — `BLAS.set_num_threads()` is process-global; if DM thermalize sets it to N and errors before restore, subsequent trajectory calls with Julia threading oversubscribe by N×. Prevention: mandatory `try/finally` save/restore at every threading entry point; existing `trajectories.jl:500-508` is the model.

2. **Race conditions in omega-loop accumulator matrices (CRIT-03)** — matrix `.+=` is not atomic; concurrent writes to shared `scratch.R` or `scratch.rho_jump` silently lose updates. Prevention: pre-allocate independent per-task accumulators indexed by chunk (not `threadid()`); merge serially after `@sync`; each task also needs its own `jump_oft`/`sandwich_tmp` scratch buffers.

3. **`threadid()`-indexed buffers unsafe from task migration (CRIT-03 sub-issue)** — Julia tasks can migrate between OS threads. Prevention: index accumulators by `enumerate(chunks)` task index, not `Threads.threadid()`; this matches the existing trajectory threading pattern.

4. **Hermitianize placement in threaded reduction (CRIT-06)** — hermitianizing each per-task R before merging vs. after produces different floating-point results (O(ntasks × dim^2 × eps)); the latter matches serial semantics. Prevention: always hermitianize after the reduction sum; use `isapprox` with realistic tolerance in serial-vs-threaded tests.

5. **Exponential fit model mismatch from premature fitting window (MOD-01)** — the trace distance curve has a multi-mode transient before settling into single-exponential decay; fitting the transient produces a biased gap. Prevention: `skip_initial = 0.2` default; R^2 > 0.95 quality gate; cross-validate against Krylov spectral gap for small systems.

## Implications for Roadmap

Based on combined research, the suggested phase structure reflects a strict dependency order where each phase is both buildable and testable in isolation, and where the pitfall-heaviest work (threading) is deferred until the prerequisite (allocation-free hot path) is in place.

### Phase 1: Per-Jump CPTP Channel Precomputation
**Rationale:** Addresses the single biggest performance bottleneck (per-step `eigen()`) without introducing threading complexity. Must come first because CRIT-02 (stale eigendecomposition) and MOD-05 (changed numerical expectations for regression tests) are best resolved in isolation before threading is added. Also a prerequisite for safe omega-loop threading — a zero-allocation hot path is needed before parallelizing it.
**Delivers:** 2–10× speedup on `run_thermalize` depending on system size; verified numerical equivalence between precomputed and recomputed paths; updated regression test baselines.
**Addresses:** TS-01 (precomputation), TS-02 (precomputed channel application).
**Avoids:** CRIT-02 (assert S eigenvalues >= -eps per jump during construction), MOD-05 (regenerate regression references), MIN-03 (warn if precomputed data exceeds 100 MB).

### Phase 2: save_every for Trace Distance
**Rationale:** Independent of threading; can be done in parallel with or immediately after Phase 1. Generates the time-series data format that mixing time estimation (Phase 4) depends on. Off-by-one bugs in the time grid are much easier to debug in a serial, single-threaded context.
**Delivers:** Configurable observation frequency in `run_thermalize`; backward-compatible (default save_every=1); `save_every` stored in `ThermalizeResults` for downstream use.
**Addresses:** TS-03 (save_every), MIN-02 (backward compatibility via default).
**Avoids:** MOD-03 (copy the exact formula from trajectories.jl; assert `length(trace_distances) == length(time_steps)`).

### Phase 3: BLAS Thread Management and Optional Omega-Loop Threading
**Rationale:** Depends on Phase 1 (allocation-free hot path makes threading safe) and Phase 2 (save_every reduces observation overhead, making threading benefit more visible). Omega-loop threading is the most pitfall-dense phase and is explicitly optional for v2.1 — the BLAS threading (TS-05) alone provides meaningful speedup for dim >= 64 with no race condition risk.
**Delivers:** Explicit BLAS thread management in `run_thermalize`; optionally thread-parallel rho_jump accumulation for large frequency grids; verified serial-vs-threaded numerical equivalence within expected FP tolerance.
**Addresses:** TS-05 (BLAS threads), DIFF-01 (omega-loop threading, optional/deferred).
**Avoids:** CRIT-01 (try/finally BLAS scoping), CRIT-03 (per-chunk accumulators, no `threadid()`), CRIT-04 (independently allocated scratch matrices, not 3D views), CRIT-05 (set BLAS threads to 1 before Julia threading), CRIT-06 (hermitianize after merge), MOD-04 (verify zero allocation in hot path), MIN-01 (threading threshold dim >= 32 && n_omega >= 20).

### Phase 4: Mixing Time Estimation
**Rationale:** Depends on Phase 2 (`save_every` produces the clean time-series input). The scientific deliverable of the milestone. Can be developed and tested independently of threading (Phase 3). Activating `staging/fitting.jl` has no risk to existing code paths.
**Delivers:** `estimate_mixing_time(result::ThermalizeResults)` post-processing function; `MixingTimeEstimate` struct with fitted gap, confidence intervals, optional extrapolated t_mix; effective rate diagnostic; fit quality warnings; LsqFit re-added to Project.toml; `fitting.jl` promoted from staging.
**Addresses:** TS-04 (mixing time estimation), DIFF-02 (effective rate), DIFF-03 (quality gates).
**Avoids:** MOD-01 (skip_initial=0.2 default, R^2 quality gate), MOD-02 (extrapolation ratio limit, propagate gap_se to t_mix uncertainty), MOD-06 (SingularException handling already in staging code; add data quality pre-check).

### Phase Ordering Rationale

The ordering is driven by three constraints from combined research:

- **Allocation prerequisite for threading:** The per-step `eigen()` in Phase 1 allocates on every thermalization step. Threading a code path with per-step allocations causes GC pressure that stalls all threads simultaneously (MOD-04). Phase 1 must come before Phase 3.
- **Time-series format prerequisite for fitting:** Phase 4's `estimate_mixing_time` takes `ThermalizeResults` with `time_steps` aligned to saved trace distances. Phase 2 establishes this format. Phase 2 must come before Phase 4.
- **Pitfall isolation:** Phases 1 and 2 are serial changes with deterministic behavior. Phase 3 introduces non-determinism and must be tested against Phase 1's verified precomputed results. Interleaving threading with numerical correctness work would make defects ambiguous.

### Research Flags

Phases likely needing deeper research during planning:
- **Phase 3 (omega-loop threading):** The specific threading threshold (dim >= 32 && n_omega >= 20 from STACK.md) is a rough estimate from BLAS benchmark data, not measured on the actual system. Requires benchmarking on target hardware before enabling by default. Also, the BohrDomain `_precompute_R` variant (needed in Phase 1) has no existing precedent — the R accumulation logic is interleaved with rho_jump in the current BohrDomain `_jump_contribution!` and must be carefully untangled.
- **Phase 4 (mixing time estimation):** The burn-in fraction (skip_initial = 0.2) and the R^2 threshold (0.95) are based on general exponential fitting practice, not system-specific calibration. Cross-validating the fitted gap against `run_krylov_spectrum` results for small systems is recommended during implementation.

Phases with standard patterns (skip research-phase):
- **Phase 1 (precomputation):** Pattern is directly proven in the trajectory engine (`_build_trajectory_workspace`). The code to adapt already exists; this is primarily a refactor and extension.
- **Phase 2 (save_every):** Direct copy of the pattern from `trajectories.jl`; no novel design decisions.

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Stack | HIGH | All threading patterns verified against Julia official docs, PSA blog, and ITensors.jl recommendations. LsqFit.jl is already tested in staging. The only MEDIUM item is the precise dim/n_omega threading threshold, which requires empirical tuning. |
| Features | HIGH | Feature set derived from exhaustive codebase reading, directly grounded in the current v2.0 code structure, existing staging code, and theoretical papers (Chen et al. 2023, Ramkumar & Soleimanifar 2024). Feature dependencies are explicit and verified. |
| Architecture | HIGH | Analysis based on full source code reading, not external docs. Data flow diagrams confirm component boundaries. The BohrDomain gap (no existing `_precompute_R` variant) is the only unproven element — approach is clear but requires careful implementation. |
| Pitfalls | HIGH | Critical pitfalls (CRIT-01 through CRIT-06) are all grounded in official Julia documentation, Julia GitHub issues, and the Julia PSA blog. Moderate pitfalls grounded in numerical analysis literature and direct codebase analysis. |

**Overall confidence:** HIGH

### Gaps to Address

- **BohrDomain `_precompute_R` implementation:** This is new code with no codebase precedent. The R accumulation loop in `_jump_contribution!` for BohrDomain (jump_workers.jl:218-282) must be read carefully to extract the rho-independent R part from the rho-dependent rho_jump part. Plan for extra validation time in Phase 1.
- **Threading threshold calibration:** The recommended `dim >= 32 && n_omega >= 20` threshold for omega-loop parallelism is an estimate. During Phase 3, benchmark on the target system (both small and large problem sizes) before hardcoding the default.
- **Regression test baseline update:** Per-jump precomputation produces results that differ from per-step recomputation at O(1e-14) level due to eigendecomposition floating-point variation. Existing regression test BSON references must be regenerated after Phase 1. Flag this explicitly in the Phase 1 task.
- **save_every final-step alignment:** Architecture research notes two options for whether the final step is always saved regardless of `save_every`. The decision must be made explicitly and documented — the current trajectory code does NOT force-save the final step.

## Sources

### Primary (HIGH confidence)
- Codebase analysis: `furnace.jl:143-223`, `jump_workers.jl:172-440`, `furnace_utensils.jl:30-200`, `trajectories.jl:60-300`, `staging/fitting.jl`, `structs.jl:291-417` — direct code inspection, all architectural claims verified
- [Julia PSA: Don't Use threadid()](https://julialang.org/blog/2023/07/PSA-dont-use-threadid/) — per-chunk accumulator pattern, task migration safety
- [Julia Multi-Threading Documentation](https://docs.julialang.org/en/v1/base/multi-threading/) — `@spawn`, `@sync`, `@threads` semantics
- [Julia GitHub Issue #44201](https://github.com/JuliaLang/julia/issues/44201) — confirms BLAS threads are process-global, no per-task control
- [ITensors.jl Multithreading](https://itensor.github.io/ITensors.jl/dev/Multithreading.html) — BLAS vs Julia thread strategy for dense linear algebra
- [LsqFit.jl Documentation](https://julianlsolvers.github.io/LsqFit.jl/latest/) — curve_fit API, Levenberg-Marquardt, confidence intervals
- Chen et al. 2023, Propositions II.3, III.1 — theoretical mixing time bound and CPTP channel construction
- Ramkumar & Soleimanifar 2024, Theorems 2.1-2.2 — spectral gap bounds
- `supplementary-informations/error_catalogue_spectral_gap_estimation.md` — Errors 1, 6 (bias, multi-mode fitting)

### Secondary (MEDIUM confidence)
- [Julia for HPC: Multithreading (ENCCS)](https://enccs.github.io/julia-for-hpc/multithreading/) — nested parallelism, BLAS thread interaction
- [BLAS Thread Count Discussion (Julia Discourse)](https://discourse.julialang.org/t/ideal-number-of-blas-threads/79197) — dim thresholds for BLAS threading crossover
- [False Sharing in Julia Threading](https://blog.jling.dev/blog/false_share/) — independent allocation vs 3D array views
- [OhMyThreads.jl Documentation](https://juliafolds2.github.io/OhMyThreads.jl/stable/) — evaluated and rejected; confirms existing `@spawn` pattern is sufficient
- [Exponential Data Fitting, SDSU Research Report](https://www.csrc.sdsu.edu/research_reports/CSRSR2009-04.pdf) — Prony problem ill-conditioning for multi-exponential fits

---
*Research completed: 2026-03-01*
*Ready for roadmap: yes*
