# Project Research Summary

**Project:** QuantumFurnace.jl v2.2 — Hamiltonian Simulation Time Counter
**Domain:** Analytical quantum algorithm cost accounting for quantum Gibbs samplers
**Researched:** 2026-03-04
**Confidence:** HIGH

## Executive Summary

This milestone adds a purely analytical cost accounting layer to QuantumFurnace.jl that computes how much Hamiltonian simulation time the quantum Gibbs sampler algorithm (Chen 2023 / Chen-Kastoryano-Gilyen 2024) would require on a real quantum computer. The computation is simple arithmetic over the QPE grid and the existing b_minus/b_plus dictionaries — no simulation is run, no matrices are constructed, no eigensolver is invoked. All required numerical infrastructure already exists in the codebase; the implementation is approximately 200-350 lines across one new file (`src/simulation_time.jl`) with zero new external dependencies.

The recommended approach follows the established QuantumFurnace pattern: a single immutable result struct (`SimulationTimeBudget`) returned by a pure function (`compute_simulation_time`) that accepts a `HamHam` plus scalar parameters, with convenience overloads for `MixingTimeEstimate` chaining. The critical design constraint is that the quantum algorithm's cost is domain-independent — EnergyDomain vs TrotterDomain is a classical simulation distinction, and the cost counter must not dispatch on domain. It must also not couple to the `Config` struct, which carries simulation-specific infrastructure irrelevant to quantum cost counting.

The primary implementation risk is conceptual: confusing "what the classical simulator computes" with "what the quantum computer pays for." Six critical pitfalls document specific, potentially order-of-magnitude-wrong mistakes: using the truncated energy grid instead of the full `2^r` QPE grid for the inner time sum, applying filter weights to the QPE time sum, miscounting the factor of 2 for the weak measurement channel, treating the B coherent term as a single sum instead of a double sum, not applying the 2x factor correctly (OFT only, not B), and using classical simulation mixing time as a proxy for the quantum algorithm's mixing time. All six pitfalls are preventable with closed-form validation: the QPE inner time sum has the exact closed form `t0 * N^2/4`, which is sigma-independent and domain-independent.

## Key Findings

### Recommended Stack

No new dependencies are required. The implementation uses only Julia's standard library (`sum`, `abs`, `ceil`, `log2`, `size`) and calls existing internal functions already present in the codebase. The `Project.toml` requires zero changes. This is possible because the cost computation is elementary grid arithmetic: sums of `|t_k|` over a discrete QPE grid, with truncated b-function dictionaries computed by already-existing code in `coherent.jl`.

**Core technologies:**
- Julia stdlib only (>= 1.9, existing requirement) — all counting operations are `sum(abs, collection)` arithmetic; no new packages
- `_create_energy_labels(r, w0)` (energy_domain.jl) — full QPE grid generation, already implemented
- `_compute_truncated_func` / `_compute_b_minus` / `_compute_b_plus*` (coherent.jl) — B term truncated dictionaries, already implemented
- `_select_b_plus_calculator` (furnace_utensils.jl) — transition weight dispatch, already implemented

See `.planning/research/STACK.md` for complete dependency analysis.

### Expected Features

The feature set is compact and well-scoped. All features are small enough to ship in one milestone with no deferrals necessary; the ordering below reflects implementation dependencies, not priority.

**Must have (table stakes):**
- `SimTimeBudget` result struct — immutable, following `MixingTimeEstimate` pattern, with full docstring
- OFT Hamiltonian simulation time — the primary cost component: sum over the full `2^r` QPE grid using Gaussian-weighted `|t_k|`
- B coherent term Hamiltonian simulation time — double-sum structure using existing truncated b-dicts (KMS only; 0.0 for GNS)
- Total step count — `ceil(mixing_time / delta)`, conservative and correct
- Per-step cost and total cost assembly — `2 * OFT_time + B_time` (KMS) or `2 * OFT_time` (GNS)
- All 3 transition weight functions — Gaussian, Metropolis (a=0), Smooth Metropolis (a>0, b>0); affects B term and effective energy count
- KMS vs GNS construction toggle — `with_coherent` trait already exists; B returns 0.0 for GNS
- QPE grid parameter reporting — `r`, `N`, `w0`, `t0`, energy range, for paper resource estimation tables

**Should have (competitive differentiators):**
- Per-component breakdown in result struct — `oft_time`, `b_time`, `per_step_time`, `n_steps`, `total_time` as separate fields
- `compute_simulation_time(::MixingTimeEstimate, ...)` overload — chains naturally from mixing time estimation workflow
- Rescaled and physical time reporting — store `rescaling_factor` in struct; physical time = rescaled time * rescaling_factor
- `HamHam` convenience overload — extracts `rescaling_factor` and `num_qubits` automatically from the struct

**Defer (v2+):**
- Trotter step counting for the quantum algorithm's internal circuit — deepens into gate-level compilation, out of scope for this milestone
- Gate-level compilation — depends on choice of Hamiltonian simulation method (Trotter order, QSP, etc.)
- Asymptotic complexity analysis — this counter provides exact numerical values, not big-O bounds

See `.planning/research/FEATURES.md` for full API surface specification and formula inventory.

### Architecture Approach

The implementation is a single new file `src/simulation_time.jl` (200-350 lines) included after `src/mixing.jl`, following the same established pattern as `fitting.jl` and `mixing.jl` — standalone post-processing modules with their own structs, no modifications to any existing struct or function. The module exports exactly two symbols: `SimulationTimeBudget` and `compute_simulation_time`. Internal helpers carry underscore prefixes per codebase convention. Zero changes to `HamHam`, `Config`, `Workspace`, or any existing result struct.

**Major components:**
1. `SimulationTimeBudget` (struct) — immutable result container with per-step breakdown, total cost, grid info, and parameter provenance
2. `compute_simulation_time` (main API) — pure function accepting `HamHam` + scalar params, dispatching internally on construction type and transition weight symbol
3. `_qpe_grid_info` (internal) — computes `N = 2^r`, `t0 = 2pi/(N*w0)`, and validates the Fourier relation
4. `_oft_hamiltonian_time` (internal) — sums `|t_k|` over full QPE grid; sigma-independent in the time sum; truncated grid for energy outer sum
5. `_b_hamiltonian_time` (internal) — double sum over truncated b-dicts using `_compute_truncated_func`; returns 0.0 for GNS
6. Transition weight dispatch — via `:gaussian`/`:metropolis`/`:smooth_metro` Symbol keyword, no coupling to `Config`

The deliberate anti-patterns to avoid: coupling to `Config` struct (carries irrelevant simulation state), domain dispatch on `EnergyDomain`/`TrotterDomain` (the quantum algorithm has no domain), adding fields to existing structs (breaks BSON serialization), and any mutable state or workspace allocation.

See `.planning/research/ARCHITECTURE.md` for build order, data flow diagram, and anti-pattern rationale.

### Critical Pitfalls

1. **Full QPE grid vs truncated simulation grid (CRIT-01)** — The classical simulator's `_truncate_energy_labels` reduces the `2^r` QPE grid to ~200-400 active points; the quantum computer's QPE circuit processes ALL `2^r` time points. Using `precomputed_data.energy_labels` (the truncated set) for the inner QPE time sum underestimates cost by a factor of ~16-25 for typical parameters. Prevention: use `_create_energy_labels(r, w0)` (full grid) for the inner QPE time sum; use the truncated grid only for counting how many energy evaluations the outer sum performs.

2. **OFT cost must NOT apply filter function weights (CRIT-02)** — The Gaussian envelope `exp(-sigma^2 * t_k^2)` appears in LCU amplitudes, not in gate cost. The inner QPE time sum is `2 * t0 * sum_{k=-N/2}^{N/2-1} |k|` with no sigma dependence at all. A computed cost that changes when sigma changes is wrong. The factor of 2 accounts for both `e^{+iHt}` and `e^{-iHt}` evolutions (forward and backward Hamiltonian simulation).

3. **B term is a double sum, not a single sum (CRIT-03)** — The coherent B term has nested outer (over `t` via b_minus) and inner (over `s` via b_plus) loops with different Hamiltonian simulation times. Total B time = `n_outer * sum_s 4|s*beta|` + `n_inner * sum_t 2|t/sigma|`. Counting only the outer loop and ignoring the inner misses a potentially dominant contribution at large beta.

4. **Factor of 2 applies to OFT only, not to B (CRIT-04)** — The dissipative channel requires both `U` (implementing the channel) and `controlled-U^dagger` (for measurement reversal), giving factor 2 for OFT. The coherent unitary `exp(-i*delta*B)` is a single application — factor 1. Per-step cost is `2*OFT_time + 1*B_time`, NOT `2*(OFT_time + B_time)`. The error grows as beta increases and B becomes more significant.

5. **Classical mixing time is not the quantum mixing time (CRIT-05)** — `estimate_mixing_time` from v2.1 reflects classical simulation artifacts: floor effects scaling as `k_energy * delta + floor_Trotter`, Bohr frequency collision effects (documented in MEMORY.md), and Trotter approximation error. The quantum algorithm's ideal mixing time is different. Accept `mixing_time` as an explicit user-provided input; do not derive it from simulation state.

6. **Transition weight type changes effective energy grid size (CRIT-06)** — Gaussian, Metropolis, and Smooth Metropolis transitions have very different support widths. The number of active energy grid points (used for the outer sum over energies) varies significantly between them. The cost counter must use `_truncate_energy_labels` (or equivalent) with the correct weight function for each type, and report `n_effective_energies` as a diagnostic in the result struct.

See `.planning/research/PITFALLS.md` for the full inventory including 6 moderate and 3 minor pitfalls with detection heuristics.

## Implications for Roadmap

Based on the research, the phase structure follows a strict linear dependency chain: each phase provides exactly the inputs the next phase requires, and each phase is responsible for a distinct class of pitfalls. All four phases are needed. None can be parallelized because the struct must exist before any computation, and OFT must be validated independently before B term complexity is added.

### Phase 1: Struct and Grid Infrastructure

**Rationale:** Zero-risk foundation before any physics. Defining the output container and grid utilities first ensures the units convention (rescaled vs physical time) and the w0/t0 Fourier relation are decided and enforced before they can infect other code. The struct's field names and docstring also force explicit documentation decisions about what each quantity means.
**Delivers:** `SimulationTimeBudget` struct with full docstring, `_qpe_grid_info(r, w0)` returning `(N, t0, energy_range)` with Fourier relation validation, `_num_qubits(ham)` helper computing `Int(log2(size(ham.data, 1)))`, and unit tests validating grid arithmetic and closed-form consistency.
**Addresses:** `SimTimeBudget` struct (table stakes), QPE grid parameter reporting (table stakes), per-component breakdown (differentiator)
**Avoids:** MOD-02 (rescaling units convention decided upfront with `rescaling_factor` stored in struct), MOD-07 (Fourier relation enforced in `_qpe_grid_info`, not left to callers), MIN-01 (clear struct naming hierarchy for per-step vs total)

### Phase 2: OFT Time Counting

**Rationale:** The OFT is the dominant cost component and has the cleanest closed-form validation (`T = 2 * t0 * N^2/4` for large N). Implementing and validating it first before B term complexity establishes confidence in the grid arithmetic approach and the transition weight dispatch pattern that the B term will reuse.
**Delivers:** `_compute_transition_weight(w, weight_type, ...)` lightweight dispatcher for `:gaussian`/`:metropolis`/`:smooth_metro`, `_oft_hamiltonian_time(r, w0, sigma, transition_weight, ...)` computing both the full-grid QPE time sum (inner) and the truncated energy count (outer), and unit tests verifying sigma-independence of the inner time sum, the factor-of-2 convention, and the closed-form `t0 * N^2/4` approximation.
**Addresses:** OFT Ham sim time (table stakes), all 3 transition weight functions (table stakes)
**Avoids:** CRIT-01 (full QPE grid for inner sum, truncated grid for outer energy count — explicitly separated), CRIT-02 (no filter weighting; sigma-independence tested), CRIT-06 (transition type correctly controls the truncated energy count)

### Phase 3: B Coherent Term Counting

**Rationale:** Builds on the transition weight dispatch from Phase 2. Comes after OFT because B is more complex (double-sum structure) and only applies to KMS. Having OFT independently validated avoids confusion about which phase introduced a bug.
**Delivers:** `_b_hamiltonian_time(r, w0, beta, sigma, transition_weight, a, b, ...)` using `_compute_truncated_func` to obtain the actual truncated b-dicts and computing the double sum explicitly, returning 0.0 unconditionally for GNS; unit tests verifying beta-dependence of the result, zero for GNS, and explicit inner vs outer sum terms.
**Addresses:** B coherent term Ham sim time (table stakes), KMS vs GNS toggle (table stakes)
**Avoids:** CRIT-03 (explicit double-sum with inner and outer loops computed separately, not conflated), CRIT-04 (factor-of-1 for B term established here as a named constant), MIN-04 (correct b_plus variant dispatch via reuse of `_select_b_plus_calculator`)

### Phase 4: Integration API and Exports

**Rationale:** Wires all independently validated components into the public-facing API. The `HamHam` and `MixingTimeEstimate` convenience overloads are pure wiring. Final integration tests validate the full pipeline and cross-check that `per_step_cost == 2 * oft_cost + b_cost` (not `2 * (oft_cost + b_cost)`).
**Delivers:** `compute_simulation_time(ham, r, delta, mixing_time; beta, ...)` main orchestrating API, `MixingTimeEstimate` convenience overload extracting `mixing_time` from the estimate, `HamHam` convenience overload extracting `rescaling_factor`, module wiring (`include("simulation_time.jl")` and exports in `src/QuantumFurnace.jl`), integration tests covering GNS vs KMS comparison, step count scaling, and physical vs rescaled time.
**Addresses:** Total step count (table stakes), total Ham sim time (table stakes), KMS vs GNS construction toggle (table stakes), chaining overload (differentiator), rescaled/physical time (differentiator), n_jumps factoring (differentiator)
**Avoids:** CRIT-04 (integration test asserts `2*OFT + B`, not `2*(OFT+B)`), CRIT-05 (mixing_time is an explicit required parameter, not derived internally), MOD-01 (`ceil` for step count; test that total cost scales linearly with `1/delta`), MOD-03 (per-step cost is per one jump, total cost reported separately), MOD-04 (GNS result has b_time=0.0; KMS result has b_time>0), MOD-05 (API accepts grid params explicitly; no domain dispatch), MOD-06 (beta required as keyword argument with no default)

### Phase Ordering Rationale

- The result struct is the output type of every computation function; it must exist first.
- OFT time sum has a trivial closed form (`t0 * N^2/4`) that enables confident independent validation — validate it before adding B term complexity.
- B term reuses the transition weight dispatch established in Phase 2; the shared code path is tested in both phases.
- Integration wiring cannot be tested until both OFT and B are implemented and individually validated.
- This 4-phase structure exactly mirrors the pitfall analysis's "Recommended Phase Ordering Based on Pitfall Dependencies" in PITFALLS.md, providing strong cross-document consistency validation.

### Research Flags

Phases with standard, well-documented patterns (no `/gsd:research-phase` warranted):
- **Phase 1 (Struct + Grid):** Entirely reuses the `FitResult`/`MixingTimeEstimate` pattern and standard QPE grid arithmetic. Fully resolved.
- **Phase 2 (OFT Time):** Formula is fully worked out with closed-form validation in both FEATURES.md and PITFALLS.md. Fully resolved.
- **Phase 3 (B Term):** Double-sum formula is documented in PITFALLS.md (CRIT-03) with explicit pseudocode. B-dict computation reuses existing internal functions. The b_plus dispatch is documented in MIN-04. Fully resolved.
- **Phase 4 (Integration):** Pure wiring following the established `include`/`export` module pattern. No research needed.

All four phases have sufficient clarity from codebase analysis and paper abstracts. The most complex physics (B term double-sum structure, LCU cost accounting) is fully resolved in the research documents.

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Stack | HIGH | All functions verified to exist by direct inspection of `src/*.jl`; `Project.toml` confirmed to require zero changes |
| Features | HIGH | Grounded in direct OFT/B/coherent code inspection plus paper abstracts; formulas are unambiguous; API surface fully specified |
| Architecture | HIGH | Mirrors `fitting.jl`/`mixing.jl` pattern exactly; build order derived from explicit dependency analysis; anti-patterns enumerated |
| Pitfalls | HIGH | Six critical pitfalls are concrete, verified against specific codebase code paths, include numerical examples (truncated grid ~200-400 pts vs full grid 4096 pts), and include detection heuristics |

**Overall confidence:** HIGH

### Gaps to Address

- **n_jumps factoring semantics:** FEATURES.md notes that whether `n_jumps` multiplies the total cost depends on the exact algorithm variant and the user's framing (per-jump vs per-step). The result struct should expose per-jump and total costs as separate fields and let the user decide whether to multiply by `n_jumps`. Document this ambiguity explicitly in the function docstring rather than silently choosing one convention.

- **Classical vs quantum mixing time gap:** CRIT-05 identifies that `estimate_mixing_time` output reflects classical simulation artifacts (floor effects, Trotter error, Bohr frequency collisions) not present in the ideal quantum algorithm. There is no automated way to quantify this gap within QuantumFurnace.jl. The function docstring must explicitly acknowledge this distinction and recommend comparing against the Krylov spectral gap (`run_krylov_spectrum`) as an independent sanity check.

- **B term Metropolis singularity validation:** The Metropolis `b_plus` function has a `1/t` singularity regulated by `eta`. The cost counter reuses `_compute_truncated_func` which already handles this, but a test verifying `length(b_plus_dict) > 0` for all three transition weight types at typical parameters would confirm correct dispatch.

## Sources

### Primary (HIGH confidence)
- `src/energy_domain.jl` — `_create_energy_labels`, `_truncate_energy_labels`, `pick_transition`, QPE grid generation
- `src/ofts.jl` — `time_oft!`, OFT prefactor formula `exp(-t^2*sigma^2 - i*omega*t)`, Gaussian envelope structure
- `src/coherent.jl` — `B_time`, `B_trotter`, `_compute_b_minus`, `_compute_b_plus*`, `_compute_truncated_func`, B double-sum structure
- `src/furnace_utensils.jl` — `_precompute_labels`, `_select_b_plus_calculator`, `_build_cptp_channel`, energy truncation logic
- `src/hamiltonian.jl` — `HamHam{T}` struct, `rescaling_factor` field, `_rescaling_and_shift_factors`
- `src/misc_tools.jl` — `validate_config!`, Fourier relation `t0 * w0 = 2*pi / N`
- `src/mixing.jl` — `MixingTimeEstimate` struct, pattern for result struct and convenience overloads
- `src/structs.jl` — `Config{S,D,C,T}`, `with_coherent` trait, `KMS`/`GNS` construction types
- `Project.toml` — current dependency list (zero changes required confirmed by inspection)

### Secondary (MEDIUM confidence)
- [Chen et al. 2023, arXiv:2311.09207](https://arxiv.org/abs/2311.09207) — GNS construction, OFT cost structure (abstract only; full text not read)
- [Chen-Kastoryano-Gilyen 2024, arXiv:2404.05998](https://arxiv.org/abs/2404.05998) — KMS construction, B term, total Ham sim time metric (abstract only; full text not read)
- [Chen et al. 2023, arXiv:2303.18224](https://arxiv.org/abs/2303.18224) — quantum thermal state preparation, original GNS (abstract only)

### Tertiary (context only)
- Project MEMORY.md — floor analysis (`floor = k_energy * delta + floor_Trotter`), Bohr frequency collision analysis, bi-exponential fitting (HIGH confidence for codebase facts; context for why classical mixing time estimates must not be used uncritically as quantum algorithm mixing times)

---
*Research completed: 2026-03-04*
*Ready for roadmap: yes*
