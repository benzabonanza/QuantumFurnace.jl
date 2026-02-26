# Phase 34: Code Deduplication - Context

**Gathered:** 2026-02-26
**Status:** Ready for planning

<domain>
## Phase Boundary

Extract shared helpers to replace copy-pasted patterns across simulation types (Thermalize, Krylov, Trajectory). Domain variants (EnergyDomain/TimeDomain/BohrDomain) stay as near-duplicates — the deduplication targets are cross-simulation-type patterns.

Prerequisite: Fix krylov_matvec sandwich convention bug (quick task) before deduplication.

</domain>

<decisions>
## Implementation Decisions

### domain_prefactor
- **Name:** `domain_prefactor` (purely scalar constant, no NUFFT/OFT involvement)
- **Dispatch:** Domain-only, no construction type needed
  - `domain_prefactor(::EnergyDomain, w0, sigma)` → `w0 / (σ√(2π))`
  - `domain_prefactor(::TimeDomain, w0, sigma, t0)` / `domain_prefactor(::TrotterDomain, w0, sigma, t0)` → `w0 * t0² * σ√(2/π) / (2π)`
  - BohrDomain: **no method** — callers use `gamma_norm_factor` directly, don't call `domain_prefactor`
- **Arguments:** Explicit fields (w0, sigma, t0), NOT Config
- **Called in:** `_precompute_data`, result stored in `precomputed_data` NamedTuple
- **Callers multiply:** `precomputed_data.domain_prefactor * gamma_norm_factor * [their own scaling]`

### Frequency iteration (SKIPPED)
- **Decision:** Keep the hermitian half-grid for-loop pattern as-is in each call site
- **Rationale:** User prefers the explicit for-loop structure; deduplication effort should target higher-level functions, not micro-patterns

### OFT unification
- **Merge:** `oft!` + `_krylov_oft!` → single `oft!` with optimized signature
- **Signature:** `oft!(out, eigenbasis::Matrix, bohr_freqs::Matrix, energy, inv_4sigma2)` — concrete types, precomputed `inv_4σ²`
- **Old JumpOp-based `oft!` signature:** Removed
- **File:** All OFT functions in `ofts.jl` (unified `oft!`, `time_oft!`, `trotter_oft!`)
- **Scope:** EnergyDomain only. Time/TrotterDomain production path uses NUFFT prefactors as-is
- **`inv_4sigma2`:** Computed at call site (`1/(4σ²)`), NOT stored in precomputed_data
- **Test utilities:** `time_oft!` and `trotter_oft!` remain as clearly-marked test/debug utilities

### Cross-simulation-type deduplication (NEW — identified during discussion)

Three shared helpers to extract into `furnace_utensils.jl`:

1. **R accumulation** (`_accumulate_R!`)
   - Same formula in 3 places: `R += Σ_ω rate²(ω) * (L†L)`
   - Thermalize: per-operator R^a (per step, hot path)
   - Krylov: summed R_total across all operators (once at construction)
   - Trajectories: per-operator R^a (once at framework build)
   - Difference: Thermalize per-step vs others precomputed; Trajectories rescales by n_jumps
   - Domain variants stay inline (Bohr/Energy/Time each have their own OFT)

2. **CPTP channel construction** (`_build_kraus_matrices!`)
   - Same formula in 3 places: `K0 = I - αR`, `S = (2α-δ)R - α²R²`, `U_res = √psd(S)`
   - Thermalize: jump_workers.jl
   - Krylov: krylov_workspace.jl
   - Trajectories: trajectories.jl

3. **Physics-convention sandwich** (shared helper, after convention bug fix)
   - After fixing krylov_matvec to use `L * ρ * L†` everywhere
   - Thermalize and Krylov eigsolve already use physics convention
   - Extract shared `_accumulate_sandwich!(out, L, rho, scalar, ws)` using `L * ρ * L†`

### Sandwich convention bug (prerequisite quick task)
- krylov_matvec.jl `_accumulate_sandwich!` uses `conj(L) * ρ * L^T`
- krylov_eigsolve.jl `_accumulate_jump_sandwich!` and Thermalize use `L * ρ * L†`
- These are mathematically different for complex L: `conj(L)ρL^T = (LρL†)^T`
- Tests pass because Pauli operators in real eigenbasis are real (σx,σz) or purely imaginary (σy), making both conventions numerically equivalent
- **Fix:** Separate quick task before Phase 34 — unify to `L * ρ * L†` everywhere
- Also need to fix the adjoint variants consistently

### Claude's Discretion
- Exact internal structure of `_accumulate_R!` (how to parameterize per-operator vs all-operator accumulation)
- Whether sandwich helper needs separate forward/adjoint variants or can be parameterized
- Order of extraction (prefactor → OFT → R → CPTP → sandwich, or different sequencing)

</decisions>

<specifics>
## Specific Ideas

- User wants simplification of "layers" — combining higher-level functions across simulation types, not micro-pattern dedup
- Domain variants (Energy/Time/Bohr) intentionally stay as near-duplicates — different OFT methods make them legitimately different
- DLL future construction type will have different NUFFT prefactor matrices than KMS — dispatch on (domain, construction) may be needed at OFT level eventually, but `domain_prefactor` itself is domain-only
- `_precompute_data` is the natural place to compute and store `domain_prefactor` since all necessary fields are available there

</specifics>

<deferred>
## Deferred Ideas

- Sandwich convention consolidation for BohrDomain 2-operator sandwich (`_accumulate_sandwich_2op!` / `_accumulate_adjoint_sandwich_2op!`) — different enough to warrant separate treatment
- `_run_chunk_*` variant consolidation (3 variants differing only in measurement flags) — minor, can be done in Phase 38 test cleanup
- `B_time`/`B_trotter` single/vector deduplication in coherent.jl — low priority
- `_pick_transition_kms`/`_pick_transition_gns` consolidation — belongs in Phase 35 or later when construction type dispatch is more mature

</deferred>

---

*Phase: 34-code-deduplication*
*Context gathered: 2026-02-26*
