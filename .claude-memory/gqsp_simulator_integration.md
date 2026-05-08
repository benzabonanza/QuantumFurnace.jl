---
name: GQSP simulator integration
description: with_gqsp option in Config replaces matrix-exp coherent step with the Jacobi-Anger Chebyshev polynomial f_d(B/α); default off, Time/TrotterDomain only.
type: project
---

# GQSP coherent step in the Julia simulator (epic qf-63j, completed 2026-04-30)

The Thermalize/Trajectory simulator can now optionally use the **GQSP polynomial
approximation** `f_d(B_a/α_a)` instead of `exp(-iδ B_a)` for the coherent step,
matching the thesis Alg. ref{alg:coh} prescription.

## Public API

Two new fields in `Config{S,D,C,T}`:

- `with_gqsp::Bool = false` — opt in to the GQSP polynomial. Default off; existing
  scripts and tests are unchanged.
- `gqsp_degree::Int = 1` — Jacobi-Anger truncation degree. d=1 is faithful to
  `O((δα)²)` and matches the splitting-error budget for typical numerics-chapter δ.

`validate_config!` rejects `with_gqsp=true` outside `KMS/DLL` constructions, outside
`Time/TrotterDomain`, or with `gqsp_degree ∉ [1, 100]`.

## Why: operator-level shortcut

By qubitization + Jacobi-Anger, the post-selected anc=|0⟩ block of the GQSP circuit is

  f_d(B/α) = J_0(δα) I + Σ_{k=1}^{d} 2(-i)^k J_k(δα) T_k(B/α)

so the simulator computes this **directly as a matrix function** — no need to build
the full block-encoding circuit. d=1 is just `J_0(δα) I − 2i J_1(δα) (B/α)` (scalar
axpy, cheaper than matrix-exp).

## Files touched

- `src/structs.jl` — `Config` schema additions
- `src/misc_tools.jl` — `validate_config!` GQSP validation block
- `src/QuantumFurnace.jl` — `using SpecialFunctions: erfc` extended to `besselj`
- `src/coherent.jl` — new `_gqsp_apply_polynomial` (d=1 fast path; d≥2 Clenshaw on
  T_k with three reused n×n scratch buffers) + `_gqsp_block_encoding_alpha`;
  `_precompute_coherent_unitary` dispatches on `config.with_gqsp` for Time/Trotter
- `src/trajectories.jl` — `_build_trajectory_workspace` per-jump U_B branch dispatches
  on `config.with_gqsp`
- `test/test_gqsp_config.jl`, `test/test_gqsp_polynomial.jl`, `test/test_gqsp_thermalize.jl`

## α formula (thesis-faithful, verified)

  α_a = γ_nf · t0_sim² · ‖b_-‖_{ℓ¹} · ‖b_+‖_{ℓ¹} · ‖A_a‖²_op

The two simulator factors `γ_nf` and `t0_sim²` compensate for `B_time`/`B_trotter`
returning a pre-scaled `B`. Triangle-inequality bound matches the helper exactly
(physics-checker confirmed). For the n=3 BETA=10/SIGMA=0.1 fixture the ratio
`opnorm(B_a)/α_a` lands in the 0.05–0.3 range — well under 1, comfortably in the
qubitization regime.

## Regime where d=1 is sufficient

d=1 truncation error: `‖f_1 − exp(-iδB)‖ ≤ 2(eδα/4)²`. For typical numerics-chapter
parameters (n=3–5, δ ∈ {1e-3, 1e-2}, βσ ∼ 1, normalized jumps) δα is well below 0.1
so the d=1 error is comfortably below the unraveling step error. **d ≥ 2 becomes
preferable** when δα approaches 1 — e.g. βσ ≳ 5 (the b_- cosh broadens linearly in
βσ), δ ≳ 0.1, very large t0, or unnormalized jumps with ‖A‖ ∼ O(1) per site.

### Per-step error budget (concrete numbers, n=3 fixture)

For the smoke-test config (n=3 KMS Time, δ=1e-3, n_jumps=9, δ_eff=9e-3): α_a ≈ 0.0185
(all 9 jumps), so δ_eff·α ≈ 1.66e-4. Per-step errors:

| Source | Per-step | Per-step value |
|--------|----------|---------------:|
| Lie-Trotter splitting (coh ⊕ diss) | O(δ_eff²) | 8e-5 |
| GQSP d=1 truncation | O((δ_eff α)²) | 3e-8 |
| GQSP d=2 truncation | O((δ_eff α)³) | 5e-12 |
| Matrix-exp truncation | machine ε | 1e-15 |

So GQSP d=1 sits **3000× below** the leading-order splitting error (the α² ratio).
**Strictly subleading**: d=1 does not degrade simulator accuracy at all in the
qubitization regime (α < 1). Diagnostic: `scripts/scratch_gqsp_trace_drift.jl`.

### Trace drift is a real-but-tiny side effect of d=1 non-unitarity

f_d(B/α) is unitary up to O((δα)^{d+1}); at d=1 the tiny non-unitarity causes
`tr(U ρ U†) = 1 + ε` per step with ε ≈ (δα)²/2. Drift accumulates **linearly** in
n_steps (not exponentially — the dissipative Kraus channel is exactly trace-preserving
by Chen Eq. 3.2, so it doesn't compound). For the smoke test (50 steps): predicted
drift `50 × 1.39e-8 ≈ 6.9e-7`, observed `6.78e-7` — 2% match. Matrix-exp baseline
and GQSP d=2 both stay at machine epsilon. **Not a bug**; the drift is below all
other error sources and disappears at d=2.

## Test coverage

- `test_gqsp_config.jl` — 8 validation tests
- `test_gqsp_polynomial.jl` — 33 polynomial-level tests (closed forms d=1,2;
  slope-(d+1) for d ∈ {1,2,3}; joint (B,α,δ)-invariance; tail-bound ratio;
  operator/circuit equivalence vs 1-ancilla qubitization; n=3 Heisenberg slope-2
  anchor with `B_time` + actual `_compute_b_minus`/`_compute_b_plus` kernels)
- `test_gqsp_thermalize.jl` — 16 integration tests (Time + Trotter smoke;
  regression vs matrix-exp baseline at small δ; higher d strictly tightens; single-
  trajectory regression at fixed seed)

Total test count after epic: **1917 (was 1860 baseline + 57 new)**.

## POC anchor

Operator-level identity `f_d(B/α) = [L_d(W)]_anc=0` was first verified end-to-end at
the **circuit level** by epic `qf-0x6` (closed). The simulator implements the
polynomial directly; the circuit-level POC scripts remain in
`scripts/scratch_gqsp_*.jl` as the verification anchor.

## Three thesis discrepancies still open (qf-45p)

The thesis chapters 1_preliminaries.tex (Fig. circ:gqsp) and 2_methods.tex
(Alg. alg:coh) currently disagree with the verified POC. The simulator implements
the **POC**, not the still-incorrect thesis text. Discrepancies tracked under epic
`qf-45p` (separate from this implementation epic).
