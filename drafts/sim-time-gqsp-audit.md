# `compute_simulation_time` audit — GQSP cost model + per-term registers (qf-e4z.18)

This note documents what `src/simulation_time.jl::compute_simulation_time` reports, what it *does not* report, and the two findings from the qf-e4z.18 audit:

1. **Per-term register triples (qf-9z0) are now surfaced in `SimulationTimeBudget`.** Old budgets only carried the dissipative register; the b_- and b_+ registers were used in the formula and discarded. Any pre-2026-05-07 saved budget is therefore ambiguous about the b_± parameters.

2. **The B-time cost model is now GQSP-aware.** The current code multiplies the per-block-encoding cost `b_per_be` by `2 · gqsp_degree` when `with_gqsp = true`, anticipating the Form-B circuit refactor tracked in [`qf-e4z.19`](#followup-form-b-refactor-qf-e4z19). Until that refactor lands, the cost model is *one line ahead of* the implementation by a factor of 1.5×.

Old `SimulationTimeBudget` BSONs in `scripts/output/` and `drafts/error-analysis/` predate both changes (single-register, GQSP-blind, pre-qf-9z0 + pre-qf-e4z.18) and **must not be reused** for plot composition.

## Per-step cost model

For `Config{Thermalize, <:Union{TimeDomain, TrotterDomain}}` (the implemented δ-channel):

```
per_step_time  = 2 · oft_time + b_time
total_time     = n_steps · per_step_time              n_steps = ⌈T / δ⌉
```

### `oft_time` — dissipative grid

The OFT integrand `Σ_{t̄ on D-grid} |t̄| · γ̃(ω(t̄))` (`_oft_hamiltonian_time`) reads the dissipative QPE register:

- `r_D = register_r_D(config)`
- `w0_D = register_w0_D(config)`
- `t0_D = 2π / (2^{r_D} · w0_D)` (Fourier relation)
- `transition_weights[i] = γ̃(ω_i)` clipped at NaN to 0 (extreme-energy `Inf · 0` cases in the erfc branches).

Closed form when `transition_weights ≡ 1`: `oft_time = t0_D · N_D² / 4` (validated in `test/test_simulation_time.jl::"OFT time — closed-form"`).

### `b_per_be` — per-block-encoding B-time

The legacy formula in `_b_hamiltonian_time` reads the **two** coherent registers separately (qf-9z0):

```
b_per_be = t0_bm · t0_bp · [4 · ‖b_-‖₁ · Σ_τ |b_+(τ)| · |τβ|
                              + 2 · ‖b_+‖₁ · Σ_t |b_-(t)| · |t/σ|]
```

with

- outer (`b_-`) register: `r_bm = register_r_b_minus(config)`, `w0_bm`, `t0_bm`.
- inner (`b_+`) register: `r_bp = register_r_b_plus(config)`, `w0_bp`, `t0_bp`.
- per (t, τ) pair: 3 inner Heisenberg evolutions (total 4|τβ|) + 2 outer time evolutions (total 2|t/σ|).

This is the cost of **one** application of the block encoding $U_{B_a}$, in Hamiltonian-simulation time units. It is GQSP-blind. Returns 0 for GNS (`with_coherent(GNS()) == false`).

### `b_time` — per-step B-time

```
b_time = with_gqsp ? 2 · gqsp_degree · b_per_be : b_per_be
```

The 2·d multiplier comes from MW2024 Theorem 6 / Eq. 46 (Form B):

- The Hamiltonian-simulation target $e^{-i\delta\alpha\cos\theta}$ is a *symmetric* Laurent polynomial in $z = e^{i\theta_\lambda}$ of bilateral degree $d$ (Jacobi-Anger truncation, MW Eq. 62).
- The underlying ordinary polynomial $P(z) = z^d L_d(z)$ has degree $D = 2d$ with shift $k = d$.
- MW Eq. 46 realises $L_d = z^{-d} P$ as **D = 2d slots: $(D-k) = d$ controlled-$W$ ($A$) + $k = d$ closed-controlled-$W^\dagger$ ($A'$ — fires on $|1\rangle$, MW Eq. 45)**, plus a closing rotation $L_0$. **Total block-encoding queries: $D = 2d$.**
- Each $W$ invocation is one application of $U_{B_a}$ (block encoding of $B_a$); the joint reflection $R_T$ inside $W = R_T · U_{B_a}$ is a Clifford and contributes nothing to Hamiltonian-simulation time. The $2d+1$ rotation triples on the QSP ancilla are also constant-cost single-qubit gates and contribute nothing.
- BS+MW Algorithm 1 produces the same $\{(\theta_j, \phi_j)\}_{j=0}^{2d}$ angles for both Form B (Eq. 46) and the equivalent Form C (Eq. 52, $D = 2d$ controlled-$W$ + uncontrolled $W^{-d}$ tail = $3d$ queries). The choice between them is purely a circuit-routing decision.

### What is NOT counted

- **PREP / un-PREP** state-preparation subroutines (gate-level cost; subdominant to controlled-Hamiltonian evolutions for $r_\pm \le 15$ per `drafts/coherent-step.md`). Future gate-complexity plots would need a separate counter.
- **Boltzmann rotations** (the controlled-$Y$ rotation that converts the OFT amplitude to a transition-weight amplitude). Same story: gate-level, subdominant.
- **GQSP ancilla rotations** ($2d + 1$ single-qubit gates). Constant gate cost per CoherentStep, not Hamiltonian-simulation time.
- **Reflection $R_T$** inside $W = R_T \cdot U_{B_a}$ — Clifford, gate-level only.

Per the spec (`supplementary-informations/the-main-numerical-plots.md` §6 (P6 precursor ii)): "I think these would be only counted if we have time to do gate complexity plots and not just Ham sim time plots." This audit confirms we are not counting them in `SimulationTimeBudget`, consistent with the spec's preference.

## Per-term register surfacing in `SimulationTimeBudget`

Old struct (pre-qf-e4z.18): one set of `(r, N, w0, t0)` and `energy_range`. The b_± per-term parameters were used in the formula and discarded.

New struct (qf-e4z.18):

| group | fields | notes |
|---|---|---|
| dissipative | `r_D, N_D, w0_D, t0_D, energy_range` | OFT grid |
| outer coherent | `r_bm, N_bm, w0_bm, t0_bm` | b_-(t) Riemann sum |
| inner coherent | `r_bp, N_bp, w0_bp, t0_bp` | b_+(τ) Riemann sum |
| GQSP | `with_gqsp::Bool, gqsp_degree::Int` | cost-model flag + multiplier base |
| cost | `oft_time, b_per_be, b_time, per_step_time, n_steps, total_time` | `b_per_be` is the new audit field |

The b_± fields are populated regardless of construction (`register_*_b_minus / _b_plus` fall back to legacy single-register fields when the per-term ones are not set; for GNS configs they record what the config has even though `b_per_be = 0`).

## Status of pre-existing budgets in `scripts/output/` and `drafts/error-analysis/`

Anything quoting Ham-sim time before today (2026-05-07) is a single-register + GQSP-blind artifact:

- Single-register: B-time formula uses `t0_outer = t0_inner = t0_D` instead of the per-term spacings. Materially wrong B-time when r_D ≠ r_b±, which is the common case — see e.g. `drafts/error-analysis/parameter-recommendations.md` §β=10 row "ε=1e-6, r_b+ = 14, r_D = 6": `t0_bp / t0_D ≈ 1/200`, so the old formula over-counts B-time by ~200×.
- GQSP-blind: missing the 2× / 3× multiplier.

These files **must not be re-used** for plot composition. P0a/P0b/P0c regenerate post-qf-e4z.18 numbers from scratch; S3/S5/S6 record the new fields per cell.

## Form B is the implementation target (qf-e4z.19, closed)

The cost-model multiplier `2 · gqsp_degree` matches the live GQSP circuit. The Python POC (`src/python/gqsp/circuit.py::gqsp_circuit(realise="L_d_form_B")`) and the v9 thesis figure both implement **Form B (MW Eq. 46)**: `d` open-controlled-$W$ slots + `d` closed-controlled-$W^\dagger$ slots (the `A'` of MW Eq. 45, fires on $|1\rangle$), no uncontrolled tail — `2d` block-encoding queries per `CoherentStep`. Same BS+MW Algorithm 1 angles drive both forms.

The mathematically equivalent **Form C (MW Eq. 52)**: all controlled-$W$ slots + uncontrolled $W^{-d}$ tail = `3d` queries. Form B and Form C produce the *same* post-selected unitary on the QSP=$|0\rangle$ block (algebraic identity from MW2024 Eqs. 49→53), so the project's classical Clenshaw evaluator `_gqsp_apply_polynomial` is form-agnostic — only the gate-level circuit (and hence the Hamiltonian-simulation cost) differs. Numerical cross-check at ≤ 1e-12 in `src/python/tests/test_gqsp.py::test_form_b_equivalent_to_form_c`.

History on the `2d` vs `3d` decision: an earlier draft committed to Form C from an empirical bug-check that used `A^† = |0\rangle\langle 0|\otimes U^\dagger + |1\rangle\langle 1|\otimes I` (open-controlled $U^\dagger$) instead of MW's `A' = |0\rangle\langle 0|\otimes I + |1\rangle\langle 1|\otimes U^\dagger` (closed-controlled $U^\dagger$). With the correct $A'$, MW Eq. 48 and the Eq. 49→53 commutation argument both go through, recovering the optimal `2d` cost. The v9 thesis text and the live Python POC use the corrected form.

## References

- `src/simulation_time.jl` — implementation.
- `test/test_simulation_time.jl` — 107 tests (struct shape, OFT closed-form, B-cost positivity, GQSP multiplier at d ∈ {1,2,3}, per-term register threading via `r_b_plus ∈ {8, 10}`).
- `drafts/qsp-subsection-v9.md` — current thesis QSP subsection (committed to Form C; targeted by qf-e4z.19 for revision).
- `supplementary-informations/Motlagh and Wiebe - 2024 - ...pdf` — Theorem 6 (Eqs. 45–53), Theorem 7, Eq. 62.
- `supplementary-informations/Berntson and Sünderhauf - 2025 - ...pdf` — independent of circuit form; Algorithm 1 (FFT-based $Q$ from $P$).
- `.claude-memory/feedback_construct_lindbladian_over_build_dense_superop.md`, `.claude-memory/per_term_registers_qf_9z0.md`, `.claude-memory/quadrature_register_recipe_qf_7xt.md` — context.
