# Physics Check Report — qf-fzj.6

## Summary

Three modules verified against the source papers (Ding–Li–Lin 2024 arXiv:2404.05998 v4, Motlagh–Wiebe / Berntson–Sünderhauf GQSP). All three implementations match the paper definitions; one notation discrepancy is flagged (the code's "Eq. 3.4 / 3.5 / 3.7 / 3.22 / 4.7" labels are **thesis** numbering, not paper numbering — paper has these at Eq. 3.3 / 3.7 / 3.10 / 3.23 / 4.7-but-only-in-CKG-recovery-section). All physics/math is correct; numbering drift is documentation-only.

---

## E1. Motlagh-Wiebe GQSP (Jacobi-Anger Chebyshev polynomial)

**Equations checked**:
- Jacobi-Anger expansion: `e^{-iδx} = J_0(δ) + 2 Σ_{k=1}^∞ (-i)^k J_k(δ) T_k(x)` for `|x| ≤ 1`
- Block-encoding norm formula: thesis Algorithm `alg:coh`
- Bessel-tail bound: `‖f_d − e^{-iδB/α}‖ = O((δα)^{d+1})`

**Code locations**:
- `src/coherent.jl:283-294` — `_gqsp_block_encoding_alpha`
- `src/coherent.jl:314-374` — `_gqsp_apply_polynomial`
- `src/coherent.jl:247-265` — `_coherent_unitary_step` (new helper from A3)
- `src/python/gqsp/circuit.py:118-148` — POC angle-extractor (BS convention)

**Verdict per check**:

| # | Check | Status | Evidence |
|---|-------|--------|----------|
| 1 | Polynomial `f_d(B/α) = J_0(δα) I + Σ_{k=1}^d 2(-i)^k J_k(δα) T_k(B/α)` matches Jacobi-Anger truncation of `exp(-iδB)` | **PASS** | Standard identity. Code at `coherent.jl:325-356`: `a0 = besselj(0, δα)`, `ak = 2 cis(-π/2 · k) besselj(k, δα)` matches `2(-i)^k J_k(δα)` since `cis(-π/2 · k) = (-i)^k`. d=1 fast path at line 328-338 yields exactly `J_0 I − 2i J_1 (B/α)`. Clenshaw recurrence implements `b_k = a_k I + 2x b_{k+1} − b_{k+2}` and final `f = a_0 I + x b_1 − b_2`, the standard Chebyshev Clenshaw form. |
| 2 | α formula `γ_nf · t0² · ‖b_-‖_{ℓ¹} · ‖b_+‖_{ℓ¹} · ‖A‖²` ensures `‖B/α‖_op ≤ 1` (faithful block encoding) | **PASS** | Triangle inequality on the inner OFT structure of `B_time`/`B_trotter`. Code at `coherent.jl:290-293` computes exactly the bound; `_compute_b_minus`/`_compute_b_plus` produce truncated Dicts whose ℓ¹ sums are `sum(abs, values(...))`. The `t0²` factor matches the trapezoidal weight from `B = ... .* t0^2` at `coherent.jl:180, 229`. |
| 3 | Slope-(d+1) Bessel-tail bound | **PASS** | `J_k(z) ≈ (z/2)^k / k!` so `‖f_d − exp(-iδB/α)‖ ≤ 2(δα/2)^{d+1}/(d+1)! · 1/(1-δα/2)`. Slope `d+1` in `δ` at fixed `α`, scales as `α^{d+1}` in `α`. Tested at `test/test_gqsp_polynomial.jl:126-145` and `:166-185` — both pass (now consolidated in commit A5). |
| 4 | BS-vs-MW: Julia simulator is convention-agnostic | **PASS** | `_gqsp_apply_polynomial` (`coherent.jl:314-374`) computes `f_d(B/α)` directly via Clenshaw on `T_k(B/α)`. There is **no `L_0`** anywhere in the Julia code — no rotations, no controlled walks, no ancilla register. The function evaluates the post-selected `anc=|0⟩` block as a closed-form polynomial, identical for BS-form and MW-native circuits (they differ in angle assignment, not in the resulting scalar polynomial). The convention only matters in `src/python/gqsp/circuit.py`. **Cleared. Do NOT recommend "fixing".** |
| 5 | Hermitian invariant before GQSP call | **CONCERN→PASS (after A3)** | Pre-A3 asymmetry: `_precompute_coherent_unitary` Thermalize TimeDomain branch passed raw `B` to `_gqsp_apply_polynomial` while matrix-exp branch wrapped in `Hermitian(B)`. The new helper `_coherent_unitary_step` calls `hermitianize!(B)` unconditionally before the GQSP/exp branch, so both paths now operate on the same hermitised input. Resolved by A3. |
| 6 | `gqsp_degree ≤ 100` cap | **OK (overkill but harmless)** | `besselj(100, δα)` for `δα ≤ 1` is `< 10^{-150}` — far below ULP. The cap is a defensive guardrail, not a physical constraint. Default `d=1` produces `O((δα)²)` error matching the `O(δ²)` Lie-Trotter splitting, so production never exceeds `d=2`. |

**Concerns**: none. The GQSP implementation is paper-faithful and physically consistent.

**Recommended fixes**: none physics-related.

---

## E2. Berntson-Sünderhauf angle solver

**Equations checked**:
- BS Eq. 1.2: `M(z) = L_0 · ∏_{j=1}^N D(z) R_j` (transpose of MW Eq. 7)
- MW Eq. 7 / Fig. 2: `∏_{j=1}^d R_j A R_0` with `L_0 = R_0` rightmost

**Code locations**:
- `src/python/gqsp/circuit.py:118-125` — explicit comment on time order (BS-form layout)
- `.claude-memory/gqsp_bs_vs_mw_convention.md` — full convention note

**Verdict**:

| # | Check | Status | Evidence |
|---|-------|--------|----------|
| 1 | Julia simulator does not use angles | **PASS** | `_gqsp_apply_polynomial` evaluates the closed-form Chebyshev polynomial. No `thetas`/`phis` extraction, no rotations, no walk operator `W`. The angle solver is exclusively in the Python POC for circuit verification. |
| 2 | Closed-form Chebyshev expansion equals post-selected anc=0 block | **PASS** | Equivalence test at `test/test_gqsp_polynomial.jl:211-234`: `f_d(B/α) ≈ [L_d(W)]_{anc=0}` for `d ∈ {1,2,3}`, `δ ∈ {1e-1, 1e-2}`. 6 cases, all pass. Polynomial equivalence holds because `(M(z))_{1,1} = ((M(z))^T)_{1,1}`. |

**Concerns**: none. Quick check passes.

**Recommended fixes**: none.

---

## E3. Ding-Li-Lin DLL

### Numbering convention

The code/audit cite **thesis** numbering, not paper numbering. The DLL paper (arXiv:2404.05998 v4) uses:

| Code/audit cite | Actual paper Eq. | Description |
|---|---|---|
| "Eq. 3.4 first form" | **3.3** | `L_a = Σ q^a(ν) e^{-βν/4} A^a_ν = ∫ f^a(t) A^a(t) dt` |
| "Eq. 3.5" coherent kernel | **3.7** | `ĝ(ν) := -(i/2) tanh(-βν/4) κ(ν)` — single argument! |
| "Eq. 3.7" coherent G | **3.10** | `G = Σ_a Σ_ν ĝ(ν) (L_a† L_a)_ν = ∫ g(t) H_L(t) dt` — single time integral! |
| "Eq. 3.19" Metropolis q | **3.20** | `q(ν) = exp(-√(1+β²ν²)/4) w(ν/S)` |
| "Eq. 3.20" Metropolis f̂ | **3.21** | `f̂(ν) = q(ν) e^{-βν/4} ≈ min{1, e^{-βν/2}}` |
| "Eq. 3.22" Gaussian | **3.22 + 3.23** | `q(ν) = exp(-(βν)²/8) w(ν/S)`; `f̂(ν) ∝ exp(-(βν+1)²/8)` |
| "Eq. 4.7" KMS skew-symm | **4.7** in Section 4 | `α_{ν,ν'} = α_{-ν',-ν} e^{-β(ν+ν')/2}` |

**Not a bug**, but a documentation drift.

### Coherent operator G time-form: paper has SINGLE integral, code has 2D

Per paper Eq. 3.10: `G = ∫ g(t) H_L(t) dt`. Code's `dll_coherent_op_time` and `_dll_coherent_from_g_tt` implement a **2D** integral `Σ_a ∬ g(t,t') A^a(t')A^a(t) dt dt'`. Mathematically equivalent — substituting `L_a = ∫ f^a(t') A^a(t') dt'` into `H_L(t) = Σ_a e^{iHt}(L_a† L_a) e^{-iHt}` gives the 2D form after re-parameterisation. The operator-ordering correction `A^a(t')A^a(t)` (paired-with-t' on LEFT) is correctly derived and implemented; legacy reference at `src/dll.jl:611-625` confirms `mul!(prod_buf, Atn, Atm)` with `nidx → t'` and `m → t`, giving `A(t')·A(t)`.

### Verdict per check (using **thesis** numbering)

| # | Check | Status |
|---|-------|--------|
| 1 | **Eq. 3.4** Lindblad operator: Bohr `freq_kernel = q · e^{-βν/4}`, time form `f(t)` is inverse FT | **PASS** |
| 2 | **Eq. 3.5** coherent kernel `ĝ(ν,ν') = (1/2i) tanh(β(ν'−ν)/4) f̂(ν) conj(f̂(ν'))` | **PASS** (equivalent to paper's 1D Eq. 3.7 via Bohr-form double sum) |
| 3 | **Eq. 3.7** time-domain G: paired-with-t' LEFT order `A^a(t')A^a(t)` | **PASS** (operator order matches across all three implementations: bohr, time-NUFFT, time-legacy) |
| 4 | **Eq. 3.19-3.20** Metropolis `q(ν) = exp(-√(1+(βν)²)/4)·w(ν/S)`, `f̂(ν) = q(ν) e^{-βν/4}` | **PASS** |
| 5 | **Eq. 3.22** Gaussian `f̂(ν) = e^{1/8} exp(-(βν+1)²/8)` | **PASS** (verified algebraically: completing the square `-(βν)²/8 - βν/4 = -(βν+1)²/8 + 1/8`) |
| 6 | **Eq. 4.7** KMS skew-symmetry α(ν,ν') = α(-ν',-ν) e^{-β(ν+ν')/2} | **PASS** (test fixture `test/test_helpers.jl:34-45::assert_kms_skew_symmetric` matches; analytic verification: rank-1 outer product of `freq_kernel` with `q(-ν) = conj(q(ν))` gives the identity) |

### Concerns

- **Documentation drift (low severity)**: code comments use thesis Eq. numbers. Recommend adding a header note per file: "Eq. 3.X (thesis) / Eq. 3.Y (paper Ding-Li-Lin 2024 v4)". No physics fix needed.

- **`S/2 ≥ max|ν_BH|` not validated** (already addressed by C3): `validate_config!` only checks `S > 0`. The `DLLMetropolisFilter` docstring was previously claiming a non-existent warning; commit C3 updated the docstring to "the caller must ensure". Severity: **low** (caller-controlled; default S=2 is correct for the test fixtures with `max|ν_BH| ≈ 0.9`).

- **No correctness fixes required.**

---

## References consulted

- Ding, Li, Lin (2024). "Efficient quantum Gibbs samplers with Kubo–Martin–Schwinger detailed balance condition." arXiv:2404.05998 v4. Pages 15-21 for Sec. 3.1-3.3, page 25 for Sec. 4.
- `.claude-memory/gqsp_bs_vs_mw_convention.md` — BS-vs-MW transpose convention note
- Audit reports `audit-{gqsp,hamiltonian,dll}.md` in `.planning/phases/52-cleanup-gqsp-ham-dll/`

---

## Final summary

- **E1 GQSP**: PASS. Jacobi-Anger/Chebyshev polynomial, α-norm formula, Bessel-tail bound, BS-vs-MW convention all match paper and POC; the new `_coherent_unitary_step` helper from A3 also resolves the pre-existing Hermitian-asymmetry between Thermalize and Trajectory paths.
- **E2 Berntson-Sünderhauf**: PASS. The Julia simulator never touches angles — it computes the closed-form Chebyshev polynomial directly, so the BS/MW transpose convention is a Python-POC-only concern.
- **E3 DLL**: PASS on all six checks. Two low-severity doc concerns: (a) thesis-vs-paper equation-number drift in code comments (already largely consistent with thesis style; physics is correct); (b) `DLLMetropolisFilter` docstring overpromise — already fixed by C3. No correctness fixes required.
