# Physics Review: QuantumFurnace.jl (qf-6z9.5)

**Overall summary**: comprehensive physics validation pass over `src/`, `test/`, and `scripts/`. The codebase is in good physical/mathematical shape — **no factor-of-2 sign or scaling errors that would invalidate downstream claims**. KMS detailed-balance conventions, Trotter splitting orders, gamma-norm closed-form sup, and DLL canonical-form decompositions are all consistent with the source papers. The main flag is a documentation drift (DLL paper Eq. numbers in code/comments are off-by-one vs. arXiv:2404.05998 v4); several `[ambiguous]` items where test fixtures don't match thesis fixtures (subtle but worth being aware of when interpreting test results).

**Counts by severity**:
- **[OK]** confirmations: 22
- **[ambiguous]** (worth noting): 6
- **[non-physical-regime]** (test fixture caveats): 3
- **[wrong-test-fixture]** (real fixture mismatch with thesis convention): 1
- **[formula-error]**, **[sign-error]**, **[wrong-scaling]**: **0**

There are **no outright sign or scaling errors**.

---

## 1. KMS Detailed-Balance Conventions

### 1.1 Skew-symmetry helper
**[OK]** `test/test_helpers.jl:63-74` — `assert_kms_skew_symmetric(α, ν_grid, β; atol=1e-12)` checks
```
α[p, q] = α[neg_idx[q], neg_idx[p]] · exp(-β(ν_p + ν_q)/2)
```
which encodes `α(ν, ν') = α(-ν', -ν) · exp(-β(ν+ν')/2)` (Ding-Li-Lin 2024 Eq. 4.7 ≡ Ramkumar-Soleimanifar Lemma 7.1). Correct.

**[OK]** Used consistently in:
- `test/test_dll_kossakowski.jl:157` (single-channel DLL Gaussian + Metropolis, β ∈ {1,5,10})
- `test/test_dll_multichannel_bohr.jl:61` (multi-channel sum)
- `test/test_dll_multichannel_filter.jl:228` (shifted-symmetric translates)
- `test/test_dll_multichannel_simulator.jl:192` (full simulator path)

### 1.2 Legacy `check_alpha_skew_symmetry`
**[OK]** `src/bohr_domain.jl:171-173` — same convention as `assert_kms_skew_symmetric`.

### 1.3 σ_β fixed-point check
**[OK]** `src/discriminant.jl:370-378` — `verify_detailed_balance` checks `‖D · vec(σ^{1/2})‖ ≈ 0`. Correct because `D(σ^{1/2}) = σ^{-1/4} L(σ) σ^{-1/4} = 0` when `L(σ) = 0`.

**[OK]** Applied consistently in `test/test_dll_dissipator.jl:55-67`, `test/test_dll_kms_db.jl:84-90`, `test/test_dm_detailed_balance.jl:30-31`.

### 1.4 CKG `create_alpha` skew-symmetry (algebraic verification)
**[OK]** `src/bohr_domain.jl:90-104` — under (ν_1,ν_2) → (-ν_2,-ν_1), only `exp(-C)` flips with `C = β(ν_1+ν_2)/4`; ratio = `exp(-2C) = exp(-β(ν_1+ν_2)/2)` ✓.

### 1.5 `create_alpha_gauss` skew-symmetry
**[OK]** `src/bohr_domain.jl:140-152`. Ratio = `exp(-w_γ·(ν_1+ν_2)/(σ²+σ_γ²)) = exp(-β(ν_1+ν_2)/2)` when `2w_γ/(σ²+σ_γ²) = β`. Matches `validate_config!` at `src/misc_tools.jl:269` ✓.

---

## 2. β / σ Scalings

### 2.1 Standard Chen choice σ = 1/β
**[OK]** Used uniformly across test fixtures: `test/test_helpers.jl:81`, `test/test_dll_*.jl`, `src/lindblad_action.jl:875`.

### 2.2 σ ≠ 1/β regression
**[OK]** `test/test_trotter_caches.jl:243-275` — explicitly tests σ=0.05, β=10 (σβ = 0.5 ≠ 1) for the per-leg natural Trotter step formula:
```
b_-(t/σ) → t0_b_minus_evol = register_t0_b_minus / σ
b_+(τβ) → t0_b_plus_evol  = β · register_t0_b_plus
```
Slope -2 in M_user verified for non-default σ. Matches `feedback_qf_d0w_shared_delta_t0.md`.

### 2.3 GNS β = 2ω_γ/σ_γ²
**[OK]** `src/misc_tools.jl:266-267` — different formula (no σ²) for GNS construction. Consistent with `_pick_alpha_gns` and `create_alpha_gns`.

---

## 3. Bohr-frequency vs continuous-ω discretisation

### 3.1 Hermitian-jump half-grid fold
**[OK]** Pattern is consistent across all forward/adjoint EnergyDomain & Time/TrotterDomain matvecs:
- `src/krylov_matvec.jl:131-152` (forward EnergyDomain)
- `src/krylov_matvec.jl:199-225` (adjoint EnergyDomain)
- `src/krylov_matvec.jl:521-547` (forward Time/TrotterDomain)
- `src/krylov_matvec.jl:594-621` (adjoint Time/TrotterDomain)

Hermitian branch: only `w_raw <= 1e-12` queued, conjugate-transpose sandwich for negative-frequency partner. Non-Hermitian branch: full ω-grid.

### 3.2 Threaded-ω work-list builder
**[OK]** `src/krylov_matvec.jl:640-659` — `_populate_lindblad_work_list!` mirrors the serial pattern. Per-thread chunks correctly use `w = abs(w_raw)` for Hermitian and signed `w = w_raw` for non-Hermitian (lines 732, 828). Consistent with `feedback_non_hermitian_jumps.md`.

### 3.3 jump_workers parallel pattern
**[OK]** `src/jump_workers.jl:68-87, 217-238, 423-465` (CKG TimeDomain, EnergyDomain, BohrDomain Thermalize). Same Hermitian-fold convention.

### 3.4 Thermalize Energy/Time uses (Aw)† for negative branch
**[OK]** `src/jump_workers.jl:441-450` — for Hermitian jumps, the negative-frequency branch uses `mul!(scratch.LdagL, scratch.jump_oft, scratch.jump_oft')` (i.e. `Aw · Aw†`) and `evolving_dm * jump_oft` (i.e. `ρ · Aw`), reflecting `A_{-w} = (A_w)†`. Consistent.

---

## 4. gamma-family normalisation (qf-etx)

### 4.1 `pick_gamma_sup(config) = 1.0` closed form
**[OK]** `src/energy_domain.jl:149-150`:
```julia
pick_gamma_sup(config::Config{<:Any, <:Any, KMS}) = 1.0
pick_gamma_sup(config::Config{<:Any, <:Any, GNS}) = 1.0
```
Used in 4 `_precompute_data` overloads (`src/furnace_utensils.jl:46, 60, 96, 192`).

### 4.2 No remaining grid-dependent `1/maximum(transition.(...))` callers
**[OK]** All `gamma_norm_factor` consumers route through `pick_gamma_sup`.

### 4.3 No `/gnf` workarounds remaining in scripts
**[OK]** `scripts/scratch_*_ref_convergence.jl`, `scripts/scratch_quadrature_n*_bench.jl`, etc. all updated with the post-qf-etx.2 comment "gnf is grid-independent so plain L_test - L_ref is the error".

### 4.4 Test coverage
**[OK]** `test/test_gamma_norm_invariance.jl:19-130` — 110 assertions cover KMS Gaussian / kinky Metro / a-regularized / smooth Metro and GNS variants, all giving `pick_gamma_sup = 1.0`.

### 4.5 Caveat for (a > 0, s > 0) case
**[wrong-test-fixture]** `test/test_helpers.jl:207-208` — `make_config` test fixture uses `a = β/30 = 1/3, s = 0.4`. Per `feedback_gamma_sup_smooth_a0.md`, this `(a > 0, s > 0)` case has continuum sup ≈ 0.74 < 1 (genuinely under-amplified by `gnf = 1.0`). The thesis fixture is `a = 0, s = 0.25` (in-principle continuum sup = 1, in practice exactly 1 in Float64). Most non-DLL tests use this off-thesis `(a, s) = (1/3, 0.4)` fixture; under-amplification is consistent across them and they're internally self-consistent, but absolute γ scales should not be cross-compared with thesis numerics until this is unified. Documented user decision per memory.

### 4.6 In-principle smooth-Metro a=0 caveat
**[ambiguous]** For `(a=0, s>0)` smooth Metropolis (thesis fixture `s=0.25`), the continuum sup = 1 only in the ω → -∞ limit. On any finite grid, sup < 1 in principle, but Float64 arithmetic hits exactly 1.0. `pick_gamma_sup = 1.0` is correct in practice but a tiny under-amplification remains in principle. Test `test/test_gamma_norm_invariance.jl:96-101` evaluates at ω = -10 to confirm sup ≥ 0.999.

---

## 5. Coherent term for non-Hermitian jumps

### 5.1 Production jump sets are all Hermitian (X, Y, Z Paulis)
**[OK]** `src/lindblad_action.jl:874-888` and `test/test_helpers.jl:140-150` use `paulis = ([X], [Y], [Z])`; `(A, A†)` pair requirement of `feedback_non_hermitian_jumps.md` automatically satisfied.

### 5.2 Single non-Hermitian jump in tests
**[OK]** `test/test_krylov_matvec.jl:475` — admissible per the third constraint in `feedback_non_hermitian_jumps.md`: only used to compare two evaluations of the same physics (dense vs Krylov matvec); does NOT claim KMS DBC.

### 5.3 CKG B = α · A^†(ω_1) A(ω_2) for Hermitian jumps
**[OK]** `src/bohr_domain.jl:13-27` correctly assembles `Σ A^†(ν_2) f A(ν_1)`. Valid because Hermitian A satisfies `A_w† = A_{-w}`.

### 5.4 No production code uses non-Hermitian jumps
**[OK]** The "extra term in CKG paper for non-Hermitian jumps" is a known limitation but not currently exercised in production.

---

## 6. DLL canonical-form / Theorem 10 multi-channel decomposition

### 6.1 `dll_multichannel_translates` factory
**[OK]** `src/dll_multichannel.jl:344-388` — builds `DLLMultiChannelFilter` via `ShiftedSymmetricFilter` channels. The √(w/2) factor preserves rank-1 outer-product KMS skew-symmetry per channel.

### 6.2 G^multi = Σ_ℓ G^(ℓ) (linearity)
**[OK]** `src/dll_multichannel.jl:457-468` (Bohr) and `:484-499` (Time). Verified at `test/test_dll_multichannel_bohr.jl:84-110`.

### 6.3 Per-channel dissipator sum (no cross-terms)
**[OK]** `src/dll_multichannel.jl:413-420` and `:436-445` return `Vector{Matrix}` (no L^(i)ρL^(j)† cross-terms). Verified at `test/test_dll_multichannel_bohr.jl:38-41`.

### 6.4 DLL Kossakowski rank-1
**[OK]** `src/dll.jl:649-682` `dll_kossakowski_bohr` returns `v · v†`. `test/test_dll_kossakowski.jl:24-36` verifies rank-1, Hermitian, PSD.

### 6.5 Multi-channel α^multi KMS skew-symmetric
**[OK]** `test/test_dll_multichannel_bohr.jl:48-63`, `test/test_dll_multichannel_simulator.jl:183-194`.

### 6.6 Documentation drift: DLL paper Eq. numbers
**[ambiguous]** `.planning/phases/52-cleanup-gqsp-ham-dll/physics-check.md:60-78` documents that code/comments cite **thesis** Eq. numbers, not paper Eq. numbers (arXiv:2404.05998 v4). Examples:
- "Eq. 3.4 first form" → paper's Eq. 3.3
- "Eq. 3.5" → paper's Eq. 3.7
- "Eq. 3.7" → paper's Eq. 3.10
- "Eq. 3.19" Metropolis q → paper's Eq. 3.20
- "Eq. 3.22" Gaussian → paper's Eq. 3.22 + 3.23
- "Eq. 4.7" KMS skew-symm → paper's Eq. 4.7 (matches)

**Not a physics bug** but a documentation hazard. Suggest a per-file header note ("Eq. 3.X (thesis) ≡ Eq. 3.Y (paper Ding-Li-Lin v4)") in `src/dll.jl`, `src/filters.jl`, `src/dll_multichannel.jl`. **Already noted in physics-check.md but not yet fixed.**

### 6.7 DLL coherent G operator-ordering A^a(t')A^a(t) (paired-with-t' LEFT)
**[OK]** `src/dll.jl:118-128` documents the operator-ordering correction. Verified at `src/dll.jl:606-625` (legacy reference). Numerical Bohr↔Time agreement at ~1e-10 op-norm validates this.

### 6.8 DLL TimeDomain `t0` matches Riemann sum weight `τ`
**[OK]** `src/dll.jl:81-103` — `weight = ft * t0`. Comment line 67 confirms.

---

## 7. Trotter splitting orders

### 7.1 Shared-δt₀ scheme (qf-d0w)
**[OK]** `src/trotter_domain.jl:99-159` — `TrottTrott(ham, t0_D, t0_b_minus, t0_b_plus, M_user)`: picks `δt₀ = min(t0_X) / M_user`, asserts `t0_X / δt₀ ∈ ℤ`, single Strang at δt₀, vector powers `eigvals_t0_X = λ_S^M_X`.

### 7.2 `make_trotter_for_config` per-leg natural step
**[OK]** `src/misc_tools.jl:228-246` — generalized for σ ≠ 1/β.

### 7.3 B_trotter slope -2 in M_user
**[OK]** `test/test_trotter_caches.jl:119-158` — error ratio between consecutive M_user doublings is > 3 (Strang slope -2 ≈ ratio 4). At M_user=8: error/norm < 1e-5.

### 7.4 Legacy single-cache regression
**[OK]** `test/test_trotter_caches.jl:160-188`: legacy `TrottTrott(ham, t0_D, M)` saturates — confirms the bug fixed by qf-d0w.

### 7.5 σ ≠ 1/β slope -2 regression
**[OK]** `test/test_trotter_caches.jl:243-275`.

### 7.6 Shared eigenbasis sanity
**[OK]** `test/test_trotter_caches.jl:103-117`.

### 7.7 Generator splitting slope +2 in δ
**[OK]** `error_analysis_qf_b4d.md` records 1.985-1.997 measured slope.

---

## 8. Quadrature register conventions (qf-7xt)

### 8.1 Per-register triples on Config
**[OK]** `src/structs.jl:131-150` — independent `(num_energy_bits_X, t0_X, w0_X)` for X ∈ {D, b_minus, b_plus}; legacy auto-promotion via `register_*_X`.

### 8.2 Fourier-relation validation per register
**[OK]** `src/misc_tools.jl:403-455`. DLL TimeDomain (`:463-468`) only requires `(r_D, t0_D)`.

### 8.3 t0 · w0 = 2π/N convention
**[OK]** `test/test_helpers.jl:101` — `T0 = 2pi / (2^NUM_ENERGY_BITS * W0)`.

### 8.4 "Channel" = Trotter+GQSP convention
**[OK]** `feedback_channel_means_trotter_gqsp.md`. Production tests at `test/test_predict_channel.jl` and `test/test_gqsp_thermalize.jl` confirm.

### 8.5 Quadrature recipe (qf-7xt)
**[OK]** Recommended bit budgets at ε=1e-6: smooth Metro `r_D=5, r_-=6, r_+=14`; kinky `r_D=12, r_-=6, r_+=14`; Gaussian `r_D=5, r_-=6, r_+=6`. Test fixtures use r=12 across the board (conservative).

---

## 9. Trap-rule t=0 L'Hôpital

### 9.1 `_compute_b_plus_metro` t=0 limit
**[OK]** `src/coherent.jl:452-462`:
```julia
if abs(t) < 1e-12
    return complex((2 - sigma^2 * beta^2 * (1 + s)) / (2 * sqrt(2) * pi^2))
```
At σβ=1, s=0: 1/(2√2 π²) ✓. Verified at `test/test_smooth_metro_eta.jl:33-50`.

### 9.2 η < t0' = 2T_+/2^r_+ convention
**[OK]** Per `trap_rule_t0_lhopital_origin.md`: choice of `η < t0'` makes the η-cutoff branch dead code. Slope-(-1) `‖B_bohr − B_time‖_op` error dominated by t=0 L'Hôpital sample. **No code bug** — explained in memory.

---

## 10. Additional checks

### 10.1 GQSP block-encoding norm α_be
**[OK]** `src/coherent.jl:343-355`. Tested at `test/test_gqsp_polynomial.jl:189-197`.

### 10.2 GQSP polynomial f_d(B/α) Jacobi-Anger
**[OK]** `src/coherent.jl:379-439` — `a_k = 2 (-i)^k J_k(δα)` matches Jacobi-Anger truncation. d=1 fast path: `f_1 = J_0 I − 2i J_1 (B/α)`.

### 10.3 CPTP weak-measurement K0, U_residual
**[OK]** `src/furnace_utensils.jl:275-293`. Matches Chen et al. Eq. 3.2.

### 10.4 Discriminant materialisation
**[OK]** `src/discriminant.jl:165-186` — equivalent to `D = (σ^{-1/4} ⊗ σ^{-1/4}) · L · (σ^{1/4} ⊗ σ^{1/4})` without forming Kronecker products. Verified at `test/test_discriminant.jl:88-92`.

### 10.5 OFT sign convention
**[OK]** `src/structs.jl:65` documents `Σ_t̄ b̄(t̄) e^{-iωt̄} A(t̄)`. NUFFT prefactor uses `isign = +1` so `out = Σ input · exp(+i · sx · tx)` — combined with `input_weights = base_weights · exp(-iω·t)` gives the right sign.

### 10.6 DLL Gaussian-type Gevrey filter
**[OK]** `src/filters.jl:130-162`. Verified at `test/test_dll_filter.jl:21-65`.

### 10.7 Hörmander mollifier
**[OK]** `src/filters.jl:206-249` — even, smooth, =1 on |x| ≤ 1/2, =0 on |x| ≥ 1.

### 10.8 Workspace structure for KrylovScratch ω-loop threading
**[OK]** `src/structs.jl:423-454` — pre-allocated `task_scratches` pool (size = nthreads).

### 10.9 Predict-trajectory Krylov (qf-ev5)
**[OK]** `src/lindblad_action.jl:595-662` (Lindbladian) and `:731-859` (channel). Channel reference at lines 712-730: trace-distance reference is Gibbs-in-Trotter-eigenbasis for TrotterDomain (qf-ev5.8 fix).

### 10.10 Bohr-dict diagonal element at ν=0
**[OK]** `src/bohr_domain.jl:159-161`.

### 10.11 Sign flip in CKG `f`
**[OK]** `src/bohr_domain.jl:44-47`. `tanh(-β(ν_1-ν_2)/4)` anti-symmetric, `α(ν_1, ν_2)` symmetric, combined with `1/(2i)`: `Σ A^†(ν_2) f A(ν_1)` is Hermitian ✓.

---

## 11. Minor observations / noting

### 11.1 `test/test_helpers.jl` test fixture vs thesis fixture
**[wrong-test-fixture]** As noted in 4.5, the default test fixture uses `(a = β/30, s = 0.4)` not the thesis `(a = 0, s = 0.25)`. Deliberate user decision per `feedback_gamma_sup_smooth_a0.md`.

### 11.2 NUFFT path & gnf=1.0 cross-grid invariance
**[OK]** `test/test_gamma_norm_invariance.jl:135-460` — 110 assertions.

### 11.3 DLL TimeDomain: no `w0_D`
**[OK]** `src/misc_tools.jl:463-468` — DLL TimeDomain only requires `(r_D, t0_D)`.

### 11.4 DLL TrotterDomain & EnergyDomain: deferred
**[OK]** `src/misc_tools.jl:362-378` — `validate_config!` rejects DLL with EnergyDomain or TrotterDomain.

### 11.5 Smooth Metropolis s=0.4 fixture vs thesis s=0.25
**[non-physical-regime]** Difference is small (~3-6% τ_mix penalty per `smooth_metropolis_s_decision_qf_3il.md`); tests internally consistent.

### 11.6 Test `test/test_dm_detailed_balance.jl` n=3 fixture uses `make_config` defaults
**[non-physical-regime]** Uses `(a = β/30, s = 0.4)` — KMS DBC structurally exact, passes regardless of `s`.

### 11.7 σ_β = 1/β default
**[non-physical-regime]** Filter band ~few · σ overlaps Bohr spectrum nontrivially. Healthy for stress-testing.

### 11.8 GQSP simulator tests use higher s
**[OK]** Tests verify GQSP polynomial truncation error vs matrix-exp baseline; independent of γ choice.

### 11.9 `test/test_jump_selection.jl` rate scaling
**[OK]** Both `:sweep` and `:random` reproduce `e^{T𝓛}` in expectation.

### 11.10 `test/test_dll_dissipator.jl` r_D choice
**[ambiguous]** Uses `_DLL_NUM_ENERGY_BITS = 12`, conservative for thesis recipe (r_D=5 for Gaussian, r_D=12 for kinky). Higher r_D ≠ wrong — just over-allocates; tests pass.

---

## 12. Summary — no observed sign errors, factor-of-2 errors, or wrong scaling laws

Across all the above checks, **no formula errors, sign errors, or scaling errors that would invalidate downstream physics claims** were found. The codebase is internally consistent with the cited papers (DLL 2024 v4, Chen-Kastoryano-Gilyén, Chen-Kastoryano-Brandao-Gilyén) and the thesis conventions.

The most pressing follow-ups (priority order):

1. **[ambiguous]** Documentation drift for DLL paper Eq. numbers (item 6.6). Suggest per-file header note in `src/dll.jl`, `src/filters.jl`, `src/dll_multichannel.jl`. Already noted in `.planning/phases/52-cleanup-gqsp-ham-dll/physics-check.md` from cleanup epic, not yet applied.
2. **[wrong-test-fixture]** The default test fixture `make_config` (`a = β/30, s = 0.4`) is NOT the thesis fixture (`a = 0, s = 0.25`). Either (a) explicitly document the discrepancy in `test/test_helpers.jl`, or (b) migrate to `a = 0, s = 0.25` and re-baseline regression numbers. The DLL test files are clean.
3. **[ambiguous]** Multiple test files don't match the qf-7xt quadrature register recipe; they over-allocate (r_D=12 across the board). Not wrong, just non-minimal. No action needed unless test wall time becomes a concern.
