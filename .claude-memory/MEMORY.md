# QuantumFurnace.jl Memory

## Project Defaults

### Thesis scope: finite-dim only — infinite-dim mentioned but never prioritised (2026-05-14)
- See [Thesis scope finite-dim only](feedback_thesis_scope_finite_dim_only.md). The thesis works exclusively with finite-dimensional Hilbert spaces (spins, fermions on a lattice with bounded local Fock dim, qubits/qudits). Infinite-dim / bosonic / continuous-variable Gibbs samplers (Becker 2026, Bose–Hubbard, oscillator chains) are out of scope. In any synthesis, frontier ranking, or chapter focus: mention infinite-dim once with an `[OUT-OF-SCOPE: infinite-dim]` tag, but never place it in a top-N "where quantum could first beat classical" list. The thesis machinery (KMS-DB construction, Trotter/energy/Bohr quadrature ladders, Krylov spectral analysis, numerics chapter) is built and tested only for finite-dim — infinite-dim needs a different functional-analytic framework the thesis does not develop.

### 1D Heisenberg β_phys∈[0.25,2.5] τ_mix is M0 not M1 (qf-e4z.22, 2026-05-15) — **REVISED by qf-e4z.23**
- See [Heisenberg 1D no-Arrhenius qf-e4z.22](heisenberg_1d_no_arrhenius_qf_e4z_22.md). Extended v4 sweep to factor-10 β_phys range, n∈{3..8}. M0 (separable power) decisively beats M1 (Arrhenius): AICc weight 0.997 (CKG) / 0.991 (DLL), Δ_AICc ≈ 10. Fit: τ_mix ≈ 7.5·n^1.23·β_phys^0.32. n-exponent matches R(n)~n^1.24 algorithm-frame rescaling, NOT intrinsic mixing scaling. β-exponent 0.32 is artefactual saturation (τ_mix flat above β_phys≈1). gap_phys stays in [2.27, 6.45] across all 36 cells — Ω(1) thermodynamic gap, no collapse. Consistent with KB16/Bardet23/BC25. CKG ~30% faster than DLL at every cell. r_D=8 needed for ≤1e-9 vs Bohr at β_alg≳50 (v5 cells); v4 cells at r_D=7 fine. **Caveats from [[heisenberg-1d-multiseed-even-odd-qf-e4z-23]]**: this fit was single-seed `find_typical`; 5-seed disorder averaging exposes a real even/odd-n parity split with even-n gap_phys collapsing at high β. The M0 fit becomes marginal (weight 0.69, Δ_AICc 1.6) and both M0/M1 are misspecified. Use the parity-resolved characterisation, not the pooled fit.

### 1D Heisenberg PBC v6 5-seed multiseed: REAL even/odd-n parity split (qf-e4z.23, 2026-05-16)
- See [Heisenberg 1D multiseed even/odd qf-e4z.23](heisenberg_1d_multiseed_even_odd_qf_e4z_23.md). 180-cell PBC sweep (n=3..8, β_phys=0.25..2.5, 5 seeds) on `build_heis_1d` multi-seed fixtures. Even-n gap_phys COLLAPSES at high β: n=8, β_phys=2.5 → 1.08; odd-n n=7 stays at 5.00 at same β. Intra-parity 5-seed half-range ~1% << inter-parity log-spread 0.11→0.89 across β. v5 single-seed scaling fit (τ_mix ~ n^1.23·β_phys^0.32, M0 weight 0.997) is **not robust** — v6 median fit gives n^1.89·β_phys^0.32 with M0 weight only 0.69. Both single-model fits misspecified; data needs parity-resolved fit. β-exponent 0.32 stays robust. Even-n τ_mix grows steeply with β at n≥6; odd-n + n=4 saturate above β_phys≈0.5. Methodology: 3 even-n + 3 odd-n is on the edge per [[feedback_more_data_points_for_scaling_claims]] — we can claim the structure exists with 1% noise floor and monotone β-growth, NOT asymptotic exponents. Cluster follow-ups: OBC counterpart + n∈{10,12,14} both filed.

### Krylov x_0 must break symmetry; gap_alg ≠ gap_phys (qf-8fr, 2026-05-13) — **EXTENDED by qf-e4z.26**
- See [Krylov x_0 symmetric-system bug qf-8fr](krylov_x0_symmetric_bug_qf_8fr.md). Two findings: (1) pre-patch `krylov_spectral_gap` seeded Arnoldi with `vec(I/d)` — symmetry-protected on clean Hamiltonians, silently returned wrong gap (3.8× overestimate on classical Ising n=4). Patched: `_krylov_default_x0` builds `vec(I/d + 1e-10·H_GUE_traceless)`. (2) Reported λ_gap is `λ_alg` in rescaled frame H_alg = H_phys/R + s·I, R∝n for extensive H — true Ω(1) physical gap shows up as ~1/n in algorithm units; conversion λ_phys = λ_alg · ham.rescaling_factor. **CORRECTION** (qf-e4z.26): the original write-up claimed disordered Heisenberg was unaffected — this is WRONG, see [[predict-trajectory-parity-bug-qf-e4z-26]]: [[Z],[Z,Z]] disorder preserves Z^⊗N parity, so the trap re-triggered in `predict_lindbladian_trajectory` (which the qf-8fr patch did NOT cover). 

### predict_lindbladian_trajectory had the parity-trap bug too (qf-e4z.26, 2026-05-16) — **SUPERSEDED by qf-e4z.27**
- See [predict_trajectory parity bug qf-e4z.26](predict_trajectory_parity_bug_qf_e4z_26.md). qf-e4z.26 1e-6-GUE-perturbed seed in `_krylov_spectral_decomposition` caught the bug at n ≤ 5 but FAILED at n=7 (gap 0.1203 vs true 0.0966). Proper fix in qf-e4z.27 below.

### qf-e4z.27 two-pass Krylov fix for predict_*_trajectory (2026-05-16)
- See [predict_trajectory two-pass qf-e4z.27](predict_trajectory_two_pass_qf_e4z_27.md). `predict_lindbladian_trajectory` / `predict_channel_trajectory` now use TWO Krylov passes: (1) single-seed `_arnoldi_factorize` from `vec(rho_0)` (trajectory tightly aligned, rho_0 ∈ span(Q) exactly, eigenvalues reflect rho_0's symmetry sector); (2) separate `krylov_spectral_gap` call with `_krylov_default_x0` + KrylovKit thick restart for the true Lindbladian gap. Returned `spectral_gap` field comes from pass 2. Spot test n=7 β_phys=2.5 seed=42 (Z+ZZ legacy): rel_err = 0.0 (was 24% pre-fix). Cost +50–100 matvecs. Dual-seed Arnoldi tested+rejected (splits Krylov budget). KrylovKit eigsolve tested+rejected (Gram-inverse breaks for KMS-DB degenerate spectra). Test (f) rtol 1e-6→1e-8; new (d) channel regression. **τ_mix unchanged** (parity-odd c_i ≡ 0 for parity-even rho_0). Sidecars from qf-e4z.{20-25} are pre-fix on `gap_arnoldi`; τ_mix data stands.

### qf-e4z.27 XZ-disorder fixtures: parity-trap-free Heisenberg (2026-05-16) — **scope narrowed by qf-e4z.28**
- See [XZ disorder fixtures qf-e4z.27](xz_disorder_fixtures_qf_e4z_27.md). New 1D Heisenberg fixture family with `[[X], [Z, Z]]` disorder at `heis_xxx_XZdisordered_periodic_n{n}_seed{seed}.bson`, parallel to legacy `[[Z], [Z, Z]]`. X-disorder anti-commutes with `P = Z^⊗N`, structurally breaking parity symmetry at Hamiltonian level (verification: `‖[P, H]‖ > 0` on X+ZZ, `= 0` exactly on Z+ZZ). Generator: `run_1d_heisenberg_XZ_disorder()` in `scripts/scratch_multiseed_disordered_fixtures.jl`. **Scope correction from qf-e4z.28**: the X-disorder strength (0.1) is too weak to eliminate the algorithmic parity trap inside plain MGS Arnoldi at krylovdim=40 from `vec(I/d)`. X+ZZ is the right *physics* fixture (no even/odd-n parity artefact in qf-e4z.23-style sweeps), but it does NOT remove the need for the qf-e4z.27 two-pass in trajectory predictors. See [[single-pass-xz-disorder-decorative-qf-e4z-28]].

### qf-e4z.28 single-pass on X+ZZ: structural fix decorative; keep two-pass (2026-05-16) — **SUPERSEDED by qf-e4z.30**
- See [Single-pass X+ZZ decorative qf-e4z.28](single_pass_xz_disorder_decorative_qf_e4z_28.md). Tested whether the qf-e4z.27 X+ZZ fixtures structurally eliminate the parity-trap inside `predict_lindbladian_trajectory` such that Pass 2 could be dropped. **They do NOT when paired with `rho_0 = I/d`.** Spot cells n=5, n=7 β_phys=2.5 seed=42: single-pass Arnoldi from `vec(I/d)` gives the SAME wrong gap on X+ZZ as Z+ZZ (rel_err 16.6% / 24.3% respectively) — the 0.1 X-disorder is below the MGS orthogonalisation threshold at krylovdim=40. Test E confirmed: 1e-10 GUE-perturbed seed FAILS identically in plain Arnoldi. **Conclusion was "keep two-pass"** — superseded by qf-e4z.30 which solved the problem from the INPUT side with `rho_0 = |+⟩⟨+|^⊗N` instead of changing the seed. Driver: `scripts/scratch_qf_e4z_28_single_pass.jl`.

### qf-e4z.29 + qf-e4z.30 canonical setup: |+⟩⟨+|^⊗N + single-pass kdim=60 (2026-05-16)
- See [Canonical τ_mix setup qf-e4z.30](canonical_taumix_setup_qf_e4z_30.md). Two empirical resolutions: (1) qf-e4z.29 at n=5 with dense `exp(L·t)·vec(rho_0)` reference: `rho_0 = |+⟩⟨+|^⊗N` (single Hadamard layer on `|0⟩^⊗N`) gives trajectory err_max = 3.9e-14 AND gap_traj = gap_dense EXACTLY at krylovdim=40 — machine precision on BOTH Z+ZZ and X+ZZ fixtures. (Compare: `rho_0 = I/d` gives correct trajectory to 1e-9 — verifying the [[predict-trajectory-two-pass-qf-e4z-27]] symmetry argument — but gap_traj is 17% off, that's the parity-even sub-spectrum.) (2) qf-e4z.30 at n=7 krylovdim sweep ∈ {40, 60, 80, 100, 120}: `|+⟩⟨+|^⊗N` reaches rel_err 2.4e-9 (X+ZZ) / 2.0e-9 (Z+ZZ) at krylovdim=60 (single-pass, 60 matvecs vs ~166 for qf-e4z.27 two-pass; 2.75× faster, 4 orders of magnitude more accurate). Fixture choice (X+ZZ vs Z+ZZ) irrelevant once `rho_0` is parity-broken. **Canonical thesis recipe** for 1D Heisenberg τ_mix / gap plots: Z+ZZ-disordered fixtures (legacy canonical) + `rho_0 = |+⟩⟨+|^⊗N` + `predict_lindbladian_trajectory` single-pass at krylovdim=60 + odd-n preferred (n ∈ {3, 5, 7, 9}); even-n kept open until larger-n data lands and we know whether the qf-e4z.23 even/odd splitting is real physics or a parity-sub-spectrum artefact from the I/d convention. krylovdim=60 verified through n=7; **must be re-verified at n ≥ 10** before signing off on plots there. X+ZZ fixture family becomes unnecessary for new work (preserved on disk for audit). Drivers: `scripts/scratch_qf_e4z_29_trajectory_vs_dense.jl`, `scripts/scratch_qf_e4z_30_kdim_sweep_n7.jl`.

### 1D fixtures: find_typical_heisenberg + [[Z],[Z,Z]], n=3..9 (qf-2kd, 2026-05-12)
- See [Fixture migration to find_typical](fixture_migration_find_typical_qf_2kd.md). Canonical 1D family is `heis_xxx_zzdisordered_periodic_n{3..9}.bson` (find_typical, batch_size=256). Legacy + clean families removed; F3 2D fixtures still find_ideal (deferred). β_phys path unchanged (raw.rescaling_factor). Smaller nu_min than find_ideal — recipe in `[[Quadrature register recipe qf-7xt]]` may shift 1–2 bits at n≥6.

### Thesis target precision: split — ε=1e-6 algorithm-level, ε≈1e-5 simulated-algorithm (revised 2026-05-04)
- See [Thesis target precision split](thesis_target_precision_1e6.md). The prior "ε=10⁻⁶ everywhere" framing is **scratched**. Two regimes: (1) **Algorithm-level** Lindbladian fixed-point analysis (Krylov, no δ scheme) → ε = 10⁻⁶. (2) **Simulated-algorithm** plots (Thermalize / Trajectory CPTP Kraus scheme) are floored at ε ≈ δ²/λ ≈ 10⁻⁵ since δ ≥ 10⁻³. Always clarify which regime when citing precision targets. See also [δ and t0 floors](project_delta_t0_floor.md).

### Krylov for fixed point / spectral gap, never full eigendecomp
- See [Krylov for fixed-point + gap](feedback_krylov_for_fixed_point_and_gap.md). Use `krylov_spectral_gap(config, ham, jumps; ...)` from `src/krylov_eigsolve.jl`; returns `fixed_point` and `spectral_gap` from matvec `eigsolve` (smallest |Re(λ)| pair). Never call `eigen()` on a dense d²×d² superoperator for fixed-point / gap analysis — only build the dense superoperator when genuinely needed (norm comparisons, KMS geometry).

### BohrDomain gap cross-check at n ≤ 6: dense `construct_lindbladian + eigvals`, NOT matrix-free Krylov (2026-05-14)
- See [BohrDomain gap dense at small n](feedback_bohr_gap_dense_at_small_n.md). Matrix-free Krylov in BohrDomain costs d⁴·n_jumps PER matvec × ~krylovdim+restarts matvecs (15+ min at n=6, 2D TFIM h=J=1). Dense `construct_lindbladian` does d⁴·n_jumps ONCE + LAPACK `geev` (~3 min at n=6, 4096×4096 complex matrix). Same gap to 8.2e-14 rel vs EnergyDomain Krylov. ~10× speedup. Reference helper `bohr_gap_dense` in `scripts/scratch_2d_tfim_gap_sweep.jl`. Forbid n > 6 (4 GB at n=7, 64 GB at n=8). Refines [[Krylov for fixed-point + gap]] for the specific Bohr-cross-check use case. **EnergyDomain still uses matrix-free Krylov at any n** (its matvec is d²·2^r_D, much cheaper).

### 2D TFIM phase diagram + ordered-phase gap scaling (2026-05-14)
- See [T_c(h) for 2D TFIM + ordered gap scaling](tc_2d_tfim_phase_diagram.md). Square-lattice T_c(h) from Hesselmann–Wessel 2016 (arXiv:1602.02096) SSE-QMC: T_c(0) = 2.269 (Onsager), T_c(0.5) ≈ 2.20, T_c(1.0) ≈ 2.07, T_c(1.5) ≈ 1.85, T_c(2.5) = 1.27369(5), T_c(3.0) = 0.2977(9), h_c = 3.04438(2). The qf-833 operating point (h=1, β_phys=0.5 → T=2.0) was MARGINALLY inside the ordered phase (~3% below T_c(1)), NOT disordered as initially labelled — pre-asymptotic finite-size flat gap, not Ω(1) thermodynamic. For clear ordered-phase sweeps: use **(h=1.5, β_phys=2.0)** primary (T/T_c≈0.27) or **(h=1.0, β_phys=2.0)** (matches qf-833 h). Ordered-phase Lindbladian gap_phys ~ exp(−γ(β,h)·√n) — stretched-exponential, same exponent as classical Glauber on 2D Ising (Lubetzky–Sly 2013); Gamarnik–Kiani–Zlokapa 2024 prove T_mix ≥ 2^{n^{1/2−o(1)}} for CKG sampler on 2D TFIM. PRE-ASYMPTOTIC at sandbox n ≤ 8 (√n range 0.83) — needs n ∈ {4,6,8,10,12} on the cluster to even start resolving the exponent. **Z₂-doublet caveat**: in symmetry-broken phase the steady-state manifold is rank-2; Krylov gap with the qf-8fr 1e-10 GUE-traceless seed locks onto the doublet-splitting eigenvalue (= correct τ_mix bottleneck), but must be reported as such, not as "spectral gap of L". Diagnose by checking |λ_3/λ_2| ≥ 5 (clean doublet) vs <5 (contamination → fire qf-833 disorder=1e-3 contingency).

### 2D TFIM ordered-vs-disordered side-by-side (qf-1jj, 2026-05-14)
- See [2D TFIM ordered vs disordered qf-1jj](2d_tfim_ordered_vs_disordered_qf_1jj.md). Two operating points at h=1: ORDERED (β_phys=2, T/T_c≈0.24) and DISORDERED (β_phys=0.25, T/T_c≈1.93). gap_phys n=4→6→8: ORD 9.96e-2 → 3.06e-3 → 1.57e-4 (collapsing ~30×/step); DIS 1.910 → 1.575 → 1.741 (flat, Ω(1)). |λ_3/λ_2|^ORD: 57 → 1100 → 18000 (clean Z₂-doublet capture, qf-8fr seed enough; disorder contingency NOT needed). |λ_3/λ_2|^DIS: 3.58 → 3.16 → 2.53 (no doublet). τ_mix(ε=1e-3) ratio ORD/DIS only 1.4→2.1, NOT exploding with the gap, because from I/d initial state the symmetric mixture has zero overlap with the Z₂-odd doublet — τ_mix tracks 1/λ_3, NOT 1/λ_2. Energy↔Bohr 1e-7 gate: PASS for DIS at machine precision, soft FAIL for ORD (rel ~1e-6 to 1e-4) traceable to abs |Δλ_2|~1e-8 within r_D=7 quadrature floor — needs r_D=8 for cluster ordered cells. Subtask B' dense-Gibbs diagnostic confirmed phase membership before the sweep ran. Driver: `scripts/scratch_2d_tfim_ordered_vs_disordered_sweep.jl` + `..._catchup.jl`. Findings: `drafts/2d-tfim-ordered-vs-disordered.md`. **Caveat**: 3 cells on a 2×L_y thin ladder are qualitative only — does NOT pin down asymptotic scaling form (see [[feedback-more-data-points-for-scaling-claims]]). Cluster plan revised: bulk-2D 3×3 + 3×4 at BOTH operating points (5+ cells across 2 geometry classes needed for asymptotic claims), r_D=8 in ordered.

### Sparse sweeps are qualitative — never fit asymptotic scaling from ≤3 cells (set 2026-05-14)
- See [more data points for scaling claims](feedback_more_data_points_for_scaling_claims.md). Three points along an axis (n, β, h, k, …) cannot discriminate among exp(-c√n), exp(-cn), exp(-cL_y), n^{-p}, finite-size convergence to a floor. Report qualitative trends + per-step ratios; do NOT extract scaling exponents γ. Bigger sweeps (≥5 cells, ≥2 geometry classes for system-size axes) needed for asymptotic claims. Came up while interpreting qf-1jj 2D TFIM ordered-phase gap collapse — initial draft fitted γ≈7.8 in λ ~ exp(-γ√n) from 3 ladder cells, but on 2×L_y ladders the asymptotic theory predicts a *finite floor* (cheapest tunnelling cut is 2 bonds, L_y-independent), not √n decay. Onsager 2D Ising σ at βJ=2 is ≈2, the fitted γ was 4× larger — geometry mismatch. Generalises to: any qf-mto β-sweep over {5,10,20}, DLL k-sweep claims, β_phys "flat" claims at canonical 3-point grid {0.25, 0.5, 1.0}.

### Coherent term ON by default; flag dissipator-only runs explicitly (2026-05-12)
- See [Coherent term on by default](feedback_coherent_term_on_by_default.md). All Lindbladian APIs default to `include_coherent=true` — the physical KMS Lindbladian. Setting `include_coherent=false` yields a dissipator-only operator (non-physical, no KMS DB). Only use it when explicitly isolating dissipator quadrature, and always call out the limitation in script header, output naming, and summary drafts — never present dissipator-only numbers as full-Lindbladian register sizing. Rule encoded in `.claude/rules/julia-code.md` (Coherent Term section) and `.claude/rules/scripts.md`.

### ε must be below the spectral gap for precise Krylov gap/τ_mix (2026-05-12)
- See [ε must be below gap](project_epsilon_must_be_below_gap.md). The "ε" label on register tables = total ‖ΔL‖_op (sum of dissipator + coherent + FINUFFT + Trotter + δ-step contributions). Krylov gap error bounded by κ·‖ΔL‖_op; relative τ_mix error blows up when gap is small. Rule: pick registers so ε ≤ 0.01·min|λ_1| for 1% τ_mix, 0.1·min|λ_1| for 10%. EnergyDomain at r_D=7 is gap-exact (coherent ≡ B_bohr, dissipator ε ≈ 10⁻⁹) for any cell we measure. TimeDomain has a β-amplified coherent floor (2.4e-5 at n=6, β_phys=2) that caps low-T accuracy. Datapoint: n=6 β_phys=2 TRUE gap=3.97e-3; TimeDomain ε=10⁻⁶ off by +2.2%, TimeDomain ε=10⁻³ off by -14.5%. Detailed write-up at `drafts/error-analysis/epsilon-vs-gap-rule.md`. **Not yet promoted to a rule** — pending the P1/S1 thesis sweeps to establish actual gap magnitudes across the grid.

### construct_lindbladian, not build_dense_superoperator, for dense Lindbladian superops
- See [construct_lindbladian over build_dense_superop](feedback_construct_lindbladian_over_build_dense_superop.md). For a `Config{Lindbladian}`, the canonical dense d²×d² superop comes from `construct_lindbladian(jumps, cfg, ham)` in `src/furnace.jl` — BLAS-vectorised in-place per-jump assembly. `build_dense_superoperator` (`src/kms_geometry.jl`) is a generic d²-matvec fallback only for matvec closures with no direct constructor (channel walks, custom diagnostics).

### Per-term register triples on Config (qf-9z0, set 2026-05-04)
- See [Per-term registers qf-9z0](per_term_registers_qf_9z0.md). Config now carries three independent triples `(num_energy_bits_X, t0_X, w0_X)` for X ∈ {D, b_minus, b_plus} — one per QPE register the algorithm uses. Use helper accessors `register_t0_X / register_w0_X / register_r_X`. Legacy `t0/w0/num_energy_bits` kwargs still work via fallback. DLL has no b_minus/b_plus split (no outer/inner coherent integration in Ding–Li–Lin).

### δ ≥ 1e-3 and t0 ≥ δ floors for simulated-algorithm plots (set 2026-05-04)
- See [δ and t0 floors](project_delta_t0_floor.md). For all plot-worthy simulations: weak-measurement timestep `δ ≥ 10⁻³` (computational floor), and `t0 ≥ δ` since K·t0 below the δ/λ floor is wasted refinement. Total simulated error ≈ `K·t0 + δ/λ_gap` with K ≈ 2.4×10⁻³ at β=10, σ=1/β. **Cumulative-trajectory scaling is `δ/λ`, NOT δ²/λ** (per-step δ² × number-of-steps 1/(λδ)). λ is the Lindbladian gap, generically closing exponentially in n at low T — floor system-size-dependent. At our n≤6, β≤20 fixtures with λ~0.1: floor ~ 10⁻². The thesis ε=10⁻⁶ target is the **algorithm-level** Lindbladian-fixed-point precision (Krylov analysis), not the simulated-algorithm output. Lever for reducing quantum cost: t0 LARGER than δ via Fourier pairing — fewer qubits.

### Krylov-spectral trajectory predictors (qf-ev5, 8/9 children closed 2026-05-04)
- See [Krylov-spectral trajectory predictors](krylov_trajectory_qf_ev5.md). New `predict_lindbladian_trajectory` (Config{Lindbladian}) and `predict_channel_trajectory` (Config{Thermalize, :sweep}) in src/lindblad_action.jl. Same NamedTuple shape as `lindblad_action_integrate` so estimate_mixing_time consumes them unchanged. Single forward Arnoldi + dense eigen(H) lift gives biorthogonal decomposition WITHOUT a separate adjoint Krylov (qf-ev5.2 cancelled). At n=3 β=10: Lindbladian 32 matvecs / 3.8e-9 vs dense; channel 36 matvecs vs run_thermalize 30k steps / 1.6e-12 abs (byte-identical). `sweep_mixing_times` has new `method=:ode|:krylov` kwarg, 41x wall-time speedup at n=3 β=10. TrotterDomain + GQSP wired (qf-ev5.8 fix in commit 28736e4: trace-distance reference must be Gibbs in trotter eigenbasis, not Hamiltonian eigenbasis). 48 new test assertions across 9 testsets in test_predict_lindbladian.jl + test_predict_channel.jl. Only qf-ev5.9 (n=11 push) deferred.

### Krylov trajectory predictor benchmarks + accuracy (qf-ev5 follow-up, 2026-05-05)
- See [Krylov trajectory bench + accuracy](krylov_trajectory_bench_qf_ev5.md). Sandbox wall times (4-BLAS, krylovdim=30, Eb=12) up to n=7: predict_lindbladian (Energy) 0.028s→224s, predict_channel (Trotter+GQSP, d=1) 0.039s→15.3s. **Trotter+GQSP channel is ~15× faster than Lindbladian (Energy)** at n=7 because TimeDomain/Trotter `_accumulate_rho_jump!` reads precomputed NUFFT prefactors (no `exp()` in hot loop). **BohrDomain channel is d⁴·n** (d² Bohr keys × d³ per-key GEMM) — never use it for scaling work. **Krylovdim=30 is overkill**: kdim sweep at n=5,6, β∈{5,10,20} shows max trace-distance diff vs kdim=60 reference is 1e-9 (Lindbladian) / 1e-12 (Channel), saturating at the KrylovKit `tol` floor — 4-9 orders below the δ algorithmic floor (~1e-3 to 1e-5). Even kdim=15 already saturates. **predict_channel wall is δ-independent** (sweep δ ∈ {1e-2..1e-5} at n=5: 1.04× wall spread); smaller δ is strictly BETTER for accuracy at zero wall cost (CPTP |μ₁−1| 8.9e-8→8.7e-14, gap converges to true L gap). The δ ≥ 1e-3 floor from `project_delta_t0_floor.md` was a `run_thermalize` step-count constraint, no longer binding for the predictor. **n=11 path on 512 GB cluster**: TrotterDomain + GQSP, drop num_energy_bits to 10 (OFT cache: 4096·d²·16 bytes is the binding memory at 256GB for n=11/Eb=12), JULIA_NUM_THREADS=64 / BLAS=1 (channel is omega-loop-threaded above OMEGA_THREAD_THRESHOLD=50), krylovdim=20, expected ~30 min wall. Scripts: `scripts/scratch_predict_trajectory_bench.jl`, `scripts/scratch_krylov_dim_convergence.jl`, `scripts/scratch_predict_channel_delta_independence.jl`. Draft: `drafts/krylov-trajectory-predictor-scaling.md`.

### "Channel" always means TrotterDomain + GQSP (set 2026-05-05)
- See [Channel = Trotter+GQSP](feedback_channel_means_trotter_gqsp.md). When user says "channel" or "the implemented channel" they always mean `Config{Thermalize, TrotterDomain}` with `with_gqsp=true, gqsp_degree=1, jump_selection=:sweep` — that's what hardware would compile to. BohrDomain is d⁴·n (15×+ slower at n=7) and never the intended target for scaling work.

### Quadrature register recipe v2 (qf-yt9, canonical, 2026-05-12)
- See [Quadrature register recipe v2](quadrature_register_recipe_v2.md). **Replaces** the qf-7xt v1 recipe after β_phys/β_alg, per-register, typical-fixture, and `default_smooth_s` refactors. **Smooth Metro `default_smooth_s` recipe at ε=1e-9 EnergyDomain reference**: r_D=6 (r_D=7 at n=6,β_phys=1.0). Coherent slope-(-1) in r_+ unchanged from qf-7xt; ε=1e-9 is NOT reachable in TimeDomain (use EnergyDomain B which is `B_energy ≡ B_bohr` closed-form). **Methodology hard rule**: S1 dissipator is the ONLY sweep using BohrDomain (matvec scales d⁵, infeasible at n≥5); everything downstream uses EnergyDomain as the 1e-9 reference. Hybrid path for S1: dense `L_b` once per cell + matrix-free `apply_lindbladian!` for L_e per r, power iter on `(L_e−L_b)^†(L_e−L_b)` reading `‖A x_final‖`. Scripts: `scripts/scratch_quad_{S1_krylov,S34_coherent,synthesis}_campaign.jl`. Summary draft: `drafts/error-analysis/quadrature-convergence-summary-v2.md`.

### Quadrature register recipe v1 (qf-7xt, superseded 2026-05-12)
- See [Quadrature register recipe v1](quadrature_register_recipe_qf_7xt.md). Historical — measured at one fixture (n=4, β=10) before β_phys/β_alg split. Superseded by qf-yt9 v2 above. Kept for reference only; cite v2 going forward.

### qf-7xt n=5 kinky feasibility — dense + BohrDomain ref (2026-05-06)
- See [qf-7xt n=5 kinky feasibility](qf_7xt_n5_kinky_feasibility.md). At n=5 the kinky $r_D$ sweep (the only filter that needs $r_D \ge 10$) is feasible: dense + BohrDomain reference, ~17 min per $\beta$, ~50 min for $\beta \in \{5, 10, 20\}$. **Always BohrDomain ref**: r-independent (~2s), measures actual quadrature error; EnergyDomain ref is 25–100× slower and only catches NUFFT precision (Energy↔Time saturate at 2.5e-13). **Dense not Krylov**: opnorm itself is 0.4s on the 1024×1024 superop; bottleneck is `construct_lindbladian` for TimeDomain. Krylov svdsolve on $(L_t - L_b)$ wins at most 2× at $r=14$ — not worth complexity. Scripts: `scripts/scratch_quadrature_n5_{bench,matvec}.jl`. Krylov becomes the right tool only at $n \ge 6$ where dense memory ($d^4 \cdot 16$ bytes) starts hurting.

### B_time slope-(-1) origin: t=0 L'Hôpital sample, not η-jump (set 2026-05-05, scratches qf-2o5/qf-oiq)
- See [trap-rule t=0 L'Hôpital origin](trap_rule_t0_lhopital_origin.md). Convention: we always choose `η < t0' = 2T_+/2^r_+`, so the η-cutoff branch of `_compute_b_plus_metro` is dead code in the discretisation. The `‖B_bohr - B_time‖_op` slope-(-1) error is dominated by `t0'·b_+(0)·K^a(0)`, the L'Hôpital sample at the single t=0 grid point that the Cauchy P.V. integral excludes. Earlier "smooth-bump η-fix" line of work (qf-2o5, qf-oiq) is scrapped — it fixes a problem that doesn't exist when η < t0'. Beads `qf-xfa` filed for a thesis explanation paragraph. **Quantitative recipe lives in [Quadrature register recipe](quadrature_register_recipe_qf_7xt.md).**

### Shared-δt₀ TrottTrott for KMS coherent (qf-d0w, closed 2026-05-05)
- See [Shared-δt₀ scheme](qf_d0w_shared_delta_t0.md). For TrotterDomain KMS use `make_trotter_for_config(ham, config)` — picks shared-δt₀ scheme (one Trotterization at the elementary δt₀ shared across (D, b_-, b_+) registers, per-register `eigvals_t0_X = λ_S .^ M_X` are vector powers). Slope -2 in M_user recovered for B_trotter inner+outer loops; legacy single-cache saturated at ‖ΔB‖ ≈ 5e-5 regardless of M. Per-leg natural Trotter steps written generally: outer b_-(t/σ) → `register_t0_b_minus / σ`; inner b_+(τβ) → `β · register_t0_b_plus`. For σ = 1/β both coincide; written generally so σ ≠ 1/β plots also work. Legacy `TrottTrott(ham, t0, M)` constructor still works (per-register fields default to nothing → byte-identical pre-qf-d0w behaviour in B_trotter).

### Kinky Metro quadrature error: always dense, never Krylov-SVD (set 2026-05-06)
- See [Kinky dense not Krylov](feedback_kinky_dense_not_krylov.md). Kinky r_ref≥14 kills Krylov-SVD (16k labels/matvec × 120 iters). Gaussian/smooth: Krylov wins at n≥5 (r_ref=8 suffices). Kinky: always dense through n=6.

### Krylov two routes: Lindbladian (EnergyDomain) vs Channel (TrotterDomain+GQSP)
- See [Krylov two routes](krylov_two_routes.md). Krylov-Lindbladian = EnergyDomain with w0 for 1e-9 Bohr match; Krylov-Channel = TrotterDomain+GQSP faithful δ-step. Kinky Metro slow in Energy (needs r_D≥12). Coherent same in both. Use Lindbladian route for generator-level analysis, Channel route for algorithmic-step simulation.

### Non-Hermitian jumps in CKG: 3 separate constraints (set 2026-05-06)
- See [Non-Hermitian jumps physics & code](feedback_non_hermitian_jumps.md). (1) **Coherent term has an extra term** in the Metropolis-CKG construction for non-Hermitian jumps — see CKG paper. Existing $B = \alpha \cdot A^\dagger(\omega_1) A(\omega_2)$ formula does not carry over. (2) **Production jump sets must contain $(A, A^\dagger)$ pairs** for KMS DB to hold. Single non-Hermitian jumps are ONLY admissible in unit tests of internal code paths (e.g., the qf-in3 threading bit-match test) where two evaluations of the same physics are compared. (3) **OFT half-grid fold** — Hermitian jumps iterate only `w_raw <= 1e-12` and reuse the conjugate-transpose sandwich for the negative-frequency partner; non-Hermitian jumps iterate the full ω-grid. The qf-in3 threading dispatch and `_build_lindblad_work_list` in `src/krylov_matvec.jl` preserve this fold.

### gamma_norm_factor is grid-independent (qf-etx, closed 2026-05-06)
- See [gamma_norm grid-independent](gamma_norm_grid_independent_qf_etx.md). `_precompute_data` now uses closed-form `pick_gamma_sup(config) = 1.0` for every standard γ family, replacing the prior grid-dependent `1.0 / maximum(transition.(energy_labels))`. Cross-grid Lindbladian comparisons no longer need the `/gnf` workaround; BohrDomain `construct_lindbladian` is byte-identical across `(r_D, w0_D)` choices. Smooth Metro Lindbladians and τ_mix shift up by ~5% post-fix; TINF-02 BSON references regenerated. 110 new regression assertions in `test/test_gamma_norm_invariance.jl`. The qf-mto FAIR.4 Λ_max table and the `no_normalization_for_metropolis_taumix_plots` rule still hold qualitatively but with ~5% shifted numbers (qf-etx.9 deferred re-sweep).

### Test coverage parity + 1e-9 cross-domain invariant (set 2026-05-12)
- See [Test coverage and 1e-9 invariant](feedback_test_coverage_and_1e9_invariant.md). Three rules when touching the test suite: (1) every NO_SANDBOX test needs a sandbox-runnable toned-down equivalent — prefer shrinking n/r over loosening tolerances; (2) if no sandbox equivalent is possible, ALWAYS print an explicit ⚠️ NO_SANDBOX TEST REQUIRED warning when changing the covered code, telling the user to run `QUANTUMFURNACE_FULL_TESTS=true`; (3) cross-domain tests must be controllable to 1e-9 (Energy↔Bohr, Time↔{Energy,Bohr}, Trotter↔Time, faithful-channel↔ideal-Lindblad). 1e-4 "good enough for a plot" is NOT a passing test — historically masked register-sizing bugs. Encoded in `.claude/rules/julia-code.md` Test Suite section.

## Completed Work

### Error analysis epic qf-b4d (completed + revised 2026-05-05; quadrature numbers SUPERSEDED by qf-7xt 2026-05-06)
- See [Error analysis qf-b4d](error_analysis_qf_b4d.md). Master synthesis (Trotter $M$, generator splitting, cross-check) is still authoritative. **The quadrature register recommendations within this note are superseded** by [Quadrature register recipe qf-7xt](quadrature_register_recipe_qf_7xt.md) — the prior "joint $r_{b\pm} = r$" numbers conflated $r_-$ and $r_+$, and the prior dissipative recommendations relied on a confounded methodology (now redone with the unified $\omega$-range / truncation principle).
- **Still authoritative findings**: L_diss Strang slope $-2$ in $M_D$; generator splitting slope $+2$ in $\delta$ (1.985–1.997); coh-diss splitting ~30× smaller than jump-wise at large $\beta$; Trotter $\beta^3$ prefactor empirically $\beta^{\sim 2.6}$; **single shared trotter.t0 caps B-Trotter at $\sim 5 \times 10^{-5}$** (qf-d0w split caches required; pending P2 task).
- **Cross-check** at $(n=3, \beta=10, \varepsilon_{\text{gen}}=10^{-4}, \text{smooth})$: asymptotic TD = $1.6 \times 10^{-3} = 16 \times \varepsilon_{\text{gen}}$ (CURRENT code with single trotter cache).
- 8 scripts in `scripts/scratch_{coherent_quadrature,coherent_quadrature_v2,coherent_quadrature_split,dissipative_quadrature,trotter_M_selection,trotter_M_selection_v2,generator_splitting,kinky_slope_check,b4d5_crosscheck}.jl`. The first three coherent-quadrature scratches and the original dissipative-quadrature scratch are HISTORICAL — the qf-7xt convergence scripts are the canonical reference.
- Open code changes: qf-d0w (split Trotter caches). qf-oiq (smooth η-regularization) was scrapped — the slope-(-1) does not come from the η-jump (see [trap-rule t=0 L'Hôpital origin](trap_rule_t0_lhopital_origin.md)).

### Smooth-Metropolis s-decision (qf-3il, closed 2026-05-04)
- See [Smooth-Metropolis s-decision qf-3il](smooth_metropolis_s_decision_qf_3il.md). s = 0.25 confirmed as evidence-backed default (already the `sweep_mixing_times` default). Trade: 4-bit energy-register savings at ε=1e-4 (16× fewer N_ω points; uniform across β=5,10,20) for 3-6% τ_mix penalty (uniform across n=3,4,5, β=5,10,20). Two commits 0be215b..b768bfd. Fills the [REF] todo after Eq. eq:smooth-metro in 2_methods.tex line 305. Quadrature analysis is n-independent → extends to all n.

### Jump-selection :sweep | :random (qf-2vo, closed 2026-05-04)
- See [Jump-selection qf-2vo](jump_selection_qf_2vo.md). 5 atomic commits 4858427..c723d4b. Default switched to :sweep (thesis-preferred Lie-Trotter); :random remains opt-in. 13 new tests; full suite 5029/5029.

### Multi-rank DLL epic (qf-7go, closed 2026-05-03)
- See [Multi-rank DLL completion](multirank_dll_complete_qf_7go.md). 7/7 sub-issues closed. ~561 new tests across 4 new test files.
- Construction: `dll_multichannel_translates(base; centers, weights)` builds DLLMultiChannelFilter via ShiftedSymmetricFilter wrappers. Symmetrised translates `q_ℓ(ν) = √(w/2)·[q_base(ν−ν_ℓ) + q_base(ν+ν_ℓ)]`.
- **KMS-DBC of multi-channel G**: `G^multi = Σ_ℓ G^(ℓ)` by linearity of Theorem 10's canonical-form map. `verify_detailed_balance` rel_norm ≤ 1e-10 for k ∈ {1,2,4} at β ∈ {1,5,10}.
- **Headline (n=3, target_ε=1e-3)**: DLL Gaussian β=20 collapse (τ_mix=264 at k=1) **fully recovers to 13.3 at k=8 — 19.9× speedup**. DLL Metropolis 5-13× speedup k=1→k=8 across β.
- **Revises §8.7 prediction**: H2 holds, but the optimal diagonal-of-α shape is recoverable inside DLL via shifted base translates — not unique to CKG.
- **Cost caveat**: per-channel ‖G^(ℓ)‖ inherits Remark 23 → ‖G^multi‖ = O(k·β^1.4). β=20 Gaussian net total-cost reduction ~2.5×, not 19.9×.
- §9 thesis draft at `drafts/dll-multirank-taumix-findings.md`. Plot at `drafts/figures/numerics/dll_multirank_taumix.{png,pdf}`.
- Library refactor: `oft_nufft_at_zero` → `oft_nufft_at_zero_list::Vector{Matrix}` in TimeDomain DLL precompute (length-1 single, length-k multi).

### Fair comparison qf-mto epic (closed 2026-05-02)
- See [Fair-comparison qf-mto epic](fair_comparison_dirichlet_qf_mto.md). 7/7 sub-issues closed. 6 atomic commits 3c11403..00ffebe (after src integration 78b36cf, 3cdd462).
- New library: `src/kms_geometry.jl` (10 public functions: kms_inner_product/norm/variance, kms_dirichlet_form, build_dense_superoperator, spectral_gap_kms, max_dirichlet_rate_kms, intrinsic_mixing_ratio, dissipator_one_to_one_norm_bound, dissipator_trace_alpha, hs_operator_norm). Reuses `materialize_discriminant`. 60 new tests, krylov_spectral_gap cross-check 2.5e-16.
- **Headline ρ_intrinsic = λ/Λ_max**: CKG sM is 1.02-1.43× ahead of DLL Metro across (n,β); not a structural tie, modest CKG advantage growing with β. DLL Gauss collapses 6× at (n=3,β=20).
- **H2 confirmed**: DLL Metro α_diag within 6% L² of CKG sM across all 12 cells. Rank decoupled — Metropolis-shaped diagonal is the driver. Predicts qf-7go.6 multi-rank DLL gives null result.
- **Bound check**: τ_pred ≈ 3·τ_meas constant factor across 32 cells (slope ~1, R²=0.93). λ is the right primary predictor.
- Surprise: CKG sM λ non-monotone in β at n ∈ {4,5} between β=5 and β=10 — not yet ruled out as fixture artefact, worth follow-up at n ∈ {6,7}.
- Thesis §8 in `drafts/ckg-vs-dll-comparison-findings.md` (203→369 lines) + 4 new figures in `drafts/figures/numerics/`.

### CKG vs DLL first comparison findings (2026-05-02)
- See [CKG vs DLL first findings](ckg_vs_dll_first_findings.md). Draft for thesis comparison subsection at `drafts/ckg-vs-dll-comparison-findings.md`.
- Headline: **CKG smooth-Metropolis ≈ DLL Metropolis** within 10–25% across n ∈ {3,4,5}, β ∈ {1..20}, with CKG slightly faster in most cells. **DLL Gaussian collapses at low T** — 8× slower than CKG at n=3, β=20.
- Mechanism: rank(α) = number of independent Lindblad channels per coupling. Standard DLL is rank-1 (single q per coupling); CKG is generically full-rank (~24 channels at n=3). DLL Gaussian's q(ν) is super-exp localised in βν → near-dephasing single channel at low T. DLL Metropolis keeps fat exp tail in βν.
- DLL paper Remark 23: wider Kossakowski → β-factor extra cost in coherent-term G block-encoding. DLL total cost: Õ(β² · S · t_mix) vs CKG: Õ(β · t_mix). So even *if* DLL beat CKG in t_mix, it would have to beat by factor β to be competitive in total cost.
- **Numerical confirmation of Remark 23**: log–log slope of ‖G‖ vs β over β ∈ {5..30} is **1.35–1.46 for DLL Metropolis** at n ∈ {3,4,5}; CKG smooth-Metro has flatter ~0.7–0.8 slopes. DLL Gauss shrinks at n=3 (slope −0.33, same collapse mechanism as τ_mix). Script: `scripts/scratch_coherent_norm_lcu_check.jl`. Absolute ‖G‖ < 1 on rescaled fixtures, but the *scaling* is the asymptotic statement and matches the paper.
- The qf-lkb closure note "DLL Gaussian τ_mix 40–60% lower than CKG/DLL Metro" is true at β ≤ 5 only — inverts at β ≥ 15.

### CKG EnergyDomain sweep epic qf-lkb.11 (closed 2026-05-02)
- See [CKG EnergyDomain sweep](ckg_energy_domain_sweep_qf_lkb_11.md).
- 5 atomic commits: dispatch relaxation, BSON loader dual-schema, scaling
  benchmark, qf-lkb.6 plot regen.
- Bohr→Energy matvec speedup: 22.5× @ n=4, scaling continues favorably.
- Smooth-Metropolis defaults locked at a=0, s=0.25 (thesis convention) in
  `sweep_mixing_times` and qf-lkb.6 plot.
- 12 new tests in `test_lindblad_action.jl` (Bohr ≈ Energy ~6e-10 rel) +
  2 alloc tests in `test_krylov_matvec.jl`.

### Lindbladian-action ODE integrator epic qf-lkb (closed 2026-05-02)
- See [Lindbladian integrator epic qf-lkb](lindbladian_integrator_epic_qf_lkb.md).
- 9/9 children + epic closed; 8 commits 8542685..710c24d (plus prerequisite f160085).
- New library API: `lindblad_action_integrate`, `discriminant_action_integrate`, `integrate_to_gibbs`, `sweep_mixing_times`; `estimate_mixing_time` gained vector-method overload.
- Matrix-free DLL `apply_lindbladian!` added (qf-lkb.9) — closes the n>5 cliff for DLL.
- Comparison plot at `drafts/figures/numerics/ckg_vs_dll_taumix.{png,pdf}`; headline: DLL Gaussian τ_mix is 40-60% lower than CKG/DLL Metropolis.
- ~485 new tests across `test_lindblad_action.jl` + `test_mixing.jl`; full DLL/CKG regression clean.

### Cleanup epic qf-fzj (GQSP + Hamiltonian + DLL, completed 2026-05-02)
- See [Cleanup epic qf-fzj](cleanup_qf_fzj.md).
- 17 atomic commits across 8 sub-issues; 17 files changed, +396 / -584 LOC (-188 net).
- Audits in `.planning/phases/52-cleanup-gqsp-ham-dll/audit-{gqsp,hamiltonian,dll}.md`;
  physics-check + verifier reports in same directory.
- Two real bugs fixed: (1) `with_gqsp + DLL` was validatable but unrunnable — now rejects;
  (2) `HamHam(...; periodic=false)` multi-term ctor was silently ignoring kwarg — now forwards.
- Hermitisation symmetry: new `_coherent_unitary_step` helper in `src/coherent.jl`
  hermitises B before both GQSP and matrix-exp branches (Thermalize path used to hit raw B
  on GQSP, matched Trajectory path's pre-existing explicit `hermitianize!`).
- New shared test fixtures in `test/test_helpers.jl`: `make_dll_n3_system(β)`,
  `assert_kms_skew_symmetric(α, ν_grid, β)`. Use these for any future Kossakowski / DLL test.
- Wall time +49s (4m42s → 5m31s) — **deliberate trade** for n=2→n=3 fixture migration in
  `test_dll_dissipator.jl` per `feedback_n3_minimum_test_size.md`.

### DLL Metropolis filter (qf-wmg, completed 2026-05-01)
- See [DLL Metropolis filter qf-wmg](dll_metropolis_filter_qf_wmg.md).
- 10 atomic commits (qf-wmg.1..10), +257 new tests, full DLL suite 1648/1648 passing.
- `DLLMetropolisFilter{T}(beta; S=2)` implements Eq. 3.19–3.20 alongside `DLLGaussianFilter`.
- BohrDomain + TimeDomain end-to-end; coherent G via 2D type-3 NUFFT over [-S, S]² (Nν=256).
- Refactor: extracted `_dll_coherent_from_g_tt` helper (filter-agnostic Steps 3+4).
- Headline qualitative win at β=10, ν=-0.45: |α_metro| = 0.95 vs |α_gauss| = 0.06 (15× larger).
- Diagnostic: `scripts/plot_kossakowski_metropolis_vs_gaussian.jl` → 2x3 PNG heatmap.

### DLL performance refactor (qf-hur, completed 2026-05-01)
- See [DLL performance refactor qf-hur](dll_performance_refactor_qf_hur.md).
- 3 phases, 5 atomic commits. Closed-form `dll_coherent_op_bohr` (Phase A); NUFFT-based `dll_lindblad_op_time` at ω=0 (Phase B); closed-form g(t,t')+2D NUFFT for `dll_coherent_op_time` (Phase C).
- DLL TimeDomain end-to-end: 167 s → 0.26 s (640×); now **3.5× faster than CKG TimeDomain** at n=5, β=5.
- DLL BohrDomain end-to-end: 18 s → 0.013 s (1380×).
- 3142 → 3171 tests (+29 new cross-checks at FINUFFT precision).

### Phase 51: DLL Lindbladian — Bohr + Time domains (qf-3i8, completed 2026-05-01)
- See [DLL implementation Phase 51](dll_implementation_phase51.md).
- Five atomic phases (DLL-1..5), 19 commits, 1169 new tests (1973 → 3142 total).
- `src/filters.jl`, `src/dll.jl` modules; `Config.filter::Union{Nothing, AbstractFilter}` opt-in.
- TrotterDomain deferred per user scope; rejected at `validate_config!`.
- Paper typo found: Eq. 3.7 third equality has `A^a(t)A^a(t')` reversed; correct order is `A^a(t')A^a(t)` (paired-with-t' LEFT).
- KMS-DBC skew-symmetry α(ν,ν') = α(-ν',-ν)·e^{-β(ν+ν')/2} (Eq. 4.7) verified at machine precision for both CKG and DLL.

### Phase 48 prerequisite: Hamiltonian families and 2D builders (qf-k1u.5, completed 2026-04-30)
- See [Phase 48 Hamiltonian families and 2D builders](hamiltonian_families_phase48.md).
- Three cached families: `heis_xxx_zzdisordered_periodic_n*.bson` (1D full disorder),
  `heis_xxx_clean_periodic_n*.bson` (1D + ε-disorder), `heis_xxz_2d_*.bson` (2D XXZ
  J_z=1.5, periodic both directions, X+Z ε-disorder).
- New 2D builders in `src/hamiltonian.jl`: `_pad_two_site_op`,
  `_construct_2d_heisenberg_base`, `find_ideal_2d_heisenberg`. Refactored
  shared optimisation kernel `_optimize_disordered_heisenberg`.
- 1973 tests (1917 existing + 56 new in `test/test_hamiltonian.jl`).
- Physics verified: F2 1D AFM e_GS/bond approaches Bethe-ansatz limit −1.7726;
  F3 2x2 (bipartite-compatible PBC) gives `<S_stag²>/n² = 0.73` confirming Néel
  correlations.

### GQSP coherent step in simulator (epic qf-63j, completed 2026-04-30)
- See [GQSP simulator integration](gqsp_simulator_integration.md) for full notes.
- New `Config` fields `with_gqsp::Bool=false`, `gqsp_degree::Int=1` opt into the
  Jacobi-Anger Chebyshev polynomial `f_d(B/α)` for the coherent step (Time/Trotter).
- 1860 → 1917 tests; α-helper, slope-(d+1), regression vs matrix-exp baseline.

### Phase 43: Bi-exponential fitting (completed 2026-03-04)
- `BiexpFitResult` struct and `fit_biexponential_decay` in `src/fitting.jl`
- `estimate_mixing_time` extended with `model=:biexp` keyword in `src/mixing.jl`
- Extrapolation via `Roots.Bisection` (no closed-form for multi-exponential)
- Result: <0.001% error on synthetic data (vs 0.13% single-exp, vs 26% on real data)
- 1273 tests pass (1246 existing + 27 new)

## Key Architecture

- Fitting: `src/fitting.jl` (LsqFit.jl v0.15.x, Levenberg-Marquardt)
- Mixing estimation: `src/mixing.jl` (post-processing wrapper)
- Module: `src/QuantumFurnace.jl` — already has `using Roots`, `using LsqFit`
- Structs: `FitResult`, `BiexpFitResult`, `MixingTimeEstimate` (all immutable)
- Tests: `test/test_fitting.jl`, `test/test_mixing.jl`
- Diagnostic scripts: `scripts/diagnose_floor_v3.jl`, `scripts/mixing_time_extrapolate_verify.jl`

## Findings from Diagnostics (2026-03-03)
- Floor scales linearly with delta (exponent ~1.0 with coherent unitary): floor ≈ 0.068 × delta
- Coherent unitary is essential — without it, floor saturates at ~0.0017
- Single-exp overestimates offset C by ~47% → 26% error in extrapolated mixing time
- Root cause: sensitivity of `t = -ln((ε-C)/A)/gap` when ε ≈ C

## TrotterDomain Floor Analysis (2026-03-04)

### Two-component floor model: `floor = k_energy × δ + floor_Trotter(constant)`
- EnergyDomain k (scales with δ): n=3→0.059, n=4→0.094, n=5→0.142
- TrotterDomain constant floor: n=3→~8e-6, n=4→~8e-5 (10× larger!), n=5→~3.5e-5
- **n=4 "anomaly" explained**: at δ=0.0001, Trotter constant dominates → effective k≈0.879
- Trotter error in OFT/R-matrix is NOT the cause (~1e-8 for all n)
- Files: `src/errors.jl`, `scripts/trotter_error_comparison.jl`, `scripts/floor_superoperator.jl`
- Required δ to reach ε=1e-4: n=3→0.001, n=4→0.0001, n=5→0.0005

### Bohr Frequency Collision Root Cause (2026-03-04)
- **Even-n periodic chains are bipartite** → sublattice operator U=diag(1,-1,1,-1,...) exists
- U anti-commutes with hopping (XX+YY), commutes with on-site Z disorder
- In single-magnon sectors (Sz=±(n/2-1)): ZZ diagonal is constant → U creates spectral symmetry
- Result: E_k(Sz=+m) + E_{d+1-k}(Sz=-m) = C (exact cross-sector pairing, spread ~5e-16)
- This gives exact Bohr frequency collisions (gap reversal palindrome)
- n=4: 12 collisions (all from Sz=±1, single-magnon), n=6: 30 (all from Sz=±2)
- Odd n: ring is frustrated (not bipartite) → no U → no collisions
- Multi-magnon sectors: ZZ diagonal varies → symmetry breaks (no collisions from Sz=0)
- Scripts: `scripts/bohr_collision_sectors.jl`, `scripts/zz_disorder_symmetry_break.jl`

### ZZ Bond Disorder Symmetry Breaking Experiment (2026-03-04)
- Adding ZZ bond disorder (Σ_q ε_q Z_q Z_{q+1}) on top of existing Z disorder breaks bipartite pairing
- **Collisions eliminated**: n=4: 12→0, n=6: 30→0. Cross-sector pair sum spread 5e-16 → 1.4e-03
- ZZ-only disorder (without Z) fails: P=∏X_i commutes with all bilinear terms → Sz=±m degenerate → nu_min=0
- Need BOTH Z (breaks P) + ZZ (breaks bipartite pairing)
- **Trotter mismatch resolved (2026-03-04)**: HamHam now supports multiple disordering terms natively

## Bi-exp Verification Results

### Z-only disordering (2026-03-04, legacy)
- n=3 (δ=0.001): 1.35% error, PASS
- n=4 (δ=0.0001): 0.78% error, PASS (fit reports converged=false, still works)
- n=5 (δ=0.0005): 0.29% error, PASS
- Scripts: `scripts/biexp_mixing_verify.jl` (EnergyDomain), `scripts/biexp_mixing_verify_trotter.jl` (TrotterDomain)

### Z+ZZ disordering (2026-03-04, symmetry-breaking)
- All n values give <2% error when extrap time is sufficient
- δ=0.0005 for all n; floor ~3-6e-5, well below target 1e-4
- n=3: 0.22%, n=4: 1.3-1.6%, n=5: 0.12% (EnergyDomain)
- **Key finding**: bi-exp fit needs data ~3× past crossing to separate decay from floor
- Previous n=4 "failure" was simply insufficient extrapolation time (T=60 vs crossing at t≈59)
- Collisions broken for all n (min Bohr gap > 0)
- No qualitative even/odd difference — n=4 anomaly eliminated
- extrap_time_map: n=3→100, n=4→150, n=5→150
- Script: `scripts/biexp_mixing_verify_zzdisordering.jl`

## Hamiltonian Construction Notes
- `find_ideal_heisenberg`: optimizes disorder for max min-gap, `disordering_terms` kwarg (default `[[Z]]`)
- `HamHam` struct: `disordering_terms::Union{Vector{Vector{Matrix{Complex{T}}}}, Nothing}` — multiple terms
- `disordering_coeffs::Union{Vector{Vector{T}}, Nothing}` — per-term per-site coefficients
- Constructors: (1) no disorder, (2) multi-term, (2b) single-term convenience, (3) from NamedTuple
- NamedTuple constructor handles both legacy `disordering_term` and new `disordering_terms` keys
- `_construct_disordering_terms` has overloads for single-term and multi-term
- `_trotterize2` iterates over all disordering terms (each gets its own Trotter layer)
- `pad_term` handles 2-site terms with periodic wrapping
- Trotter: even-n uses 2 bond groups (clean), odd-n needs 3rd group for wrapping bond

## Thesis
- [Thesis Structure and Status](thesis_structure.md) — MSc thesis structure, completed/empty sections, original contributions (Props 5,7,9,10)
- [Thesis Numerics Plan](thesis_numerics_plan.md) — Full plan for Ch5 numerical data: plots, cross-checks of all analytical results, priority order
- [GQSP angle-finding decision](gqsp_angle_finding_decision.md) — Final: use GQSP (Motlagh-Wiebe) + Berntson-Sünderhauf FFT for Q, not optimization
- [Trotter sign convention mismatch](thesis_sign_convention_trotter.md) — 1_preliminaries.tex mixes $\ee^{+\ii Ht}$ (QPE) with $\ee^{-\ii Ht}$ (Trotter subsection); 2_methods.tex uses $+$. Unify in a later pass
- [TODO: bring R_b reflection into 2_methods.tex GQSP](todo_methods_gqsp_reflection.md) — prelim uses W = R_b U_B explicitly; methods chapter still treats U_B alone as the walk. Also resolve sin/cos convention.
- [GQSP thesis discrepancies (verified by POC)](gqsp_thesis_discrepancies.md) — three concrete fixes for circuits in 1_preliminaries.tex (Fig. circ:gqsp) and 2_methods.tex (Alg. alg:coh): L_0 placement reversed, slot-reassignment claim wrong, Alg. 1 needs asymmetric prep. Full draft in drafts/gqsp-circuit-thesis-discrepancies.md.
- [GQSP figure follows BS (transpose) convention, not MW-native](gqsp_bs_vs_mw_convention.md) — L_0 placed LAST in circ:gqsp / alg:coh because POC angles are BS-form; MW Eq. 52 / Fig. 2 has L_0 first (transpose-related). Don't "fix" by flipping.
- [Thesis Review chapter — quantum vs classical Gibbs sampling](thesis_review_chapter.md) — beads epic `qf-yk4`; framing: where could quantum Gibbs sampling first beat classical? Regime map + /lit-review corpus.

## Conventions
- [β_phys vs β_alg distinction (qf-6vr)](beta_phys_beta_alg_convention.md) — Hamiltonian fixtures store a rescaled spectrum, so every β has two meanings: `β_phys` (against `H_phys`) and `β_alg = β_phys · ham.rescaling_factor` (against the stored spectrum). `cfg.beta` stays β_alg; `cfg.beta_phys` is the new explicit kwarg. Sweep harness takes `beta_phys_values`; sidecar emits the triple `:beta_phys`, `:beta_alg`, `:rescaling_factor`. **Canonical β_phys grid: `{0.25, 0.5, 1.0}`** (2026-05-11; replaces legacy `{5, 10, 20}` β_alg). 0.25 = smallest β with non-trivial thermal contrast (S/log d ≈ 0.80 uniformly across n); 1.0 = practical ceiling (σ tightens to 0.01 at n=10; `default_smooth_s` would force s≈25 at n=11 — kept at fixed s=0.25 pending audit of large-s smooth-Metro behaviour).
- [gap_phys vs gap_alg relation (qf-8fr)](gap_phys_alg_relation_qf_8fr.md) — every gap / τ_mix reported by `krylov_spectral_gap`, `predict_lindbladian_trajectory.spectral_gap`, etc. is in **algorithm time units**. **Physical gap = R · alg gap** (R = `ham.rescaling_factor`), i.e., `λ_gap_phys = R · λ_gap_alg`, `τ_mix_phys = τ_mix_alg / R`. INDEPENDENT of the β_phys→β_alg input correction: that one fixes the Gibbs state inside L; this one is just unit-relabelling at readout. For extensive 1D H, R ∝ n, so the raw code output looks like 1/n decay even when the physical gap is Ω(1). **Cross-validated KB16 / Bardet23 / BC25 / RSA24 Ω(1)-gap theorems in qf-8fr**: slopes of log λ_gap_phys vs log n in [-0.2, +0.2] for classical 1D Ising at β_phys=0.5 (slope -0.05) and disordered Heisenberg at β_phys ∈ {0.25, 0.5, 1.0} (slopes -0.11, -0.20, -0.15), n=4..8, after dropping n=3 PBC outlier. Use alg-gap for gate-complexity statements (algorithm runs ~1/λ_gap_alg steps); phys-gap for theorem comparisons. Full draft: `drafts/phys-vs-alg-conventions.md`, evidence: `drafts/classical-ising-crosscheck.md`.
- [No normalisation for CKG vs DLL Metro τ_mix plots](no_normalization_for_metropolis_taumix_plots.md) — Λ_max coincidentally agrees to ~1%, production-normalised plots are already at the fair-comparison scale. Long-form: `drafts/no-normalization-metropolis-taumix.md`.
- [PREP gates inside block encodings](feedback_prep_convention.md) — figures draw `U_A`/`U_{B_a}`/`U_{diss}` as monolithic boxes (PREP inside); externalisation is gate-cost optimisation mentioned only in §Cost analysis text and a one-liner in `circ:gqsp` caption.
- [Circuits folder stays .tex-only](feedback_circuits_folder.md) — never leave .aux/.log/.pdf/.DS_Store in `drafts/circuits/`
- [Thesis .tex are read-only; drafts are readable Markdown](feedback_thesis_tex_readonly.md) — never Edit/Write `.tex` in `supplementary-informations/`; drafts go to `drafts/<slug>.md` as Markdown with inline `$math$`, not LaTeX dumps. Canonical example: `drafts/quantum-circuits-basics.md`
- [Plot PNGs may fail Read with API 400](feedback_image_read_fails.md) — verify generated plots via stdout diagnostics or ask the user; do not retry the same PNG after a "Could not process image" error
- [DLL implementation: consult the paper before coding](feedback_dll_paper_first.md) — for epic `qf-3i8`, always read the DLL paper formula before writing each piece; do not infer from CKG analogues
- [Random-vs-sweep is dissipative-only](feedback_jump_selection_dissipative_only.md) — coherent/GQSP fixes are tracked separately; do not conflate scopes when changing per-step structure
- [Tests should use n=3 minimum, not n=2](feedback_n3_minimum_test_size.md) — 2-qubit toys hide bugs that surface at n=3; default to `heis_disordered_periodic_n3.bson`
- [Targeted tests during integration](feedback_targeted_tests_during_integration.md) — between integration commits, run only relevant test files (not the full ~3-min suite); run full suite only at end of phase
- [Always include β=10 in numerical tests](feedback_beta_test_values.md) — β=1 hides errors; β=5 borderline; β=10 reliably exposes them. Sweep β ∈ {1, 5, 10} for any β-dependent verification
- [Kossakowski matrices must be checked for KMS skew-symmetry](feedback_kossakowski_skew_symmetry_check.md) — α(ν,ν') = α(-ν',-ν) e^{-β(ν+ν')/2}; structural witness of detailed balance at matrix level (Ding–Li–Lin 2024 Eq. 4.7)
- [Fast/NUFFT path is the default; legacy slow paths are references only](feedback_fast_path_default.md) — when extending DLL operators to new filters, reuse the existing 2D-NUFFT factorisation (`dll_coherent_op_time` Step 3+4); never propose the legacy `O(Nt²·n³)` Riemann sum as a default path
- [`test_helpers.jl` silently shadows script-level constants](feedback_test_helpers_const_shadow.md) — `include("test/test_helpers.jl")` rebinds `SIGMA, BETA, NUM_QUBITS, ...` without warning. Don't reuse those names; use prefixes like `CKG_SIGMA`. Sanity-check by printing the constant before use.

## Plotting
- [Thesis colour palette](reference_thesis_colors.md) — named colours (pinegreen, bordeaux, dustyplum, deepplum, aubergine, slateblue, sage, ochre, terracotta, dustyteal, mustard); refer to by name in plotting requests
- [Thesis gradient palettes](reference_gradient_palettes.md) — cold (mint→navy), warm (cream→mulberry), diverging (teal↔purple) gradients for heatmaps and continuous-scale plots; prefer over `:inferno`/`:viridis`

## Memory Persistence
- Memory lives in `.claude-memory/` at repo root, symlinked to `~/.claude/projects/.../memory/`
- Setup via `SessionStart` hook in `.claude/settings.json` (runs once per session)
- Migrated from `.planning/memory/` + `PreToolUse` hook on 2026-03-05
