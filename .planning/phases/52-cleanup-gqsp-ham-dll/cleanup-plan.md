# Cleanup epic qf-fzj — execution plan

Synthesised from audits at `audit-{gqsp,hamiltonian,dll}.md` (qf-fzj.1).
Each numbered item below becomes one atomic commit.

---

## Phase A — qf-fzj.2 (GQSP)

### A1. Drop dead `_get_truncated_indices(::Real)` overload
- File: `src/coherent.jl:385-393`
- Two byte-identical overloads; only the `Complex` one is exercised.
- Action: drop the `Real` overload, or replace both with one `AbstractVector{<:Number}` signature.
- Tests: existing GQSP polynomial tests + truncated-func tests cover the live path.

### A2. Tighten `validate_config!` to reject `with_gqsp ∧ DLL`
- File: `src/misc_tools.jl:170-183`
- DLL precompute returns `(filter, time_labels, t0, oft_nufft_at_zero)` — no `b_minus`/`b_plus`/`gamma_norm_factor`.
- `_precompute_coherent_unitary` would crash. Either reject upfront or extend the implementation; cleanup-only path is rejection.
- Tests: add `test_gqsp_config.jl` case asserting `validate_config!` rejects `(DLL, TimeDomain, with_gqsp=true)`.
- Note: this is a **public-facing API change** (a config that was previously "valid in validation but unrunnable" now errors at validation). Acceptable — better to fail fast.

### A3. Extract `_coherent_unitary_step` helper for the 3× near-clone
- Files: `src/coherent.jl:108-127`, `src/trajectories.jl:101-110`.
- New helper signature (in `coherent.jl`):
  ```julia
  _coherent_unitary_step(B, b_minus, b_plus, t0_sim, gamma_norm_factor,
                         delta_eff, with_gqsp::Bool, gqsp_degree::Int) -> Matrix
  ```
- Make Thermalize-path GQSP branches go through `Hermitian(B)` to match the matrix-exp branch and the Trajectory path's explicit `hermitianize!`.
- Three call sites become a one-liner each.
- Tests: existing regression vs matrix-exp + slope-(d+1) cover both branches.

### A4. Cache `b_minus_l1`, `b_plus_l1` next to truncated dicts (optional)
- Tiny perf win, no behaviour change. Skip if it adds surface area.

### A5. Test consolidation in `test_gqsp_polynomial.jl` (slope-(d+1) + tail-bound)
- Soft duplicates at lines 126 and 166. Merge into one parameterised testset over `axis ∈ (:delta, :alpha)`.

### A6. Test consolidation in `test_gqsp_thermalize.jl` (smoke + regression)
- Lines 44/66 and 55/77 are weak duplicates. Fold smoke into regression.

---

## Phase B — qf-fzj.3 (Hamiltonian)

### B1. Fix ctor (2) `periodic` kwarg drop (BUG)
- File: `src/hamiltonian.jl:127` — `_construct_base_ham(terms, coeffs, num_qubits)` should be `_construct_base_ham(terms, coeffs, num_qubits; periodic=periodic)`.
- Tests: add `test_hamiltonian.jl` case checking `HamHam(...; periodic=false)` builds an open-chain base.

### B2. Remove dead single-term `_construct_disordering_terms` overload
- File: `src/hamiltonian.jl:528-541` — never called.
- If a single-term call site ever returns, callers wrap as `[term]`/`[coeffs]`.

### B3. Collapse ctor (1) into thin wrapper around ctor (2)
- File: `src/hamiltonian.jl:60-107` → 5-10 line wrapper that calls (2) with empty `disordering_terms`/`disordering_coeffs`.
- Saves ~30 LOC of boilerplate. Keep ctor (2b) as already-thin wrapper.
- BSON wire format: ctor (1) currently writes `disordering_terms = nothing`; the wrapper would write `Vector{Vector{Matrix}}[]`. Loader at `_unpack_disordering_fields` already accepts both. Verify cached BSONs still load.

### B4. Dedupe `_load_hamiltonian_bson` between src and test
- Files: `src/misc_tools.jl:27-74` and `test/test_helpers.jl:21-67` — ~50 lines duplicated.
- Make `test/test_helpers.jl` call `QuantumFurnace._load_hamiltonian_bson` directly.

### B5. Fix `load_hamiltonian` filename pattern OR add a F1/F2/F3 loader
- File: `src/misc_tools.jl:13` — currently builds `heis_xxx_disordered_periodic_n*.bson` which doesn't match any new family.
- Three families need different patterns: `heis_xxx_zzdisordered_periodic_n*`, `heis_xxx_clean_periodic_n*`, `heis_xxz_2d_*`.
- Action: parametrise filename construction on a "family" symbol or string suffix.
- Out of scope if the user prefers minimal changes — flag for follow-up.

### B6. Drop n=2 test cases violating n=3 minimum
- `test_hamiltonian.jl:19, 53, 64, 117` — n=2 placements + (2,2) lattices.
- Re-base on n=3 / (3,1) / (1,3) / (3,3). Per `feedback_n3_minimum_test_size.md`.

### B7. Add direct ctor coverage for ctors (1), (2), (2b)
- Currently only ctor (3) (NamedTuple) is exercised. Add minimal regression for each.

---

## Phase C — qf-fzj.4 (DLL)

### C1. Demote `dll_coherent_kernel_bohr` from exports
- File: `src/QuantumFurnace.jl:92`
- Used only by test-only `dll_coherent_op_time_legacy`. Either rename to `_dll_coherent_kernel_bohr` (private) or drop from export list.
- Verify no `drafts/` / `scripts/` references the public name.

### C2. Document `dll_lindblad_op_time` as reference path
- File: `src/dll.jl:73` (docstring) — add "production amortises this via FINUFFT in `_precompute_data`."
- No code change; doc-only.

### C3. Fix `DLLMetropolisFilter` docstring vs `validate_config!` drift
- File: `src/filters.jl:293-298` claims `validate_config!` warns on `S/2 < max|ν_BH|`. `src/misc_tools.jl:194-203` only checks `S > 0`.
- Action: update docstring to "caller must ensure S/2 ≥ max|ν_BH|; otherwise asymptote silently breaks". (Adding the actual check requires Hamiltonian access in `_precompute_data`; defer.)

### C4. Migrate `test_dll_dissipator.jl` (a-d) from n=2 toy to n=3
- File: `test/test_dll_dissipator.jl:5-12` (toy fixture).
- Per `feedback_n3_minimum_test_size.md`. The toy was needed pre-G-wiring; now redundant.
- Verify n=3 fixture passes the same tolerances.

### C5. Consolidate Gaussian/Metropolis test pairs into parameterised testsets
- Files: `test_dll_coherent.jl`, `test_dll_kms_db.jl`, `test_dll_kossakowski.jl`.
- Pattern: `for filter in (DLLGaussianFilter(β), DLLMetropolisFilter(β))` over the property-checking testsets.
- Estimated win: ~40% LOC reduction in those three files.
- Per audit pairings:
  - `coherent (a)/(i)` — Hermiticity Bohr
  - `coherent (b)/(j)` — Hermiticity Time
  - `coherent (c)/(k)` — Bohr ↔ Time on G
  - `coherent (d)/(l)` — keep separate (different intent: Gaussian asserts ≤1, Metropolis asserts O(1))
  - `kms_db (b)/(h)` — Time KMS-DB up to quadrature
  - `kossakowski (a-c)/(h2-h5)` — shape/Hermitian/PSD/skew-sym

### C6. Consolidate Bohr ↔ Time agreement testsets
- Four near-identical tests at `coherent.c`, `coherent.k`, `dissipator.c`, `kms_db.i`.
- Single parameterised testset `@testset "Bohr ↔ Time agreement (DLL)"` over (filter, fixture, β, what, tol) matrix.

### C7. Shared KMS skew-symmetry helper
- New `test/test_helpers.jl::assert_kms_skew_symmetric(α, ν_grid; β, atol)` — shared by `test_dll_kossakowski.jl` Gaussian + Metropolis variants and any CKG counterpart.
- Per `feedback_kossakowski_skew_symmetry_check.md`.

### C8. Inline `_time_oft_prefactor_dll`
- File: `src/filters.jl:171` — single call site at `:182` in `filter_time_cutoff(::DLLGaussianFilter)`. Inline.

### C9. Rename `_dll_coherent_op_time_legacy` test-only doc note
- File: `src/dll.jl:542-550` — add "TEST-ONLY reference; never on a production path." per `feedback_fast_path_default.md`.

### C10. Restrict legacy reference filter type to DLL filters
- File: `src/dll.jl:546` — `dll_coherent_op_time_legacy(... filter::AbstractFilter ...)` could silently accept CKG `GaussianFilter` and produce nonsense.
- Either change to `Union{DLLGaussianFilter, DLLMetropolisFilter}` or assert `filter isa AbstractDLLFilter` (need to introduce a marker if pursuing the latter).

### C11. Make `dll_coherent_op_bohr` `beta` arg defensive or remove
- File: `src/dll.jl:171-174` (docstring), `src/coherent.jl:59` (call site).
- Either remove `beta` from signature (read from filter), or assert `beta == filter.beta` defensively.

---

## Phase D — qf-fzj.5 (test-suite consolidation)

### D1. Add `test/fixtures.jl` with shared n=3 fixture loader
- Single canonical `make_n3_disordered_heis(β)` returning a HamHam.
- Replaces per-test re-load.

### D2. Move shared `assert_kms_skew_symmetric` into `test/test_helpers.jl`
- Used by C7 above.

### D3. Record full-suite wall time before/after
- Pre-cleanup baseline: see `baseline_run.jl` output.
- Post-cleanup: re-run after C7.

---

## Phase E — qf-fzj.6 (physics check)

Run `physics-checker` agent on three papers:

### E1. Motlagh–Wiebe GQSP
- Polynomial choice (Jacobi-Anger Chebyshev), α-helper, slope-(d+1).
- BS-form vs MW-native L_0 (memory: deliberate; DO NOT flip).
- R_b reflection (memory: TODO in 2_methods.tex).

### E2. Berntson–Sünderhauf FFT angle solver
- Sign / convention. L_0 placement matches BS.

### E3. Ding–Li–Lin DLL
- Eq. 3.7 paired-with-t' LEFT order (paper typo: third equality).
- Eq. 3.19–3.20 Metropolis filter.
- Eq. 4.7 Kossakowski KMS skew-symmetry.
- Filter-vs-coherent factorisation `g(t,t')` shared kernel.

Deliverable: `physics-check-report.md` with discrepancies + severity.

---

## Phase F — qf-fzj.7 (integration pass)

`code-integrator` agent: ensure cross-module consistency post-cleanup.

### F1. Naming sweep
- `with_X` for Config opt-in flags. Already consistent (`with_gqsp`, `with_coherent`).
- `_xxx_bohr` / `_xxx_time` suffix for domain variants. Already consistent for DLL.

### F2. Export-list audit
- Drop `dll_coherent_kernel_bohr` (per C1).
- Confirm only public-facing names exported.

### F3. Config field grouping
- `with_gqsp` / `gqsp_degree` (GQSP block) → already adjacent at `structs.jl:122-123`.
- `filter` (DLL block) → at `:126`. Confirm comment block separates the two.

---

## Phase G — qf-fzj.8 (final verification)

### G1. `code-verifier` agent reviews cumulative diff vs epic-start commit
### G2. Full `Pkg.test()` (or per-file fallback if sandbox OOMs)
### G3. Wall-time delta record
### G4. MEMORY.md "Completed Work" entry referencing the epic

---

## Sequencing rules
- Phase A, B, C are independent — could run in parallel, but per `feedback_targeted_tests_during_integration.md`, run only the relevant `test_*.jl` between commits.
- Phase D depends on C (KMS fixture).
- Phase E (physics-check) depends on A, B, C completing — papers vs final cleaned code.
- Phase F (integration) depends on D and E.
- Phase G (verification) is terminal.

## Atomic commit count estimate
- A: 6 commits
- B: 7 commits
- C: 11 commits
- D: 3 commits
- E: physics-check report (1 file, no commits unless E flags fixes)
- F: 1-2 commits (integration adjustments)
- G: 1-2 commits (memory note + any verifier-flagged fixes)

**Total: ~30 atomic commits.**
