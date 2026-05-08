# QuantumFurnace.jl src/ Code Review (qf-6z9.4)

Date: 2026-05-06. Reviewer: code-verifier agent. Scope: every file under
`src/` (excluding `src/staging/` and `src/python/`).

## Overall summary

The code is in **good shape overall**. Recent epics (qf-in3 threading, qf-9z0
per-register triples, qf-d0w shared-δt₀ Trotter, qf-etx grid-independent gamma
norm, qf-7go multi-rank DLL, qf-7xt quadrature recipes) are cleanly wired and
the per-domain dispatch tables are coherent. The ω-loop threaded variants
mirror their serial counterparts faithfully. No matrix-free path silently falls
back into a dense build, no exported symbol is undefined, no `include()`d file
is missing, no obvious off-by-one in indices, and the test suite passes on
sub-targets I could spot-run (`test_register_validation`,
`test_register_independence`, `test_trotter_caches`).

That said, I found a handful of small dead-code, perf, and correctness items
worth cleaning up. Counts by severity:

- **bug**: 2 (one latent, one cosmetic)
- **dead-code**: 5
- **perf**: 4
- **test-gap**: 3
- **nit**: 3

No BLOCKERs. Most findings are individually small atomic-commit fixes.

---

## Findings

### Bugs

- **[bug]** `src/trotter_domain.jl:225` — `for (term, coupling) in (groups.noncommuting[1], groups.noncommuting[2])` is a malformed iteration. The right-hand side is a 2-tuple of vectors, not a paired iterator; iterating yields `groups.noncommuting[1]` (a `Vector{Vector{Matrix{CT}}}`) on the first pass and `groups.noncommuting[2]` (a `Vector{T}`) on the second pass — destructuring those as `(term, coupling)` is undefined / incorrect. **Latent**: the entire branch sits behind `if length(groups.noncommuting[1]) != 0`, and `_does_term_differ_at_both_sites` (line 244) only ever flags non-commuting terms when there are 2-site terms that differ at exactly one site. The default XYZ-Heisenberg + Z-disordering Hamiltonians used everywhere have all 2-site terms in the "kinda commuting" bucket (XX, YY, ZZ all differ at both sites), so the branch never fires in tests. But any future Hamiltonian that introduces a single asymmetric 2-site term (say XY or XZ) would walk into this corrupted code path. **Fix**: change to `for (term, coupling) in zip(groups.noncommuting[1], groups.noncommuting[2])`. Also retype `sequence_2site_not_commuting = []` (line 223) to `Matrix{ComplexF64}[]` for type stability. Add a regression test that constructs an asymmetric 2-site Hamiltonian and verifies `_trotterize2` matches `exp(-iHt)` to the expected Trotter error.

- **[bug]** `src/qi_tools.jl:162` — `frobenius_norm(A::Matrix{<:Complex})` returns `sqrt(sum(abs.(eigvals(A)).^2))`, which equals `‖A‖_F` only for Hermitian / normal `A`. For a general non-Hermitian matrix the Frobenius norm is `sqrt(sum(σ_i^2))` (singular values), not `sqrt(sum(|λ_i|^2))`. **Cosmetic** because the function is exported but never called inside `src/` or `test/`; still, an exported function with that name silently giving wrong values for non-normal inputs is a footgun. **Fix**: replace with `LinearAlgebra.norm(A)` (which is the Frobenius norm for matrices) or `sqrt(real(tr(A' * A)))`, or — since it's unused — drop the export and the function.

### Dead code

- **[dead-code]** `src/krylov_eigsolve.jl:309` — `A_nu2_dag = zeros(T, dim, dim)` is allocated and `fill!(A_nu2_dag, 0)`-style logic was clearly intended (cf. the comment on line 308: "Allocate A_nu2_dag buffer (one per call, acceptable for Bohr)"), but the function (`_accumulate_jump_sandwich!`, BohrDomain, lines 296–332) never reads or writes the buffer. The actual scatter pattern for `rho * A_nu2_dag` runs straight into `sc.sandwich_tmp` (lines 317–325), bypassing `A_nu2_dag` entirely. Net effect: a `dim×dim` complex matrix is allocated per `apply_delta_channel!` call for the BohrDomain path, then discarded. **Fix**: delete the line. Saves 16·d² bytes per matvec.

- **[dead-code]** `src/energy_domain.jl:210` — `return energy_labels[start_index:end_index]` is unreachable; the function returned on line 208. **Fix**: delete line 210.

- **[dead-code]** `src/structs.jl:222` (`JumpOp.orthogonal`) — the `orthogonal` field is only consulted by the deprecated `time_oft!` and `trotter_oft!` (`src/ofts.jl:50,103`), which are themselves only kept "for tests" (errors.jl + test_dm_scaling.jl). All mainline matvec / channel / trajectory code keys off `JumpOp.hermitian` instead. The field is essentially redundant with `hermitian` for real-valued Pauli jumps (X, Y, Z) and unused everywhere else. Either retire the field (and the deprecated `*_oft!` functions / `OFTCaches`) or document explicitly that it's only for legacy `time_oft!`/`trotter_oft!`.

- **[dead-code]** `src/ofts.jl:21–133` (`time_oft!`, `trotter_oft!`) — explicitly labeled "Depricated but used for tests". They duplicate logic that lives in `_prepare_oft_nufft_prefactors` (NUFFT path) and could be retired. Their existence keeps `OFTCaches` and `JumpOp.orthogonal` alive, both of which are otherwise dead. If the legacy tests stay, mark the functions with the standard `@deprecate` or at minimum spell the misspelling correctly ("Deprecated").

- **[dead-code]** `src/qi_tools.jl:162` — `frobenius_norm` (see [bug] above). Unused; export should be dropped.

### Performance

- **[perf]** `src/krylov_matvec.jl:306` and `src/krylov_matvec.jl:356` — `A_nu2_dag = zeros(T, dim, dim)` allocates a `dim×dim` complex matrix on every Krylov matvec for the BohrDomain Lindbladian (`apply_lindbladian!` and `apply_adjoint_lindbladian!`). Unlike the dead allocation in `_accumulate_jump_sandwich!` above, here the buffer **is** used (lines 312–317, 319, etc.). With KrylovKit doing ~30 matvecs per restart × ~5 restarts on the typical sweep, that's ~150 unnecessary `dim²` allocations per spectral-gap call. **Fix**: add an `A_nu2_dag::Matrix{Complex{T}}` field to `KrylovScratch` and reuse it (alongside the existing `sandwich_tmp` and `sandwich_out` slots). The BohrDomain workspace constructor is the only one that needs to populate it — the Energy/Time/Trotter paths don't use a dense `A_nu2_dag`.

- **[perf]** `src/jump_workers.jl:814,898,977` and `src/trajectories.jl:443,526,603` — the threaded `_accumulate_rho_jump_*` and `_precompute_R_threaded_*` helpers all do `task_scratches = [ThermalizeScratch(CT, dim) for _ in 1:length(chunks)]` per call. For `_accumulate_rho_jump_*` this sits inside `run_thermalize`'s outer δ-step loop — at `n=5, β=10` that's `n_steps × n_jumps ≈ 5000 × 5 = 25 000` ThermalizeScratch allocations per simulation, each carrying 7 dim² matrices (≈ 7 × 1024 × 16 bytes ≈ 112 kB). Net allocation churn ~3 GB. **Fix**: mirror the qf-in3.4 pattern in `KrylovScratch` — pre-allocate a `task_scratches::Vector{ThermalizeScratch{CT}}` pool on the `Workspace` (or pass it in via the precomputed_data NamedTuple) and reuse it. The Krylov path already does this pre-allocation correctly (`KrylovScratch(CT, dim; num_threads=…)` in `src/structs.jl:444`). The Thermalize path is the one that still pays per-call.

- **[perf]** `src/jump_workers.jl:803,888` (and the parallel `_precompute_R` paths) — `half_indices = [i for i in eachindex(energy_labels) if energy_labels[i] <= 1e-12]` allocates a fresh Vector{Int} per call. Could be cached on the workspace alongside `work_list`. Minor (~few hundred ints) but it's per-step.

- **[perf]** `src/lindblad_action.jl:383,386–393` (`_arnoldi_factorize`) — modified Gram–Schmidt does `w = f(Q[:, j])` (allocates) then `H[i, j] = dot(view(Q, :, i), w); w .-= H[i, j] .* view(Q, :, i)` (more allocations from views). This is called once per `_krylov_spectral_decomposition` call so the impact is modest, but `mul!` and pre-allocated Gram-Schmidt scratch would clean it up. Document as "minor follow-up; current Arnoldi runs in <1 s for n ≤ 6".

### Test gaps

- **[test-gap]** `src/trotter_domain.jl:225` non-commuting branch — see the `[bug]` above. No test currently exercises a Hamiltonian whose `group_hamiltonian_terms` returns a non-empty `noncommuting` bucket. A targeted unit test (e.g. a 2-site Hamiltonian with both XX and XZ terms) would catch the malformed iteration.

- **[test-gap]** Multi-channel DLL filter validation messages — `validate_config!` (`src/misc_tools.jl:336–359`) emits per-channel error strings for `DLLMultiChannelFilter`. I didn't find a regression test that asserts the specific error string contents (e.g. "channel ℓ.beta=… does not match Config.beta=…"). `test_dll_multichannel_filter.jl` covers construction; the `validate_config!` error path may be uncovered.

- **[test-gap]** GQSP polynomial branch `d ≥ 2` (`src/coherent.jl:405–438`) — the `d=1` fast-path is exercised by `test_gqsp_polynomial.jl` and `test_gqsp_thermalize.jl`. The Clenshaw recurrence for `d ≥ 2` is implemented but I didn't see a numerical test that bumps `gqsp_degree=2` or `=3` and compares against the Bessel/Jacobi-Anger expansion or the `d=1` reference. Worth a single assertion that says "at `d=1`, `f₁(B/α)` agrees with the matrix-free implementation for a small fixed `B`; at `d=2`, the residual to `exp(-iδB)` shrinks by `O((δα)²)` over `d=1`".

### Nits

- **[nit]** `src/coherent.jl:21` — `Depricated` typo in the source comment of `src/ofts.jl:21`. (Already flagged under dead-code.)

- **[nit]** `src/jump_workers.jl:571` (`_partition_range`) and `src/trajectories.jl:725` (`_partition_trajectories`) implement the same partitioning algorithm twice. `_partition_trajectories` could just call `_partition_range`; the docstring already says they share the algorithm.

- **[nit]** `src/qi_tools.jl:167–183` and `src/qi_tools.jl:185–201` — two `is_density_matrix` methods that differ only in tolerance (`digits=15` vs `digits=13`). The non-`Hermitian{Complex{T}}` overload is the catch-all; the more specific overload keeps a tighter tolerance. Worth a docstring comment that explains why the looser tolerance is needed in the specific-T branch (presumably for `Float32`?).

- **[nit]** `src/structs.jl:551` — `scratch` field of `Workspace` is untyped (`Any`). The hot-path consumers cast it via `ws.scratch::KrylovScratch{T}` etc., which is fine for type stability at the call site, but the struct itself reports `scratch::Any`. This is a deliberate trade (one struct, four scratch shapes, parametric on simulation type already). Document the rationale next to the field, or — bolder — make `Workspace` parametric on a fifth `scratch` type.

---

## Aside: things I checked and they're fine

- All `export`ed symbols in `src/QuantumFurnace.jl` are defined somewhere in
  the module. All `include()`d files exist.
- The qf-in3 threading bit-match invariants in `test_threading.jl` cover both
  Hermitian and non-Hermitian jump paths and are passing (the recently fixed
  stale-test bugs at lines 89, 188 are unrelated).
- The qf-9z0 per-register triple plumbing is consistent across
  `furnace_utensils.jl::_precompute_data` (all five domain-specific overloads),
  `coherent.jl::_precompute_coherent_unitary`, `coherent.jl::B_time` and
  `coherent.jl::B_trotter` (with the explicit `t0_outer / t0_inner` pair),
  `simulation_time.jl::compute_simulation_time` (b_minus / b_plus separately),
  and `make_trotter_for_config`. Legacy single-register fallback is preserved
  symmetrically.
- The qf-d0w shared-δt₀ Trotter constructor (`src/trotter_domain.jl:99–139`)
  enforces the integer-M condition with `tol = 1e-9 · max(1, M[i])` and
  reports the offending leg in the error message — clean.
- The qf-etx `gamma_norm_factor = 1.0 / pick_gamma_sup(config) ≡ 1.0` change
  is wired correctly: every `_precompute_data` overload that sets
  `gamma_norm_factor` uses the new closed form, and the DLL overloads
  (`furnace_utensils.jl:115–175`) correctly do **not** set it (DLL has no γ
  rates). `pick_gamma_sup` is only defined for KMS / GNS, never queried in
  DLL paths.
- The qf-7go multi-channel DLL operator-level overloads in
  `src/dll_multichannel.jl` correctly return `Vector{Matrix}` for
  `dll_lindblad_op_*` (so the dissipator path can sum per-channel
  `L^(ℓ) ρ (L^(ℓ))†` without cross terms) and a single `Matrix` for
  `dll_coherent_op_*` (where `G^multi = Σ_ℓ G^(ℓ)` is linear, no cross
  terms). The simulator coupling at `_accumulate_dll_bohr_dissipator!`
  (`dll_multichannel.jl:508`) iterates over the per-channel operators
  correctly.
- `pick_transition` 2-arg form (`src/energy_domain.jl:7–44`) is allocation-
  free on the hot path, by intent and by inspection. Smooth-Metropolis
  branches handle `s == 0` and `a == 0` consistently with the closure
  `_pick_transition_kms` (line 47).
- `apply_delta_channel!` correctly handles aliasing between `rho` and
  `rho_eff` (using `sc.sandwich_out` as the staging buffer when
  `U_coherent !== nothing`, and the input `rho` directly otherwise; the
  `K0 * rho_eff * K0'` and `U_res * rho_eff * U_res'` GEMMs both write to
  `sc.rho_out` via `sc.sandwich_tmp`, so no aliasing surprise).
- `Workspace{KrylovSpectrum, BohrDomain, DLL}` is a separate constructor in
  `src/krylov_workspace.jl:295–349` that doesn't call `_accumulate_R_total!`
  (which has no DLL dispatch); instead it uses `_accumulate_R_total_dll!`
  which builds the per-jump `L_a = dll_lindblad_op_bohr(...)` matrices
  in-place and accumulates `R = Σ_a L_a' L_a` correctly. The matrix-free
  matvec at `apply_lindbladian!`/`apply_adjoint_lindbladian!` for
  `Config{Lindbladian, BohrDomain, DLL}` (krylov_matvec.jl:399, 440) reuses
  the existing sandwich helpers.

## Recommended commit order

If picking up only some of these, I'd order them:

1. `[bug]` trotter_domain.jl noncommuting iteration (one-liner + test)
2. `[dead-code]` krylov_eigsolve.jl:309 unused A_nu2_dag (one-liner)
3. `[dead-code]` energy_domain.jl:210 unreachable return (one-liner)
4. `[perf]` move BohrDomain `A_nu2_dag` into `KrylovScratch` (medium; touches
   structs.jl + krylov_workspace.jl + krylov_matvec.jl)
5. `[perf]` move ThermalizeScratch task pool out of the hot loops (medium;
   touches structs.jl + jump_workers.jl + trajectories.jl)
6. `[bug]` `frobenius_norm` (one-liner + drop export)
7. `[dead-code]` retire deprecated `time_oft!`/`trotter_oft!` (medium;
   touches ofts.jl + structs.jl OFTCaches + JumpOp.orthogonal + a couple
   of regression tests in `test_dm_scaling.jl`, `test_dm_detailed_balance.jl`)

Items 1–6 are independent and should each fit under a single atomic commit
of <30 lines.
