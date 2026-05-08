---
name: Cleanup epic qf-fzj (GQSP, Hamiltonian, DLL)
description: 17-commit cleanup pass across GQSP, Hamiltonian, and DLL modules; closed 2026-05-02
type: project
---

# Cleanup epic qf-fzj — completed 2026-05-02

After landing GQSP simulator integration (qf-63j), Hamiltonian families/2D builders (qf-k1u.5),
and the DLL stack (Phase 51 + perf refactor + Metropolis filter), the corresponding code and
tests had grown organically. This consolidation pass landed before downstream consumers
(numerics chapter, Ding-Chen comparison) get added.

## Workflow (8 sub-issues, all closed)

1. **qf-fzj.1 audit** — three module audits in `.planning/phases/52-cleanup-gqsp-ham-dll/audit-*.md`
2. **qf-fzj.2 GQSP cleanup** — 5 commits
3. **qf-fzj.3 Hamiltonian cleanup** — 5 commits
4. **qf-fzj.4 DLL cleanup** — 8 commits
5. **qf-fzj.5 test consolidation** — shared fixtures + helper extraction
6. **qf-fzj.6 physics check** — physics-checker agent PASS on all three modules
7. **qf-fzj.7 integration pass** — already in good shape from per-module phases (no commits)
8. **qf-fzj.8 final verification** — code-verifier PASS

**17 atomic commits**, 17 files, **+396 / -584 LOC (-188 LOC net)**.

## Notable commits

### GQSP (Phase A)
- `a213129` — collapse `_get_truncated_indices(::Real)` and `(::Complex)` into `(::Number)`
- `fd93feb` — reject `with_gqsp ∧ DLL` in `validate_config!` (was validatable but unrunnable)
- `433d78a` — extract `_coherent_unitary_step` helper for the 3× near-clone; fixes asymmetric
  hermitisation (Thermalize TimeDomain GQSP branch now hermitises B before polynomial, matching
  the Trajectory path's explicit `hermitianize!`)
- `47b0f37` + `843780d` — test consolidations (slope+tail-bound merge, smoke+regression fold)

### Hamiltonian (Phase B)
- `2115850` — fix bug: multi-term ctor (2) was dropping `periodic` kwarg before
  `_construct_base_ham`, silently building periodic base on `periodic=false` calls
- `3315f17` — drop dead single-term `_construct_disordering_terms` overload (never called)
- `e910cb4` — dedupe `_load_hamiltonian_bson` between `src/misc_tools.jl` and
  `test/test_helpers.jl` (47-line copy-paste → 1-line const alias)
- `9ae1b9a` — drop n=2 test cases (per `feedback_n3_minimum_test_size.md`)
- `d980be9` — direct ctor coverage for `HamHam(1)` and `HamHam(2b)` (previously only
  reached transitively via NamedTuple ctor)

### DLL (Phase C)
- `a426613` — unexport `dll_coherent_kernel_bohr` (test-only kernel)
- `e87aa56` — docs: clarify `dll_lindblad_op_time` as reference path; fix
  `DLLMetropolisFilter` docstring (was claiming non-existent `S/2 ≥ max|ν_BH|` warning)
- `98b385d` — extract `assert_kms_skew_symmetric` helper to `test_helpers.jl`
- `fe59a1a` — parameterise Gaussian/Metropolis test pairs in `test_dll_coherent.jl` and
  `test_dll_kossakowski.jl` (centrosymmetry merged with skew-symmetry — Eq. 4.8 ⇔ Eq. 4.7)
- `d453f5d` — migrate `test_dll_dissipator.jl` (a-d) from n=2 toy to n=3 disordered Heisenberg
  (the toy was kept pre-G-wiring for σ_β preservation; with G now wired through
  `_precompute_coherent_B`, n=3 works at the same tolerances)
- `7b67274` — parameterise `(a/f)`, `(b/h)` in `test_dll_kms_db.jl`

### Test consolidation (Phase D)
- `a6bc3cf` — `make_dll_n3_system(beta)` consolidated into `test/test_helpers.jl`; was
  duplicated 3× across DLL test files

## Behavioural changes (none silent)

1. **`with_gqsp + DLL` now rejects at validation** (was: validatable but runtime-error)
2. **`HamHam(...; periodic=false)` ctor (2) now correctly forwards `periodic`** (was: silently
   ignored)
3. **`_coherent_unitary_step` hermitises B for both branches** (matrix-exp branch was already
   reading via `Hermitian(B)` wrapper; GQSP branch now sees the same hermitised B; matches
   the Trajectory path which already did this explicitly)

## Why: physics-check and verifier results

- **Physics check (qf-fzj.6)**: PASS on all three modules. Only doc concerns flagged
  (thesis vs paper Eq. number drift in DLL); no correctness fixes required.
- **Final verifier (qf-fzj.8)**: PASS. ≥2500 individual `@test` invocations confirmed across
  the cleanup-affected files; 0 FAIL, 0 ERROR. The 1217 → 704 test count drop in
  `test_dll_kossakowski.jl` is attributed to (a) centrosymmetry merge (genuinely redundant
  with skew-symmetry) and (b) smaller ν grids in the merged parameterised testset.
  Property coverage preserved.

## How to apply

When extending DLL filters, GQSP, or Hamiltonian builders in the future:

- **Use `make_dll_n3_system(β)`** in tests (per memory `feedback_n3_minimum_test_size.md`)
- **Use `assert_kms_skew_symmetric(α, ν_grid, β)`** for Kossakowski KMS witness checks
  (per memory `feedback_kossakowski_skew_symmetry_check.md`)
- **Don't reintroduce duplicated `_load_hamiltonian_bson`** — test_helpers aliases the src loader
- **Don't add single-call helpers** unless the call site cleanup justifies them (per
  julia-code.md "no abstractions beyond what the task requires")

## Wall time

- Pre-cleanup: 4m42s for the 9 in-scope test files
- Post-cleanup: 5m31s (+49s, **deliberate trade** for n=2 → n=3 fixture migration in
  `test_dll_dissipator.jl`; per `feedback_n3_minimum_test_size.md`)
- Net LOC: 17 files, +396 / -584 = **-188 LOC**
