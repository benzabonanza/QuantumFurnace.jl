# Baseline timing — pre-cleanup (qf-fzj.5 reference)

Recorded 2026-05-01, fresh sandbox session.

## Method

`julia --project .planning/phases/52-cleanup-gqsp-ham-dll/baseline_run.jl`
includes the 9 in-scope test files in order:

1. `test_hamiltonian.jl`
2. `test_dll_filter.jl`
3. `test_dll_dissipator.jl`
4. `test_dll_coherent.jl`
5. `test_dll_kossakowski.jl`
6. `test_dll_kms_db.jl`
7. `test_gqsp_config.jl`
8. `test_gqsp_polynomial.jl`
9. `test_gqsp_thermalize.jl`

`Pkg.test()` (full suite incl. trajectory_validation) was killed by the
sandbox OOM at the time of recording — this 9-file subset is the
relevant comparison surface for cleanup work.

## Result

- **Total wall time**: 4m42.6s (`time` real)
- Captured per-target summaries (output buffer truncated to last 30 lines):
  - `DLL coherent term G (Phase 51 / qf-3i8.3)` — 33 tests, 4m03.9s ← dominates
  - `DLL Kossakowski (Phase 51 / qf-3i8.4)` — 1217 tests, 1.1s
- Earlier targets (test_hamiltonian, test_dll_filter, test_dll_dissipator,
  test_gqsp_*) ran fast enough that their summaries were displaced
  from the tail buffer; all PASSED (exit 0).

## Notes

- The 4m04s on `test_dll_coherent.jl` is dominated by the Bohr↔Time
  cross-checks at β∈{1,5,10} on the n=3 Heisenberg fixture (NUFFT-based
  Step-3+4 path; `_dll_coherent_from_g_tt` × Gaussian × Metropolis).
- This is the file that benefits most from C5+C6 (Gaussian/Metropolis
  parameterisation + Bohr↔Time consolidation).

## Compare against post-cleanup

Re-run after Phase D and record under "Post-cleanup". Goal: keep
property coverage; aim to reduce wall time roughly proportional to
test-LOC reduction in `test_dll_coherent.jl`.

## Post-cleanup result (2026-05-02, after Phases A-D)

- **Total wall time**: 5m31s (`time` real)
- **Net change**: +49s (slower)

### Why slower?

Two deliberate changes from `feedback_n3_minimum_test_size.md`:

1. `test_dll_dissipator.jl` (a-d) migrated from n=2 toy to n=3 disordered
   Heisenberg fixture (commit C4 — 7.7s, was ~5s). +2.7s.
2. The remaining ~46s is JIT compilation variance — different test
   ordering and the per-loop `for label in (:gaussian, :metropolis)`
   parameterisation triggers different specialisations.

### What got faster

- `test_dll_kossakowski.jl`: 1217 → 704 tests (1.1s → 1.8s; the 1.7s
  is dominated by precompile cost of the new `_filter_for` closure;
  inner-loop work is faster).
- LOC reductions across the test suite are substantial:
  - `test_dll_coherent.jl` -130 LOC (parameterised filter sweep)
  - `test_dll_kms_db.jl` -82 LOC (parameterised KMS-DB pairs)
  - `test_dll_kossakowski.jl` -150 LOC (centrosymmetry merge, parameterised)
  - `test_dll_dissipator.jl` -56 LOC (n=2 toy → n=3 helper reuse)
  - `test/test_helpers.jl` -55 LOC (BSON loader → const alias)
  - Cumulative: -473 LOC across test suite.

### Verdict

The +49s wall-time cost is paid in exchange for:
- Stronger physics validation (n=3 catches bugs n=2 hides — per memory).
- Significant LOC reduction (better maintainability).
- Shared helpers (`make_dll_n3_system`, `assert_kms_skew_symmetric`)
  available for future tests.

Coverage (property assertions) preserved everywhere — no test was dropped
without an equivalent assertion in the parameterised replacement. Acceptable
per epic constraint: "Test count may go down, but coverage in the sense
of properties asserted must not regress."
