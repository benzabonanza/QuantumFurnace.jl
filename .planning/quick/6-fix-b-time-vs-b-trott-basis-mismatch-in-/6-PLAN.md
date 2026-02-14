---
phase: quick-6
plan: 01
type: execute
wave: 1
depends_on: []
files_modified:
  - src/coherent.jl
  - test/test_dm_scaling.jl
autonomous: true
must_haves:
  truths:
    - "B_trotter() produces results in a consistent basis (Trotter eigenbasis) by transforming jump.in_eigenbasis to Trotter eigenbasis before mixing with Trotter time evolution operators"
    - "DMTST-05 coherent term consistency test passes with B_bohr ~ B_trott to within 1e-5"
    - "All existing tests continue to pass (no regressions)"
  artifacts:
    - path: "src/coherent.jl"
      provides: "Fixed B_trotter() single-jump and multi-jump functions with basis-correct jump operators"
      contains: "trafo_from_eigen_to_trotter"
    - path: "test/test_dm_scaling.jl"
      provides: "Tightened DMTST-05 threshold from 0.02 to 1e-5"
      contains: "dist_bohr_trott < 1e-5"
  key_links:
    - from: "src/coherent.jl B_trotter()"
      to: "trotter.trafo_from_eigen_to_trotter"
      via: "basis transformation of jump.in_eigenbasis before matrix multiply"
      pattern: "trafo_from_eigen_to_trotter.*jump\\.in_eigenbasis"
---

<objective>
Fix basis mismatch bug in B_trotter() functions in src/coherent.jl.

Purpose: B_trotter() mixes jump.in_eigenbasis (Hamiltonian eigenbasis) with Trotter time evolution
operators (Trotter eigenbasis), producing incorrect results. The fix transforms the jump operator
to Trotter eigenbasis before the computation, matching the pattern used in the DMTST-06 OFT fix
(quick task 3). After fixing, tighten the DMTST-05 test threshold from 0.02 to 1e-5.

Output: Corrected B_trotter() functions and tightened test threshold.
</objective>

<context>
@src/coherent.jl
@test/test_dm_scaling.jl
@src/trotter_domain.jl
</context>

<tasks>

<task type="auto">
  <name>Task 1: Fix basis mismatch in both B_trotter() functions and tighten DMTST-05 threshold</name>
  <files>src/coherent.jl, test/test_dm_scaling.jl</files>
  <action>
In src/coherent.jl, fix both B_trotter() functions (single-jump at ~line 214 and multi-jump at ~line 243).

The bug: `jump.in_eigenbasis` is in the Hamiltonian eigenbasis, but the Trotter time evolution
`Diagonal(trotter.eigvals_t0 .^ n)` is diagonal in the Trotter eigenbasis. These are different
bases, so multiplying them together is a basis mismatch.

For the SINGLE-JUMP B_trotter() (~line 214):
- At the start of the function (after `trotter_time_evolution` definition), add:
  ```julia
  U = trotter.trafo_from_eigen_to_trotter
  jump_in_trotter = U * jump.in_eigenbasis * U'
  ```
- Replace `jump.in_eigenbasis` with `jump_in_trotter` (2 occurrences in the b_plus_summand loop at line 230):
  - `jump.in_eigenbasis'` becomes `jump_in_trotter'`
  - `jump.in_eigenbasis` becomes `jump_in_trotter`

For the MULTI-JUMP B_trotter() (~line 243):
- At the start of the function (after `trotter_time_evolution` definition), add:
  ```julia
  U = trotter.trafo_from_eigen_to_trotter
  ```
- Inside the inner loop over `jump_a in jumps`, transform each jump before use:
  ```julia
  jump_a_trotter = U * jump_a.in_eigenbasis * U'
  ```
  Then use `jump_a_trotter'` and `jump_a_trotter` instead of `jump_a.in_eigenbasis'` and `jump_a.in_eigenbasis`.

Do NOT modify B_time() functions -- they are correct (H-eigenbasis time evolution with H-eigenbasis jumps).

Do NOT modify coherent_term_trotter() (~line 303) -- it has the same bug pattern but is currently
unused (all call sites are commented out). Fixing dead code risks introducing regressions. If it
becomes active later, it should be fixed at that time.

In test/test_dm_scaling.jl, in the DMTST-05 testset (~line 132):
- Change `@test dist_bohr_trott < 0.02` to `@test dist_bohr_trott < 1e-5`
- This matches the pattern from quick task 5 where TrotterDomain thresholds were tightened.
  </action>
  <verify>
Run the full test suite from the project root:
```
cd /Users/bence/code/QuantumFurnace.jl && julia --project=. -e 'using Pkg; Pkg.test()'
```
All tests must pass, including:
- DMTST-05 with the tightened threshold (dist_bohr_trott < 1e-5)
- DMTST-06 (OFT consistency, should be unaffected)
- All other test files
  </verify>
  <done>
Both B_trotter() functions transform jump.in_eigenbasis to Trotter eigenbasis before computation.
DMTST-05 passes with threshold 1e-5 (down from 0.02). Full test suite green.
  </done>
</task>

</tasks>

<verification>
- DMTST-05 `dist_bohr_trott` value drops from ~0.011 to ~1e-8 or smaller
- DMTST-05 `dist_bohr_trott < 1e-5` passes
- DMTST-06 OFT tests unaffected (they already handle basis transformation at the call site)
- No regressions in any test file
</verification>

<success_criteria>
- B_trotter() single-jump and multi-jump both use jump operators transformed to Trotter eigenbasis
- DMTST-05 threshold tightened to 1e-5 and test passes
- Full test suite passes with no regressions
</success_criteria>

<output>
After completion, create `.planning/quick/6-fix-b-time-vs-b-trott-basis-mismatch-in-/6-SUMMARY.md`
</output>
