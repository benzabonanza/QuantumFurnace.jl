---
phase: 21-fix-test-errors-after-removing-transform
plan: 01
type: execute
wave: 1
depends_on: []
files_modified:
  - src/QuantumFurnace.jl
  - test/test_allocation.jl
autonomous: true
must_haves:
  truths:
    - "All 539 tests pass (zero failures, zero errors)"
    - "transform_jumps_to_basis is fully removed from exports and test imports"
    - "B_trotter allocation test uses pre-built TEST_TROTTER_JUMPS from test_helpers.jl"
  artifacts:
    - path: "src/QuantumFurnace.jl"
      provides: "Module exports without transform_jumps_to_basis"
    - path: "test/test_allocation.jl"
      provides: "Allocation tests using TEST_TROTTER_JUMPS"
  key_links: []
---

<objective>
Remove all remaining references to `transform_jumps_to_basis` after its deletion from `src/qi_tools.jl`.

Purpose: The function was removed but two files still reference it -- the module export list and one test file -- causing 1 test failure (Aqua undefined exports) and 1 test error (missing function call).

Output: Clean test suite with all 539 tests passing.
</objective>

<execution_context>
@./.claude/get-shit-done/workflows/execute-plan.md
@./.claude/get-shit-done/templates/summary.md
</execution_context>

<context>
@src/QuantumFurnace.jl
@test/test_allocation.jl
@test/test_helpers.jl
</context>

<tasks>

<task type="auto">
  <name>Task 1: Remove transform_jumps_to_basis from exports and fix allocation test</name>
  <files>src/QuantumFurnace.jl, test/test_allocation.jl</files>
  <action>
1. In `src/QuantumFurnace.jl` line 56: remove `, transform_jumps_to_basis` from the export statement. The line should become:
   ```
   hermitianize!
   ```
   (just `hermitianize!` ending the qi_tools export block).

2. In `test/test_allocation.jl` line 18: remove `transform_jumps_to_basis,` from the `using QuantumFurnace:` import. The import block (lines 16-19) should become:
   ```julia
   using QuantumFurnace: B_bohr, B_time, B_trotter,
                         _precompute_data, _jump_contribution!,
                         KrausScratch,
                         TrajectoryWorkspace
   ```

3. In `test/test_allocation.jl` lines 86-87: replace the `transform_jumps_to_basis` call with the pre-built `TEST_TROTTER_JUMPS` constant (already defined in `test/test_helpers.jl` line 217, constructed via `make_test_system(; trotter=TEST_TROTTER).jumps` which uses `trotter.eigvecs` as the basis -- identical transform). Change:
   ```julia
   trotter_jumps = transform_jumps_to_basis(TEST_JUMPS, TEST_TROTTER.eigvecs)
   jump = trotter_jumps[1]
   ```
   to:
   ```julia
   trotter_jumps = TEST_TROTTER_JUMPS
   jump = trotter_jumps[1]
   ```

4. Remove the now-stale comment on line 85 ("Need Trotter-basis jumps (callers now transform before calling B_trotter)") or update it to: `# Use pre-built Trotter-basis jumps from test_helpers.jl`.
  </action>
  <verify>
Run the full test suite:
```bash
cd /Users/bence/code/QuantumFurnace.jl && julia --project=. -e 'using Pkg; Pkg.test()'
```
Expected: all 539 tests pass, 0 failures, 0 errors.
  </verify>
  <done>
- `transform_jumps_to_basis` does not appear in `src/QuantumFurnace.jl` exports
- `transform_jumps_to_basis` does not appear in `test/test_allocation.jl`
- B_trotter allocation test uses `TEST_TROTTER_JUMPS` and passes
- Full test suite passes with zero failures and zero errors
  </done>
</task>

</tasks>

<verification>
```bash
# Verify no remaining references to transform_jumps_to_basis in source or test code
grep -r "transform_jumps_to_basis" src/ test/

# Run full test suite
cd /Users/bence/code/QuantumFurnace.jl && julia --project=. -e 'using Pkg; Pkg.test()'
```
</verification>

<success_criteria>
- grep returns no matches for transform_jumps_to_basis in src/ or test/
- All 539 tests pass with zero failures and zero errors
</success_criteria>

<output>
After completion, create `.planning/quick/21-fix-test-errors-after-removing-transform/21-SUMMARY.md`
</output>
