---
phase: 15-remove-unused-build-common-fields-helper
plan: 01
type: execute
wave: 1
depends_on: []
files_modified: [src/structs.jl]
autonomous: true

must_haves:
  truths:
    - "_build_common_fields function no longer exists in codebase"
    - "All existing tests still pass (no behavioral change)"
  artifacts:
    - path: "src/structs.jl"
      provides: "Config structs without dead helper function"
  key_links: []
---

<objective>
Remove the unused `_build_common_fields()` helper function and its docstring from `src/structs.jl` (lines 56-81).

Purpose: This function was created during phase 08-01 as a shared constructor helper for config types, but it is never called anywhere. It cannot actually simplify construction because `@kwdef` structs auto-generate keyword constructors and GNS manual constructors need explicit kwargs. Removing it eliminates dead code.

Output: Cleaner `src/structs.jl` with no dead helper function.
</objective>

<execution_context>
@./.claude/get-shit-done/workflows/execute-plan.md
@./.claude/get-shit-done/templates/summary.md
</execution_context>

<context>
@src/structs.jl
</context>

<tasks>

<task type="auto">
  <name>Task 1: Remove _build_common_fields() dead code</name>
  <files>src/structs.jl</files>
  <action>
    Remove lines 56-81 from src/structs.jl. These lines contain the docstring (lines 56-61) and the function definition (lines 63-81) for `_build_common_fields()`.

    After removal, there will be two consecutive blank lines (the original line 55 after `AbstractThermalizeConfig` and the original line 82 before the `LiouvConfig` comment). Collapse these to a single blank line so the file reads cleanly:

    ```
    abstract type AbstractThermalizeConfig{D<:AbstractDomain} <: AbstractConfig{D} end

    # Let's keep this structure...
    ```

    Do NOT modify any other code. No other files reference this function in actual code.
  </action>
  <verify>
    1. `grep -n "_build_common_fields" src/structs.jl` returns no results
    2. Run the full test suite to confirm no behavioral change (use whatever test runner the project uses -- look for Julia test commands)
  </verify>
  <done>
    The `_build_common_fields()` function and its docstring are completely removed from `src/structs.jl`. All tests pass unchanged.
  </done>
</task>

</tasks>

<verification>
- `grep -rn "_build_common_fields" src/` returns zero results
- Test suite passes with no failures
</verification>

<success_criteria>
- Dead function removed from src/structs.jl
- No references to _build_common_fields remain in source code (planning docs are fine)
- All tests pass
</success_criteria>

<output>
After completion, create `.planning/quick/15-remove-unused-build-common-fields-helper/15-SUMMARY.md`
</output>
