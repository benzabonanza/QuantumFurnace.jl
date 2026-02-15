---
phase: 16-defer-lindbladianworkspace-construction
plan: 16
type: execute
wave: 1
depends_on: []
files_modified: [src/structs.jl]
autonomous: true
must_haves:
  truths:
    - "LindbladianWorkspace cannot be constructed without an explicit type parameter"
    - "All existing code continues to work (the convenience constructor was unused)"
  artifacts:
    - path: "src/structs.jl"
      provides: "LindbladianWorkspace{T} struct without Float64 default constructor"
  key_links: []
---

<objective>
Remove the eager default constructor `LindbladianWorkspace(dim::Int) = LindbladianWorkspace{Float64}(dim)` from structs.jl.

Purpose: This convenience constructor allows `LindbladianWorkspace` to be constructed without specifying a type parameter, defaulting to Float64. The user wants to ensure LindbladianWorkspace is only constructed explicitly when needed for Lindbladian simulations, never eagerly or accidentally. The only actual construction site (furnace.jl:69) already uses `LindbladianWorkspace{T}(dim)` with an explicit type parameter, so this convenience constructor is dead code.

Output: structs.jl with the convenience constructor removed.
</objective>

<execution_context>
@./.claude/get-shit-done/workflows/execute-plan.md
@./.claude/get-shit-done/templates/summary.md
</execution_context>

<context>
@.planning/STATE.md
@src/structs.jl
@src/furnace.jl
</context>

<tasks>

<task type="auto">
  <name>Task 1: Remove LindbladianWorkspace default constructor</name>
  <files>src/structs.jl</files>
  <action>
Delete line 54 from src/structs.jl:
```
LindbladianWorkspace(dim::Int) = LindbladianWorkspace{Float64}(dim)
```

This is the only line to change. The convenience constructor is not used anywhere in src/ or test/. The sole construction site in furnace.jl:69 already uses `LindbladianWorkspace{T}(dim)` with an explicit type parameter derived from `eltype(hamiltonian.eigvals)`.

Do NOT modify any other files. Do NOT change the `LindbladianWorkspace{T}(dim)` inner constructor -- that stays.
  </action>
  <verify>
1. `grep -n "LindbladianWorkspace(dim" src/structs.jl` returns no matches (the convenience constructor is gone)
2. `grep -n "LindbladianWorkspace{T}" src/structs.jl` still shows the inner constructor (line ~44)
3. Run full test suite: `cd /Users/bence/code/QuantumFurnace.jl && julia --project -e 'using Pkg; Pkg.test()'` -- all tests pass
  </verify>
  <done>
The convenience constructor `LindbladianWorkspace(dim::Int)` no longer exists. The parameterized constructor `LindbladianWorkspace{T}(dim)` remains. All tests pass.
  </done>
</task>

</tasks>

<verification>
- `LindbladianWorkspace(dim)` (without type parameter) is no longer callable from within the package
- `LindbladianWorkspace{Float64}(dim)` (with explicit type parameter) still works
- All 224+ tests pass
</verification>

<success_criteria>
Line 54 of structs.jl (`LindbladianWorkspace(dim::Int) = LindbladianWorkspace{Float64}(dim)`) is deleted. No other changes needed. Full test suite passes.
</success_criteria>

<output>
After completion, create `.planning/quick/16-defer-lindbladianworkspace-construction-/16-SUMMARY.md`
</output>
