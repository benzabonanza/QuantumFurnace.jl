---
phase: 13-unify-residual-cholesky
plan: 01
type: execute
wave: 1
depends_on: []
files_modified: [src/jump_workers.jl]
autonomous: true

must_haves:
  truths:
    - "apply_cptp_channel! uses eigendecomposition with clamped eigenvalues instead of Cholesky"
    - "No cholesky! call remains in jump_workers.jl"
    - "No eps_shift diagonal hack remains in jump_workers.jl"
    - "All existing tests pass with identical or improved numerical results"
  artifacts:
    - path: "src/jump_workers.jl"
      provides: "Eigendecomposition-based U_residual computation in apply_cptp_channel!"
      contains: "eigen"
  key_links:
    - from: "apply_cptp_channel!"
      to: "eigen(Hermitian(scratch.tmp2))"
      via: "eigendecomposition with clamped negative values"
      pattern: "eig\\.values .= max"
---

<objective>
Replace the fragile Cholesky-based U_residual computation in apply_cptp_channel! (jump_workers.jl)
with the numerically robust eigendecomposition approach already used in trajectories.jl.

Purpose: The current `cholesky!(Hermitian(scratch.tmp2), check=false)` silently produces NaN/garbage
when the matrix S is not PSD despite the eps_shift heuristic. The eigendecomposition approach handles
all negative eigenvalue cases by clamping them to zero, which is mathematically correct (nearest PSD
square root).

Output: Updated apply_cptp_channel! function using eigendecomposition.
</objective>

<execution_context>
@./.claude/get-shit-done/workflows/execute-plan.md
@./.claude/get-shit-done/templates/summary.md
</execution_context>

<context>
@src/jump_workers.jl (lines 160-216: apply_cptp_channel! function)
@src/trajectories.jl (lines 112-117: eigendecomposition approach to copy)
@src/qi_tools.jl (hermitianize! helper)
</context>

<tasks>

<task type="auto">
  <name>Task 1: Replace Cholesky with eigendecomposition in apply_cptp_channel!</name>
  <files>src/jump_workers.jl</files>
  <action>
In src/jump_workers.jl, modify the apply_cptp_channel! function (lines 174-216).

**Replace lines 195-202** (the eps_shift + cholesky block):
```julia
    # Guard against tiny negative eigenvalues from roundoff (S is O(delta^2))
    eps_shift = 10 * eps(Float64)
    @inbounds for i in 1:dim
        scratch.tmp2[i,i] += eps_shift
    end

    cholesky_S = cholesky!(Hermitian(scratch.tmp2), check=false)
    U_residual = cholesky_S.U
```

**With the eigendecomposition approach** (matching trajectories.jl pattern):
```julia
    # PSD guard: clamp negative eigenvalues to zero (more robust than Cholesky + eps shift)
    hermitianize!(scratch.tmp2)
    S_herm = Hermitian(scratch.tmp2)
    eig = eigen(S_herm)
    eig.values .= max.(eig.values, 0.0)
    U_residual = Matrix{ComplexF64}(Diagonal(sqrt.(eig.values)) * eig.vectors')
```

The hermitianize! call on scratch.tmp2 is needed because scratch.tmp2 holds S (computed from
arithmetic on R and R^2), which may have tiny imaginary asymmetries from floating-point arithmetic
even though R itself was hermitianized before entry. This matches the trajectories.jl pattern.

**Also update the docstring** (line 168): change
  `U_residual = cholesky(S + eps*I).U`
to
  `U_residual = sqrt_psd(S)  (eigendecomposition with clamped eigenvalues)`

Everything else in the function remains unchanged -- the K0 computation, the rho_next sandwich
products, the final hermitianize! and copyto! are all untouched.
  </action>
  <verify>
Run the full test suite:
```bash
cd /Users/bence/code/QuantumFurnace.jl && julia --project=. -e 'using Pkg; Pkg.test()'
```

Additionally verify the code changes:
- `grep -c 'cholesky' src/jump_workers.jl` should return 0
- `grep -c 'eps_shift' src/jump_workers.jl` should return 0
- `grep -c 'eigen' src/jump_workers.jl` should return 1 (inside apply_cptp_channel!)
- `grep -c 'hermitianize!(scratch.tmp2)' src/jump_workers.jl` should return 1
  </verify>
  <done>
apply_cptp_channel! uses eigendecomposition with clamped eigenvalues. No cholesky! or eps_shift
references remain in jump_workers.jl. All tests pass. The DM simulator and trajectory simulator
now use the same numerically robust approach for computing U_residual.
  </done>
</task>

</tasks>

<verification>
1. All existing tests pass (the test suite exercises all 3 domain-specific jump_contribution! methods
   that call apply_cptp_channel!)
2. No cholesky or eps_shift references remain in jump_workers.jl
3. The eigendecomposition pattern in jump_workers.jl matches the one in trajectories.jl
</verification>

<success_criteria>
- apply_cptp_channel! uses eigendecomposition (eigen + clamp + sqrt) instead of Cholesky
- Zero references to cholesky! or eps_shift in jump_workers.jl
- Full test suite passes
</success_criteria>

<output>
After completion, create `.planning/quick/13-unify-residual-cholesky-computation-comp/13-SUMMARY.md`
</output>
