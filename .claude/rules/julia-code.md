---
paths:
  - "src/**/*.jl"
  - "test/**/*.jl"
  - "experiments/**/*.jl"
  - "simulations/**/*.jl"
  - "hamiltonians/**/*.jl"
---

# Julia Code Rules

## Coding Conventions

- Parametric types: `T` consistently, no hardcoded `Float64` in generic code
- Immutable structs for result types
- Numerical tests use `atol`/`rtol`, never exact equality for floats
- Library code lives in `src/`; read nearby code before writing — match existing patterns
- Use existing deps: `LinearAlgebra`, `SparseArrays`, `LsqFit`, `Roots` — don't add new deps without justification
- Performance: preallocate, `@views` for slices, avoid allocations in hot loops
- No `Any` in struct fields or hot-path return types

## Architecture

- **Domain hierarchy**: `BohrDomain` -> `EnergyDomain` -> `TimeDomain` -> `TrotterDomain`
- **Core structs**: `HamHam` (Hamiltonian + disorder), `TrottTrott` (Trotterized), `Workspace`
- **Fitting**: `FitResult`, `BiexpFitResult`, `MixingTimeEstimate` — all immutable, via LsqFit.jl
- **Simulation**: Krylov (sparse, large n) vs dense (exact, small n)
- **Disordering**: `find_ideal_heisenberg` with `disordering_terms` kwarg; Z-only or Z+ZZ for symmetry breaking

## Atomic Commits

Every code change MUST be an atomic git commit. One logical change per commit.

- Each commit: single, self-contained change (one new function, one bug fix, one struct, one test file)
- Tests must pass after each commit
- Commit message format: `<type>(<scope>): <description>`
- Stage only the specific files — never `git add .` or `git add -A`
- Include `Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>`
- `src/` + `test/` in one commit is fine if the test is specifically for that change

## Physics Sanity Rule

Whenever code makes a decision based on physical or mathematical reasoning — parameter choices, expected scalings, algorithm selection, convergence criteria — spawn the `physics-checker` agent to verify. Flag any `# PHYSICS CHECK:` comments for review.
