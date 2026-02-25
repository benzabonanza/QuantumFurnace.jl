# Coding Conventions

**Analysis Date:** 2026-02-25

## Naming Patterns

**Files:**
- `snake_case.jl` for all source files: `krylov_matvec.jl`, `gap_estimation.jl`, `jump_workers.jl`
- Module entry point: `src/QuantumFurnace.jl`
- One logical concern per file (domains, structs, trajectories, fitting, etc.)

**Functions:**
- `snake_case` for all functions: `run_lindbladian`, `build_trajectoryframework`, `estimate_spectral_gap`
- Mutating functions end in `!`: `hermitianize!`, `apply_lindbladian!`, `step_along_trajectory!`, `_kron!`
- Internal/private functions prefixed with `_`: `_precompute_data`, `_jump_contribution!`, `_gibbs_in_eigen`, `_load_hamiltonian_bson`
- Inline hot-path functions use `@inline`: `_krylov_oft!`

**Types / Structs:**
- `PascalCase` for all types: `HamHam`, `TrottTrott`, `LiouvConfig`, `JumpOp`, `KrylovWorkspace`
- Quirky abbreviations are used for domain-specific types: `HamHam` (Hamiltonian), `TrottTrott` (Trotter)
- Abstract types prefixed with `Abstract`: `AbstractConfig`, `AbstractDomain`, `AbstractLiouvConfig`
- Result structs suffixed with `Result`: `LindbladianResult`, `SpectralGapResult`, `FitResult`, `KrylovGapResult`
- Config structs suffixed with `Config`: `LiouvConfig`, `ThermalizeConfig`, `LiouvConfigGNS`
- Workspace structs suffixed with `Workspace`: `KrylovWorkspace`, `TrajectoryWorkspace`, `LindbladianWorkspace`
- GNS variant structs suffixed with `GNS`: `LiouvConfigGNS`, `ThermalizeConfigGNS`

**Variables:**
- `snake_case` for local variables: `num_qubits`, `jump_normalization`, `eigvals_vec`
- Physics/math variables use short conventional names: `rho`, `psi`, `beta`, `sigma`, `delta`, `eigvals`, `eigvecs`
- Scratch/temporary matrices use descriptive names: `jump_tmp`, `jump_conj`, `jump_dag_jump`, `rho_acc`
- Constants are `UPPER_SNAKE_CASE`: `TOL_EXACT`, `NUM_QUBITS`, `BETA`, `TEST_DELTA`

**Constants in source code:**
- Module-level constants for index aliases: `const _IDX_A = 1`, `const _IDX_GAP = 2`

## Code Style

**Formatting:**
- No autoformatter is configured (no `.prettierrc`, no JuliaFormatter config detected)
- Standard Julia style: 4-space indentation
- Spaces around operators and after commas
- Keyword arguments use `=` without spaces in function signatures

**Linting:**
- Aqua.jl is used for package quality checks (`test/test_aqua.jl`): `Aqua.test_all(QuantumFurnace; ambiguities=false, piracies=false)`

## Type Annotations and Parametric Programming

**Type parameters on structs:**
- Structs are parametric on element type `T<:AbstractFloat`: `HamHam{T}`, `LindbladianWorkspace{T}`
- Complex matrices use `Complex{T}` parametrically: `Matrix{Complex{T}}`
- Abstract type hierarchies used for dispatch: `AbstractDomain`, `AbstractConfig{D,T}`

**Function signatures:**
- Type-constrained signatures for performance-critical functions: `where {T<:Complex}`
- Use `Union{X, Nothing}` for optional fields: `Union{Nothing, Matrix{T}}`, `Union{TrottTrott, Nothing}`
- Keyword arguments with default `nothing` for optional parameters

**Type annotations in structs:**
- All struct fields explicitly typed; no untyped fields except where `Any` is required (e.g. `precomputed_data::Any` for NamedTuple variability)
- `@kwdef` macro used heavily to allow keyword-argument construction: `@kwdef struct LiouvConfig{D,T} ...`

## Import Organization

**Order in source files:**
1. Standard library imports (`using LinearAlgebra`, `using Random`, `using Printf`)
2. Third-party package imports (`using BSON`, `using Arpack`, `using KrylovKit`)
3. No relative imports — files included via `include()` calls in `src/QuantumFurnace.jl`

**In test files:**
1. `using Test` first
2. Then any additional packages needed
3. No need to re-import `QuantumFurnace` (done in `runtests.jl`)
4. Comment stating `# test_helpers.jl is already included by runtests.jl`

**Path Aliases:**
- Not used; all cross-file references go through the module namespace or are `include()`d

## Error Handling

**Patterns:**
- Use `error(message)` for domain invariant violations: `error("GNS configs must have with_coherent=false")`
- Use `@assert` for runtime precondition checks: `@assert delta < 1.0 "delta = $(delta) >= 1.0: too large for CPTP channel"`
- Use `throw(ArgumentError(...))` for invalid input arguments: `throw(ArgumentError("Expected $(Complex{T}) term data..."))`
- Constructors enforce invariants via outer constructor methods (not inner)
- No exception types beyond base Julia: no custom exception types

## Logging and Debug Output

**Framework:** `@printf` / `println` (no structured logging library)

**Patterns:**
- Use `@printf("Done.\n")` for progress in long computations
- Use `@info` in tests for diagnostic values: `@info "DMTST-01: Bohr fixed point trace distance to Gibbs" dist`
- Commented-out debug prints are left in-place (e.g. `# @printf("Gaussian\n")`, `# @printf("Smooth Metro\n")`)
- `ProgressMeter` is imported but usage is in long-running scripts

## Comments

**When to Comment:**
- All public functions get docstrings in `"""..."""` format with description, parameters, returns
- Internal helpers get brief docstrings explaining purpose
- Section delimiters use `# ---` or `# ===` banners with descriptive text
- Physics/math derivations are commented inline: formula references, convention notes
- Struct fields documented with `#` inline comments listing buffer purpose
- Phase/plan references included as comments: `# Phase 27: Core Matvec Infrastructure`, `# TVAL-01`

**Docstring format:**
```julia
"""
    function_name(arg1, arg2) -> ReturnType

Brief description.

More details if needed.

# Arguments
- `arg1`: Description
"""
```

**Inline annotation style:**
```julia
# Compute: out[i,j] = eigenbasis[i,j] * exp(-(energy - bohr_freqs[i,j])^2 * inv_4sigma2)
```

## Function Design

**Size:** Functions are kept focused; large files are split by logical concern. Hot-path functions are kept minimal.

**Parameters:**
- Positional arguments for primary data, keyword arguments for options with defaults
- Workspace/scratch objects passed explicitly rather than created inside hot-path functions
- Config structs bundle related parameters rather than long positional argument lists

**Return Values:**
- Functions return named tuples `(; field1, field2)` for lightweight multi-value returns
- Heavy results use dedicated Result structs (`LindbladianResult`, `KrylovGapResult`, etc.)
- Mutating functions (`!`) return `nothing` by convention (except `hermitianize!` which returns `A`)
- Non-mutating hot-path functions return their output matrix

## Module Design

**Exports:** All public API is explicitly listed in `src/QuantumFurnace.jl` under grouped `export` statements with comments (Types, Domains, Trajectory, Results, etc.)

**Internal naming convention:** Internal symbols prefixed with `_` are not exported but can be accessed as `QuantumFurnace._precompute_data(...)` in tests when needed.

**Barrel Files:** Not used; a single module file (`src/QuantumFurnace.jl`) includes all source files via `include()`.

---

*Convention analysis: 2026-02-25*
