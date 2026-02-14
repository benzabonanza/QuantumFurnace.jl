# Coding Conventions

**Analysis Date:** 2026-02-13

## Naming Patterns

**Files:**
- Lowercase with underscores: `hamiltonian.jl`, `jump_workers.jl`, `coherent.jl`
- One main concept per file in most cases
- Grouped logically by domain (e.g., `bohr_domain.jl`, `energy_domain.jl`, `time_domain.jl`)
- Test files follow pattern: `{concept}_test.jl` (e.g., `trajectory_test.jl`, `ham_test.jl`)
- Main module file: `QuantumFurnace.jl`

**Functions:**
- Snake_case with trailing function indicators: `run_lindbladian`, `construct_lindbladian`, `precompute_data`
- Action verbs prefix: `create_*`, `compute_*`, `build_*`, `jump_contribution!`
- Abbreviations common in domain: `oft!` (Operator Fourier Transform), `psi`, `rho`, `dm` (density matrix)
- Query functions: `is_density_matrix`, `pick_transition`
- Mutating functions use `!` suffix: `vectorize_liouv_diss_and_add!`, `oft!`, `kron!`, `rmul!`

**Variables:**
- Greek letters used extensively: `beta`, `sigma`, `omega`, `nu`, `eta`, `gamma`, `delta`
- Math-inspired short names: `H` (Hamiltonian), `A` (jump operator), `B` (coherent term), `L` (Lindblad operator), `U` (unitary), `R` (dissipation matrix), `K` (Kraus operator)
- Subscripts indicated by underscore: `jump_oft`, `b_plus`, `b_minus`, `f_minus`, `f_plus`
- Loop indices: `i`, `j`, `k`, `l` (short loops), `site`, `pauli`, `trajectory` (descriptive loops)
- Configuration abbreviations: `nqb` (num_qubits), `w0`, `t0` (time/energy units)

**Types:**
- PascalCase: `HamHam`, `TrottTrott`, `JumpOp`, `LiouvConfig`, `ThermalizeConfig`
- Enum-like domain types: `BohrDomain`, `EnergyDomain`, `TimeDomain`, `TrotterDomain`
- Abstract types prefix with `Abstract`: `AbstractDomain`, `AbstractConfig`, `AbstractLiouvConfig`
- Result types: `HotAlgorithmResults`, `HotSpectralResults`
- Workspace/cache types: `LindbladianWorkspace`, `TrajectoryWorkspace`, `OFTCaches`, `KrausScratch`

## Code Style

**Formatting:**
- Indentation: 4 spaces
- Line length: Generally follows Julia convention (~92 characters)
- No explicit formatter configured (no `.prettierrc` or equivalent)
- Spacing around operators: `a + b`, `a * b`, `a / b`

**Linting:**
- No explicit linter configured
- Follows Julia style guidelines implicitly

**Documentation:**
- Triple-quoted docstrings for public functions: `""" ... """`
- Docstrings placed immediately before function definition
- Docstrings include: description, `# Fields`, `# Arguments`, implementation notes
- Examples from `src/structs.jl`:
  ```julia
  """
      LiouvConfig

      A configuration object that holds all the parameters...

      # Fields
      - `num_qubits`: The number of system qubits.
      - `with_coherent`: Option to add coherent term.
      ...
  """
  ```
- Internal functions use short comments instead of docstrings

## Import Organization

**Order:**
1. Package imports first: `using Pkg`, `using Base`, `using Printf`
2. External packages: `using LinearAlgebra`, `using SparseArrays`, `using BSON`, `using Arpack`
3. Specialized packages: `using FINUFFT`, `using QuadGK`, `using Optim`
4. Utilities: `using Random`, `using ProgressMeter`, `using Distributed`

**Module Export:**
- Module defines exports at top after imports: `export AbstractConfig, LiouvConfig, HamHam, ...`
- Grouped by functionality
- Public API clearly demarcated with `# --- Public API ---` comment
- Implementation includes listed with `include()` statements after exports

**Example from `src/QuantumFurnace.jl`:**
```julia
module QuantumFurnace
using Pkg
using Printf
using LinearAlgebra
...
# --- Public API ---
export AbstractConfig, AbstractLiouvConfig, ...
# --- Internal Implementation ---
include("constants.jl")
include("hamiltonian.jl")
...
end
```

## Error Handling

**Patterns:**
- `@assert` with error message for invariant checks: `@assert ishermitian(rescaled_hamiltonian) "..."`
- `throw(ArgumentError("..."))` for invalid inputs: Used in `is_density_matrix`, `validate_config!`
- Guard clauses with `===` checks: `trotter === nothing && error("A Trotter object must be provided")`
- Configuration validation via `validate_config!(config)` before major operations

**Example from `src/furnace.jl`:**
```julia
if config.domain isa TrotterDomain
    trotter === nothing && error("A Trotter object must be provided for the TrotterDomain")
    trotter
else
    hamiltonian
end
```

## Logging

**Framework:** `Printf.@printf` for output (no external logging library)

**Patterns:**
- Progress output: `@printf("Done.\n")`, `@printf("Worst quadrature error for the energy integral: %s\n", energy_error)`
- Debugging: `display()` and `norm()` calls in test/debug files
- BenchmarkTools: `@time` and `@btime` macros for performance measurement

**Example from test files:**
```julia
kraus_jumps = @time precompute_kraus_jumps(config.domain, jumps, hamiltonian, config, precomputed_data)
psi = @btime evolve_along_trajectory(psi0, fw, total_time)
```

## Comments

**When to Comment:**
- Physics/algorithm explanations: Comments describe mathematical concepts
- Non-obvious choices: Why a particular approximation or representation is used
- Performance notes: Explains in-place operations and memory optimization
- TODO/FIXME: Marked with `#!` or `#*` for visibility

**JSDoc/TSDoc:**
- Julia docstrings use `"""` format
- Structured fields described in `# Fields` section
- Parameters in `# Arguments` section (when needed)

**Style:**
- Comments prefixed with `#` or `#*` (for section markers)
- Complex algorithms documented with step-by-step comments
- Example from `src/qi_tools.jl`:
  ```julia
  """
      Computes C .+= alpha .* kron(A, B) completely in-place, without allocating
      the result of the Kronecker product. Speed.
  """
  ```

## Function Design

**Size:**
- Most functions 30-100 lines (balance between modularity and readability)
- Helper functions factored out for reusability
- Some domain-specific functions can be longer (100-200 lines) due to algorithm complexity

**Parameters:**
- Use type annotations for dispatch: `jump::JumpOp`, `hamiltonian::HamHam`
- Domain dispatch via type: Methods specialized by `AbstractDomain` subtype
- Keyword arguments for optional parameters: `hermitian_check = false`, `do_adjoint::Bool=false`
- Union types for multiple accepted types: `Union{TrottTrott, Nothing}`, `Union{HamHam, TrottTrott}`

**Return Values:**
- Explicit return statements for non-trivial returns
- In-place functions return modified argument or `nothing`
- Structured return types for complex results: `HotSpectralResults(data=..., fixed_point=...)`

**Example from `src/hamiltonian.jl`:**
```julia
function HamHam(terms::Vector{Vector{Matrix{ComplexF64}}}, coeffs::Vector{Float64}, num_qubits::Int64;
    periodic::Bool = true, hermitian_check = false)
    # Implementation
    return hamiltonian
end
```

## Module Design

**Exports:**
- All public types exported explicitly: Types, main functions
- Private implementation functions not exported
- Two-step structure: public API section, then internal includes

**Barrel Files:**
- Main module file `src/QuantumFurnace.jl` acts as barrel, importing and re-exporting
- Each domain (bohr, energy, time, trotter) has dedicated file
- Utility files: `constants.jl`, `qi_tools.jl`, `misc_tools.jl`

**Example structure:**
```julia
# src/QuantumFurnace.jl
export AbstractConfig, LiouvConfig, HamHam, run_lindbladian, ...
include("constants.jl")
include("structs.jl")
include("hamiltonian.jl")
include("bohr_domain.jl")
```

## Dispatch Patterns

**Multiple Dispatch:**
- Heavy use of dispatch on domain type: `jump_contribution!(::BohrDomain, ...)`, `jump_contribution!(::EnergyDomain, ...)`
- Dispatch on config type: `create_alpha(...)` functions have multiple implementations
- Dispatch on struct type: `HamHam` vs `TrottTrott` for basis selection

**Example from `src/jump_workers.jl`:**
```julia
function jump_contribution!(
    L_target::AbstractMatrix{ComplexF64},
    ::BohrDomain,
    jump::JumpOp,
    hamiltonian::HamHam,
    ...
)
    # Bohr-domain specific implementation
end

function jump_contribution!(
    L_target::AbstractMatrix{ComplexF64},
    ::EnergyDomain,
    jump::JumpOp,
    ...
)
    # Energy-domain specific implementation
end
```

## Type Annotations

**Usage:**
- Function arguments typed: `jump::JumpOp`, `beta::Float64`, `dim::Int`
- Return type annotations sometimes used but not consistently: `-> Matrix{ComplexF64}`
- Type-parameterized structs: `struct TrajectoryWorkspace{T}`
- Union types for optional parameters: `Union{Matrix{T}, Nothing}`

**Example from `src/trajectories.jl`:**
```julia
struct TrajectoryWorkspace{T}
    jump_oft::Matrix{T}
    psi_tmp::Vector{T}
    Rpsi::Vector{T}
end

function TrajectoryWorkspace(::Type{T}, dim::Int) where {T}
    return TrajectoryWorkspace{T}(
        zeros(T, dim, dim),
        zeros(T, dim),
        zeros(T, dim),
    )
end
```

---

*Convention analysis: 2026-02-13*
