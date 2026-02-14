# Testing Patterns

**Analysis Date:** 2026-02-13

## Test Framework

**Runner:**
- Built-in Julia `Test` module (from `Project.toml` extras)
- Not using external test frameworks like Pytest or Jest

**Run Commands:**
```bash
julia --project -e "using Pkg; Pkg.test()"
```

**Note:** Currently, test suite appears to use interactive evaluation rather than automated `@test` assertions (see Test File Organization below).

## Test File Organization

**Location:**
- Dedicated `/test` directory: `test/` folder at project root
- Interactive test files, not automated test suites

**Naming:**
- Pattern: `{concept}_test.jl`
- Examples:
  - `test/B_test.jl`
  - `test/trajectory_test.jl`
  - `test/ham_test.jl`
  - `test/time_tests.jl`
  - `test/kossakowski_test.jl`
  - `test/trott_test.jl`
  - `test/log_sobolev_test.jl`

**File Count:** 7 test files

## Test Structure

**Pattern:**
Test files are interactive Julia scripts with setup, computation, and display/verification rather than automated test suites. They follow a consistent structure:

```julia
using Revise

if !isdefined(Main, :QuantumFurnace)
    includet("../src/QuantumFurnace.jl")
end

using .QuantumFurnace
using LinearAlgebra, Random

#* Config
num_qubits = 4
dim = 2^num_qubits
beta = 10.
...

#* Hamiltonian
hamiltonian = load_hamiltonian("heis", num_qubits)
hamiltonian = finalize_hamham(hamiltonian, beta)

#* Jumps
jump_paulis = [[X], [Y], [Z]]
...

#* Computation
result = compute_something()
display(result)

# Verification (manual inspection)
norm(expected - actual)
```

**Key Characteristics:**

1. **Module Loading:** Always includes conditional loading with Revise for development
   ```julia
   using Revise
   if !isdefined(Main, :QuantumFurnace)
       includet("../src/QuantumFurnace.jl")
   end
   using .QuantumFurnace
   ```

2. **Section Organization:** Configuration, setup, computation marked with `#*` comment header
   ```julia
   #* Config
   #* Hamiltonian
   #* Jumps
   #* Computation
   ```

3. **Display vs Assert:** Results displayed with `display()` or computed norms rather than `@test` assertions
   ```julia
   norm(B_bohr - B_t)           # Manual verification
   display(norm(...))           # Display for inspection
   ```

4. **Parameter Exploration:** Tests often contain commented-out alternative configurations
   ```julia
   # beta = 5  # or 10, 30
   # a = 0.0   # Alternative parameter
   ```

## Benchmarking

**Included Imports:**
- `BenchmarkTools` in some tests for performance measurement

**Patterns:**
```julia
using BenchmarkTools

# Time single evaluation
kraus_jumps = @time precompute_kraus_jumps(config.domain, jumps, hamiltonian, config, precomputed_data)

# Benchmark repeated evaluation
psi = @btime evolve_along_trajectory(psi0, fw, total_time)
```

**Used in:**
- `test/trajectory_test.jl`
- `test/kossakowski_test.jl`

## Test Fixtures and Data

**Test Data:**
- Hamiltonians loaded from BSON files in `hamiltonians/` directory
- Pattern in tests:
  ```julia
  hamiltonian = load_hamiltonian("heis", num_qubits)
  ```

**Configuration Objects:**
Tests create full `LiouvConfig` or `ThermalizeConfig` objects to define test scenarios:

```julia
config = LiouvConfig(
    num_qubits = num_qubits,
    with_coherent = with_coherent,
    with_linear_combination = with_linear_combination,
    domain = domain,
    beta = beta,
    sigma = sigma,
    gaussian_parameters = (w_gamma, sigma_gamma),
    a = a,
    b = b,
    num_energy_bits = num_energy_bits,
    w0 = w0,
    t0 = t0,
    eta = eta,
    num_trotter_steps_per_t0 = num_trotter_steps_per_t0
)
```

**Hamiltonian Creation:**
```julia
# Method 1: Load from file
hamiltonian = load_hamiltonian("heis", num_qubits)

# Method 2: Create programmatically
hamiltonian_terms = [[X, X], [Y, Y], [Z, Z]]
hamiltonian_coeffs = fill(1.0, length(hamiltonian_terms))
hamiltonian = HamHam(hamiltonian_terms, hamiltonian_coeffs, num_qubits)
```

**Jump Operator Creation:**
```julia
jump_paulis = [[X], [Y], [Z]]
jump_sites = 1:num_qubits
jumps::Vector{JumpOp} = []
for pauli in jump_paulis
    for site in jump_sites
        jump_op = pad_term(pauli, num_qubits, site) ./ jump_normalization
        basis_unitary = hamiltonian.eigvecs
        jump_op_in_eigenbasis = basis_unitary' * jump_op * basis_unitary

        jump = JumpOp(jump_op, jump_op_in_eigenbasis, orthogonal, hermitian)
        push!(jumps, jump)
    end
end
```

## Test Types

**What Tests Cover:**

1. **Configuration & Setup Tests** (e.g., `ham_test.jl`, `trott_test.jl`)
   - Verify Hamiltonian creation and comparison
   - Test configuration object building
   - Check helper function outputs

2. **Domain Transformation Tests** (e.g., `time_tests.jl`)
   - Verify consistency across domains (Bohr, Energy, Time, Trotter)
   - Test quadrature errors and OFT (Operator Fourier Transform)
   - Compare results from different approximation levels

3. **Algorithm Verification Tests** (e.g., `trajectory_test.jl`, `B_test.jl`)
   - Run thermalization/trajectory evolution
   - Compare coherent term B across implementations
   - Verify jump operator contributions

4. **Mathematical Validation** (e.g., `kossakowski_test.jl`)
   - Compare alpha functions across parameter sets
   - Verify symmetry properties (skew-symmetry checks)

## Verification Patterns

**Manual Verification:**
- Norm differences: `norm(expected - actual)` for approximate equality
- Zero checks: Verify small norms indicate correctness
- Display with `display()` for manual inspection

**Example from `test/B_test.jl`:**
```julia
B_bohr = coherent_bohr(hamiltonian, jump, config)
rmul!(B_bohr, precomputed_data.gamma_norm_factor)

B_t = B_time(jump, hamiltonian, precomputed_data.b_minus,
             precomputed_data.b_plus, config.t0, config.beta, config.sigma)
rmul!(B_t, precomputed_data.gamma_norm_factor)

norm(B_bohr - B_t)  # Should be small (near zero)
```

## Missing Automated Test Suite

**Current State:**
- No `@testset` blocks in test files
- No formal `@test` assertions
- Tests are interactive/exploratory rather than automated
- Results are visually inspected or norms computed manually

**CI Configuration:**
- Documentation workflow exists: `.github/workflows/Documentation.yml`
- No dedicated test CI workflow configured
- Tests are currently manual/interactive only

## Test Coverage

**Untested Areas:**
- No automated regression testing
- Edge cases and error conditions not systematically tested
- Configuration validation (`validate_config!`) used but not exhaustively tested
- Error handling paths not covered by automated tests

**Tested Areas (via interactive tests):**
- Core functionality: Liouvillian construction, B term computation
- Domain transformations: Consistency checks across approximation levels
- Hamiltonian operations: Creation, finalization, Gibbs state computation
- Jump operators: Creation, basis transformations

## Test Fixtures Location

**Data Files:**
- `hamiltonians/` - BSON files with precomputed Hamiltonian data
- Pattern: `{type}_disordered_periodic_n{num_qubits}.bson`
- Example: `heis_disordered_periodic_n4.bson`

## Common Test Patterns

**Configuration Setup:**
```julia
# Physical parameters
num_qubits = 4
dim = 2^num_qubits
beta = 10.0
sigma = 1 / beta

# Domain and approximation parameters
domain = TimeDomain()
num_energy_bits = 12
w0 = 0.05
t0 = 2pi / (2^num_energy_bits * w0)
num_trotter_steps_per_t0 = 10

# Linear combination parameters
a = beta / 30.0
b = 0.4
eta = 0.0

# Create configuration
config = LiouvConfig(
    num_qubits = num_qubits,
    with_coherent = true,
    with_linear_combination = true,
    domain = domain,
    beta = beta,
    sigma = sigma,
    a = a,
    b = b,
    num_energy_bits = num_energy_bits,
    w0 = w0,
    t0 = t0,
    eta = eta,
    num_trotter_steps_per_t0 = num_trotter_steps_per_t0
)
```

**Precomputation Pattern:**
```julia
precomputed_data = precompute_data(config.domain, config)
time_oft_caches = OFTCaches(dim)
```

**Comparison Pattern:**
```julia
# Compute both versions
result1 = approach_A(input)
result2 = approach_B(input)

# Verify consistency
error = norm(result1 - result2)
println("Difference: $error")  # Should be near machine epsilon
```

## Debugging Utilities

**Revise Integration:**
Tests use `Revise.jl` for interactive development:
```julia
using Revise
if !isdefined(Main, :QuantumFurnace)
    includet("../src/QuantumFurnace.jl")
end
```

**BenchmarkTools:**
For performance measurement:
```julia
using BenchmarkTools
result = @btime expensive_function(args)
```

**Plotting (when used):**
```julia
using Plots
# display(plot(...))  # Often commented out in test files
```

---

*Testing analysis: 2026-02-13*
