# Phase 11: Allocation Optimization - Research

**Researched:** 2026-02-15
**Domain:** Julia performance optimization -- eliminating heap allocations in hot-path numerical kernels
**Confidence:** HIGH

## Summary

Phase 11 targets four distinct allocation hotspots in the core simulation paths: sparse matrix construction in `B_bohr` (`bohr_domain.jl`), `Diagonal` wrapper creation in `B_time`/`B_trotter` (`coherent.jl`), redundant basis transforms in `B_trotter` multi-jump variant, and the `abs.(filter(...))` intermediate in the Time/Trotter thermalize `_jump_contribution!` (`jump_workers.jl`). All four are well-scoped, independently fixable, and affect functions that are called either once during precomputation (B variants) or once per thermalization step (jump_workers). The codebase already demonstrates the correct allocation-free patterns in adjacent code (e.g., `_precompute_R` in `trajectories.jl` uses the half-grid `continue` pattern instead of `filter`+`abs`).

The existing test suite has 224 tests across 7 test files, providing strong regression coverage. No `@allocated` tests exist yet -- this phase will introduce function-level allocation tests for each optimized hot path.

**Primary recommendation:** Fix each hotspot independently following the patterns already established in the codebase, then add `@allocated` tests to lock in the zero-allocation guarantees.

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

#### Sparse matrix strategy (coherent_bohr)
- Replace per-iteration `spzeros` + scatter with index-based accumulation: loop over `bohr_dict` indices directly, accumulate contributions element-wise into B without building any intermediate matrix
- The `f_A_nu_1` full-matrix broadcast (`f(bohr_freqs, nu_2) * jump.in_eigenbasis`) is acceptable -- element-wise broadcast into a pre-allocated dense matrix is cheap
- Apply the same index-based approach consistently to both single-jump and multi-jump `B_bohr` variants
- Do NOT precompute/cache per-frequency A_nu matrices -- Bohr frequency count explodes at larger system sizes (12 qubits for DM/trajectory)

#### Diagonal elimination (B_time / B_trotter)
- Replace `Diagonal(exp.(...))` and `Diagonal(eigvals .^ n)` wrappers with pre-allocated vector buffers + element-wise broadcasting (`.*`)
- No `Diagonal` objects should be created in any loop iteration
- Apply to all 4 variants: B_time single-jump, B_time multi-jump, B_trotter single-jump, B_trotter multi-jump

#### Redundant basis transforms (B_trotter)
- Both single-jump and multi-jump `B_trotter` variants currently compute `trotter.eigvecs' * jump.data * trotter.eigvecs` -- this is redundant because `JumpOp.in_eigenbasis` already stores the jump in Trotter eigenbasis when using TrotterDomain
- Fix: use `jump.in_eigenbasis` directly instead of recomputing the basis transform
- This applies to both single-jump (`jump_in_trotter = ...`) and multi-jump (`jump_a_trotter = ...` inside the inner loop)

#### Filter intermediate (jump_workers)
- The `abs.(filter(w -> w < 1e-12, energy_labels))` allocation may already be eliminated from earlier refactoring
- Researcher should verify whether this hotspot still exists in current code
- If still present: inline the condition in the loop body (skip + abs), keeping current Hermitian branch logic
- Note: `abs()` for energies may be obsolete since energy labels are now returned as a symmetrized grid around 0

#### Struct mutation scope
- Prefer function-local pre-allocated buffers over adding new struct fields
- Slightly open to struct changes ONLY if truly needed AND doesn't overcomplicate existing struct design
- No new structs
- Do NOT add domain-specific fields to cross-domain structs (structs have good structure across domains now)
- Keep in mind scale targets: DM/trajectory should run for 12 qubits, Lindbladian constructor for ~8 qubits

#### Verification
- Function-level `@allocated` tests for each optimized hot path: B_bohr, B_time, B_trotter, _jump_contribution! (if filter fix needed)
- Existing 224 tests continue to verify correctness
- No entry-point level allocation tests needed

### Claude's Discretion
- Exact implementation of index-based accumulation loop structure
- Choice of broadcasting vs manual loops for element-wise diagonal scaling
- How to structure @allocated test assertions (exact zero vs threshold)
- Whether f_A_nu_1 buffer needs pre-allocation or can rely on in-place broadcast

### Deferred Ideas (OUT OF SCOPE)
- **Peak memory estimation / allocation transparency** -- ability to see how many large objects (Lindbladians, DM-sized matrices) are simultaneously allocated during simulation. Should become its own phase for memory profiling/budgeting, enabling users to estimate RAM requirements for cluster simulations.
</user_constraints>

## Standard Stack

### Core
| Library | Purpose | Why Standard |
|---------|---------|--------------|
| Julia `Base` | `@allocated` macro for measuring allocations | Built-in, no dependency needed |
| `LinearAlgebra` | `mul!`, `rmul!`, broadcasting operators | Already used throughout codebase |
| `SparseArrays` | `spzeros`, `sparse` (to be eliminated) | Currently imported; elimination target |
| `Test` | `@testset`, `@test` for allocation tests | Already the test framework |

### Supporting
No additional libraries needed. All optimizations use existing Julia Base and LinearAlgebra primitives.

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| `@allocated` | `BenchmarkTools.@ballocated` | More statistically robust but adds dependency; `@allocated` is sufficient for zero-allocation assertions after warmup |
| Manual loops for diagonal scaling | `Diagonal` with pre-allocated diag vector | Still creates a wrapper struct; direct `.*` broadcasting is simpler and guaranteed zero-alloc |

## Architecture Patterns

### Relevant File Structure
```
src/
  bohr_domain.jl       # B_bohr single-jump and multi-jump (ALLOC-01)
  coherent.jl          # B_time (2 variants) and B_trotter (2 variants) (ALLOC-02, ALLOC-04)
  jump_workers.jl      # _jump_contribution! Time/Trotter thermalize (ALLOC-03)
  structs.jl           # JumpOp, workspaces (no changes needed)
test/
  test_helpers.jl      # Test system factories
  runtests.jl          # Test runner
```

### Pattern 1: Index-Based Accumulation for B_bohr (ALLOC-01)

**What:** Replace `spzeros` + index scatter + sparse transpose with direct index-based accumulation into pre-allocated dense B matrix.

**Current code** (`src/bohr_domain.jl` lines 12-19):
```julia
B = zeros(CT, dim, dim)
f_A_nu_1 = zeros(CT, dim, dim)
for nu_2 in unique_freqs
    A_nu_2 = spzeros(CT, dim, dim)                    # ALLOCATION per iteration
    indices = hamiltonian.bohr_dict[nu_2]
    A_nu_2[indices] .= jump.in_eigenbasis[indices]     # scatter
    @. f_A_nu_1 = f(hamiltonian.bohr_freqs, nu_2) * jump.in_eigenbasis
    B .+= A_nu_2' * f_A_nu_1                           # sparse*dense = allocation
end
```

**Optimized approach:** The operation `B += A_nu_2' * f_A_nu_1` where `A_nu_2` is sparse with nonzeros only at `indices` can be rewritten as:
```julia
# For each (i,j) in indices (where bohr_freqs[i,j] == nu_2):
#   A_nu_2[i,j] = jump.in_eigenbasis[i,j]
#   (A_nu_2')[j,i] = conj(jump.in_eigenbasis[i,j])
#   B[j,:] += conj(jump.in_eigenbasis[i,j]) * f_A_nu_1[i,:]
# This is a rank-1-like update per index pair.
```

The accumulation directly loops over `indices` from `bohr_dict[nu_2]`, extracting `(i,j)` pairs via `CartesianIndex`, and accumulates into B row-by-row. No sparse matrix is ever created.

**Key detail:** `bohr_dict` values are `Vector{CartesianIndex{2}}`. Each `CartesianIndex{2}` gives `(i,j)`. The operation `B .+= A_nu_2' * f_A_nu_1` when `A_nu_2` has only entries at those indices reduces to:
```julia
for idx in indices
    i, j = idx[1], idx[2]
    val = conj(jump.in_eigenbasis[i, j])
    @inbounds for col in 1:dim
        B[j, col] += val * f_A_nu_1[i, col]
    end
end
```

This pattern is already used in the Kraus BohrDomain `_jump_contribution!` in `jump_workers.jl` lines 253-301 (the `bohr_is`/`bohr_js` loop). The B_bohr code should follow the same pattern.

### Pattern 2: Vector Buffer Diagonal Replacement (ALLOC-02)

**What:** Replace `Diagonal(exp.(...))` closures with in-place element-wise operations on pre-allocated vectors.

**Current code** (`src/coherent.jl` line 164):
```julia
diag_time_evolve(t) = Diagonal(exp.(1im * hamiltonian.eigvals * t))
# ... later in loop:
U .= diag_time_evolve(s * beta)       # creates Diagonal, then broadcasts into matrix
```

**Why this allocates:** `exp.(1im * hamiltonian.eigvals * t)` allocates a new vector. `Diagonal(...)` wraps it (small allocation). Then `U .= Diagonal(...)` broadcasts the Diagonal into the dense matrix U, which works but the intermediate vector and Diagonal wrapper were allocated.

**Optimized approach:**
```julia
diag_buf = Vector{CT}(undef, dim)    # pre-allocate once outside loop
# ... in loop:
@. diag_buf = exp(1im * hamiltonian.eigvals * (s * beta))
# Then use diag_buf as a column-scaling vector in matrix operations
```

For the matrix products like `U * jump.in_eigenbasis' * U_minus_2 * jump.in_eigenbasis * U`, when `U` is diagonal, the product `U * M` is equivalent to `diag_buf .* M` (broadcasting each row), and `M * U` is equivalent to `M .* diag_buf'` (broadcasting each column). The full expression can be rewritten using element-wise scaling with the diagonal vectors.

**Critical detail for B_time/B_trotter:** The current code uses `zeros(CT, dim)` for `U` and `U_minus_2` -- but `dim = size(hamiltonian.data)` returns a tuple `(rows, cols)`, not an integer. This means `zeros(CT, dim)` actually creates a **matrix** of size `(rows, cols)`, not a vector. So `U .= diag_time_evolve(...)` broadcasts a Diagonal into a full dense matrix. The optimized version should work with vectors directly and use element-wise column/row scaling.

Actually, re-examining: `size(hamiltonian.data)` returns `(16, 16)` for 4 qubits, and `zeros(CT, (16,16))` creates a 16x16 matrix. So `U` is already a dense matrix receiving the broadcast of the Diagonal. The fix is to avoid creating the Diagonal entirely and instead store the diagonal values in a vector, then use element-wise operations.

### Pattern 3: Redundant Basis Transform Elimination (ALLOC-04)

**What:** `B_trotter` single-jump computes `jump_in_trotter = trotter.eigvecs' * jump.data * trotter.eigvecs` (line 223), and `B_trotter` multi-jump computes `jump_a_trotter = trotter.eigvecs' * jump_a.data * trotter.eigvecs` inside the inner loop (line 264). Both are redundant.

**Why redundant:** When using TrotterDomain, the `JumpOp.in_eigenbasis` field is already populated with the jump operator in the Trotter eigenbasis. This happens at:
- `simulations/main_thermalize.jl` lines 113-114: `basis_unitary = trotter.eigvecs; jump_op_in_eigenbasis = basis_unitary' * jump_op * basis_unitary`
- `simulations/main_liouv.jl` lines 114-115: same pattern
- `furnace.jl` lines 80-81: `transform_jumps_to_basis(jumps, trotter.eigvecs)` (creates new JumpOps with in_eigenbasis = trotter basis)

**Fix:** Replace `trotter.eigvecs' * jump.data * trotter.eigvecs` with `jump.in_eigenbasis` directly. This eliminates two matrix multiplications per call (single-jump) or two matrix multiplications per jump per b_plus iteration (multi-jump).

**IMPORTANT CAVEAT:** The test helpers (`test_helpers.jl`) construct JumpOps with `in_eigenbasis = hamiltonian.eigvecs' * jump_op * hamiltonian.eigvecs` (Hamiltonian eigenbasis), NOT Trotter eigenbasis. This means for TrotterDomain tests, the B_trotter function cannot simply use `jump.in_eigenbasis` unless the jump operators are first transformed to Trotter basis (which `construct_lindbladian` and `run_thermalization` already do via `transform_jumps_to_basis`). The B_trotter functions are called from `_precompute_coherent_total_B` and `_precompute_coherent_unitary_terms`, which pass the ORIGINAL jumps (not the Trotter-transformed ones).

Looking at the call chain:
1. `construct_lindbladian` calls `_precompute_coherent_total_B(jumps, ham_or_trott, ...)` with original jumps
2. `_precompute_coherent_total_B` calls `B_trotter(jumps, ham_or_trott, ...)` with original jumps
3. Then separately does `jumps_for_diss = transform_jumps_to_basis(jumps, trotter.eigvecs)` for dissipative part

So `B_trotter` currently receives jumps in Hamiltonian eigenbasis and manually transforms them. The fix requires either:
- (a) Transform jumps to Trotter basis BEFORE calling B_trotter, then use `jump.in_eigenbasis` inside, OR
- (b) Pass the Trotter-basis jumps to B_trotter

The cleanest approach: modify the callers (`_precompute_coherent_total_B`, `_precompute_coherent_unitary_terms`) to pass Trotter-basis jumps when domain is TrotterDomain. The `transform_jumps_to_basis` call already exists; it just needs to happen earlier or be shared.

### Pattern 4: Half-Grid Continue Pattern for filter (ALLOC-03)

**What:** Replace `abs.(filter(w -> w < 1e-12, energy_labels))` with inline loop condition.

**Current code** (`src/jump_workers.jl` line 415):
```julia
energies = jump.hermitian ? abs.(filter(w -> w < 1e-12, energy_labels)) : energy_labels
for w in energies
    ...
end
```

**This allocates** two vectors: `filter(...)` creates a new filtered vector, then `abs.(...)` creates another.

**Optimized approach** (already in `_precompute_R` at `trajectories.jl` lines 226-229):
```julia
for w_raw in energy_labels
    if jump.hermitian
        w_raw > 1e-12 && continue
        w = abs(w_raw)
    else
        w = w_raw
    end
    ...
end
```

Or more compactly:
```julia
if jump.hermitian
    for w_raw in energy_labels
        w_raw > 1e-12 && continue
        w = abs(w_raw)
        # ... body ...
        if w > 1e-12
            # ... negative frequency partner ...
        end
    end
else
    for w in energy_labels
        # ... body ...
    end
end
```

This matches the pattern already used in the EnergyDomain `_jump_contribution!` (lines 339-341, 370) and the `_precompute_R` function.

**Verification note:** The `abs()` call IS still needed even with symmetrized grid because the energy labels contain values like `-0.05, 0.0, 0.05, ...` and the `_prefactor_view` dict uses the original energy label values as keys. So for `w_raw = -0.05`, we need `w = abs(-0.05) = 0.05` which maps to the same NUFFT prefactor entry as the positive label `0.05`. The symmetrized grid guarantees both `+0.05` and `-0.05` exist in the energy_labels, and the NUFFT prefactor dict has entries for all energy labels, so `_prefactor_view(oft_nufft_prefactors, 0.05)` will succeed.

### Anti-Patterns to Avoid
- **Diagonal wrapper in hot loops:** `Diagonal(exp.(...))` allocates on every call. Use pre-allocated vectors and `.*` broadcasting.
- **filter() + abs.() instead of continue:** Creates two temporary vectors. Use inline loop conditions.
- **spzeros + scatter for sparse sub-matrix extraction:** Allocates a new sparse matrix per frequency. Use index-based accumulation.
- **Recomputing basis transforms available in struct fields:** Two matrix multiplications per call that duplicate existing data.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Diagonal matrix operations | Custom diagonal struct | Pre-allocated vectors + `.*` broadcasting | Julia's `.*` is already optimized; Diagonal wrapper adds allocation overhead |
| Sparse sub-matrix extraction | Custom sparse indexing | Direct index iteration over `bohr_dict` values | `CartesianIndex` iteration is allocation-free and cache-friendly |
| Allocation measurement | Custom timing/memory tracking | `@allocated` macro | Built-in, precise, composable with `@test` |

**Key insight:** The codebase already has the correct allocation-free patterns in newer code (trajectories.jl, kraus thermalize BohrDomain). The optimization is about bringing older code up to the same standard.

## Common Pitfalls

### Pitfall 1: @allocated Warmup Requirement
**What goes wrong:** `@allocated` returns non-zero on first call due to JIT compilation, method caching, and type inference.
**Why it happens:** Julia's JIT compiles code on first execution, which itself allocates.
**How to avoid:** Always call the function once before measuring with `@allocated`. Pattern:
```julia
# Warmup call (triggers compilation)
my_function(args...)
# Measurement call
allocs = @allocated my_function(args...)
@test allocs == 0
```
**Warning signs:** Tests that pass on re-run but fail on fresh `Pkg.test()` invocation.

### Pitfall 2: Broadcasting Into Pre-Existing Arrays vs Creating New Arrays
**What goes wrong:** `x = exp.(v)` allocates a new array, while `x .= exp.(v)` or `@. x = exp(v)` broadcasts into existing `x`.
**Why it happens:** Julia's dot-broadcasting syntax fuses when using `.=` or `@.`, but plain `=` with dotted calls still allocates the RHS.
**How to avoid:** Always use `.=` or `@.` macro when the target already exists. When introducing vector buffers for diagonal values, use `@. diag_buf = exp(1im * eigvals * t)`.
**Warning signs:** `@allocated` still showing allocations after seemingly in-place operations.

### Pitfall 3: Matrix Multiplication With Diagonal Vectors
**What goes wrong:** Replacing `Diagonal(v) * M` with `v .* M` seems straightforward but the broadcasting semantics differ. `Diagonal(v) * M` scales rows of M by v[i]. `v .* M` when v is a column vector scales the same way only if shapes align.
**Why it happens:** Julia broadcasting dimension alignment.
**How to avoid:** For row scaling (left multiply by diagonal): `v .* M` where v is a column vector. For column scaling (right multiply by diagonal): `M .* v'` (transpose v to a row vector).
**Warning signs:** Incorrect numerical results, dimension mismatch errors.

### Pitfall 4: Confusing Hamiltonian vs Trotter Eigenbasis
**What goes wrong:** Using `jump.in_eigenbasis` in B_trotter when it contains the Hamiltonian eigenbasis transform, not the Trotter eigenbasis transform.
**Why it happens:** The test helpers always store Hamiltonian eigenbasis in `in_eigenbasis`. The simulation entry points (main_*.jl) correctly set it to Trotter basis for TrotterDomain, but the test helpers do not.
**How to avoid:** The callers of B_trotter must ensure jumps are in Trotter basis before calling. Use `transform_jumps_to_basis(jumps, trotter.eigvecs)` before passing to B_trotter. This is what `construct_lindbladian` and `run_thermalization` already do for the dissipative part; the coherent precomputation must do the same.
**Warning signs:** B_trotter results change after optimization. Must maintain existing `trotter.eigvecs' * jump.data * trotter.eigvecs` behavior.

### Pitfall 5: The `dim = size(matrix)` Tuple vs Integer
**What goes wrong:** In `B_time` and `B_trotter`, `dim = size(hamiltonian.data)` or `dim = size(trotter.eigvecs)` returns a tuple `(rows, cols)`, not an integer. Then `zeros(CT, dim)` creates a matrix, not a vector.
**Why it happens:** `size()` without a dimension argument returns a tuple.
**How to avoid:** Use `size(matrix, 1)` to get an integer dimension. When creating vector buffers, use `Vector{CT}(undef, d)` where `d` is `size(matrix, 1)`.
**Warning signs:** Type errors when trying to use vector operations on what is actually a matrix.

## Code Examples

### Example 1: Index-Based B_bohr Accumulation (for single-jump variant)
```julia
function B_bohr(hamiltonian::HamHam{T}, jump::JumpOp, config) where {T}
    dim = size(hamiltonian.data, 1)
    CT = Complex{T}
    unique_freqs = keys(hamiltonian.bohr_dict)
    f = _pick_f(config)

    B = zeros(CT, dim, dim)
    f_A_nu_1 = zeros(CT, dim, dim)
    for nu_2 in unique_freqs
        @. f_A_nu_1 = f(hamiltonian.bohr_freqs, nu_2) * jump.in_eigenbasis

        indices = hamiltonian.bohr_dict[nu_2]
        # B += A_nu_2' * f_A_nu_1
        # A_nu_2 has nonzeros only at indices, so A_nu_2'[j,i] = conj(A[i,j])
        @inbounds for idx in indices
            i, j = idx[1], idx[2]
            val = conj(jump.in_eigenbasis[i, j])
            @inbounds for col in 1:dim
                B[j, col] += val * f_A_nu_1[i, col]
            end
        end
    end
    return B
end
```

### Example 2: Vector-Based Diagonal Replacement (B_time single-jump)
```julia
function B_time(jump::JumpOp, hamiltonian::HamHam, b_minus, b_plus, t0, beta, sigma)
    d = size(hamiltonian.data, 1)
    CT = Complex{eltype(hamiltonian.eigvals)}

    # Pre-allocated diagonal vector buffers
    diag_u = Vector{CT}(undef, d)
    diag_u2 = Vector{CT}(undef, d)

    # Pre-allocated matrix buffers for intermediate results
    b_plus_summand = zeros(CT, d, d)
    tmp1 = zeros(CT, d, d)
    tmp2 = zeros(CT, d, d)
    B = zeros(CT, d, d)

    eigvals = hamiltonian.eigvals
    jump_eig = jump.in_eigenbasis

    for s in keys(b_plus)
        t_s = s * beta
        @. diag_u = exp(1im * eigvals * t_s)
        @. diag_u2 = exp(-2im * eigvals * t_s)

        # U * jump' * U_minus_2 * jump * U
        # = diag(u) * jump' * diag(u2) * jump * diag(u)
        # Step 1: scale columns of jump by u -> tmp1 = jump * diag(u)
        # Step 2: scale rows by u2 -> diag(u2) * tmp1
        # etc.
        # ... (exact loop structure is Claude's discretion)

        b_plus_summand .+= b_plus[s] * (result)
    end

    for t in keys(b_minus)
        @. diag_u = exp(1im * eigvals * (t / sigma))
        # B += b_minus[t] * diag(u)' * b_plus_summand * diag(u)
        # = b_minus[t] * conj(diag_u) .* b_plus_summand .* diag_u'
        # (element-wise: B[i,j] += b_t * conj(u[i]) * bps[i,j] * u[j])
        # ... (exact loop structure is Claude's discretion)
    end

    return B .* t0^2
end
```

### Example 3: @allocated Test Pattern
```julia
@testset "Allocation: B_bohr" begin
    config = make_liouv_config(BohrDomain())
    jump = TEST_JUMPS[1]

    # Warmup (triggers JIT compilation)
    B = QuantumFurnace.B_bohr(TEST_HAM, jump, config)

    # Measure allocations on second call
    allocs = @allocated QuantumFurnace.B_bohr(TEST_HAM, jump, config)

    # B_bohr still allocates its return value (B matrix) and f_A_nu_1 buffer.
    # But it should NOT allocate sparse matrices per frequency.
    # Strategy: compare against known minimum (return matrix + one buffer).
    @test allocs <= 2 * sizeof(ComplexF64) * DIM^2 + 1024  # output + buffer + overhead
end
```

**Note on @allocated thresholds:** Functions that return newly-allocated matrices (like `B_bohr` returning `B`) will always show some allocation for the return value. The test should verify that per-iteration allocations are eliminated, not that the function is zero-alloc overall. Options:
1. Test that allocations scale as O(dim^2) not O(num_freqs * dim^2)
2. Extract the inner loop into a testable helper and verify that helper is zero-alloc
3. Use a threshold: `allocs <= C * dim^2` where C accounts for return value + pre-allocated buffers

### Example 4: Half-Grid Pattern (already in codebase)
From `src/trajectories.jl` lines 226-244 -- this is the reference pattern:
```julia
for w_raw in energy_labels
    w_raw > 1e-12 && continue
    w = abs(w_raw)

    nufft_prefactor_matrix = _prefactor_view(oft_nufft_prefactors, w)
    @. scratch.jump_oft = jump.in_eigenbasis * nufft_prefactor_matrix

    rate2_pos = base_prefactor * transition(w)
    mul!(scratch.LdagL, scratch.jump_oft', scratch.jump_oft)
    @. scratch.R += rate2_pos * scratch.LdagL

    if w > 1e-12
        rate2_neg = base_prefactor * transition(-w)
        mul!(scratch.LdagL, scratch.jump_oft, scratch.jump_oft')
        @. scratch.R += rate2_neg * scratch.LdagL
    end
end
```

## State of the Art

| Old Approach | Current Approach | Where | Impact |
|--------------|------------------|-------|--------|
| `spzeros` + scatter per Bohr freq | Index-based accumulation (done in Kraus path) | `jump_workers.jl` Kraus BohrDomain | Eliminates O(num_freqs) sparse allocations per jump |
| `filter` + `abs` vector creation | Half-grid `continue` pattern (done in `_precompute_R`) | `trajectories.jl` | Eliminates 2 vector allocations per step |
| `Diagonal(exp.(...))` per loop iter | Not yet optimized | `coherent.jl` B_time/B_trotter | Will eliminate vector + Diagonal wrapper per iter |

**Key insight:** The codebase has already solved 2 of 4 allocation patterns in newer code. Phase 11 propagates these solutions to the remaining older code paths.

## Specific Findings by Hotspot

### ALLOC-01: B_bohr Sparse Matrix (bohr_domain.jl)

**Location:** `src/bohr_domain.jl` lines 13-18 (single-jump) and 35-39 (multi-jump)

**Current allocations per inner loop iteration:**
1. `spzeros(CT, dim, dim)` -- allocates a sparse matrix
2. `A_nu_2[indices] .= ...` -- may trigger sparse structural modification
3. `A_nu_2' * f_A_nu_1` -- sparse-dense multiply creates a new dense matrix

**Number of iterations:** `length(keys(hamiltonian.bohr_dict))` = number of unique Bohr frequencies. For 4 qubits: ~16 unique frequencies. For 8 qubits: ~256. For 12 qubits: ~4096. This allocation scales badly.

**Fix complexity:** Moderate. The index-based loop pattern exists in `jump_workers.jl` (Kraus BohrDomain). Adapt it for the simpler B_bohr computation.

### ALLOC-02: Diagonal Wrappers (coherent.jl)

**Location:** `src/coherent.jl` lines 164, 191 (B_time) and 220, 252 (B_trotter)

**Current allocations per inner loop iteration:**
1. `Diagonal(exp.(...))` or `Diagonal(eigvals .^ n)` -- allocates vector + Diagonal wrapper
2. The `.=` broadcast into pre-allocated U then performs correctly, but the RHS intermediate is allocated

**Number of iterations:** `length(b_plus)` + `length(b_minus)` calls to the closure. With truncated b_plus/b_minus, typically 10-50 iterations each.

**Additional issue in B_time/B_trotter:** The `dim = size(...)` returns a tuple, and matrix operations like `U * jump.in_eigenbasis' * U_minus_2 * ...` create multiple intermediate matrices because `*` is left-associative and each binary `*` allocates a result. However, fixing ALL intermediate matrix allocations is a larger refactoring. The phase scope is specifically about eliminating the Diagonal wrapper allocations.

**Fix complexity:** Moderate. Replace the closure with in-place vector computation, then restructure the matrix products to use the diagonal vectors directly via element-wise operations. The product `diag(u) * M` for column vector u and matrix M is equivalent to `u .* M` (row-wise scaling), and `M * diag(u)` is `M .* u'` (column-wise scaling).

### ALLOC-03: filter+abs in jump_workers (jump_workers.jl)

**Location:** `src/jump_workers.jl` line 415

**Status:** CONFIRMED STILL PRESENT in current code.

**Current allocations:**
1. `filter(w -> w < 1e-12, energy_labels)` -- allocates a new vector
2. `abs.(...)` -- allocates another vector

These allocations happen once per `_jump_contribution!` call in the Time/Trotter thermalize path, which is called once per thermalization step.

**Fix complexity:** Trivial. Replace with the half-grid `continue` pattern already used in `_precompute_R` (trajectories.jl) and the Liouvillian Time/Trotter `_jump_contribution!` (jump_workers.jl lines 111-132). Exact same logic, just inline the condition.

### ALLOC-04: Redundant Basis Transform (coherent.jl)

**Location:**
- `src/coherent.jl` line 223: `jump_in_trotter = trotter.eigvecs' * jump.data * trotter.eigvecs` (single-jump)
- `src/coherent.jl` line 264: `jump_a_trotter = trotter.eigvecs' * jump_a.data * trotter.eigvecs` (multi-jump, inside inner loop)

**Current allocations:**
1. `trotter.eigvecs' * jump.data` -- allocates dim x dim matrix
2. `(...) * trotter.eigvecs` -- allocates another dim x dim matrix
3. In multi-jump: this happens per jump per b_plus iteration

**Approach:** The callers (`_precompute_coherent_total_B`, `_precompute_coherent_unitary_terms` in `coherent.jl`) should transform jumps to Trotter basis before calling B_trotter. Then B_trotter uses `jump.in_eigenbasis` directly.

**Key finding:** The callers pass ORIGINAL jumps (Hamiltonian eigenbasis) to B_trotter. The fix requires:
1. In `_precompute_coherent_total_B` and `_precompute_coherent_unitary_terms`: when domain is TrotterDomain, transform jumps via `transform_jumps_to_basis(jumps, trotter.eigvecs)` BEFORE calling B_trotter
2. In B_trotter: replace `trotter.eigvecs' * jump.data * trotter.eigvecs` with `jump.in_eigenbasis`

This is safe because `transform_jumps_to_basis` creates new JumpOps where `in_eigenbasis = trotter.eigvecs' * jump.data * trotter.eigvecs`, which is exactly what B_trotter currently computes.

**Fix complexity:** Low-moderate. Requires a small caller-side change (transform before call) and a one-line replacement inside B_trotter.

## Open Questions

1. **Matrix chain intermediate allocations in B_time/B_trotter**
   - What we know: The expression `U * jump.in_eigenbasis' * U_minus_2 * jump.in_eigenbasis * U` creates 3-4 intermediate matrices from left-to-right evaluation, even after eliminating the Diagonal allocation.
   - What's unclear: Whether eliminating these intermediates is in scope for this phase. The phase description says "eliminate Diagonal wrappers" specifically.
   - Recommendation: Focus on Diagonal elimination as specified. The intermediate matrix allocations from `*` chains are a separate optimization (could use `mul!` with scratch buffers) but are not in the phase requirements. The Diagonal fix alone still provides meaningful improvement.

2. **@allocated test threshold strategy**
   - What we know: Functions that return newly-allocated matrices will always show some allocation. `@allocated` includes the return value allocation.
   - What's unclear: Whether to test zero-alloc on extracted inner loops or use a dimensional threshold on the full function.
   - Recommendation: For B_bohr, B_time, B_trotter: test full function with a threshold of `C * dim^2 * sizeof(ComplexF64)` where C accounts for return value + pre-allocated buffers (2-4 matrices). For _jump_contribution! (Time/Trotter thermalize): test that the function is zero-alloc since it operates in-place via scratch buffers and the filter fix eliminates the only per-call allocation.

3. **B_trotter callers and jump basis transform timing**
   - What we know: `_precompute_coherent_total_B` and `_precompute_coherent_unitary_terms` need to transform jumps before calling B_trotter. `construct_lindbladian` already does `transform_jumps_to_basis` for dissipative but after coherent precomputation.
   - What's unclear: Whether to hoist the transform earlier (before coherent precomputation) and share it, or duplicate the transform call.
   - Recommendation: Transform jumps at the start of `_precompute_coherent_total_B`/`_precompute_coherent_unitary_terms` when domain is TrotterDomain. The cost of `transform_jumps_to_basis` is one-time (12 jumps * 2 matrix muls each = 24 muls), negligible compared to the per-iteration savings.

## Sources

### Primary (HIGH confidence)
- Direct codebase analysis of all source files in `src/` -- every finding verified by reading the actual code
- `src/bohr_domain.jl`: B_bohr implementation (lines 1-44)
- `src/coherent.jl`: B_time (lines 160-214), B_trotter (lines 216-277)
- `src/jump_workers.jl`: _jump_contribution! variants (lines 1-452)
- `src/trajectories.jl`: _precompute_R half-grid pattern (lines 209-260)
- `src/qi_tools.jl`: transform_jumps_to_basis (lines 18-25)
- `src/structs.jl`: JumpOp struct (lines 228-233)
- `src/furnace.jl`: construct_lindbladian and run_thermalization call chains
- `src/nufft.jl`: NUFFTPrefactors and _prefactor_view (lines 1-96)
- `test/test_helpers.jl`: test system factories and config factories

### Secondary (MEDIUM confidence)
- Julia documentation on `@allocated` macro behavior (warmup requirements, JIT interaction)
- Julia broadcasting semantics for `Diagonal` types

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- no external libraries needed; all Julia Base
- Architecture: HIGH -- all patterns verified by reading existing codebase code
- Pitfalls: HIGH -- based on direct code analysis and Julia language semantics
- Implementation details: HIGH -- four discrete, well-scoped changes with clear before/after patterns

**Research date:** 2026-02-15
**Valid until:** 2026-03-15 (stable -- Julia language features are mature)
