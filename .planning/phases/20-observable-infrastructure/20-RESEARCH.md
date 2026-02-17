# Phase 20: Observable Infrastructure - Research

**Researched:** 2026-02-17
**Domain:** Quantum observable construction and basis transformation for spectral gap estimation
**Confidence:** HIGH

## Summary

Phase 20 adds two new observable builders to `convergence.jl`: `build_total_magnetization` (M_z = sum(Z_i)/n) and `build_gap_estimation_observables` (H + M_z bundle). Both must work in Hamiltonian eigenbasis and Trotter eigenbasis via a `trotter` keyword argument. The existing `build_convergence_observables` and `build_convergence_observables_trotter` functions provide a complete, working template for how observables are constructed, basis-transformed, and tested in this codebase.

This is a straightforward additive phase. No new dependencies are needed. The existing patterns in `convergence.jl` (lines 13-77) cover observable construction, eigenbasis transformation, and return format. The test patterns in `test_convergence.jl` (testsets 2-5) cover Hermiticity verification, Gibbs trace checks, and cross-basis validation. The main technical content is: (1) constructing M_z = sum(Z_i)/n in the computational basis, (2) transforming to the correct eigenbasis, and (3) bundling H + M_z into the gap estimation function.

**Primary recommendation:** Follow `build_convergence_observables` pattern exactly, unifying Hamiltonian and Trotter paths via a `trotter` keyword argument rather than separate function names.

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions
- Start with only H (energy) and M_z (total magnetization per site) for gap estimation
- No ZZ correlations in the gap estimation bundle -- keep it lean
- Rationale: choose observables that overlap with the first excited mode (slowest-decaying) to capture the spectral gap
- If H + M_z don't give a good gap estimate (cross-checked against exact Liouvillian), reconsider adding more observables later
- Per-site magnetization: M_z = sum(Z_i) / n (not total sum)
- Normalizes across system sizes, making amplitude comparison easier for n=4 vs n=6 cross-validation
- Decay rate (which determines spectral gap) is unaffected by normalization
- Both `build_total_magnetization` and `build_gap_estimation_observables` must have Trotter variants
- The Trotter spectral gap (algorithmic evolution) is the primary quantity of interest for the paper
- Pattern: `build_total_magnetization(ham, n; trotter=trotter)` and `build_gap_estimation_observables(ham, n; trotter=trotter)`
- Keep existing `(observables::Vector{Matrix{ComplexF64}}, names::Vector{String})` tuple pattern
- Consistent with `build_convergence_observables` -- downstream code already handles this

### Claude's Discretion
- Whether `build_gap_estimation_observables` internally calls `build_total_magnetization` or constructs directly
- Test system sizes and tolerance values for regression tests
- Internal helper organization

### Deferred Ideas (OUT OF SCOPE)
- ZZ correlations and other observables for gap estimation -- revisit if H + M_z insufficient (after Phase 24 cross-validation)
</user_constraints>

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| LinearAlgebra (stdlib) | Julia 1.11+ | Matrix operations, `diagm`, `tr`, `Hermitian`, `eigvals` | Already used throughout codebase |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| SparseArrays (stdlib) | Julia 1.11+ | `pad_term` returns sparse matrices | Already used via `pad_term` |

### Alternatives Considered
None -- this phase adds no new dependencies. All building blocks already exist in the codebase.

## Architecture Patterns

### Recommended Project Structure
```
src/
  convergence.jl    # ADD new functions here (same file as existing observable builders)
test/
  test_convergence.jl  # ADD new testsets here (same file as existing observable tests)
```

No new files. Both functions belong in `convergence.jl` alongside `build_convergence_observables`. Tests go into `test_convergence.jl` following the existing testset numbering.

### Pattern 1: Unified Trotter Keyword (User Decision)

**What:** Single function with optional `trotter` keyword instead of separate `_trotter` suffix functions.
**When to use:** For the new `build_total_magnetization` and `build_gap_estimation_observables` functions.
**Why:** The CONTEXT.md specifies the API pattern `build_total_magnetization(ham, n; trotter=trotter)`.
**Note:** This differs from the existing v1.2 pattern (`build_convergence_observables` vs `build_convergence_observables_trotter` as separate functions). The user's decision to use a keyword argument is locked.

**Example:**
```julia
function build_total_magnetization(hamiltonian::HamHam, num_qubits::Int;
                                    trotter::Union{TrottTrott, Nothing}=nothing)
    # Select basis transformation matrix
    V = trotter !== nothing ? trotter.eigvecs : hamiltonian.eigvecs

    # Build M_z in computational basis: sum(Z_i) / n
    Mz_comp = zeros(ComplexF64, 2^num_qubits, 2^num_qubits)
    for i in 1:num_qubits
        Mz_comp .+= Matrix{ComplexF64}(pad_term([Z], num_qubits, i))
    end
    Mz_comp ./= num_qubits  # Per-site normalization (locked decision)

    # Transform to eigenbasis: O_eigen = V' * O_comp * V
    Mz_eigen = Matrix{ComplexF64}(V' * Mz_comp * V)

    observables = [Mz_eigen]
    names = ["Mz"]

    return observables, names
end
```

### Pattern 2: Energy Observable in Eigenbasis vs Trotter Basis

**What:** H is diagonal in Hamiltonian eigenbasis (eigenvalues on diagonal) but NOT diagonal in Trotter eigenbasis (requires full basis transform).
**When to use:** For the H observable inside `build_gap_estimation_observables`.
**Critical detail:** This is already implemented in the existing code and must be replicated exactly.

```julia
# In Hamiltonian eigenbasis: H is diagonal
H_eigen = Matrix{ComplexF64}(diagm(ComplexF64.(hamiltonian.eigvals)))

# In Trotter eigenbasis: H is NOT diagonal, requires full transform
V_T = trotter.eigvecs
H_trotter = Matrix{ComplexF64}(V_T' * hamiltonian.data * V_T)
```

Source: `convergence.jl` lines 41-42 (eigenbasis) and lines 72-73 (Trotter basis).

### Pattern 3: Gap Estimation Bundle Calls Magnetization

**What:** `build_gap_estimation_observables` internally calls `build_total_magnetization` for M_z, then constructs H separately, and concatenates.
**Why:** Avoids duplicating M_z construction logic. If M_z definition changes later, it changes in one place.
**Recommendation (Claude's Discretion):** Use internal call.

```julia
function build_gap_estimation_observables(hamiltonian::HamHam, num_qubits::Int;
                                           trotter::Union{TrottTrott, Nothing}=nothing)
    # Get M_z from the dedicated builder
    mz_obs, mz_names = build_total_magnetization(hamiltonian, num_qubits; trotter=trotter)

    # Build H in the correct basis
    if trotter !== nothing
        V_T = trotter.eigvecs
        H = Matrix{ComplexF64}(V_T' * hamiltonian.data * V_T)
    else
        H = Matrix{ComplexF64}(diagm(ComplexF64.(hamiltonian.eigvals)))
    end

    observables = [H; mz_obs]
    names = ["H"; mz_names]

    return observables, names
end
```

### Pattern 4: Return Format (Locked Decision)

**What:** Return `(observables::Vector{Matrix{ComplexF64}}, names::Vector{String})` tuple.
**Source:** `convergence.jl` lines 19, 46, 55, 77 -- all existing builders use this exact pattern.
**Why locked:** Downstream code (`run_trajectories_convergence`, `run_trajectories_adaptive`) already accepts this format.

### Anti-Patterns to Avoid
- **Separate `_trotter` suffix functions for the new API:** The user locked the keyword argument pattern. Do NOT create `build_total_magnetization_trotter` as a separate function.
- **Including ZZ correlations in gap estimation:** Explicitly deferred. The gap bundle is H + M_z ONLY.
- **Using total sum instead of per-site normalization for M_z:** The user locked M_z = sum(Z_i) / n, not sum(Z_i).
- **Constructing M_z from `hamiltonian.data`:** M_z is sum(Z_i)/n, which is a Pauli sum -- it has nothing to do with the Hamiltonian matrix. Build it from `pad_term([Z], ...)`.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Pauli tensor products | Custom kron logic | `pad_term([Z], num_qubits, site)` | Already handles periodic boundaries, sparse ops, site padding |
| Basis transformation | Manual index manipulation | `V' * O_comp * V` | Standard unitary transformation, used everywhere in codebase |
| Gibbs state in Trotter basis | Re-derive transformation | `_gibbs_in_trotter_basis(ham, trotter)` | Already exists in `convergence.jl` line 89 |
| Observable Gibbs values | Manual trace computation | `_compute_gibbs_observable_values(gibbs, observables)` | Already exists in `convergence.jl` line 99 |

**Key insight:** Every building block needed already exists. The only new code is the M_z construction (summing Z_i and dividing by n) and the thin bundle function.

## Common Pitfalls

### Pitfall 1: Wrong Basis for Observables
**What goes wrong:** Observable built in computational basis is used directly with `run_trajectories` output, which returns rho in eigenbasis. `tr(rho * O)` gives wrong expectation values.
**Why it happens:** Forgetting the basis transform step `V' * O_comp * V`.
**How to avoid:** Always apply the basis transform. The `trotter !== nothing` branch selects `trotter.eigvecs`, the default branch selects `hamiltonian.eigvecs`.
**Warning signs:** `tr(gibbs * Mz)` does not match sum of individual `<Z_i>_gibbs`.
**Source:** This exact bug class is documented in the docstring of `build_convergence_observables` (convergence.jl line 23: "CRITICAL: All observables are in the Hamiltonian eigenbasis").

### Pitfall 2: H Diagonal Only in Hamiltonian Eigenbasis
**What goes wrong:** Using `diagm(hamiltonian.eigvals)` for H in the Trotter basis.
**Why it happens:** H is diagonal in its own eigenbasis, but NOT in the Trotter eigenbasis.
**How to avoid:** For Trotter: `V_T' * hamiltonian.data * V_T`. For eigenbasis: `diagm(eigvals)`.
**Warning signs:** `tr(gibbs_trotter * H_trotter)` gives the wrong energy.
**Source:** Already correctly handled in `build_convergence_observables_trotter` (convergence.jl line 72).

### Pitfall 3: Forgetting Per-Site Normalization
**What goes wrong:** M_z = sum(Z_i) instead of sum(Z_i) / n.
**Why it happens:** Physicist instinct to use total magnetization.
**How to avoid:** Divide by `num_qubits` after the sum. This is a locked decision.
**Warning signs:** Amplitude of M_z observable scales with system size during cross-validation.

### Pitfall 4: `pad_term` Returns Sparse Matrix
**What goes wrong:** Type mismatch when accumulating into dense matrix.
**Why it happens:** `pad_term` returns `SparseMatrixCSC{ComplexF64}`.
**How to avoid:** Wrap with `Matrix{ComplexF64}(pad_term(...))` or accumulate with `.+=` into pre-allocated dense matrix.
**Warning signs:** Compilation errors or unexpected sparse matrix types in output.
**Source:** `pad_term` in `misc_tools.jl` line 326: return type is `SparseMatrixCSC{ComplexF64}`.

### Pitfall 5: Export Statement
**What goes wrong:** New public functions are not exported, cannot be called by users.
**Why it happens:** Forgetting to add to the `export` block in `QuantumFurnace.jl`.
**How to avoid:** Add `build_total_magnetization` and `build_gap_estimation_observables` to the convergence tracking export block (line 47 of `QuantumFurnace.jl`).

## Code Examples

Verified patterns from the existing codebase:

### Constructing a Single-Site Observable in Eigenbasis
```julia
# Source: convergence.jl lines 31-37
V = hamiltonian.eigvecs
ZZ_comp = Matrix{ComplexF64}(pad_term([Z, Z], num_qubits, i; periodic=true))
ZZ_eigen = V' * ZZ_comp * V
```

### Constructing the Energy Observable
```julia
# Source: convergence.jl lines 41-42 (eigenbasis, diagonal)
H_eigen = Matrix{ComplexF64}(diagm(ComplexF64.(hamiltonian.eigvals)))

# Source: convergence.jl lines 72-73 (Trotter basis, NOT diagonal)
H_trotter = V_T' * hamiltonian.data * V_T
```

### Gibbs Trace Verification (Test Pattern)
```julia
# Source: test_convergence.jl lines 81-84
boltz = exp.(-BETA .* TEST_HAM.eigvals)
gibbs_energy_analytical = sum(TEST_HAM.eigvals .* boltz) / sum(boltz)
gibbs_energy_from_obs = real(tr(Matrix(TEST_GIBBS) * H_eigen))
@test isapprox(gibbs_energy_from_obs, gibbs_energy_analytical; atol=1e-10)
```

### Analytical Gibbs Magnetization for Test Verification
```julia
# Compute sum_i <Z_i>_gibbs / n analytically (computational basis)
# 1. Build Z_i in computational basis
# 2. Compute gibbs_comp = V * gibbs_eigen * V' (transform back)
# 3. Sum tr(gibbs_comp * Z_i) for each site, divide by n
gibbs_comp = hamiltonian.eigvecs * Matrix(hamiltonian.gibbs) * hamiltonian.eigvecs'
mz_analytical = sum(
    real(tr(gibbs_comp * Matrix{ComplexF64}(pad_term([Z], num_qubits, i))))
    for i in 1:num_qubits
) / num_qubits
```

### Trotter Gibbs Trace Verification (Test Pattern)
```julia
# Source: test_convergence.jl lines 113-125, convergence.jl lines 89-91
gibbs_trotter = QuantumFurnace._gibbs_in_trotter_basis(TEST_HAM, TEST_TROTTER)
# Verify: trace 1, Hermitian, positive eigenvalues
@test isapprox(real(tr(gibbs_trotter)), 1.0; atol=1e-12)
# Then: tr(gibbs_trotter * Mz_trotter) should equal mz_analytical
```

## Verification Approach for M_z

The key verification is: `tr(gibbs * M_z) = sum_i <Z_i>_gibbs / n`.

**In Hamiltonian eigenbasis:**
1. Build M_z in computational basis: `Mz_comp = sum(pad_term([Z], n, i) for i in 1:n) / n`
2. Transform: `Mz_eigen = V' * Mz_comp * V`
3. Verify: `real(tr(Matrix(gibbs) * Mz_eigen))` equals analytical value

**Analytical value computation:**
- `gibbs_comp = V * gibbs_eigen * V'` (Gibbs in computational basis)
- `analytical_mz = sum(real(tr(gibbs_comp * pad_term([Z], n, i))) for i in 1:n) / n`
- These must agree to machine precision (~1e-12)

**In Trotter eigenbasis:**
- Same logic but with `V_T = trotter.eigvecs` and `gibbs_trotter = _gibbs_in_trotter_basis(ham, trotter)`
- The analytical value should be the same (it's a physical observable, basis-independent)

## Recommended Test Sizes and Tolerances (Claude's Discretion)

| System | Qubits | Dim | Purpose |
|--------|--------|-----|---------|
| TEST_HAM (4-qubit) | 4 | 16 | Primary regression (matches existing convergence tests) |
| SMALL_HAM (3-qubit) | 3 | 8 | Fast sanity check |

| Test Type | Tolerance | Rationale |
|-----------|-----------|-----------|
| Hermiticity check | `atol=1e-12` | Machine precision, matching existing testset 2 |
| Gibbs trace check | `atol=1e-10` | Matching existing testset 2 line 84 |
| Trotter trace check | `atol=1e-10` | Same as eigenbasis check |
| Cross-basis consistency | `atol=1e-10` | Physical observable must agree across bases |

## Recommended Internal Organization (Claude's Discretion)

**Recommendation:** `build_gap_estimation_observables` should internally call `build_total_magnetization`.

**Rationale:**
1. DRY: M_z logic in one place
2. `build_total_magnetization` is independently useful (exported separately)
3. If M_z definition changes, it propagates automatically
4. The bundle function becomes a thin combinator: get M_z, get H, concatenate

## State of the Art

| Old Approach (v1.2) | Current Approach (Phase 20) | Impact |
|---------------------|----------------------------|--------|
| Separate `_trotter` suffix functions | Single function with `trotter` keyword | Cleaner API, one entry point |
| ZZ + H for convergence | H + M_z for gap estimation | Different observable set tuned for spectral gap capture |
| No per-site normalization | M_z = sum(Z_i) / n | System-size-independent amplitude |

## Open Questions

1. **Should existing `build_convergence_observables` be retrofitted to use the `trotter` keyword pattern?**
   - What we know: The existing pair works fine with separate functions. The new functions use keyword pattern (locked decision).
   - What's unclear: Whether to unify the old API to match the new one.
   - Recommendation: Do NOT change existing API in this phase. That would be a refactoring concern for a future phase. Keep backward compatibility.

2. **Note on OBS-03 vs CONTEXT.md tension**
   - The original OBS-03 requirement says "H, M_z, all ZZ correlations" but CONTEXT.md locks "No ZZ correlations in the gap estimation bundle."
   - Resolution: CONTEXT.md (user decisions) overrides the original requirement text. The bundle is H + M_z only.
   - The requirement OBS-03 should be considered satisfied by H + M_z, per user decision.

## Sources

### Primary (HIGH confidence)
- `src/convergence.jl` -- existing `build_convergence_observables` and `build_convergence_observables_trotter` functions (lines 13-77)
- `src/convergence.jl` -- `_gibbs_in_trotter_basis` and `_compute_gibbs_observable_values` helpers (lines 89-101)
- `test/test_convergence.jl` -- testsets 2-5 covering observable construction and Gibbs verification
- `src/hamiltonian.jl` -- `HamHam` struct with `eigvals`, `eigvecs`, `gibbs`, `data` fields
- `src/trotter_domain.jl` -- `TrottTrott` struct with `eigvecs` field
- `src/constants.jl` -- Pauli `Z` matrix definition
- `src/misc_tools.jl` -- `pad_term` function for tensor product padding
- `src/QuantumFurnace.jl` -- export list (line 47) for convergence tracking functions
- `test/test_helpers.jl` -- shared test fixtures (`TEST_HAM`, `TEST_TROTTER`, `TEST_GIBBS`, `NUM_QUBITS`, `DIM`, `BETA`)

### Secondary (MEDIUM confidence)
- None needed -- all information from codebase inspection

### Tertiary (LOW confidence)
- None

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- no new dependencies, all building blocks exist in codebase
- Architecture: HIGH -- follows established patterns from v1.2 convergence tracking
- Pitfalls: HIGH -- all pitfalls derived from existing code patterns and documented gotchas in the codebase
- Code examples: HIGH -- all examples copied from verified, working source code

**Research date:** 2026-02-17
**Valid until:** No expiry -- codebase-internal research, no external dependency version concerns
