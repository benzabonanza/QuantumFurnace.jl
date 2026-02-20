---
phase: 34-extend-exact-diagnostics-to-support-trot
plan: 01
type: execute
wave: 1
depends_on: []
files_modified:
  - src/diagnostics.jl
  - test/test_diagnostics.jl
autonomous: true
must_haves:
  truths:
    - "run_exact_diagnostics works with a TrotterDomain Lindbladian and returns valid ExactDiagnosticsResult"
    - "DIAG-02 fixed point trace distance for TrotterDomain is finite (non-zero, since Trotter fixed point differs from Gibbs)"
    - "DIAG-05 overlap coefficients use Trotter-basis observables and Trotter-basis Gibbs correctly"
    - "DIAG-06 Sz labels use Trotter eigenvectors for basis transform when trotter object is provided"
    - "Existing BohrDomain diagnostics tests continue to pass unchanged"
  artifacts:
    - path: "src/diagnostics.jl"
      provides: "run_exact_diagnostics_from_config convenience wrapper, updated compute_sz_labels accepting eigvecs"
    - path: "test/test_diagnostics.jl"
      provides: "TrotterDomain diagnostics test coverage"
  key_links:
    - from: "src/diagnostics.jl"
      to: "src/furnace.jl"
      via: "construct_lindbladian call in convenience wrapper"
      pattern: "construct_lindbladian"
    - from: "test/test_diagnostics.jl"
      to: "src/diagnostics.jl"
      via: "run_exact_diagnostics with TrotterDomain L"
      pattern: "run_exact_diagnostics.*SMALL_TROTTER"
---

<objective>
Extend the exact diagnostics infrastructure (DIAG-01 through DIAG-06 and run_exact_diagnostics)
to support TrotterDomain Lindbladians, so spectral analysis can be performed on Trotterized
Lindbladian simulations -- not just BohrDomain ones.

Purpose: The diagnostics currently only work end-to-end with BohrDomain/Hamiltonian-eigenbasis
Lindbladians. TrotterDomain Lindbladians operate in the Trotter eigenbasis, requiring different
Gibbs state projection, observable basis transforms, and Sz label computation. This is needed
for Phase 27+ which will do fitting validation on TrotterDomain data.

Output: Updated src/diagnostics.jl with TrotterDomain support, comprehensive tests.
</objective>

<execution_context>
@./.claude/get-shit-done/workflows/execute-plan.md
@./.claude/get-shit-done/templates/summary.md
</execution_context>

<context>
@src/diagnostics.jl
@src/furnace.jl
@src/structs.jl
@src/trotter_domain.jl
@src/hamiltonian.jl
@test/test_diagnostics.jl
@test/test_helpers.jl
</context>

<tasks>

<task type="auto">
  <name>Task 1: Add TrotterDomain support to diagnostics core functions</name>
  <files>src/diagnostics.jl</files>
  <action>
Make three targeted changes to src/diagnostics.jl:

1. **Update `compute_sz_labels` signature** to accept an eigenvector matrix instead of requiring
   a HamHam object. Add a new method that accepts `eigvecs::Matrix{<:Complex}` alongside the
   existing `hamiltonian::HamHam` method (which becomes a thin wrapper calling the new method
   with `hamiltonian.eigvecs`). The new signature:

   ```julia
   function compute_sz_labels(eigen_result::EigenDecompositionResult, eigvecs::Matrix{<:Complex},
                               n_qubits::Int; n_modes::Int=20)
   ```

   This computes `Sz_eigen = eigvecs' * Sz_comp * eigvecs` using the provided eigvecs (which can
   be either `hamiltonian.eigvecs` or `trotter.eigvecs`). The existing `compute_sz_labels(eigen_result, hamiltonian; n_modes)` method should be preserved as a convenience wrapper that extracts `hamiltonian.eigvecs` and `Int(log2(size(hamiltonian.data, 1)))` and delegates.

2. **Add `run_exact_diagnostics_from_config` convenience function** that takes a `LiouvConfig`
   (or any `AbstractLiouvConfig`), a `HamHam`, and optionally a `TrottTrott`, and:
   - Constructs the dense Lindbladian via `construct_lindbladian(jumps, config, hamiltonian; trotter=trotter)`
   - For TrotterDomain: computes Gibbs in Trotter basis as `Hermitian(trotter.eigvecs' * hamiltonian.eigvecs * Matrix(hamiltonian.gibbs) * hamiltonian.eigvecs' * trotter.eigvecs)` (matching the pattern in `run_thermalization` at furnace.jl line 107)
   - For non-TrotterDomain: uses `hamiltonian.gibbs` directly
   - Selects `eigvecs_for_sz = trotter !== nothing && config.domain isa TrotterDomain ? trotter.eigvecs : hamiltonian.eigvecs`
   - Calls `run_exact_diagnostics(L_dense, hamiltonian, gibbs; eigvecs_for_sz=eigvecs_for_sz, ...)`

   However, since this requires jumps which are external to the diagnostics module, a simpler
   approach is better: **just add an `eigvecs_for_sz` keyword to `run_exact_diagnostics`** so
   callers can pass `trotter.eigvecs` when using TrotterDomain. The function already takes L
   as a pre-built matrix, so the caller already handles Lindbladian construction.

   Add this keyword argument to `run_exact_diagnostics`:
   ```julia
   eigvecs_for_sz::Union{Nothing, Matrix{<:Complex}}=nothing
   ```
   When `nothing` (default), uses `hamiltonian.eigvecs` (preserving backward compatibility).
   When provided, passes it to `compute_sz_labels` for Trotter-basis Sz labeling.

3. **Update the default observables and initial states in `run_exact_diagnostics`** to use
   `eigvecs_for_sz` when transforming computational-basis operators to the working basis.
   Specifically, add an `eigvecs` keyword that defaults to `nothing`, and when provided,
   use it instead of `hamiltonian.eigvecs` for:
   - Building default Z1 observable: `V' * Z1_comp * V` where `V = eigvecs` or `hamiltonian.eigvecs`
   - Building default initial states: `psi0_eigen = V' * psi0_comp` where V is the same
   - H diagonal in eigenbasis should still use `hamiltonian.eigvals` regardless (energy eigenvalues don't change)
   - BUT for TrotterDomain: H is NOT diagonal in the Trotter basis. So H_eigen should be `V' * hamiltonian.data * V` when V = trotter.eigvecs. Actually, simpler: just use `Matrix{ComplexF64}(V' * hamiltonian.data * V)` in all cases (it's diagonal when V = hamiltonian.eigvecs, non-diagonal when V = trotter.eigvecs -- both are correct).

   Actually, since the Lindbladian L is in the Trotter eigenbasis, ALL operators that interact
   with L must be in the Trotter eigenbasis. So the eigvecs_for_sz keyword should really be a
   general `basis_eigvecs` keyword that determines the working basis for everything:
   observables, initial states, AND Sz labels.

   **Final design:** Add a single keyword `basis_eigvecs::Union{Nothing, Matrix{<:Complex}}=nothing`
   to `run_exact_diagnostics`. When `nothing`, defaults to `hamiltonian.eigvecs`. This matrix
   is used for:
   - Default observable construction (Z1, H transform to working basis)
   - Default initial state construction (|0>^n, |+>^n transform to working basis)
   - Sz label computation (passed to compute_sz_labels)

   This keeps the API simple: one keyword controls the basis, backward compatible when omitted.
  </action>
  <verify>
  The module loads without errors:
  ```
  julia -e 'using Pkg; Pkg.activate("."); include("src/QuantumFurnace.jl")'
  ```
  Check that the existing `compute_sz_labels(eigen_result, hamiltonian; n_modes)` signature still works.
  </verify>
  <done>
  1. `compute_sz_labels` has a new method accepting `(eigen_result, eigvecs, n_qubits; n_modes)` and the old HamHam-based method delegates to it
  2. `run_exact_diagnostics` accepts `basis_eigvecs` keyword; when `nothing`, uses `hamiltonian.eigvecs`
  3. Default observables, initial states, and Sz labels all use the `basis_eigvecs` matrix
  4. Existing callers with no `basis_eigvecs` keyword see identical behavior
  </done>
</task>

<task type="auto">
  <name>Task 2: Add TrotterDomain diagnostics tests</name>
  <files>test/test_diagnostics.jl</files>
  <action>
Add a new `@testset "Diagnostics TrotterDomain"` block at the end of test/test_diagnostics.jl.
This testset validates that all six DIAG functions work correctly when the Lindbladian is
constructed for TrotterDomain.

The test should:

1. **Build TrotterDomain Lindbladian:**
   ```julia
   config_trott = make_small_liouv_config(TrotterDomain())
   L_trott_sparse = construct_lindbladian(SMALL_TROTTER_JUMPS, config_trott, SMALL_HAM; trotter=SMALL_TROTTER)
   L_trott = Matrix{ComplexF64}(L_trott_sparse)
   ```

2. **Compute Gibbs in Trotter basis:**
   ```julia
   gibbs_trott = Hermitian(SMALL_TROTTER.eigvecs' * SMALL_HAM.eigvecs * Matrix(SMALL_GIBBS) * SMALL_HAM.eigvecs' * SMALL_TROTTER.eigvecs)
   ```

3. **Test DIAG-01 on TrotterDomain L:**
   - `extract_leading_eigendata(L_trott; n_modes=10)` returns valid EigenDecompositionResult
   - Eigenvalues sorted by |Re|, first near zero, spectral gap > 0
   - Biorthonormality holds

4. **Test DIAG-02 with Trotter-basis Gibbs:**
   - `compute_fixed_point_distance(eigen_trott, gibbs_trott)` returns FixedPointResult
   - Trace distance should be non-trivially larger than BohrDomain (Trotter error shifts fixed point)
   - But still finite and reasonable (e.g., < 0.5 for 3-qubit system with 10 Trotter steps)
   - Fixed point is valid density matrix (Hermitian, trace 1, non-negative eigenvalues)

5. **Test DIAG-03/04 with Trotter-basis Gibbs:**
   - `compute_anti_hermitian_defect(L_trott, gibbs_trott)` returns DefectResult
   - All fields are valid (A_norm >= 0, H_gap > 0, ratio consistent)

6. **Test DIAG-05 with Trotter-basis observables:**
   - Build Z1 in Trotter basis: `Z1_trott = Matrix{ComplexF64}(SMALL_TROTTER.eigvecs' * Z1_comp * SMALL_TROTTER.eigvecs)`
   - Build rho0 in Trotter basis: maximally mixed `I(dim)/dim` (same in any basis)
   - `compute_overlap_coefficients(eigen_trott, [Z1_trott], ["Z1"], rho0, gibbs_trott; n_modes=10)`
   - Coefficients have correct dimensions, c_1 near zero

7. **Test DIAG-06 with Trotter eigenvectors:**
   - `compute_sz_labels(eigen_trott, SMALL_TROTTER.eigvecs, 3; n_modes=10)` returns valid labels
   - Each label has valid purity in [0, 1], non-empty sector_weights
   - Steady-state mode (k=1) should still have delta_sz ~ 0

8. **Test run_exact_diagnostics bundle with basis_eigvecs:**
   ```julia
   result = run_exact_diagnostics(L_trott, SMALL_HAM, gibbs_trott;
       n_modes=10, basis_eigvecs=SMALL_TROTTER.eigvecs)
   ```
   - Returns valid ExactDiagnosticsResult
   - All sub-results are the correct types
   - Default observables and initial states were built in Trotter basis (overlaps c_1 near zero)
   - Sz labels use Trotter eigenvectors

9. **Verify backward compatibility:** The existing BohrDomain test suite at the top of the file
   must continue to pass without modification (no `basis_eigvecs` argument = uses hamiltonian.eigvecs).
  </action>
  <verify>
  Run the full diagnostics test suite:
  ```
  cd /Users/bence/code/QuantumFurnace.jl && julia --project -e 'using Pkg; Pkg.test(test_args=["diagnostics"])'
  ```
  If test filtering is not available, run:
  ```
  cd /Users/bence/code/QuantumFurnace.jl && julia --project test/test_diagnostics.jl
  ```
  All existing BohrDomain tests pass, all new TrotterDomain tests pass.
  </verify>
  <done>
  1. TrotterDomain Lindbladian eigendecomposition returns valid results (sorted eigenvalues, biorthonormality)
  2. TrotterDomain fixed point has larger trace distance to Gibbs than BohrDomain (Trotter error) but is a valid density matrix
  3. TrotterDomain overlap coefficients computed in Trotter basis have c_1 near zero
  4. TrotterDomain Sz labels computed with trotter.eigvecs return valid sector assignments
  5. run_exact_diagnostics with basis_eigvecs=trotter.eigvecs produces complete ExactDiagnosticsResult
  6. All pre-existing BohrDomain diagnostics tests pass unchanged
  </done>
</task>

</tasks>

<verification>
- All existing tests pass: `julia --project -e 'using Pkg; Pkg.test()'`
- New TrotterDomain diagnostics tests pass
- No regressions in any test file
- The `basis_eigvecs` keyword is backward-compatible (omitting it gives identical results to before)
</verification>

<success_criteria>
- run_exact_diagnostics produces valid ExactDiagnosticsResult for TrotterDomain Lindbladians
- compute_sz_labels works with both HamHam and raw eigvecs arguments
- TrotterDomain fixed point distance is larger than BohrDomain (expected due to Trotter error) but finite
- All 6 DIAG functions produce valid outputs on TrotterDomain L
- Zero regressions in existing test suite
</success_criteria>

<output>
After completion, create `.planning/quick/34-extend-exact-diagnostics-to-support-trot/34-SUMMARY.md`
</output>
