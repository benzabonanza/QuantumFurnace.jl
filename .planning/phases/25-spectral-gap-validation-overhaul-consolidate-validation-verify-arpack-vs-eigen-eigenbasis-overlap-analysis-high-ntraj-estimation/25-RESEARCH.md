# Phase 25: Spectral Gap Validation Overhaul - Research

**Researched:** 2026-02-18
**Domain:** Lindbladian spectral gap estimation cleanup, ARPACK verification, eigenbasis overlap analysis
**Confidence:** HIGH

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

#### Consolidation -- fresh start
- Delete ALL existing validation scripts (experiments/validate_*.jl, experiments/run_sweep.jl, etc.)
- Delete cross_validate_gap function and CrossValidationResult struct
- Delete overlapping/redundant observable builder functions -- audit what exists, keep only what's needed
- Read quick-fix summaries (quick-22, quick-23, quick-24) for context on what went wrong, but treat them as potentially containing errors and contradictions -- that's why this phase exists
- Goal: minimal function count, logically organized

#### Observable set
- Single function: `build_preset_trajectory_observables()` (or similar name)
- Fixed set of 5 observables: H (energy), M_z (total magnetization per site), XX_avg, YY_avg, ZZ_avg (per-bond averaged 2-site correlations)
- Per-bond averaged versions, NOT individual bond pairs
- Must be correctly transformed to the simulation basis

#### Eigenbasis overlap analysis
- Separate exported diagnostic function (not embedded in validation script)
- Decomposes each observable into the Lindbladian's eigenbasis
- Reports overlap of each observable with the slowest decaying mode (first excited eigenmode of L)
- Larger overlap = better observable for gap estimation
- This is the key diagnostic: if gap estimation fails, overlap analysis explains whether the chosen observables can see the gap at all

#### ARPACK vs eigen verification
- Check that ARPACK eigs method in run_lindbladian delivers the same spectral gap as eigen() for n=4
- This is a sanity check -- if they disagree, something is fundamentally wrong
- If they agree, proceed with confidence that the exact gap reference is correct

#### Estimation protocol
- System sizes: n=4 and n=6 (Heisenberg chain)
- Trajectories: 20,000 (high enough to rule out statistical noise)
- Parameters: beta=10, delta=0.01
- Target accuracy: relative error < 1e-2 (1%) between fitted and exact gap
- If target not met: must explain WHY with evidence, suggest followup tests
  - Possible reason: insufficient overlap with first excited mode
  - Possible reason: discrete-step Kraus effect (known from quick-22)
  - Must provide concrete diagnostic output, not just "it didn't work"

#### Validation script
- One script that runs everything: exact gap computation, trajectory estimation, comparison
- Prints clear pass/fail with gap values, relative error, per-observable fit quality
- Calls the eigenbasis overlap diagnostic to show overlap table
- Lives in experiments/ directory

### Claude's Discretion
- Exact function signatures and names (within the spirit of minimal API)
- How to structure the overlap computation internally
- Whether to use KMS or GNS detailed balance (whatever the codebase currently defaults to)
- Test structure (integration test vs script-only validation)

### Deferred Ideas (OUT OF SCOPE)
None -- discussion stayed within phase scope
</user_constraints>

## Summary

This phase is a clean-slate overhaul of the spectral gap validation infrastructure built across phases 20-24 and quick-fixes 22-24. The existing code works but is bloated: there are 4 experiment scripts, 4 observable builder functions with overlapping responsibility, a `CrossValidationResult` struct that duplicates simple arithmetic, and scattered validation logic across files. The user suspects hallucination may have crept in during rapid quick-fix iteration -- and indeed the code shows evidence of layered fixes (quick-22 fixed a fundamental delta double-counting bug, quick-23 changed selection logic, quick-24 added observables) rather than a coherent design.

The overhaul deletes all validation scripts and the `CrossValidationResult`/`cross_validate_gap` API surface. It replaces the 4 observable builder functions (`build_convergence_observables`, `build_convergence_observables_trotter`, `build_total_magnetization`, `build_gap_estimation_observables`) with a single `build_preset_trajectory_observables()`. It adds a new eigenbasis overlap diagnostic function (an exported analysis tool, not just a script). It verifies ARPACK vs eigen for the exact gap reference, then runs 20k-trajectory estimation for n=4 and n=6 at beta=10, targeting <1% relative error.

**Primary recommendation:** Execute deletion first, rebuild second. Delete all 4 experiment scripts and `CrossValidationResult`/`cross_validate_gap` from `gap_estimation.jl`. Consolidate observable builders into one function. Add eigenbasis overlap analysis as a new exported function. Write one unified validation script. The core machinery (`estimate_spectral_gap`, `fit_exponential_decay`, `run_observable_trajectories`, `_select_best_observable`) is sound and stays.

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| LsqFit.jl | (existing dep) | Levenberg-Marquardt exponential fitting | Already locked in phase 21, no change |
| Arpack.jl | (existing dep) | Sparse eigenvalue computation (eigs) | Used by run_lindbladian for spectral gap |
| LinearAlgebra | stdlib | Dense eigen(), matrix operations | Julia stdlib, no dep management |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| Printf | stdlib | Formatted output in validation scripts | Script output only |

### Alternatives Considered
None -- this phase uses only existing dependencies. No new packages needed.

## Architecture Patterns

### Current File Layout (relevant files)
```
src/
  convergence.jl      # Observable builders + convergence runners
  fitting.jl          # fit_exponential_decay, FitResult
  gap_estimation.jl   # estimate_spectral_gap, SpectralGapResult, CrossValidationResult, cross_validate_gap
  furnace.jl          # run_lindbladian (ARPACK eigs), construct_lindbladian
  QuantumFurnace.jl   # Module file with exports

experiments/
  validate_gap_estimation.jl   # n=4 + n=6 cross-validation (TimeDomain)
  run_gap_validation.jl        # n=6 only, 10k traj (TrotterDomain)
  run_sweep.jl                 # 27-experiment KMS-vs-GNS sweep (unrelated to gap)
  eigenmode_decomposition.jl   # Eigenmode analysis diagnostic (TrotterDomain, n=4)

test/
  test_gap_estimation.jl       # Tests for estimate_spectral_gap + cross_validate_gap
  test_convergence.jl          # Tests for observable builders + convergence runners
```

### Target File Layout (after Phase 25)
```
src/
  convergence.jl      # ONLY: run_trajectories_convergence, run_trajectories_adaptive, helpers
  fitting.jl           # UNCHANGED: fit_exponential_decay, FitResult
  gap_estimation.jl    # estimate_spectral_gap, SpectralGapResult, _select_best_observable
                       # + NEW: eigenbasis_overlap_analysis (the new diagnostic function)
                       # - DELETED: CrossValidationResult, cross_validate_gap
  furnace.jl           # UNCHANGED: run_lindbladian, construct_lindbladian
  QuantumFurnace.jl    # Updated exports

experiments/
  validate_spectral_gap.jl    # NEW: single unified validation script
  - DELETED: validate_gap_estimation.jl, run_gap_validation.jl, run_sweep.jl, eigenmode_decomposition.jl

test/
  test_gap_estimation.jl       # Updated: remove cross_validate_gap tests, add overlap analysis tests
  test_convergence.jl          # Updated: remove old builder tests, update for renamed function
```

### Pattern 1: Observable Builder Consolidation

**What:** Replace 4 overlapping observable builders with 1 unified function.

**Current state (4 functions):**
1. `build_convergence_observables(ham, n)` -- returns per-bond ZZ_ij + H in eigenbasis
2. `build_convergence_observables_trotter(ham, trotter, n)` -- same in Trotter basis
3. `build_total_magnetization(ham, n; trotter=nothing)` -- returns Mz in either basis
4. `build_gap_estimation_observables(ham, n; trotter=nothing)` -- returns H + Mz + XX_avg + YY_avg + ZZ_avg in either basis

**Target state (1 function):**
- `build_preset_trajectory_observables(ham, n; trotter=nothing)` -- returns the same 5 observables as current `build_gap_estimation_observables`, with the `trotter` keyword controlling basis

**Key insight:** `build_convergence_observables` and `build_convergence_observables_trotter` are used ONLY by:
- `experiments/run_sweep.jl` (being deleted)
- `test/test_convergence.jl` (tests for convergence runners use them as fixture inputs)

The convergence runner tests need observables as fixture input. They can use the new `build_preset_trajectory_observables` or a simple inline construction. The convergence runners themselves (`run_trajectories_convergence`, `run_trajectories_adaptive`) take observables as parameters and never call any builder internally.

**Risk:** `build_convergence_observables` is exported in `QuantumFurnace.jl` line 48. Removing it is a public API break. But since this is pre-release code (not published on General registry), this is acceptable per the user's directive.

### Pattern 2: ARPACK vs eigen Verification

**What:** The existing `run_lindbladian` uses ARPACK `eigs` with a shift-invert strategy to find the 2 eigenvalues nearest to zero. The validation scripts use `eigen()` (dense) for exact reference. This phase verifies they agree.

**Current ARPACK call in `furnace.jl` (lines 16-17):**
```julia
shift = 1e-9 * (1 + 1im)
eigvals_near_zero, eigvecs_near_zero = eigs(liouv, nev=2, sigma=shift, tol=1e-12)
```

**Current eigen() approach in `validate_gap_estimation.jl` (lines 111-135):**
```julia
eig = eigen(L)
sorted_idx = sortperm(abs.(real.(eig.values)))
spectral_gap = eig.values[sorted_idx[2]]
```

**Verification approach:** For n=4 (dim=16, Liouvillian is 256x256), both methods are tractable. Construct L, run both, compare `abs(real(spectral_gap))`. They should agree to ~1e-10 or better.

**Confidence: HIGH** -- this is straightforward linear algebra verification. Both methods are well-understood. The only subtlety is that ARPACK uses shift-invert mode with a small complex shift to avoid the exact singularity at eigenvalue 0.

### Pattern 3: Eigenbasis Overlap Analysis

**What:** A diagnostic function that decomposes each observable into the Lindbladian eigenbasis and measures how strongly each observable couples to the first excited mode (the gap mode).

**Mathematical formulation:**

Given Liouvillian L with eigendecomposition L = V diag(lambda) V^{-1}:
- Let rho_0 = |psi_0><psi_0| be the initial density matrix
- Let alpha = V^{-1} vec(rho_0) be the initial state expansion in eigenmodes
- For observable O, the time-dependent expectation value is:
  <O>(t) = sum_k c_k exp(lambda_k * t)
  where c_k = [vec(O)^dagger * v_k] * alpha_k

The "overlap" of observable O with mode k is |c_k|. The gap mode is k=1 (second eigenvalue after steady state). An observable with zero overlap with the gap mode cannot see the spectral gap, regardless of trajectory count.

**The key diagnostic metric:** For each observable, report:
1. |c_1| (absolute overlap with gap mode)
2. |c_1| / sum_k |c_k| (relative overlap -- fraction of signal in gap mode)
3. The effective initial decay rate: -f'(0)/f(0) compared to the exact gap

This function already exists as a script (`experiments/eigenmode_decomposition.jl`) but needs to be converted into an exported library function.

**Recommended signature:**
```julia
function eigenbasis_overlap_analysis(
    L::Matrix{<:Complex},       # Full Liouvillian matrix
    observables::Vector{<:Matrix{<:Complex}},
    observable_names::Vector{String},
    rho0::Matrix{<:Complex};    # Initial density matrix
) -> NamedTuple or struct
```

**Discretion area: struct vs NamedTuple for return type.** Recommendation: Use a lightweight NamedTuple or a new struct `EigenbasisOverlapResult`. A struct is better for Aqua compliance and future-proofing.

### Pattern 4: Validation Script Structure

**What:** One script that runs:
1. ARPACK vs eigen verification (n=4)
2. Exact gap computation (n=4, n=6 via eigen)
3. Eigenbasis overlap analysis (n=4, n=6)
4. Trajectory estimation with 20k trajectories (n=4, n=6)
5. Comparison with pass/fail judgment

**Domain choice (Claude's discretion):** The existing validation scripts use two different domains:
- `validate_gap_estimation.jl` uses **TimeDomain** with `with_coherent=false`
- `run_gap_validation.jl` uses **TrotterDomain** with `with_coherent=true`

**Recommendation: Use TimeDomain with with_coherent=false** for the canonical validation. Reasons:
1. TimeDomain is simpler (no Trotter error on top of everything else)
2. `with_coherent=false` removes the coherent term complexity
3. The eigenmode decomposition script used TrotterDomain, but that adds Trotter approximation error which confuses the gap estimation analysis
4. If TimeDomain validation passes, TrotterDomain can be tested separately

**However:** The test infrastructure uses a disordered Heisenberg chain (loaded from BSON), while validation scripts create a clean periodic Heisenberg chain via `HamHam([[X,X],[Y,Y],[Z,Z]], [1,1,1], n, beta; periodic=true)`. The validation script should use the same clean Heisenberg chain pattern (not the test BSON fixtures), since the analysis is about gap estimation methodology, not a specific disorder realization.

### Anti-Patterns to Avoid

- **Don't recreate CrossValidationResult:** The relative error `|fitted - exact| / exact` is a single line of arithmetic. Wrapping it in a struct with 7 fields is over-engineering.
- **Don't mix domains in one script:** Use a single domain (TimeDomain) for clean results.
- **Don't use TrotterDomain for validation unless necessary:** Trotter error adds a confounding variable. Save it for a separate follow-up.
- **Don't keep dead observable builders:** If `build_convergence_observables` is only needed by deleted scripts and by test fixtures that can use the new function, remove it.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Exponential fitting | Custom optimizer | `fit_exponential_decay` (existing) | Already handles LsqFit, CI, bounds |
| Observable selection | New selection logic | `_select_best_observable` (existing) | Already has smallest-gap criterion |
| Trajectory simulation | Custom loop | `estimate_spectral_gap` (existing) | Composes trajectory runner + fitting |
| Eigendecomposition | Custom eigensolver | `eigen()` from LinearAlgebra | Exact for small systems |

**Key insight:** The core estimation pipeline (`build_observables` -> `run_observable_trajectories` -> `fit_exponential_decay` -> `_select_best_observable`) is sound. This phase only replaces the validation/diagnostic scaffolding around it, not the pipeline itself.

## Common Pitfalls

### Pitfall 1: Removing exports that are used in tests
**What goes wrong:** Deleting `build_convergence_observables` from exports but forgetting to update `test_convergence.jl` which uses it extensively as a fixture input.
**Why it happens:** `test_convergence.jl` uses `build_convergence_observables` in 8 testsets as the observable fixture for `run_trajectories_convergence` and `run_trajectories_adaptive`.
**How to avoid:** When removing `build_convergence_observables`, update all convergence tests to use the new `build_preset_trajectory_observables` instead. The convergence runners take observables as parameters -- they don't care about the builder function.
**Warning signs:** Tests fail with `UndefVarError: build_convergence_observables`.

### Pitfall 2: Breaking test assertions on observable count/names
**What goes wrong:** `test_convergence.jl` asserts `length(observables) == NUM_QUBITS + 1` and `names == ["ZZ_12", "ZZ_23", "ZZ_34", "ZZ_41", "H"]`. The new function returns 5 observables with names `["H", "Mz", "XX_avg", "YY_avg", "ZZ_avg"]`.
**Why it happens:** The old convergence builders returned per-bond ZZ pairs (4 of them for n=4) plus H. The new function returns averaged correlations.
**How to avoid:** Update all test assertions that check observable count or names.
**Warning signs:** Test failures with `5 != 5` or name mismatch assertions.

### Pitfall 3: Liouvillian basis mismatch
**What goes wrong:** The Liouvillian from `construct_lindbladian` is in a specific eigenbasis (Hamiltonian eigenbasis for TimeDomain, Trotter eigenbasis for TrotterDomain). The observables and initial state must be in the SAME basis.
**Why it happens:** This was the root cause of the quick-20 fix and has been a persistent source of confusion.
**How to avoid:** For TimeDomain: observables use `hamiltonian.eigvecs` for basis transform, initial state is in eigenbasis (psi0[end]=1 gives highest-energy eigenstate). For TrotterDomain: use `trotter.eigvecs`.
**Warning signs:** Fitted gaps that are orders of magnitude off, or overlaps that are all near zero.

### Pitfall 4: Initial state choice
**What goes wrong:** Using ground state (psi0[1]=1) at high beta produces near-zero decay signal because the ground state is already close to the Gibbs state.
**Why it happens:** At beta=10, the Gibbs state population is exponentially concentrated on the ground state. The ground state *is* essentially the steady state.
**How to avoid:** Use the excited state (psi0[end]=1 in eigenbasis). This is far from Gibbs, producing large decay amplitudes. This was already fixed in quick-24.
**Warning signs:** All fitted gaps are near zero or fitting fails to converge.

### Pitfall 5: ARPACK shift-invert subtlety
**What goes wrong:** ARPACK's `eigs(L, nev=2, which=:SM)` may fail or converge to wrong eigenvalues for the Liouvillian because the steady-state eigenvalue is exactly 0.
**Why it happens:** `:SM` looks for smallest magnitude, but the zero eigenvalue causes numerical issues.
**How to avoid:** The existing code uses shift-invert mode: `eigs(L, nev=2, sigma=1e-9*(1+1im), tol=1e-12)`. This finds eigenvalues nearest to the shift (near zero) via (L - sigma*I)^{-1}. Keep this approach.
**Warning signs:** ARPACK fails to converge or returns eigenvalues far from zero.

### Pitfall 6: Removing CrossValidationResult breaks test_gap_estimation.jl
**What goes wrong:** `test_gap_estimation.jl` has an entire "Cross-Validation" testset (lines 206-305) testing `cross_validate_gap` and `CrossValidationResult`. Deleting the source without updating the tests causes test failures.
**Why it happens:** Tight coupling between source and test code.
**How to avoid:** Delete the Cross-Validation testset entirely when removing the source code. The validation logic moves to the script.

## Code Examples

### Example 1: Current `run_lindbladian` ARPACK usage (from src/furnace.jl)
```julia
# Source: src/furnace.jl lines 14-22
shift = 1e-9 * (1 + 1im)
eigvals_near_zero, eigvecs_near_zero = eigs(liouv, nev=2, sigma=shift, tol=1e-12)
sorted_permutation_eigen = sortperm(abs.(real.(eigvals_near_zero)))
ss_index = sorted_permutation_eigen[1]   # Smallest
gap_index = sorted_permutation_eigen[2]  # Second smallest
spectral_gap = eigvals_near_zero[gap_index]
```

### Example 2: Current `eigen()` dense approach (from experiments/validate_gap_estimation.jl)
```julia
# Source: experiments/validate_gap_estimation.jl lines 114-135
eig = eigen(L)
sorted_idx = sortperm(abs.(real.(eig.values)))
ss_idx = sorted_idx[1]   # eigenvalue closest to 0
gap_idx = sorted_idx[2]  # second eigenvalue (spectral gap)
spectral_gap = eig.values[gap_idx]
```

### Example 3: Eigenbasis overlap computation (from experiments/eigenmode_decomposition.jl)
```julia
# Source: experiments/eigenmode_decomposition.jl lines 96-123
F = eigen(L)
lambda = F.values
V = F.vectors
perm = sortperm(abs.(real.(lambda)))
lambda = lambda[perm]
V = V[:, perm]

# Initial state expansion
rho0_vec = reshape(rho0, :)
alpha = V \ rho0_vec

# Overlap coefficients for observable O
O_vec = reshape(O, :)
c = zeros(ComplexF64, length(lambda))
for k in 1:length(lambda)
    o_k = dot(O_vec, V[:, k])
    c[k] = o_k * alpha[k]
end
```

### Example 4: Recommended overlap analysis function structure
```julia
# Recommended new function in gap_estimation.jl
struct OverlapAnalysisResult
    eigenvalues::Vector{ComplexF64}          # Sorted Liouvillian eigenvalues
    exact_gap::Float64                        # abs(real(eigenvalues[2]))
    observable_names::Vector{String}
    overlap_coefficients::Matrix{ComplexF64}  # n_obs x n_modes
    gap_mode_overlap::Vector{Float64}         # |c_1| for each observable
    relative_gap_overlap::Vector{Float64}     # |c_1| / sum |c_k| for each obs
end

function eigenbasis_overlap_analysis(
    L::Matrix{<:Complex},
    observables::Vector{<:Matrix{<:Complex}},
    observable_names::Vector{String},
    rho0::Matrix{<:Complex},
)
    # Full eigendecomposition (dense, only for small systems)
    F = eigen(L)
    perm = sortperm(abs.(real.(F.values)))
    lambda = F.values[perm]
    V = F.vectors[:, perm]

    # Initial state expansion
    alpha = V \ reshape(rho0, :)

    # Compute overlap coefficients
    n_obs = length(observables)
    n_modes = length(lambda)
    coeffs = zeros(ComplexF64, n_obs, n_modes)
    for (i, O) in enumerate(observables)
        O_vec = reshape(O, :)
        for k in 1:n_modes
            coeffs[i, k] = dot(O_vec, V[:, k]) * alpha[k]
        end
    end

    exact_gap = abs(real(lambda[2]))
    gap_overlap = [abs(coeffs[i, 2]) for i in 1:n_obs]
    rel_overlap = [abs(coeffs[i, 2]) / sum(abs.(coeffs[i, 2:end])) for i in 1:n_obs]

    return OverlapAnalysisResult(lambda, exact_gap, observable_names, coeffs,
                                  gap_overlap, rel_overlap)
end
```

## Inventory of Code to Delete

### Files to DELETE entirely:
1. `experiments/validate_gap_estimation.jl` -- replaced by new unified script
2. `experiments/run_gap_validation.jl` -- replaced by new unified script
3. `experiments/run_sweep.jl` -- KMS-vs-GNS sweep, out of scope for gap validation
4. `experiments/eigenmode_decomposition.jl` -- absorbed into library function

### Functions to DELETE from source:
1. `cross_validate_gap(::SpectralGapResult, ::LindbladianResult)` -- from `gap_estimation.jl`
2. `cross_validate_gap(::SpectralGapResult, ::Complex)` -- from `gap_estimation.jl`
3. `CrossValidationResult` struct -- from `gap_estimation.jl`
4. `build_convergence_observables(::HamHam, ::Int)` -- from `convergence.jl`
5. `build_convergence_observables_trotter(::HamHam, ::TrottTrott, ::Int)` -- from `convergence.jl`
6. `build_total_magnetization(::HamHam, ::Int; trotter)` -- from `convergence.jl`

### Functions to RENAME:
1. `build_gap_estimation_observables` -> `build_preset_trajectory_observables` (or similar)

### Exports to UPDATE in `QuantumFurnace.jl`:
1. Remove: `CrossValidationResult`, `cross_validate_gap`
2. Remove: `build_convergence_observables`, `build_convergence_observables_trotter`, `build_total_magnetization`, `build_gap_estimation_observables`
3. Add: `build_preset_trajectory_observables` (new name)
4. Add: `eigenbasis_overlap_analysis` (new function), `OverlapAnalysisResult` (new struct)

### Tests to UPDATE:
1. `test_gap_estimation.jl`: Remove entire "Cross-Validation" testset (lines 200-305). Remove `_make_test_spectral_gap_result` helper if only used there. Update observable name assertions if renamed. Add eigenbasis_overlap_analysis tests.
2. `test_convergence.jl`: Replace all `build_convergence_observables(TEST_HAM, NUM_QUBITS)` calls with `build_preset_trajectory_observables(TEST_HAM, NUM_QUBITS)`. Update assertions on observable count (was NUM_QUBITS+1=5 per-bond, now 5 preset), names (was ZZ_ij+H, now H+Mz+XX_avg+YY_avg+ZZ_avg). Delete testsets 2 and 3 (old builder tests). Testset 21+22 (existing gap estimation observable tests) cover the same logic -- keep or merge.

## Entanglement Analysis

### Functions that STAY (core pipeline, untouched):
- `estimate_spectral_gap` -- stays, just uses the new observable builder name internally
- `_select_best_observable` -- stays
- `fit_exponential_decay` / `FitResult` -- stays
- `run_observable_trajectories` / `ObservableTrajectoryResult` -- stays
- `run_trajectories_convergence` / `run_trajectories_adaptive` -- stays
- `run_lindbladian` / `construct_lindbladian` -- stays
- `SpectralGapResult` -- stays
- `ConvergenceData` -- stays
- `_gibbs_in_trotter_basis`, `_compute_gibbs_observable_values`, `_windowed_relative_change` -- stays

### Internal reference in estimate_spectral_gap:
Line 158-159 in `gap_estimation.jl` calls `build_gap_estimation_observables` by name. This must be updated to the new function name.

### Convergence runner tests:
`test_convergence.jl` testsets 8-18 use `build_convergence_observables` as fixture input for `run_trajectories_convergence` and `run_trajectories_adaptive`. These just need the observables as matrix inputs -- any builder function works. The new function returns different observables (H, Mz, XX_avg, YY_avg, ZZ_avg instead of ZZ_12, ZZ_23, ZZ_34, ZZ_41, H), so assertions about `observable_names`, `observable_values` matrix sizes, and energy convergence need updating.

**Critical detail:** The convergence integration test (testset 8) asserts `conv_data.trace_distances[end] < conv_data.trace_distances[1]` and energy convergence. These are physics-level assertions about the trajectory dynamics, not about the observable builder. They should still pass with different observables.

## ARPACK vs eigen: What to Expect

For n=4, the Liouvillian is 256x256. Both `eigs(L, nev=2, sigma=shift, tol=1e-12)` and `eigen(L)` should produce the same spectral gap to within ~1e-10.

Potential discrepancy sources:
1. **ARPACK tol**: The shift-invert tolerance is 1e-12, which is very tight. Should match eigen to ~1e-11.
2. **Shift perturbation**: The shift `1e-9 * (1+1im)` perturbs the matrix slightly. For the gap mode (eigenvalue ~-0.1), this is negligible.
3. **Eigenvalue sorting**: ARPACK sorts by proximity to sigma, while eigen sorts by... whatever Julia's LAPACK binding uses. Both need manual re-sorting by `abs(real(eigenvalue))`.

**Expected outcome:** Agreement to ~1e-10 or better. If they disagree significantly, it indicates a bug in the ARPACK call setup.

## Eigenbasis Overlap: What the Existing Script Found

The `eigenmode_decomposition.jl` script (n=4, beta=5, TrotterDomain) found:
- Some observables have sign cancellation in their gap-mode coefficients
- The effective initial decay rate can be above or below the true gap depending on the observable
- ZZ_avg showed sign cancellation (positive + negative coefficients near the gap mode), explaining its tendency to under-estimate the gap

This analysis needs to be repeated with the **validation parameters** (beta=10, TimeDomain) for the Phase 25 results to be self-consistent.

## Open Questions

1. **What is the actual residual factor at beta=10 with 20k trajectories?**
   - Prior quick-fix results: factor ~0.87x (n=4) and ~0.92x (n=6) at beta=10 with only 1000 trajectories
   - With 20k trajectories, statistical noise should be negligible, isolating the systematic bias
   - If relative error exceeds 1e-2, the eigenbasis overlap analysis must explain why
   - Recommendation: Run the script and report the result; the RESEARCH cannot predict this

2. **Should `build_convergence_observables` be preserved for backward compatibility?**
   - It's an exported function, so removing it is an API break
   - But it's only used in `run_sweep.jl` (being deleted) and test fixtures
   - The user explicitly said "delete overlapping/redundant observable builder functions"
   - Recommendation: Delete it. This is pre-release code. The new function serves the same purpose.

3. **Should the overlap analysis function take a dense Liouvillian or construct it internally?**
   - Taking a pre-constructed L is simpler and more composable
   - The caller already knows how to build L (via `construct_lindbladian`)
   - Recommendation: Take L as input; don't couple the overlap function to Lindbladian construction

## Sources

### Primary (HIGH confidence)
- Direct codebase inspection of all source files, test files, and experiment scripts
- Quick-fix summaries 22, 23, 24 (read and cross-referenced against actual code)
- CONTEXT.md from user discussion session

### Codebase Evidence
- `src/furnace.jl` lines 14-17: ARPACK `eigs` call with shift-invert
- `src/convergence.jl` lines 14-162: All 4 observable builder functions
- `src/gap_estimation.jl` lines 207-314: CrossValidationResult + cross_validate_gap (to be deleted)
- `src/gap_estimation.jl` lines 141-205: estimate_spectral_gap (to be kept, line 158 references build_gap_estimation_observables)
- `experiments/eigenmode_decomposition.jl` lines 96-123: Eigenbasis overlap math (to be absorbed into library)
- `test/test_convergence.jl` lines 54-91: Tests for old observable builders (to be updated)
- `test/test_gap_estimation.jl` lines 206-305: Cross-validation tests (to be deleted)

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- no new libraries, only existing deps
- Architecture: HIGH -- direct inspection of all affected files, clear dependency graph
- Pitfalls: HIGH -- derived from actual bugs found in quick-22/23/24 and code inspection
- Overlap analysis math: HIGH -- verified against existing eigenmode_decomposition.jl implementation

**Research date:** 2026-02-18
**Valid until:** 2026-03-18 (stable domain, no external dependencies changing)
