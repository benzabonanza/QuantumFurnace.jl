# Phase 24: Cross-Validation - Research

**Researched:** 2026-02-17
**Domain:** Cross-validation of trajectory-fitted spectral gap against exact Liouvillian eigenvalues; validation script for Heisenberg chains
**Confidence:** HIGH

## Summary

Phase 24 closes the v1.3 milestone by establishing quantitative trust in the trajectory-based gap estimation method (Phase 23). The core task is threefold: (1) implement a `cross_validate_gap` function that compares a trajectory-fitted `SpectralGapResult` against the exact Liouvillian spectral gap from `LindbladianResult`, using `abs(real(spectral_gap))` as the exact reference (locked decision), (2) add a warning when the exact eigenvalue has significant imaginary part (`|Im/Re| > 0.1`), indicating oscillatory decay incompatible with pure exponential fitting, and (3) create a validation script demonstrating agreement for n=4 and n=6 Heisenberg chains.

All building blocks exist. The `estimate_spectral_gap` function (Phase 23) returns `SpectralGapResult` with gap, CI, and per-observable fits. The `run_lindbladian` function computes exact eigenvalues and returns `LindbladianResult` with `spectral_gap::Complex{T}`. The Hamiltonians for n=3, 4, 5, 6 are available as pre-computed BSON files, and `HamHam` can also be constructed from scratch for uniform chains. The `construct_lindbladian` function builds the full Liouvillian matrix -- feasible for n=4 (256x256) and n=6 (4096x4096) but not n=8 (65536x65536, ~32 GB). No new dependencies are needed.

The main design challenge is where `cross_validate_gap` lives and what it returns. It should be a thin comparison function in `gap_estimation.jl` (co-located with `estimate_spectral_gap` and `SpectralGapResult`), returning a lightweight result struct (`CrossValidationResult`) with relative error, absolute error, whether the fitted gap falls within CI, and the imaginary part warning. The validation script goes in `experiments/` (matching the `run_sweep.jl` pattern for standalone scripts), not in the test suite (it takes minutes to run for n=6).

**Primary recommendation:** Add `cross_validate_gap(estimated::SpectralGapResult, exact::LindbladianResult)` to `gap_estimation.jl` with a `CrossValidationResult` return type. Add unit tests for the function logic. Create `experiments/validate_gap_estimation.jl` as a standalone script that runs n=4 and n=6 cross-validation.

## Standard Stack

### Core (already in project)
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| LinearAlgebra | stdlib | `eigen()` for Liouvillian spectral decomposition | Already used throughout |
| LsqFit.jl | 0.15 | Behind `fit_exponential_decay` (via Phase 21) | Already added |
| Arpack | -- | Behind `run_lindbladian` (sparse eigensolve) | Already used, but for n<=6 dense `eigen()` via `construct_lindbladian` is preferred |

### Supporting (already in project)
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| Printf | stdlib | Formatted output in validation script | Already used in `run_sweep.jl` |
| Dates | stdlib | Timestamps in script output | Already used |
| BSON | -- | Hamiltonian loading via `load_hamiltonian` | Already used |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| `run_lindbladian` for exact gap | `construct_lindbladian` + `eigen()` directly | `run_lindbladian` uses Arpack's `eigs(nev=2)` which is iterative. For n<=6, full dense `eigen()` on `construct_lindbladian` output is simpler and more reliable. The validation script should use `construct_lindbladian` + `eigen()` directly (matching `test_dm_detailed_balance.jl` pattern). |
| `CrossValidationResult` struct | Just return a NamedTuple | Struct is better for export, documentation, and Aqua compliance. NamedTuples are fine for internal use but this is public API. |
| Separate `cross_validation.jl` file | Add to existing `gap_estimation.jl` | `gap_estimation.jl` is currently 192 lines. Adding ~80-100 lines of cross-validation keeps it under 300 lines. The functions are tightly related (both operate on `SpectralGapResult`). Separate file would be overengineering. |

**Installation:** No new packages needed.

## Architecture Patterns

### Recommended Project Structure
```
src/
  gap_estimation.jl     # Phase 23: SpectralGapResult, estimate_spectral_gap
                        # Phase 24 ADD: CrossValidationResult, cross_validate_gap
  QuantumFurnace.jl     # Add exports: CrossValidationResult, cross_validate_gap

test/
  test_gap_estimation.jl  # Phase 24 ADD: tests for cross_validate_gap

experiments/
  validate_gap_estimation.jl  # Phase 24 NEW: standalone validation script for n=4,6
```

### Pattern 1: cross_validate_gap Function
**What:** A comparison function that takes trajectory-fitted results and exact Liouvillian results, computes relative error, and warns about imaginary eigenvalue components.
**When to use:** After running both `estimate_spectral_gap` and exact Liouvillian eigendecomposition on the same system.
**Example:**
```julia
# Source: codebase analysis (furnace.jl line 22, gap_estimation.jl line 37)
struct CrossValidationResult
    fitted_gap::Float64              # from SpectralGapResult.gap
    exact_gap::Float64               # abs(real(exact_result.spectral_gap))
    relative_error::Float64          # |fitted - exact| / exact
    absolute_error::Float64          # |fitted - exact|
    within_ci::Bool                  # exact_gap in [gap_ci[1], gap_ci[2]]
    imaginary_ratio::Float64         # |Im(spectral_gap)| / |Re(spectral_gap)|
    imaginary_warning::Bool          # true if imaginary_ratio > 0.1
end

function cross_validate_gap(
    estimated::SpectralGapResult,
    exact_result::LindbladianResult,
)
    exact_eigenvalue = exact_result.spectral_gap  # Complex{T}
    exact_gap = abs(real(exact_eigenvalue))        # LOCKED DECISION

    fitted_gap = estimated.gap
    relative_error = abs(fitted_gap - exact_gap) / exact_gap
    absolute_error = abs(fitted_gap - exact_gap)
    within_ci = estimated.gap_ci[1] <= exact_gap <= estimated.gap_ci[2]

    im_ratio = abs(imag(exact_eigenvalue)) / abs(real(exact_eigenvalue))
    im_warning = im_ratio > 0.1

    if im_warning
        @warn "Exact eigenvalue has significant imaginary part" im_ratio exact_eigenvalue
    end

    return CrossValidationResult(
        fitted_gap, exact_gap, relative_error, absolute_error,
        within_ci, im_ratio, im_warning,
    )
end
```

### Pattern 2: Exact Gap Extraction via construct_lindbladian + eigen
**What:** For small systems (n<=6), build the full Liouvillian and do dense eigendecomposition. Do NOT use `run_lindbladian` (which uses Arpack's iterative solver and prints output).
**When to use:** In the validation script, and optionally in tests.
**Example:**
```julia
# Source: test_dm_detailed_balance.jl lines 29-36 (established pattern)
liouv_config = make_small_liouv_config(TimeDomain(); with_coherent=false)
L = construct_lindbladian(SMALL_JUMPS, liouv_config, SMALL_HAM)
eig = eigen(L)

# Sort by |Re(lambda)| to find steady state (index 1) and gap mode (index 2)
sorted_idx = sortperm(abs.(real.(eig.values)))
ss_idx = sorted_idx[1]      # eigenvalue closest to 0 (steady state)
gap_idx = sorted_idx[2]     # second eigenvalue (spectral gap)

exact_gap_eigenvalue = eig.values[gap_idx]
exact_gap = abs(real(exact_gap_eigenvalue))
```

### Pattern 3: Validation Script Structure
**What:** A standalone script (not in test suite) that runs full gap estimation + exact eigendecomposition for n=4 and n=6.
**When to use:** The script runs for minutes (especially n=6 with many trajectories). It belongs in `experiments/`, not `test/`.
**Example structure:**
```julia
# experiments/validate_gap_estimation.jl
# Run via: julia --project=. experiments/validate_gap_estimation.jl

using QuantumFurnace
using LinearAlgebra
using Printf

function validate_system(num_qubits::Int, beta::Float64; ntraj::Int, seed::Int)
    # 1. Build Hamiltonian (uniform Heisenberg XXX chain)
    terms = Vector{Vector{Matrix{ComplexF64}}}([[X, X], [Y, Y], [Z, Z]])
    coeffs = [1.0, 1.0, 1.0]
    hamiltonian = HamHam(terms, coeffs, num_qubits, beta; periodic=true)

    # 2. Build jump operators
    jumps = _build_jumps(hamiltonian, num_qubits)

    # 3. Exact gap via Liouvillian
    liouv_config = LiouvConfig(...)
    L = construct_lindbladian(jumps, liouv_config, hamiltonian)
    exact_result = _extract_liouvillian_result(L, hamiltonian)

    # 4. Trajectory gap estimation
    config = ThermalizeConfig(...)
    psi0 = zeros(ComplexF64, 2^num_qubits); psi0[1] = 1.0
    estimated = estimate_spectral_gap(jumps, config, psi0, hamiltonian; ntraj, seed)

    # 5. Cross-validate
    cv = cross_validate_gap(estimated, exact_result)
    _print_results(cv, num_qubits, beta)
    return cv
end
```

### Pattern 4: LindbladianResult Construction from Dense Eigen
**What:** The validation script needs a `LindbladianResult` to pass to `cross_validate_gap`, but should use dense `eigen()` rather than `run_lindbladian` (which uses Arpack and prints output).
**Options:**
  a. Construct `LindbladianResult` directly from `eigen()` output (the struct fields are `liouvillian`, `fixed_point`, `gap_mode`, `spectral_gap`).
  b. Have `cross_validate_gap` accept `Complex{Float64}` directly instead of `LindbladianResult`.

**Recommendation:** Option (b) is more flexible -- `cross_validate_gap(estimated::SpectralGapResult, exact_eigenvalue::Complex)`. However, the success criteria say `cross_validate_gap(estimated, exact_result)`, implying it takes `LindbladianResult`. Use option (a): construct `LindbladianResult` from dense eigen. This is straightforward since `LindbladianResult` is `@kwdef` with 4 fields.

Actually, reviewing the success criteria more carefully: "User can call `cross_validate_gap(estimated, exact_result)` comparing fitted gap against `abs(real(exact_result.spectral_gap))`." The `exact_result` must have a `.spectral_gap` field, which `LindbladianResult` has. But the function could also just take a `Complex` eigenvalue. The most reusable approach: accept either `LindbladianResult` or `Complex` via method dispatch.

### Anti-Patterns to Avoid
- **Using `run_lindbladian` in validation script:** It uses Arpack's iterative solver (less reliable for small systems), prints output, and returns via `eigs(nev=2, sigma=shift)` which can miss eigenvalues. Dense `eigen()` is exact for n<=6.
- **Including validation script in `test/runtests.jl`:** The n=6 case takes minutes. It should be a standalone script in `experiments/`.
- **Using `abs(spectral_gap)` instead of `abs(real(spectral_gap))`:** Locked decision. For complex eigenvalues, `abs()` would incorporate the imaginary part into the gap, but the trajectory-based exponential decay only measures the real part of the decay rate.
- **Failing silently when imaginary part is large:** The whole point of VAL-02 is to warn the user. If `|Im/Re| > 0.1`, the pure exponential fit `A * exp(-gap * t) + C` cannot capture the oscillatory component `cos(Im * t)`, so the gap estimate is inherently less reliable.
- **Putting CrossValidationResult in structs.jl:** Like `SpectralGapResult`, it should be co-located with its producer in `gap_estimation.jl`.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Gap estimation from trajectories | Manual orchestration | `estimate_spectral_gap(jumps, config, psi0, ham; ...)` | Phase 23 -- tested and working |
| Observable construction | Manual matrix building | `build_gap_estimation_observables(ham, n; trotter)` | Phase 20 -- basis transforms correct |
| Exponential fitting | Manual LsqFit calls | `fit_exponential_decay(times, values; skip_initial)` | Phase 21 -- auto-init, bounds, CI |
| Liouvillian construction | Custom matrix assembly | `construct_lindbladian(jumps, config, ham)` | Existing -- tested across all domains |
| Eigendecomposition | Custom spectral analysis | `eigen(L)` from LinearAlgebra stdlib | Standard, exact for dense matrices |
| Hamiltonian construction | Manual matrix building | `HamHam(terms, coeffs, n, beta; periodic=true)` | Existing constructor handles rescaling, eigvals, Gibbs |
| Jump operator creation | Manual construction | Follow `run_sweep.jl` `build_trotter_system` or `test_helpers.jl` `make_test_system` pattern | Established pattern with normalization |

**Key insight:** Phase 24 adds ~100-150 lines of new code: `CrossValidationResult` struct (~10 lines), `cross_validate_gap` function (~30 lines), tests (~50 lines), and the validation script (~100 lines). Everything else is composition of existing building blocks.

## Common Pitfalls

### Pitfall 1: abs(spectral_gap) vs abs(real(spectral_gap))
**What goes wrong:** Using `abs(spectral_gap)` instead of `abs(real(spectral_gap))` gives a different (larger) exact gap when the eigenvalue has imaginary part.
**Why it happens:** The eigenvalue is complex (e.g. `-0.5 + 0.3im`). `abs()` gives `sqrt(0.25 + 0.09) = 0.583`, while `abs(real())` gives `0.5`. The trajectory exponential fit only captures the `0.5` decay rate.
**How to avoid:** This is a locked decision: always use `abs(real(spectral_gap))`. The `cross_validate_gap` function enforces this.
**Warning signs:** Relative error is consistently positive (fitted gap < exact gap) when imaginary part is non-zero.

### Pitfall 2: Eigenvalue Selection for Spectral Gap
**What goes wrong:** Selecting the wrong eigenvalue as the spectral gap (e.g., picking a degenerate zero eigenvalue or sorting by `abs(lambda)` instead of `abs(real(lambda))`).
**Why it happens:** The Liouvillian has one eigenvalue at 0 (steady state). The spectral gap is the second-smallest `|Re(lambda)|`. But some eigenvalues might have very small real part and large imaginary part.
**How to avoid:** Sort by `abs(real(eig.values))`. The steady state is the one closest to 0. The gap eigenvalue is the next one. This matches the pattern in `run_lindbladian` (`sortperm(abs.(real.(eigvals_near_zero)))`) and `test_dm_detailed_balance.jl` (`argmin(abs.(real.(eig.values)))`).
**Warning signs:** The "exact gap" is unreasonably large or doesn't match expectations from DM evolution timescale.

### Pitfall 3: Config Mismatch Between Trajectory and Liouvillian
**What goes wrong:** The trajectory estimation uses a `ThermalizeConfig` (with delta, mixing_time) while the Liouvillian uses `LiouvConfig`. If physics parameters differ (beta, sigma, a, b, domain), the gaps will not match.
**Why it happens:** Two separate config objects are constructed. Easy to have different `num_energy_bits`, `w0`, `a`, `b`, etc.
**How to avoid:** In the validation script, derive both configs from the same physical parameters. Use a helper function that creates matched `LiouvConfig` and `ThermalizeConfig` from a shared parameter set.
**Warning signs:** Relative error is large and doesn't decrease with more trajectories.

### Pitfall 4: n=6 Liouvillian Construction Time and Memory
**What goes wrong:** n=6 gives a 4096x4096 Liouvillian (16.8M entries, ~128 MB as ComplexF64). `construct_lindbladian` allocates this plus workspace matrices. `eigen()` on this matrix takes ~10-30 seconds.
**Why it happens:** Dense matrix eigendecomposition is O(n^3) where n=4096.
**How to avoid:** This is expected and acceptable. The validation script should print progress messages ("Constructing Liouvillian for n=6...") so the user knows it's working. n=8 (65536x65536, ~32 GB) is NOT feasible and should NOT be attempted.
**Warning signs:** If the script tries n=8, it will either OOM or take hours.

### Pitfall 5: Domain Approximation Error Masking Gap Agreement
**What goes wrong:** For TimeDomain or TrotterDomain, the Liouvillian has domain approximation errors. The exact gap from the Liouvillian is the gap of the approximated system, not the exact Bohr-domain system. However, the trajectory simulation operates in the same approximated domain, so they should still agree.
**Why it happens:** Confusion between "exact gap of the Bohr Lindbladian" and "exact gap of the TimeDomain Lindbladian." Both are exact eigenvalues of their respective Liouvillians.
**How to avoid:** Use the same domain for both trajectory estimation and Liouvillian construction. Document that cross-validation confirms the trajectory method agrees with the Liouvillian within statistical error, not that it recovers the "true" theoretical gap.
**Warning signs:** Using BohrDomain for exact eigenvalues but TimeDomain for trajectories gives unexplained discrepancy.

### Pitfall 6: Insufficient Trajectories for n=6
**What goes wrong:** With n=6 (dim=64), observable means have more noise per trajectory. The fitted gap has wide confidence intervals that make cross-validation inconclusive.
**Why it happens:** Signal-to-noise ratio decreases with system size. The spectral gap may also be smaller for larger systems, making the decay slower and harder to fit.
**How to avoid:** Use more trajectories for n=6 than n=4 (e.g., ntraj=5000 for n=6 vs ntraj=1000 for n=4). Use longer total_time to capture more of the decay curve. Consider `skip_initial` to exclude early multi-exponential transient.
**Warning signs:** R-squared < 0.7 or gap_ci spans zero.

### Pitfall 7: Mixing Time Too Short
**What goes wrong:** If total_time is much shorter than the gap inverse (1/gap), the observable barely decays, and the exponential fit extracts noise rather than signal.
**Why it happens:** The gap is unknown before estimation. Using `mixing_time` from config might be too short.
**How to avoid:** For the validation script, set total_time to be large enough (e.g., 5-10x expected 1/gap). For the Heisenberg chain at beta=10, typical gaps are O(0.01-0.1) in rescaled units, so total_time=50-100 is appropriate.
**Warning signs:** Fitted gap is near zero with large standard error.

### Pitfall 8: Hamiltonian Type in Validation Script
**What goes wrong:** Tests use pre-computed disordered Hamiltonians (via `_load_test_hamiltonian`), but the validation script description says "Heisenberg chains." Using disordered Hamiltonians might have different gap behavior than uniform chains.
**Why it happens:** The existing BSON files are disordered. The `run_sweep.jl` script builds uniform chains via `HamHam(terms, coeffs, n, beta; periodic=true)`.
**How to avoid:** The validation script should build uniform Heisenberg chains (matching `run_sweep.jl` pattern), not load disordered ones. This gives reproducible, disorder-independent results. For tests, the existing SMALL system (disordered) is fine since we're testing the function logic, not physical agreement.
**Warning signs:** Gap values vary between runs if disorder realization changes.

## Code Examples

Verified patterns from the existing codebase:

### Building Uniform Heisenberg Hamiltonian (from run_sweep.jl)
```julia
# Source: experiments/run_sweep.jl lines 48-52
terms = Vector{Vector{Matrix{ComplexF64}}}([[X, X], [Y, Y], [Z, Z]])
coeffs = [1.0, 1.0, 1.0]
hamiltonian = HamHam(terms, coeffs, num_qubits, beta; periodic=true)
```

### Building Jump Operators (from run_sweep.jl)
```julia
# Source: experiments/run_sweep.jl lines 59-76
function build_jumps(hamiltonian::HamHam, num_qubits::Int)
    jump_paulis = [[X], [Y], [Z]]
    n_jumps = length(jump_paulis) * num_qubits
    norm_factor = sqrt(n_jumps)

    jumps = JumpOp[]
    for pauli in jump_paulis
        for site in 1:num_qubits
            op = Matrix(pad_term(pauli, num_qubits, site)) ./ norm_factor
            in_eigen = hamiltonian.eigvecs' * op * hamiltonian.eigvecs
            push!(jumps, JumpOp(op, in_eigen, op == transpose(op), op == op'))
        end
    end
    return jumps
end
```

### Exact Eigendecomposition Pattern (from test_dm_detailed_balance.jl)
```julia
# Source: test/test_dm_detailed_balance.jl lines 29-36, 53-74
liouv_config = make_liouv_config(domain)
L = construct_lindbladian(jumps, liouv_config, hamiltonian)

eig = eigen(L)
ss_idx = argmin(abs.(real.(eig.values)))
# Second smallest |Re(lambda)| is the spectral gap eigenvalue
sorted_by_real = sortperm(abs.(real.(eig.values)))
gap_idx = sorted_by_real[2]
spectral_gap = eig.values[gap_idx]
```

### Constructing LindbladianResult from Dense Eigen
```julia
# Source: furnace.jl lines 18-38 (adapted for dense eigen)
dim = 2^num_qubits
L = construct_lindbladian(jumps, liouv_config, hamiltonian)
eig = eigen(L)
sorted_idx = sortperm(abs.(real.(eig.values)))

ss_idx = sorted_idx[1]
gap_idx = sorted_idx[2]

ss_vec = eig.vectors[:, ss_idx]
ss_dm = reshape(ss_vec, dim, dim)
hermitianize!(ss_dm)
ss_dm ./= tr(ss_dm)

gap_vec = eig.vectors[:, gap_idx]
gap_mode = reshape(gap_vec, dim, dim)

exact_result = LindbladianResult(
    liouvillian = L,
    fixed_point = ss_dm,
    gap_mode = gap_mode,
    spectral_gap = eig.values[gap_idx],
)
```

### Using @warn for Imaginary Part Warning
```julia
# Source: src/time_domain.jl line 12, src/trajectories.jl line 916 (existing @warn pattern)
if im_warning
    @warn "Exact eigenvalue has significant imaginary part (|Im/Re| = $(round(im_ratio; digits=3))). " *
          "Pure exponential fit may not capture oscillatory decay." exact_eigenvalue=exact_eigenvalue
end
```

### Config Factory Pattern (from test_helpers.jl)
```julia
# Source: test/test_helpers.jl lines 227-243, 249-270
# For validation script: create matched LiouvConfig and ThermalizeConfig
function make_matched_configs(num_qubits, beta, sigma; domain, with_coherent, delta, mixing_time)
    shared = (
        num_qubits = num_qubits,
        with_coherent = with_coherent,
        with_linear_combination = true,
        domain = domain,
        beta = beta,
        sigma = sigma,
        a = beta / 30.0,
        b = 0.4,
        num_energy_bits = 12,
        w0 = 0.05,
        t0 = 2pi / (2^12 * 0.05),
        num_trotter_steps_per_t0 = 10,
    )
    liouv = LiouvConfig(; shared...)
    therm = ThermalizeConfig(; shared..., mixing_time=mixing_time, delta=delta)
    return liouv, therm
end
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Manual comparison of trajectory results vs eigenvalues in notebooks | Programmatic `cross_validate_gap` function | Phase 24 | Reproducible, documented comparison |
| Visual inspection of decay curves | Quantitative relative error + CI overlap check | Phase 24 | Trust metric for publication |
| No warning about oscillatory eigenvalues | Automatic warning when `\|Im/Re\| > 0.1` | Phase 24 | Users know when pure exponential fit is unreliable |

**Important context:** The existing `main_liouv.jl` simulation script already uses `abs(real(liouv_result.spectral_gap))` for printing the gap (line 130). This confirms the locked decision about using `abs(real(...))`.

## Open Questions

1. **Should `cross_validate_gap` accept a bare `Complex` eigenvalue in addition to `LindbladianResult`?**
   - What we know: The success criteria say `cross_validate_gap(estimated, exact_result)` where `exact_result` has `.spectral_gap`. This implies `LindbladianResult`.
   - What's unclear: Whether users might want to cross-validate against an eigenvalue from a different source (e.g., manual Arpack call).
   - Recommendation: Implement two methods: one taking `LindbladianResult`, one taking `Complex`. The `LindbladianResult` method just extracts `.spectral_gap` and forwards. This adds ~3 lines and makes the API more flexible.

2. **What domain should the validation script use for cross-validation?**
   - What we know: TimeDomain (without coherent) is the simplest non-trivial domain. EnergyDomain is slightly more accurate. BohrDomain is exact but not trajectory-compatible.
   - What's unclear: Whether to validate multiple domains or just one.
   - Recommendation: Use TimeDomain with `with_coherent=false` for primary validation. It matches the test helpers pattern (`make_small_thermalize_config(TimeDomain(); with_coherent=false)`) and has a known domain approximation error that is small for 12-bit energy grid. Optionally show a second run with EnergyDomain or BohrDomain for comparison.

3. **What beta value for the validation script?**
   - What we know: Tests use beta=10. The run_sweep uses beta=5, 10, 20. Higher beta gives smaller gaps (slower decay) requiring more time and trajectories.
   - What's unclear: Which beta gives the cleanest cross-validation demonstration.
   - Recommendation: Use beta=10 (matching test infrastructure defaults). Document that beta=5 gives faster convergence (larger gap) while beta=20 is harder.

4. **What tolerance for "agreement" in VAL-03?**
   - What we know: The success criterion says "fitted gap within confidence interval of exact gap (or within documented tolerance)."
   - What's unclear: What if the CI is very wide but includes the exact gap? Or CI is tight but misses by a small amount?
   - Recommendation: Primary check: `within_ci == true`. Secondary check: relative_error < 0.3 (30%). Document that tighter agreement requires more trajectories. The validation script should report both metrics.

5. **Should `CrossValidationResult` be exported?**
   - What we know: `SpectralGapResult` is exported. `FitResult` is exported.
   - Recommendation: Yes, export both `CrossValidationResult` and `cross_validate_gap` for user access.

## Sources

### Primary (HIGH confidence)
- `src/gap_estimation.jl` -- `SpectralGapResult` struct (lines 37-50), `estimate_spectral_gap` function (lines 127-191). Verified by reading source code.
- `src/structs.jl` -- `LindbladianResult` struct (lines 258-264), showing `spectral_gap::Complex{T}` field. Verified by reading source code.
- `src/furnace.jl` -- `run_lindbladian` function (lines 1-39), showing eigenvalue extraction pattern with `sortperm(abs.(real.(...)))`. Verified by reading source code.
- `simulations/main_liouv.jl` -- Line 130: `abs(real(liouv_result.spectral_gap))` confirming locked decision. Verified by reading source code.
- `test/test_dm_detailed_balance.jl` -- Dense eigendecomposition pattern with `eigen(L)` and `argmin(abs.(real.(...)))`. Verified by reading source code.
- `experiments/run_sweep.jl` -- Hamiltonian/jump construction pattern for uniform Heisenberg chains. Verified by reading source code.
- `test/test_helpers.jl` -- Config factory functions, SMALL system constants, test infrastructure. Verified by reading source code.
- `.planning/REQUIREMENTS.md` -- VAL-01, VAL-02, VAL-03 requirement definitions. Verified by reading source.
- `.planning/ROADMAP.md` -- Phase 24 success criteria. Verified by reading source.
- `.planning/STATE.md` -- Locked decisions including `abs(real(spectral_gap))`. Verified by reading source.

### Secondary (MEDIUM confidence)
- System size analysis: n=4 gives 256x256 Liouvillian (trivial), n=6 gives 4096x4096 (feasible, ~10-30s for eigen), n=8 gives 65536x65536 (infeasible, ~32 GB). Computed from 2^n dimensions.

### Tertiary (LOW confidence)
- None needed -- all Phase 24 work is internal composition of verified building blocks.

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- no new dependencies, all building blocks exist and are tested
- Architecture: HIGH -- `cross_validate_gap` is a thin comparison function; `CrossValidationResult` follows `SpectralGapResult` pattern
- Pitfalls: HIGH -- all pitfalls derived from actual code patterns and physics understanding (eigenvalue selection, domain mismatch, imaginary components)
- Code examples: HIGH -- all examples use verified function signatures from source code

**Research date:** 2026-02-17
**Valid until:** No expiry -- codebase-internal research, no external dependency version concerns
