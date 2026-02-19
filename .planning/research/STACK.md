# Stack Research: Spectral Gap Refinement Diagnostics

**Domain:** Advanced spectral analysis, bootstrap statistics, two-exponential fitting, similarity transforms, diagnostic visualization for quantum Lindbladian simulation
**Researched:** 2026-02-19
**Confidence:** HIGH (existing deps verified from codebase; new capabilities analyzed against official docs and Julia ecosystem)

## Scope

This stack research covers ONLY the additions needed for v1.4 Spectral Gap Refinement -- the diagnostic and improved estimation capabilities described in the spectral-gap-refinements-instructions.md (Tasks 1.1--5.2).

Specifically:
1. Anti-Hermitian defect computation (similarity transform with rho^{+/-1/4})
2. Effective rate plot lambda_eff(t) computation
3. Bootstrap resampling over trajectories for error bars
4. Two-exponential fitting with robust initialization (Prony method)
5. Richardson extrapolation for delta-convergence
6. Automatic fitting window selection (SNR-based t_max, stability-based t_min)
7. Symmetry sector labeling on Lindbladian eigenvectors
8. Dashboard/summary figure generation

It does NOT re-research the existing stack. The v1.3 STACK.md covered LsqFit.jl (already in [deps]), Arpack.jl, single-exponential fitting, and log-linear initial guess. All of those remain valid and unchanged.

## Verdict: Zero New Production Dependencies

The existing production dependency set is sufficient for ALL v1.4 features. The key insight: every new capability (two-exponential fitting, bootstrap, Richardson extrapolation, matrix fourth roots, symmetry labeling) can be implemented with LinearAlgebra + LsqFit + existing stdlib, with plotting handled via the already-present extras dependency on Plots.jl. No new `[deps]` entries needed.

---

## Existing Stack (Relevant to v1.4 Features)

### Already in [deps] -- Used Directly by New Features

| Existing Dep | v1.4 Role | Sufficient? |
|---|---|---|
| **LsqFit.jl** (0.15) | Two-exponential fitting: same `curve_fit` with 5-param model `c1*exp(-g1*t) + c2*exp(-g2*t)` + bounds. Confidence intervals via `confint`. | YES -- no API changes needed. Two-exponential is just a different model function passed to the same `curve_fit`. |
| **LinearAlgebra** (stdlib) | Matrix fourth root via eigendecomposition: `F = eigen(Hermitian(rho)); rho_quarter = F.vectors * Diagonal(F.values .^ 0.25) * F.vectors'`. Also: `opnorm()` for anti-Hermitian defect, `eigen()` for full Lindbladian diagonalization (n<=6), `Hermitian()` wrapper. | YES -- `eigen(Hermitian(...))` guarantees real eigenvalues for the Gibbs state; raising diagonal to 1/4 power is trivial. |
| **Arpack.jl** (0.5.4) | Leading 20-30 eigenvalues via shift-invert `eigs(L, nev=30, sigma=shift)`. Already used for 2 eigenvalues in `run_lindbladian`; just increase `nev`. | YES -- shift-invert mode already working for Lindbladian. Increasing nev from 2 to 20-30 is the only change. |
| **Optim.jl** (1) | NOT directly used for fitting (LsqFit handles that), but available if custom loss functions needed for fitting window optimization. | Not needed for v1.4. Keep as-is. |
| **SparseArrays** (stdlib) | Lindbladian stored as dense (dim^2 x dim^2) for n<=6. Sparse not needed at these sizes. | No change. |

### Already in [extras] -- Used in Scripts/Tests

| Existing Extra | v1.4 Role | Sufficient? |
|---|---|---|
| **Plots.jl** (1) | Dashboard figures (Task 5.1). Multi-panel layout via `plot(p1, p2, ..., layout=(rows, cols))`. Supports PNG, PDF, SVG output. Already in `[extras]` and `[compat]` with bounds. | YES for diagnostic scripts. See plotting section below for detailed analysis. |
| **StatsBase** (0.34) | `mean`, `std`, `sample` (with replacement) for bootstrap resampling. `sample(1:N, N; replace=true)` gives bootstrap indices. | YES -- `StatsBase.sample` with `replace=true` is the standard Julia bootstrap primitive. |
| **HypothesisTests** (0.11) | Statistical validation of gap estimates. | YES -- no change. |
| **StableRNGs** (1) | Reproducible bootstrap seeds. | YES -- no change. |
| **Statistics** (stdlib) | `mean`, `std`, `var` for bootstrap statistics, SNR computation. Already resolved in Manifest. | YES -- stdlib, always available. |

---

## Feature-by-Feature Stack Analysis

### 1. Anti-Hermitian Defect (Task 1.2): Matrix Fourth Root

**What's needed:** Compute rho^{1/4} and rho^{-1/4} for the Gibbs state to form the similarity-transformed generator D = rho^{-1/4} * L[rho^{1/4} * (.) * rho^{1/4}] * rho^{-1/4}.

**Stack decision: Use LinearAlgebra.eigen() on Hermitian Gibbs state. No new dependency.**

**Confidence: HIGH** -- Verified from Julia docs and codebase analysis.

**Rationale:** The Gibbs state rho_beta is diagonal in the energy eigenbasis (already computed in `HamHam.gibbs`). The project stores it as a `Hermitian{Complex{T}, Matrix{Complex{T}}}`. Two paths exist:

**Path A (preferred): Exploit diagonal structure in eigenbasis.**
The `HamHam` struct stores `eigvals` and `eigvecs`. The Gibbs state in the eigenbasis is `Diagonal(exp.(-beta .* eigvals) ./ Z)`. Fourth root is simply `Diagonal((exp.(-beta .* eigvals) ./ Z) .^ 0.25)`. Transform back to computational basis if needed via `eigvecs * D * eigvecs'`. This is O(dim^2), exact, and numerically stable (all Gibbs eigenvalues are strictly positive).

```julia
# In energy eigenbasis: rho is diagonal
gibbs_eigs = exp.(-beta .* hamiltonian.eigvals)
gibbs_eigs ./= sum(gibbs_eigs)
rho_quarter_diag = Diagonal(gibbs_eigs .^ 0.25)
rho_inv_quarter_diag = Diagonal(gibbs_eigs .^ (-0.25))

# In computational basis (if needed for Lindbladian that's in computational basis):
V = hamiltonian.eigvecs
rho_quarter = V * rho_quarter_diag * V'
rho_inv_quarter = V * rho_inv_quarter_diag * V'
```

**Path B (for TrotterDomain): General eigendecomposition.**
In TrotterDomain, the Gibbs state is NOT diagonal in the Trotter eigenbasis. The Lindbladian and Gibbs state are both in the Trotter basis. Use:

```julia
F = eigen(Hermitian(Matrix(gibbs_in_trotter_basis)))
rho_quarter = F.vectors * Diagonal(F.values .^ 0.25) * F.vectors'
rho_inv_quarter = F.vectors * Diagonal(F.values .^ (-0.25)) * F.vectors'
```

This is safe because `Hermitian()` guarantees real eigenvalues. The Gibbs state is positive definite (all eigenvalues > 0), so the quarter-power is well-defined and real.

**Numerical concern: near-zero eigenvalues.**
For high beta (low temperature), some Gibbs eigenvalues approach zero, making rho^{-1/4} diverge. At n=6 with beta=1 (the current test case), the smallest Gibbs eigenvalue is ~exp(-beta * E_max) / Z, which is well above machine epsilon. For beta >> 1, a regularization floor `max(eigval, epsilon)` may be needed, but this is a v1.4+ concern. **No new package needed.**

### 2. Effective Rate Plot lambda_eff(t) (Task 2.1): Pure Arithmetic

**What's needed:** Compute lambda_eff(t) = -ln|Delta(t+tau)/Delta(t)| / tau from the trajectory-averaged signal.

**Stack decision: Pure Julia arithmetic. No dependency at all.**

**Confidence: HIGH** -- This is a 15-line function using only division, abs, log, and array indexing.

```julia
function effective_rate(Delta::Vector{Float64}, dt::Float64; lag::Int=3)
    n = length(Delta)
    tau = lag * dt
    t_vals = [(i-1)*dt for i in 1:(n-lag)]
    lambda_eff = Vector{Float64}(undef, n-lag)
    for i in 1:(n-lag)
        d_now = Delta[i]
        d_later = Delta[i + lag]
        if d_now != 0.0 && d_later != 0.0 && sign(d_now) == sign(d_later)
            lambda_eff[i] = -log(abs(d_later / d_now)) / tau
        else
            lambda_eff[i] = NaN
        end
    end
    return t_vals, lambda_eff
end
```

### 3. Bootstrap Resampling (Task 2.2): StatsBase.sample + Custom Loop

**What's needed:** Resample N_traj trajectories with replacement, recompute averaged signal, recompute lambda_eff, collect statistics.

**Stack decision: Use StatsBase.sample (already in test extras) for index sampling. No new dependency.**

**Confidence: HIGH** -- StatsBase.sample with `replace=true` is the standard Julia primitive for bootstrap.

**Critical design consideration:** The current `run_observable_trajectories` does NOT store per-trajectory measurements -- it accumulates into `mean_data` and divides by `ntraj`. For bootstrap, we need per-trajectory data OR we need a modified trajectory runner.

**Two approaches (architecture decision, not a stack decision):**

**Approach A: Store per-trajectory measurements.** Add a `measurements_per_traj::Array{Float64, 3}` field (n_obs x n_saves x n_traj) to ObservableTrajectoryResult. Memory: for n=6, 8 observables, 500 save points, 10000 trajectories = 8 * 500 * 10000 * 8 bytes = 320 MB. Feasible for n<=6 validation but NOT for production use at larger n.

**Approach B: Batch-level bootstrap.** Run K independent trajectory batches (e.g., K=200 batches of ntraj/K trajectories each). Store per-batch means. Bootstrap resample at the batch level. Memory: 8 * 500 * 200 * 8 bytes = 6.4 MB. This is the approach used in lattice QCD (jackknife/bootstrap over configurations).

**Recommended: Approach B (batch-level bootstrap).** The trajectory runner already supports deterministic per-trajectory seeding via `Xoshiro(seed + traj_id)`. Run batches sequentially, store per-batch mean arrays, then bootstrap over batch indices using `StatsBase.sample(1:K, K; replace=true)`.

```julia
using StatsBase: sample
# After collecting batch_means::Vector{Matrix{Float64}}  (K matrices of size n_obs x n_saves)
n_boot = 200
boot_gaps = Vector{Float64}(undef, n_boot)
for b in 1:n_boot
    idx = sample(1:K, K; replace=true)
    boot_mean = mean(batch_means[idx])  # requires Statistics.mean
    # compute lambda_eff or fit from boot_mean
    boot_gaps[b] = ...
end
gap_se = std(boot_gaps)
```

**Why NOT Bootstrap.jl:** Bootstrap.jl provides infrastructure for standard bootstrap workflows, but our use case is non-standard (resampling trajectory batches, then re-fitting). The custom loop above is 10 lines and gives full control. Adding Bootstrap.jl would be dependency bloat for no benefit.

### 4. Two-Exponential Fitting (Task 3.1): LsqFit.jl with 5-Parameter Model

**What's needed:** Fit `Delta(t) = c1 * exp(-g1 * t) + c2 * exp(-g2 * t)` with constraint `0 < g1 < g2`.

**Stack decision: LsqFit.jl (already in [deps]). No new dependency.**

**Confidence: HIGH** -- Verified that LsqFit supports arbitrary model functions with bounded parameters.

```julia
# Two-exponential model
_two_exp_model(t, p) = @. p[1] * exp(-p[2] * t) + p[3] * exp(-p[4] * t) + p[5]

# Parameters: [c1, g1, c2, g2, offset]
# Bounds: g1 > 0, g2 > g1 (enforce g2 > some_minimum)
lower = [-Inf, 0.0, -Inf, 0.0, -Inf]
upper = [Inf, Inf, Inf, Inf, Inf]

fit = curve_fit(_two_exp_model, times, data, p0; lower=lower, upper=upper)
```

**The challenge is initialization, not the fitting library.** Two-exponential fits are notoriously sensitive to initial conditions. The instructions describe two initialization strategies:

**Strategy A: Effective-rate-informed initialization.**
Read plateau value from lambda_eff plot as g1_init, early-time value as g2_init. Simple, requires lambda_eff to be computed first.

**Strategy B: Prony two-point method.**
Pick two time points, form ratios, solve quadratic for decay constants. More automated but sensitive to noise in the chosen points.

Both strategies are ~20 lines of pure Julia arithmetic. No external package needed. The Prony method in particular does NOT need SignalDecomposition.jl or any signal processing package -- it's a 2x2 linear algebra problem:

```julia
function _prony_two_exp_init(t, Delta; t_frac_a=0.2, t_frac_b=0.5)
    n = length(t)
    ia = max(1, round(Int, t_frac_a * n))
    ib = max(ia+2, round(Int, t_frac_b * n))
    im = div(ia + ib, 2)
    dt_half = t[im] - t[ia]

    r1 = Delta[im] / Delta[ia]
    r2 = Delta[ib] / Delta[im]

    # z1 + z2 = r1 + r2,  z1 * z2 = r1 * r2  (Prony relations)
    S = r1 + r2
    P = r1 * r2
    disc = S^2 - 4*P
    if disc < 0
        # Complex roots: fall back to single-exponential init
        return nothing
    end
    z1 = (S + sqrt(disc)) / 2
    z2 = (S - sqrt(disc)) / 2

    # Convert to decay rates
    g1 = -log(max(abs(z1), 1e-10)) / dt_half
    g2 = -log(max(abs(z2), 1e-10)) / dt_half
    if g1 > g2; g1, g2 = g2, g1; end  # ensure g1 < g2

    # Amplitudes from linear system
    # Delta[ia] = c1*exp(-g1*t[ia]) + c2*exp(-g2*t[ia])
    # Delta[ib] = c1*exp(-g1*t[ib]) + c2*exp(-g2*t[ib])
    E = [exp(-g1*t[ia]) exp(-g2*t[ia]); exp(-g1*t[ib]) exp(-g2*t[ib])]
    c = E \ [Delta[ia]; Delta[ib]]

    return [c[1], g1, c[2], g2, 0.0]  # offset = 0 initially
end
```

### 5. Richardson Extrapolation (Task 3.2): Pure Arithmetic

**What's needed:** Combine gap estimates at delta and delta/2: `gap_rich = 2*gap(delta/2) - gap(delta)`.

**Stack decision: One-line formula. No dependency.**

**Confidence: HIGH** -- This is literally one line of arithmetic.

```julia
gap_richardson = 2 * gap_half_delta - gap_delta
sigma_richardson = sqrt(4 * sigma_half_delta^2 + sigma_delta^2)
```

**Why NOT Richardson.jl:** The Richardson.jl package (v1.4.0, Dec 2020) is designed for adaptive extrapolation of scalar functions with automatic convergence detection. Our use case is a one-shot linear extrapolation from two or three points -- literally `2*f(h/2) - f(h)`. Adding a package dependency for this would be absurd. Richardson.jl would be useful if we were doing higher-order Richardson with automatic order detection, but the instructions explicitly describe the simple two-point formula.

### 6. Automatic Fitting Window Selection (Task 3.3): SNR + Stability

**What's needed:**
- t_max: last time where SNR(t) = |Delta(t)| / sigma_Delta(t) > 3
- t_min: stability test -- sweep t_min, look for plateau in fitted g1

**Stack decision: Pure Julia loops and array operations. Statistics.std for variance. No new dependency.**

**Confidence: HIGH** -- This is iterative fitting with varying parameters, all using existing LsqFit infrastructure.

SNR computation requires per-timepoint variance from the bootstrap batches. This feeds directly into the batch-level data collection discussed in section 3 above.

### 7. Symmetry Sector Labeling (Task 4.1): LinearAlgebra + Eigenbasis Indexing

**What's needed:** Compute S_z quantum numbers for each energy eigenstate, then label each Lindbladian eigenvector by the Delta_S_z sectors it has support in.

**Stack decision: Pure Julia. Pauli Z matrix already exists in the codebase (see `src/hamiltonian.jl` exports `Z`). Total S_z = sum of Z_i/2 for each eigenstate.**

**Confidence: HIGH** -- The only operation is computing Z_total = sum(pad_term(Z, i, n)) in the energy eigenbasis, then checking which (i,j) blocks of each Lindbladian eigenvector (reshaped to dim x dim) have non-zero weight.

```julia
# Compute S_z for each energy eigenstate
Z_total = sum(pad_term([Z], i, n_qubits) for i in 1:n_qubits)  # existing pad_term
Sz_values = [real(hamiltonian.eigvecs[:, k]' * Z_total * hamiltonian.eigvecs[:, k]) / 2
             for k in 1:dim]

# For each Lindbladian eigenvector reshaped to dim x dim:
# entry (i,j) corresponds to |E_i><E_j| with Delta_Sz = Sz[i] - Sz[j]
function label_symmetry_sector(eigvec, dim, Sz_values; threshold=1e-6)
    R = reshape(eigvec, dim, dim)
    delta_sz_weights = Dict{Float64, Float64}()
    for i in 1:dim, j in 1:dim
        w = abs2(R[i, j])
        w < threshold && continue
        dsz = round(Sz_values[i] - Sz_values[j], digits=6)
        delta_sz_weights[dsz] = get(delta_sz_weights, dsz, 0.0) + w
    end
    return delta_sz_weights
end
```

### 8. Dashboard/Summary Figures (Task 5.1): Plots.jl

**What's needed:** 7-panel diagnostic figure with eigenvalue spectrum, defect metrics, overlap coefficients, effective rate plot, delta-convergence, fit overlay, and t_min stability.

**Stack decision: Use Plots.jl (already in [extras] with compat bounds). Generate figures in diagnostic scripts, not in production src/.**

**Confidence: HIGH** -- Plots.jl is already a project dependency (extras section), already has compat bounds (`Plots = "1"`), and the docs Project.toml already uses it.

**Why Plots.jl and not CairoMakie:**

| Factor | Plots.jl | CairoMakie |
|--------|----------|------------|
| Already in Project.toml | YES (extras) | NO |
| Multi-panel layout | `plot(p1, p2, ..., layout=(r,c))` -- simple | `Figure()`/`Axis` -- more flexible but more verbose |
| Publication quality | Good with GR backend (default) | Better (vector graphics focus) |
| Dependency weight | GR backend is lightweight | CairoMakie 0.15.8 brings Makie 0.24.8, Cairo, FreeType, GeometryBasics + many transitive deps |
| Learning curve for team | Low (already used in docs) | Higher (scene graph model) |
| Julia compat | Plots 1.x works with Julia 1.11+ | CairoMakie 0.15.8 requires Julia 1.3+ (fine) but Makie ecosystem adds ~30 packages |

**Verdict: Use Plots.jl.** The diagnostic figures are for internal validation and paper drafts, not for interactive dashboards. Plots.jl's `plot(; layout=...)` handles the 7-panel figure adequately. If publication-quality vector graphics are needed later, CairoMakie can be considered, but adding ~30 new transitive dependencies for slightly nicer fonts is not justified at this stage.

**Plots.jl usage pattern for dashboard:**

```julia
using Plots

# Panel A: Eigenvalue spectrum
p1 = scatter(real.(eigs), imag.(eigs); xlabel="Re(lambda)", ylabel="Im(lambda)",
             marker_z=delta_sz_labels, title="Lindbladian Spectrum")

# Panel D: Effective rate plot
p4 = plot(t_eff, lambda_eff; ribbon=sigma_eff, xlabel="t", ylabel="lambda_eff(t)",
          title="Effective Rate")
hline!([exact_gap]; linestyle=:dash, label="Exact gap")

# Combine into dashboard
dashboard = plot(p1, p2, p3, p4, p5, p6, p7; layout=(4, 2), size=(1200, 1600))
savefig(dashboard, "diagnostic_dashboard.png")
```

---

## What Arpack.jl Can Do for Extended Eigenvalue Computation

The existing `run_lindbladian` uses `eigs(liouv, nev=2, sigma=shift, tol=1e-12)` to get the two smallest-magnitude eigenvalues. For v1.4 Task 1.1, we need the leading 20-30 eigenvalues and their eigenvectors.

**Stack decision: Keep Arpack.jl. Increase nev to 20-30. No change to dependency.**

**Confidence: HIGH** -- Arpack's shift-invert mode supports arbitrary nev up to dim-1.

```julia
# Current (v1.3):
eigvals, eigvecs = eigs(liouv, nev=2, sigma=shift, tol=1e-12)

# New (v1.4):
eigvals, eigvecs = eigs(liouv, nev=30, sigma=shift, tol=1e-12)
# Returns 30 eigenvalues nearest to shift, with their eigenvectors.
```

For n=6, the Lindbladian is 4096 x 4096 (dense). Arpack with shift-invert at this size takes ~seconds and is well within capabilities. No need for KrylovKit.jl or any other eigenvalue solver.

**Why NOT switch to KrylovKit.jl:**
KrylovKit.jl (v0.8+) does NOT support shift-invert mode. Its `eigsolve` only finds extremal eigenvalues (largest or smallest magnitude/real part). For finding eigenvalues near zero in a non-Hermitian Lindbladian, shift-invert is essential. KrylovKit's documentation explicitly warns: "since no (shift-and)-invert is used, this will only be successful if you somehow know that eigenvalues close to zero are also close to the periphery of the spectrum." For our Lindbladian, eigenvalues near zero are interior eigenvalues, NOT extremal. Arpack's Fortran-based shift-invert is the right tool here.

KrylovKit.jl would be relevant for n=8+ (dim=65536) where the Lindbladian is too large to factorize for shift-invert. That's deferred to a future milestone (Task 5.3 is explicitly deferred in the instructions).

---

## Dependencies to NOT Add

| Package | Why Evaluated | Why NOT Adding |
|---------|---------------|----------------|
| **Bootstrap.jl** | Bootstrap resampling for error bars on lambda_eff and gap estimates | Our bootstrap is non-standard (batch-level trajectory resampling + re-fitting). A 10-line loop with `StatsBase.sample` gives full control. Bootstrap.jl's API is designed for standard statistic-on-sample workflows. |
| **Richardson.jl** (v1.4.0) | Richardson extrapolation for delta-convergence | One-line formula: `2*gap(delta/2) - gap(delta)`. Adding a dependency for this is absurd overhead. Richardson.jl's adaptive polynomial extrapolation is overkill for a fixed two-point linear combination. |
| **KrylovKit.jl** | Alternative sparse eigenvalue solver | No shift-invert mode. Cannot find eigenvalues near zero for non-Hermitian Lindbladian. Arpack is the right tool for n<=6. |
| **CairoMakie.jl** (0.15.8) | Publication-quality vector graphics | Adds ~30 transitive dependencies (Makie 0.24.8, Cairo, FreeType, GeometryBasics, etc.). Plots.jl already in extras and sufficient for diagnostic figures. Consider for paper polish milestone only. |
| **SignalDecomposition.jl** | Prony/ESPRIT method for multi-exponential decomposition | Immature Julia package. Prony two-point initialization is ~20 lines of arithmetic. No external package needed. |
| **Measurements.jl** | Uncertainty propagation | Heavyweight for propagating bootstrap error through Richardson formula. Manual formula `sigma_rich = sqrt(4*sigma1^2 + sigma2^2)` is simpler and more transparent. |
| **Statistics.jl** (stdlib) as production dep | mean, std for bootstrap/SNR | Already available as stdlib. Use in scripts/tests only (already works via `using Statistics`). Not needed as `[deps]` entry since diagnostic functions can take pre-computed statistics as arguments. |

---

## Recommended Stack Changes Summary

### Production Dependencies [deps]: NO CHANGES

The existing dependency set is sufficient:

| Package | Status | Role in v1.4 |
|---------|--------|--------------|
| LsqFit (0.15) | Already present | Two-exponential fitting (new model function, same API) |
| LinearAlgebra (stdlib) | Already present | Matrix fourth root, opnorm for defect, full eigendecomposition |
| Arpack (0.5.4) | Already present | 20-30 eigenvalue extraction (increase nev) |
| SparseArrays (stdlib) | Already present | Lindbladian storage (no change) |

### Test/Script Dependencies [extras]: NO CHANGES

| Package | Status | Role in v1.4 |
|---------|--------|--------------|
| Plots (1) | Already present | Diagnostic dashboard figures |
| StatsBase (0.34) | Already present | Bootstrap index sampling via `sample(replace=true)` |
| HypothesisTests (0.11) | Already present | Gap estimate statistical validation |
| StableRNGs (1) | Already present | Reproducible bootstrap seeds |

### Project.toml Changes: NONE

```toml
# No changes to [deps], [compat], [extras], or [targets]
# Everything needed is already declared.
```

---

## Integration Points with Existing Code

### New Files (src/)

| New File | Dependencies Used | Purpose |
|----------|-------------------|---------|
| `src/diagnostics.jl` | LinearAlgebra (eigen, opnorm, Hermitian, Diagonal) | Anti-Hermitian defect, symmetry sector labeling, exact reference data |
| `src/effective_rate.jl` | None (pure arithmetic) | lambda_eff(t) computation, SNR computation |
| `src/two_exp_fit.jl` | LsqFit (curve_fit, confint, stderror, coef) | Two-exponential model definition, Prony initialization, fitting window selection |
| `src/richardson.jl` | None (pure arithmetic) | Richardson extrapolation formula and error propagation |

### New Files (scripts/ or experiments/)

| New File | Dependencies Used | Purpose |
|----------|-------------------|---------|
| `experiments/diagnostics/run_diagnostics.jl` | Plots, StatsBase, Statistics | Full diagnostic pipeline: run trajectories, compute diagnostics, generate dashboard |
| `experiments/diagnostics/bootstrap_analysis.jl` | StatsBase (sample), Statistics (mean, std) | Batch-level bootstrap for error bars |

### Modified Files (src/)

| File | Change | Stack Impact |
|------|--------|--------------|
| `src/QuantumFurnace.jl` | Add `include` for new files, add exports for new public functions | No new `using` statements needed |
| `src/furnace.jl` | Modify `run_lindbladian` to support `nev` parameter (default 2, allow 20-30) | Same Arpack API, just parameterize `nev` |
| `src/fitting.jl` | Keep existing `fit_exponential_decay` unchanged. New two-exponential functions go in separate file. | No modification to existing LsqFit usage |

---

## Numerical Considerations

### Matrix Fourth Root Stability

For the anti-Hermitian defect computation, rho^{-1/4} involves inverting the fourth root of potentially small eigenvalues:

| Scenario | Smallest Gibbs Eigenvalue | rho^{-1/4} Magnitude | Risk |
|----------|---------------------------|----------------------|------|
| n=4, beta=1 | ~0.04 | ~2.2 | None |
| n=6, beta=1 | ~0.003 | ~4.3 | None |
| n=6, beta=5 | ~1e-10 | ~5600 | Moderate (amplifies numerical errors in Lindbladian) |
| n=6, beta=10 | ~1e-22 | ~1e5 | HIGH (similarity transform numerically meaningless) |

**Recommendation:** Compute the condition number of the similarity transform (max/min of rho^{-1/4} diagonal entries). If > 1e8, warn that the defect metric is unreliable. For the Heisenberg chain at beta=1, this is not an issue.

### Two-Exponential Fit Conditioning

The condition of a two-exponential fit depends on the ratio g2/g1. If g2/g1 < 2, the two exponentials are nearly indistinguishable and the Jacobian becomes ill-conditioned, leading to:
- Non-convergence of Levenberg-Marquardt
- Huge confidence intervals
- Rate estimates that swap (g1 and g2 exchange roles)

**Mitigation built into LsqFit:** The `maxIter` parameter (default 1000) prevents infinite loops. The `confint` function returns Inf when the Jacobian is rank-deficient (the codebase already handles this for single-exponential in `fitting.jl` lines 201-208). The same `try/catch LinearAlgebra.SingularException` pattern should be applied to two-exponential fits.

### Bootstrap Sample Size

For K batches, the number of distinct bootstrap samples is `binomial(2K-1, K)`. With K=200 and N_boot=200, each bootstrap resample gives an independent gap estimate. The standard error on the bootstrap SE itself is ~SE/sqrt(2*N_boot), so N_boot=200 gives ~5% precision on the error estimate. Sufficient for diagnostic purposes.

---

## Alternatives Considered

| Category | Recommended | Alternative | Why Not |
|----------|-------------|-------------|---------|
| Eigenvalues | Arpack.jl (shift-invert) | KrylovKit.jl | No shift-invert; cannot target interior eigenvalues |
| Eigenvalues | Arpack.jl | ArnoldiMethod.jl | Pure Julia but also no shift-invert mode |
| Two-exp fitting | LsqFit.jl | Optim.jl + manual Jacobian | LsqFit handles Jacobian, CIs, and convergence tracking. Reimplementing is error-prone. |
| Two-exp fitting | LsqFit.jl | Custom IRLS/Prony standalone | Prony only for initialization; LsqFit for refinement is standard practice |
| Bootstrap | StatsBase.sample + loop | Bootstrap.jl | Non-standard use case (batch resampling + refitting); 10-line loop is clearer |
| Richardson | Manual formula | Richardson.jl | One-line computation; package adds dependency for no benefit |
| Plotting | Plots.jl (GR backend) | CairoMakie.jl | ~30 extra transitive deps; Plots.jl already in project, sufficient quality |
| Plotting | Plots.jl | UnicodePlots.jl | No file output; terminal-only; not suitable for paper figures |
| Matrix power | eigen + diagonal ^ 0.25 | General `A^0.25` operator | Exploiting known diagonal structure in eigenbasis is more efficient and numerically stable |

---

## Version Compatibility

All v1.4 features use existing packages at their current versions. No version bumps needed.

| Package | Current Version | v1.4 API Used | Compatibility |
|---------|----------------|---------------|---------------|
| LsqFit | 0.15 | `curve_fit` with 5-param model, `confint`, `stderror`, `coef` | Same API as single-exponential. Verified from v0.15 docs. |
| Arpack | 0.5.4 | `eigs(L; nev=30, sigma=shift)` | `nev` parameter supported since Arpack inception. |
| LinearAlgebra | stdlib (Julia 1.11) | `eigen(Hermitian(...))`, `opnorm`, `Diagonal`, `norm` | All standard stdlib functions. |
| Plots | 1 | `plot(layout=...)`, `scatter`, `hline!`, `savefig` | Standard Plots.jl API. |
| StatsBase | 0.34 | `sample(1:N, N; replace=true)` | Core StatsBase function, stable since v0.30+. |
| Statistics | stdlib | `mean`, `std`, `var` | Standard stdlib functions. |

**No version conflicts.** No new packages to resolve. No compat bound changes needed.

---

## Sources

- [Julia LinearAlgebra documentation](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/) -- eigen, Hermitian, opnorm, Diagonal. HIGH confidence.
- [LsqFit.jl official documentation](https://julianlsolvers.github.io/LsqFit.jl/latest/) -- curve_fit API, bounds, confidence intervals. HIGH confidence.
- [LsqFit.jl Getting Started](https://julianlsolvers.github.io/LsqFit.jl/latest/getting_started/) -- Multi-parameter models, weighted fitting. HIGH confidence.
- [KrylovKit.jl eigenvalue problems](https://jutho.github.io/KrylovKit.jl/stable/man/eig/) -- Confirmed NO shift-invert mode. Quote: "since no (shift-and)-invert is used, this will only be successful if you somehow know that eigenvalues close to zero are also close to the periphery of the spectrum." HIGH confidence.
- [Arpack.jl documentation](https://arpack.julialinearalgebra.org/latest/) -- eigs API with sigma (shift-invert) and nev parameters. HIGH confidence.
- [Arpack.jl GitHub](https://github.com/JuliaLinearAlgebra/Arpack.jl) -- v0.5.4, Julia wrapper for arpack-ng Fortran library. HIGH confidence.
- [Richardson.jl GitHub](https://github.com/JuliaMath/Richardson.jl) -- v1.4.0 (Dec 2020). Only dependency: LinearAlgebra. Evaluated and rejected (overkill for one-line formula). HIGH confidence.
- [Bootstrap.jl GitHub](https://github.com/juliangehring/Bootstrap.jl) -- Evaluated and rejected (non-standard use case). HIGH confidence.
- [Plots.jl vs Makie.jl discussion](https://discourse.julialang.org/t/what-are-the-biggest-differences-between-makie-jl-and-plots-jl/76643) -- Ecosystem comparison. MEDIUM confidence (community discussion, not official).
- [CairoMakie.jl in Makie monorepo](https://github.com/MakieOrg/Makie.jl/tree/master/CairoMakie) -- v0.15.8, requires Makie 0.24.8. Evaluated and deferred. HIGH confidence.
- [Makie.jl releases](https://github.com/JuliaPlots/Makie.jl/releases) -- v0.24.8 (Dec 5, 2024). HIGH confidence.
- [StatsBase.jl documentation](https://juliastats.org/StatsBase.jl/stable/) -- sample function with replace keyword. HIGH confidence.
- [Prony's method Wikipedia](https://en.wikipedia.org/wiki/Prony%27s_method) -- Mathematical foundation for two-exponential initialization. HIGH confidence (textbook method).
- Eigendecomposition-based matrix power: f(A) = V * f(Lambda) * V^{-1} for diagonalizable A. Standard numerical linear algebra (Golub & Van Loan). HIGH confidence.
- [Julia Discourse: sparse eigenvalues](https://discourse.julialang.org/t/suggestions-needed-diagonalizing-large-hermitian-sparse-matrix/96580) -- Comparison of Arpack, KrylovKit, ArnoldiMethod for sparse eigenproblems. MEDIUM confidence.

---
*Stack research for: QuantumFurnace.jl v1.4 Spectral Gap Refinement Diagnostics*
*Researched: 2026-02-19*
