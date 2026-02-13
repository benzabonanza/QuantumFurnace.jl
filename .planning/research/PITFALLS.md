# Pitfalls Research

**Domain:** Quantum Lindbladian simulation (Julia package), trajectory unraveling, Python/Qiskit interop
**Researched:** 2026-02-13
**Confidence:** HIGH (codebase-grounded findings) / MEDIUM (interop and ecosystem patterns)

## Critical Pitfalls

### Pitfall 1: Trajectory Fallback Logic Silently Masks Probability Normalization Bugs

**What goes wrong:**
In `step_along_trajectory!` (both `TimeDomain/TrotterDomain` and `EnergyDomain` variants in `/Users/bence/code/QuantumFurnace.jl/src/trajectories.jl`), the dissipative jump branch has a `!chosen` fallback at lines 564-568 and 722-725 that applies the *last computed candidate state* when the cumulative probability scan fails to hit the target. This fallback is meant to catch floating-point rounding, but it also catches genuine probability normalization bugs. If `p_nojump + p_res + p_jump_total` does not sum to the expected total weight because `R` was computed inconsistently with the per-jump rates, the fallback absorbs the error and the simulation appears to work while producing statistically wrong steady states.

**Why it happens:**
The total jump probability `p_jump_total = delta * expR` is computed from the precomputed `R` matrix, while individual jump probabilities are computed on-the-fly from `base_prefactor * transition(w) * ||A_w psi||^2`. If the base_prefactor formula in `step_along_trajectory!` drifts from the formula in `precompute_R()`, the probabilities will not sum correctly, and the fallback silently picks an arbitrary jump instead of flagging the inconsistency.

**How to avoid:**
- Add a completeness assertion during framework construction: `|p_nojump + p_res + p_jump_total - 1.0| < epsilon` for a known test state.
- During development, add a counter for how often the fallback fires. If it fires more than ~0.1% of steps, it indicates a normalization bug, not rounding.
- Write a dedicated test that runs a trajectory for a 1-2 qubit system and checks that the trajectory-averaged density matrix converges to the same Gibbs state as the density matrix evolution path (`run_thermalization`).

**Warning signs:**
- The fallback branch fires frequently (add a debug counter).
- Trajectory-averaged rho deviates from DM-evolved rho by more than statistical sampling error.
- Changing `delta` changes the steady state (it should only change convergence speed).

**Phase to address:**
Testing/Validation phase -- immediately when adding trajectory tests.

---

### Pitfall 2: Cholesky Residual S Matrix Can Silently Go Non-PSD, Breaking CPTP Guarantee

**What goes wrong:**
The residual Kraus operator `U_residual` is computed via Cholesky factorization of `S = (2*alpha - delta)*R - alpha^2 * R^2` (lines 70-86 in `trajectories.jl`, and identically in `jump_workers.jl` lines 464-485). This Cholesky uses `check=false`, meaning if `S` has negative eigenvalues (beyond the `10*eps(Float64)` shift), the factorization silently produces garbage. The resulting CPTP channel is no longer trace-preserving, causing density matrices to develop negative eigenvalues or traces that drift from 1.0.

**Why it happens:**
`S` is `O(delta^2)`, so its eigenvalues are tiny. When `R` has eigenvalues near 1 (which happens for strong coupling or large jump operator norms), `alpha^2 * R^2` can dominate `(2*alpha - delta)*R`, making `S` indefinite. The `10*eps(Float64)` diagonal shift (~2.2e-15) is often too small to save this. Additionally, accumulated floating-point error in the `R` Hermitianization step can leave asymmetric residuals that push eigenvalues negative.

**How to avoid:**
- After computing `S`, check `minimum(real(eigvals(Hermitian(S))))` in tests and assert it is non-negative (or only negligibly negative, e.g., > -1e-14).
- Consider an adaptive epsilon shift: `eps_shift = max(10*eps(Float64), abs(minimum_eigval_of_S) + 1e-14)` computed from a single `eigmin` call (cheap for small matrices).
- For large-scale use, replace the `check=false` Cholesky with a pivoted Cholesky that truncates negative directions, or fall back to eigendecomposition-based square root.
- Add a trace-preservation test: `|tr(K0' * K0) + delta * tr(R) + tr(S) - dim| < epsilon`.

**Warning signs:**
- `issuccess(cholesky_S)` returns `false` (not currently checked because of `check=false`).
- Density matrix trace drifts from 1.0 over many steps.
- Negative eigenvalues appear in the evolved density matrix.
- NaN or Inf in trajectory state vectors.

**Phase to address:**
Core stabilization phase -- before trajectory validation since it affects all downstream results.

---

### Pitfall 3: Hardcoded `1e-12` Tolerance Used for Two Incompatible Purposes

**What goes wrong:**
The value `1e-12` appears 58 times across the source code in `/Users/bence/code/QuantumFurnace.jl/src/`. It serves as: (a) a frequency cutoff for the Hermitian half-grid optimization (`w_raw > 1e-12 && continue`), (b) numerical integration tolerance for `quadgk`, (c) NUFFT accuracy parameter, (d) a diagonal regularization shift in the LSI optimization, and (e) a truncation threshold for coherent term dictionaries. These are fundamentally different quantities -- a frequency resolution depends on `nu_min` (the smallest Bohr spacing), while an integration tolerance depends on the integrand's smoothness. Using the same magic number for all means that changing system parameters (larger beta, smaller energy gaps) can silently break one use while being fine for another.

**Why it happens:**
Quick prototyping. When all tests are on 4-qubit Heisenberg chains with beta=10, a single tolerance works. But `nu_min` scales as `O(2^{-n})` for n qubits, so for 10-12 qubits the smallest Bohr frequency could be close to or below `1e-12`, causing the half-grid optimization to misclassify zero-frequency contributions.

**How to avoid:**
- Define named constants in `constants.jl`: `FREQ_ZERO_TOL = nu_min / 10`, `INTEGRATION_TOL = 1e-12`, `NUFFT_EPS = 1e-12`, `EIGENVALUE_FLOOR = 1e-14`.
- The frequency zero tolerance should be derived from `hamiltonian.nu_min`, not hardcoded.
- Run a sweep test: for a range of qubit counts (2, 4, 6, 8), verify that the half-grid optimization gives the same Lindbladian as the full-grid computation.

**Warning signs:**
- Tests pass at 4 qubits but fail at 8+ qubits.
- Unexplained asymmetry in the Lindbladian (some Bohr frequency contributions are dropped).
- `nu_min` of the Hamiltonian is smaller than `1e-12` (check explicitly during validation).

**Phase to address:**
Refactoring phase -- before scaling tests to larger qubit counts.

---

### Pitfall 4: Density Matrix Thermalization and Trajectory Paths Use Different Delta Scaling Conventions

**What goes wrong:**
In `run_thermalization` (`furnace.jl` lines 80-145), when `rescale_by_inv_prob = true`, the coherent unitary is scaled by `delta_scale = 1.0 / p_jump` (line 109), and the jump weight scaling uses `gamma_norm_factor / jump_prob` (line 161 of `jump_workers.jl`). In the trajectory code (`trajectories.jl`), `build_trajectoryframework` does NOT apply this per-jump probability rescaling at all -- it applies all jumps deterministically in each step. This means the density matrix path and the trajectory path implement subtly different CPTP channels. Comparing their outputs will show disagreement that looks like a bug in the trajectory code but is actually a convention mismatch.

**Why it happens:**
The DM thermalization randomly selects one jump per step (line 119: `idx = rand(rng, 1:length(jumps))`), requiring a 1/p_jump rescaling to maintain the correct average Lindbladian. Trajectories apply the full Kraus map (all jumps) in each step. These are mathematically equivalent only if the scaling conventions are matched. The trajectory code was recently refactored (per commit history) and this reconciliation was likely not yet completed.

**How to avoid:**
- Before comparing DM and trajectory outputs, verify that both paths produce the same single-step channel on a known density matrix (identity/dim or a random state). Compute `E(rho) = K0 rho K0' + sum_jumps + residual` from both paths and check they agree.
- Document which convention each path uses in the function docstrings.
- Add a `verify_channel_equivalence(fw_dm, fw_traj, test_rho)` test helper.

**Warning signs:**
- DM thermalization converges to Gibbs state but trajectory average does not (or converges to a different state).
- Changing `ntraj` does not reduce the discrepancy (indicating systematic bias, not sampling noise).

**Phase to address:**
Trajectory validation phase.

---

### Pitfall 5: `unsafe_wrap` Pointer Arithmetic in Log-Sobolev Optimization Can Corrupt Memory

**What goes wrong:**
In `log_sobolev.jl` (lines 44-45, 148-149), `unsafe_wrap(Array, pointer(A_flat), (d, d))` and `unsafe_wrap(Array, pointer(A_flat, n_params * sizeof(Float64) + 1), (d, d))` create array views into a flat `Vector{Float64}` by raw pointer offset. The second call computes the offset as `n_params * sizeof(Float64) + 1`, but `pointer(A_flat, k)` uses *element* indexing, not byte indexing. If `n_params` is the number of Float64 elements, then `pointer(A_flat, n_params + 1)` is the correct offset, and the `* sizeof(Float64)` factor is wrong -- it overshoots by 8x, reading garbage memory.

**Why it happens:**
Confusion between Julia's `pointer(array, i)` which uses 1-based element indices, and C-style byte offsets. The `d` variable on line 44 is also undefined in the current code (should be `dim`) -- the closure captures `dim` from the outer scope, but `d` appears to be a typo or leftover. If `d != dim`, the wrap creates misshapen arrays, causing silent data corruption or crashes.

**How to avoid:**
- Replace `unsafe_wrap` with `reshape(view(...))` which is safe and equally fast for BLAS operations in modern Julia (1.9+).
- If `unsafe_wrap` must stay, add `@assert n_params == dim^2` and use `pointer(A_flat, n_params + 1)` without the `sizeof` multiplication.
- Fix the `d` vs `dim` inconsistency.
- Add a unit test: round-trip `A_flat -> A_real, A_imag -> reconstruct -> compare`.

**Warning signs:**
- Segmentation faults during LSI optimization.
- Gradient values that are NaN or unreasonably large.
- LSI optimization returning different results across runs with the same seed.
- `d` is referenced but not defined in the local scope (would be a runtime NameError if the closure does not capture it).

**Phase to address:**
Core stabilization -- this is a memory safety issue that should be fixed before any testing relies on LSI results.

---

### Pitfall 6: Coherent Term B Has an Unresolved `trotter` Keyword Bug in Trajectory Framework

**What goes wrong:**
In `build_trajectoryframework` (`trajectories.jl` line 53), the call `precompute_coherent_total_B(jumps, ham_or_trott, config, precomputed_data; trotter=trotter)` passes a keyword argument `trotter=trotter`, but the variable `trotter` is not defined in the function scope. The function signature has no `trotter` parameter. The `precompute_coherent_total_B` function in `coherent.jl` (line 14-39) also does not accept a `trotter` keyword. This means that if `config.with_coherent == true`, the trajectory framework will throw an `UndefVarError` at construction time.

**Why it happens:**
The trajectory framework was refactored from an older codebase where `trotter` was passed as a separate argument. The refactoring consolidated `ham_or_trott` but did not clean up the `trotter=trotter` keyword in this call path.

**How to avoid:**
- Run the trajectory code with `with_coherent = true` in a test to trigger the error immediately.
- Remove the `; trotter=trotter` keyword from line 53 (the `ham_or_trott` already carries the correct object).
- Add a CI test that exercises all four domain types with `with_coherent = true`.

**Warning signs:**
- Running trajectories with `with_coherent = true` crashes immediately.
- All trajectory tests only test `with_coherent = false` (hiding the bug).

**Phase to address:**
Testing phase -- a basic smoke test would catch this immediately.

---

## Technical Debt Patterns

| Shortcut | Immediate Benefit | Long-term Cost | When Acceptable |
|----------|-------------------|----------------|-----------------|
| ~200 lines of commented-out code in `trajectories.jl` and `jump_workers.jl` | Preserves old implementations for reference | Confuses contributors, makes diffs noisy, hides the actual API surface | Never in a published package; delete and rely on git history |
| `check=false` on Cholesky factorization | Avoids crash when S is marginally non-PSD | Silently produces wrong Kraus operators, breaks CPTP guarantee | Only with a runtime eigenvalue guard checking S is PSD first |
| Duplicated config structs (`LiouvConfig`/`LiouvConfigGNS`/`ThermalizeConfig`/`ThermalizeConfigGNS`) sharing 90% identical fields | Quick to add GNS variant | 4x maintenance burden for any field change, invitation for inconsistency | Pre-alpha only; refactor to parametric struct or composition |
| Immutable `HamHam` struct requiring `finalize_hamham` to add Gibbs/Bohr data | Avoids mutable state | Two-step construction is error-prone (forgetting to finalize); old tests show `hamiltonian.bohr_freqs = ...` attempts on immutable fields | Acceptable if enforced by a single constructor or builder pattern |
| Using `includet` (Revise) in test files instead of proper `using` | Fast iteration during development | Tests are not runnable via `Pkg.test()`, CI cannot execute them | During active development only; must be converted for CI |

## Integration Gotchas

| Integration | Common Mistake | Correct Approach |
|-------------|----------------|------------------|
| Qiskit circuit generation from Julia parameters | Generating circuits in Python that use different Hamiltonian rescaling conventions than the Julia code (Julia rescales spectrum to [0, 0.45], Python code in `quantum.py` uses raw parameters) | Export the rescaling factor and shift from Julia's `HamHam` and apply them in the Python Trotter step generator; add a cross-validation test computing the Trotter unitary in both languages |
| PythonCall.jl / juliacall for Julia-Python bridge | Using PyCall.jl (deprecated path) or PyJulia (known JIT startup issues); data format conversion between numpy arrays and Julia arrays silently transposing (row-major vs column-major) | Use PythonCall.jl (actively maintained); explicitly test array layout by round-tripping a known non-symmetric matrix and checking element order |
| BSON serialization for Hamiltonians | Saving/loading `HamHam` structs where the struct definition has changed between versions (added/removed fields) | Version the serialized data; add a `version` field to BSON files; use JSON or HDF5 for cross-language data exchange |
| FINUFFT library (called from Julia) | FINUFFT accuracy parameter `eps=1e-12` may be tighter than needed, causing slow transforms at high qubit counts where matrices are large | Profile FINUFFT time vs. accuracy; `eps=1e-8` may suffice for energy-domain calculations where the Gaussian filter provides natural smoothing |

## Performance Traps

| Trap | Symptoms | Prevention | When It Breaks |
|------|----------|------------|----------------|
| Dense `dim^2 x dim^2` Liouvillian matrix for `run_lindbladian` | Memory usage ~`dim^4 * 16 bytes` = 64 GB for 10 qubits, 16 TB for 12 qubits | Use `LinearMaps` (already imported) for iterative eigensolvers; never materialize the full Liouvillian above 8-9 qubits | At 10 qubits (dim=1024, Liouvillian is 1M x 1M) |
| Per-step full-grid loop in `step_along_trajectory!` iterates over all energy labels | Each trajectory step costs O(num_energy_labels * dim^2); with 2^12 energy labels and dim=256, that is ~268M flops per step | Pre-sort jumps by weight and use early-exit more aggressively; for large systems, consider sparse jump representations | At 8+ qubits with num_energy_bits >= 12 |
| `push!`-based JumpOp vector construction in test files | Array reallocation on each push, type instability if vector not typed | Pre-allocate `Vector{JumpOp}(undef, n)` and fill; annotate types | Minor; only affects setup time, not simulation |
| Shared `KrausScratch` workspace not thread-safe | Silent data races if Distributed/threads are enabled (code has `@distributed` and `@threads` references commented out) | One scratch per thread/worker; use `threadid()` indexing as the commented-out `apply_lindbladian!` partially implements | When re-enabling parallel execution |

## Numerical Validation Pitfalls

### Pitfall: Testing Only Against the Gibbs State Without Checking the Approach Dynamics

**What goes wrong:**
A common mistake is checking only that the final evolved state is close to the Gibbs state. This passes even if the evolution takes a wrong path (e.g., briefly leaving the PSD cone then returning, or thermalizing to a wrong state and then slowly correcting). The Lindbladian's spectral gap determines the convergence rate; if the trajectory converges at a different rate, the unraveling is wrong even if the endpoint happens to match.

**Prevention:**
Compare the trace distance decay curve `d(rho(t), gibbs)` from the DM path against the trajectory-averaged curve. They should match within `O(1/sqrt(ntraj))` at every time point, not just the final one. Also verify that the decay rate matches `exp(-gap * t)` from the spectral analysis.

### Pitfall: Detailed Balance Verification Only at Bohr Level

**What goes wrong:**
The Lindbladian satisfies detailed balance (KMS condition) by construction at the Bohr domain level. But each successive approximation (Energy, Time, Trotter) introduces errors. Testing detailed balance only in BohrDomain gives false confidence that lower domains also satisfy it.

**Prevention:**
For each domain, check that `L(gibbs) ~ 0` (the Gibbs state is approximately stationary). The deviation `||L(gibbs)||` quantifies the approximation error of that domain and should be tracked as a regression metric.

## "Looks Done But Isn't" Checklist

- [ ] **Trajectory framework:** Missing `precompute_R` implementation for `BohrDomain` -- only `EnergyDomain` and `TimeDomain/TrotterDomain` variants exist. Attempting trajectory simulation in BohrDomain will fail at construction.
- [ ] **Coherent term in trajectories:** `with_coherent = true` path has an undefined variable (`trotter` on line 53 of `trajectories.jl`). Currently untested and will crash.
- [ ] **`evolve_and_measure_along_trajectory`:** Commented out (lines 730-786 of `trajectories.jl`). The functionality is partially reimplemented inside `run_trajectories` but the standalone function referenced in `export` does not exist as callable code.
- [ ] **Error computation functions:** `compute_time_oft_quadrature_error`, `compute_trotter_oft_quadrature_error`, `compute_time_B_quadrature_error`, `compute_trotter_B_quadrature_error` are defined as empty functions in `errors.jl` (lines 40-50). They are placeholders.
- [ ] **Package test suite:** Test files use `includet` (Revise-based) and direct `include` instead of `using QuantumFurnace`. They are interactive scripts, not `Pkg.test()`-compatible test suites. No `test/runtests.jl` exists.
- [ ] **LSI optimization:** Uses `d` variable that may not be in scope (should be `dim`). No test validates the gradient computation against finite differences.
- [ ] **Export list:** `precompute_data` and `verify_completeness` are exported but `verify_completeness` is entirely commented out.

## Recovery Strategies

| Pitfall | Recovery Cost | Recovery Steps |
|---------|---------------|----------------|
| Trajectory fallback masking normalization bugs | LOW | Add probability sum assertion; if it fails, audit `base_prefactor` formula in each `step_along_trajectory!` variant against `precompute_R` |
| Cholesky on non-PSD S matrix | MEDIUM | Replace `cholesky!(... check=false)` with eigendecomposition-based square root; zero out negative eigenvalues; add PSD assertion |
| Hardcoded tolerances | LOW | Search-and-replace with named constants; parameterize frequency tolerance by `nu_min` |
| Delta scaling convention mismatch | MEDIUM | Write a channel equivalence test; reconcile the two paths to use the same normalization; document the convention |
| `unsafe_wrap` memory corruption | LOW | Replace with `reshape(view(...))`, fix `d` -> `dim` |
| Coherent term `trotter` bug | LOW | Remove `; trotter=trotter` keyword from `build_trajectoryframework` line 53 |

## Pitfall-to-Phase Mapping

| Pitfall | Prevention Phase | Verification |
|---------|------------------|--------------|
| Trajectory fallback masks bugs | Testing & Validation | Probability sum assertion passes for 1-4 qubit systems across all domains |
| Cholesky non-PSD S | Core Stabilization | `minimum(eigvals(S)) > -1e-14` for all test configurations |
| Hardcoded 1e-12 | Refactoring | Named constants used everywhere; sweep test across qubit counts shows consistent Lindbladians |
| Delta scaling mismatch | Trajectory Validation | DM and trajectory paths produce identical single-step channels on identity/dim |
| `unsafe_wrap` memory safety | Core Stabilization | Gradient finite-difference test passes for LSI; no segfaults |
| Coherent term `trotter` bug | Smoke Testing | `with_coherent=true` trajectory runs without error for all domains |
| No `runtests.jl` | Testing Infrastructure | `Pkg.test("QuantumFurnace")` runs and passes |
| Commented-out code clutter | Refactoring | Source files have < 5% commented-out code |
| Config struct duplication | Refactoring (optional) | Single parametric config type or documented reason for duplication |
| Python convention mismatch (rescaling) | Python Interop | Cross-validation: Julia Trotter unitary == Python Trotter unitary for 4-qubit test case |
| Full Liouvillian materialization at scale | Performance | LinearMaps-based Lindbladian for 10+ qubits; memory stays under 512 GB |

## Sources

- Codebase analysis: `/Users/bence/code/QuantumFurnace.jl/src/trajectories.jl`, `jump_workers.jl`, `log_sobolev.jl`, `coherent.jl`, `furnace.jl`, `constants.jl`, `errors.jl`, `structs.jl`
- [Quantum-Trajectory-Inspired Lindbladian Simulation (PRX Quantum)](https://doi.org/10.1103/ssrs-8x32)
- [Large-scale stochastic simulation of open quantum systems](https://arxiv.org/html/2501.17913)
- [Kraus is king: High-order CPTP low rank method for Lindblad master equation](https://www.sciencedirect.com/science/article/abs/pii/S0021999125003195)
- [PythonCall.jl FAQ & Troubleshooting](https://juliapy.github.io/PythonCall.jl/stable/faq/)
- [Bridging Worlds: Julia and Python Interoperability](https://arxiv.org/html/2404.18170v1)
- [Julia Package Testing Best Practices](https://blog.glcs.io/package-testing)
- [Efficient Quantum Gibbs Sampling with Local Circuits (2025)](https://arxiv.org/abs/2506.04321)

---
*Pitfalls research for: QuantumFurnace.jl -- Lindbladian quantum Gibbs sampling simulation*
*Researched: 2026-02-13*
