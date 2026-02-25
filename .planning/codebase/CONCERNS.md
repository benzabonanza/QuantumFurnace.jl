# Codebase Concerns

**Analysis Date:** 2026-02-25

## Tech Debt

**`linearmaps_liouv.jl` is entirely dead commented-out code:**
- Issue: The entire file is a commented-out `create_liouvillian_map` function plus a large inline test script. The file carries a `#TODO: Finish this with precomputed kraus operators.` header comment that has never been acted on.
- Files: `src/linearmaps_liouv.jl` (214 lines, all commented out)
- Impact: The `LinearMaps` dependency in `Project.toml` may exist solely for this unfinished feature. The file adds noise during navigation and signals unfinished design work.
- Fix approach: Either implement the `LinearMap` wrapper using the current `apply_lindbladian!` infrastructure, or delete the file and remove `LinearMaps` from `Project.toml` if no other code uses it.

**`log_sobolev_manopt.jl` is an empty file:**
- Issue: The file has 1 line (empty). It is included in `QuantumFurnace.jl` at line 125 and contributes nothing.
- Files: `src/log_sobolev_manopt.jl`
- Impact: Misleads readers into expecting a Manopt-based LSI implementation. Module include of an empty file is harmless but signals unfinished work.
- Fix approach: Either implement the planned Manopt variant of `compute_LSI_alpha2`, or remove the file and the corresponding `include` in `src/QuantumFurnace.jl:125`.

**`errors.jl` is a stub placeholder:**
- Issue: Contains only a comment `# Error computation utilities (placeholder for Phase 10 API cleanup)` with no actual code. Still included in the module at `src/QuantumFurnace.jl:112`.
- Files: `src/errors.jl`
- Impact: Zero functional impact, but occupies a module slot and implies API surface cleanup was planned and not completed.
- Fix approach: Populate with actual error helpers (e.g. `compute_trotter_error` could be moved here from `hamiltonian.jl`) or remove the file and its include.

**`#TODO: Rewrite compute_LSI_alpha2 with apply_lindbladian!()`:**
- Issue: `src/log_sobolev.jl:17` has an explicit TODO to rewrite this function using the now-available `apply_lindbladian!` (Phase 27/32 infrastructure). The current implementation builds a dense `dim^2 x dim^2` liouvillian matrix from `LindbladianResult.liouvillian` and passes it directly, which prohibits use with large systems or the matrix-free Krylov path.
- Files: `src/log_sobolev.jl`
- Impact: LSI computation is not usable for any system where building the full Liouvillian is infeasible. Coupling to the dense matrix path also means it lags behind performance improvements in the matvec layer.
- Fix approach: Replace the `L_mat` multiply at `src/log_sobolev.jl:87` with `apply_lindbladian!` / `apply_adjoint_lindbladian!` calls using a `KrylovWorkspace`.

**`#TODO: test it` in `jump_workers.jl`:**
- Issue: A `#TODO: test it; set BLAS threads to 1, let julia threads be more.` comment at `src/jump_workers.jl:464` marks a section of threading interaction code that has never been validated with this specific BLAS threading configuration.
- Files: `src/jump_workers.jl`
- Impact: Potential silently incorrect behaviour when running multi-threaded. The parallel path in `run_thermalization` sets BLAS threads at `src/trajectories.jl:544`, but it is unclear if the jump_workers path respects this consistently.
- Fix approach: Add a targeted threading test similar to `test/test_threading.jl` that exercises `run_thermalization` under multiple Julia threads with BLAS thread count explicitly verified.

**Distributed computing is imported but barely used:**
- Issue: `using Distributed` and `using SharedArrays` are in `src/QuantumFurnace.jl` at lines 12 and 20. The only active usage is in `src/nufft.jl:26,51-53` for `SharedArray` when `nprocs() > 1`. The `@distributed` accumulation loop in `src/furnace.jl:58-60` is commented out.
- Files: `src/QuantumFurnace.jl`, `src/furnace.jl`, `src/nufft.jl`
- Impact: Adds load-time cost for two packages on every import. The multi-process Lindbladian construction is never exercised. `SharedArray` in NUFFT prefactors adds path complexity that is largely untested.
- Fix approach: Remove `using Distributed` and `using SharedArrays` from the main module and guard `SharedArray` use in `nufft.jl` behind a conditional `@require` or dedicated extension. Delete the commented-out `@distributed` loop in `furnace.jl`.

**`OFTCaches` struct is declared obsolete but retained:**
- Issue: `src/structs.jl:346` has comment `# Became obsolete with NUFFTCaches. But used for debugging.` The struct is still exported-adjacent (referenced in `ofts.jl` by `time_oft!` which itself is labelled "Deprecated but used for tests").
- Files: `src/structs.jl`, `src/ofts.jl`
- Impact: Obsolete types pollute the namespace and may confuse new contributors. The `time_oft!` function at `src/ofts.jl:8` depends on it and is only kept for tests.
- Fix approach: Confirm no production path uses `OFTCaches`/`time_oft!`. Remove both and update the tests that depend on them to use the NUFFT path instead.

**Disabled dimension-checking assert in `qi_tools.jl`:**
- Issue: `src/qi_tools.jl:31` has `# @assert size(C) == (m_a * m_b, n_a * n_b)` commented out with the note "can be removed for performance". The assertion guards `kron_add!`, a core linear algebra kernel. Without it, dimension mismatches silently corrupt results.
- Files: `src/qi_tools.jl`
- Impact: Any caller that passes incorrectly sized `C` will get silent memory corruption or wrong results. The performance gain is trivial for the off-hot-path construction code that calls this.
- Fix approach: Restore the `@assert` or gate it behind a `@debug` check.

**`bohr_domain.jl` modified without debug verification:**
- Issue: First line is `#! Changed it slightly for speed without debugging`, indicating a performance-motivated change made without regression testing. The exclamation mark notation is not a standard Julia comment convention.
- Files: `src/bohr_domain.jl`
- Impact: Unknown — the change may be correct, but the comment signals the developer was uncertain and intended to verify later. If the Bohr domain function is incorrect, convergence tests using `BohrDomain` will silently produce wrong results.
- Fix approach: Verify correctness of `B_bohr` against the pre-change reference implementation and add a regression test that pins the output values.

**Config struct fields are heavily `Union{T, Nothing}` with no cross-field enforcement at construction time:**
- Issue: `LiouvConfig`, `ThermalizeConfig`, and GNS variants all declare `a`, `b`, `num_energy_bits`, `t0`, `w0`, `eta`, `num_trotter_steps_per_t0` as `Union{T, Nothing}`. Domain-specific required fields are only validated by `validate_config!` at runtime, not at struct construction.
- Files: `src/structs.jl` (lines 80-87, 107-113, 151-157, 182-188)
- Impact: Users can construct an invalid `ThermalizeConfig{TrotterDomain}` without `t0` or `num_trotter_steps_per_t0` and only receive an error when `validate_config!` is explicitly called (inside `run_trajectories`). Structs created outside the main execution path may carry undetectable invalid state indefinitely.
- Fix approach: Add inner constructors per domain type (e.g. `ThermalizeConfig{TrotterDomain}` requiring non-nothing `t0`, `num_trotter_steps_per_t0`) or use a builder/factory pattern.

**Prefactor computation is duplicated between `jump_workers.jl`, `trajectories.jl`, and `krylov_workspace.jl`:**
- Issue: The `EnergyDomain` prefactor formula `config.w0 / (config.sigma * sqrt(2 * pi)) * gamma_norm_factor` and the Time/Trotter formula `config.w0 * config.t0^2 * (config.sigma * sqrt(2 / pi)) / (2 * pi) * gamma_norm_factor` appear identically in at least three files.
- Files: `src/jump_workers.jl`, `src/trajectories.jl` (lines 186-189), `src/krylov_workspace.jl`, `src/krylov_matvec.jl` (lines 170, 483, 553)
- Impact: Any change to the prefactor physics requires coordinated updates across all occurrences. Prior bugs in this area have caused sign and scaling errors.
- Fix approach: Extract `_compute_prefactor(config, gamma_norm_factor)` as a single canonical function in a shared location (e.g. `misc_tools.jl`) and call it from all sites.

## Known Bugs

**Unreachable code after early return in `energy_domain.jl`:**
- Symptoms: `src/energy_domain.jl:164` has `return energy_labels[start_index:end_index]` that can never execute because the early return at line 162 (`return energy_labels[abs.(energy_labels) .<= sym_limit]`) always fires first.
- Files: `src/energy_domain.jl` (lines 160-164)
- Trigger: Every call to the energy truncation function follows the symmetrized path at line 162. The line 164 return is dead code.
- Workaround: The symmetrized return at line 162 appears to be the intended behaviour. The unreachable line is a copy-paste artifact.
- Fix approach: Delete line 164 (`return energy_labels[start_index:end_index]`).

**`log(sigma_eigen)` applied to an `Eigen` struct instead of values:**
- Symptoms: `src/log_sobolev.jl:38` uses `log(sigma_eigen)` where `sigma_eigen` is an `Eigen` decomposition struct. The intent is to compute `log` of eigenvalues. This relies on Julia's `log` broadcast behaviour over `Eigen` structs, which may not be the element-wise log of `.values`.
- Files: `src/log_sobolev.jl` (line 38)
- Trigger: Any call to `compute_LSI_alpha2`.
- Workaround: Not confirmed as a bug — may work due to Julia's matrix function dispatch on `Eigen` — but the usage pattern is non-obvious and fragile.
- Fix approach: Replace with the explicit `log.(sigma_eigen.values)` to make intent clear and avoid relying on undocumented `Eigen` dispatch.

## Security Considerations

**`unsafe_wrap` and raw `pointer` arithmetic in LSI optimizer:**
- Risk: `src/log_sobolev.jl` lines 44-45 and 148-149 use `unsafe_wrap(Array, pointer(A_flat, offset * sizeof(Float64) + 1), ...)` to partition a flat Float64 vector into real and imaginary views. If `n_params` is computed incorrectly (it is computed as `length(A_flat) / 2`, a Float64 division), or the optimizer reallocates `A_flat`, the pointer is dangling and may corrupt memory.
- Files: `src/log_sobolev.jl` (lines 42-45, 148-149)
- Current mitigation: The division `length(A_flat) / 2` will always be integer-valued given the construction, but Julia does not guarantee this statically.
- Recommendations: Replace `unsafe_wrap` + `pointer` with a safe `view` into the flat vector using integer indexing: `reshape(view(A_flat, 1:dim^2), dim, dim)` and `reshape(view(A_flat, dim^2+1:end), dim, dim)`. The comment at line 43 itself acknowledges this equivalence but dismisses it for BLAS compatibility. Test whether `Array(view(...))` suffices.

## Performance Bottlenecks

**BohrDomain matvec allocates per-matvec:**
- Problem: `src/krylov_matvec.jl:361` allocates `A_nu2_dag = zeros(T, dim, dim)` on every call to `apply_lindbladian!` for `BohrDomain`. The comment says "one allocation per matvec -- acceptable for Bohr". For Krylov iterations this is called thousands of times.
- Files: `src/krylov_matvec.jl` (lines 360-361, 417-418)
- Cause: `KrylovWorkspace` does not include an `A_nu2_dag` scratch buffer for the Bohr two-operator sandwich path.
- Improvement path: Add `bohr_A_nu2_dag::Union{Nothing, Matrix{T}}` to `KrylovWorkspace` and allocate it in the `BohrDomain` constructor branch. Then replace `zeros(T, dim, dim)` with `fill!(ws.bohr_A_nu2_dag, 0)`.

**`apply_lindbladian!` for Time/Trotter recomputes `prefactor` on every call:**
- Problem: `src/krylov_matvec.jl:483` and `553` compute `prefactor = config.w0 * config.t0^2 * (config.sigma * sqrt(2 / pi)) / (2 * pi) * gamma_norm_factor` inside the hot matvec function. This is a config-level constant that never changes between calls.
- Files: `src/krylov_matvec.jl` (lines 483, 553)
- Cause: Prefactor not precomputed in `KrylovWorkspace` construction.
- Improvement path: Compute and cache `scaled_prefactor` in `KrylovWorkspace` at construction time (analogous to `TrajectoryFramework.scaled_prefactor` at `src/trajectories.jl:95`).

**`construct_lindbladian` builds a dense `dim^2 x dim^2` matrix:**
- Problem: `src/furnace.jl` constructs a full Liouvillian matrix. At `n=8` qubits, `dim^2 = 65536`, making the matrix `65536^2 * 16 bytes ≈ 70 GB`. This makes `run_lindbladian` / `run_thermalization` via dense Liouvillian infeasible beyond ~6 qubits.
- Files: `src/furnace.jl`
- Cause: Dense construction is the original design; Krylov path was added later but is separate.
- Improvement path: `krylov_spectral_gap` (via `apply_lindbladian!`) is the scalable path. Add a warning when `dim^2 > 4096` (n > 6) that dense construction will be memory-intensive.

## Fragile Areas

**`TrajectoryFramework.precomputed_data::Any`:**
- Files: `src/trajectories.jl:84`
- Why fragile: The `precomputed_data` field is typed as `Any`. Any access to its fields in the hot `step_along_trajectory!` loop that doesn't go through the concrete-typed cached fields (`scaled_prefactor`, `sigma`, `transition`, `energy_labels`, `oft_nufft_prefactors`) will trigger dynamic dispatch and boxing. The comment at line 95 explicitly warns about this.
- Safe modification: When adding new hot-path computations that need precomputed data, always add a concrete-typed field to `TrajectoryFramework` rather than accessing `precomputed_data` directly in the step loop.
- Test coverage: `test/test_allocation.jl` exists and pins allocations, but coverage depends on the specific test cases included.

**`krylov_spectral_gap` memory warning is advisory only:**
- Files: `src/krylov_eigsolve.jl:60-78`
- Why fragile: The `_check_krylov_memory` function issues `@warn` but does not error. A caller on a memory-constrained machine will proceed to an OOM crash with no graceful handling. The 80% threshold is a heuristic and does not account for other in-flight allocations.
- Safe modification: Always check the memory estimate before calling `krylov_spectral_gap` with large `krylovdim`. Consider making the threshold configurable and adding an `abort_on_oom::Bool` kwarg.
- Test coverage: Not tested for OOM conditions.

**Normalization violation warnings in trajectory step are advisory-only:**
- Files: `src/trajectories.jl:919-920`, `src/trajectories.jl:1057-1058`
- Why fragile: Both `step_along_trajectory!` methods emit `@warn` when `abs(total_weight - 1.0) > 1e-6` but continue execution by sampling from the unnormalized distribution. The trajectory is not aborted or corrected, so accumulated error from repeated violations will silently corrupt the long-time average.
- Safe modification: Track violation counts in `TrajectoryWorkspace` and surface them in `TrajectoryResult`. Consider aborting after a configurable number of consecutive violations above a stricter threshold.
- Test coverage: `test/test_trajectory_fixes.jl` tests some CPTP properties but does not stress-test normalization degradation under edge-case parameters.

**`_select_best_observable` falls back silently to diagnostic-only results:**
- Files: `src/gap_estimation.jl:71-98`
- Why fragile: The tertiary fallback at line 94 selects the observable with highest `R-squared` regardless of whether the fit converged or the gap is physically meaningful. The returned `SpectralGapResult.gap` can be from a non-converged fit with negative `gap` or near-zero amplitude, which looks like a real result to the caller.
- Safe modification: Add a `fit_quality::Symbol` field to `SpectralGapResult` (`:good`, `:fallback1`, `:fallback2`) and document that `:fallback2` results should not be trusted as gap estimates.
- Test coverage: `test/test_gap_estimation.jl` covers the happy path; no tests exercise the fallback 2 case explicitly.

## Scaling Limits

**Dense Liouvillian construction (`construct_lindbladian`):**
- Current capacity: Feasible up to ~6 qubits (dim=64, dim^2=4096, matrix ~128 MB).
- Limit: 8 qubits (dim=256, dim^2=65536) requires ~70 GB just for the matrix; this is beyond typical workstation RAM.
- Scaling path: Use `krylov_spectral_gap` with `apply_lindbladian!` for n>6; the matrix-free path scales to much larger systems constrained only by Krylov subspace size.

**Krylov memory for large systems:**
- Current capacity: For n=8 qubits with default `krylovdim=30`, estimated ~3 GB per solve.
- Limit: n=10 qubits with `krylovdim=30` requires ~50 GB.
- Scaling path: Reduce `krylovdim`, increase `maxiter`, or use restart strategies. The `_eigsolve_with_retry` function increases `krylovdim` on failure which makes OOM worse; retry logic should cap absolute memory not just `krylovdim`.

## Dependencies at Risk

**`Arpack` for spectral gap estimation:**
- Risk: `Arpack.jl` wraps the legacy Fortran ARPACK library. The package is in maintenance mode and has known issues with eigenvalue ordering. The codebase now uses `KrylovKit` as the primary eigensolver (Phases 27-30); `Arpack` is retained for the dense Liouvillian path in `src/gap_estimation.jl` and `src/furnace.jl`.
- Impact: If `Arpack` stalls or produces incorrect eigenvalues for a given system (a documented upstream failure mode), the dense gap estimation path silently returns wrong results.
- Migration plan: Replace `Arpack.eigs` calls with `KrylovKit.eigsolve` using the `apply_lindbladian!` matvec for the dense path as well. `Arpack` can then be removed from `Project.toml`.

**`BSON` for data persistence:**
- Risk: `BSON.jl` serialization is tied to Julia's internal struct layout. Field reordering, type renames, or Julia version changes can break deserialization of existing experiment files. `src/misc_tools.jl:35-57` already contains a hard-coded "legacy HamHam field order (14 fields)" workaround for exactly this class of breakage.
- Impact: Existing `*.bson` files in `hamiltonians/` and `results/` may become unreadable after struct changes without a migration path.
- Migration plan: Adopt a schema-versioned serialization format (e.g. `JLD2.jl` with explicit field names, or `JSON3.jl` for the config/metadata layer). The `_experiment_to_dict` / `_dict_to_experiment` pattern in `src/results.jl` is already a step in this direction.

## Missing Critical Features

**No GPU support:**
- Problem: All matrix operations use CPU BLAS. For large systems (n=8+), GPU-accelerated GEMM could reduce matvec time by 10-100x.
- Blocks: Scaling benchmarks above n=8 within practical wall-clock time.

**`validate_config!` is not called at struct construction:**
- Problem: Invalid configs (e.g. `ThermalizeConfig` for `TrotterDomain` without `t0`) can be constructed without error and passed around. `validate_config!` is only called inside `run_trajectories`, `construct_lindbladian`, and `run_thermalization`. Code that builds configs for later use (e.g. in scripts) can easily miss validation.
- Blocks: Reliable config error reporting; users see errors deep in the call stack.

## Test Coverage Gaps

**`B_bohr` and Bohr domain coherent path:**
- What's not tested: There is no test that verifies the Bohr domain B operator computation against an analytic reference for a non-trivial system. The `#! Changed it slightly for speed` comment in `src/bohr_domain.jl` makes this gap particularly concerning.
- Files: `src/bohr_domain.jl`, `src/coherent.jl`
- Risk: Undetected incorrect B operator silently shifts the steady state away from the Gibbs state.
- Priority: High

**`compute_LSI_alpha2` (log-Sobolev computation):**
- What's not tested: No test file exercises `compute_LSI_alpha2` end-to-end. The function is complex, uses `unsafe_wrap`, and has a known TODO to rewrite it.
- Files: `src/log_sobolev.jl`
- Risk: Silent wrong output or memory corruption during optimization.
- Priority: High

**`linearmaps_liouv.jl` / LinearMap path:**
- What's not tested: The entire file is commented out; no test verifies that a LinearMap-based Liouvillian produces the same spectrum as the dense or Krylov paths.
- Files: `src/linearmaps_liouv.jl`
- Risk: If this path is ever uncommented and shipped, it is completely untested.
- Priority: Low (path is currently unreachable)

**Distributed / SharedArray NUFFT path:**
- What's not tested: The `use_shared_array=true` branch in `src/nufft.jl:51-53` is only active when `nprocs() > 1`, which is never true in the standard test suite.
- Files: `src/nufft.jl`
- Risk: SharedArray-based prefactor computation may produce wrong results or hang without test coverage.
- Priority: Medium

**Gap estimation fallback 2:**
- What's not tested: The tertiary fallback in `_select_best_observable` (highest R-squared regardless of convergence) is not exercised by any test.
- Files: `src/gap_estimation.jl`
- Risk: A bad gap estimate passes silently to downstream consumers.
- Priority: Medium

---

*Concerns audit: 2026-02-25*
