# Codebase Concerns

**Analysis Date:** 2026-02-13

## Tech Debt

**Unfinished Linear Map Implementation:**
- Issue: `construct_linbdladian_map()` in `src/linearmaps_liouv.jl` is incomplete and uncommented
- Files: `src/linearmaps_liouv.jl` (lines 2-14)
- Impact: Linear map construction for Liouvillian evolution is not usable; large block of test/debug code (lines 16-215) left commented out
- Fix approach: Complete the function signature validation, remove obsolete debug code, or deprecate if no longer needed

**Incomplete Coherent Term Implementation:**
- Issue: `compute_b_minus()` and related functions in `src/coherent.jl` have TODO comment "Reintroduce sigmas here"
- Files: `src/coherent.jl` (line 473)
- Impact: The sigma parameter is not being used in `compute_b_minus()` despite being passed as argument; this may affect correctness of coherent term calculations
- Fix approach: Either remove sigma parameter or properly integrate it into calculation; verify against mathematical derivation

**Incomplete Log-Sobolev Implementation:**
- Issue: `src/log_sobolev.jl` contains TODO at line 17: "Rewrite this with apply_lindbladian!()"
- Files: `src/log_sobolev.jl` (line 17)
- Impact: Function uses outdated approach; refactoring would improve consistency with rest of codebase
- Fix approach: Refactor to use `apply_lindbladian!()` pattern consistent with `src/furnace.jl`

**Threading/BLAS Configuration Not Tested:**
- Issue: TODO at `src/jump_workers.jl` line 640: "test it; set BLAS threads to 1, let julia threads be more"
- Files: `src/jump_workers.jl` (line 640)
- Impact: Threading configuration in computationally heavy loop not validated; potential performance loss or unexpected behavior
- Fix approach: Add performance benchmarks comparing BLAS=1 vs multi-threaded configurations; document optimal settings

## Fragile Areas

**Large Module: jump_workers.jl (976 lines)**
- Files: `src/jump_workers.jl`
- Why fragile: Multiple dispatch implementations for different domains (Bohr, Energy, Time, Trotter) with overlapping logic; contains ~200 lines of old commented-out implementations
- Safe modification: Add comprehensive tests for each domain variant; use integration tests to catch cross-domain regressions
- Test coverage: Existing tests in `test/` appear to focus on individual domains; need cross-domain consistency checks

**Large Module: trajectories.jl (893 lines)**
- Files: `src/trajectories.jl`
- Why fragile: Recent major refactoring (commits 54a1d3c through aa3308e) introduced new `TrajectoryFramework` struct and changed interface; trajectory evolution logic uses in-place mutations
- Safe modification: Add regression tests for trajectory final states vs. Lindbladian evolution; validate against known analytical results for simple systems
- Test coverage: `test/trajectory_test.jl` exists but appears to be timing/benchmark focused, not validation focused

**Unsafe Pointer Operations in log_sobolev.jl:**
- Files: `src/log_sobolev.jl` (lines 44-45, 148-149)
- Why fragile: Uses `unsafe_wrap()` to reinterpret real/imaginary parts of complex arrays; pointer arithmetic with `sizeof(Float64)` offsets
- Safe modification: Add assertions to validate memory layout assumptions; consider using safe reshape alternatives if possible; document pointer arithmetic clearly
- Test coverage: No unit tests for unsafe pointer wrapping; behavior depends on Julia's internal representation

**Recent Trajectory Refactoring (Not Yet Run):**
- Files: `src/trajectories.jl`, `src/coherent.jl`, `src/furnace.jl`
- Why fragile: Commits indicate "trajectory updates without test runs yet" (435d9da) and "Working through trajectory update. Definitely not yet done." (aa3308e); trajectory interface was refactored from `evolve_along_trajectory(psi0, fw, total_time)` to new `TrajectoryFramework` struct
- Safe modification: Run `test/trajectory_test.jl` thoroughly; compare outputs against previous known-good results; validate against Lindbladian predictions
- Test coverage: High risk - recent structural changes not yet verified

## Known Bugs / Issues Under Investigation

**Potential Memory Allocation Issue in Old Code:**
- Issue: Commented FIXME at `src/jump_workers.jl` line 562: "Shouldnt this be out of the loop"
- Files: `src/jump_workers.jl` (lines 558-574 - commented block)
- Trigger: This appears to be from obsolete NUFFT-based Kraus code; unclear if still relevant
- Symptoms: Potential unnecessary allocations inside loop if this code pattern exists elsewhere
- Workaround: Current active code does not use this pattern; verify no similar patterns in production paths

**Trajectory Last-Resort Jump Selection:**
- Issue: Comment at `src/trajectories.jl` line 890: "If somehow we haven't picked any jumps, then use the last one"
- Files: `src/trajectories.jl` (line 890)
- Trigger: Edge case where no jump probabilities exceed random threshold
- Impact: Silent fallback behavior may mask numerical issues or configuration errors
- Recommendation: Make this explicit with warning/error; validate that this path is never hit in normal operation

## Numerical Stability Concerns

**Tolerance Hardcoding:**
- Issue: Multiple hardcoded tolerances (1e-12) scattered throughout codebase without centralized configuration
- Files: `src/coherent.jl` (multiple lines: 367, 371, 385, 386, 387, 392, 406, 407, 408, 413, 427, 428, 429, 434, 461, 485, 501, 508, 510, 513, 515, 518, 520), `src/log_sobolev.jl` (line 52)
- Impact: Inconsistent tolerances across different algorithms; difficult to debug numerical failures; tolerance not appropriate for all problem scales
- Improvement path: Create configuration constants in `src/QuantumFurnace.jl` or config structs; document tolerance rationale for each use case

**Numerical Symmetrization with Epsilon Shift:**
- Issue: Manual symmetrization and epsilon shift at `src/trajectories.jl` lines 75-80
- Files: `src/trajectories.jl` (lines 75-80)
- Impact: Comment states "matches thermalization code logic" but unclear if epsilon value (10 * eps(Float64)) is appropriate; may introduce bias in matrix structure
- Improvement path: Document epsilon choice; validate that shift doesn't affect convergence behavior; consider using Hermitian() constructor instead

**Potential Division by Zero in Normalization:**
- Issue: `src/log_sobolev.jl` line 60: `scale = 1.0 / sqrt(norm_val)` with no zero-check
- Files: `src/log_sobolev.jl` (line 60)
- Trigger: If `norm_val` becomes zero (unlikely but possible with pathological initial conditions)
- Impact: NaN/Inf propagation; silent failure in optimization
- Recommendation: Add assertion or early return; handle zero norm explicitly

**Cholesky Factorization Without Guarantee:**
- Issue: `src/trajectories.jl` line 87: `cholesky!(Hermitian(scratch.tmp2), check=false)` with `check=false`
- Files: `src/trajectories.jl` (line 87)
- Impact: If matrix S is not positive definite (due to numerical errors), Cholesky will fail silently; may result in NaN values in U_residual
- Improvement path: Use `check=true` to catch failures; add fallback using eigendecomposition or SVD if matrix is near-singular

## Performance Bottlenecks

**Commented-Out Old Implementation Not Removed:**
- Issue: `src/jump_workers.jl` contains ~200 lines of commented "Slow and old" implementations (lines 734-948)
- Files: `src/jump_workers.jl` (lines 734-948)
- Impact: Code bloat; confusing for maintenance; no clarity on why old code is retained (historical reference?)
- Improvement path: Either move to separate archive file or remove entirely with git history preserved; document migration path if needed for reference

**Obsolete Caching Struct Still Present:**
- Issue: `OFTCaches` struct in `src/structs.jl` marked as "Became obsolete with NUFFTCaches. But used for debugging."
- Files: `src/structs.jl` (lines 9-21)
- Impact: Maintenance burden; suggests incomplete migration to NUFFT approach
- Fix approach: Either document debugging use cases or remove; verify no production code depends on it

**Quadrature Tolerance Setting Across Integration Calls:**
- Issue: All `quadgk()` calls use fixed `atol=1e-12, rtol=1e-12` without considering problem scale
- Files: `src/coherent.jl` (multiple locations), `src/errors.jl` (line 20)
- Impact: Unnecessarily strict tolerances for some problems (expensive); insufficiently strict for others (inaccurate)
- Improvement path: Make tolerance parameters configurable; scale based on physical parameters (beta, sigma, etc.)

## Security Considerations

**No Input Validation for Quantum State:**
- Issue: State vector normalization assumed but not enforced in trajectory evolution
- Files: `src/trajectories.jl` (line 262 has comment "ensure normalized input (optional, but makes probabilities consistent)")
- Risk: Denormalized input could lead to incorrect jump probabilities; no error thrown
- Recommendation: Add assertion for `norm(psi) ≈ 1.0`; consider adding option to auto-normalize

**Pointer Manipulation Without Bounds Checking:**
- Issue: `unsafe_wrap()` in `src/log_sobolev.jl` assumes correct memory layout without validation
- Files: `src/log_sobolev.jl` (lines 44-45, 148-149)
- Risk: If array layout changes or `sizeof(Float64)` assumption violated, silent memory corruption
- Recommendation: Add docstring explaining memory layout assumptions; consider safer alternatives like `reinterpret`

## Scaling Limits

**Matrix Dimension Scaling for Trajectory Workspace:**
- Issue: Trajectory workspace allocates O(dim²) matrices for every trajectory
- Files: `src/trajectories.jl` (lines 4-16)
- Current capacity: Tested on 4-qubit systems (dim=16); `jump_oft` and related matrices scale as dim²
- Limit: For 10+ qubit systems (dim ≥ 1024), these buffers become memory-intensive; may exceed cache
- Scaling path: Consider sparse representations for jump operators; use FFT-based OFT implementations for larger systems; benchmark memory usage vs. system size

**Bohr Frequency Matrix Scaling:**
- Issue: Bohr frequency calculations create dim×dim matrices (line 39 in `src/trotter_domain.jl`)
- Files: `src/trotter_domain.jl` (line 39), `src/bohr_domain.jl` (various)
- Current capacity: O(dim²) storage for frequency grids
- Limit: For 12+ qubit systems, frequency dictionary/matrix becomes large
- Scaling path: Implement sparse Bohr frequency storage; use approximate frequency binning for large systems

**OFT Precomputation Memory:**
- Issue: Jump operator Fourier transforms precomputed for all energy/time labels
- Files: `src/jump_workers.jl` (precompute_kraus_jumps)
- Current capacity: For 12-bit energy grid, creates ~4000 Kraus operators (one per label per jump per domain state)
- Limit: For 14+ bits, memory becomes prohibitive
- Scaling path: Implement on-the-fly OFT computation with caching; use adaptive grid refinement

## Test Coverage Gaps

**Trajectory Evolution Not Validated Against Lindbladian:**
- What's not tested: Direct comparison between trajectory final state (averaged) vs. Lindbladian steady state
- Files: `test/trajectory_test.jl` (runs @btime benchmark but doesn't validate output)
- Risk: Trajectory updates could silently diverge from correct dynamics; wouldn't be caught by existing tests
- Priority: High - core functionality

**Domain Consistency Not Systematically Tested:**
- What's not tested: Same physical evolution compared across domains (Bohr, Energy, Time, Trotter)
- Files: `src/jump_workers.jl` has 4 domain implementations; tests exist but scattered
- Risk: Domain-specific bugs only manifest when switching domains
- Priority: High - affects correctness

**Coherent Term Validation Missing:**
- What's not tested: Coherent term B against independent calculation or reference
- Files: `src/coherent.jl` (multiple B calculation functions)
- Risk: Todo at line 473 suggests sigma handling incomplete; no test validates correctness
- Priority: Medium - affects accuracy

**Unsafe Pointer Operations Not Tested:**
- What's not tested: `unsafe_wrap()` in log Sobolev code; memory layout assumptions
- Files: `src/log_sobolev.jl` (lines 44-45, 148-149)
- Risk: Silent memory corruption on different Julia versions or configurations
- Priority: Medium - affects reliability

**Large System Scaling Not Tested:**
- What's not tested: Trajectory/Lindbladian evolution for 8+ qubit systems
- Files: `test/trajectory_test.jl` uses 4 qubits; `test/ham_test.jl` varies
- Risk: Memory/performance issues only manifest at scale; numerical precision issues hidden
- Priority: Medium - affects usability for realistic systems

## Missing Critical Features

**No Configuration Validation/Defaults:**
- Problem: Config structures passed to functions without comprehensive validation
- Blocks: Users must know all valid configurations; no helpful error messages for typos/invalid combinations
- Recommendation: Add `validate_config!()` pattern with specific error messages for each config type

**No Convergence Diagnostics for Optimization:**
- Problem: `src/log_sobolev.jl` optimization runs but no convergence monitoring exposed
- Blocks: Users cannot diagnose optimization failures or suboptimal solutions
- Recommendation: Add convergence metrics; expose iteration history; add early stopping criteria

**No API Documentation for Public Functions:**
- Problem: Most functions lack docstrings specifying: preconditions, postconditions, error conditions
- Blocks: External use of library; internal developers must read implementation to understand contracts
- Recommendation: Add comprehensive docstrings to all public functions in `src/QuantumFurnace.jl`

**No Error Recovery or Graceful Degradation:**
- Problem: Numerical failures (NaN, Inf, non-positive-definite matrices) cause silent failures
- Blocks: Users unaware of computation failures; results silently invalid
- Recommendation: Add error checking after critical operations; propagate errors with context

---

*Concerns audit: 2026-02-13*
