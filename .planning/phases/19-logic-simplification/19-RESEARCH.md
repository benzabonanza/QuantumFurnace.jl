# Phase 19: Logic Simplification - Research

**Researched:** 2026-02-16
**Domain:** Internal refactoring of Julia codebase (QuantumFurnace.jl)
**Confidence:** HIGH

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

#### Call chain flattening
- Current chain is too deep: `run_experiment() -> run_trajectories_adaptive() -> run_trajectories() -> _evolve_along_trajectory!() -> step_along_trajectory!()` -- 5 levels before actual computation
- Reduce indirection so the path from experiment entry point to trajectory stepping is shorter and more readable
- API signatures, structs, and internal organization can all change freely -- no backward compatibility constraint

#### Jump basis construction
- Eliminate the `jumps_for_diss_raw` / `transform_jumps_to_basis` / `convert` code block that appears in multiple places
- Fix at the source: when constructing `JumpOp`, set `in_eigenbasis` to the correct basis immediately
- For TrotterDomain: use `trotter.eigvecs` for the basis transform when building the JumpOp
- For other domains: use `hamiltonian.eigvecs` as before
- This is safe because TrotterDomain computation is entirely done in Trotter basis
- The downstream `transform_jumps_to_basis` call and `convert(Vector{JumpOp{Matrix{CT}}}, ...)` become unnecessary

#### Result struct simplification
- Current structs (ExperimentResult, TrajectoryResult, ConvergenceData) feel overly complex
- Simplify to distinct result types per simulation method: Lindbladian results, DM simulator results, Trajectory results
- Each type should be clean and self-contained rather than one flexible type with optional fields

### Claude's Discretion
- How to flatten the call chain -- which layers to merge, inline, or keep separate
- Internal architecture of simplified result structs
- Whether `step_along_trajectory!` remains its own function or gets merged
- How to handle ConvergenceData in the new struct hierarchy

### Deferred Ideas (OUT OF SCOPE)
None -- discussion stayed within phase scope
</user_constraints>

## Summary

Phase 19 is a pure internal simplification of three entangled complexities that accumulated during the v1.0-v1.2 development cycle. No new capabilities are added; all existing tests must continue to pass (539 total). The three targets are: (1) a 5-level deep trajectory call chain that makes debugging and comprehension difficult, (2) a redundant jump basis transform pattern duplicated across 6 call sites, and (3) a result struct hierarchy with too many optional fields and overlapping responsibilities.

The research confirms that all three targets are well-scoped and safe to refactor. The jump basis fix is the simplest (user provided exact code). The call chain flattening requires careful analysis of what each layer does to determine the right merge points. The result struct simplification touches serialization code that needs careful migration.

**Primary recommendation:** Execute the three simplifications in dependency order -- jump basis first (enables cleaner call chain code), then call chain flattening (changes how results are produced), then result structs (changes what is returned). Each can be independently verified with the existing test suite.

## Standard Stack

Not applicable -- this is an internal refactoring phase. No new libraries or dependencies are needed. All changes are within the existing Julia codebase.

## Architecture Patterns

### Current Call Chain (Target 1)

The trajectory execution path from experiment entry to actual physics has 5 layers:

```
LAYER 1: run_experiment()                [experiments/run_sweep.jl:84]
  - Builds config, Gibbs state, observables, initial state
  - Calls run_trajectories_adaptive()
  - Wraps result in ExperimentResult, saves to BSON

LAYER 2: run_trajectories_adaptive()     [src/convergence.jl:330]
  - Adaptive batching loop (convergence detection)
  - For each batch, calls run_trajectories()
  - Accumulates rho, computes trace distance per batch
  - Returns (TrajectoryResult, ConvergenceData)

LAYER 3: run_trajectories()              [src/trajectories.jl:472]
  - Builds TrajectoryFramework (one-time precompute)
  - Validates config, allocates workspaces
  - Dispatches serial vs threaded path
  - For each trajectory: calls _evolve_along_trajectory!()
  - Returns TrajectoryResult

LAYER 4: _evolve_along_trajectory!()     [src/trajectories.jl:329]
  - Computes num_steps = ceil(total_time / delta)
  - Normalizes input state
  - Loops: calls step_along_trajectory!() per step

LAYER 5: step_along_trajectory!()        [src/trajectories.jl:648/786]
  - Actual physics: coherent unitary, branching (no-jump/residual/jump)
  - Domain-dispatched (EnergyDomain vs Time/Trotter)
```

**Key inefficiency in the adaptive path:** Layer 2 calls Layer 3 per batch, which rebuilds `TrajectoryFramework` every time (including `_precompute_data`, `build_trajectoryframework`, `validate_config!`, `_print_press`). The framework is immutable and identical across batches -- this is wasted work.

### Recommended Flattened Architecture

**Keep `step_along_trajectory!` as-is** -- it is the hot loop, well-optimized, and its function boundary provides domain dispatch. Merging it would create an unreadable function.

**Merge layers 3+4** -- `_evolve_along_trajectory!` is a trivial loop wrapper (13 lines) that adds no meaningful abstraction. Its body can be inlined into the chunk runners (`_run_chunk_no_obs!`, `_run_chunk_with_obs!`).

**Lift framework building out of `run_trajectories`** -- Create a lower-level function that takes a pre-built `TrajectoryFramework` and just runs N trajectories. The adaptive runner builds the framework once and reuses it across batches.

**Result: 3 effective layers instead of 5:**
```
LAYER 1: run_trajectories_adaptive() / run_trajectories()  [public API]
  - One-time setup: config validation, framework building, workspace allocation
  - Batching/convergence logic (adaptive) or direct run (non-adaptive)

LAYER 2: _run_batch!() or _run_chunk!()  [internal]
  - Given pre-built framework, run N trajectories
  - Step loop inlined from _evolve_along_trajectory!

LAYER 3: step_along_trajectory!()  [internal, hot path]
  - Unchanged: domain-dispatched single-step physics
```

### Current Jump Basis Pattern (Target 2)

The same `transform_jumps_to_basis` pattern appears in **6 locations**:

| Location | File:Line | Pattern |
|----------|-----------|---------|
| Lindbladian construction | `furnace.jl:80-84` | `jumps_for_diss = if TrotterDomain; transform_jumps_to_basis(jumps, trotter.eigvecs)` |
| DM thermalization | `furnace.jl:127-131` | Same pattern |
| Trajectory framework | `trajectories.jl:92-98` | Plus `convert(Vector{JumpOp{Matrix{CT}}}, ...)` |
| Coherent B (total) | `coherent.jl:30` | `trotter_jumps = transform_jumps_to_basis(jumps, ham_or_trott.eigvecs)` |
| Coherent unitaries | `coherent.jl:85` | Same |
| Coherent terms | `coherent.jl:141` | Same |

**What `transform_jumps_to_basis` does** (from `qi_tools.jl:23-25`):
```julia
function transform_jumps_to_basis(jumps::AbstractVector{<:JumpOp}, eigvecs::AbstractMatrix)
    return JumpOp[JumpOp(j.data, eigvecs' * j.data * eigvecs, j.orthogonal, j.hermitian) for j in jumps]
end
```
It replaces `in_eigenbasis` with `eigvecs' * data * eigvecs`. Currently, JumpOps are always constructed with `hamiltonian.eigvecs' * data * hamiltonian.eigvecs` as `in_eigenbasis`, and then **re-transformed** for TrotterDomain.

**The fix:** When building JumpOps for TrotterDomain, set `in_eigenbasis = trotter.eigvecs' * data * trotter.eigvecs` directly. User provided the exact code:
```julia
basis_unitary = (domain isa TrotterDomain) ? trotter.eigvecs : hamiltonian.eigvecs
jump_op_in_eigenbasis = basis_unitary' * jump_op * basis_unitary
```

**IMPORTANT caveat for coherent.jl:** The coherent term functions (`B_trotter`) expect `jump.in_eigenbasis` to be in Trotter basis, and currently get it via `transform_jumps_to_basis`. After the fix, `jump.in_eigenbasis` WILL already be in Trotter basis for TrotterDomain, so the `transform_jumps_to_basis` calls in coherent.jl become unnecessary and must be removed.

However, there is a subtlety in `_precompute_coherent_total_B` (trajectories.jl:106): it passes the ORIGINAL jumps (not `jumps_for_diss`) to the coherent B computation. The comment says "precompute_coherent_total_B handles its own Trotter basis transform internally". After this fix, passing the original jumps would work correctly because their `in_eigenbasis` is already in Trotter basis.

### Current JumpOp Construction Sites

All JumpOp construction sites that need updating:

| Site | File | Current `in_eigenbasis` |
|------|------|------------------------|
| Test helpers (4-qubit) | `test/test_helpers.jl:130-133` | `hamiltonian.eigvecs' * op * hamiltonian.eigvecs` |
| Test helpers (3-qubit) | `test/test_helpers.jl:180-183` | Same |
| Sweep script | `experiments/run_sweep.jl:70-71` | `hamiltonian.eigvecs' * op * hamiltonian.eigvecs` |
| DM scaling test | `test/test_dm_scaling.jl:165,233` | Manual trotter transform for specific tests |
| Old tests/simulations | `simulations/main_*.jl`, `test/old_tests/*.jl` | Various (legacy) |

The test helpers and sweep script are the active sites. The `test_dm_scaling.jl` cases do a manual Trotter transform on top of the energy-basis JumpOp -- these would change because the JumpOp would already be in the correct basis.

### Current Result Struct Hierarchy (Target 3)

```
AbstractConfig
  AbstractLiouvConfig        -> HotSpectralResults (Lindbladian eigendecomposition)
  AbstractThermalizeConfig   -> HotAlgorithmResults (DM step-by-step simulation)
                             -> TrajectoryResult (trajectory averaging)
                             -> ExperimentResult (wraps TrajectoryResult + metadata)
                             -> ConvergenceData (batch convergence metrics)
```

**Current structs and their fields:**

```julia
TrajectoryResult{T}
    rho_mean::Matrix{T}
    n_trajectories::Int
    seed::Int
    times::Union{Nothing, Vector{Float64}}           # nothing if no observables
    measurements_mean::Union{Nothing, Matrix{Float64}} # nothing if no observables

ExperimentResult{C<:AbstractConfig, T<:AbstractFloat}
    config::C
    trajectory_result::TrajectoryResult{Complex{T}}
    hamiltonian_params::Dict{Symbol, Any}
    metadata::Dict{Symbol, Any}

ConvergenceData  (10 fields -- 6 from Phase 16, 4 from Phase 17)
    batch_sizes::Vector{Int}
    cumulative_n_traj::Vector{Int}
    trace_distances::Vector{Float64}
    observable_names::Vector{String}
    observable_values::Matrix{Float64}
    observable_gibbs_values::Vector{Float64}
    converged::Bool
    final_relative_change::Float64
    consecutive_stable_batches::Int
    total_batches::Int

HotSpectralResults{D, T}   (Lindbladian path)
    data::Matrix{Complex{T}}
    fixed_point::Matrix{Complex{T}}
    gap_mode::Matrix{Complex{T}}
    spectral_gap::Complex{T}
    hamiltonian::HamHam{T}
    trotter::Union{TrottTrott{T}, Nothing}
    config::AbstractLiouvConfig{D,T}

HotAlgorithmResults{D, T}   (DM simulation path)
    evolved_dm::Matrix{Complex{T}}
    distances_to_gibbs::Vector{T}
    time_steps::Vector{T}
    hamiltonian::HamHam{T}
    trotter::Union{TrottTrott{T}, Nothing}
    config::AbstractThermalizeConfig{D,T}
```

**Pain points:**
1. `TrajectoryResult` has `Union{Nothing, ...}` fields -- callers must always null-check
2. `ExperimentResult` wraps `TrajectoryResult` + loosely-typed Dict metadata
3. `ConvergenceData` is returned as a separate value alongside `TrajectoryResult` (not embedded)
4. The adaptive/convergence runners return a tuple `(TrajectoryResult, ConvergenceData)` instead of a single coherent result
5. `HotAlgorithmResults` and `HotSpectralResults` carry full `hamiltonian` and `config` -- heavy for no benefit

### Recommended Result Struct Architecture

**Three clean result types, one per simulation method:**

```julia
# Lindbladian eigendecomposition result
struct LindbladianResult{T}
    fixed_point::Matrix{Complex{T}}
    spectral_gap::Complex{T}
    gap_mode::Matrix{Complex{T}}
    liouvillian::Matrix{Complex{T}}  # optional, could be dropped
end

# DM (density matrix) step-by-step simulation result
struct DMSimulationResult{T}
    final_dm::Matrix{Complex{T}}
    trace_distances::Vector{T}
    time_steps::Vector{T}
end

# Trajectory simulation result (the main workhorse)
struct TrajectoryResult{T}
    rho_mean::Matrix{T}
    n_trajectories::Int
    seed::Int
    convergence::Union{Nothing, ConvergenceData}  # present if run with convergence tracking
end
```

**Key changes:**
- Drop `times`/`measurements_mean` from `TrajectoryResult` (observable time series is a separate concern, rarely used)
- Embed `ConvergenceData` inside `TrajectoryResult` instead of returning a tuple
- `ExperimentResult` stays as the persistence wrapper (it has good serialization infra)
- `HotAlgorithmResults` -> `DMSimulationResult` (drops hamiltonian/config baggage)
- `HotSpectralResults` -> `LindbladianResult` (drops hamiltonian/config baggage)

**Serialization impact:** The `_trajectory_to_dict` / `_dict_to_experiment` functions in `results.jl` need updating but the Dict-based BSON pattern remains sound. ConvergenceData already has `_convergence_to_dict` / `_dict_to_convergence`. The merge just nests one inside the other.

**Note on `ConvergenceData` itself:** ConvergenceData with 10 fields is dense but each field has clear purpose and documented provenance (Phase 16 vs Phase 17). Recommend keeping it as-is internally, just embedding it in TrajectoryResult rather than returning separately.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| BSON serialization | Custom binary format | Existing Dict-based BSON pattern | Already proven across 80+ round-trip tests |
| Test validation | New test suite | Existing 539 tests unchanged | The whole point is "same behavior, cleaner code" |

**Key insight:** This is a refactoring phase. Every change must be provably behavior-preserving. The existing test suite is the primary validation tool.

## Common Pitfalls

### Pitfall 1: Breaking the Coherent Term Path
**What goes wrong:** After changing JumpOp construction for TrotterDomain, the coherent.jl functions still call `transform_jumps_to_basis`, double-transforming the jumps.
**Why it happens:** The 6 call sites are spread across 3 files and easy to miss one.
**How to avoid:** Search for ALL `transform_jumps_to_basis` calls. After the fix, this function should have zero callers in the TrotterDomain path. The function itself may still be needed for the `transform_jumps_to_basis` export (public API), but internal callers should be eliminated.
**Warning signs:** TrotterDomain results diverge from pre-refactor baselines while Energy/TimeDomain remain correct.

### Pitfall 2: Framework Rebuilt Per Batch in Adaptive Path
**What goes wrong:** When flattening the call chain, the naive approach is to keep calling `run_trajectories` per batch. This rebuilds the framework every time.
**Why it happens:** `run_trajectories` is a monolithic function that does setup + execution.
**How to avoid:** Split `run_trajectories` into setup (framework building) and execution (running N trajectories with a pre-built framework). The adaptive runner calls setup once and execution per batch.
**Warning signs:** Adaptive runs are slower than expected due to redundant precomputation.

### Pitfall 3: Serialization Forward/Backward Compatibility
**What goes wrong:** Changing result structs breaks loading of previously saved BSON files.
**Why it happens:** BSON files store struct layout, and Dict conversion functions assume specific field names.
**How to avoid:** Keep the Dict-based serialization layer as an adapter. Old files produce old Dicts; new code can handle both via `get()` with defaults.
**Warning signs:** `load_experiment` throws on existing BSON files in experiments/ directory.

### Pitfall 4: Seed Determinism Across Refactored Call Chain
**What goes wrong:** Changing the call chain structure changes which RNG seeds are used for which trajectories, breaking reproducibility guarantees.
**Why it happens:** The current code uses `seed + traj_id` for per-trajectory seeding. If the trajectory ID assignment changes (e.g., different batching), results change.
**How to avoid:** Preserve the exact same `seed + n_total + traj_offset` pattern used in adaptive mode. The trajectory ID is what matters, not the call structure.
**Warning signs:** Determinism tests fail (test_threading.jl, test_convergence.jl).

### Pitfall 5: Removing `times`/`measurements_mean` from TrajectoryResult
**What goes wrong:** Some callers (tests, tutorials) pass `observables` to `run_trajectories` and expect time-resolved measurements back.
**Why it happens:** Observable time series was an early feature, now mostly superseded by convergence tracking.
**How to avoid:** Check all callers of `run_trajectories` with `observables != nothing`. If any are in active use, either keep the feature or migrate them.
**Warning signs:** Tests that pass `observables` to `run_trajectories` fail.

**Current callers with observables:**
- `test/test_workspace_independence.jl` -- no observables (passes nothing)
- `test/test_threading.jl` -- some tests pass observables for measurement verification
- `test/test_dm_scaling.jl` -- no observables
- Active use in `run_trajectories_convergence` / `run_trajectories_adaptive` -- but these use `observables=nothing` when calling `run_trajectories` internally

Grepping reveals `run_trajectories(...; observables=...)` in `test_threading.jl` lines 105-107. This needs to be preserved or migrated.

## Code Examples

### Jump Basis Fix (User-Provided)

**Before** (current pattern in test_helpers.jl):
```julia
jump_in_eigen = hamiltonian.eigvecs' * jump_op * hamiltonian.eigvecs
push!(jumps, JumpOp(jump_op, jump_in_eigen, orthogonal, herm))
```

**After** (domain-aware construction):
```julia
basis_unitary = (domain isa TrotterDomain) ? trotter.eigvecs : hamiltonian.eigvecs
jump_in_eigenbasis = basis_unitary' * jump_op * basis_unitary
push!(jumps, JumpOp(jump_op, jump_in_eigenbasis, orthogonal, herm))
```

### Call Chain Flattening: Lifting Framework Building

**Before** (adaptive calls run_trajectories per batch):
```julia
# convergence.jl -- run_trajectories_adaptive
for batch_idx in 1:max_batches
    result = run_trajectories(jumps, config, psi0, hamiltonian;  # rebuilds framework!
        trotter=trotter, ntraj=batch_size, seed=batch_seed)
    ...
end
```

**After** (framework built once, reused):
```julia
# Shared setup
validate_config!(config)
ham_or_trott = _pick_ham_or_trott(config, hamiltonian, trotter)
precomputed_data = _precompute_data(config, ham_or_trott)
fw = build_trajectoryframework(jumps, ham_or_trott, config, precomputed_data, scratch, delta)

# Batch loop reuses fw
for batch_idx in 1:max_batches
    rho_batch = _run_batch!(fw, psi0, batch_size, batch_seed)  # no rebuild
    ...
end
```

### Result Struct: Embedding ConvergenceData

**Before** (returns tuple):
```julia
function run_trajectories_adaptive(...) -> (TrajectoryResult, ConvergenceData)
```

**After** (single return):
```julia
function run_trajectories_adaptive(...) -> TrajectoryResult
# where TrajectoryResult.convergence::Union{Nothing, ConvergenceData} is populated
```

## State of the Art

Not applicable -- internal refactoring, no external technology changes.

## Open Questions

1. **Should `transform_jumps_to_basis` be removed from the public API?**
   - What we know: After the fix, no internal code will call it for TrotterDomain. It may still be useful for external callers or debugging.
   - What's unclear: Are there external users of this function?
   - Recommendation: Keep as public API but mark as rarely needed. Remove from internal call sites.

2. **What to do with observable time-series in `run_trajectories`?**
   - What we know: The `observables`/`save_every` path exists and is tested (test_threading.jl). Convergence tracking uses batch-level observables instead.
   - What's unclear: Whether any active workflow depends on per-step observable time series.
   - Recommendation: Keep the observable path in `run_trajectories` for now. It costs nothing if not used (the `if observables === nothing` branch). Removing it is a separate cleanup if desired.

3. **Should `run_trajectories_convergence` (fixed-batch) be kept or merged into `run_trajectories_adaptive`?**
   - What we know: `run_trajectories_convergence` is the fixed-batch variant (no stopping criterion). `run_trajectories_adaptive` subsumes it (set n_max = batch_size * n_batches, convergence_threshold = 0).
   - What's unclear: Whether any workflow specifically depends on the fixed-batch API.
   - Recommendation: Keep both for now. Merging is optional and could be a follow-up.

## Sources

### Primary (HIGH confidence)
- Direct codebase inspection of all source files in `src/` (the definitive source of truth for internal refactoring)
- `src/trajectories.jl` -- call chain layers 3-5, framework building, trajectory execution
- `src/convergence.jl` -- call chain layer 2, convergence/adaptive runners
- `src/furnace.jl` -- Lindbladian and DM simulation paths, `jumps_for_diss` pattern
- `src/coherent.jl` -- 3 additional `transform_jumps_to_basis` call sites
- `src/qi_tools.jl` -- `transform_jumps_to_basis` definition
- `src/structs.jl` -- all struct definitions
- `src/results.jl` -- serialization/deserialization code
- `test/test_helpers.jl` -- JumpOp construction patterns
- `experiments/run_sweep.jl` -- experiment entry point, JumpOp construction

### Secondary (MEDIUM confidence)
- Phase 12-18 ROADMAP entries confirming which features are active vs legacy
- STATE.md accumulated decisions documenting design rationale

## Metadata

**Confidence breakdown:**
- Call chain analysis: HIGH -- complete traceable path through source code
- Jump basis transform: HIGH -- all 6 call sites identified, fix validated against user-provided code
- Result struct hierarchy: HIGH -- all struct definitions and serialization code reviewed
- Pitfalls: HIGH -- based on direct code analysis, not speculation

**Research date:** 2026-02-16
**Valid until:** Indefinite (internal codebase analysis, no external dependency versioning concerns)
