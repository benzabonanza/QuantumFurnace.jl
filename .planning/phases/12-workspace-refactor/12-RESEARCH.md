# Phase 12: Workspace Refactor - Research

**Researched:** 2026-02-15
**Domain:** Julia struct refactoring for thread-safe trajectory stepping
**Confidence:** HIGH

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

#### Workspace scope
- Density matrix accumulator lives inside the workspace (not separate)
- `TrajectoryWorkspace{T}` matches framework's type parameterization for type safety
- Claude's discretion on which scratch arrays go into workspace vs stay in framework (goal: framework becomes read-only during stepping)
- Claude's discretion on constructor design (from-framework vs independent)

#### User-facing API
- Workspace is internal -- not exported or exposed to users
- High-level API (thermalize functions) keeps the same external signature; workspace/RNG managed inside
- Trajectory run returns a `TrajectoryResult` struct containing averaged density matrix + step count (minimal -- no forward-looking fields for convergence etc.)
- Result struct extended by later phases as needed

#### Backward compatibility
- Clean break: old function signatures removed, all call sites updated to new signatures
- No deprecation wrappers -- direct migration
- Test logic and assertions stay the same; call sites updated to new signatures
- `TrajectoryFramework` struct modified in place (mutable fields move to workspace) -- no duplicate types

#### RNG contract
- Low-level `step_along_trajectory!` takes `AbstractRNG` as explicit argument
- High-level API takes a seed integer; RNG created internally from seed (Xoshiro default)
- If no seed provided, auto-generate random seed from system entropy
- Seed stored in `TrajectoryResult` -- every run is reproducible after the fact

### Claude's Discretion

- Which scratch arrays go into workspace vs stay in framework (goal: framework becomes read-only during stepping)
- Constructor design (from-framework vs independent)

### Deferred Ideas (OUT OF SCOPE)

None -- discussion stayed within phase scope

</user_constraints>

## Summary

Phase 12 separates the mutable `TrajectoryWorkspace` from the immutable `TrajectoryFramework`, makes `step_along_trajectory!` accept an explicit workspace and `AbstractRNG` argument, and introduces a `TrajectoryResult` struct for returning trajectory outcomes. This is a surgical internal restructuring: the high-level API (`run_trajectories`) keeps the same external signature, but internally allocates workspace and RNG per call.

The current `TrajectoryFramework` struct embeds a mutable `ws::TrajectoryWorkspace{T}` field (line 46 of `trajectories.jl`). The `step_along_trajectory!` function accesses workspace buffers via `fw.ws` (lines 464, 611) and calls the global RNG via bare `rand()` (lines 483, 517, 627, 661). Both patterns block thread-safe execution. The fix is well-scoped: remove `ws` from the framework struct, add `ws::TrajectoryWorkspace` and `rng::AbstractRNG` parameters to the stepping functions, and update all internal and test call sites.

**Primary recommendation:** Move `ws` out of `TrajectoryFramework`, expand `TrajectoryWorkspace` to include `rho_acc` (density matrix accumulator), add `ws` and `rng` parameters to `step_along_trajectory!` and its internal callers, introduce `TrajectoryResult` as a minimal return struct, and update all call sites. The framework struct itself is modified in-place (field removed, not a new type).

## Standard Stack

### Core

No new dependencies required. This phase uses only existing Julia stdlib and project dependencies:

| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| Random (stdlib) | Julia 1.x | `AbstractRNG`, `Xoshiro`, `rand(rng, ...)` | Julia's built-in RNG hierarchy; `Xoshiro` is the default PRNG since Julia 1.7, fast and well-tested |
| LinearAlgebra (stdlib) | Julia 1.x | `mul!`, `dot`, matrix operations in workspace | Already used extensively throughout `trajectories.jl` |

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| `Random.Xoshiro` | `StableRNGs.StableRNG` | StableRNG is version-stable across Julia releases but slower; Xoshiro is faster and sufficient since we store the seed for reproducibility. StableRNG stays in test deps only. |
| Bare `rand()` | `rand(rng, ...)` | Explicit RNG is the correct pattern for reproducible, thread-safe code. Zero overhead difference. |

## Architecture Patterns

### Current State: What Exists

```
TrajectoryFramework{T,D}    (struct, trajectories.jl:29-47)
  domain::D                  # read-only
  jumps::Vector{JumpOp}      # read-only
  ham_or_trott::...          # read-only
  config::...                # read-only
  precomputed_data::Any      # read-only
  per_operator::Vector{...}  # read-only
  n_jumps::Int               # read-only
  delta::Float64             # read-only
  delta_eff::Float64         # read-only
  alpha::Float64             # read-only
  ws::TrajectoryWorkspace{T} # *** MUTABLE -- the problem ***

TrajectoryWorkspace{T}       (struct, trajectories.jl:4-8)
  jump_oft::Matrix{T}        # scratch buffer, mutated every step
  psi_tmp::Vector{T}         # scratch buffer, mutated every step
  Rpsi::Vector{T}            # scratch buffer, mutated every step
```

**Call chain:**
```
run_trajectories(...)                        # high-level API (trajectories.jl:330)
  -> build_trajectoryframework(...)          # builds fw with embedded ws
  -> _evolve_along_trajectory!(psi, fw, t)   # inner loop (trajectories.jl:283)
     -> step_along_trajectory!(psi, fw)      # per-step (trajectories.jl:459, 606)
        -> fw.ws.{jump_oft, psi_tmp, Rpsi}   # accesses mutable workspace
        -> rand(1:fw.n_jumps)                # global RNG
        -> rand() * total_weight             # global RNG
```

### Target State: After Refactoring

```
TrajectoryFramework{T,D}    (struct, modified in-place)
  domain::D
  jumps::Vector{JumpOp}
  ham_or_trott::...
  config::...
  precomputed_data::Any
  per_operator::Vector{...}
  n_jumps::Int
  delta::Float64
  delta_eff::Float64
  alpha::Float64
  # ws field REMOVED

TrajectoryWorkspace{T}       (struct, expanded)
  jump_oft::Matrix{T}        # scratch: OFT buffer
  psi_tmp::Vector{T}         # scratch: matrix-vector result
  Rpsi::Vector{T}            # scratch: R*psi cache
  rho_acc::Matrix{T}         # density matrix accumulator (NEW)

TrajectoryResult{T}          (struct, NEW)
  rho_mean::Matrix{T}        # averaged density matrix
  n_trajectories::Int        # number of trajectories run
  seed::Int                  # RNG seed used (for reproducibility)
```

**New call chain:**
```
run_trajectories(...)                        # same external signature + optional seed kwarg
  -> build_trajectoryframework(...)          # builds fw WITHOUT ws
  -> seed = (user-provided or rand(UInt64))  # capture seed
  -> rng = Xoshiro(seed)                     # create RNG from seed
  -> ws = TrajectoryWorkspace(CT, dim)       # allocate workspace separately
  -> _evolve_along_trajectory!(psi, fw, ws, rng, t)
     -> step_along_trajectory!(psi, fw, ws, rng)
        -> ws.{jump_oft, psi_tmp, Rpsi}     # explicit workspace
        -> rand(rng, 1:fw.n_jumps)           # explicit RNG
        -> rand(rng) * total_weight          # explicit RNG
  -> TrajectoryResult(rho_mean, ntraj, seed) # structured return
```

### Pattern 1: Workspace Allocation -- From-Framework Constructor

**What:** Convenience constructor that creates a workspace sized to match a framework.
**When to use:** When the caller has a framework but needs a fresh workspace.
**Recommendation:** Provide both a from-framework and a standalone constructor.

```julia
# Standalone (type + dimension)
function TrajectoryWorkspace(::Type{T}, dim::Int) where {T}
    TrajectoryWorkspace{T}(
        zeros(T, dim, dim),  # jump_oft
        zeros(T, dim),       # psi_tmp
        zeros(T, dim),       # Rpsi
        zeros(T, dim, dim),  # rho_acc (NEW)
    )
end

# From-framework convenience
function TrajectoryWorkspace(fw::TrajectoryFramework{T}) where {T}
    dim = size(fw.per_operator[1].R, 1)
    TrajectoryWorkspace(T, dim)
end
```

**Rationale:** The standalone constructor is needed for Phase 13 (multi-threading) where each thread creates its own workspace independently. The from-framework constructor is convenient for the single-threaded path and tests.

### Pattern 2: Explicit RNG Threading

**What:** Every `rand()` call in the step function uses an explicit `rng::AbstractRNG` argument.
**When to use:** Always -- this is the new mandatory pattern for all randomness in stepping.

```julia
# BEFORE (current):
a = rand(1:fw.n_jumps)
r = rand() * total_weight

# AFTER (new):
a = rand(rng, 1:fw.n_jumps)
r = rand(rng) * total_weight
```

There are exactly **4 `rand()` call sites** in `step_along_trajectory!`:
- `trajectories.jl:483` -- `rand(1:fw.n_jumps)` (Time/Trotter variant)
- `trajectories.jl:517` -- `rand() * total_weight` (Time/Trotter variant)
- `trajectories.jl:627` -- `rand(1:fw.n_jumps)` (EnergyDomain variant)
- `trajectories.jl:661` -- `rand(rng) * total_weight` (EnergyDomain variant)

### Pattern 3: Seed Capture for Reproducibility

**What:** High-level API accepts optional `seed::Union{Int,Nothing}` kwarg. If `nothing`, generates a random seed from system entropy. The seed is always stored in the result.

```julia
function run_trajectories(...; seed::Union{Int,Nothing}=nothing, ...)
    actual_seed = seed === nothing ? rand(Random.RandomDevice(), UInt64) : UInt64(seed)
    rng = Random.Xoshiro(actual_seed)
    # ... run trajectories with rng ...
    return TrajectoryResult(rho_mean, ntraj, Int(actual_seed))
end
```

**Rationale:** This matches the user's decision: "If no seed provided, auto-generate random seed from system entropy. Seed stored in TrajectoryResult -- every run is reproducible after the fact."

### Anti-Patterns to Avoid

- **Copying the entire framework per thread:** The framework contains large read-only matrices (per-operator Kraus data). Only the workspace needs to be per-thread. Never clone the framework for threading purposes.
- **Using `Random.seed!()` for reproducibility:** `Random.seed!()` modifies the global/task-local RNG state and is not composable. Always pass explicit `rng` objects.
- **Adding new mutable fields to TrajectoryFramework:** The whole point is to make the framework immutable/read-only during stepping. Any new mutable state goes in the workspace.
- **Exporting TrajectoryWorkspace:** The user decision locks this as internal. Only `TrajectoryFramework`, `build_trajectoryframework`, and `step_along_trajectory!` are exported.

## Inventory of Changes

### What Moves Into TrajectoryWorkspace

Current `TrajectoryWorkspace` fields (all stay):
- `jump_oft::Matrix{T}` -- scratch buffer for OFT computation during dissipative jump sampling
- `psi_tmp::Vector{T}` -- scratch buffer for matrix-vector products (K0*psi, U_B*psi, etc.)
- `Rpsi::Vector{T}` -- scratch buffer for R*psi (reused for U_residual*psi)

New field added:
- `rho_acc::Matrix{T}` -- density matrix accumulator (currently allocated separately as `rho_mean` in `run_trajectories`)

**Goal check: Framework becomes read-only during stepping?** YES. After removing `ws`, every field of `TrajectoryFramework` is either immutable (Int, Float64, struct values) or a reference to data that is never mutated during stepping (per_operator matrices, precomputed_data, jumps, config). The framework is fully read-only.

### What Stays in TrajectoryFramework

All current fields EXCEPT `ws`:
- `domain`, `jumps`, `ham_or_trott`, `config`, `precomputed_data` -- input references
- `per_operator` -- precomputed Kraus matrices (read-only during stepping)
- `n_jumps`, `delta`, `delta_eff`, `alpha` -- scalar parameters

### Source File Changes

**`src/trajectories.jl`** -- Primary target:

| Function | Current Signature | New Signature | Change |
|----------|------------------|---------------|--------|
| `TrajectoryWorkspace` struct | 3 fields | 4 fields (+rho_acc) | Add density matrix accumulator |
| `TrajectoryWorkspace(T, dim)` | Constructor | Constructor | Add rho_acc allocation |
| `TrajectoryFramework` struct | Has `ws` field | No `ws` field | Remove field |
| `build_trajectoryframework(...)` | Returns fw with ws | Returns fw without ws | Remove ws creation from builder |
| `step_along_trajectory!(psi, fw)` [Time/Trotter] | 2 args | `step_along_trajectory!(psi, fw, ws, rng)` 4 args | Add ws + rng params, replace `fw.ws` with `ws`, replace `rand()` with `rand(rng, ...)` |
| `step_along_trajectory!(psi, fw)` [EnergyDomain] | 2 args | `step_along_trajectory!(psi, fw, ws, rng)` 4 args | Same changes |
| `_evolve_along_trajectory!(psi, fw, t)` | 3 args | `_evolve_along_trajectory!(psi, fw, ws, rng, t)` 5 args | Pass ws + rng through |
| `_accumulate_density_matrix!(rho, psi)` | Separate rho arg | Use `ws.rho_acc` | Caller passes `ws.rho_acc` |
| `run_trajectories(...)` | Returns NamedTuple | Returns `TrajectoryResult` | Major rework of internal logic |

**`src/QuantumFurnace.jl`** -- Export list:

| Change | Detail |
|--------|--------|
| Remove `TrajectoryWorkspace` from exports (if present) | Currently NOT exported -- no change needed |
| Keep `TrajectoryFramework`, `build_trajectoryframework`, `step_along_trajectory!` exported | Signature changes but names stay |
| Add `TrajectoryResult` to exports | New struct that users receive from `run_trajectories` |

### Test File Changes

All test call sites use the pattern `step_along_trajectory!(psi, fw)` and must be updated to `step_along_trajectory!(psi, fw, ws, rng)`.

| File | Call Sites | What Changes |
|------|------------|--------------|
| `test/test_trajectory_fixes.jl` | 3 calls (lines 26, 53, 120) | Create local `ws` + `rng`, pass to step function |
| `test/test_regression.jl` | 2 calls (lines 82, 146) | Create local `ws` + `rng`, pass to step function |
| `test/trajectory_validation/run_trajectory_validation.jl` | 1 call (line 63) | Create local `ws` + `rng`, pass to step function |
| `test/trajectory_validation/run_convergence_tests.jl` | 2 calls (lines 66, 85) | Create local `ws` + `rng`, pass to step function |

**Test logic and assertions are unchanged** -- only the call site signatures change. Tests that currently use `Random.seed!(42)` before calling step will instead create `rng = Xoshiro(42)` and pass it explicitly. The RNG stream is identical for `Xoshiro(42)` as for `Random.seed!(42)` followed by `rand()` (both use the same Xoshiro algorithm with the same seed).

### run_trajectories Return Value Change

**Current return:** Named tuples of varying shapes depending on code path:
```julia
# Path 1 (ntraj==1, no observables, no store_states):
(framework = fw, psi = psi, rho_mean = rho_mean)

# Path 2 (ntraj>1, no observables):
(framework = fw, states = states, rho_mean = rho_mean)

# Path 3 (with observables):
(framework = fw, times = times, measurements_mean = mean_data, rho_mean = rho_mean)
```

**New return:** `TrajectoryResult` struct (minimal per user decision):
```julia
struct TrajectoryResult{T}
    rho_mean::Matrix{T}        # averaged density matrix (Hermitian, trace 1)
    n_trajectories::Int        # number of trajectories run
    seed::Int                  # RNG seed (for reproducibility)
end
```

The observables/measurements path in `run_trajectories` is preserved but the return value wraps data differently. The `framework` is no longer included in the return (it was only there for debugging; the workspace is now separate). States storage (`store_states`) can be kept as an optional field or removed if not needed -- recommend keeping it as a keyword argument that returns a richer result type or attaches states to the result.

**Decision point (Claude's discretion):** For `run_trajectories` with `observables`, the return must include times and measurement data. Recommend extending `TrajectoryResult` with optional fields rather than returning a different type:

```julia
struct TrajectoryResult{T}
    rho_mean::Matrix{T}
    n_trajectories::Int
    seed::Int
    # Optional measurement data (nothing if no observables)
    times::Union{Nothing, Vector{Float64}}
    measurements_mean::Union{Nothing, Matrix{Float64}}
end
```

This keeps a single return type while accommodating both use cases.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Random number generation | Custom RNG | `Random.Xoshiro(seed)` | Julia's built-in Xoshiro is fast, well-tested, and supports the `AbstractRNG` interface needed for explicit RNG passing |
| System entropy for seed generation | Reading `/dev/urandom` | `rand(Random.RandomDevice(), UInt64)` | Julia's `RandomDevice` is the standard way to get system entropy |
| Density matrix Hermitianization | Manual `(A + A') / 2` | Existing `hermitianize!(A)` helper | Already extracted in Phase 7 (DRY-01) |

## Common Pitfalls

### Pitfall 1: Breaking RNG Stream Compatibility

**What goes wrong:** After switching from `rand()` to `rand(rng, ...)`, the RNG stream may differ from the old global-RNG path even with the same seed. Tests that relied on `Random.seed!(42)` followed by exact numerical comparisons will fail.

**Why it happens:** `Random.seed!(42)` seeds the global `TaskLocalRNG`, while `rng = Xoshiro(42)` creates a new independent Xoshiro instance. Both use the same algorithm, but the initial states are set up differently by `seed!` vs the constructor. In Julia 1.7+, `Random.seed!(42)` and `Random.Xoshiro(42)` produce the same stream, but this must be verified.

**How to avoid:** The trajectory regression tests (`test_regression.jl`) compare averaged density matrices with `atol=0.05` tolerance. This tolerance is large enough to absorb any minor RNG stream differences. The trajectory fix tests (`test_trajectory_fixes.jl`) test structural properties (normalization, CPTP completeness) not exact values. No exact-value RNG tests exist, so this pitfall is low risk.

**Warning signs:** A test that passes with `Random.seed!(42)` but fails with `rng = Xoshiro(42)` indicates an RNG stream difference.

### Pitfall 2: Forgetting to Pass Workspace Through Internal Functions

**What goes wrong:** `_evolve_along_trajectory!` calls `step_along_trajectory!` in a loop. If the workspace is threaded through `step_along_trajectory!` but `_evolve_along_trajectory!` still tries to use `fw.ws`, compilation fails (field no longer exists).

**Why it happens:** The call chain is 3 levels deep (`run_trajectories` -> `_evolve_along_trajectory!` -> `step_along_trajectory!`). It's easy to update the leaf function but forget the intermediate.

**How to avoid:** Start from the bottom (step function), add the parameter, then follow compiler errors upward. Julia's type system will catch any missed `fw.ws` access at the first call.

**Warning signs:** `FieldError: type TrajectoryFramework has no field ws` at compile/call time.

### Pitfall 3: Measurement Path Using Wrong Workspace Buffer

**What goes wrong:** The `run_trajectories` observables path currently reuses `fw.ws.psi_tmp` as a measurement scratch buffer (line 414: `tmp_meas = fw.ws.psi_tmp`). After the refactor, this buffer lives in the workspace. If the measurement path creates its own temporary instead of using `ws.psi_tmp`, allocations increase.

**Why it happens:** The measurement path is a separate code branch in `run_trajectories` that accesses `fw.ws` directly.

**How to avoid:** When updating the observables path, use `ws.psi_tmp` as the measurement scratch buffer. This buffer is safe to reuse between steps because `_accumulate_measurements!` is called after `step_along_trajectory!` completes.

**Warning signs:** Increased allocations in the observables path. The allocation test (`test_allocation.jl`) does not currently cover this path, but new allocation regressions would be visible in benchmarking.

### Pitfall 4: TrajectoryResult Not Capturing Seed When Auto-Generated

**What goes wrong:** If the seed auto-generation happens but the seed value is not stored before creating the RNG, the result contains a stale or wrong seed. The user cannot reproduce the run.

**Why it happens:** The seed capture must happen in the right order: generate seed, store it, create RNG from it, use RNG, return seed in result.

**How to avoid:**
```julia
actual_seed = seed === nothing ? rand(Random.RandomDevice(), UInt64) : UInt64(seed)
rng = Random.Xoshiro(actual_seed)
# ... use rng ...
return TrajectoryResult(..., seed=Int(actual_seed))
```
Always capture the seed in a local variable before creating the RNG.

**Warning signs:** Two runs with `seed=nothing` that should differ produce the same result (seed wasn't regenerated), or a run's result has `seed=0` or some default value.

## Code Examples

Verified patterns from the existing codebase:

### Current step_along_trajectory! Signature (to be replaced)

```julia
# Source: trajectories.jl:459-462
function step_along_trajectory!(
    psi::Vector{<:Complex},
    fw::TrajectoryFramework{<:Complex,D},
) where {D<:Union{TimeDomain,TrotterDomain}}
    ws = fw.ws          # <-- accesses embedded workspace
    # ...
    a = rand(1:fw.n_jumps)  # <-- global RNG
    r = rand() * total_weight  # <-- global RNG
```

### New step_along_trajectory! Signature

```julia
function step_along_trajectory!(
    psi::Vector{<:Complex},
    fw::TrajectoryFramework{<:Complex,D},
    ws::TrajectoryWorkspace{<:Complex},
    rng::AbstractRNG,
) where {D<:Union{TimeDomain,TrotterDomain}}
    # ws passed explicitly (no fw.ws access)
    # ...
    a = rand(rng, 1:fw.n_jumps)  # explicit RNG
    r = rand(rng) * total_weight  # explicit RNG
```

### Workspace Creation in run_trajectories

```julia
# Current (inside build_trajectoryframework):
ws = TrajectoryWorkspace(CT, dim)
return TrajectoryFramework(..., ws)

# New (inside run_trajectories, after build):
fw = build_trajectoryframework(jumps, ham_or_trott, config, precomputed_data, scratch, delta)
ws = TrajectoryWorkspace(CT, dim)  # separate allocation
actual_seed = seed === nothing ? rand(Random.RandomDevice(), UInt64) : UInt64(seed)
rng = Random.Xoshiro(actual_seed)
```

### Test Update Pattern

```julia
# BEFORE:
Random.seed!(42)
step_along_trajectory!(psi, fw)

# AFTER:
ws = TrajectoryWorkspace(fw)
rng = Random.Xoshiro(42)
step_along_trajectory!(psi, fw, ws, rng)
```

### Independence Test (Success Criterion 2)

```julia
# Two independent workspaces stepping from the same framework
ws1 = TrajectoryWorkspace(fw)
ws2 = TrajectoryWorkspace(fw)
rng1 = Xoshiro(100)
rng2 = Xoshiro(200)

psi1 = copy(psi0)
psi2 = copy(psi0)

step_along_trajectory!(psi1, fw, ws1, rng1)
step_along_trajectory!(psi2, fw, ws2, rng2)

# psi1 and psi2 evolved independently -- different RNG seeds produce different outcomes
# ws1 and ws2 contain different scratch data -- no cross-contamination
# fw is unchanged (read-only)
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Workspace embedded in framework | Workspace passed explicitly | This phase | Enables thread-safe execution |
| Global `rand()` in stepping | Explicit `rng::AbstractRNG` | This phase | Enables reproducible parallel execution |
| NamedTuple return from run_trajectories | `TrajectoryResult` struct | This phase | Structured, extensible result type |
| No seed capture | Seed always captured in result | This phase | Every run is reproducible after the fact |

## Open Questions

1. **Observable path complexity**
   - What we know: `run_trajectories` with `observables !== nothing` has a separate code branch that saves measurements at intervals. This path uses `fw.ws.psi_tmp` as a scratch buffer for measurement accumulation.
   - What's unclear: Whether the `TrajectoryResult` struct should include measurement data fields or whether the observable path should return a different/extended result type.
   - Recommendation: Include optional `times` and `measurements_mean` fields in `TrajectoryResult` (set to `nothing` when no observables). This keeps a single return type. Later phases can extend the struct.

2. **`store_states` option in run_trajectories**
   - What we know: There is a `store_states::Bool=false` kwarg that stores per-trajectory final state vectors. This is a debugging feature.
   - What's unclear: Whether this should be preserved in the new API.
   - Recommendation: Keep the kwarg but do NOT add states to `TrajectoryResult`. Instead, when `store_states=true`, return a richer NamedTuple or extended struct. This avoids bloating the minimal result struct. Lower priority -- can be deferred or dropped.

3. **run_thermalization interaction**
   - What we know: `run_thermalization` in `furnace.jl` is the DM-based (density matrix) evolution function. It does NOT use `TrajectoryFramework` or `TrajectoryWorkspace`. It has its own `KrausScratch` workspace and already takes an `rng::AbstractRNG` parameter.
   - What's unclear: Nothing -- this function is unaffected by the workspace refactor.
   - Recommendation: No changes to `run_thermalization`. It already follows the explicit-RNG pattern.

## Sources

### Primary (HIGH confidence)
- **Direct codebase analysis:** `src/trajectories.jl` (740 lines), `src/structs.jl` (312 lines), `src/furnace.jl` (172 lines), `src/QuantumFurnace.jl` (95 lines), `src/kraus.jl` (15 lines), `src/furnace_utensils.jl` (136 lines)
- **Test files:** `test/test_trajectory_fixes.jl` (125 lines), `test/test_regression.jl` (155 lines), `test/test_helpers.jl` (326 lines), `test/test_allocation.jl` (131 lines), `test/runtests.jl` (17 lines)
- **Validation scripts:** `test/trajectory_validation/run_trajectory_validation.jl`, `test/trajectory_validation/run_convergence_tests.jl`
- **Prior milestone research:** `.planning/research/ARCHITECTURE.md`, `.planning/research/STACK.md`, `.planning/research/FEATURES.md`, `.planning/research/PITFALLS.md`
- **Phase context:** `.planning/phases/12-workspace-refactor/12-CONTEXT.md`

### Secondary (MEDIUM confidence)
- Julia Random documentation -- `AbstractRNG`, `Xoshiro`, `RandomDevice` APIs confirmed
- Julia threading documentation -- `TaskLocalRNG` behavior for future Phase 13 compatibility

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- no new dependencies, all Julia stdlib
- Architecture: HIGH -- direct codebase analysis, every affected line identified
- Pitfalls: HIGH -- patterns well-understood from v1.0/v1.1 experience and prior milestone research
- Test impact: HIGH -- all call sites enumerated, test logic unchanged

**Research date:** 2026-02-15
**Valid until:** Indefinite (internal refactoring of a stable codebase; no external dependency drift)
