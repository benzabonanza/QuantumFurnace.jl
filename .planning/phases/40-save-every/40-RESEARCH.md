# Phase 40: Save Every - Research

**Researched:** 2026-03-01
**Domain:** Observation frequency control for `run_thermalize` trace distance computation
**Confidence:** HIGH

## Summary

Phase 40 adds a `save_every` keyword argument to `run_thermalize` that controls how often the trace distance to the Gibbs state is computed and stored. Currently, trace distance is computed on every single step (line 222 of `furnace.jl`), which involves an eigendecomposition of the difference matrix (`trace_distance_h` calls `eigvals`). For long simulations (e.g., 4-qubit system with `mixing_time=1200`, `delta=0.005` yields 240,000 steps), this per-step eigendecomposition is a significant fraction of total runtime.

The trajectory engine already implements an identical `save_every` pattern for observable measurements: `step % save_every == 0` gating (trajectories.jl lines 483, 529, 718, 835), with `save_every::Int = 1` as the default and `@assert save_every >= 1` for validation. This phase ports the same pattern to the DM thermalization loop, applying it to trace distance computation rather than observable measurement.

The key design decisions are: (1) how to build the `time_steps` and `trace_distances` arrays with the new stride, (2) whether to always record the initial (step 0) and final values regardless of save_every, (3) how the convergence cutoff interacts with save_every (cutoff check can only trigger on steps where trace distance is computed), and (4) the serialization/deserialization implications (ThermalizeResults is a struct, so its field types cannot change).

**Primary recommendation:** Add `save_every::Int = 1` keyword argument to `run_thermalize` with `step % save_every == 0` gating on trace distance computation. Always record the initial state (step 0). Build `time_steps` to match `trace_distances` exactly. Default `save_every=1` preserves identical behavior to current code.

## Standard Stack

### Core

No new libraries required. This is a pure logic change within the existing Julia codebase.

| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| LinearAlgebra (stdlib) | Julia 1.x | `eigvals` inside `trace_distance_h` | Already used, no change |

### Supporting

No new supporting libraries needed.

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| `save_every::Int` (step count) | Time-based saving (save every N seconds) | Step-based is simpler, matches trajectory engine pattern, and is deterministic |
| Always record initial + final | Only record at stride boundaries | Always recording initial/final is more useful for users (can always see start and end state) |
| Modulo check `step % save_every == 0` | Pre-computed save indices set | Modulo is simpler, matches trajectory engine, negligible overhead for integer modulo |

## Architecture Patterns

### Current Architecture (Before This Phase)

```
run_thermalize (furnace.jl, lines 195-243)
  |-- num_steps = ceil(mixing_time / delta)
  |-- trace_distances = [trace_distance_h(initial_dm, gibbs)]    # step 0
  |-- for step in 1:num_steps
  |     |-- [coherent + rho_jump + channel application]
  |     |-- dist = trace_distance_h(evolving_dm, gibbs)          # EVERY step
  |     |-- push!(trace_distances, dist)
  |     |-- @printf("Dist to Gibbs: %s\n", dist)
  |     |-- if dist < convergence_cutoff: break
  |-- time_steps = collect(0.0:delta:(num_steps * delta))
  |-- return ThermalizeResults(config, final_dm, trace_distances, time_steps, metadata)
```

### Target Architecture (After This Phase)

```
run_thermalize (furnace.jl)
  |-- num_steps = ceil(mixing_time / delta)
  |-- trace_distances = [trace_distance_h(initial_dm, gibbs)]    # step 0 (always)
  |-- for step in 1:num_steps
  |     |-- [coherent + rho_jump + channel application]           # EVERY step (physics unchanged)
  |     |-- if step % save_every == 0 || step == num_steps        # GATED
  |     |     |-- dist = trace_distance_h(evolving_dm, gibbs)
  |     |     |-- push!(trace_distances, dist)
  |     |     |-- @printf("Dist to Gibbs: %s\n", dist)
  |     |     |-- if dist < convergence_cutoff: break
  |-- time_steps = [0.0, save_every*delta, 2*save_every*delta, ...]  # matched to trace_distances
  |-- return ThermalizeResults(config, final_dm, trace_distances, time_steps, metadata)
```

### Pattern: Trajectory Engine save_every (Reference Implementation)

The trajectory engine uses the following pattern for save_every gating (trajectories.jl):

```julia
# From run_trajectories (trajectories.jl lines 668-674):
num_steps = ceil(Int, total_time / delta_step)
num_saves = div(num_steps, save_every) + 1   # +1 for initial state
times = Vector{Float64}(undef, num_saves)
@inbounds for s in 1:num_saves
    times[s] = (s - 1) * save_every * delta_step
end

# From inner loop (trajectories.jl lines 715-721):
save_idx = 1
for step in 1:num_steps
    step_along_trajectory!(psi, ws, rng)
    if step % save_every == 0
        save_idx += 1
        _accumulate_measurements!(mean_data, save_idx, psi, observables, tmp_meas)
    end
end
```

Key features:
- `+1` in `num_saves` accounts for the initial state (step 0)
- `save_idx` starts at 1 (initial state already stored at index 1)
- Only steps where `step % save_every == 0` trigger computation
- `times` array is pre-computed to match saved indices

### Design Decision: Always Record Final Step

The trajectory engine does NOT force-record the final step because it uses fixed-size pre-allocated arrays. For `run_thermalize`, we use dynamic `push!` so we can easily force-record the final step if it does not fall on a `save_every` boundary. However, the success criteria says "lengths match ceil(n_steps / save_every) (plus or minus 1 for initial/final step handling)", so we should keep the behavior simple and predictable.

**Recommended approach:** Record at step 0 (initial), then at every step where `step % save_every == 0`. Also always record the last step if the simulation ends early due to convergence cutoff (to capture the final state). Do NOT force-record the absolute last step of a full run that does not coincide with a `save_every` boundary -- this keeps array lengths deterministic: `1 + div(num_steps, save_every)` entries (or fewer if convergence cutoff triggers).

Actually, on reflection, recording the final step is important for usability. Users want to see the last trace distance. The simplest approach:
- Always record step 0
- Record on `step % save_every == 0`
- If the loop completes (no early convergence) and the final step was NOT on a save_every boundary, also record it
- If convergence cutoff triggers, that step is already being checked (it is a save_every boundary), so it gets recorded naturally

But this makes the array length less predictable. The cleanest approach, matching the trajectory engine, is:
- Always record step 0
- Record on `step % save_every == 0`
- Convergence cutoff only checked at save points
- No special final-step recording

This gives deterministic array lengths: `1 + div(actual_steps_taken, save_every)`.

### ThermalizeResults Struct (No Change Needed)

```julia
struct ThermalizeResults{T<:AbstractFloat} <: AbstractResults
    config::Config
    final_dm::Matrix{Complex{T}}
    trace_distances::Vector{T}
    time_steps::Vector{T}
    metadata::Dict{Symbol, Any}
end
```

The struct does not need any changes. `trace_distances` and `time_steps` are already `Vector{T}` -- they just get shorter when `save_every > 1`. The `save_every` value can be stored in `metadata` for provenance.

### Anti-Patterns to Avoid

- **Changing the physics:** The density matrix evolution must happen on EVERY step regardless of `save_every`. Only the observation (trace distance computation and printing) is gated. If someone accidentally gates the DM evolution, the simulation becomes wrong.
- **Breaking backward compatibility:** Default `save_every=1` must produce bit-identical results to the current code. The `trace_distances` and `time_steps` arrays must match exactly.
- **Off-by-one in time_steps:** The current code uses `collect(0.0:delta:(num_steps * delta))`, which always has `num_steps + 1` entries (matching `length(trace_distances)` because trace_distances starts with the initial state). With `save_every > 1`, time_steps must still match trace_distances element-by-element.
- **Convergence cutoff only on save points:** The convergence check currently happens every step. With `save_every > 1`, it can only trigger on save points. This is acceptable -- the user chose to save less frequently, so they accept coarser convergence detection.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Save frequency gating | Custom iterator or index set | `step % save_every == 0` | Simple, proven pattern from trajectory engine, zero allocation |
| Time array construction | Manual accumulation | `collect(0.0:(save_every * delta):(num_steps * delta))` or equivalent | Matches stride to save frequency |

**Key insight:** This is a simple conditional gating of existing logic. There is no complex machinery to build.

## Common Pitfalls

### Pitfall 1: Off-by-One in time_steps vs trace_distances Length

**What goes wrong:** `time_steps` and `trace_distances` have different lengths, causing downstream plotting or analysis code to crash.
**Why it happens:** `trace_distances` includes the initial state (step 0) as its first element, plus one element per save point. If `time_steps` is constructed differently, lengths diverge.
**How to avoid:** Build `time_steps` to have exactly `length(trace_distances)` entries. The initial entry is `0.0`, then each subsequent entry corresponds to the step where a save occurred. Use the same counting logic for both arrays.
**Warning signs:** `@assert length(time_steps) == length(trace_distances)` fails.

### Pitfall 2: Early Convergence Cutoff Misalignment

**What goes wrong:** Simulation converges between save points, but convergence is not detected until the next save point. Or worse, the `num_steps` variable is mutated after break but time_steps construction uses the original value.
**Why it happens:** Currently `num_steps = step` is set before `break`. With `save_every > 1`, the last recorded step number differs from the total step count.
**How to avoid:** Track the number of recorded save points separately from the total step count. Use the recorded save points to construct time_steps. Better: build time_steps from the actual indices where saves occurred, or from `length(trace_distances)`.
**Warning signs:** time_steps has entries beyond the actual simulation duration.

### Pitfall 3: Printf Noise with Large save_every

**What goes wrong:** With `save_every=1`, every step prints. With `save_every=100`, only every 100th step prints. This is the desired behavior, but the `@printf` is currently tied to the trace distance computation. If someone separates them, they might print without computing or compute without printing.
**Why it happens:** The `@printf` and `push!` are currently independent statements after `trace_distance_h`. They should remain coupled inside the same `if` block.
**How to avoid:** Keep the `@printf`, `push!`, and convergence check inside the same `if step % save_every == 0` block.
**Warning signs:** Inconsistent printing frequency vs actual save frequency.

### Pitfall 4: Metadata Does Not Record save_every

**What goes wrong:** A ThermalizeResults loaded from BSON has `save_every=10` trace distances, but the user does not know what save_every was used. They assume every-step data and misinterpret results.
**Why it happens:** `save_every` is not stored in config or metadata.
**How to avoid:** Store `save_every` in the metadata dict: `metadata[:save_every] = save_every`.
**Warning signs:** User confusion about time resolution of trace distance data.

### Pitfall 5: Serialization Backward Compatibility

**What goes wrong:** Old ThermalizeResults BSON files (with every-step data) fail to load because the loader expects a `save_every` field that does not exist.
**Why it happens:** Adding `save_every` to the struct or to a required field in the serialization format.
**How to avoid:** Do NOT change the ThermalizeResults struct. Store `save_every` only in the metadata dict (which uses `get` with defaults for backward compatibility). Old BSON files will simply not have `metadata[:save_every]`, which can default to 1.
**Warning signs:** `load_result` crashes on pre-Phase-40 BSON files.

## Code Examples

Verified patterns from the existing codebase:

### Trajectory Engine save_every Pattern (trajectories.jl)

```julia
# Source: trajectories.jl lines 641-646 (run_trajectories signature)
function run_trajectories(jumps, config, psi0, hamiltonian;
    ...,
    save_every::Int = 1,
    ...
)
    @assert save_every >= 1
    # ...
    num_saves = div(num_steps, save_every) + 1
    times = Vector{Float64}(undef, num_saves)
    @inbounds for s in 1:num_saves
        times[s] = (s - 1) * save_every * delta_step
    end
end
```

### Current run_thermalize Hot Loop (furnace.jl lines 195-231)

```julia
# Source: furnace.jl lines 195-231 (current implementation)
num_steps = Int(ceil(config.mixing_time / config.delta))
convergence_cutoff = 1e-5
trace_distances = [trace_distance_h(Hermitian(evolving_dm), gibbs)]

for step in 1:num_steps
    # ... physics (coherent, rho_jump, channel) ...

    dist = trace_distance_h(Hermitian(evolving_dm), gibbs)
    push!(trace_distances, dist)
    @printf("Dist to Gibbs: %s\n", dist)
    if dist < convergence_cutoff
        num_steps = step
        break
    end
end

time_steps = collect(0.0:config.delta:(num_steps * config.delta))
```

### Proposed Refactored Hot Loop (New)

```julia
# Proposed implementation for run_thermalize with save_every:
@assert save_every >= 1 "save_every must be >= 1"

num_steps = Int(ceil(config.mixing_time / config.delta))
convergence_cutoff = 1e-5
trace_distances = [trace_distance_h(Hermitian(evolving_dm), gibbs)]
recorded_steps = Int[0]  # track which steps were recorded

for step in 1:num_steps
    idx = rand(rng, 1:length(jumps))
    jump = jumps[idx]

    # Physics (ALWAYS runs, regardless of save_every)
    _apply_coherent_unitary!(
        evolving_dm,
        coherent_unitaries === nothing ? nothing : coherent_unitaries[idx],
        scratch,
    )
    _accumulate_rho_jump!(
        scratch, evolving_dm, jump, ham_or_trott, config, precomputed_data;
        jump_weight_scaling = jump_weight_scaling,
    )
    _apply_precomputed_channel!(
        evolving_dm, K0s[idx], U_residuals[idx], scratch,
    )

    # Observation (GATED by save_every)
    if step % save_every == 0
        dist = trace_distance_h(Hermitian(evolving_dm), gibbs)
        push!(trace_distances, dist)
        push!(recorded_steps, step)
        @printf("Dist to Gibbs: %s\n", dist)
        if dist < convergence_cutoff
            break
        end
    end
end

time_steps = T.(recorded_steps .* config.delta)
```

### Alternative: Simpler time_steps Construction

```julia
# If we don't track recorded_steps explicitly, we can compute time_steps from
# the length of trace_distances and save_every:
n_saves = length(trace_distances)  # includes initial state
time_steps = T[i * save_every * config.delta for i in 0:(n_saves - 1)]
```

This is cleaner but only works if no early convergence truncation changes the stride pattern. Since convergence can only trigger on save_every boundaries, the stride IS regular, so this works.

However, the `recorded_steps` approach is more robust and explicit. Both are acceptable.

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Compute trace distance every step | Still every step (Phase 40 adds control) | This phase (40) | Eliminates unnecessary eigendecompositions for long runs |
| No save_every in DM thermalization | Trajectory engine has save_every | Phase 22 (v1.3) | DM path now gets the same feature |

**Already implemented in trajectory engine:**
- `save_every::Int = 1` keyword with `@assert save_every >= 1`
- `step % save_every == 0` gating
- `num_saves = div(num_steps, save_every) + 1` array sizing
- `times[s] = (s - 1) * save_every * delta_step` time construction

## Implementation Strategy

### Phase Decomposition

This is a small, focused change. A single plan is likely sufficient:

**Plan 40-01: Add save_every to run_thermalize**
1. Add `save_every::Int = 1` keyword argument to `run_thermalize` signature
2. Add `@assert save_every >= 1` validation
3. Gate the trace distance computation, push, printf, and convergence check inside `if step % save_every == 0`
4. Build `time_steps` to match `trace_distances` (accounting for save_every stride)
5. Store `save_every` in metadata
6. Update docstring
7. Tests:
   - `save_every=1` produces identical results to current behavior (backward compatibility)
   - `save_every=10` produces correct trace_distances and time_steps lengths
   - `time_steps` and `trace_distances` have same length
   - `time_steps[2] - time_steps[1] == save_every * delta` (correct stride)
   - Convergence cutoff works correctly with save_every > 1
   - Serialization round-trip works (ThermalizeResults with shorter arrays)

### Touched Files

| File | Change | Scope |
|------|--------|-------|
| `src/furnace.jl` | Add `save_every` kwarg, gate trace distance computation | ~15 lines modified |
| `test/test_results.jl` | Add test for save_every metadata in round-trip | ~10 lines added |
| New test file or section | Behavioral tests for save_every | ~40-60 lines added |

### What Does NOT Change

- `ThermalizeResults` struct (unchanged)
- `_thermalize_to_dict` / `_dict_to_thermalize_results` (unchanged -- arrays are just shorter)
- `_write_result_companion_txt` (unchanged -- already handles variable-length arrays)
- `Config` struct (save_every is a runtime parameter, not a config parameter)
- Physics/evolution code (coherent, rho_jump, channel application)
- Convergence cutoff value (still 1e-5)
- `save_result` / `load_result` (unchanged)

## Open Questions

1. **Should save_every be stored in Config or as a keyword argument?**
   - What we know: The trajectory engine uses it as a keyword argument, not a Config field. Config holds physics parameters (beta, sigma, delta, mixing_time), not observation parameters.
   - What's unclear: Whether anyone would want save_every to be serialized as part of the config for reproducibility.
   - Recommendation: Keyword argument (not in Config), stored in metadata dict. This matches trajectory engine convention and avoids changing the Config struct.

2. **Should we always record the final step even if it does not fall on a save_every boundary?**
   - What we know: The trajectory engine does NOT force-record the final step. Success criteria says "ceil(n_steps / save_every) (plus or minus 1 for initial/final step handling)".
   - What's unclear: Whether users expect to always see the final trace distance.
   - Recommendation: Do NOT force-record the final step. Keep behavior simple and deterministic: initial + every save_every-th step. The final DM is always returned in `final_dm` regardless. If the user needs the final trace distance, they can compute it from `final_dm` and the Gibbs state. This also makes array length exactly `1 + div(actual_steps, save_every)` (or fewer if convergence cutoff triggers within a save window).

3. **Should the convergence cutoff check happen on non-save steps?**
   - What we know: Currently the cutoff check is coupled to trace distance computation. Computing trace distance just for the cutoff check would negate the performance benefit.
   - What's unclear: Whether users expect convergence detection to be as fine-grained as before.
   - Recommendation: Only check convergence on save steps. This is the natural consequence -- convergence detection has `save_every`-step granularity. Document this behavior.

## Sources

### Primary (HIGH confidence)
- `src/furnace.jl` lines 120-243 -- `run_thermalize` implementation (current state post-Phase 39)
- `src/trajectories.jl` lines 618-723 -- `run_trajectories` with `save_every` implementation
- `src/trajectories.jl` lines 734-838 -- `run_observable_trajectories` with `save_every`
- `src/trajectories.jl` lines 1210-1308 -- `run_trajectory` unified API with `save_every`
- `src/structs.jl` lines 200-207 -- `ThermalizeResults` struct definition
- `src/results.jl` lines 274-283 -- `_thermalize_to_dict` serialization
- `src/results.jl` lines 356-365 -- `_dict_to_thermalize_results` deserialization
- `src/qi_tools.jl` lines 132-136 -- `trace_distance_h` implementation
- `test/test_results.jl` -- ThermalizeResults serialization tests
- `test/test_helpers.jl` -- `make_config`, `make_test_system` test factories

### Secondary (MEDIUM confidence)
- `.planning/REQUIREMENTS.md` -- SAVE-01, SAVE-02 requirement definitions
- `.planning/ROADMAP.md` -- Phase 40 placement and dependencies

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- no new libraries, pure logic change
- Architecture: HIGH -- pattern is directly copied from trajectory engine which is production-validated
- Pitfalls: HIGH -- all identified from direct code reading; off-by-one and backward compatibility concerns are well-understood
- Implementation: HIGH -- small, focused change with clear reference implementation

**Research date:** 2026-03-01
**Valid until:** 2026-04-01 (stable -- internal logic change, no external dependencies)
