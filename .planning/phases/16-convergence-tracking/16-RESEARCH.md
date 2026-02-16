# Phase 16: Convergence Tracking - Research

**Researched:** 2026-02-16
**Domain:** Batch-level convergence monitoring for quantum trajectory sampling (trace distance, per-observable tracking)
**Confidence:** HIGH

## Summary

Phase 16 adds convergence monitoring to the trajectory sampling engine in QuantumFurnace.jl. The core requirement is: after running batches of trajectories, report how the running average density matrix evolves toward the Gibbs state, and how specific physical observables converge toward their thermal equilibrium values. This gives users visibility into whether a simulation has run long enough, and produces the convergence curves needed for the KMS-vs-GNS comparison paper (Phase 18).

The implementation requires three coordinated additions: (1) a `ConvergenceData` struct that stores scalar convergence metrics (trace distance to Gibbs, per-observable values) at batch checkpoints, (2) a new `run_trajectories_batched` function (or equivalent) that runs trajectory batches, computes the running average density matrix after each batch, measures trace distance and observable values, and stores them in the convergence data, and (3) helper functions that construct the nearest-neighbor Z_iZ_{i+1} and energy <H> observable matrices from the existing Hamiltonian data.

The existing codebase already contains all the mathematical building blocks. `trace_distance_h` in `qi_tools.jl` computes trace distance between Hermitian matrices. `_accumulate_density_matrix!` and `_accumulate_measurements!` in `trajectories.jl` handle density matrix and observable accumulation. The `run_thermalization` function in `furnace.jl` already tracks `distances_to_gibbs::Vector{Float64}` at every step -- Phase 16 adapts this pattern to the trajectory (pure-state unraveling) path with batch-level granularity. The `HamHam` struct contains `gibbs::Hermitian` and `data::Matrix` (the Hamiltonian matrix), both needed for CONV-01 and CONV-03.

**Primary recommendation:** Create a `ConvergenceData` struct holding `(batch_sizes, trace_distances, observable_names, observable_values, observable_gibbs_values)`. Implement `run_trajectories_convergence` that wraps the existing `run_trajectories` infrastructure to run in batches, compute running-average rho after each batch, and fill the convergence data. Use `trace_distance_h(Hermitian(rho_running), gibbs)` for CONV-01. Build Z_iZ_{i+1} and H observable matrices using existing `pad_term` and `HamHam.data`. This phase does NOT implement adaptive stopping (that is Phase 17, CONV-04/CONV-05).

## Standard Stack

### Core

| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| `LinearAlgebra` (stdlib) | Julia 1.11+ | `Hermitian`, `eigvals`, `tr`, `dot`, `mul!` | All matrix operations for trace distance and observable expectation values |
| `Random` (stdlib) | Julia 1.11+ | `Xoshiro(seed)` for per-trajectory RNG | Deterministic batch seeding, same pattern as Phase 13 |
| Existing `qi_tools.jl` | N/A | `trace_distance_h`, `gibbs_state`, `hermitianize!` | Already verified, used in `run_thermalization` |
| Existing `trajectories.jl` | N/A | `run_trajectories`, `_run_chunk_no_obs!`, `_partition_trajectories` | The batch runner wraps these existing functions |
| Existing `misc_tools.jl` | N/A | `pad_term` | Builds observable matrices (Z_iZ_{i+1}) from Pauli operators |

### Supporting

| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| `constants.jl` | N/A | `X`, `Y`, `Z` Pauli matrices | Building Z_iZ_{i+1} observable matrices |
| Existing `results.jl` | N/A | `ExperimentResult`, BSON serialization | ConvergenceData should be serializable alongside ExperimentResult |

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Scalar metrics at checkpoints | Full rho snapshots at checkpoints | Full snapshots enable post-hoc analysis but O(n_batches * dim^2) memory; Pitfall 11 from PITFALLS.md warns against this; scalars are O(n_batches) |
| Trace distance to Gibbs | Fidelity to Gibbs | Trace distance is the standard metric in the Chen papers and is already used by `run_thermalization`; fidelity is an option but adds no value |
| Separate `run_trajectories_convergence` | Adding convergence tracking to existing `run_trajectories` | Separate function avoids API bloat on the existing function; the existing function's observable-tracking path serves a different purpose (time-resolved measurements within a single trajectory evolution, not batch-level convergence) |

## Architecture Patterns

### Current State (After Phase 13/15)

```
run_trajectories(jumps, config, psi0, ham; ntraj, seed, observables, ...)
  -> TrajectoryResult(rho_mean, n_trajectories, seed, times, measurements_mean)

ExperimentResult{C,T}
  config::C
  trajectory_result::TrajectoryResult{Complex{T}}
  hamiltonian_params::Dict{Symbol, Any}
  metadata::Dict{Symbol, Any}
```

The current `run_trajectories` runs ALL ntraj trajectories in one call and returns a single `TrajectoryResult`. There is no concept of batches, checkpoints, or convergence curves.

### Target State (After Phase 16)

```
ConvergenceData
  batch_sizes::Vector{Int}                    # [1000, 1000, 1000, ...]
  cumulative_n_traj::Vector{Int}              # [1000, 2000, 3000, ...]
  trace_distances::Vector{Float64}            # td to Gibbs at each checkpoint
  observable_names::Vector{String}            # ["ZZ_12", "ZZ_23", ..., "H"]
  observable_values::Matrix{Float64}          # n_obs x n_batches
  observable_gibbs_values::Vector{Float64}    # <O_i>_gibbs for reference

run_trajectories_convergence(jumps, config, psi0, ham;
    gibbs, observables, observable_names,
    batch_size, n_batches, seed, ...)
  -> (TrajectoryResult, ConvergenceData)
```

### Pattern 1: Batch-Accumulate-Measure Loop

**What:** Run trajectories in fixed-size batches. After each batch, compute the running average density matrix, measure trace distance and observables, record the scalars.

**When to use:** Always for convergence tracking. This is the core Phase 16 pattern.

**Example:**
```julia
function run_trajectories_convergence(
    jumps, config, psi0, hamiltonian;
    gibbs::Hermitian,
    observables::Vector{<:Matrix{<:Complex}},
    observable_names::Vector{String},
    batch_size::Int = 1000,
    n_batches::Int = 10,
    seed::Union{Int,Nothing} = nothing,
    trotter = nothing,
    total_time = config.mixing_time,
    delta = config.delta,
)
    actual_seed = seed === nothing ? Int(rand(Random.RandomDevice(), UInt64) >> 1) : seed

    CT = eltype(psi0)
    dim = length(psi0)
    rho_acc = zeros(CT, dim, dim)
    n_total = 0

    # Pre-allocate convergence data storage
    trace_dists = Float64[]
    obs_values = Matrix{Float64}(undef, length(observables), 0)
    cum_n_traj = Int[]
    batch_sizes_vec = Int[]

    # Compute Gibbs observable values for reference
    obs_gibbs = [real(tr(Matrix(gibbs) * obs)) for obs in observables]

    for batch_idx in 1:n_batches
        # Seed offset: batch_idx * batch_size ensures non-overlapping seeds
        batch_seed = actual_seed + n_total

        # Run batch using existing run_trajectories
        result = run_trajectories(
            jumps, config, psi0, hamiltonian;
            trotter = trotter,
            total_time = total_time,
            delta = delta,
            ntraj = batch_size,
            seed = batch_seed,
        )

        # Accumulate (running sum, not average)
        rho_acc .+= result.rho_mean .* batch_size
        n_total += batch_size

        # Compute running average
        rho_running = rho_acc ./ n_total
        hermitianize!(rho_running)

        # Measure trace distance
        td = trace_distance_h(Hermitian(rho_running), gibbs)
        push!(trace_dists, td)

        # Measure observable values: tr(rho * O)
        batch_obs = [real(tr(rho_running * obs)) for obs in observables]
        obs_values = hcat(obs_values, batch_obs)

        push!(cum_n_traj, n_total)
        push!(batch_sizes_vec, batch_size)
    end

    rho_final = rho_acc ./ n_total
    hermitianize!(rho_final)

    conv_data = ConvergenceData(
        batch_sizes_vec, cum_n_traj, trace_dists,
        observable_names, obs_values, obs_gibbs,
    )

    traj_result = TrajectoryResult(rho_final, n_total, actual_seed, nothing, nothing)

    return traj_result, conv_data
end
```

### Pattern 2: Observable Matrix Construction

**What:** Build the Z_iZ_{i+1} correlation matrices and the Hamiltonian energy observable from existing infrastructure.

**When to use:** Before calling the convergence tracking function. These matrices are passed as the `observables` argument.

**Example:**
```julia
function build_convergence_observables(hamiltonian::HamHam, num_qubits::Int)
    observables = Matrix{ComplexF64}[]
    names = String[]

    # Nearest-neighbor Z_iZ_{i+1} correlations (periodic)
    for i in 1:num_qubits
        ZZ = Matrix(pad_term([Z, Z], num_qubits, i; periodic=true))
        push!(observables, ComplexF64.(ZZ))
        j = (i % num_qubits) + 1
        push!(names, "ZZ_$(i)$(j)")
    end

    # Energy <H>: the Hamiltonian matrix itself (in computational basis)
    push!(observables, hamiltonian.data)
    push!(names, "H")

    return observables, names
end
```

### Pattern 3: Seed Management Across Batches

**What:** Each batch gets a seed offset so that trajectory IDs never overlap across batches. The first batch uses seeds `actual_seed + 1` through `actual_seed + batch_size`, the second batch uses `actual_seed + batch_size + 1` through `actual_seed + 2*batch_size`, etc.

**When to use:** Always. This ensures the entire batched run is deterministic and reproducible with the same master seed.

**Key insight:** The existing `run_trajectories` seeds each trajectory with `master_seed + traj_id` where `traj_id` ranges from 1 to ntraj. By passing `seed = actual_seed + n_total_so_far` to each batch call, the trajectories across all batches form a contiguous, non-overlapping seed sequence.

### Anti-Patterns to Avoid

- **Storing full rho at every checkpoint:** O(n_batches * dim^2) memory. For n=8 with 100 batches, this is 100 * 256^2 * 16 bytes = 100 MB. Use scalar metrics instead (Pitfall 11 from PITFALLS.md).
- **Computing trace distance with non-Hermitian rho:** Always wrap in `Hermitian()` before calling `trace_distance_h`. Numerical accumulation introduces small non-Hermitian artifacts that cause `eigvals` to return complex values.
- **Using `run_trajectories`' observable path for convergence:** The existing `observables` parameter in `run_trajectories` tracks time-resolved measurements WITHIN each trajectory (recording <O> at intermediate steps during the evolution). This is different from convergence tracking, which monitors how the BATCH-averaged result changes as more trajectories are added. Do not conflate the two.
- **Reusing the same seed across batches:** Each batch must get different seeds. Using the same seed would produce identical trajectories, making the running average never change after the first batch.
- **Computing observable values from individual psi states:** Observable values at batch checkpoints should be computed from `tr(rho_running * O)`, NOT from averaging `<psi|O|psi>` across trajectories. While mathematically equivalent in expectation, the density-matrix-based computation is numerically more stable and directly uses the running average rho that is already being accumulated.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Trace distance computation | Custom matrix norm | `trace_distance_h(Hermitian(rho), gibbs)` from `qi_tools.jl` | Already handles Hermitianization, eigendecomposition, absolute value summation |
| Gibbs state construction | Manual exp(-beta*H)/Z | `hamiltonian.gibbs` field from `HamHam` struct | Pre-computed at HamHam construction time, already in eigenbasis |
| Observable matrix construction | Manual kron products | `pad_term([Z, Z], num_qubits, site)` from `misc_tools.jl` | Handles periodic boundary conditions, sparse construction, correct padding |
| Trajectory batching | Custom trajectory loop | Wrap existing `run_trajectories` with batch_size as ntraj | The existing function handles serial/threaded dispatch, workspace management, BLAS thread control |
| Hermitianization | `(A + A')/2` inline | `hermitianize!(A)` from `qi_tools.jl` | Consistent, in-place, already used throughout codebase |

**Key insight:** Phase 16 does not need any new numerical algorithms. Every computation (trace distance, observable expectation, density matrix averaging, Pauli matrix construction) already exists in the codebase. The new code is orchestration: running batches, computing metrics after each batch, and storing the results in a structured format.

## Common Pitfalls

### Pitfall 1: Observable Matrices Must Be in Computational Basis

**What goes wrong:** The `hamiltonian.data` matrix is in computational basis, but `hamiltonian.gibbs` is in eigenbasis. If observables are constructed in one basis and rho_running is in another, `tr(rho * O)` gives wrong results.

**Why it happens:** `run_trajectories` returns `rho_mean` in the basis of the state vector `psi`. For Hamiltonian eigenbasis simulations, `psi` lives in the eigenbasis. For TrotterDomain, `psi` lives in the Trotter eigenbasis. The observable matrices must match.

**How to avoid:**
- For non-Trotter domains (Energy, Time, Bohr): `psi` is in the Hamiltonian eigenbasis. Observable matrices must also be in the eigenbasis: `O_eigen = V' * O_comp * V` where `V = hamiltonian.eigvecs`.
- For TrotterDomain: `psi` is in the Trotter eigenbasis. Transform observables accordingly.
- The Gibbs state `hamiltonian.gibbs` is already in eigenbasis, so `trace_distance_h(Hermitian(rho_running), gibbs)` works directly when rho_running is also in eigenbasis.

**Warning signs:** Observable values that do not match the Gibbs expectation even after many trajectories. Trace distance that plateaus at a value inconsistent with known convergence behavior.

### Pitfall 2: Eigenbasis Gibbs vs Computational Basis Gibbs

**What goes wrong:** `hamiltonian.gibbs` is diagonal in the eigenbasis (constructed by `_gibbs_in_eigen`). If `rho_running` from `run_trajectories` is in the eigenbasis (which it is for EnergyDomain), the trace distance is computed correctly. But `gibbs_state(hamiltonian, beta)` returns the Gibbs state in the *computational* basis. Using the wrong Gibbs state function for comparison gives incorrect trace distances.

**Why it happens:** Two Gibbs state functions exist: `gibbs_state` (computational basis) and `gibbs_state_in_eigen` (eigenbasis), plus `hamiltonian.gibbs` (eigenbasis, Hermitian-wrapped).

**How to avoid:** Always use `hamiltonian.gibbs` for the reference Gibbs state when comparing with trajectory results. This is what `run_thermalization` uses. The trajectory states evolve in the eigenbasis (or Trotter eigenbasis), and `hamiltonian.gibbs` is in the eigenbasis by construction.

**Warning signs:** Trace distance does not decrease toward zero even with many trajectories and long evolution time.

### Pitfall 3: Running Average Numerical Stability

**What goes wrong:** Accumulating `rho_acc += rho_batch * batch_size` then dividing by `n_total` is mathematically correct but can lose precision if `n_total` becomes very large (>10^6). The absolute values of `rho_acc` grow proportionally to `n_total`, and the final division amplifies any rounding errors.

**Why it happens:** Standard floating-point accumulation error. For n=8 (dim=256) with 100,000 trajectories, each element of `rho_acc` is the sum of ~100,000 complex numbers, each ~1/256 in magnitude. The sum is ~390, and the final average is ~0.004. The relative error from accumulation is ~sqrt(100,000) * eps / 0.004 ~ 1e-11, which is negligible.

**How to avoid:** For the target trajectory counts (up to ~100,000), this is not a practical concern. If ever needed, use compensated summation (Kahan). For Phase 16, standard accumulation is sufficient.

**Warning signs:** Trace distance oscillates or increases at very high trajectory counts (>10^6).

### Pitfall 4: Batch Size Too Small Produces Noisy Convergence Curves

**What goes wrong:** With batch_size=10 and dim=256, each batch's `rho_batch` is extremely noisy (rank at most 10, vs full-rank dim=256). The trace distance curve is jagged and hard to interpret.

**Why it happens:** `rho_batch` for M trajectories has rank at most M. For M << dim, it is far from the true Gibbs state regardless of convergence. The running average improves, but the per-batch noise makes the convergence curve noisy.

**How to avoid:** Use batch_size >= dim for interpretable convergence curves. For n=4 (dim=16), batch_size=100 is fine. For n=8 (dim=256), batch_size=500-1000 is recommended. This also aligns with the adaptive sampling requirements in Phase 17 (CONV-04: batches of reasonable size for convergence detection).

**Warning signs:** Convergence curve is jagged with batch-to-batch oscillations larger than the trend.

### Pitfall 5: TrotterDomain Requires Gibbs State in Trotter Basis

**What goes wrong:** For TrotterDomain, `psi` evolves in the Trotter eigenbasis, not the Hamiltonian eigenbasis. `hamiltonian.gibbs` is in the Hamiltonian eigenbasis. Computing `trace_distance_h(Hermitian(rho_running), hamiltonian.gibbs)` compares matrices in different bases.

**Why it happens:** The Trotter eigenbasis differs from the Hamiltonian eigenbasis by a unitary rotation (`trotter.eigvecs`). The current `run_thermalization` handles this correctly (line 116-118 of furnace.jl): `gibbs_trotter = V_T' * V_H * gibbs_H * V_H' * V_T`.

**How to avoid:** When using TrotterDomain, transform the Gibbs state to the Trotter basis before computing trace distance. Also transform observable matrices to the Trotter basis. This basis-matching must be done at the convergence tracking setup, not inside the hot loop.

**Warning signs:** Trace distance starts at ~0 or does not decrease for TrotterDomain runs.

## Code Examples

### ConvergenceData Struct

```julia
# Source: Codebase analysis of HotAlgorithmResults pattern + FEATURES.md data architecture
"""
    ConvergenceData

Stores convergence metrics at batch checkpoints during trajectory sampling.
Scalars only (no density matrix snapshots) to keep memory O(n_batches).
"""
struct ConvergenceData
    batch_sizes::Vector{Int}
    cumulative_n_traj::Vector{Int}
    trace_distances::Vector{Float64}
    observable_names::Vector{String}
    observable_values::Matrix{Float64}     # n_obs x n_checkpoints
    observable_gibbs_values::Vector{Float64}
end
```

### Building Z_iZ_{i+1} Observables in Eigenbasis

```julia
# Source: test_helpers.jl pattern for jump operator construction + pad_term from misc_tools.jl
function build_zz_observables_eigenbasis(hamiltonian::HamHam, num_qubits::Int)
    V = hamiltonian.eigvecs
    observables = Matrix{ComplexF64}[]
    names = String[]

    for i in 1:num_qubits
        # ZZ in computational basis
        ZZ_comp = Matrix(pad_term([Z, Z], num_qubits, i; periodic=true))
        # Transform to eigenbasis: O_eigen = V' * O_comp * V
        ZZ_eigen = V' * ZZ_comp * V
        push!(observables, ZZ_eigen)
        j = (i % num_qubits) + 1
        push!(names, "ZZ_$(i)$(j)")
    end

    return observables, names
end
```

### Building Energy Observable

```julia
# Source: HamHam struct, hamiltonian.data is in computational basis
function build_energy_observable_eigenbasis(hamiltonian::HamHam)
    # H in eigenbasis is diagonal: diag(eigvals)
    # But we need it as a full matrix for tr(rho * H)
    dim = length(hamiltonian.eigvals)
    H_eigen = zeros(ComplexF64, dim, dim)
    for i in 1:dim
        H_eigen[i, i] = hamiltonian.eigvals[i]
    end
    return H_eigen
end
```

### Computing Gibbs Observable Values

```julia
# Source: qi_tools.jl gibbs_state pattern
function compute_gibbs_observable_values(gibbs::Hermitian, observables::Vector{<:Matrix})
    return [real(tr(Matrix(gibbs) * obs)) for obs in observables]
end
```

### ConvergenceData Dict Serialization (for BSON)

```julia
# Source: results.jl Dict-based serialization pattern
function _convergence_to_dict(conv::ConvergenceData)
    return Dict{Symbol, Any}(
        :batch_sizes           => conv.batch_sizes,
        :cumulative_n_traj     => conv.cumulative_n_traj,
        :trace_distances       => conv.trace_distances,
        :observable_names      => conv.observable_names,
        :observable_values     => conv.observable_values,
        :observable_gibbs_values => conv.observable_gibbs_values,
    )
end

function _dict_to_convergence(d::Dict)
    return ConvergenceData(
        d[:batch_sizes],
        d[:cumulative_n_traj],
        d[:trace_distances],
        d[:observable_names],
        d[:observable_values],
        get(d, :observable_gibbs_values, Float64[]),
    )
end
```

### Trace Distance Computation at Checkpoint

```julia
# Source: qi_tools.jl trace_distance_h + furnace.jl run_thermalization pattern (line 143)
function _checkpoint_trace_distance(rho_running::Matrix{<:Complex}, gibbs::Hermitian)
    hermitianize!(rho_running)
    return trace_distance_h(Hermitian(rho_running), gibbs)
end
```

### TrotterDomain Gibbs Basis Transform

```julia
# Source: furnace.jl run_thermalization (line 116-118)
function gibbs_in_trotter_basis(hamiltonian::HamHam, trotter::TrottTrott)
    return Hermitian(
        trotter.eigvecs' * hamiltonian.eigvecs * hamiltonian.gibbs *
        hamiltonian.eigvecs' * trotter.eigvecs
    )
end
```

## State of the Art

| Old Approach (run_thermalization) | New Approach (Phase 16) | Difference |
|---|---|---|
| Per-step trace distance in DM evolution | Per-batch trace distance in trajectory sampling | DM evolves one rho; trajectory batches average many |psi><psi| |
| `distances_to_gibbs::Vector{Float64}` stored per step | `ConvergenceData` struct with batches, observables, names | Richer data structure supporting per-observable tracking |
| No per-observable tracking in DM path | Z_iZ_{i+1} and <H> tracked at batch checkpoints | Required for paper (CONV-02, CONV-03) |
| Convergence cutoff hardcoded (1e-5) | No adaptive stopping in Phase 16 (deferred to Phase 17) | Clean separation: Phase 16 = monitoring, Phase 17 = adaptive |
| Single TrajectoryResult return | (TrajectoryResult, ConvergenceData) tuple return | Convergence data alongside final result |

**Scope boundary:**
- Phase 16 tracks and records convergence data (CONV-01, CONV-02, CONV-03)
- Phase 17 implements adaptive stopping based on convergence data (CONV-04, CONV-05)
- Phase 15 implemented ExperimentResult + BSON serialization (DATA-01, DATA-02)

## Open Questions

1. **Should ConvergenceData be embedded in ExperimentResult?**
   - What we know: ExperimentResult already has `trajectory_result::TrajectoryResult` and `metadata::Dict`. ConvergenceData could be stored as an additional field or serialized into the metadata Dict.
   - What's unclear: Whether to add a new field to ExperimentResult (breaking its parametric signature) or store convergence data as a separate BSON dict entry.
   - Recommendation: Store ConvergenceData as an optional field in ExperimentResult, OR serialize it as a nested Dict inside the existing metadata. The Dict approach avoids changing ExperimentResult's type signature and is forward-compatible. Use the Dict approach initially.

2. **Should run_trajectories_convergence support multi-threading internally?**
   - What we know: Each batch call to `run_trajectories` already handles multi-threading (Phase 13). The convergence wrapper just calls `run_trajectories` in a loop.
   - What's unclear: Whether the trace distance computation (eigendecomposition of dim x dim matrix) between batches is a bottleneck.
   - Recommendation: Let multi-threading happen inside each `run_trajectories` call (already works). The between-batch bookkeeping (trace distance, observable measurements) is cheap relative to the batch computation. For n=8 (dim=256), `eigvals` on a 256x256 matrix takes ~1ms, negligible vs. batch runtime.

3. **Observable names: strings or symbols?**
   - What we know: The existing metadata Dict uses `Symbol` keys. Observable names are used for display and identification.
   - What's unclear: Whether to use `String` or `Symbol` for observable names.
   - Recommendation: Use `String`. Observable names like "ZZ_12" are display-oriented and may contain characters that are awkward as Symbols. Strings also serialize naturally to BSON.

4. **Should the convergence function also track local magnetization <Z_i>?**
   - What we know: The requirements (CONV-01, CONV-02, CONV-03) specify trace distance, Z_iZ_{i+1}, and <H>. The FEATURES.md analysis also recommends <Z_i> as useful due to disorder breaking SU(2) symmetry. However, the REQUIREMENTS.md explicitly lists "Local magnetization <Z_i> tracking" as out of scope for this milestone.
   - Recommendation: Do NOT include <Z_i> in Phase 16. The observable construction is generic (any Matrix can be passed), so users can add <Z_i> later without code changes. Focus on the three required metrics.

## Sources

### Primary (HIGH confidence)
- **QuantumFurnace.jl codebase** -- Direct analysis of `src/trajectories.jl` (914 lines: run_trajectories, _run_chunk_no_obs!, _accumulate_density_matrix!, _accumulate_measurements!), `src/qi_tools.jl` (230 lines: trace_distance_h, gibbs_state, hermitianize!), `src/furnace.jl` (172 lines: run_thermalization with distances_to_gibbs tracking), `src/structs.jl` (317 lines: TrajectoryResult, HotAlgorithmResults with convergence fields), `src/results.jl` (396 lines: ExperimentResult with Dict serialization), `src/hamiltonian.jl` (332 lines: HamHam with gibbs field), `src/constants.jl` (Pauli matrices), `src/misc_tools.jl` (pad_term for observable construction)
- **Phase 13 RESEARCH.md** -- Multi-threading architecture, per-task workspace, deterministic seeding, BLAS management -- all patterns reused by Phase 16 batch runner
- **REQUIREMENTS.md** -- CONV-01/02/03 requirements definition, explicit out-of-scope for <Z_i>
- **.planning/research/FEATURES.md** -- Observable selection analysis, adaptive sampling protocol, data architecture design, convergence data structure
- **.planning/research/PITFALLS.md** -- Pitfall 11 (convergence tracking memory), Pitfall 12 (batch size effects), Pitfall 7 (premature convergence), Pitfall 5 (sigma comparison), Pitfall 6 (fixed-point vs mixing rate separation)

### Secondary (MEDIUM confidence)
- **test_helpers.jl** -- Test system construction patterns (make_test_system, make_thermalize_config, jump operator construction) -- directly reusable for Phase 16 test setup
- **HotAlgorithmResults pattern** -- The existing `distances_to_gibbs::Vector{Float64}` field demonstrates the project's established pattern for tracking convergence scalars, not density matrix snapshots

### Tertiary (LOW confidence)
- None -- all findings verified against codebase analysis

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- No new dependencies; all functionality exists in codebase (trace_distance_h, pad_term, run_trajectories, HamHam.gibbs)
- Architecture: HIGH -- ConvergenceData struct + batch wrapper is a straightforward orchestration layer over existing proven components
- Pitfalls: HIGH -- All pitfalls identified from codebase analysis and prior research (PITFALLS.md Pitfall 11, 12); basis-matching pitfall identified from furnace.jl TrotterDomain handling
- Code examples: HIGH -- All examples derived from existing codebase patterns (results.jl Dict serialization, test_helpers.jl observable construction, furnace.jl convergence tracking)

**Research date:** 2026-02-16
**Valid until:** 60 days (stable domain; no external dependency drift; all code is internal)
