# Phase 36: API and Results - Context

**Gathered:** 2026-02-27
**Status:** Ready for planning

<domain>
## Phase Boundary

Define 4 clean public entry points (`run_lindblad`, `run_thermalize`, `run_krylov_spectrum`, `run_trajectory`) each returning a typed Result struct with optional BSON save capability. Consolidate the current 7+ entry points and inconsistent result types into a uniform API. Replace `ExperimentResult` with 4 typed Result structs.

</domain>

<decisions>
## Implementation Decisions

### Entry point signatures
- Keep separate positional args: `run_*(jumps, config, hamiltonian; kwargs...)` — Config stays focused on simulation parameters, not problem definition
- Hamiltonian and trotter are both positional args at the same level (both are eigenbasis objects, used interchangeably as `ham_or_trott`): `run_*(jumps, config, hamiltonian, trotter; kwargs...)`
- Trotter defaults to `nothing` when not used
- Trajectory-specific kwargs (ntraj, total_time, delta, seed, save_every, observables) stay as flat keyword args, no options struct

### Trajectory API consolidation
- All 4 current trajectory functions (`run_trajectories`, `run_observable_trajectories`, `run_trajectories_convergence`, `run_trajectories_adaptive`) collapse into a single `run_trajectory`
- `convergence=true` keyword triggers convergence tracking; Gibbs reference comes from `hamiltonian.gibbs` (already available in HamHam struct)
- `adaptive=true` keyword (requires `convergence=true`) enables adaptive early stopping with additional kwargs (convergence_threshold, patience, etc.)
- DM reconstruction logic: no observables = reconstruct DM at save points (mid-run); with observables = track observables mid-run, reconstruct DM once at the end as final average

### Result struct contents
- 4 typed Result structs: `LindbladResults`, `ThermalizeResults`, `KrylovSpectrumResults`, `TrajectoryResults`
- `ExperimentResult` is removed entirely — replaced by the 4 typed structs
- `AbstractResults` as common supertype if needed for dispatch
- Each Result stores the full `Config{S,D,C,T}` for reproducibility
- Metadata in all results: git hash, timestamp, wall_time_seconds, thread count (no Julia version)
- `LindbladResults`: spectral data only (eigenvalues, fixed_point, gap_mode, spectral_gap) — no full Liouvillian matrix
- `ThermalizeResults`: includes trace distances over time for convergence-to-Gibbs plotting
- `TrajectoryResults`: convergence data and observable time series are `Union{Nothing, ...}` — only populated when those features were used

### Save/load workflow
- Standalone `save_result(result, path)` function — not a keyword on run_* functions
- Every `save_result` call always creates a companion `.txt` file with: date, config summary (n_qubits, beta, domain, construction), key results
- `load_result(path)` auto-detects the result type from a tag stored in the BSON — returns the correct typed Result
- Clean break with old `ExperimentResult` .bson files — no backward compatibility needed

### Claude's Discretion
- Internal serialization format (Dict-based vs direct struct BSON)
- Companion .txt content formatting and level of detail per result type
- Whether `AbstractResults` is actually needed or if duck typing suffices
- Auto-filename generation patterns (if any)
- How metadata is captured internally (_capture_metadata helper pattern)

</decisions>

<specifics>
## Specific Ideas

- Hamiltonian and trotter are both eigenbasis objects — in the codebase they're often used as `ham_or_trott` since most of the time only one is needed and the simulation works in its eigenbasis
- ThermalizeResults should support convergence plotting (trace distances over time to Gibbs state)
- TrajectoryResults should support convergence plotting when convergence tracking was used
- `hamiltonian.gibbs` provides the Gibbs reference for convergence — no need to pass it as a separate argument

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 36-api-and-results*
*Context gathered: 2026-02-27*
