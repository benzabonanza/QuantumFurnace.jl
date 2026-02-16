# Phase 18: KMS-vs-GNS Experiments - Context

**Gathered:** 2026-02-16
**Status:** Ready for planning

<domain>
## Phase Boundary

Run a parameter sweep comparing KMS (exact detailed balance) vs GNS (approximate detailed balance) convergence across system sizes and inverse temperatures. Produce paper-ready data as BSON files. The full grid is n=4,6,8 x beta=5,10,20 x {KMS, GNS@sigma=1/beta, GNS@sigma=0.5/beta} = 27 experiments using TrotterDomain.

</domain>

<decisions>
## Implementation Decisions

### Hamiltonian model
- 1D Heisenberg chain (XXX) with periodic boundaries, J=1.0 uniform coupling
- Experiment function accepts a pre-built Hamiltonian object (not constructed internally) — supports future Hamiltonian models
- System sizes n=4,6,8 are qubit counts (Hilbert space dims 16, 64, 256)
- Trajectory-based simulation uses local jump operators (TrotterDomain), no full Lindbladian construction — n=8 is feasible

### Simulation parameters
- delta = 0.01 for all 27 experiments (fixed, same as Phase 14 baseline)
- Adaptive sampling with shared n_max = 10000 cap across all experiments
- Rationale for n_max: trajectory sampling error should be subdominant to the delta=0.01 systematic error (~0.03 trace distance)
- Mixing time, batch size, convergence thresholds: Claude's discretion (see below)

### Data organization
- Output directory: `experiments/` in project root
- Descriptive file names: e.g., `kms_n4_beta5.bson`, `gns_sigma1_n6_beta10.bson`
- No summary/index file — individual BSON files contain all metadata, glob the directory
- `experiments/` is gitignored — data is regenerable

### Execution strategy
- Per-experiment function: `run_experiment()` called in a loop over the parameter grid
- Sweep script lives in `experiments/run_sweep.jl` — standalone script using the package, not part of package API
- Failure handling: skip failed experiments, log the error, save what completed, continue with remaining — print warning at end listing failures
- No resume capability — always rerun all experiments (overwrite existing files)

### Claude's Discretion
- Mixing time value (scaling with beta or fixed)
- Batch size for adaptive sampling
- Convergence thresholds (relative change threshold, window size)
- Exact file naming format details
- How `run_experiment()` is structured internally (kwargs, return type)
- Progress/logging output during sweep

</decisions>

<specifics>
## Specific Ideas

- "Trajectory sampling error should always be a bit better than the delta algorithmic error" — n_max=10000 chosen to ensure statistical error is subdominant to Trotter error
- Make sure the experiment API works with future Hamiltonian objects too (accept pre-built Hamiltonians, don't hardcode Heisenberg)
- Quick-20 established baselines: GNS-to-Gibbs gap at sigma=0.1 is ~0.08 for TrotterDomain

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 18-kms-vs-gns-experiments*
*Context gathered: 2026-02-16*
