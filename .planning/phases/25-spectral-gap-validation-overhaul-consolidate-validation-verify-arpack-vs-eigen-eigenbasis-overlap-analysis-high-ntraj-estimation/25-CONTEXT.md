# Phase 25: Spectral Gap Validation Overhaul - Context

**Gathered:** 2026-02-18
**Status:** Ready for planning

<domain>
## Phase Boundary

Clean-slate rebuild of spectral gap validation. Delete all existing validation code (scripts, CrossValidationResult, cross_validate_gap, overlapping observable builders). Replace with minimal, correct infrastructure: one preset observable builder, one eigenbasis overlap diagnostic, one validation script. Verify ARPACK eigs agrees with eigen() for exact gap. Run high-statistics trajectory estimation (20k trajectories) for n=4 and n=6 at beta=10, delta=0.01. Target: relative error < 1e-2 between trajectory-estimated and exact spectral gap. If not achievable, explain WHY with evidence (e.g., observable-eigenmode overlap analysis).

</domain>

<decisions>
## Implementation Decisions

### Consolidation — fresh start
- Delete ALL existing validation scripts (experiments/validate_*.jl, experiments/run_sweep.jl, etc.)
- Delete cross_validate_gap function and CrossValidationResult struct
- Delete overlapping/redundant observable builder functions — audit what exists, keep only what's needed
- Read quick-fix summaries (quick-22, quick-23, quick-24) for context on what went wrong, but treat them as potentially containing errors and contradictions — that's why this phase exists
- Goal: minimal function count, logically organized

### Observable set
- Single function: `build_preset_trajectory_observables()` (or similar name)
- Fixed set of 5 observables: H (energy), M_z (total magnetization per site), XX_avg, YY_avg, ZZ_avg (per-bond averaged 2-site correlations)
- Per-bond averaged versions, NOT individual bond pairs
- Must be correctly transformed to the simulation basis

### Eigenbasis overlap analysis
- Separate exported diagnostic function (not embedded in validation script)
- Decomposes each observable into the Lindbladian's eigenbasis
- Reports overlap of each observable with the slowest decaying mode (first excited eigenmode of L)
- Larger overlap = better observable for gap estimation
- This is the key diagnostic: if gap estimation fails, overlap analysis explains whether the chosen observables can see the gap at all

### ARPACK vs eigen verification
- Check that ARPACK eigs method in run_lindbladian delivers the same spectral gap as eigen() for n=4
- This is a sanity check — if they disagree, something is fundamentally wrong
- If they agree, proceed with confidence that the exact gap reference is correct

### Estimation protocol
- System sizes: n=4 and n=6 (Heisenberg chain)
- Trajectories: 20,000 (high enough to rule out statistical noise)
- Parameters: beta=10, delta=0.01
- Target accuracy: relative error < 1e-2 (1%) between fitted and exact gap
- If target not met: must explain WHY with evidence, suggest followup tests
  - Possible reason: insufficient overlap with first excited mode
  - Possible reason: discrete-step Kraus effect (known from quick-22)
  - Must provide concrete diagnostic output, not just "it didn't work"

### Validation script
- One script that runs everything: exact gap computation, trajectory estimation, comparison
- Prints clear pass/fail with gap values, relative error, per-observable fit quality
- Calls the eigenbasis overlap diagnostic to show overlap table
- Lives in experiments/ directory

### Claude's Discretion
- Exact function signatures and names (within the spirit of minimal API)
- How to structure the overlap computation internally
- Whether to use KMS or GNS detailed balance (whatever the codebase currently defaults to)
- Test structure (integration test vs script-only validation)

</decisions>

<specifics>
## Specific Ideas

- "I think some hallucination has been going on and we could not settle on what the actual solution is" — prior quick fixes may have introduced contradictory logic. Start fresh, trust only the math.
- The overlap with the worst decaying mode is the key physics insight: an observable that doesn't overlap with the first excited eigenmode of L cannot estimate its decay rate, no matter how many trajectories you run.
- beta=10, delta=0.01 are the fixed validation parameters — not a sweep, just these values.

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 25-spectral-gap-validation-overhaul*
*Context gathered: 2026-02-18*
