# Phase 17: Adaptive Sampling - Context

**Gathered:** 2026-02-16
**Status:** Ready for planning

<domain>
## Phase Boundary

Convergence-driven trajectory batching with automatic stopping and hard cap. The system runs trajectory batches and monitors trace distance to determine when to stop, so the user doesn't need to specify a fixed trajectory count. Returns the same result structure as fixed-count mode.

</domain>

<decisions>
## Implementation Decisions

### Convergence criteria
- **Stop trigger: trace distance only** — observables (correlations, energy) are tracked but do NOT gate stopping
- **Relative change computed via windowed average** — compare average of last K=3 batches vs previous K=3 batches
- **3 consecutive stable checks required** — relative change must be below threshold for 3 consecutive batch checkpoints
- **Minimum 5 batches before stopping can trigger** — prevents premature stopping from lucky early batches

### Non-convergence handling
- **Return with `converged=false` flag** — no warning logged, no error thrown; user checks the flag programmatically
- **Include diagnostics** — final relative change value, number of consecutive stable batches reached, total batches run; helps user decide next steps
- **Extend existing ConvergenceData struct** — add `converged::Bool` and stop diagnostics fields directly to ConvergenceData (no new wrapper type)
- **Return immediately when converged** — no extra confirmation batches after convergence criteria are met

### Batch sizing strategy
- **Fixed batch size throughout** — every batch has the same number of trajectories
- **Default 200 trajectories per batch** — matches Phase 16 convergence tracking structure
- **Batch size configurable** — user can pass batch_size kwarg, default 200
- **Wraps `run_trajectories_convergence`** — adaptive function calls the existing Phase 16 function in a loop with increasing trajectory counts, reusing existing batch infrastructure

### Default parameters (all configurable)
- **N_max = 20,000 trajectories** (100 batches of 200) — generous ceiling for large systems at high beta
- **Relative change threshold = 0.01** (1%) — user can tighten (0.001) or loosen (0.05)
- **Patience = 3 consecutive stable batches** — user can increase for more confidence
- **Minimum batches = 5** — user can adjust burn-in period
- **Batch size = 200** — user can adjust per-batch trajectory count

### Claude's Discretion
- How to accumulate windowed averages efficiently
- Internal loop structure for wrapping run_trajectories_convergence
- Exact field names and types for diagnostics in ConvergenceData
- How to handle the seed progression across adaptive batches

</decisions>

<specifics>
## Specific Ideas

- All adaptive parameters are configurable with sensible defaults — the API should feel like "just works" out of the box but allows expert tuning
- Wrapping run_trajectories_convergence means the adaptive layer is thin — most complexity lives in Phase 16's existing batch runner

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 17-adaptive-sampling*
*Context gathered: 2026-02-16*
