# Phase 31: Scaling Benchmarks - Context

**Gathered:** 2026-02-24
**Status:** Ready for planning

<domain>
## Phase Boundary

Empirical timing and memory measurements at n=3-7 for Krylov spectral gap estimation, establishing the 4^n scaling law and producing resource estimates for n=10,12 production runs. This is a standalone benchmark script, not part of the test suite.

</domain>

<decisions>
## Implementation Decisions

### Benchmark scope
- Domains: EnergyDomain and TrotterDomain only, both with KMS balance (coherent term enabled)
- System sizes: n=3,4,5,6 guaranteed; n=7 conditional on extrapolated time fitting under 15 minutes
- Repetitions: 3-5 runs with median for n=3,4,5,6; single run for n=7 if attempted
- Execution: Standalone script in `simulations/` folder (fourth simulation kind alongside Lindbladian construction, DM, trajectory)

### Output & reporting
- Print summary tables to stdout during execution
- Save markdown report to `results/` folder (top-level, same level as `simulations/`)
- Report includes extrapolation with time + memory + feasibility notes for n=10,12
- Separate sections for EnergyDomain and TrotterDomain (not side-by-side)

### Breakdown granularity
- Total eigsolve wall-clock time only (no per-matvec breakdown)
- Peak memory allocation measured at each system size
- Memory guard thresholds to be calibrated from empirical data (current guards were set blindly)

### Extrapolation approach
- Power-law fit: t = a * b^n, expect b ~ 4
- Hard assertion that b is in [3.5, 4.5] — script errors if scaling is wrong
- Separate fits for EnergyDomain and TrotterDomain
- Report predicted time, memory, and feasibility notes for n=10 and n=12

### Krylovdim selection
- Target Krylov precision: 1e-8 absolute (well below 1e-6 physical errors from Trotter/quadrature)
- No krylovdim sweep — use heuristic values that increase with n
- At n=6: probe krylovdim=50 against dense reference for 1e-8 precision; if insufficient, try krylovdim=100 (memory permitting). Chosen krylovdim feeds into timing data and extrapolation
- Be memory-conscious at n=6+ (laptop RAM constraints)
- Only include krylovdim-vs-n recommendation table if research supports it with prior experience from the literature

### Claude's Discretion
- Exact krylovdim heuristic for n=3,4,5 (small systems, less sensitive)
- Memory measurement approach (peak allocation tracking method)
- Warmup strategy (JIT compilation handling)
- Report formatting and section structure
- How to handle n=7 time extrapolation logic

</decisions>

<specifics>
## Specific Ideas

- The user sees four simulation kinds in the project: Lindbladian construction, DM simulator, trajectory simulator, and Krylov spectral gap estimation — this benchmark is the fourth
- Krylov precision must be a few orders below physical errors (1e-8 vs 1e-6) — this is a design principle, not just a number
- The krylovdim probe at n=6 is critical because it determines the timing data that feeds into extrapolation — wrong krylovdim means wrong predictions
- Memory guards from Phase 29 were set blindly and should be calibrated from the empirical measurements

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 31-scaling-benchmarks*
*Context gathered: 2026-02-24*
