# Phase 44: Struct and Grid Infrastructure - Context

**Gathered:** 2026-03-04
**Status:** Ready for planning

<domain>
## Phase Boundary

Result container (`SimulationTimeBudget`) and QPE grid utilities that all subsequent cost computations (phases 45-47) build upon. This phase delivers the data structures and grid arithmetic — no actual cost computation yet.

</domain>

<decisions>
## Implementation Decisions

### Struct metadata scope
- Store full physical context so the struct alone reproduces a paper table row
- Fields include: oft_time, b_time, per_step_time, n_steps, total_time (computed costs)
- Grid info: r, N, w0, t0, energy_range
- Physical parameters: beta, sigma, delta, construction type, n_qubits, rescaling_factor, mixing_time
- Filter info: filter type (:gaussian, :metropolis, :smooth_metropolis) and parameters (sigma for Gaussian, a/b for smooth Metropolis)
- Scalars only from HamHam (n_qubits, rescaling_factor) — no HamHam reference stored
- Struct is immutable

### Display/printing format
- Two-method Julia convention: compact `show` (one-liner) + verbose `show(io, MIME"text/plain", ...)` (multi-line table)
- Verbose display includes formula breakdown: `Per step: 2460.0 (2x1200.0 + 60.0)` making cost structure transparent
- Times shown in raw scientific notation (no human-friendly suffixes — physicists read scientific notation)
- Verbose display includes filter info (e.g., `Filter: Gaussian(sigma=0.5)`)

### Claude's Discretion
- Whether to include sub-breakdowns (raw QPE time sum, number of b-dict terms) beyond the main time totals
- Grid return type for `_qpe_grid_info` (NamedTuple vs struct vs individual values)
- Validation strictness on construction (which physical constraints to enforce)
- Exact compact one-liner format

</decisions>

<specifics>
## Specific Ideas

- The struct should be self-contained enough that someone can look at it and understand the full cost breakdown without needing the original inputs
- Formula transparency: showing `2xOFT + B` in the display helps verify correctness at a glance

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 44-struct-and-grid-infrastructure*
*Context gathered: 2026-03-04*
