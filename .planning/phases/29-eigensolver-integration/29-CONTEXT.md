# Phase 29: Eigensolver Integration - Context

**Gathered:** 2026-02-24
**Status:** Ready for planning

<domain>
## Phase Boundary

Single-call API (`krylov_spectral_gap()`) wrapping KrylovKit eigsolve to compute spectral gaps matrix-free, with both Lindbladian and CPTP channel targeting. All four domains (Energy, Time, Trotter, Bohr) supported. Includes pre-flight memory estimation and convergence retry logic.

</domain>

<decisions>
## Implementation Decisions

### API surface design
- Single function `krylov_spectral_gap(config, ham)` — config type dispatches the path
- `LindbladConfig` → Lindbladian eigsolve with `:LR` targeting (largest real part)
- `ThermalizeConfig` → CPTP channel eigsolve with `:LM` targeting (largest magnitude), using delta from config
- All four domains supported from day one (matvec already exists from Phase 28)
- Workspace allocated internally — no user-facing workspace kwarg
- Conservative defaults: `krylovdim=30`, `howmany=4`, `tol=1e-10` — user tunes down for speed
- KrylovKit kwargs passed through as keyword arguments for full control

### Result struct contents
- `KrylovGapResult` stores top 2 eigenvalues and eigenvectors (fixed-point + gap mode)
- Eigenvectors reshaped to dim x dim matrices (not raw vectors)
- Unified `gap` field — always represents Lindbladian spectral gap regardless of path
- Channel path auto-converts channel eigenvalues back to Lindbladian gap internally
- Convergence info and matvec count included as metadata
- **RESEARCH FLAG:** Exact channel→Lindbladian gap conversion formula needs careful verification during planning (log vs linear approximation, complex eigenvalue handling)

### CPTP channel strategy
- `apply_delta_channel!` computes E(rho) = (I + delta*L)(rho) using delta from ThermalizeConfig
- Delta used as-is from config — no validation or auto-selection
- Channel matvec lives in same Krylov module alongside Lindbladian matvec
- `:LM` targeting for channel (eigenvalues cluster near 1)
- `:LR` targeting for Lindbladian (steady state at 0, gap from lambda_2)

### Failure and guard behavior
- Pre-flight memory estimate: `krylovdim * 4^n * 16 * 1.5` bytes
- Threshold: 80% of `Sys.free_memory()` — adapts to the machine
- Memory guard issues `@warn` and proceeds (not a hard error)
- Convergence retry strategy: Claude's discretion (research KrylovKit behavior to pick between doubling vs incremental krylovdim increase)
- Final failure: simple error message — "KrylovKit failed to converge: n_converged/howmany eigenvalues after N retries"

### Claude's Discretion
- Convergence retry strategy (krylovdim increase amount and number of retries)
- Internal workspace allocation details
- KrylovKit algorithm selection (Arnoldi vs other iterative methods)
- Exact channel→Lindbladian gap conversion formula (after research)

</decisions>

<specifics>
## Specific Ideas

- Config type dispatch: LindbladConfig → Lindbladian path, ThermalizeConfig → channel path — no extra kwargs needed, Julia's type system handles it
- "I want auto-convert internally but I wanna make sure you are computing the right thing" — research the exact conversion formula carefully before implementing
- Top 2 eigenpairs give fixed-point state, gap mode, and spectral gap — that's all we need

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 29-eigensolver-integration*
*Context gathered: 2026-02-24*
