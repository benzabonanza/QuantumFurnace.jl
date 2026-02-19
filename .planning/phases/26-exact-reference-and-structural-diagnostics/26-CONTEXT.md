# Phase 26: Exact Reference and Structural Diagnostics - Context

**Gathered:** 2026-02-19
**Status:** Ready for planning

<domain>
## Phase Boundary

Build exact Lindbladian spectral data and structural diagnostics (eigenvalues, fixed point, anti-Hermitian defect, symmetry sectors, observable overlaps) as ground truth for validating trajectory-based gap estimates at n=4,6. This phase produces reference data consumed by Phases 27-30.

</domain>

<decisions>
## Implementation Decisions

### Defect ratio interpretation
- Advisory warning only — report ||A||/lambda_gap(H) and flag if high, but do NOT gate or block downstream fitting
- Threshold value to be determined from supplementary information papers and numerical evidence at n=4,6 (Claude's discretion)
- KMS similarity transform uses a fixed internal spectrum truncation threshold (no user-facing parameter) — chosen to be robust for our n=4,6 cases
- Report scalar ratio only (no full spectrum decomposition of anti-Hermitian part A) — full eigen() infeasible above n=6
- Imaginary-to-real eigenvalue ratios |Im(lambda_k)/Re(lambda_k)| also reported for leading modes (per supplementary info Task 1.2)

### Observable overlap scope
- Canonical observable set from supplementary info + Mz_stagg:
  - Z_1 (Pauli-Z on first site)
  - X_1 (Pauli-X on first site)
  - Z_1 * Z_{floor(n/2)} (two-point correlator)
  - H (the Hamiltonian itself)
  - Random traceless Hermitian matrix (control, fixed seed for reproducibility)
  - Mz_stagg (staggered magnetization, user addition)
- Replace existing v1.3 observable set (XX_avg, YY_avg, ZZ_avg, XZ_stagg, etc.) with this new canonical set in build_preset_observables
- Compute c_k for best-coupling observables across all three initial states:
  - |0>^n (all spins up)
  - |+>^n (all spins in X-plus state)
  - I/2^n (maximally mixed state)
- 20 leading modes for overlap coefficient computation
- c_k = Tr[O R_k] * Tr[L_k^dagger (rho_0 - rho_beta)] per supplementary info formula

### Symmetry sector reporting
- General Sz-based analysis (works for any Hamiltonian that conserves total Sz: Heisenberg, Ising, XXZ, etc.)
- Delta_Sz labels assigned to each Lindbladian eigenvector based on density matrix support structure

### Eigenvalue extraction defaults
- Default: 20 leading eigenvalues via Arpack shift-invert (nearest to zero)
- Store eigenvalues as Vector{ComplexF64} (natural representation, caller decomposes)
- Extract both left AND right eigenvectors upfront in DIAG-01 (both needed for overlap computation)
- API: standalone functions for each diagnostic PLUS a `run_exact_diagnostics()` bundle that calls them all and returns a single result struct

### Claude's Discretion
- Defect ratio warning threshold (determined from supplementary papers and numerical evidence)
- Near-degeneracy detection and multiplet grouping approach for symmetry sectors
- Handling of mixed-sector eigenvectors (dominant sector + purity fraction vs threshold-based labeling)
- Dashboard visualization choices (color-by-sector, marker shapes for the spectrum plot)
- Fixed internal truncation threshold for KMS similarity transform

</decisions>

<specifics>
## Specific Ideas

- The supplementary information (spectral-gap-refinements-instructions.md) provides detailed formulas for all computations:
  - KMS similarity transform: D(rho, L) = rho^{-1/4} L[rho^{1/4} (.) rho^{1/4}] rho^{-1/4}
  - Overlap formula: c_k = Tr[O R_k] * Tr[L_k^dagger (rho_0 - rho_beta)]
  - Symmetry labeling via Delta_Sz = Sz(E_i) - Sz(E_j) of eigenvector support
- Fixed point trace distance should be ~0 for BohrDomain, ~1e-8 for TrotterDomain (per STATE.md)
- n=6 spectral gap is approximately 0.11 (per supplementary info)
- Full eigen() only feasible up to n=6 (4^6 = 4096 Lindbladian dimension); n=8 deferred to v2

</specifics>

<deferred>
## Deferred Ideas

- n=8 sparse Lindbladian diagonalization via KrylovKit.jl — explicitly deferred to v2 (EVAL-01)
- Damped-oscillation fit model c*exp(-gamma*t)*cos(omega*t+phi) — deferred to v2 (EVAL-02)
- GEVP / matrix pencil methods for multi-eigenvalue extraction — deferred to v2 (EVAL-03)

</deferred>

---

*Phase: 26-exact-reference-and-structural-diagnostics*
*Context gathered: 2026-02-19*
