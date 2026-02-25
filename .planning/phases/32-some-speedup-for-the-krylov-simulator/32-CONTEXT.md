# Phase 32: Krylov Simulator Speedup - Context

**Gathered:** 2026-02-25
**Status:** Ready for planning

<domain>
## Phase Boundary

Precompute aggregate Lindbladian matrices (R_total, effective Hamiltonian) at workspace construction time to reduce per-matvec GEMM count. Clean up legacy Euler channel code. No new capabilities -- purely performance optimization and dead code removal.

</domain>

<decisions>
## Implementation Decisions

### Precomputed effective Hamiltonian for Lindbladian matvec
- Currently each dissipator term computes L'L (1 GEMM) and anticommutator (2 GEMMs) per term, totaling 5N GEMMs for N = jumps x energy_labels
- Since Krylov always applies the FULL Lindbladian (unlike trajectory which samples individual jumps), precompute R_total = sum_i scalar_i * L_i'L_i at workspace construction
- Define G_left = i*B^T - 0.5*R_total^T and G_right = -i*B^T - 0.5*R_total^T (the non-Hermitian effective Hamiltonian in kron convention)
- Per-matvec becomes: rho_out = G_left*rho + rho*G_right + sandwiches (2 GEMMs constant + 2 GEMMs per sandwich term)
- This reduces from 5N GEMMs to 2 + 2N GEMMs -- roughly 2.5x fewer BLAS calls
- R_total accumulation already exists (`_accumulate_R_total!` in krylov_workspace.jl) for the channel path -- reuse for Lindbladian path
- Apply to all domains: EnergyDomain, TimeDomain, TrotterDomain, BohrDomain (Bohr has different dissipator structure but same R_total precomputation applies)
- Similarly precompute for adjoint: adjoint just swaps G_left and G_right

### Delete legacy Euler apply_delta_channel!
- Remove the 5-argument `apply_delta_channel!(ws, rho, delta, config_liouv, hamiltonian)` Euler approximation (E(rho) = rho + delta*L(rho))
- It is faulty: doesn't retain correct O(delta^2) error properties of Chen's CPTP algorithm
- We only ever need the faithful Chen version (4-argument form) which precomputes R_total, K0, U_residual
- Also delete the corresponding test: "apply_delta_channel! legacy Euler" testset in test_krylov_eigsolve.jl

### Energy label count is small (truncated)
- Energy labels are truncated by `_truncate_energy_labels` from 2^num_energy_bits to just a few dozen relevant labels (based on transition function cutoff 1e-12)
- This means the inner loop (over energy_labels) is ~20-50 iterations, not 2^n
- The optimization still matters: saving 3 GEMMs per iteration x ~30 iterations x ~100 matvecs is significant

### Claude's Discretion
- Whether to store G_left and G_right as separate fields or compute G_right = -conj(G_left) on the fly
- Exact scratch matrix reuse strategy for the sandwich-only loop
- Whether BohrDomain needs a separate precomputed matrix or can share the same R_total pattern

</decisions>

<specifics>
## Specific Ideas

- R_total accumulation helpers already exist in krylov_workspace.jl for the channel path; reuse the same `_accumulate_R_total!` functions
- The kron convention used throughout: coherent is i*B^T*rho - i*rho*B^T, anticommutator uses (L'L)^T, sandwich is conj(L)*rho*L^T
- For the adjoint Lindbladian, G_left_adj = G_right and G_right_adj = G_left (coherent sign flips, anticommutator same)
- Keep existing round-trip correctness tests -- they validate the optimized matvec produces identical results to the current implementation

</specifics>

<deferred>
## Deferred Ideas

None -- discussion stayed within phase scope

</deferred>

---

*Phase: 32-some-speedup-for-the-krylov-simulator*
*Context gathered: 2026-02-25*
