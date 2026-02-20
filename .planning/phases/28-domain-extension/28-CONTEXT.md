# Phase 28: Domain Extension - Context

**Gathered:** 2026-02-20
**Status:** Ready for planning

<domain>
## Phase Boundary

Matrix-free Lindbladian action (`apply_lindbladian!` and `apply_adjoint_lindbladian!`) for TimeDomain, TrotterDomain, and BohrDomain. Extends Phase 27's EnergyDomain matvec to all four domains. Round-trip tested against dense `construct_lindbladian()` at n=4 for each domain with both KMS and GNS balance.

</domain>

<decisions>
## Implementation Decisions

### BohrDomain dissipator structure
- New `_accumulate_dissipator_2op!(out, A, B, rho, scalar, ws)` helper — separate function, keeps existing single-op helpers untouched
- BohrDomain uses two-operator dissipator D(A, B, rho) = A rho B' - 0.5(B'A rho + rho B'A) where A and B are different matrices
- Entrywise alpha computation in the hot path: `@. alpha(bohr_freqs, nu_2) * eigenbasis` — cannot precompute because too many Bohr frequencies to store
- BohrDomain iterates over `bohr_dict` keys (bucket iteration), fundamentally different loop from Energy/Time/Trotter
- A_nu2_dag scratch: do NOT add a new workspace field for Bohr — BohrDomain is not the primary code path. Reuse existing scratch matrices creatively (Claude's discretion on which to repurpose)

### Adjoint scope
- Forward + adjoint for all three new domains — full parity with EnergyDomain
- For Time/Trotter adjoint: same NUFFT prefactors, just swap dissipator sandwich (coherent sign flip + sandwich swap, no prefactor changes)
- For BohrDomain adjoint: swap arguments at call site — `_accumulate_dissipator_2op!(out, B, A, rho, scalar, ws)` with A and B swapped. No separate adjoint 2op function needed (researcher should verify the math)

### Plan structure / domain ordering
- Plan 1: Time + Trotter together (mechanical, structurally similar to EnergyDomain — swap `_krylov_oft!` for NUFFT prefactor multiply, change scalar prefactor formula)
- Plan 2: Bohr separately (new 2op helper, different loop structure, bucket iteration)
- EnergyDomain (Phase 27) is the template — Time/Trotter are the closest structural analogs

### Round-trip testing
- Per-domain round-trip against dense `construct_lindbladian()` at n=4 (no cross-domain tests)
- Both KMS and GNS balance for all three domains
- Allocation regression tests for Time/Trotter only — Bohr may have unavoidable allocations from entrywise alpha computation

### Workspace and API
- Uniform matvec signature: `apply_lindbladian!(ws, rho, config, hamiltonian)` — workspace has everything, same API across all domains
- TrotterDomain: trotter data already in `ws.precomputed_data` from construction, no extra argument needed
- BohrDomain: access `hamiltonian.bohr_dict` and `hamiltonian.bohr_freqs` at call time (not stored in workspace)
- Time/Trotter OFT: `@. ws.jump_oft = ws.jump_eigenbases[k] * nufft_prefactor_matrix` — use concrete-typed eigenbases to avoid boxing, trust existing `_prefactor_view` for zero-alloc access

### Claude's Discretion
- Dense vs sparse for A_nu2_dag in BohrDomain (whichever avoids allocation best while keeping code clear)
- Which existing scratch matrix to repurpose for BohrDomain's second operator
- Internal helper decomposition within each domain's apply_lindbladian! method
- Test matrix generation strategy for round-trip validation

</decisions>

<specifics>
## Specific Ideas

- "Follow the pattern of DM simulator codes" — domain-dispatched config types, same structural patterns as existing `_jump_contribution!` methods
- EnergyDomain from Phase 27 is the direct template for Time/Trotter — only the OFT computation and prefactor formula change
- Time/Trotter share a single `_jump_contribution!` in the dense code (`Union{TimeDomain, TrotterDomain}` dispatch) — the Krylov versions should mirror this shared structure
- Bohr is secondary to Energy/Time/Trotter in importance

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 28-domain-extension*
*Context gathered: 2026-02-20*
