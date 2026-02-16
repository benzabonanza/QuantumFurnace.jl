# Quick Task 20 Summary: Debug GNS TrotterDomain 0.83 gap

## Root Cause
**Basis mismatch in test_gns_trajectory.jl.** The TrotterDomain Liouvillian is built in the Trotter eigenbasis (jump operators transformed via `transform_jumps_to_basis`), so its fixed point is also in the Trotter eigenbasis. But `SMALL_GIBBS` (from `hamiltonian.gibbs`) is in the energy eigenbasis. Comparing trace distance across different bases gave the spurious 0.83 gap.

## Diagnostic Evidence
- Trotter quasi-energies match exact eigenvalues to ~1e-9 (10 Trotter steps)
- Trotter Bohr frequency max difference: 0.9 (due to reversed eigenvector ordering from `eigen` on unitary)
- After basis transformation: TrotterDomain gap = 0.0807, matching EnergyDomain = 0.0807 (diff 7e-9)

## Changes Made

### test/test_gns_trajectory.jl
1. **GNS-01 (line 35-42):** Added `U_t2e = SMALL_HAM.eigvecs' * SMALL_TROTTER.eigvecs` basis change before comparing fixed point to `SMALL_GIBBS`. Restored sanity bound from `< 1.0` to `< 0.5`.
2. **GNS-01 CPTP (line 44):** Fixed `_precompute_data(config, SMALL_HAM)` → `_precompute_data(config, SMALL_TROTTER)` for consistency.
3. **GNS-02 (line 91-94):** Added same Trotter-to-energy basis change for Gibbs distance baseline.

## Verification
- All 284 tests pass
- TrotterDomain GNS-to-Gibbs gap: 0.0807 (was 0.83)
- BohrDomain GNS-to-Gibbs gap: 0.0352 (unchanged)
- Trajectory convergence: 0.038 trace distance (unchanged)
