# Quick Task 20: Debug GNS TrotterDomain 0.83 gap

## Task
Investigate and fix the suspiciously large 0.83 GNS-to-Gibbs gap for TrotterDomain (expected ~0.08).

## Tasks

### Task 1: Diagnose basis mismatch in TrotterDomain GNS test
- **Root cause:** `construct_lindbladian` builds the Liouvillian in the Trotter eigenbasis (because jump operators are transformed to Trotter eigenbasis for element-wise OFT). The steady-state density matrix extracted from the Liouvillian eigendecomposition is therefore in the Trotter eigenbasis. But `SMALL_GIBBS` is in the energy eigenbasis. Comparing matrices in different bases gives meaningless trace distances.
- **Evidence:** Transforming the fixed point to energy eigenbasis via `U = ham.eigvecs' * trotter.eigvecs` gives gap=0.0807, matching EnergyDomain exactly (within 7e-9 Trotter error).

### Task 2: Fix GNS-01 test (line 36)
- Add Trotter-to-energy basis change before comparing to `SMALL_GIBBS`
- Restore sanity bound to `gap < 0.5`

### Task 3: Fix GNS-02 test (line 89)
- Add same basis change for the Gibbs distance baseline logging

### Task 4: Fix CPTP test (line 44)
- Pass `SMALL_TROTTER` instead of `SMALL_HAM` to `_precompute_data` for consistency
