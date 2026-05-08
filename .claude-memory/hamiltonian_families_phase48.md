---
name: Phase 48 Hamiltonian families and 2D builders
description: Three families of cached Heisenberg Hamiltonians for the thesis numerics chapter (epic qf-k1u, task qf-k1u.5), including the new 2D builders.
type: project
---

# Phase 48 Hamiltonian Families (qf-k1u.5, completed 2026-04-30)

Cached BSON Hamiltonians under `hamiltonians/` for the thesis numerics chapter:

  F1. `heis_xxx_zzdisordered_periodic_n{3..10}.bson` — 1D XXX
      with full `[[Z], [Z,Z]]` disorder (strength 1.0, batch_size=200).
      Replaces the legacy Z-only `heis_disordered_periodic_n*.bson`,
      which had bipartite-pairing Bohr collisions for even n.
  F2. `heis_xxx_clean_periodic_n{3..10}.bson` — 1D XXX with weak
      `[[Z]]` ε-disorder (strength 1e-2, batch_size=500). Per-bond
      ground-state energy verified against the Bethe-ansatz limit
      `4 (1/4 − ln 2) ≈ −1.7726` (Pauli normalization). Even-n values
      hover above the limit, odd-n values below — classic finite-size
      frustration of the XXX AFM ring.
  F3. `heis_xxz_2d_{Lx}x{Ly}_n{n}.bson` for `(Lx,Ly,n) ∈
      {(2,2,4), (2,3,6), (3,3,9), (2,5,10)}` — 2D anisotropic XXZ
      `[J_x, J_y, J_z] = [1, 1, 1.5]` on a square lattice with
      periodic BCs in both directions, plus weak `[[X], [Z]]`
      ε-disorder (strength 1e-2, batch_size=500). The bulk model has
      a finite-T Néel transition (Z₂ Ising symmetry breaking — allowed
      by Mermin–Wagner because it's discrete, unlike the SU(2)-
      symmetric XXX which has no finite-T transition).

## Disorder choice rationale (2D)

Random Z-only disorder *commutes with Sz*, leaving the spectrum
block-diagonal in Sz sectors and producing astronomically small
cross-Sz accidental gaps (~1e-9). Random `[[X], [Z]]` (transverse +
longitudinal field) breaks Sz, lattice translations, and the
spin-flip symmetry ∏ X_i, achieving nu_min ≈ 1e-4 at small
lattices. Random X-only disorder has its own residual symmetry
(∏ X_i still commutes with H + ε ΣX) → essential 2-fold doublets,
nu_min collapses to machine epsilon.

## Builders (`src/hamiltonian.jl`)

  - `_pad_two_site_op(term, num_qubits, q1, q2)` — places a 2-site
    Pauli at non-adjacent qubit indices. Required for 2D right-
    neighbour bonds where the qubit-index step is `Ly`, not 1.
  - `_construct_2d_heisenberg_base(Lx, Ly, terms, coeffs;
    periodic_x, periodic_y)` — returns dense Hermitian 2D base
    Hamiltonian. Skips x-direction bonds when `Lx==1` (and
    similarly y), avoiding self-bonds. Lx==2 (or Ly==2) with PBC
    double-counts each wrap-around bond; this matches the existing
    1D `n=2` PBC convention (the `_construct_base_ham` loop does
    the same).
  - `find_ideal_2d_heisenberg(Lx, Ly, coeffs; ...)` — wraps the
    shared `_optimize_disordered_heisenberg` inner kernel.
  - `find_ideal_heisenberg(num_qubits, coeffs; ...)` — gained a
    `disorder_strength::Float64` keyword. Default 1.0 preserves
    pre-refactor behaviour; ε-disorder use cases pass 1e-2 or 1e-3.

## Acceptance numbers (achieved)

  F1 nu_min by n: n=3 → 2.2e-2, n=4 → 9.7e-3, n=5 → 2.0e-3, n=6 →
  4.4e-4, n=7 → 1.5e-4, n=8 → 2.5e-5, n=9 → 6.1e-6, n=10 → 2.5e-6.

  F2 nu_min by n: n=3 → 5.5e-4, n=4 → 2.0e-4, n=5 → 7.7e-5, n=6 →
  4.9e-7, n=7 → 1.6e-6, n=8 → 1.6e-7, n=9 → 9.0e-8, n=10 → 3.3e-8.

  F2 e_GS/bond converges toward Bethe limit -1.7726:
  n=3 → -1.00, n=4 → -2.00, n=5 → -1.50, n=6 → -1.87, n=7 → -1.63,
  n=8 → -1.83, n=9 → -1.69, n=10 → -1.81.

  F3 ⟨S_stag²⟩/n² (2D Néel correlator): 2x2 → 0.73 (clearly ordered;
  bipartite-PBC compatible), 2x3, 3x3, 2x5 → smaller (geometric
  frustration from PBC with odd Ly creates non-bipartite wrap).
  Bulk physics gives long-range Néel; small frustrated PBC lattices
  give partial correlator and a non-trivial Gibbs state — exactly
  the "less obvious" Gibbs target the thesis simulator should learn
  to thermalize to.

## File-size scaling

BSON file size for find_ideal_* (NamedTuple + saved with @save) is
about 0.4× of the legacy 14-field HamHam serialisation (no nothing
fields for bohr_freqs/bohr_dict/gibbs):

  n=10: legacy 50 MB → new family ~20 MB (n=10 NamedTuple).

## Test fixtures still use legacy files (migration tried and reverted)

`test/test_helpers.jl::make_test_system` loads
`heis_disordered_periodic_n{3,4,5}.bson` (n=5 from
`test/test_discriminant.jl`'s loop). I migrated to the new family-1
files + regenerated `test/reference/*.bson`, but one Krylov eigsolve
test failed on the new Hamiltonian: the channel→Lindbladian
spectral-gap relative error was 3.8e-3, just over the 2e-3
calibrated tolerance. The original fixture's tolerance budget is
narrowly tuned to the old eigenstructure. Reverted the migration —
keep legacy n=3, n=4, n=5 BSONs in place for tests. Removed only
the unused legacy n=6..10 files. A future migration would also need
to re-tune the rtol in `test/test_krylov_eigsolve.jl:168`.

## Diagnostic scripts still load legacy files

The following scripts load `heis_disordered_periodic_n*.bson` via
the legacy 14-field BSON.parse pattern:

  - `scripts/biexp_mixing_verify.jl`
  - `scripts/biexp_mixing_verify_trotter.jl`
  - `scripts/mixing_time_extrapolate_verify.jl`
  - `scripts/simulation_time_budget.jl`
  - `scripts/scratch_gqsp_thermalization_check.jl`
  - `scripts/scratch_gqsp_trace_drift.jl`

They keep working with the retained legacy n=3..5 files (and any
larger n the user regenerates locally). Migration to the new format
is a one-line replacement per script:
`hamiltonian = HamHam(BSON.load(path)[:hamiltonian], beta)`.

## Soft nu_min target

The simulator does not require `nu_min ≥ 1e-4`. The generation script
warns when it falls below this but does not abort; consumers (M1/M2/E1
plot scripts) just need finer simulator parameters when nu_min is
small. Closely-spaced Bohr frequencies become distinct dictionary
keys (exact float equality), so the algorithm handles them correctly
provided the OFT grid resolves them.
