---
name: CKG EnergyDomain production sweep (qf-lkb.11)
description: Switch CKG KMS Lindbladian sweeps from BohrDomain (O(d⁴)) to EnergyDomain (O(N·d²) at fixed N=4096), unblocking n≥5 mixing-time computations
type: project
originSessionId: 55d42cf9-ca55-414e-aede-9b64d6027ec7
---
## qf-lkb.11 — Speed up CKG KMS BohrDomain Lindbladian matvec (closed 2026-05-02)

### Why
CKG BohrDomain matvec was the bottleneck for n≥5 mixing-time sweeps:
n=4 wall ≈ 100s, n=5 ≈ 33min/cell extrapolating d^4. The matvec calls
`_pick_alpha(config, …)` ≈ d²·|B_H| ≈ d⁴ times per call, with 3× erfc per
evaluation in the smooth-Metropolis branch.

CKG EnergyDomain expresses the same KMS-DB Lindbladian in OFT-smeared form,
where the alpha factorisation `α(ν1,ν2) = G(ν1+ν2) · exp(-(ν1-ν2)²/(8σ²))`
becomes a Riemann sum over a fixed N=2^num_energy_bits=4096 grid. Per matvec
scales as O(num_jumps · N · d²), no erfc on the hot path.

### Headline numbers (laptop, 4 BLAS, 4 Julia threads)
- Bohr→Energy speedup: 6× at n=3, 22.5× at n=4 (matvec wall).
- EnergyDomain matvec scaling: 0.32 → 1.84 → 13.8 → 137 → 322 ms (n=3..7), ~d^2.5.
- Integration wall (4 betas): n=3 ~5s, n=4 ~3.5s, n=5 ~27s, n=6 ~170s, n=7 ~24min.
- Dense Liouvillian agreement Bohr ≈ Energy: ~6e-10 rel at n=3,4,5 β=10
  (FP-accumulation floor, well below quadrature error ~exp(-78)).

### Decisions locked
1. **Skipped Path A** (BohrDomain alpha caching): user direction —
   memory blow-up at n>5 (`O(d⁴)` cache). Path B is the only path.
2. **Skipped Task 11.3** (extract `_apply_lindbladian_kraus_body!`):
   user direction in description — Bohr/Energy/DLL have genuinely
   different inner-product structure; routing through CKG EnergyDomain
   (the natural Kraus form) is the production path.
3. **Smooth-Metropolis defaults locked at a=0, s=0.25** (was a=β/30, s=0.4)
   — the thesis-numerics convention. Applied uniformly in
   `sweep_mixing_times`, scaling benchmark, and qf-lkb.6 plot regen.

### API additions
- `integrate_to_gibbs` accepts `Config{Lindbladian, <:Union{BohrDomain, EnergyDomain}}`.
- `sweep_mixing_times` gains `domain::AbstractDomain=BohrDomain()`,
  `a::Real=0.0`, `s::Real=0.25`, `hamiltonian_filename::Function`
  kwargs. Result NamedTuple gains `:domain` field. Sidecar filename
  embeds domain tag (Bohr|Energy) to prevent cross-domain cache collisions.
- `_load_hamiltonian_bson` now handles both schemas: legacy 14-field
  HamHam-typed BSON and the new 11-field NamedTuple-typed BSON
  (heis_xxx_zzdisordered, heis_xxx_clean, heis_xxx_2d families).

### Test coverage added
- `test_lindblad_action.jl (q0)`: dense Bohr ≈ Energy at n∈{3,4,5}, β=10,
  rel < 1e-9 (achieved ~6e-10).
- `test_lindblad_action.jl (q1)/(q2)`: EnergyDomain CKG end-to-end mixing
  at n=3, β=10, both mode=:L and mode=:K (chi² Lyapunov check on K).
- `test_krylov_matvec.jl 'EnergyDomain CKG @ a=0, s=0.25 (qf-lkb.11.4)'`:
  zero-byte allocation regression on the production smooth-Metropolis config.

### Files changed (5 atomic commits)
- `src/lindblad_action.jl`: dispatch relaxation, kwargs, domain-tagged sidecars.
- `src/misc_tools.jl`: dual-schema BSON loader.
- `test/test_lindblad_action.jl`: 12 new tests (q0, q1, q2 bundles).
- `test/test_krylov_matvec.jl`: 2 new allocation tests.
- `scripts/numerics_ckg_vs_dll_taumix_comparison.jl`: CKG → EnergyDomain,
  N_VALUES extended {3,4} → {3,4,5}, all sweeps a=0 s=0.25.
- `scripts/scratch_bench_ckg_energy_scaling.jl`: new scaling benchmark.

### Why: see project memory `dll_performance_refactor_qf_hur.md` for the
prior DLL refactor that hit the same crossover problem from the DLL side
(n>5 cliff in DLL apply, closed by qf-lkb.9 matrix-free path).

### How to apply: when a downstream sweep needs CKG at n≥5, pass
`domain=EnergyDomain()` to `sweep_mixing_times`. Keep `a=0, s=0.25` as
defaults unless explicitly studying transition-weight shape (in which case
override).
