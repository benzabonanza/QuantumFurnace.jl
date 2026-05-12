# Dissipator Quadrature at β_phys = 2 — Smooth Metropolis (default_smooth_s)

**Scope caveat**: These numbers measure the **dissipator term only** —
the matvec ran with `include_coherent = false`, so the coherent
commutator $i[B, \rho]$ is masked. Full-Lindbladian register sizing
needs the coherent registers `r_b_minus, r_b_plus` sized separately
(see [v2 recipe](../../.claude-memory/quadrature_register_recipe_v2.md)).
The `r_D = 7` cutoff below is **NOT** the full-Lindbladian register
recipe at β_phys = 2.

**Task**: qf-w2u. How much harder is the EnergyDomain dissipator quadrature
at β_phys = 2 vs β_phys = 1 for smooth Metropolis (`default_smooth_s`)?

**Method**: matrix-free on both sides. Test
`L_E(r_test)` against reference `L_E(r_ref = 12)` via custom power iteration
on `(L_E(r) − L_E(r_ref))^† (L_E(r) − L_E(r_ref))` reading the norm as
`‖A x_final‖`. No BohrDomain dense build anywhere — at n=6 the dense `L_b`
costs 256 MB and ~98 s/cell. The reference at `r_ref = 12` sits ~5 bits
above the saturation point of `‖L_E(r) − L_B‖_op`, so `L_E(r_ref) ≡ L_B`
within machine matrix-element precision.

Both `apply_lindbladian!` calls per matvec use the ω-loop threading
(`OMEGA_THREAD_THRESHOLD = 10`, src/krylov_matvec.jl:122) — `JULIA_NUM_THREADS=8`,
`BLAS=1`.

## Headline

**β_phys = 2 costs at most one extra bit on `r_D` vs β_phys = 1, and r_D = 7
suffices for ε = 10⁻⁹ at every n ∈ {3, 4, 5, 6}.** The smooth Metropolis
filter narrows by 4× (kink width ∝ 1/(β·s) with `s = (0.05/σ)² ∝ β²`),
but EnergyDomain Riemann-sum convergence is super-exponential enough that
the same r_D bucket (6–7) still hits machine precision.

## β_phys = 2 cutoffs

| n | β_alg | σ | s = default_smooth_s | ω_range | r_D@1e-3 | r_D@1e-6 | r_D@1e-9 |
|---|-------|---|----------------------|---------|----------|----------|----------|
| 3 | 43.72  | 0.0229  | 4.78  | 1.27 | 6 | 6 | **6** |
| 4 | 77.38  | 0.0129  | 14.97 | 1.11 | 6 | 7 | **7** |
| 5 | 88.50  | 0.0113  | 19.58 | 1.08 | 6 | 7 | **7** |
| 6 | 112.22 | 0.00891 | 31.48 | 1.04 | 6 | 7 | **7** |

Full per-r error sweeps (operator-norm of `L_E(r) − L_E(r_ref=12)`):

| n | r=3 | r=4 | r=5 | r=6 | r=7 | r=8 |
|---|-----|-----|-----|-----|-----|-----|
| 3 | 6.1e-1 | 1.2e-1 | 1.8e-3 | 3.3e-11 | **5.8e-16** | — |
| 4 | 7.6e-1 | 2.7e-1 | 4.1e-2 | 8.9e-6  | **2.1e-16** | — |
| 5 | 6.0e-1 | 2.4e-1 | 5.1e-2 | 9.5e-5  | **2.3e-16** | — |
| 6 | 4.2e-1 | 2.5e-1 | 1.1e-1 | 8.9e-4  | 2.8e-11 | **8.6e-17** |

The floor is the matrix-element precision (~ε_machine · ‖L‖_op ≈ 10⁻¹⁶);
that we hit it at r = 7 (or r = 8 at n=6) confirms `r_ref = 12` is far in
excess of any practical resolution.

## β_phys = 2 vs β_phys = 1 (smooth_def_s, side-by-side)

`r_D@1e-9` from the v2 recipe (β_phys ∈ {0.25, 0.5, 1.0}, dense BohrDomain
reference) against the new β_phys = 2 row:

| n | β_phys=0.25 | β_phys=0.5 | β_phys=1.0 | **β_phys=2.0** | Δ vs β_phys=1.0 |
|---|-------------|------------|------------|----------------|-----------------|
| 3 | 7 | 6 | 6 | **6** | 0 |
| 4 | 6 | 6 | 6 | **7** | +1 |
| 5 | 6 | 6 | 6 | **7** | +1 |
| 6 | 6 | 6 | 7 | **7** | 0 |

So the worst-case penalty across n is one extra bit. r_D = 7 is a safe
default for `default_smooth_s` from β_phys = 0.25 all the way to β_phys = 2.

What *does* shift: the intermediate per-r error is 2–5 orders of magnitude
worse at β_phys = 2 (e.g. n=6, r=6: 8.9e-4 at β_phys = 2 vs 1.2e-8 at
β_phys = 1). Convergence has steeper slope and starts further from zero,
so r below the cutoff is much less useful — but the final cutoff is
basically unchanged.

## Wall time (apples-to-apples, both matrix-free, JULIA_NUM_THREADS=8)

| n | cell wall, β_phys=1.0 | cell wall, β_phys=2.0 |
|---|------------------------|------------------------|
| 3 | (3.3s with dense L_b) | 6.2s   |
| 4 | —                      | 9.9s   |
| 5 | —                      | 42.0s  |
| 6 | **612s** [matrix-free both sides] | **292s** [matrix-free both sides] |

n=6 β_phys=2 is roughly half the wall of n=6 β_phys=1 in this
all-matrix-free setup — counter-intuitive but real. The Krylov power
iteration converges in ~half the iterations at β_phys=2 (15 ↔ 30) because
the spectral gap of `A^†A` opens up with the larger β_alg. So if anything,
β_phys = 2 is *cheaper* to qualify than β_phys = 1 at fixed n. For
context, the dense-L_b reference run at n=6, β_phys=1.0 in the v2
campaign was 113s (98s of which was the one-time L_b build), so the
all-matrix-free path adds substantial wall when both sides are big — but
unlike the dense path, it scales to n=7+ where 4 GB of dense L_b is
out of reach in the 3.8 GB sandbox.

## What's *actually* changing as β_phys grows

- β_alg = β_phys · rescaling_factor doubles → σ = 1/β_alg halves
- `s = (0.05/σ)² ∝ β_alg²` quadruples (smoothing parameter)
- Kink width in ν ≈ 1/(β·s) ∝ 1/β_alg³ shrinks 8×
- ω_range = 2(‖H‖ + 8σ) is essentially unchanged (‖H‖ dominates;
  rescaled-spectrum [0, 0.45] is n-independent)

So in principle resolving the kink needs `2^r ≫ ω_range · β_alg · s ∝ β_alg³`,
adding ~3 bits per β_phys doubling. Empirically we see only 0–1 bits.
The smooth-Metropolis filter must have effectively
band-limited spectral content much smaller than its kink width — the
trapezoid/Riemann sum on a smooth periodic function converges
super-polynomially once you're past the band-limit.

## Practical recommendation

For β_phys up to 2 and smooth_def_s, n ∈ {3..6}: **r_D = 7** uniformly
qualifies EnergyDomain as a 10⁻⁹ proxy for BohrDomain. Drop to r_D = 6 at
β_phys ≤ 1 if you want to save the bit — but at this point the cache cost
of r=7 over r=6 is fractions of a second so there is no reason not to
keep it pinned.

## Files

- Script: `scripts/scratch_quad_S1_bphys2_krylov.jl`
- BSON sidecars (β_phys=2): `scripts/output/quad_redo_S1_bphys2_krylov/n{3..6}_bphys2.0_smooth_def_s.bson`
- BSON sidecar (β_phys=1 apples-to-apples for n=6): `scripts/output/quad_redo_S1_bphys2_krylov/_n6_bphys1_reference/n6_bphys1.0_smooth_def_s.bson`
- Campaign summary: `scripts/output/quad_redo_S1_bphys2_krylov/_campaign_summary.bson`
- Beads: qf-w2u (extension of qf-yt9)
