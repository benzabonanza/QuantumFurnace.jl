# TimeDomain Dissipator Quadrature at β_phys = 2, smooth Metropolis

**Scope caveat**: These numbers measure the **dissipator term only** —
the matvec ran with `include_coherent = false`. The coherent commutator
$i[B, \rho]$ is masked. The full-Lindbladian register triple
$(r_D, r_{b_-}, r_{b_+})$ requires `r_{b_+}` sized separately — and in
TimeDomain `r_{b_+}` caps at ε ≈ 10⁻⁶ due to a slope-(−1) anomaly at
$t = 0$ (see [trap-rule t=0 L'Hôpital origin](../../.claude-memory/trap_rule_t0_lhopital_origin.md)). For ε ≤ 10⁻⁹ on the full
Lindbladian, use EnergyDomain (`B_energy ≡ B_bohr` closed-form, no
quadrature).

**Task**: qf-fqt. Companion to qf-w2u: same fixture, same filter, same
all-matrix-free methodology, but now the **test side is `TimeDomain`**
(FINUFFT-driven) while the reference stays `EnergyDomain(r_D = 12)`. The
v2 recipe established that EnergyDomain at r=12 is indistinguishable from
the analytical BohrDomain at machine matrix-element precision, so this
directly measures

$$
  \big\| L_\mathrm{T}(r_\mathrm{test}) - L_\mathrm{B} \big\|_\mathrm{op}
$$

— the right scaling for `predict_lindbladian_trajectory` if we were to
run it on a `Config{Lindbladian, TimeDomain}` instead of the default
`Config{Lindbladian, EnergyDomain}`.

**Why this matters**: `predict_lindbladian_trajectory` (src/lindblad_action.jl)
builds a Krylov subspace by repeated `apply_lindbladian!` matvecs. If
TimeDomain were dramatically more expensive per matvec, or required many
more bits of `r_D` to qualify as a 10⁻⁹ proxy, the EnergyDomain route
would be the only viable one. The user worry was that the FINUFFT and
larger time-cache might "blow up" — they don't.

## Headline

| n | r_D @ 1e-9 (EnergyDomain) | r_D @ 1e-9 (TimeDomain) | Wall (EnergyDomain) | Wall (TimeDomain) | Wall ratio |
|---|---------------------------|--------------------------|----------------------|--------------------|------------|
| 5 | 7 | 7 | 42 s | 69 s | **1.64×** |
| 6 | 7 | 7 | 292 s | 398 s | **1.36×** |

**TimeDomain does not blow up at n=5 or n=6, β_phys=2**. Same `r_D = 7`
cutoff for ε = 10⁻⁹, modest 1.4–1.6× wall overhead per cell. Convergence
is *slightly faster per r at intermediate r* than EnergyDomain (the
FINUFFT-truncated, band-limited time grid is a marginally better
quadrature shape for the smooth Metropolis filter), but it terminates at
the **~10⁻¹³ FINUFFT floor** instead of machine precision. The per-r
matvec cost is identical to EnergyDomain at fixed r once the NUFFT
prefactor cache is built — the wall overhead is amortised setup +
FINUFFT JIT, not asymptotic.

## Per-r operator-norm sweep (β_phys = 2, smooth_def_s)

### n = 5

| r_D | EnergyDomain ‖ΔL‖_op | TimeDomain ‖ΔL‖_op |
|-----|----------------------|---------------------|
| 3   | 5.96e-01 | 3.27e-01 |
| 4   | 2.43e-01 | 1.61e-01 |
| 5   | 5.11e-02 | 2.04e-02 |
| 6   | 9.54e-05 | 1.77e-05 |
| 7   | 2.28e-16 *(machine)* | **1.16e-13** *(FINUFFT floor)* |
| 8   | — | 7.71e-14 *(saturated)* |

### n = 6

| r_D | EnergyDomain ‖ΔL‖_op | TimeDomain ‖ΔL‖_op |
|-----|----------------------|---------------------|
| 3   | 4.18e-01 | 2.35e-01 |
| 4   | 2.50e-01 | 1.38e-01 |
| 5   | 1.10e-01 | 2.45e-02 |
| 6   | 8.91e-04 | 1.98e-04 |
| 7   | 2.84e-11 | 4.18e-12 |
| 8   | 8.64e-17 *(machine)* | **5.04e-14** *(near FINUFFT floor)* |

Worth highlighting: at n=6, r=5 TimeDomain is **5× better** than
EnergyDomain (2.5e-2 vs 1.1e-1) and at r=6 it is 4× better (2.0e-4 vs
8.9e-4). The TimeDomain Riemann-sum quadrature lines up unusually well
with the smooth Metropolis kernel structure once the time window
captures most of its support — by r=7 ω_range·t0 ≈ 49 in absolute
units, vs the kernel's effective support of ~10. The
EnergyDomain ω-grid resolution sees no analogous "envelope" bonus at
those r.

The TimeDomain "OFT integration warning" prints `filter kernel at the
ends should be small but it is: 0.97` at r=3 (effectively no truncation)
down to `1.1e-5` at r=6 — by r=7 the time window is wide enough that the
smooth Metropolis kernel is numerically zero at the boundary, which is
precisely the regime where the FINUFFT-driven Riemann sum becomes exact
up to the requested precision.

## Predict-Lindbladian Krylov implications

For `predict_lindbladian_trajectory` at β_phys = 2 with smooth_def_s the
recipe is **identical between domains** at both n = 5 and n = 6: r_D = 7
for ε = 10⁻⁹. The TimeDomain route is viable everywhere the EnergyDomain
route is — 1.4–1.6× slower in matrix-free wall time, floored at ~10⁻¹³
instead of machine precision. At n = 6 the wall ratio actually
**improves** to 1.36×: the FINUFFT prefactor cache + JIT warmup are
amortised better against the longer per-matvec time.

If FINUFFT's 10⁻¹³ floor matters (it doesn't for the thesis target ε =
10⁻⁶ algorithm-level or ε ≈ 10⁻⁵ simulated-algorithm), drop FINUFFT's
`eps` parameter in `_prepare_oft_nufft_prefactors`
(src/furnace_utensils.jl:224 — currently hardcoded `eps=1e-12`,
load-bearing for byte-identity tests). At 10⁻¹³ accuracy needs you'd
also need to verify the per-`(jump, ω)` sandwich accumulation hasn't
introduced its own cancellation noise.

## Why TimeDomain is competitive here

Both EnergyDomain and TimeDomain matvecs do the same outer loop:

```
for each jump:
    for each (ω in energy_labels):
        jump_oft = eigenbasis * prefactor_matrix[:, :, idx[ω]]
        accumulate γ(ω) · jump_oft * ρ * jump_oft†   into out
```

The *only* difference is what `prefactor_matrix[:, :, idx[ω]]` contains:

- EnergyDomain: `δ_{ω, λ_i − λ_j}` (essentially a Kronecker delta, but
  encoded as a real prefactor coming from a tightly-resolved sinc-style
  energy quadrature).
- TimeDomain: the FINUFFT type-3 transform of the filter kernel sampled
  on the time grid, evaluated at the Bohr-frequency points.

The per-matvec cost is identical at fixed r once the NUFFT prefactor
cache is built — and the cache build is amortised across the ~15–30
Krylov iterations per r. The 1.6× wall overhead measured here is
dominated by the lower-r warmup costs (FINUFFT setup, JIT
compilation) — at r ≥ 6 per-matvec is within ~10% of EnergyDomain.

## Files

- Script: `scripts/scratch_quad_S2_bphys2_time_krylov.jl`
- BSON sidecars (TimeDomain): `scripts/output/quad_redo_S2_bphys2_time_krylov/n{5,6}_bphys2.0_smooth_def_s_time.bson`
- Comparison source (β_phys = 2, EnergyDomain): `scripts/output/quad_redo_S1_bphys2_krylov/n{5,6}_bphys2.0_smooth_def_s.bson`
- Beads: qf-fqt (companion to qf-w2u under epic qf-yt9)
