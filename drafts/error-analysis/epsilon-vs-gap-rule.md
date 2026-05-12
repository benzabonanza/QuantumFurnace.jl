# ε vs spectral gap — when Krylov gap/τ_mix estimates are reliable

**Topic**: The "ε" in any register/quadrature sizing is shorthand for
$\|\mathcal L_\text{disc} - \mathcal L_\infty\|_\text{op}$ — the sum of
*all* quadrature errors. For a Krylov gap or τ_mix predictor to deliver
results precise to relative tolerance $\delta_\tau$, that operator-norm
error must satisfy $\varepsilon \lesssim \delta_\tau \cdot |\lambda_1|$.
This document records the supporting evidence and a single concrete
data point that nails it numerically.

**Origin**: qf-1au speed bench (TimeDomain `predict_lindbladian_trajectory`
at ε=10⁻³ vs ε=10⁻⁶) showed a 16% τ_mix shift at n=6, β_phys=2 between
the two register sizings, which the qf-f45 split coherent sweep traced
to the β-amplified TimeDomain coherent-term floor.

## The principle

For any Krylov-based gap predictor:

$$
  \frac{\big|\delta\lambda_1\big|}{|\lambda_1|} \;\le\;
  \kappa\,\frac{\|\Delta\mathcal L\|_\text{op}}{|\lambda_1|},
$$

with $\kappa$ the eigenvector condition number (modest in practice for
KMS-DB Lindbladians but not strictly 1, since these operators are
non-normal). The Krylov decomposition itself reads the spectrum of the
operator we built to machine precision (`eigen` on the small Hessenberg
+ `krylovdim` ≥ 2 captures the gap). So *all* of the gap error is in
the operator-norm distance from the continuum, not in the eigensolve.

That operator-norm distance is the *sum* of every quadrature
contribution (dissipator, coherent, FINUFFT precision, Trotter $M$ for
the Trotter path, δ-step for the channel path). The label "ε" in any
single register table is just a stand-in for this sum.

So the operational rule:

> **Pick register/discretisation precision relative to the smallest gap
> you intend to resolve.** Don't pick "ε" first and hope it's enough.

| Target τ_mix accuracy | Required ε ≈ |
|------------------------|---------------|
| 10% | $\,10^{-1} \cdot \min|\lambda_1|$ |
| 1%  | $\,10^{-2} \cdot \min|\lambda_1|$ |
| 0.1%| $\,10^{-3} \cdot \min|\lambda_1|$ |

## Verified data point: n=6, β_phys=2, smooth_def_s

Full Lindbladian with coherent term ON, `predict_lindbladian_trajectory`
via `_krylov_spectral_decomposition`, krylovdim = 30. The "TRUE" gap
comes from EnergyDomain with B_energy ≡ B_bohr exact (no coherent
quadrature) and r_D = 8 (dissipator quadrature ≤ 10⁻¹¹ ≪ gap):

| Setup | r_D | r_b_minus | r_b_plus | actual ‖ΔL‖_op | gap | shift vs TRUE | τ_mix |
|-------|-----|-----------|----------|----------------|-----|---------------|--------|
| **EnergyDomain (TRUE)** | 8 | — | — | ~10⁻¹¹ | **3.974×10⁻³** | — | **251.6** |
| TimeDomain "ε=10⁻⁶" reg | 7 | 7 | 17 | ≈ 2.4×10⁻⁵ (coherent floor) | 4.060×10⁻³ | **+2.17%** | 246.3 |
| TimeDomain "ε=10⁻³" reg | 6 | 6 | 7  | ≈ 1.2×10⁻³ (coherent r_+ tail) | 3.398×10⁻³ | **−14.48%** | 294.3 |

Sources:
- `scripts/output/quad_redo_S34_coherent_bphys2/n6_bphys2.0_smooth_def_s.bson` — TimeDomain B_time vs B_bohr split sweep
- `scripts/output/predict_lindbladian_time_bench/n6_bphys2.0_smooth_def_s_time.bson` — TimeDomain ε=10⁻⁶ bench
- `scripts/output/predict_lindbladian_time_bench_eps1e3/n6_bphys2.0_smooth_def_s_time.bson` — TimeDomain ε=10⁻³ bench
- `scripts/scratch_predict_lindbladian_energy_spotcheck.jl` — EnergyDomain TRUE gap

### Reading the row-by-row consistency

- TimeDomain ε=10⁻³: ‖ΔB‖ ≈ 1.15×10⁻³, observed Δgap = −5.76×10⁻⁴.
  Ratio |Δgap|/‖ΔB‖ ≈ 0.5 — within the Bauer-Fike-with-κ bound (and
  the sign of the gap shift matches: under-resolved B → spectrum
  compressed toward 0 → smaller gap).
- TimeDomain ε=10⁻⁶: ‖ΔB‖ ≈ 2.4×10⁻⁵ (at the TimeDomain coherent
  floor), observed Δgap = +8.6×10⁻⁵. Ratio ≈ 3.6 — same order, picks
  up some κ.

So the perturbation bound is realised at order-1 constants in practice.

## Why TimeDomain hits a floor here

The B_time coherent quadrature for smooth_def_s has a slope-(−1) tail
in 1/$N_+$ from the $t = 0$ L'Hôpital sample of the Cauchy P.V.
integrand. The residual scales linearly in $\beta_\text{alg}$ from the
$1/t$ envelope, and as $\beta_\text{alg}^2$ from the smooth_def_s
parameter $s = (0.05/\sigma)^2$. So between β_phys = 1 and β_phys = 2
(β_alg doubling at fixed n=6: 56 → 112) the floor scales by roughly
$2 \times 4 = 8$, with empirical observation:

| β_phys | β_alg (n=6) | TimeDomain coherent floor ‖ΔB‖ |
|---------|--------------|--------------------------------|
| 1.0    | 56  | ~2×10⁻⁶ (v2 recipe) |
| 2.0    | 112 | **~2.4×10⁻⁵** (this work) |

Pushing r_b_plus past ~13 buys nothing at β_phys = 2 (verified up to
r_b_plus = 20). The floor is *physical to the discretisation*, not a
register issue.

In EnergyDomain the coherent term is the analytical $B_\text{bohr}$
closed-form — there is no quadrature at all. The only contribution
to $\|\Delta\mathcal L\|_\text{op}$ is the dissipator, which converges
super-exponentially in $r_D$: at n=6, β_phys=2 the EnergyDomain
dissipator gives ε_diss ≈ 8.9×10⁻⁴ at r_D = 6 but **≈ 2.8×10⁻¹¹ at
r_D = 7** (qf-w2u). One extra bit ⇒ 8 orders of magnitude. So
EnergyDomain at r_D = 7 essentially eliminates all spectral error at
any cell we've measured.

## Operational guidance

1. **For ideal Lindbladian τ_mix plots** (thesis P1 / S1): use
   `Config{Lindbladian, EnergyDomain, KMS}` with r_D ≥ 7 across the
   whole grid. Coherent is automatically exact. The only knob is r_D,
   and 7 reaches 10⁻⁹ to 10⁻¹¹ on the dissipator at every measured
   cell — i.e., $\varepsilon \ll$ any conceivable gap in the
   parameter range.

2. **For implemented-channel τ_mix plots** (compare δ-step CPTP
   channel mixing vs ideal Lindbladian mixing): the channel has its
   own (δ, gqsp_degree, M_user) error budget that contributes to
   $\|\Delta\Phi_\delta\|$ on top of any register quadrature. The
   relevant question is whether the channel's *intended* error budget
   (e.g. δ = 10⁻³, δ²/λ ≈ 10⁻³ accumulated) is smaller than the gap
   you want to resolve. Same rule, same threshold formula.

3. **First pass for any new sweep**: do a cheap probe at coarse
   registers, read off the *minimum* gap across cells, then pick the
   production register sizing to satisfy ε ≤ 10⁻² · min(gap). If
   that's incompatible with the cheap-register choice (most relevant
   at low T / large n where the gap closes), switch domain
   accordingly.

4. **The τ_mix you measure is the τ_mix of whatever operator you
   built**, not of the continuum L_∞. The two coincide only when
   $\|\Delta\mathcal L\|_\text{op} \ll |\lambda_1|$.

## Pointers

- v2 recipe register tables: `drafts/error-analysis/quadrature-convergence-summary-v2.md`
- t=0 L'Hôpital origin of the TimeDomain coherent floor: [[trap_rule_t0_lhopital_origin]] in memory
- TimeDomain/EnergyDomain Krylov routes: [[krylov_two_routes]] in memory
- include_coherent rule: [[feedback_coherent_term_on_by_default]] in memory
- Beads: qf-1au, qf-f45 (closed); under epic qf-yt9 thread
