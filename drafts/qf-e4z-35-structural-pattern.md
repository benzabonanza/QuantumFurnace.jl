# qf-e4z.35 σ-sweep — structural mixing pattern (n=3..7)

**Question:** at what σ value does the smooth-Metropolis CKG sampler reach
its **best structural mixing** — once we factor out the σ-dependence of
the generator's overall rate scale?

**Setup:** the qf-e4z.35 σ-sweep produces per-cell `mixing_time_phys`,
`gap_phys`, `‖L‖_HS_phys`, and `d_{1→1}_phys` (Kossakowski opnorm bound
on the dissipator). At fixed (n, β_phys), normalising `τ_mix` by any
L-norm `N` measures mixing in "L-time-units" — `τ_mix · N` removes pure
rate-scale changes, leaving the structural component.

Analysis script: `scripts/scratch_qf_e4z_35_structural_pattern.jl`.
Full tabular output: `drafts/qf-e4z-35-structural-pattern-n3-7.txt`.

## c* (argmin of τ_mix · ‖L‖) tables

Both L-norms agree to within ±1 grid step everywhere — the structural
optimum is robust under norm choice.

**τ_mix · ‖L‖_HS:**

| β_phys \\ n | n=3 | n=4 | n=5 | n=6 | n=7 |
|------------|-----|-----|-----|-----|-----|
| 0.25 (hot) | 0.25 | 0.75 | 0.75 | 0.75 | 0.75 |
| 0.5        | 0.25 | 1.50 | 1.00 | 1.50 | 1.50 |
| 1.0 (cold) | 0.25 | 2.00 | 1.50 | 2.00 | 2.00 |

**τ_mix · d_{1→1}:**

| β_phys \\ n | n=3 | n=4 | n=5 | n=6 | n=7 |
|------------|-----|-----|-----|-----|-----|
| 0.25 (hot) | 0.25 | 0.75 | 0.75 | 0.75 | 0.75 |
| 0.5        | 0.25 | 1.50 | 1.00 | 1.50 | 1.50 |
| 1.0 (cold) | 0.25 | 2.00 | 1.50 | 1.50 | 1.50 |

## Findings

1. **n=3 is anomalous** — every β_phys cell prefers c*=0.25. Likely because
   the d=8 spectrum is too small / dense for a wider smoothing kernel
   to be helpful: every nontrivial Bohr frequency is already
   well-separated, so the narrowest kink wins.

2. **For n ≥ 4: optimal σ grows with β_phys** (colder ⇒ wider σ is
   structurally better):
   - β=0.25:  c* ≈ 0.75
   - β=0.5:   c* ≈ 1.0–1.5
   - β=1.0:   c* ≈ 1.5–2.0 (often at the upper edge of our c-grid — the
     true optimum may lie at c > 2 for β=1)
   The optimum is essentially **n-independent** within each β bin —
   this is a temperature effect, not a finite-size effect.

3. **Structural spread is modest:** the max/min of τ_mix·‖L‖ across the
   c-grid is **1.25–1.7×** at n ≥ 4. Picking the right σ buys roughly
   25–70% structural improvement, not orders of magnitude.

4. **The canonical σ = 1/β_alg (c=1) is well-placed at β=0.5 but
   sub-optimal at both temperature extremes** — slightly too wide at
   β=0.25 (optimum c=0.75) and too narrow at β=1.0 (optimum c=1.5–2.0).

## Why the two normalisations agree

`‖L‖_HS` measures the full generator (dissipator + coherent commutator),
while `d_{1→1}` is a Wolf–Pérez-García style bound on the dissipator
only. They agree on c* because (a) the coherent part is σ-independent
(it's the Hamiltonian commutator i[H + B, ·] with B coming from a
σ-dependent integrand that integrates out to a near-σ-independent
operator), and (b) the dissipator's σ-dependence dominates the L-norm
σ-scaling at the (n, β) cells we sweep.

## Different metrics, different answers

- **τ_mix (absolute mixing time):** prefers c=0.25 everywhere because
  smaller σ ⇒ larger ‖L‖ ⇒ faster wall-clock mixing. This is the
  rate-scale effect we want to factor out.
- **τ_mix · ‖L‖_HS / τ_mix · d_{1→1} (structural):** see above —
  optimum migrates from c=0.75 (hot) to c=1.5–2.0 (cold) for n ≥ 4.
- **τ_mix · gap (rate-resolved):** picks c*=2.0 almost uniformly. This
  measures the |c_2| coefficient of the |+⟩⟨+|^⊗N initial state on the
  slowest mode — at wider σ, that overlap shrinks, so τ_mix per gap
  drops. Different question from "is L structurally better".
- **gap / ‖L‖_HS (relative spectral gap):** peaks at c ≈ 0.25–0.75 at
  n ≥ 4 — smaller σ ⇒ slow mode more separated from the bulk in
  relative terms — but this isn't what determines τ_mix (because |c_2|
  matters for the initial state, not just the gap).

So **for a thesis claim of "best structural mixing"**, the L-norm
normalisations (HS or d_{1→1}) are the right metric — they answer the
question "given a unit-norm budget on L, which σ gives the fastest
trace-distance decay?". The answer:

> For the canonical 1D Heisenberg PBC fixtures at canonical β_phys ≤ 1
> and n ≥ 4: the structural optimum is a wider-than-canonical σ that
> grows with β_phys (from c ≈ 0.75 at the hottest cell to c ≈ 2 at the
> coldest), nearly system-size independent.

## What's next (when n=8 resumes)

The n=8 cells will pin down whether the n-independence of c* persists
into the larger system. The β=1 optimum at c=2.0 is at the grid edge —
a follow-up sweep extending the c-grid (e.g. c ∈ {2, 3, 4}) at β=1
would tell us whether the cold-cell optimum sits inside the grid or
keeps growing.
