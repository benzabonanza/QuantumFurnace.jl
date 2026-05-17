---
name: heisenberg-1d-even-odd-mechanism
description: Physics behind the even/odd-n split in 1D Heisenberg PBC + 0.1 Z+ZZ disorder. Even-n is frustration-free (unique singlet GS, Gibbs ≈ pure state at β≥1); odd-n has PBC-frustration multiplet at low energy that stays populated even at β=2.5.
metadata:
  type: project
---

# Even/odd-N mechanism for our 1D Heisenberg + 0.1 disorder PBC chain (2026-05-17)

The qualitative even/odd-n split documented in [[heisenberg-1d-multiseed-even-odd-qf-e4z-23]] and [[qf-e4z-31-v6-plus-parity-trap-quantified]] has a physical mechanism we worked out by directly inspecting Gibbs occupations at β_phys=2.5 (seed=42):

## Eigenstate structure

| n | parity | first 3 ΔE_phys above GS | p_0 at β_phys=2.5 | effective rank |
|---|---|---|---|---|
| 6 | EVEN | 2.69, 2.76, 2.87 | **99.7 %** | **1.01** |
| 8 | EVEN | 2.04, 2.10, 2.22 | **98.5 %** | **1.03** |
| 5 | ODD  | 0.063, 0.101, 0.157 | 30 % | 3.92 |
| 7 | ODD  | 0.012, 0.029, 0.143 | 28 % | 3.93 |

## Mechanism

**Even N (unfrustrated, frustration-free):** PBC AFM XXX with even N admits an alternating-singlet covering of the ring. Ground state is a unique total-singlet (S^z=0). Lowest excitation is Ω(1) above GS (~2 in physical units). At β_phys ≥ 1 the Gibbs state is essentially pure GS (p_0 > 80 %), so the relaxation problem reduces to "how fast does the slowest L mode bring you to the unique GS?" — the slowest mode is the GS↔first-excited coherence at a fixed Ω(1) Bohr frequency. KMS detailed balance forces this rate to scale roughly with the downhill `γ(−ΔE)` value, giving `gap_phys ≈ 1` in physical units. Empirically gap_phys drifts down slowly with N (1.27 at n=6, 1.08 at n=8) — consistent with gapless clean XXX spinons + weak disorder lifting (see [[qf-e4z-33-kdim-convergence]]). Whether even-N gap closes (1/N or sub-poly Griffiths) or saturates at a disorder floor is an open empirical question; the n ∈ {10, 12, 14} cluster sweep is needed to discriminate.

**Odd N (frustrated):** PBC AFM XXX with odd N **cannot** form a perfect alternating-singlet covering. S^z conservation forces minimum |S^z| = 1/2 → doublet (or larger) low-energy multiplet. Numerically we see 3–4 near-degenerate states within ΔE ≲ 0.15 in physical units, then a real Ω(1) gap to the next band. At β_phys = 2.5 the Gibbs is spread across this multiplet (p_0 ~ 28 %, eff_rank ~ 4 — NOT a pure state).

The Lindbladian gap on odd N is set by **intra-multiplet relaxation**: Bohr frequencies within the multiplet are tiny (0.01–0.15 phys), below the smooth-Metro kink kernel width σ·√s ≈ 1/(β_alg √4), so KMS rates are unsuppressed → relaxation within the multiplet is fast → `gap_phys ≈ 4–5`, roughly N-independent. The doublet structure is symmetry-protected at every odd N, so this Ω(1) gap should saturate (not close).

## Connection to the algorithmic parity-trap fix

The mechanism explains why [[qf-e4z-31-v6-plus-parity-trap-quantified]] only had to correct odd-N gaps (v6 → v6_plus shifted odd-n by 14–21 %, even-n unchanged):

* **Even N**: slowest mode is GS↔first-excited coherence, which sits in the parity-EVEN sub-spectrum (Gibbs is pure singlet, both states parity-even). `vec(I/d)` from v6 lives there → correct gap.
* **Odd N**: slowest mode is the `|↑⟩⟨↓|` intra-doublet coherence, which is parity-ODD (the two doublet states have opposite parity on odd N). `vec(I/d)` from v6 was trapped in parity-EVEN → reported the *parity-EVEN sub-spectrum* slowest mode (~5.0 phys) instead of the true parity-ODD slowest mode (~4.1 phys). `|+⟩⟨+|^⊗N` (qf-e4z.30 canonical) breaks parity, captures the true mode.

## Bibliographic context

No rigorous theorem in our literature corpus guarantees Ω(1) gap for our XXX at β=2.5. The closest results all require structural assumptions we don't satisfy:

* KB16 / Bardet 2023 / Kochanowski 2024 / Capel 2021 / Capel 2025 / Bergamaschi-Chen-Liu 2024: all **commuting** (ours is non-commuting).
* Rouzé 2024: noncommuting but only for **β < 1/(615^D · J) ≈ 0.0016** in D=1.
* Tong-Zhan 2024 / Smid 2025: **fermionic** weakly-interacting (ours under JW gives strong interactions).
* Ramkumar-Soleimanifar 2024: **non-local 1-design jumps** + β‖H‖ = O(1).

The "low-T × noncommuting × natural-local × spin" regime is the explicit open "punchline cell" of the literature synthesis. Our XXX + 0.1 disorder + PBC at β=2.5 sits squarely in that cell; the physics-motivated heuristic above gives the most honest expectation.

## Thesis presentation

When we plot 1D Heisenberg τ_mix / gap_phys sweeps, **always plot even-N and odd-N separately** and explain the mechanism in the figure / accompanying text. Pooling even and odd into a single scaling fit is wrong (the qf-e4z.23 v6 pooled M0 fit had only 0.69 model weight precisely because of this misspecification). Per [[feedback-more-data-points-for-scaling-claims]] we will not extract asymptotic scaling exponents from 3-point sweeps — describe the trend qualitatively.
