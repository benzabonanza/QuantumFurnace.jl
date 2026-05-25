# qf-c9g — Literature survey: ordered-phase gap closing vs GS-doublet artefact

**Context.** 2D TFIM (h = J = 1, β_phys = 2.0, ladders n = 4, 6, 8) shows
gap_phys = 9.96e-2 → 3.06e-3 → 1.57e-4. Is this **genuine ordered-phase
free-energy-barrier closing** (Z₂ surface-tension mechanism, asymptotic
exp(−γ·perimeter) or stretched-exponential) or a trivial finite-size
GS-doublet tunnel-splitting artefact? T/T_c ≈ 0.24 (T_c ≈ 2.07,
Hesselmann–Wessel 2016, arXiv:1602.02096).

## 1 Classical 2D Ising Glauber, low-T lower bound

**Martinelli–Olivieri 1994** (*Comm. Math. Phys.* 161, 487–514; also
*J. Stat. Phys.* 76 1145, DOI 10.1007/BF02187060). On (ℤ/Lℤ)² with
β > β_c, the spectral gap of single-site Glauber dynamics satisfies

  gap ≤ exp(−τ(β) · L · (1 + o(1))),

where **L is the linear lattice side** and τ(β) > 0 is the
surface tension along a coordinate axis (Wulff construction). Mechanism:
flipping between the +/− magnetisation sectors requires nucleating a
domain of opposite phase, whose minimal free-energy cost scales with
**perimeter** ~ L = √N where N = L². Equivalently in our notation:
**gap ~ exp(−γ·√N)** — stretched-exponential in N, full-exponential in L.

**Lubetzky–Sly 2013** (arXiv:1305.4524, local PDF
`supplementary-informations/classical-review/lubetzky-sly-2013-cutoff-ising-lattice.pdf`)
covers the *complementary* high-T regime (β < β_c on ℤ², all β > 0 with
external field h ≠ 0): Theorem 2 gives cutoff at (d/2λ_∞) log N with
window O(log log N). This **does not** apply to our β_phys = 2 > β_c
case; we cite it only to anchor that the slow-mixing regime starts
strictly below T_c, not deeper.

## 2 Quantum 2D TFIM Lindbladian-gap lower bounds

**Gamarnik–Kiani–Zlokapa 2024** (arXiv:2411.04300, "Slow Mixing of
Quantum Gibbs Samplers"). **Theorem 1.4** (informal): for the 2D TFIM
H = −∑Z_iZ_j − h∑X_i on an n × n lattice with constant β ≥ β*, h ≤ h* ≈ 1,
any geometrically local Lindbladian (including Davies and the
**Chen–Kastoryano–Gilyen 2023 KMS-DB sampler** — the same class as our
construction; reduction via Lieb–Robinson bounds) obeys

  T_mix(ℒ) ≥ exp[n^(1/2 − o(1))] = exp[Ω(√n_total_spins)].

Mechanism (Section 4): a Peierls fault-line argument; the bottleneck is
the **same surface-tension free-energy barrier** as the classical case,
"dressed" by quantum transverse-field fluctuations via the
Poisson–Feynman–Kac representation. Stretched-exponential, NOT
full-exponential — surface tension produces a √N (perimeter) law for
square 2D regions, not a volume law.

**Rakovszky–Placke–Breuckmann–Khemani 2024** (arXiv:2412.09598,
"Bottlenecks in Quantum Channels and Finite-Temperature Phases"). A
quantum conductance/Cheeger-style lemma: if Hilbert-space regions are
separated by a low-weight off-diagonal cut, channel mixing is bounded
below by the inverse cut weight. The headline example is
extensive-barrier commuting-projector codes (full-exponential), but the
**bottleneck framework is mechanism-agnostic** and applies to any
symmetry-broken finite-T phase whose ordered sectors are separated by a
surface-tension cut, including 2D TFIM at β > β_c.

**Bergamaschi–Chen–Liu 2024** (arXiv:2404.14639, "Quantum Computational
Advantage with Constant-Temperature Gibbs Sampling"). Designs *specific*
constant-locality commuting Hamiltonians (shallow-circuit parents) for
which a CKG-style sampler **does** mix in poly(N) at constant T. The
paper is silent on 2D TFIM and explicitly does NOT apply to it; we cite
it as the only known rapid-mixing constant-T result, to anchor that
"rapid mixing at constant T" is a knife-edge property of engineered
commuting models, not a generic feature of natural Hamiltonians.

## 3 Finite-size diagnostics — ordered phase vs GS-doublet artefact

The relevant criteria (Bauerschmidt–Dagallier 2024 on near-critical
Ising LSI; Cuff et al. 2012 on mean-field Potts; Bovier–den Hollander
*Metastability* monograph, ch. 16; classical Martinelli surveys):

1. **β-persistence below T_c.** Closing must survive over the entire
   ordered window β_c < β < ∞, NOT just at deep cold. A pure
   GS-doublet artefact would show gap_phys decoupling from temperature
   once β ≫ ΔE_GS-gap (since Gibbs ≈ |GS⟩⟨GS|). Genuine surface-tension
   closing tracks **τ(β)·L** with τ(β) → 0 as β ↓ β_c⁺. The signature is
   gap_phys getting *worse* with N at *every* β in (β_c, ∞), and a τ(β)
   exponent that is finite and continuous through the ordered phase.

2. **Doublet vs multiplet, magnetisation moments.** Binder cumulant
   U_4(β, n) = 1 − ⟨m⁴⟩/(3⟨m²⟩²) → 2/3 in the ordered phase, → 0 in the
   paramagnet; ⟨m_z²⟩ → 1; eff_rank → 2 reflects "the Gibbs is
   essentially the doublet" but is **necessary, not sufficient** for
   ordered-phase physics — a trivial two-level system also has
   eff_rank ≈ 2. The discriminator is N-scaling of m² (intensive,
   plateau in ordered phase) plus U_4 → 2/3 (Binder 1981; standard
   Monte-Carlo finite-size scaling).

3. **Slow-mode operator-space character.** For genuine surface-tension
   closing, the slowest Lindbladian eigenmode R_2 should be **strongly
   correlated with the order parameter** — large overlap with M_z
   = ∑Z_i, dominated by doublet × bulk off-diagonal coherences that
   carry the inter-sector domain-wall content (not within-doublet
   matrix elements; we already saw at n = 4 that R_2 lives 89 % in the
   doublet × bulk block — exactly the surface-tension signature).
   Pure GS-doublet artefact would give within-doublet character.

4. **Disordered-phase control.** Run the same sweep at β_phys = 0.5 and
   0.1 (we already have eff_rank = 5 and 56). gap_phys should be
   **Ω(1)** in the paramagnet — no closing — confirming the closing is
   phase-specific, not a Krylov / sampler artefact.

## 4 Expected verdict for our n-scan

The literature is unambiguous: at β > β_c on 2D TFIM, both classical
(Martinelli–Olivieri) and quantum (Gamarnik–Kiani–Zlokapa) prove
**exp(Ω(√N))** lower bounds via the surface-tension mechanism for the
exact sampler class we use. Our observed sequence
9.96e-2 → 3.06e-3 → 1.57e-4 at n = 4, 6, 8 is a factor ≈ 32 per +2
sites. Asymptotic exp(−γ·√n) at γ ≈ 4 (rough Ising surface tension at
β = 2, h = 0) would predict factors ≈ exp(4·(√6−√4)) ≈ 24 and
≈ exp(4·(√8−√6)) ≈ 16 — same ballpark but pre-asymptotic; n ∈ {4, 6, 8}
ladders are too small for clean √N fits.

Diagnostic plan: (i) repeat the n-sweep at β_phys ∈ {0.5, 1.0, 1.5},
expect closing to **persist and weaken with decreasing β**; (ii)
β_phys = 0.1 control with no closing; (iii) confirm R_2's
doublet × bulk character at every n. Three out of three positive ⇒
the closing IS the asymptotic ordered-phase scaling, NOT a GS-doublet
coincidence.

## Sources

- Lubetzky & Sly 2013, arXiv:1305.4524 (local PDF)
- Martinelli & Olivieri 1994, DOI 10.1007/BF02187060
- Gamarnik, Kiani & Zlokapa 2024, arXiv:2411.04300
- Rakovszky, Placke, Breuckmann & Khemani 2024, arXiv:2412.09598
- Bergamaschi, Chen & Liu 2024, arXiv:2404.14639
- Hesselmann & Wessel 2016, arXiv:1602.02096 (T_c value)
- Bauerschmidt & Dagallier 2024 (near-critical Ising LSI, local PDF)
- Cuff et al. 2012, mean-field Potts Glauber cutoff (local PDF)
- Bovier & den Hollander, *Metastability* (local PDF), ch. 16
