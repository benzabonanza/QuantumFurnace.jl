---
name: Thesis Numerics Plan
description: Plan for generating Ch5 numerical data — plots, cross-checks of analytical results (Props 5,7,9,10, Kossakowski, quadrature, Trotter)
type: project
---

# Numerical Simulation Plan for Chapter 5

Goal: generate all plot data and cross-check analytical statements from the core Chapter 5 (Dissipative Quantum Gibbs Sampling). Results feed into thesis figures and validate the original contributions.

## A. Transition Weight Comparison Plots

1. **γ(ω) function plots**: Gaussian γ_G, Metropolis γ_M, smooth Metropolis γ_M^(s) for various s values (s=0.1, 0.5, 1, 2). Show kink at ω = -βσ²/2 disappearing with s > 0.
2. **Effective temperature**: plot β_eff vs σ (eq 5.15), confirming β_eff < β for σ > 0 and β_eff → β as σ → 0.

## B. Kossakowski Matrix Analysis

3. **α_{ν1,ν2} heatmaps**: for Gaussian, Metropolis, smooth Metropolis weights. Show off-diagonal structure and skew-symmetry (eq 5.22).
4. **Diagonal KMS condition check**: verify α_{ν,ν}/α_{-ν,-ν} = e^{-βν} numerically for all three γ choices.
5. **Corollary 8 (GNS↔KMS shift)**: verify α^GNS_{ν1,ν2} = α^KMS_{ν1-βσ²/2, ν2-βσ²/2} numerically.
6. **Positive semidefiniteness**: check eigenvalues of α matrix are non-negative for all cases.

## C. Spectral Gap & Mixing Time

7. **Proposition 5 (Metropolis optimality)**: for CKG Lindbladians, compute Gap(L^{γ_M}) and Gap(L^{γ}) for various γ choices; verify γ_M always gives largest gap.
8. **Spectral gap vs s**: sweep s ∈ [0, 2], plot gap of CKG Lindbladian. Confirm gap decreases with s but only mildly for small s.
9. **Spectral gap vs β**: temperature dependence of mixing for different transition weights.
10. **Mixing time vs system size n**: for n=3,4,5 (and 6 if feasible).

## D. Smooth Metropolis Gevrey Properties (Proposition 7)

11. **Derivative bound verification**: compute |(γ_M^(s))^{(n)}(ω)| numerically (finite differences or symbolic) and compare against the Gevrey-1/2 bound C1 · C2^n · n^{n/2} for various n and s.
12. **Gevrey constants**: verify C1=1, C2=1/(σ√(es)) match the numerical envelope.

## E. Quadrature Error Analysis

13. **Table 5.1 reproduction**: for each γ choice, compute actual quadrature error |α_{ν1,ν2} - ᾱ_{ν1,ν2}| as function of ω0 (grid spacing). Confirm:
    - Gaussian: O(e^{-c_G/ω0²})
    - Smooth Metropolis: O(e^{-c_s/ω0²}) with c_s ≤ c_G
    - Metropolis: O(β ω0²)
14. **Estimating qubits r vs 1/ε**: polylog for Gaussian and smooth Metropolis, polynomial for Metropolis (eq 5.60 vs the Metropolis scaling).
15. **Coherent term B quadrature**: outer (b-) and inner (b+) discretization errors separately.

## F. Time-Domain Kernels

16. **b-(t) plot**: universal kernel (eq 5.40), verify shape and ℓ1 norm bound.
17. **b+(t) plots**: for Gaussian (eq 5.41) and smooth Metropolis (eq 5.43). Show regularization near t=0 and effect of η.
18. **Subnormalization ‖b-‖₁ · ‖b+‖₁**: verify O(β log(1/ε)) scaling for Metropolis-like cases.

## G. Trotterization Error (Propositions 9, 10)

19. **Proposition 9 (Trotter for L_diss)**: compute ‖L̄_diss - L̃_diss‖_{1→1} for various M (Trotter steps), verify 1/M² scaling with Strang splitting and the √15 · α̃_comm / (M² σ³) prefactor.
20. **Proposition 10 (Trotter for B)**: compute ‖B̄ - B̃‖ for various M, verify the two-term structure (outer σ⁻³ and inner β³ contributions).
21. **Commutator-scaling constant α̃_comm^(2)**: compute explicitly for Heisenberg chain Hamiltonians at n=3,4,5.
22. **Palindromic vs non-palindromic (Remark 11)**: compare anti-Hermitian part of discriminant for Strang (palindromic) vs Lie-Trotter (non-palindromic).

## H. Detailed Balance Verification

23. **Anti-Hermitian part ‖A(ρ,L)‖_{2→2}**: the main DB violation metric (Definition 1). Compute for:
    - Continuous L (should be 0 for CKG)
    - Discretized L̄ (quadrature-limited)
    - Trotterized L̃ (quadrature + Trotter)
24. **Fixed-point accuracy (Corollary 2)**: verify ‖ρ_fix(L) - ρ_β‖₁ ≤ 20 · t_mix · ε numerically.

## I. Generator Splitting

25. **Jump-wise splitting error**: verify ‖Φ_A - e^{δL}‖ = O(δ² log²(β‖H‖)) (eq 5.86).
26. **BCH anti-symmetric structure**: confirm leading splitting error is KMS anti-self-adjoint (eq 5.89), so spectral gap agrees to O(δ³).
27. **Coherent-dissipative splitting**: confirm fixed-point shift of O(t_mix · δ) (eq 5.92).

## Priority Order
1. Start with B, C (Kossakowski + spectral gap) — these are the most impactful for the thesis narrative.
2. Then E (quadrature errors) — validates the smooth Metropolis advantage.
3. Then G, H (Trotter + DB) — validates implementation correctness.
4. Then A, D, F (plots, Gevrey verification, kernels) — supporting figures.
5. Finally I (generator splitting) — less critical but good to have.

## Notes
- Do NOT assume that the thesis text is true until you have verified it. It can happen that something is wrong or wrongly written that we would cross-check with numerical results. Numerics and Analytics should go hand in hand.
- All simulations use QuantumFurnace.jl infrastructure (HamHam, OFT, Lindbladian construction)
- System sizes: n=3,4,5 primary; n=6 if computationally feasible
- Heisenberg XXZ chain with Z and Z+ZZ disordering (established in prior work)
- Temperature range: β ∈ [5, 10, 20] typical; higher β more interesting but harder
- Scripts should go in `scripts/` directory, following existing naming conventions
