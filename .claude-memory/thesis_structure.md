---
name: Thesis Structure and Status
description: Bence's MSc thesis on quantum Gibbs sampling — structure, completed/incomplete sections, and original contributions
type: project
---

# Thesis: "Quantum Gibbs Sampling: Theory, Algorithms, and Numerical Analysis"
Author: Bence Temesi. LaTeX-generated, 92 pages as of 2026-03-30. Located at `supplementary-informations/thesis.pdf`.

## Structure

### Part I: Preliminaries (Ch 1-4) — mostly written
- **Ch 1**: Quantum Formalism and Notation (states, measurements) — written
- **Ch 2**: Quantum Dynamics (closed/open systems, Kraus, Choi, Lindbladian, microscopic derivation) — written
- **Ch 3**: Detailed Balanced Markov Dynamics — written, substantial (~40 pages)
  - 3.1 Classical Markov Semigroups (Metropolis-Hastings, mixing time, spectral bounds, log-Sobolev, continuous-time)
  - 3.2 Quantum Markov Semigroups (QDB, GNS vs KMS, Davies generator, quantum spectral bounds, QLSI)
- **Ch 4**: Quantum Circuits — **PLACEHOLDER** (all subsections just titles: QPE, Trotterization, LCU, QSP, weak-measurement)

### Part II: Gibbs Sampling on a Quantum Computer (Ch 5) — main chapter, ~42 pages, mostly written
- Approximate GNS Lindbladian (CKBG)
- KMS Lindbladian (CKG) — derivation of coherent term G
- Designing transition weights (Gaussian, Metropolis, Smooth Metropolis)
- **Proposition 5**: Optimality of Metropolis weights for CKG (original proof)
- **Proposition 7**: Smooth Metropolis is Gevrey-1/2 (original result)
- Improved Kossakowski matrices (Corollary 8: shift between CKG/CKBG)
- Time-domain representations (b+, b- kernels, regularization)
- Quadrature errors for Ldiss (Table 5.1 convergence hierarchy: Gaussian vs Smooth Metropolis vs Metropolis)
- Quadrature errors for B (outer/inner integrals, nested structure)
- Hamiltonian Trotterization (Proposition 9: Trotter for Ldiss, Proposition 10: Trotter for B)
- Generator splitting (jump-wise Lie-Trotter + coherent-dissipative)
- Remark 11: Palindromic product formulas preserve OFT adjoint identity
- **Algorithm** — **PLACEHOLDER** (just title)
- **Parent Hamiltonian** — **PLACEHOLDER** (just title)

### Part III: Review — **EMPTY PLACEHOLDER**
### Part IV: QuantumFurnace.jl — **EMPTY PLACEHOLDER**
### Appendix A — **EMPTY PLACEHOLDER** ("Quantum this quantum that, la la la")

## Status / TODOs visible in text
- Many margin annotations: "ref", "fig", "plot", "cites" throughout
- Abstract: "Your abstract text here..."
- Acknowledgment: empty
- Ch 1: TODO note about "correlation and entanglement and distances"
- Ch 4: All subsections empty
- Ch 5 Algorithm section: empty
- Ch 5 Parent Hamiltonian section: empty
- Parts III, IV: empty
- Appendix: empty

## Key references
- [CKBG23] = Chen et al. 2023 "Quantum Thermal State Preparation" (approximate GNS)
- [CKG23] = Chen et al. 2023 "An efficient and exact noncommutative quantum Gibbs sampler" (exact KMS)
- [DLL25] = Ding, Li, Lin 2025 "Efficient quantum Gibbs samplers with KMS detailed balance"
- [Chi+21] = Childs et al. 2021 "Theory of Trotter error with commutator scaling"

**Why:** This is the #1 reference for understanding what Bence is working on. The thesis IS the project.
**How to apply:** When implementing numerical analysis features in QuantumFurnace.jl, align with the theoretical framework and notation in the thesis. The smooth Metropolis family (Prop 5, 7), Gevrey analysis, and Trotter error bounds (Props 9, 10) are original contributions.
