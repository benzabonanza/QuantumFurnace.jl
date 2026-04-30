---
name: GQSP angle-finding: Berntson-Sünderhauf (final decision)
description: Final thesis/numerics decision — use GQSP (Motlagh-Wiebe framework) but compute the complementary polynomial Q via Berntson-Sünderhauf's FFT-based algorithm instead of Motlagh-Wiebe's optimization.
type: project
originSessionId: dee4ce01-d734-4d17-b307-e0de75b30e47
---
**Decision (2026-04-15):** For implementing e^{−iδB} via block-encoded GQSP in the thesis and upcoming numerics, use the **Motlagh-Wiebe GQSP framework** (arXiv:2401.10321 / arXiv:2308.01501) but replace their optimization-based complementary-polynomial step with **Berntson & Sünderhauf's FFT-based algorithm** (arXiv:2406.04246).

**Why:**
- GQSP (complex-polynomial convention) avoids the factor-of-2 overhead that a standard QSVT/LCU decomposition of cos/sin would incur — this is Berry et al.'s "doubling efficiency" result we already cite.
- Berntson-Sünderhauf is substantially easier to implement than Motlagh-Wiebe's nonlinear optimization: no solver, no initialization, no convergence tuning. Just FFTs + a discrete Hilbert transform + rigorous a priori error bounds.
- Endorsed by Skelton's "Hitchhiker's Guide to QSP pre-processing" (arXiv:2501.05977) — she recommends Berntson-Sünderhauf for complex polynomials (i.e. the GQSP case). Colleague = Skelton, so citation coherence is natural.
- Pipeline is end-to-end rigorous: Jacobi-Anger truncation error (Bessel tail) + Berntson-Sünderhauf error bound + Motlagh-Wiebe exact recursion for angles.

**How to apply:**
- Thesis citations for the Hamiltonian-simulation subroutine: Motlagh-Wiebe (GQSP framework) + Berntson-Sünderhauf (complementary polynomial) + Skelton (pre-processing recommendation). Keep Berry et al. for the doubling argument.
- Numerics pipeline (Julia):
  1. P = Jacobi-Anger expansion of e^{−iδx}; coefficients are Bessel J_k(δ) via `SpecialFunctions.jl`. Truncate at d = O(δ + log(1/ε)).
  2. Q via Berntson-Sünderhauf: FFT of |P|² on unit circle → log(1 − |P|²) → Hilbert transform (FFT-based) → exponentiate → inverse FFT. Watch conditioning when |P| ≈ 1 (rescale P by (1 − η) if needed).
  3. Angles: Motlagh-Wiebe's closed-form recursion (their Alg. 1) given both P and Q.
- Do **not** use QSPPACK / Lin Lin's optimization pipeline — wrong convention (QSVT, not GQSP) and pays the factor-of-2.
- **No Yamamoto-Yoshioka Prony cross-check** — dropped for time. Rely on Berntson-Sünderhauf's a priori error bounds + verifying that the GQSP circuit reproduces e^{−iδB} to target precision on test Hamiltonians.

**Key references (in supplementary-informations/):**
- Motlagh-Wiebe 2024 — Generalized Quantum Signal Processing (already present)
- Berry et al. — Doubling Efficiency of Hamiltonian Simulation via GQSP (already present)
- Berntson-Sünderhauf 2024 — Complementary polynomials in QSP (arXiv:2406.04246) — **need to add**
- Skelton 2025 — Hitchhiker's Guide to QSP pre-processing (arXiv:2501.05977) — **need to add**
