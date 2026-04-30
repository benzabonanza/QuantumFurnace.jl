# Implementing the Boltzmann rotation
<!-- Insert between: end of "Boltzmann filter and weak measurement" enumeration (line 33) and "Reverse pass" \medskip (line 35) of drafts/dissipative-step.md -->
<!-- Resolves: [REF to text] placeholder at line 1582 of 2_methods.tex -->
<!-- Depends on: Designing transition weights subsection (pp. 51--57), Proposition 7 (smooth Metropolis is Gevrey-1/2) -->

\medskip

**Implementing the Boltzmann rotation.**  The Boltzmann filter (step 7) requires a circuit that, controlled on the $r$-qubit energy register $|\bar\omega\rangle$, rotates $q_\gamma$ by the angle $\theta(\bar\omega) = 2\arcsin\sqrt{1 - \gamma(\bar\omega)}$.  This is a function-controlled rotation: the classical function $\gamma$ must be "compiled" into a quantum circuit acting on the joint $\Omega \otimes q_\gamma$ space.  The cost depends on the regularity of $\gamma$ and on the register size $r$.

For *large* $r$ (the fault-tolerant regime), the standard approach is to use QSP/QSVT \cite{gilyen2019quantum} to implement a polynomial approximation of $\sqrt{\gamma(\bar\omega)}$ as a function of the block-encoded normalised energy $\bar\omega / (\omega_0 N/2)$.  The required polynomial degree $d_\gamma$ depends on the regularity of $\gamma$:
\begin{itemize}
\item For the Gaussian $\gamma_G$ (analytic) or the smooth Metropolis $\gamma_M^{(s)}$ with $s > 0$ (Gevrey-$1/2$, Proposition 7): $d_\gamma = \mathcal{O}(\mathrm{polylog}(1/\varepsilon_\gamma))$, since Gevrey-$1/2$ functions admit super-algebraically converging polynomial approximations.
\item For the kinky Metropolis $\gamma_M$ ($s = 0$, only Lipschitz at $\omega = -\beta\sigma^2/2$): $d_\gamma = \mathcal{O}(\mathrm{poly}(1/\varepsilon_\gamma))$, since Lipschitz functions require polynomially many terms.  As noted in the Designing transition weights subsection, a comparator circuit can split the domain at the kink and apply QSP on each smooth piece separately \cite[footnote 33]{chen2023quantum}, but this adds circuit complexity --- and is one of the motivations for the smooth Metropolis family.
\end{itemize}

For *moderate* $r$ (say $r \leq 12$, which covers the parameter regime of our quadrature analysis in Table 5.1), one can bypass QSP entirely and use uniformly controlled $R_Y$ rotations \cite{mottonen2004transformation,shende2006synthesis} --- the same decomposition technique used for the state preparations $|f\rangle$, $|b_-\rangle$, $|b_+\rangle$ in the CoherentStep.  The $N = 2^r$ rotation angles $\theta(\bar\omega)$ are precomputed classically, and the resulting uniformly controlled rotation decomposes into $\mathcal{O}(2^r)$ CNOT + single-qubit gates.  At $r = 10$--$12$ this amounts to a few thousand gates, subdominant to the $\mathcal{O}(M)$ Trotter steps in the controlled Hamiltonian evolutions.  This is the approach we use in practice.

In either case, the Boltzmann rotation is not a bottleneck: for smooth $\gamma$ and large $r$, QSP gives polylogarithmic overhead; for moderate $r$, the brute-force decomposition is already cheap.  This is another instance where the regularity of $\gamma$ pays off --- a practical advantage of the smooth Metropolis family $\gamma_M^{(s>0)}$ introduced in the Designing transition weights subsection.

---
## Writing Notes
<!-- These notes are for the author, not for the thesis -->

- **Cross-references to verify in TeX**: Proposition 7 (smooth Metropolis Gevrey-1/2), Table 5.1 (quadrature parameter choices giving $r$ values), equation (5.34) for the smooth Metropolis definition, equation (5.31) for the kinky Metropolis.  The phrase "Designing transition weights subsection" should become a `\ref{}` to the actual subsection label.
- **Citations to verify**: `\cite{gilyen2019quantum}` = QSVT paper (Gilyen, Su, Low, Wiebe 2019), `\cite{chen2023quantum}` = CKBG23 (footnote 33 for the comparator circuit suggestion), `\cite{mottonen2004transformation}` = Mottonen et al. 2004 (uniformly controlled rotations), `\cite{shende2006synthesis}` = Shende, Bullock, Markov 2006 (same technique, independent derivation).
- **Review item addressed**: This insertion resolves review item 3.3 (P1: "Explain controlled-$\gamma(\bar\omega)$ rotation and its cost") and the `[REF to text]` placeholder at line 1582 of `2_methods.tex` (review item 4.1, P0).
- **Consistency with CoherentStep**: The state preparation discussion in `drafts/coherent-step.md` (paragraph starting "At these moderate sizes...") uses the same Mottonen et al. decomposition and the same $\mathcal{O}(2^r)$ gate count argument.  The parallel is deliberate and should be maintained.
- **[CHECK]**: The claim that $r \leq 12$ covers the practical regime should be cross-checked against the actual quadrature parameter choices in Table 5.1 and equation (5.60).  The dissipative step's $r$ is set by the OFT discretisation, which may differ from the coherent step's $r_\pm$.
- **What was NOT included** (by design): The detailed QSP polynomial degree analysis (which polynomial exactly, how to reduce $\sqrt{\gamma}$ to a QSP problem on the block-encoded energy) is deferred to Chapter 4 (Quantum Circuits), which is currently a placeholder.  The paragraph here is deliberately brief, following the existing draft's convention of "focus on structural features specific to our implementation" rather than reproducing CKBG23's circuit-level details.
