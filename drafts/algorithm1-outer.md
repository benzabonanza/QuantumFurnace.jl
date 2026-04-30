# Algorithm 1: CKG-GibbsSampler (Outer Structure)
<!-- Suggested placement: Chapter 5, Section 5.X "Algorithm", immediately after the Generator splitting subsection -->
<!-- Depends on: CKG Lindbladian (5.30), Generator splitting (5.82)--(5.93), mixing time (3.51)/(3.55) -->

---

**Algorithm.**
We now assemble the ingredients developed in the preceding subsections into a complete quantum algorithm for Gibbs state preparation.
The CKG Lindbladian~\eqref{eq:ckg-L} defines the target dynamics; the generator splitting analysis tells us how to decompose it into implementable pieces; what remains is to specify the outer control flow, the step-size budget, and the interplay between the two levels of splitting.
The resulting procedure, Algorithm~\ref{alg:main}, combines three independently developed techniques:
the \emph{weak-measurement simulation scheme} of \cite{chen2023quantum}, which implements the dissipative channels $e^{\delta\mathcal{L}_{a,\mathrm{diss}}}$ at first-order accuracy in $\delta$;
the \emph{coherent correction term} $B_a$ of \cite{chen2023efficient}, whose Hamiltonian evolution upgrades the approximate GNS detailed balance of the CKBG Lindbladian to exact KMS detailed balance at the continuous level;
and \emph{Generalised Quantum Signal Processing} (GQSP) from \cite{motlagh2024generalized}, which provides a near-optimal method to realise the unitary $e^{-\mathrm{i}\delta B_a}$ from a block encoding of $B_a$.
The end result is a quantum channel that, after $L$ iterations of the outer loop, produces a state $\varepsilon$-close to $\rho_\beta$ in trace distance.

We now walk through Algorithm~\ref{alg:main} line by line.

\paragraph{Total number of steps (line 1).}
The number of outer-loop iterations is set to
\begin{equation*}
    L = \left\lceil \frac{t_{\mathrm{mix}}(\mathcal{L})\,\log(2/\varepsilon)}{\delta} \right\rceil.
\end{equation*}
The logic is as follows.
From the mixing time definition~(3.51) and the spectral bound~(3.55), after evolving for a total Lindbladian time $t = t_{\mathrm{mix}}(\mathcal{L})\,\log(2/\varepsilon)$ we are guaranteed
$\|e^{t\mathcal{L}}(\rho_0) - \rho_\beta\|_1 \leq \varepsilon$
for any initial state $\rho_0$.
Each outer-loop iteration performs a sequential sweep over all $M_\mathcal{A} := |\mathcal{A}|$ jumps, and, as shown in~(5.82), the composition
$\Phi_\mathcal{A} := e^{\delta\mathcal{L}_{M_\mathcal{A}}} \circ \cdots \circ e^{\delta\mathcal{L}_1}$
approximates $e^{\delta\mathcal{L}}$ to $\mathcal{O}(\delta^2)$.
Crucially, $\delta$ here is the step size that appears in \emph{each} individual factor $e^{\delta\mathcal{L}_a}$, but the net effect of one full sweep is an advancement of the Lindbladian time by $\delta$, not by $M_\mathcal{A}\delta$.
This is the natural reading of the first-order Lie--Trotter product: the same parameter $\delta$ enters each factor, and the product approximates $e^{\delta\sum_a \mathcal{L}_a} = e^{\delta\mathcal{L}}$.
Hence $L$ sweeps advance the total simulated time by $L\delta \geq t_{\mathrm{mix}}(\mathcal{L})\,\log(2/\varepsilon)$, as required.

\paragraph{Initial state (line 2).}
The system register $S$ is initialised in an arbitrary state $\rho_0$.
This is justified by the primitivity of $\mathcal{L}$: since the Gibbs state $\rho_\beta$ is the unique fixed point and the semigroup is ergodic, every initial state converges to $\rho_\beta$ at the same asymptotic rate.
The mixing time bound~(3.55) already accounts for the worst-case initial condition through the factor $1/\sqrt{\lambda_{\min}(\rho_\beta)}$, which enters the $\log(2/\varepsilon)$ overhead.
In practice, choosing an initial state with reasonable overlap with $\rho_\beta$ -- such as the maximally mixed state $\mathbb{1}/d$ -- can reduce the transient, but the asymptotic scaling is not affected.

\paragraph{Outer and inner loops (lines 3--7).}
The outer loop runs $L$ iterations.
Within each iteration, the inner loop sweeps sequentially over the $M_\mathcal{A}$ jumps $a = 1, \ldots, M_\mathcal{A}$.
For each jump $a$, two subroutines are called in order:

\begin{enumerate}
    \item \textsc{CoherentStep} (line 5): implements the unitary $e^{-\mathrm{i}\delta B_a}$, where $B_a$ is the per-jump coherent Hamiltonian from the decomposition $\mathcal{L}_a = \mathcal{L}_{a,\mathrm{coh}} + \mathcal{L}_{a,\mathrm{diss}}$ with $\mathcal{L}_{a,\mathrm{coh}}(\cdot) = -\mathrm{i}[B_a, \,\cdot\,]$.
    \item \textsc{DissipativeStep} (line 6): implements the channel $e^{\delta\mathcal{L}_{a,\mathrm{diss}}}(\rho) + \mathcal{O}(\delta^2)$ via the weak-measurement scheme of \cite{chen2023quantum}.
\end{enumerate}

This ordering -- coherent first, then dissipative -- corresponds to the first-order coherent--dissipative splitting from~(5.90):
\begin{equation*}
    e^{\delta\mathcal{L}_{a,\mathrm{coh}}}\, e^{\delta\mathcal{L}_{a,\mathrm{diss}}}
    = \exp\!\left(\delta\mathcal{L}_a + \frac{\delta^2}{2}[\mathcal{L}_{a,\mathrm{coh}}, \mathcal{L}_{a,\mathrm{diss}}] + \mathcal{O}(\delta^3)\right).
\end{equation*}
At first order, the two possible orderings (coherent--dissipative or dissipative--coherent) produce the same leading-order generator $\delta\mathcal{L}_a$; they differ only in the sign of the $\mathcal{O}(\delta^2)$ commutator correction.
Since the weak-measurement subroutine for $\mathcal{L}_{a,\mathrm{diss}}$ is itself only first-order accurate in $\delta$, upgrading to a Strang splitting~(5.93) would not improve the end-to-end error.
A first-order splitting is therefore the simplest, yet sufficient, choice.

The inner loop over $a$ realises the jump-wise Lie--Trotter splitting~(5.82):
\begin{equation*}
    \Phi_\mathcal{A} := e^{\delta\mathcal{L}_{M_\mathcal{A}}} \circ \cdots \circ e^{\delta\mathcal{L}_2} \circ e^{\delta\mathcal{L}_1},
\end{equation*}
with per-step approximation error $\|\Phi_\mathcal{A} - e^{\delta\mathcal{L}}\|_{1\to 1} = \mathcal{O}\!\left(\delta^2 \log^2(\beta\|H\|)\right)$ from~(5.86).
The $\log^2(\beta\|H\|)$ factor originates from the norm of the individual generators $\mathcal{L}_a$, which is dominated by the coherent term $B_a$ and its Metropolis-like regularisation; for purely Gaussian transition weights this factor reduces to $\mathcal{O}(\delta^2)$.

\paragraph{Jump-wise decomposition: sequential sweep vs.\ random sampling.}
It is worth clarifying the two strategies for the inner loop that are mentioned in Algorithm~\ref{alg:main}.
As described in the Generator splitting subsection, one can either sweep sequentially over all jumps (lines 4--7 as written), or sample a single jump $a$ uniformly at random from $\mathcal{A}$ and rescale the step size by $M_\mathcal{A}$.
Asymptotically both strategies give the same complexity, but they break the KMS detailed balance condition in qualitatively different ways.

In the sequential case, the Lie--Trotter product~(5.82) breaks time-reversal symmetry: the KMS adjoint reverses the ordering of the factors, so $\Phi_\mathcal{A}^* \neq \Phi_\mathcal{A}$ as shown in~(5.87).
There is some structure to this error however, that is worth exploring.
By Baker--Campbell--Hausdorff, the effective generators of the forward and reverse orderings differ only in the sign of the leading commutator term~(5.88), which is KMS anti-self-adjoint:
$[\mathcal{L}_a, \mathcal{L}_b]^* = [\mathcal{L}_b^*, \mathcal{L}_a^*] = -[\mathcal{L}_a, \mathcal{L}_b]$.
Hence the DB violation is a purely antisymmetric perturbation, and by standard perturbation theory it shifts the eigenvalues of $\mathcal{L}$ only along the imaginary axis~(5.89).
The spectral gap -- and with it the mixing time -- is therefore preserved up to $\mathcal{O}(\delta^3)$.
Furthermore, the Gibbs state remains the exact fixed point of $\Phi_\mathcal{A}$, since each individual channel $e^{\delta\mathcal{L}_a}$ fixes $\rho_\beta$.

In the random-sampling case, each step applies $e^{\delta'\mathcal{L}_a}$ for a uniformly random $a$ (with $\delta' = M_\mathcal{A}\delta$), so the expected channel is $\frac{1}{M_\mathcal{A}}\sum_a e^{\delta'\mathcal{L}_a}$, which is a convex combination of KMS detailed balanced channels and thus preserves KMS DB on average.
The trade-off is that the leading-order errors now do perturb the spectral gap, rather than being purely imaginary.
Since these errors are of the same order as those introduced by the weak-measurement scheme, the choice between the two strategies makes no significant difference to the asymptotic scaling.
We adopt the sequential sweep throughout for concreteness.

\paragraph{The role of the coherent step.}
Without the \textsc{CoherentStep} -- that is, setting $B_a = 0$ for all $a$ -- the algorithm reduces to the CKBG Gibbs sampler of \cite{chen2023quantum}, which simulates the Lindbladian $\mathcal{L}_\mathrm{CKBG}$~(5.6) with only approximate GNS detailed balance.
The coherent correction term $B_a$, derived in~(5.28)--(5.29) as the unique Hermitian operator that upgrades the Kossakowski-matrix symmetry from GNS to KMS, is the essential new ingredient of the CKG construction \cite{chen2023efficient}.
At the continuous level, $\mathcal{L}_\mathrm{CKG}$~(5.30) satisfies exact KMS detailed balance by construction.
The two levels of splitting -- jump-wise and coherent--dissipative -- re-introduce controlled violations: the jump-wise splitting yields an antisymmetric DB error that preserves the fixed point~(5.87)--(5.89), while the coherent--dissipative splitting shifts the fixed point by $\mathcal{O}(t_\mathrm{mix}\,\delta)$~(5.92).
Both are of the same order as the baseline weak-measurement error, so they do not degrade the overall accuracy.

\medskip
The detailed circuit-level implementations of \textsc{CoherentStep} (Algorithm~\ref{alg:coh}) and \textsc{DissipativeStep} (Algorithm~\ref{alg:diss}) are described in the paragraphs that follow.

---
## Writing Notes
<!-- These notes are for the author, not for the thesis -->
- The equation reference style (3.51), (3.55), (5.30), etc. uses the **thesis PDF numbering** as of the current version (92 pages). If equation numbers shift after edits to earlier sections, all references here will need updating.
- **Verified TeX labels**: `eq:ckg-L` = CKG Lindbladian (5.30); `eq:quantum-spectral-bound` = spectral bound (3.55); `eq:gen-split-approx-error` = commutator error bound (5.83). The Lie-Trotter product formula (5.82) does **not** have a `\label{}` in the current TeX source -- one should be added (e.g. `eq:lie-trotter-jump`) when integrating this prose.
- The quantum mixing time definition (3.51) does not have a visible label in the TeX excerpt; [CHECK the label in `1_preliminaries.tex` around line 724].
- [CHECK] The statement "the net effect of one full sweep is an advancement of the Lindbladian time by $\delta$, not by $M_\mathcal{A}\delta$" -- this is the correct interpretation consistent with the pseudocode. Each factor $e^{\delta\mathcal{L}_a}$ uses the same $\delta$, and the product approximates $e^{\delta\mathcal{L}}$. Confirm this is the intended convention by checking how $\delta$ is chosen in the error budget section when it is written.
- [Plots needed: none for this subsection]
- [Connection points]: This subsection follows immediately after the Generator splitting subsection (which ends with the bold header "Algorithm." and no prose). The forward references to Algorithms 2 and 3 connect to the CoherentStep and DissipativeStep prose drafts (not yet written).
- The phrase "There is some structure to this error however, that is worth exploring" is a direct stylistic callback to p. 82 of the thesis.
- The phrase "the simplest, yet sufficient, choice" echoes p. 82 ("The simplest (yet sufficient) way to split the generator sequentially...").
- Citation keys used: `chen2023quantum` (CKBG23), `chen2023efficient` (CKG23), `motlagh2024generalized` (MW24).
- The `\eqref{}` reference on line 9 (`\eqref{eq:ckg-L}`) is the only one that uses a TeX label directly; the rest use hardcoded equation numbers like (5.82). When converting to TeX, all hardcoded numbers should be replaced with `\eqref{}` references. The corresponding labels that need to be added or verified: (5.82) needs a label, (5.86) = `eq:gen-split-explicit-bound` [CHECK], (5.87)--(5.89) need labels, (5.90) needs a label, (5.92) needs a label, (5.93) needs a label, (5.6) = `eq:L-ckbg`, (5.28)--(5.29) = check around `eq:G-bohr-domain`.
