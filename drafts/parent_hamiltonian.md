# Parent Hamiltonian
<!-- Suggested placement: Chapter 2 (Methods), final subsection (after Algorithm) -->
<!-- Depends on: Preliminaries discriminant eq:discriminant; quantum spectral bound eq:quantum-spectral-bound; KMS Lindbladian subsection -->

\subsection{Parent Hamiltonian.}
Beyond its role as a design principle for the algorithm, KMS detailed balance carries a powerful analytical payoff: it reduces the spectral analysis of the Lindbladian to a Hermitian eigenvalue problem. Via the discriminant transformation introduced in \eqref{eq:discriminant}, the non-Hermitian generator $\mathcal{L}$ is mapped to a positive semidefinite superoperator --- the \textit{parent Hamiltonian} --- with the same spectrum. We now specialise this construction to the CKG Lindbladian and discuss its implications.

Recall that for a Lindbladian $\mathcal{L}$ satisfying $\rho_\beta$-KMS detailed balance, the discriminant
\begin{equation}
    \mathcal{D}(\rho_\beta, \mathcal{L}) = \rho_\beta^{-1/4}\,\mathcal{L}\!\left(\rho_\beta^{1/4}\,(\cdot)\,\rho_\beta^{1/4}\right)\rho_\beta^{-1/4}
\end{equation}
is a Hermitian superoperator that shares the spectrum of $\mathcal{L}$, which lies on the non-positive real line. We define the \textit{parent Hamiltonian} as
\begin{equation}\label{eq:parent-hamiltonian}
    \mathcal{H}_\beta \;:=\; -\mathcal{D}(\rho_\beta, \mathcal{L}) \;\geq\; 0,
\end{equation}
the sign-flipped discriminant, so that $\mathcal{H}_\beta$ is positive semidefinite. Vectorizing $\mathcal{H}_\beta$ onto the doubled Hilbert space $\mathcal{H}\otimes\mathcal{H}$, its ground state is the \textit{purified Gibbs state} (or canonical purification / thermal field double)
\begin{equation}
    |\sqrt{\rho_\beta}\rangle\!\rangle \;:=\; \frac{1}{\sqrt{\tr[\ee^{-\beta H}]}}\sum_i \ee^{-\beta E_i/2}\,|\psi_i\rangle\otimes|\psi_i^*\rangle,
\end{equation}
where $\{|\psi_i\rangle\}$ are the energy eigenstates and $|\psi_i^*\rangle$ denotes the entrywise complex conjugate. Since the discriminant preserves the spectrum, the spectral gap of the parent Hamiltonian equals that of the Lindbladian, $\lambda_\text{gap}(\mathcal{H}_\beta) = \lambda(\mathcal{L})$. Thus, the mixing time question --- how fast does $\ee^{t\mathcal{L}}(\rho_0)$ converge to $\rho_\beta$? --- reduces to bounding the spectral gap of a Hermitian positive semidefinite operator. The quantum spectral bound \eqref{eq:quantum-spectral-bound} directly applies:
\begin{equation}\label{eq:parent-ham-gap-mixing}
    t_\text{mix}(\varepsilon) \;\leq\; \frac{1}{\lambda_\text{gap}(\mathcal{H}_\beta)} \ln \left( \frac{1}{\varepsilon\sqrt{\lambda_{\min}(\rho_\beta)}} \right).
\end{equation}
By reformulating the problem as a Hermitian spectral gap, we gain access to the full toolkit of Hamiltonian spectral theory: variational (min-max) characterisations \cite[Ch.~III]{bhatia1997matrix}, perturbation theory for gapped self-adjoint operators \cite{kato1995perturbation}, and the detectability lemma \cite{aharonov2009detectability,anshu2016simple} --- none of which are directly available for the original non-Hermitian Lindbladian.

\smallskip

The key structural property, proven by Chen, Kastoryano and Gily\'{e}n \cite[Proposition~I.1]{chen2023efficient}, is that this parent Hamiltonian is \textit{frustration-free}: the CKG Lindbladian decomposes as $\mathcal{L} = \sum_{a\in\mathcal{A}} \mathcal{L}^a_\beta$, where each per-jump term $\mathcal{L}^a_\beta$ individually satisfies $\mathcal{L}^a_\beta(\rho_\beta) = 0$. Under the discriminant transformation this carries over to the parent Hamiltonian,
\begin{equation}
    \mathcal{H}_\beta = \sum_{a\in\mathcal{A}} \mathcal{H}_\beta^a, \qquad \mathcal{H}_\beta^a \geq 0, \qquad \mathcal{H}_\beta^a\,|\sqrt{\rho_\beta}\rangle\!\rangle = 0 \quad\text{for each}\;\; a\in\mathcal{A}.
\end{equation}
That is, every local term $\mathcal{H}_\beta^a$ is positive semidefinite and has the purified Gibbs state in its kernel. For lattice Hamiltonians with local jumps $A^a$, the parent Hamiltonian inherits the quasi-locality of the Lindbladian: each term $\mathcal{H}_\beta^a$ acts on a patch of radius
\begin{equation}\label{eq:parent-ham-locality-radius}
    r \;=\; \mathcal{O}\!\left(v_\text{LR}\,\beta\,\sqrt{\log(1/\varepsilon)}\,\right)
\end{equation}
around the support of $A^a$, up to tails of operator norm at most $\varepsilon$, where $v_\text{LR}$ is the Lieb--Robinson velocity of $H$ \cite[Fig.~2]{chen2023efficient}. The $\sqrt{\log(1/\varepsilon)}$ factor --- rather than a plain $\log(1/\varepsilon)$ --- originates from the Gaussian profile of the operator Fourier transform: the time-domain truncation of the Gaussian-weighted kernel $f(t)$ of width $\sigma_t = \Theta(\beta)$ required to achieve tail error $\varepsilon$ is $T = \Theta(\beta\sqrt{\log(1/\varepsilon)})$ \cite[App.~C]{chen2023quantum}, and Lieb--Robinson then converts this time cutoff into the spatial radius \eqref{eq:parent-ham-locality-radius}. Frustration-freeness has three concrete analytical payoffs:

\begin{enumerate}
    \item \textit{Spectral gap stability.} Rouz\'{e}, Stilck~Fran\c{c}a and Alhambra \cite{rouze2025efficient} proved that for 1D lattice Hamiltonians, the CKG parent Hamiltonian satisfies local topological quantum order (LTQO) at $\beta = 0$ where the gap is explicitly computable. Using the Michalakis--Zwolak stability theorem for frustration-free gapped Hamiltonians \cite{michalakis2013stability}, they showed that the gap persists at all sufficiently small inverse temperatures $\beta$, establishing the first unconditional rapid mixing result for a KMS detailed balanced Lindbladian with noncommuting $H$. This was recently strengthened by Bergamaschi and Chen \cite{bergamaschi2025fast}, who showed that 1D spin chains mix rapidly at \textit{all} finite temperatures ($\beta < \infty$), not just in the high-temperature regime. Without frustration-freeness, the stability theorem simply does not apply.

    \item \textit{Detectability lemma and gap amplification.} Since the ground-state projector factorises into local projectors (one per kernel of $\mathcal{H}_\beta^a$), the detectability lemma of Anshu, Arad and Vidick \cite{anshu2016simple} provides a way to amplify the spectral gap via coarse-graining: grouping the local terms into $\mathcal{O}(1)$ layers of commuting projectors yields a product of projections whose spectral gap lower-bounds that of $\mathcal{H}_\beta$ up to constant factors. This is a key ingredient in the Rouz\'{e}--Stilck~Fran\c{c}a--Alhambra proof. Very recently, Fang, Lu, Tong and Zhao \cite{fang2026detectability} leveraged this structure \textit{algorithmically}, using the detectability lemma applied to the parent Hamiltonian to prepare the Gibbs state without ever simulating the Lindbladian. Their construction yields an $\mathcal{O}(|\mathcal{A}|)$-factor gate saving over simulation-based methods in general, and --- when the parent Hamiltonian is strictly local (e.g.\ for commuting $H$) --- a quadratic improvement in the spectral-gap dependence via quantum singular value transformation.

    \item \textit{Gap-to-MLSI upgrade and rapid mixing.} The bound \eqref{eq:parent-ham-gap-mixing} incurs a $\log(1/\lambda_{\min}(\rho_\beta)) = \mathcal{O}(\beta\|H\|)$ factor, so the spectral gap alone yields a mixing time polynomial in system size. Upgrading this to the genuinely \textit{rapid mixing} regime $t_\text{mix}(\varepsilon) = \tilde{\mathcal{O}}(\log N \,\log(1/\varepsilon))$ requires a modified logarithmic Sobolev inequality (MLSI) with a system-size-independent constant, and the frustration-free decomposition of $\mathcal{H}_\beta$ is the key ingredient enabling this conversion. Kochanowski, Alhambra, Capel and Rouz\'{e} \cite{kochanowski2024rapid} recently proved such a gap-to-MLSI equivalence for 1D Davies dynamics of commuting Hamiltonians at all positive temperatures, strictly strengthening earlier 1D rapid-mixing results; the analogous statement in the noncommuting CKG setting remains open.
\end{enumerate}

\smallskip

Beyond proving mixing time bounds, the parent Hamiltonian structure yields numerical advantages that fall into a three-level hierarchy, each level strictly extending the previous one.

At the base level, Lindblad dynamics can be simulated iteratively without any special structure: Minganti and Huybrechts \cite{minganti2022arnoldi} apply Arnoldi iteration directly to the non-Hermitian Liouvillian $\mathcal{L}$, maintaining a full upper-Hessenberg matrix and re-orthogonalising against all accumulated basis vectors at every step. Non-normality of $\mathcal{L}$ further complicates the relationship between its spectral gap and the iteration count.

\textit{Hermiticity}, obtained by the discriminant transformation $\mathcal{L} \mapsto \mathcal{H}_\beta$, collapses Arnoldi to the Lanczos iteration. The recurrence becomes three-term, per-step storage and re-orthogonalisation reduce to $\mathcal{O}(N)$, and convergence of the Ritz values is now controlled directly by the spectral gap $\lambda_\text{gap}(\mathcal{H}_\beta)$ and the overlap of the starting vector with the ground space.

\textit{Frustration-freeness} adds a structural layer on top, via the \textit{detectability lemma} \cite{aharonov2009detectability,anshu2016simple}. Because the purified Gibbs state lies in the common kernel of the positive semidefinite local terms $\{\mathcal{H}_\beta^a\}$, the lemma turns this collection into an explicit operator --- the product of the corresponding local kernel projectors --- whose repeated action drives any state towards $|\sqrt{\rho_\beta}\rangle\!\rangle$ at a rate set by $\lambda_\text{gap}(\mathcal{H}_\beta)$, with no prior estimate of the ground-state energy required. On the quantum side, Fang, Lu, Tong and Zhao \cite{fang2026detectability} exploit precisely this operator to prepare the Gibbs state on a quantum register via singular value transformation, bypassing simulation of the Lindbladian entirely. On the classical side, the same lemma is the analytical backbone of the rigorous polynomial-time ground-state algorithms for 1D gapped local Hamiltonians of Landau, Vazirani and Vidick \cite{landau2015polynomial} and Arad, Landau, Vazirani and Vidick \cite{arad2017rigorous}.

---
## Writing Notes
- [New bibliography entries needed]:
  - `michalakis2013stability` -- Michalakis & Zwolak, "Stability of frustration-free Hamiltonians", Comm. Math. Phys. 322, 277--302 (2013), arXiv:1109.1588.
  - `bergamaschi2025fast` -- Bergamaschi & Chen, "Fast Mixing of Quantum Spin Chains at All Temperatures", arXiv:2510.08533 (2025).
  - `anshu2016simple` -- Anshu, Arad, Vidick, "Simple proof of the detectability lemma and spectral gap amplification", Phys. Rev. B 93, 205142 (2016), arXiv:1602.01166.
  - `kochanowski2024rapid` -- Kochanowski, Alhambra, Capel, Rouz\'{e}, "Rapid thermalization of dissipative many-body dynamics of commuting Hamiltonians", arXiv:2404.16780 (2024).
  - `minganti2022arnoldi` -- Minganti & Huybrechts, "Arnoldi-Lindblad time evolution", Quantum 6, 649 (2022), arXiv:2104.00442.
  - `landau2015polynomial` -- Landau, Vazirani, Vidick, "A polynomial time algorithm for the ground state of one-dimensional gapped local Hamiltonians", Nature Physics 11, 566--569 (2015), arXiv:1307.5143.
  - `arad2017rigorous` -- Arad, Landau, Vazirani, Vidick, "Rigorous RG algorithms and area laws for low energy eigenstates in 1D", Comm. Math. Phys. 356, 65--105 (2017), arXiv:1602.08828.
  - `bhatia1997matrix` -- R. Bhatia, *Matrix Analysis*, Graduate Texts in Mathematics 169, Springer (1997). Cited for the min-max (Courant--Fischer) principle, Chapter III "Variational Principles for Eigenvalues" (Corollary III.1.2). Check whether already in thesis bibliography for other matrix-analysis results.
  - `kato1995perturbation` -- T. Kato, *Perturbation Theory for Linear Operators*, Classics in Mathematics, Springer (1995; reprint of the 1980 2nd ed.). Cited for analytic perturbation of isolated eigenvalues of self-adjoint operators under a spectral gap assumption.
  - `aharonov2009detectability` -- Aharonov, Arad, Landau, Vazirani, "The detectability lemma and quantum gap amplification", Proc. STOC 2009, pp. 417--426, arXiv:0811.3412. Original detectability lemma.
  - `fang2026detectability` -- D. Fang, J. Lu, Y. Tong, C. Zhao, "Quantum Gibbs sampling through the detectability lemma", arXiv:2604.07214 (2026).
- [Existing bib keys to verify]: `rouze2025efficient`, `chen2023efficient` -- should match what's already in the bibliography.
- [Equation cross-references]: `\eqref{eq:discriminant}` and `\eqref{eq:quantum-spectral-bound}` refer back to Preliminaries. Verify numbering after integration.
- [Length]: ~1.5 pages of LaTeX, within target.
