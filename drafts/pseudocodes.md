Algorithm 1:

```
\begin{algorithm}[t]

\caption{\textsc{CKG-GibbsSampler}}

\label{alg:main}

\begin{algorithmic}[1]

\Require

System Hamiltonian $H = \sum_{k=1}^{K} H_k$ on $n$ qubits;

jump set $\{A_a\}_{a\in\mathcal{A}}$ with $\bigl\lVert\sum_a A_a^\dagger A_a\bigr\rVert \le 1$;

inverse temperature $\beta$;

Gaussian filter width $\sigma$ for filter function $f$ \eqref{eq:gaussian-filter-t};

transition weight $\gamma$ (e.g.\ $\gamma_M^{(s)}$ from~\eqref{eq:smooth-metro});

functions $b_\pm$ for the coherent term \eqref{eq:b_plus-s-eta} and \eqref{eq:b_minus};

step size $\delta$;

target trace-distance accuracy $\varepsilon$

\Ensure System register $S$ in state

$\rho \approx \rho_\beta = e^{-\beta H}/\mathrm{tr}(e^{-\beta H})$

\Statex

\State $L \gets \bigl\lceil t_{\mathrm{mix}}(\mathcal{L})\,

\log(2/\varepsilon)\,/\,\delta \bigr\rceil$

\Comment{total Lindbladian steps}

\State Initialise $S$ in an arbitrary state $\rho_0$

\For{$\ell = 1, \ldots, L$}

\For{$a = 1, \ldots, M_{\mathcal{A}}$}

\Comment{sequential sweep over jumps; or sample $a$ uniformly}

\State \Call{CoherentStep}{$S, a, \delta, b_\pm$}

\Comment{$e^{-\mathrm{i}\delta B_a}$, Alg.~\ref{alg:coh}}

\State \Call{DissipativeStep}{$S, a, \delta, f, \gamma$}

\Comment{$e^{\delta\mathcal{L}_{a,\mathrm{diss}}}$, Alg.~\ref{alg:diss}}

\EndFor

\EndFor

\State \Return $S$

\end{algorithmic}

\end{algorithm}
```


Algorithm 1A, (CoherentStep):
```
\begin{algorithm}[t]

\caption{\textsc{CoherentStep}: Hamiltonian simulation of

$e^{-\mathrm{i}\delta B_a}$ via block encoding and QSP}

\label{alg:coh}

\begin{algorithmic}[1]

\Require System register $S$, jump index $a$, step size $\delta$

\Ensure $S$ evolved by $e^{-\mathrm{i}\delta B_a}$

up to QSP and Trotter precision of $\mathcal{O}(\delta^2)$

\Statex

\State \textbf{Ancillas:}

outer time register $T_-$ ($r_-$ qubits),

inner time register $T_+$ ($r_+$ qubits),

QSP signal qubit $q_{\mathrm{QSP}}$

\Statex

\Statex\textit{--- Block encoding $U_{B_a}$ of

$B_a/\alpha$~\eqref{eq:B-discr}, with

$\alpha := \lVert b_-\rVert_1\,\lVert b_+\rVert_1$ ---}

\Statex\textit{Forward pass:}

\State \textsc{StatePrep}$(T_-)$:

$|0\rangle \to |b_-\rangle :=

\displaystyle\sum_{\bar{t}}

\sqrt{\frac{t_0\,|b_-(\bar{t})|}{\alpha_-}}\;

\mathrm{sgn}\bigl(b_-(\bar{t})\bigr)\,|\bar{t}\rangle$

\Comment{outer kernel~\eqref{eq:b_minus}; $\alpha_- = \lVert b_-\rVert_1$}

\State \textsc{Ctrl-Trot}$(T_-, S)$: apply

$S_2^{(M_-)}(-\bar{t}/\sigma)$ on $S$, controlled on $T_-$

\Comment{outer: $\approx e^{-\mathrm{i}H\bar{t}/\sigma}$}

\State \textsc{StatePrep}$(T_+)$:

$|0\rangle \to |b_+\rangle :=

\displaystyle\sum_{\bar{t}'}

\sqrt{\frac{t_0'\,|b_+(\bar{t}')|}{\alpha_+}}\;

\mathrm{sgn}\bigl(b_+(\bar{t}')\bigr)\,|\bar{t}'\rangle$

\Comment{inner kernel~\eqref{eq:b_plus-s-eta}; $\alpha_+ = \lVert b_+\rVert_1$}

\Statex

\Statex \textit{Inner Heisenberg evolution}

$$\tilde{A}_a^\dagger(\beta\bar{t}')\,\tilde{A}_a(-\beta\bar{t}')

= S_2^{(M_+)}(\beta\bar{t}')\,A_a^\dagger\,

S_2^{(M_+)}(-2\beta\bar{t}')\,A_a\,

S_2^{(M_+)}(\beta\bar{t}'),$$

\textit{controlled on $T_+$}:

\State \textsc{Ctrl-Trott}$(T_+, S)$: apply

$S_2^{(M_+)}(+\beta\bar{t}')$ on $S$

\Comment{$\approx e^{+\mathrm{i}\beta H\bar{t}'}$}

\State Apply $A_a^\dagger$ on $S$

\State \textsc{Ctrl-Trott}$(T_+, S)$: apply

$S_2^{(M_+)}(-2\beta\bar{t}')$ on $S$

\Comment{$\approx e^{-2\mathrm{i}\beta H\bar{t}'}$}

\State Apply $A_a$ on $S$

\State \textsc{Ctrl-Trott}$(T_+, S)$: apply

$S_2^{(M_+)}(+\beta\bar{t}')$ on $S$

\Comment{$\approx e^{+\mathrm{i}\beta H\bar{t}'}$}

\Statex

\State \textsc{StatePrep}$^\dagger(T_+)$:

$|b_+\rangle \to |0\rangle$

\Comment{reflect / uncompute inner register}

\State \textsc{Ctrl-Trott}$(T_-, S)$: apply

$S_2^{(M_-)}(+\bar{t}/\sigma)$ on $S$

\Comment{undo outer evolution}

\State \textsc{StatePrep}$^\dagger(T_-)$:

$|b_-\rangle \to |0\rangle$

\Comment{reflect / uncompute outer register}

\Statex \Comment{Block encoding:

$\langle 0|_{T_-}\!\langle 0|_{T_+}\;

U_{B_a}\;

|0\rangle_{T_-}\!|0\rangle_{T_+}

\;=\; B_a\,/\,\alpha$

acting on $S$}

\Statex

\Statex \textit{--- QSP: implement $e^{-\mathrm{i}\delta B_a}$

from the block encoding ---}

\State Apply QSP~\cite{motlagh2024generalized} with precomputed rotation angles

$\{(\theta_k,\phi_k)\}_{k=0}^{d}$ for the Laurent polynomial $P(e^{i\theta}) \approx e^{-\mathrm{i}\delta\alpha\,\sin\theta}$

\Statex \Comment{Degree

$d = \mathcal{O}\!\bigl(\delta\,\alpha

+ \log(1/\varepsilon_{\mathrm{QSP}})/\log\log(1/\varepsilon_{\mathrm{QSP}})\bigr)$

calls to $U_{B_a}$}

\State Measure and reset ancillas $q_{\mathrm{QSP}}$, $T_\pm$

\end{algorithmic}

\end{algorithm}
```


Algorithm 1B, DissipativeStep:
```
\begin{algorithm}[t]

\caption{\textsc{DissipativeStep}: weak-measurement simulation of

$e^{\delta\mathcal{L}_{a,\mathrm{diss}}}$}

\label{alg:diss}

\begin{algorithmic}[1]

\Require System register $S$, jump index $a$, step size $\delta$, filter function $f$, transition weight $\gamma$

\Ensure $S$ evolved by

$e^{\delta\mathcal{L}_{a,\mathrm{diss}}}(\rho) + \mathcal{O}(\delta^2)$

\Statex

\State \textbf{Ancillas:}

frequency register $\Omega$ ($r$ qubits),

Boltzmann qubit $q_\gamma$,

weak-measurement qubit $q_\delta$

\Comment{all initialised to $|0\rangle$}

\Statex

\Statex \hfill \textit{--- Forward: Trotterized OFT~\eqref{eq:trotter-oft} ---}

\State \textsc{StatePrep}$(\Omega)$:

$|0\rangle^{\otimes r} \;\longrightarrow\;

|f\rangle := \sum_{\bar{t}\in S_{t_0}^{[N]}}

\bar{f}(\bar{t})\,|\bar{t}\rangle$

\Comment{Gaussian filter amplitudes~\eqref{eq:gaussian-filter-t}}

\State \textsc{Ctrl-Trott}$^{-}(\Omega, S)$:

$|\bar{t}\rangle\,|\psi\rangle \;\longrightarrow\;

|\bar{t}\rangle\; S_p^{(M)}(\bar{t})^\dagger \,|\psi\rangle$

\Comment{$\approx e^{-\mathrm{i}H\bar{t}}$}

\State Apply $A_a$ on $S$

\Comment{jump operator (e.g.\ single-site Pauli)}

\State \textsc{Ctrl-Trott}$^{+}(\Omega, S)$:

$|\bar{t}\rangle\,|\psi'\rangle \;\longrightarrow\;

|\bar{t}\rangle\; S_p^{(M)}(\bar{t})\,|\psi'\rangle$

\Comment{$\approx e^{+\mathrm{i}H\bar{t}}$}

\State \textsc{QFT}$(\Omega)$

\Comment{time $\to$ energy: register now encodes $|\bar{\omega}\rangle$}

\Statex \Comment{Joint state:

$\displaystyle\sum_{\bar\omega}

|\bar\omega\rangle \otimes

\tilde{A}_a(\bar\omega)\,|\psi\rangle$}

\Statex

\Statex \textit{--- Boltzmann filter ---}

\State Controlled on $|\bar\omega\rangle$, rotate $q_\gamma$:\;

$R_Y\!\Bigl(2\arcsin\!\sqrt{1-\gamma(\bar\omega)}\Bigr)$

\Comment{see [REF to text]}

\Statex \Comment{$|0\rangle_{q_\gamma} \to

\sqrt{\gamma(\bar\omega)}\,|0\rangle

+ \sqrt{1 - \gamma(\bar\omega)}\,|1\rangle$}

\Statex

\Statex \textit{--- Weak measurement ---}

\State Controlled on $q_\gamma = |0\rangle$, rotate $q_\delta$:\;

$R_Y\!\bigl(2\arcsin\!\sqrt{\delta}\bigr)$

\Comment{step-size acceptance}

\Statex \Comment{In the $q_\gamma\!=\!|0\rangle$ branch:

$|0\rangle_{q_\delta} \to

\sqrt{1-\delta}\,|0\rangle + \sqrt{\delta}\,|1\rangle$;

\quad net jump probability $= \gamma(\bar\omega)\,\delta$}

\Statex

\Statex \textit{--- Reverse: uncompute $U$ (controlled on $q_\delta = |0\rangle$) ---}

\Statex \Comment{No-jump branch ($q_\delta\!=\!|0\rangle$): full uncomputation;

jump branch ($q_\delta\!=\!|1\rangle$): $\Omega, q_\gamma, S$ remain entangled}

\State Undo Boltzmann rotation on $q_\gamma$

\Comment{undo step~7}

\State \textsc{QFT}$^\dagger(\Omega)$

\Comment{undo step~6}

\State \textsc{Ctrl-Trott}$^{-}(\Omega, S)$

\Comment{undo step~5}

\State Apply $A_a^\dagger$ on $S$

\Comment{undo step~4; for Pauli jumps $A_a^\dagger = A_a$}

\State \textsc{Ctrl-Trott}$^{+}(\Omega, S)$

\Comment{undo step~3}

\State \textsc{StatePrep}$^\dagger(\Omega)$:

$|f\rangle \to |0\rangle^{\otimes r}$

\Comment{undo step~2}

\Statex

\Statex \textit{--- Discard ancillas ---}

\State Measure and reset $q_\delta$, $q_\gamma$, and $\Omega$

\Comment{discard outcomes}

\Statex \Comment{Net channel on $S$:

$\;\rho \;\mapsto\; e^{\delta\mathcal{L}_{a,\mathrm{diss}}}(\rho)

+ \mathcal{O}(\delta^2)$}

\end{algorithmic}

\end{algorithm}
```


# Tasks for Claude: 
- Research the best practices for writing pseudocodes, e.g. should it contain scaling, should it and what kind of commets should it have, what inputs should be written out explicitly, is it correct in its current form or we need more arguments that are missing or less,  etc.
- Make the DissipativeStep a bit more compact by introducing $U$, the block encoding of the jumps just like in Chen et al papers. This was the undoing step is just a controlled $U^\dagger$ controlled by the weak measurement ancilla being 0. $U$ should probably not fully be defined within the pseudocode but in some way it should be made clear even without reading the text around the pseudocode what $U$ is really. Maybe by a comment that would suggest what the effect on the system qubits was or something. Should be clean. 
- Another rigorous cross check about the maths that the pseudocodes we write here, combined, really gives what we want, what CKG and CKBG papers suggest. (E.g. I think previously it was wrongly stated by you, but we apply an RY with 1-\gamma onto the Boltzmann qubit, and all controls assume the control to be on zero, there and also for the weak measurement if I am not mistaken.)The Dissipative part is indeed the same, though maybe written with Trotterized evolutions like now. The coherent step was never really written out in detail like this, neither the full block encoding, nor using a QSP for it. But in CKG paper we see how part of it is meant to be implemented in Figure 4. or Proposition III.1 in the CKG paper. That part should be the same, but we definitely want to make it clear (that they didnt) how the full block encoding looks like for the nested integral.
  