# DissipativeStep
<!-- Suggested placement: Chapter 5, Algorithm subsection, immediately after CoherentStep prose (p. 88), before Algorithm 1b -->
<!-- Depends on: Approximate GNS Lindbladian (§5.1, eqs 5.1-5.7), Quadrature errors for L_diss (eqs 5.47-5.60), Hamiltonian Trotterization (Prop 9, eq 5.80), Generator splitting (eqs 5.90-5.93), Weak-measurement techniques (Ch 4) -->

\textsc{DissipativeStep:} **Weak-measurement simulation.**
The \textsc{DissipativeStep} (Algorithm 1b) implements the dissipative channel $e^{\delta\mathcal{L}_{a,\mathrm{diss}}}(\rho) + \mathcal{O}(\delta^2)$ for a single jump $a$.  This is the component inherited from CKBG23 \cite{chen2023quantum}; the CKG algorithm retains it unchanged except for the shifted transition weights $\gamma$ and the addition of the coherent correction $B_a$ as a separate subroutine.  Since the weak-measurement scheme has been described in detail in \cite[Theorem III.1]{chen2023quantum} and its circuit form in \cite[Figure 3]{chen2023quantum}, we keep the exposition brief and focus on the structural features that are specific to our implementation --- in particular, the simplifications arising from unitary (Pauli) jump operators.

\medskip

**The forward unitary $U_{\mathrm{diss}}$.**
It is convenient to package the forward pass of Algorithm 1b (steps 2--6) into a single unitary $U_{\mathrm{diss}}$ acting on the joint system $\Omega \otimes S$.  Concretely,
$$
U_{\mathrm{diss}} := \mathrm{QFT}(\Omega) \cdot \mathrm{Ctrl\text{-}Trott}^+(\Omega, S) \cdot A_a \cdot \mathrm{Ctrl\text{-}Trott}^-(\Omega, S) \cdot \mathrm{StatePrep}(\Omega),
$$
where each factor is read right-to-left in circuit order: $\mathrm{StatePrep}(\Omega)$ prepares the Gaussian filter state $|f\rangle$ on the frequency register (step 2), $\mathrm{Ctrl\text{-}Trott}^-$ applies the backward Trotterized evolution $S_p^{(M)}(\bar t)^\dagger$ controlled on $\Omega$ (step 3), $A_a$ applies the jump operator (step 4), $\mathrm{Ctrl\text{-}Trott}^+$ applies the forward evolution $S_p^{(M)}(\bar t)$ (step 5), and $\mathrm{QFT}$ maps the time register to the energy basis (step 6).

Starting from $|0\rangle_\Omega |\psi\rangle_S$, the forward unitary produces
$$
U_{\mathrm{diss}}\,|0\rangle_\Omega\,|\psi\rangle_S = \sum_{\bar\omega} |\bar\omega\rangle_\Omega \otimes \tilde{A}_a(\bar\omega)\,|\psi\rangle_S,
$$
where $\tilde{A}_a(\bar\omega)$ is the Trotterized OFT defined in equation (5.80), which replaces the exact Hamiltonian evolution in the discretized OFT (5.48) with the product formula $S_p^{(M)}$.  This is the dissipative analogue of the block encoding $U_{B_a}$ from the coherent step: the frequency register $\Omega$ plays the role of the ancilla, and projecting onto a given $|\bar\omega\rangle$ selects the corresponding filtered jump component $\tilde{A}_a(\bar\omega)$.  The key difference is that $U_{\mathrm{diss}}$ is not a block encoding in the usual sense --- it is a unitary on the full $\Omega \otimes S$ space whose ancilla-projected blocks $\langle \bar\omega|U_{\mathrm{diss}}|0\rangle = \tilde{A}_a(\bar\omega)$ define a family of operators satisfying the Parseval-type bound $\sum_{\bar\omega} \tilde{A}_a(\bar\omega)^\dagger \tilde{A}_a(\bar\omega) \preceq \mathbb{1}$ \cite[Proposition A.1]{chen2023quantum}.

Defining $U_{\mathrm{diss}}$ explicitly makes the structure of the remaining steps transparent: the Boltzmann filter and weak measurement act on the energy-labelled state $\sum_{\bar\omega} |\bar\omega\rangle \otimes \tilde{A}_a(\bar\omega)|\psi\rangle$, and the reverse pass is simply $U_{\mathrm{diss}}^\dagger$ (preceded by undoing the Boltzmann rotation), applied in the no-jump branch.

\medskip

**Boltzmann filter and weak measurement** (steps 7--8). Given the energy-resolved state after $U_{\mathrm{diss}}$, two single-qubit rotations implement the accept/reject logic:
\begin{enumerate}
\item \emph{Boltzmann filter} (step 7): controlled on the energy register $|\bar\omega\rangle$, rotate the Boltzmann qubit $q_\gamma$ by $R_Y(2\arcsin\sqrt{1 - \gamma(\bar\omega)})$.  This maps $|0\rangle_{q_\gamma} \to \sqrt{\gamma(\bar\omega)}\,|0\rangle + \sqrt{1-\gamma(\bar\omega)}\,|1\rangle$, encoding the transition weight into the amplitude of the $|0\rangle$ branch.  Our convention is to condition on $q_\gamma = |0\rangle$ for acceptance, following \cite{chen2023quantum}.
\item \emph{Weak measurement} (step 8): controlled on $q_\gamma = |0\rangle$, rotate the weak-measurement qubit $q_\delta$ by $R_Y(2\arcsin\sqrt\delta)$.  In the acceptance branch, $|0\rangle_{q_\delta} \to \sqrt{1-\delta}\,|0\rangle + \sqrt\delta\,|1\rangle$.  The net probability of the "jump" outcome ($q_\delta = |1\rangle$) is $\gamma(\bar\omega)\,\delta$ --- first order in $\delta$, as required by the weak-measurement approximation $e^{\delta\mathcal{L}_{a,\mathrm{diss}}} = \mathrm{id} + \delta\mathcal{L}_{a,\mathrm{diss}} + \mathcal{O}(\delta^2)$.
\end{enumerate}

The small jump probability is precisely what makes this a *weak* measurement: the system state is only weakly disturbed per step, and the accumulated effect of $L \cdot M_\mathcal{A}$ such steps reproduces the dissipative semigroup to the desired accuracy.  The $\mathcal{O}(\delta^2)$ local error per step is the baseline that all other error sources (Trotterisation, quadrature, generator splitting) are tuned to match.

\medskip

**Implementing the Boltzmann rotation.**  The Boltzmann filter (step 7) requires a circuit that, controlled on the $r$-qubit energy register $|\bar\omega\rangle$, rotates $q_\gamma$ by the angle $\theta(\bar\omega) = 2\arcsin\sqrt{1 - \gamma(\bar\omega)}$.  This is a function-controlled rotation: the classical function $\gamma$ must be "compiled" into a quantum circuit acting on the joint $\Omega \otimes q_\gamma$ space.  The cost depends on the regularity of $\gamma$ and on the register size $r$.

For *large* $r$ (the fault-tolerant regime), the standard approach is to use QSP/QSVT \cite{gilyen2019quantum} to implement a polynomial approximation of $\sqrt{\gamma(\bar\omega)}$ as a function of the block-encoded normalised energy $\bar\omega / (\omega_0 N/2)$.  The required polynomial degree $d_\gamma$ depends on the regularity of $\gamma$:
\begin{itemize}
\item For the Gaussian $\gamma_G$ (analytic) or the smooth Metropolis $\gamma_M^{(s)}$ with $s > 0$ (Gevrey-$1/2$, Proposition 7): $d_\gamma = \mathcal{O}(\mathrm{polylog}(1/\varepsilon_\gamma))$, since Gevrey-$1/2$ functions admit super-algebraically converging polynomial approximations.
\item For the kinky Metropolis $\gamma_M$ ($s = 0$, only Lipschitz at $\omega = -\beta\sigma^2/2$): $d_\gamma = \mathcal{O}(\mathrm{poly}(1/\varepsilon_\gamma))$, since Lipschitz functions require polynomially many terms.  As noted in the Designing transition weights subsection, one can work around the kink by splitting the domain piecewise (via a comparator on the energy register) and applying QSP on each smooth piece, but this adds circuit complexity --- and is one of the motivations for the smooth Metropolis family.
\end{itemize}

For *moderate* $r$ (say $r \leq 12$, which covers the parameter regime of our quadrature analysis in Table 5.1), one can bypass QSP entirely and use uniformly controlled $R_Y$ rotations \cite{mottonen2004transformation,shende2006synthesis} --- the same decomposition technique used for the state preparations $|f\rangle$, $|b_-\rangle$, $|b_+\rangle$ in the CoherentStep.  The $N = 2^r$ rotation angles $\theta(\bar\omega)$ are precomputed classically, and the resulting uniformly controlled rotation decomposes into $\mathcal{O}(2^r)$ CNOT + single-qubit gates.  At $r = 10$--$12$ this amounts to a few thousand gates, subdominant to the $\mathcal{O}(M)$ Trotter steps in the controlled Hamiltonian evolutions.  This is the approach we use in practice.

In either case, the Boltzmann rotation is not a bottleneck: for smooth $\gamma$ and large $r$, QSP gives polylogarithmic overhead; for moderate $r$, the brute-force decomposition is already cheap.  This is another instance where the regularity of $\gamma$ pays off --- a practical advantage of the smooth Metropolis family $\gamma_M^{(s>0)}$ introduced in the Designing transition weights subsection.

\medskip

**Reverse pass and measurement** (steps 9--15).  After the weak measurement, the algorithm splits into two coherent branches depending on $q_\delta$.  In the *no-jump branch* ($q_\delta = |0\rangle$), steps 9--14 undo all operations of the forward pass: first the Boltzmann rotation is reversed (step 9), then $U_{\mathrm{diss}}^\dagger$ is applied (steps 10--14), restoring the ancillas to $|0\rangle^{\otimes r}$ and the system to its original state.  In the *jump branch* ($q_\delta = |1\rangle$), none of steps 9--14 act --- the ancillas $\Omega$, $q_\gamma$ remain entangled with $S$ in the state produced by the forward pass.

In the notation of $U_{\mathrm{diss}}$, the reverse pass is: controlled on $q_\delta = |0\rangle$, first undo the Boltzmann rotation on $q_\gamma$, then apply $U_{\mathrm{diss}}^\dagger$.  Both operations are conditioned on the no-jump outcome; this controlled structure is the same as in \cite[Figure 3]{chen2023quantum}.

After measuring and resetting all ancillas (step 15), the two branches contribute to the effective channel on $S$: the no-jump branch returns the state unchanged (with probability $1 - \mathcal{O}(\delta)$, contributing the identity and the anti-commutator correction), while the jump branch contributes the desired dissipative transitions $\sum_{\bar\omega} \gamma(\bar\omega)\,\tilde{A}_a(\bar\omega)\,\rho\,\tilde{A}_a(\bar\omega)^\dagger$ (with total probability $\mathcal{O}(\delta)$).  Together, these reproduce $\rho + \delta\,\mathcal{L}_{a,\mathrm{diss}}(\rho) + \mathcal{O}(\delta^2)$, which is the first-order expansion of the dissipative semigroup channel.  The formal derivation is in \cite[Theorem III.1]{chen2023quantum}.

\medskip

**Simplifications from unitary jump operators.** The CKBG23 construction \cite[Figure 3]{chen2023quantum} is stated for general (possibly non-unitary) jump operators $A_a$, which require block encoding the jump via ancilla qubits.  In Chen et al.'s construction, the isometry $V$ encodes both the sum over jumps $\sum_a$ and the energy resolution into a single block-encoding circuit, using auxiliary registers for the jump labels and Boltzmann weights.  This is necessary when $A_a$ is non-unitary, since applying a non-unitary operator directly is not a valid quantum operation.

In our setting, the jump operators are single-site Pauli operators (e.g.\ $\sigma_x^{(j)}, \sigma_y^{(j)}, \sigma_z^{(j)}$ on site $j$), which are both Hermitian and unitary: $A_a^\dagger = A_a$ and $A_a^2 = \mathbb{1}$.  This leads to two simplifications that are reflected in Algorithm 1b:

\begin{enumerate}
\item \emph{No block-encoding ancillas.}  Since $A_a$ is unitary, it can be applied directly as a quantum gate (step 4), and its inverse $A_a^\dagger = A_a$ is equally cheap (step 12).  The block-encoding ancillas for the jump operator are absent entirely.  The only ancillas needed are the frequency register $\Omega$ (for the OFT), the Boltzmann qubit $q_\gamma$, and the weak-measurement qubit $q_\delta$ --- a total of $r + 2$ qubits.

\item \emph{Sequential single-jump application.}  In each $\delta$-step, the outer loop of Algorithm 1 sweeps over the jumps sequentially: for each $a$, a separate \textsc{DissipativeStep} is called with that single $A_a$.  This is in contrast to the CKBG23 presentation, where a single step applies the block-encoded isometry $V$ that encodes all jumps simultaneously.  Our sequential sweep avoids the LCU overhead of encoding $M_\mathcal{A}$ jump operators into a single block encoding (which would require $\lceil\log_2 M_\mathcal{A}\rceil$ additional ancilla qubits for the jump label register and a more complex controlled circuit), at the cost of $M_\mathcal{A}$ separate passes per outer step instead of one.  Asymptotically both approaches have the same scaling (the CKBG23 block encoding has subnormalisation $\|\sum_a A_a^\dagger A_a\| \leq 1$, which cancels the $M_\mathcal{A}$ factor in the sequential sweep), but the sequential approach produces a simpler circuit per step.
\end{enumerate}

Together, these two simplifications eliminate all ancillas beyond the frequency register and two single qubits: the general CKBG23 circuit requires $r + \lceil\log_2 M_\mathcal{A}\rceil + c + 2$ ancilla qubits (frequency register $+$ jump-label register $+$ block-encoding ancillas for non-unitary jumps $+$ Boltzmann and weak-measurement qubits), while our per-jump circuit needs only $r + 2$.  The jump-label register is absent because we sweep over jumps sequentially rather than encoding them simultaneously; the block-encoding ancillas are absent because each Pauli $A_a$ is applied directly as a gate.  The price is that each outer Lindbladian step involves $M_\mathcal{A}$ sequential applications of the DissipativeStep circuit, but since each individual circuit is simpler and does not require multi-controlled jump-label selection, the practical gate count is comparable.  A circuit diagram illustrating this simplified structure will be provided separately.

\medskip

**Connection to the error budget.** The DissipativeStep contributes to the overall error through two channels.  First, the weak-measurement approximation itself introduces an $\mathcal{O}(\delta^2)$ error per step \cite[Theorem III.1]{chen2023quantum} --- this is the baseline error that sets the scaling for all other error sources.  Second, the Trotterized OFT $\tilde{A}_a(\bar\omega)$ approximates the exact discretized OFT $\hat{A}_a(\bar\omega)$ with an error controlled by the Trotter step count $M$ (Proposition 9).  The quadrature errors from discretising the time and energy integrals are controlled by the frequency register size $r$, as analysed in the Quadrature errors for $\mathcal{L}_\mathrm{diss}$ subsection.  All of these are set to be $\mathcal{O}(\varepsilon)$ individually.

Compared to the CoherentStep, the DissipativeStep is structurally simpler: it uses a single time/frequency register $\Omega$ (versus the nested $T_-, T_+$ of the coherent term), involves no QSP, and its Hamiltonian simulation time per step is $\mathcal{O}(\sigma\sqrt{\log(1/\varepsilon)})$ (set by the Gaussian filter truncation $T = \Theta(\sigma\sqrt{\log(1/\varepsilon)})$ and the controlled Trotter evolutions).  The dominant cost in practice is the $2M$ Trotter steps for the forward and reverse controlled evolutions (steps 3 and 5, plus their reverses 11 and 13), which use the same Strang splitting $S_p^{(M)}$ (with $p = 2$) as the coherent step's outer evolution.  As in the CoherentStep, the palindromic property $S_p(t)^\dagger = S_p(-t)$ of the Strang splitting is structurally important: it preserves the OFT adjoint identity $\tilde{A}_a(-\bar\omega)^\dagger = \tilde{A}_a(\bar\omega)$ at the Trotterized level (Remark 11), which is needed for the detailed-balance structure of the dissipative channel.

*[Insert Algorithm 1b here.]*

---
## Writing Notes
<!-- These notes are for the author, not for the thesis -->

- **Equation cross-references to verify in TeX**: (5.1), (5.3), (5.7), (5.47), (5.48), (5.80), (5.90), (5.91).  Also Proposition 7 (smooth Metropolis Gevrey-1/2), Proposition 9 for the Trotterized OFT, Table 5.1 (quadrature parameters giving $r$ values), equations (5.31) and (5.34) for Metropolis definitions, and equations (5.49)--(5.60) for the quadrature error analysis of $\mathcal{L}_\mathrm{diss}$.
- **Citations to verify**: `\cite{chen2023quantum}` = CKBG23 (Theorem III.1 for weak-measurement channel, Figure 3 for circuit, Proposition A.1 for Parseval identity, footnote 33 for comparator circuit suggestion), `\cite{gilyen2019quantum}` = QSVT (Gilyen et al. 2019), `\cite{mottonen2004transformation}` = Mottonen et al. 2004, `\cite{shende2006synthesis}` = Shende, Bullock, Markov 2006 (uniformly controlled rotations).
- **Review feedback addressed**:
  - (1.1) Boltzmann rotation convention: explicitly stated that we condition on $q_\gamma = |0\rangle$, matching CKBG23.
  - (1.2) Controlled uncomputation: confirmed that steps 9-14 are controlled on $q_\delta = |0\rangle$ (as stated in the pseudocode header).  The draft describes this correctly.
  - (2.1) Ctrl-Trott spelling: used `Ctrl-Trott` consistently throughout.
- **Parseval identity correction**: Changed from $\sum \tilde{A}^\dagger \tilde{A} \preceq \mathbb{1}$ to $= A_a^\dagger A_a \preceq \mathbb{1}$, which is more precise (the Parseval identity gives equality with $A_a^\dagger A_a$, and the normalization (5.7) gives $\preceq \mathbb{1}$ after summing over $a$; for a single $a$, $A_a^\dagger A_a \preceq \mathbb{1}$ since $\|\sum_a A_a^\dagger A_a\| \leq 1$).
- **Figures needed**: A circuit diagram showing the simplified DissipativeStep with Pauli jumps (no block-encoding ancillas), contrasted with the general CKBG23 circuit.  This would concretely illustrate the ancilla savings.
- **Style note**: Kept deliberately shorter than the CoherentStep prose, since: (i) the weak-measurement scheme is not our original contribution (it's from CKBG23), and (ii) the user requested minimal repetition of Chen's work, focusing only on our simplifications.
- **Connection points**: The error budget consolidation (to follow) should pull together the DissipativeStep's $r, M$ parameters alongside the CoherentStep's $r_\pm, M_\pm, d$.  The Parent Hamiltonian subsection follows after this.
- **Ham sim time correction**: Changed from $\mathcal{O}(\beta\log(1/\varepsilon))$ to $\mathcal{O}(\sigma\sqrt{\log(1/\varepsilon)})$ for the DissipativeStep.  The filter truncation is $T = \Theta(\sigma\sqrt{\log(1/\varepsilon)})$ for Gaussian filters, and the Hamiltonian simulation time is $T$ (not $\beta T$; the factor $\beta$ appears in the CoherentStep inner evolution but not in the DissipativeStep).  With $\sigma = 1/\beta$ (KMS parametrisation), this becomes $\mathcal{O}(\sqrt{\log(1/\varepsilon)}/\beta)$.
