# Algorithm Subsection -- Draft Prose

> **Target location**: `supplementary-informations/2_methods.tex`, lines 1362--1557.
> **Style**: first-person plural ("we"), run-in subsection headers in bold, matching Ch 5 conventions.
> **Notation**: follows established thesis notation ($\mathcal{L}$, $\beta$, $\delta$, $\sigma$, $\gamma$, $B_a$, etc.).
> **References**: CKBG23 = `\cite{chen2023quantum}`, CKG23 = `\cite{chen2023efficient}`, MW24 = `\cite{motlagh2024generalized}`, Chi+21 = `\cite{childs2021theory}`, LCU = `\cite{childs2012hamiltonian}` or `\cite{berry2015simulating}`.

---

## Prose to appear before and around Algorithm 1

We now assemble the subroutines developed in the preceding subsections into a complete quantum algorithm for Gibbs state preparation via the CKG Lindbladian (5.30).  The algorithm simulates the continuous-time Markov semigroup $e^{t\mathcal{L}}$ by repeated application of small-step channels that approximate $e^{\delta\mathcal{L}_a}$ for each jump $a \in \mathcal{A}$.  Its structure draws on three key ingredients from the literature: the *weak-measurement simulation scheme* of Chen, Kastoryano, Brand\~ao, and Gily\'en \cite{chen2023quantum}, which implements the dissipative channels $e^{\delta\mathcal{L}_{a,\mathrm{diss}}}$ at first-order accuracy in $\delta$; the *coherent correction term* $B_a$ introduced by Chen, Kastoryano, and Gily\'en \cite{chen2023efficient}, whose Hamiltonian simulation restores exact KMS detailed balance at the continuous level; and *Generalized Quantum Signal Processing* (GQSP) from Motlagh and Wiebe \cite{motlagh2024generalized}, which provides an efficient method to realise $e^{-i\delta B_a}$ from a block encoding of $B_a$.

The outer loop of the algorithm (Algorithm 1) iterates $L = \lceil t_{\mathrm{mix}}(\mathcal{L})\log(2/\varepsilon)/\delta\rceil$ Lindbladian steps, each of which sweeps over the full jump set $\mathcal{A}$.  Within each sweep, every jump $a$ is processed by first applying a CoherentStep (Algorithm 3) that implements $e^{-i\delta B_a}$, followed by a DissipativeStep (Algorithm 2) that implements $e^{\delta\mathcal{L}_{a,\mathrm{diss}}}(\rho) + \mathcal{O}(\delta^2)$ via weak measurement.  This ordering corresponds to the first-order coherent--dissipative splitting analysed in the preceding subsection (see equations (5.90)--(5.91)), whose $\mathcal{O}(\delta^2)$ error per step is already matched by the weak-measurement error and therefore does not require upgrading to a Strang splitting.  The jump-wise sweep itself is a first-order Lie--Trotter splitting of the full generator $\mathcal{L} = \sum_a \mathcal{L}_a$ as discussed in (5.82)--(5.86), which preserves the Gibbs state as fixed point exactly.

\medskip

*[Insert Algorithm 1 here.]*

\medskip

A note on the jump-wise decomposition: as discussed in the Generator splitting subsection, one may either sweep sequentially over the jumps (as written in Algorithm 1) or sample each jump $a$ uniformly at random from $\mathcal{A}$ and rescale by $M_\mathcal{A} := |\mathcal{A}|$.  In the sequential case, $\delta$ denotes the step size per individual jump and the effective Lindbladian time advanced per outer step is $M_\mathcal{A}\delta$.  In the random-sampling case, the same $\delta$ plays the role of the full step size and KMS detailed balance is preserved on average.  The two strategies have the same asymptotic complexity; we adopt the sequential sweep throughout for concreteness.

---

## Prose to appear before and around Algorithm 2

**Dissipative step: weak-measurement simulation.**  The DissipativeStep (Algorithm 2) implements a single channel application of $e^{\delta\mathcal{L}_{a,\mathrm{diss}}}(\rho) + \mathcal{O}(\delta^2)$ using the weak-measurement technique introduced in \cite[Theorem III.1]{chen2023quantum}.  The underlying principle is a quantum analogue of the Zeno effect: for a small step size $\delta$, the probability of a "jump" occurring -- i.e., the ancilla qubit $q_\delta$ being measured in the $|1\rangle$ state -- is only $\mathcal{O}(\delta)$.  With the complementary probability $1 - \mathcal{O}(\delta)$, no jump occurs and the system state is approximately unperturbed (the forward and reverse unitaries cancel).  After tracing out all ancilla registers, the effective channel on the system is
$$
\rho \;\mapsto\; \rho + \delta\,\mathcal{L}_{a,\mathrm{diss}}(\rho) + \mathcal{O}(\delta^2),
$$
which approximates $e^{\delta\mathcal{L}_{a,\mathrm{diss}}}(\rho)$ to first order.  The error is $\mathcal{O}(\delta^2)$ per step; over the $L = \mathcal{O}(t_{\mathrm{mix}}/\delta)$ total steps, the accumulated error is $\mathcal{O}(t_{\mathrm{mix}}\cdot\delta)$.  Setting $\delta = \varepsilon/t_{\mathrm{mix}}$ then yields $L = \mathcal{O}(t_{\mathrm{mix}}^2/\varepsilon)$ total channel applications.

The circuit proceeds in three phases: a forward pass that constructs the Trotterized OFT (5.80), a Boltzmann filter and weak-measurement rotation, and a reverse pass that uncomputes the ancilla entanglement.  We describe each in turn.

\medskip

**Forward pass: Trotterized OFT.**  The goal of the forward pass is to create the joint state $\sum_{\bar\omega}|\bar\omega\rangle\otimes\tilde{A}_a(\bar\omega)|\psi\rangle$ on the frequency register $\Omega$ and the system register $S$, where $\tilde{A}_a(\bar\omega)$ is the Trotterized discretised OFT from equation (5.80).  This is achieved in four sub-steps:

1. *State preparation* (step 2): The frequency register $\Omega$ (consisting of $r$ qubits) is initialised in the state $|f\rangle := \sum_{\bar t \in S_{t_0}^{[N]}} \bar f(\bar t)\,|\bar t\rangle$, where $\bar f(\bar t) := \sqrt{t_0}\,f(\bar t)$ absorbs the Riemann-sum prefactor into the discretised filter amplitudes (cf.\ equation (5.48)).  For the Gaussian filter (5.3), the amplitudes $\bar f(\bar t)$ form a discretised Gaussian distribution.  Efficient quantum state preparation of such structured distributions has been studied in \cite{mcArdle2022quantum}; the cost is $\mathcal{O}(\mathrm{poly}(\log N,\,\log(1/\varepsilon)))$ gates, where $N = 2^r$ is the register size.  For the circuit-level details, we refer to Section 4 (Quantum Circuits).

2. *Controlled Hamiltonian simulation* (steps 3 and 5): Controlled on the time register $|\bar t\rangle$, we apply the Trotterized evolutions $S_p^{(M)}(\bar t)^\dagger \approx e^{-iH\bar t}$ and $S_p^{(M)}(\bar t) \approx e^{+iH\bar t}$ to the system.  The product formula order $p$ and the number of Trotter steps $M$ are chosen according to the error analysis of Proposition 9.

3. *Jump operator* (step 4): Between the two controlled evolutions, the jump operator $A_a$ is applied to the system register.  For the canonical choice of single-site Pauli jumps, this is a single one- or two-qubit gate.

4. *Quantum Fourier Transform* (step 6): The QFT on $\Omega$ maps the time-domain register to the frequency domain, yielding $|\bar\omega\rangle$ labels.  The combined effect of steps 2--6 produces the joint state $\sum_{\bar\omega}|\bar\omega\rangle\otimes\tilde{A}_a(\bar\omega)\,|\psi\rangle$, which is the discretised and Trotterized analogue of the continuous OFT (5.1).

\medskip

**Boltzmann filter.**  Step 7 implements the transition weight $\gamma(\bar\omega)$ by performing a controlled rotation on the Boltzmann ancilla qubit $q_\gamma$.  Controlled on the frequency register $|\bar\omega\rangle$, we apply $R_Y(2\arcsin\sqrt{1-\gamma(\bar\omega)})$ to $q_\gamma$, producing
$$
|0\rangle_{q_\gamma} \;\longrightarrow\; \sqrt{\gamma(\bar\omega)}\,|0\rangle + \sqrt{1-\gamma(\bar\omega)}\,|1\rangle.
$$
(Our convention is that the $|0\rangle$ outcome of $q_\gamma$ corresponds to "acceptance" of the Boltzmann weight, which is the same as the convention in \cite[Figure 4]{chen2023quantum}.)

The controlled rotation $R_Y(2\arcsin\sqrt{1-\gamma(\bar\omega)})$ is parametrised by the classical function $\gamma$ evaluated at each grid point $\bar\omega$.  A naive decomposition into $N = 2^r$ individually controlled rotations would cost $\mathcal{O}(N)$ gates -- exponentially many in $r$.  An efficient implementation depends on the smoothness of $\gamma$:

- *Gaussian weight $\gamma_G$*:  Since $\gamma_G(\omega) = \exp(-((\omega+\omega_\gamma)^2)/(2\sigma_\gamma^2))$ is a smooth Gaussian, one can use controlled Hamiltonian simulation for the $Y$-rotation (implementing $\sum_{\bar\omega}e^{i\bar\omega Y/(2\|H\|)}\otimes|\bar\omega\rangle\langle\bar\omega|$) followed by QSVT to compose the polynomial $\arcsin\sqrt{1-\gamma_G(\bar\omega)}$.  The total cost is $\mathcal{O}(\mathrm{poly}(\log N, \log(1/\varepsilon)))$ as noted in \cite[footnote 33]{chen2023quantum}.

- *Metropolis weight $\gamma_M$* ($s = 0$):  The kink at $\omega = -\beta\sigma^2/2$ renders $\gamma_M$ merely Lipschitz (not $C^1$), which necessitates either a higher-degree polynomial approximation in QSVT or a piecewise implementation using a comparator circuit to distinguish the two sides of the kink.  Either approach incurs additional polynomial overhead, consistent with the broader theme that the kink degrades all aspects of the implementation (see the Quadrature errors for $\mathcal{L}_\mathrm{diss}$ subsection).

- *Smooth Metropolis weight $\gamma_M^{(s)}$* ($s > 0$):  Being Gevrey-$1/2$ by Proposition 7, $\gamma_M^{(s)}$ admits efficient polynomial approximation via QSVT, with a degree that scales polylogarithmically in $1/\varepsilon$.  The cost is comparable to the Gaussian case up to constants that depend on $s$.  This is the principal circuit-level advantage of the smooth Metropolis family: it combines the broad acceptance window of the Metropolis function with the efficient implementability of the Gaussian.

For circuit-level details of the QSVT-based controlled rotation, we refer to Section 4 (Quantum Signal Processing).

\medskip

**Weak measurement.**  In step 8, a second controlled rotation on the weak-measurement ancilla $q_\delta$ encodes the step size: controlled on $q_\gamma = |0\rangle$ (i.e., conditional on the Boltzmann acceptance), we apply $R_Y(2\arcsin\sqrt{\delta})$, producing
$$
|0\rangle_{q_\delta} \;\longrightarrow\; \sqrt{1-\delta}\,|0\rangle + \sqrt{\delta}\,|1\rangle
$$
in the accepted branch.  The net probability of a "jump" -- i.e., both $q_\gamma = |0\rangle$ and $q_\delta = |1\rangle$ -- is $\gamma(\bar\omega)\cdot\delta$, which is precisely the discretised Lindbladian transition rate for frequency $\bar\omega$.

\medskip

**Reverse pass and uncomputation.**  Steps 9--14 undo the forward circuit in reverse order: the Boltzmann rotation, QFT, controlled Hamiltonian evolutions, jump operator, and state preparation are all applied in their inverse form.  Crucially, these uncomputation steps are applied *unconditionally* -- they act on both the no-jump branch ($q_\delta = |0\rangle$) and the jump branch ($q_\delta = |1\rangle$).

This requires some explanation, as it differs from the original presentation in \cite[Theorem III.1]{chen2023quantum} where the reverse unitary is applied controlled on $q_\delta = |0\rangle$.  In the unconditional variant:
- In the *no-jump branch* ($q_\delta = |0\rangle$, probability $\approx 1 - \mathcal{O}(\delta)$): the full reverse circuit $U^\dagger U = I$ gives exact uncomputation, returning the system to its original state and the ancillas to $|0\rangle$.
- In the *jump branch* ($q_\delta = |1\rangle$, probability $\mathcal{O}(\delta)$): the reverse circuit does not cleanly uncompute -- the ancilla registers $\Omega$ and $q_\gamma$ remain entangled with the system.

The key observation is that after tracing out (measuring and discarding) the ancillas in step 15, the effective channel on the system is the same in both the controlled and unconditional variants, up to $\mathcal{O}(\delta^2)$ corrections.  This follows because the two implementations differ only in the jump branch, which has probability $\mathcal{O}(\delta)$, and within that branch the difference in the system's reduced state is $\mathcal{O}(\delta)$ (from the failure of the uncomputation to cancel).  Hence the total channel differs by $\mathcal{O}(\delta)\times\mathcal{O}(\delta) = \mathcal{O}(\delta^2)$, which is within the existing error budget.  The unconditional variant has the practical advantage of halving the controlled-operation overhead (no need for a multi-controlled reverse pass).

\medskip

*[Insert Algorithm 2 here.]*

---

## Prose to appear before and around Algorithm 3

**Coherent step: block encoding and Hamiltonian simulation via GQSP.**  The CoherentStep (Algorithm 3) implements the unitary evolution $e^{-i\delta B_a}$ generated by the per-jump coherent term $B_a$.  This is the component that distinguishes the CKG algorithm from the CKBG construction: without it (i.e., setting $B_a = 0$ for all $a$), the algorithm reduces to the CKBG Gibbs sampler \cite{chen2023quantum}, which satisfies only approximate GNS detailed balance.  The coherent term $B_a$ was derived in the KMS Lindbladian subsection (equation (5.29)) precisely to cancel the anti-Hermitian residual in the quantum discriminant, upgrading the Lindbladian from approximate GNS to exact KMS detailed balance.  Of course, the discretisation and Trotterisation errors analysed in the preceding subsections will re-introduce controlled violations of this detailed balance.

The algorithm proceeds in two stages: first we construct a *block encoding* $U_{B_a}$ satisfying $\langle 0|_{T_-}\langle 0|_{T_+}\,U_{B_a}\,|0\rangle_{T_-}|0\rangle_{T_+} = B_a/\alpha$, where $\alpha = \|b_-\|_1\,\|b_+\|_1$ is the subnormalisation factor; then we use GQSP to implement $e^{-i\delta B_a}$ from $U_{B_a}$.

\medskip

**Block encoding via Linear Combination of Unitaries.**  The block encoding is constructed from the nested double-integral representation of $B_a$ (equation (5.61), equivalently (5.39)):
$$
B_a = \sum_a \int_{-\infty}^{\infty} b_-(t)\,e^{-iHt/\sigma}\left[\int_{-\infty}^{\infty} b_+(t')\,A^{a\dagger}(\beta t')\,A^a(-\beta t')\,\dd t'\right]e^{iHt/\sigma}\,\dd t.
$$
The factorisation into outer ($b_-$) and inner ($b_+$) kernels, which arose from the Kossakowski matrix splitting $\alpha_{\nu_1,\nu_2} = \hat f_+(\nu_+)\hat f_-(\nu_-)$ (equation (5.38)), is essential for the circuit structure.  The LCU framework \cite{childs2012hamiltonian,berry2015simulating} encodes this double integral as a state-preparation / controlled-evolution / un-preparation sandwich using two ancilla registers $T_-$ ($r_-$ qubits) and $T_+$ ($r_+$ qubits).

Steps 2 and 4 prepare the outer and inner time registers:
$$
|0\rangle_{T_-} \;\longrightarrow\; |b_-\rangle := \sum_{\bar t} \sqrt{\frac{t_0\,|b_-(\bar t)|}{\alpha_-}}\;\mathrm{sgn}(b_-(\bar t))\,|\bar t\rangle, \qquad \alpha_- = \|b_-\|_1,
$$
$$
|0\rangle_{T_+} \;\longrightarrow\; |b_+\rangle := \sum_{\bar t'} \sqrt{\frac{t_0'\,|b_+(\bar t')|}{\alpha_+}}\;\mathrm{sgn}(b_+(\bar t'))\,|\bar t'\rangle, \qquad \alpha_+ = \|b_+\|_1.
$$
Here the absolute values go into the amplitudes and the signs (or, more precisely, the phases) go into the phases of the quantum state.  For the outer kernel $b_-$ (equation (5.40)), which is real-valued, $\mathrm{sgn}(b_-(\bar t)) \in \{-1, +1\}$ is the ordinary sign function.  For the inner kernel $b_+^{(s,\eta)}$ (equation (5.43)), which is complex-valued in general, the notation $\mathrm{sgn}(b_+(\bar t'))$ should be understood as the unit-modulus phase $b_+(\bar t')/|b_+(\bar t')|$.  For the Gaussian case (equation (5.41)), $b_+$ is also complex (it involves $e^{-2it\beta\omega_\gamma}$), and the same convention applies.

The state preparation costs for $|b_-\rangle$ and $|b_+\rangle$ are of the same type as the StatePrep in the dissipative step.  Since $b_-$ is a smooth, rapidly decaying function with effective support $\mathcal{O}(\beta\sigma)$, and $b_+^{(s,\eta)}$ has effective support $\mathcal{O}(1/(\sigma\beta\sqrt{1+s}))$ (with the $\eta$-regularised form from (5.43)), both admit efficient preparation circuits of cost $\mathcal{O}(\mathrm{poly}(r_\pm, \log(1/\varepsilon)))$ using the methods of \cite{mcArdle2022quantum}.  We defer the circuit-level details to Section 4.

Between the two state preparations, controlled Hamiltonian evolutions realise the operator-valued integrands.  Steps 3 and 11 provide the outer conjugation $e^{\pm iH\bar t/\sigma}$, while steps 5--9 implement the inner Heisenberg-picture expression
$$
\tilde{A}^{a\dagger}_a(\beta\bar t')\,\tilde{A}_a(-\beta\bar t') = S_2^{(M_+)}(\beta\bar t')\,A_a^\dagger\,S_2^{(M_+)}(-2\beta\bar t')\,A_a\,S_2^{(M_+)}(\beta\bar t'),
$$
which is controlled on the inner register $T_+$.  The Strang splitting ($p = 2$) is used here to maintain the palindromic structure required by Remark 11 for the preservation of the OFT adjoint identity.

After uncomputing the state preparations (steps 10 and 12), the net effect on the system register is precisely the block-encoded action $B_a/\alpha$, with $\alpha = \alpha_-\,\alpha_+ = \|b_-\|_1\,\|b_+\|_1$.

\medskip

**Subnormalisation $\alpha$ and its dependence on $\gamma$.**  The subnormalisation $\alpha = \|b_-\|_1\,\|b_+\|_1$ is a critical parameter: it enters both the QSP degree (see below) and the success probability of the block encoding.  Its value depends on the choice of transition weight:

- *Gaussian $\gamma_G$*: From (5.40) and (5.41), $\|b_-\|_1 \leq (\pi/\sqrt{2})\,e^{\beta^2\sigma^2/8}$ and $\|b_+\|_1 = (\beta/\omega_\gamma)\sqrt{\sigma_\gamma/(2\pi)}$.  With the natural KMS parametrisation $\sigma = \omega_\gamma = \sigma_\gamma = 1/\beta$, both norms are $\mathcal{O}(1)$, giving $\alpha = \mathcal{O}(1)$.

- *Smooth Metropolis $\gamma_M^{(s)}$*: The outer norm $\|b_-\|_1$ is universal (independent of $\gamma$) and remains $\mathcal{O}(1)$ for $\sigma = 1/\beta$.  The inner norm $\|b_+^{(s,\eta)}\|_1$ acquires a logarithmic divergence from the $\eta$-regularisation (equation (5.46)):
$$
\|b_-\|_1\,\|b_+^{(s,\eta)}\|_1 = \mathcal{O}\!\left(\log\frac{\beta\|H\|\,\|\sum_a A^{a\dagger}A^a\|}{\varepsilon}\right).
$$

This logarithmic overhead is the price paid for the broader Metropolis acceptance window and the $\eta$-regularisation of the $1/t$ singularity in $b_+$.  The resulting Hamiltonian simulation time for the coherent term is $\mathcal{O}(\beta\log(1/\varepsilon))$ for the Metropolis-like transitions, compared to $\mathcal{O}(\beta)$ for the Gaussian case.

\medskip

**Hamiltonian simulation via GQSP.**  Given the block encoding $U_{B_a}$ with $\langle 0|U_{B_a}|0\rangle = B_a/\alpha$, we wish to implement the unitary $e^{-i\delta B_a}$.  If $\lambda$ is an eigenvalue of $B_a/\alpha$ (so $|\lambda| \leq 1$), then the eigenvalues of $U_{B_a}$ are $e^{\pm i\theta}$ where $\sin\theta = \lambda$.  Thus, implementing $e^{-i\delta B_a} = e^{-i\delta\alpha\sin\theta}$ on the eigenspaces of $U_{B_a}$ reduces to finding a Laurent polynomial $P(e^{i\theta}) \approx e^{-i\delta\alpha\sin\theta}$ on the unit circle.

We use *Generalised Quantum Signal Processing* (GQSP) of Motlagh and Wiebe \cite{motlagh2024generalized} rather than standard QSP \cite{low2017optimal,low2019hamiltonian} for three reasons:
1. GQSP imposes no parity constraints on the polynomial $P$ -- important since $e^{-i\delta\alpha\sin\theta}$ has indefinite parity.
2. GQSP only requires controlled-$U_{B_a}$ (not $U_{B_a}^\dagger$), simplifying the circuit.
3. GQSP handles Laurent polynomials on the unit circle directly, which is the natural setting for our block encoding.

The degree of the Laurent polynomial required for $\varepsilon_\mathrm{QSP}$-approximation of $e^{-i\delta\alpha\sin\theta}$ is given by the GQSP Hamiltonian simulation bound \cite[Section V]{motlagh2024generalized}:
$$
d = \mathcal{O}\!\left(\delta\,\alpha + \frac{\log(1/\varepsilon_\mathrm{QSP})}{\log\log(1/\varepsilon_\mathrm{QSP})}\right).
$$
The first term $\delta\alpha$ is the dominant contribution when the simulation time $\delta\alpha$ is large, representing the intrinsic cost of the Hamiltonian simulation.  The second term captures the cost of achieving high precision and is nearly linear in $\log(1/\varepsilon_\mathrm{QSP})$.

Each call to $U_{B_a}$ within the GQSP protocol requires one full forward-and-reverse pass through the block encoding (steps 2--12), involving $\mathcal{O}(M_- + M_+)$ Trotter steps for the controlled Hamiltonian simulations plus the two state preparation subroutines.

\medskip

*[Insert Algorithm 3 here.]*

---

## Error budget and total complexity

We now consolidate the error sources developed individually in the preceding subsections into a complete error budget for the CKG Gibbs sampler.  Let $\varepsilon > 0$ be the target trace-distance accuracy $\|\rho_\mathrm{out} - \rho_\beta\|_1 \leq \varepsilon$.

\medskip

**Error sources.**  The total error receives contributions from five distinct mechanisms, each of which was analysed in its own subsection:

1. *Weak-measurement error* ($\mathcal{O}(\delta^2)$ per step):  The first-order simulation of $e^{\delta\mathcal{L}_{a,\mathrm{diss}}}$ via the weak-measurement scheme incurs an $\mathcal{O}(\delta^2)$ local error per channel application.  Over $L \cdot M_\mathcal{A}$ total applications, the accumulated error is $\mathcal{O}(L\cdot M_\mathcal{A}\cdot\delta^2) = \mathcal{O}(t_\mathrm{mix}\cdot M_\mathcal{A}\cdot\delta)$.  This is the dominant error source and sets the baseline scaling.

2. *Generator splitting -- jump-wise* ($\mathcal{O}(\delta^2\log^2(\beta\|H\|))$ per outer step, equations (5.82)--(5.86)):  The sequential Lie--Trotter splitting of $\mathcal{L} = \sum_a\mathcal{L}_a$ preserves the Gibbs state as fixed point exactly; the $\mathcal{O}(\delta^2)$ error per outer step accumulates to $\mathcal{O}(t_\mathrm{mix}\cdot\delta\cdot\log^2(\beta\|H\|))$ over the full evolution.

3. *Generator splitting -- coherent-dissipative* ($\mathcal{O}(\delta^2)$ per step, equations (5.90)--(5.93)):  The first-order splitting of $\mathcal{L}_a = \mathcal{L}_{a,\mathrm{coh}} + \mathcal{L}_{a,\mathrm{diss}}$ introduces a fixed-point bias of $\mathcal{O}(t_\mathrm{mix}\cdot\delta)$ over the full evolution.  As argued in the Generator splitting subsection, this is within the weak-measurement error budget.

4. *Quadrature / discretisation errors* (Table 5.1, equations (5.53), (5.60), (5.63)--(5.78)):  Truncation and aliasing errors from replacing the continuous integrals by finite sums on the frequency and time registers.  For the dissipative part, these are controlled to $\varepsilon_\mathrm{quad}$ by choosing $r = \mathcal{O}(\log(|\mathcal{A}|(\|H\|+1)(\sigma+1/\sigma)(\beta+1)) + s^{-1}\log^2(1/\varepsilon_\mathrm{quad}))$ qubits in the case of the smooth Metropolis $\gamma_M^{(s)}$ (or the same expression without the $1/s$ factor for the Gaussian $\gamma_G$).  For the coherent term $B$, the outer and inner integrals require $r_-$ and $r_+$ qubits respectively (equations (5.68) and (5.78)).

5. *Trotterisation errors* (Propositions 9 and 10):  Replacing exact Hamiltonian evolutions $e^{\pm iH\bar t}$ with product formulas $S_p^{(M)}(\bar t)$ introduces errors controlled by $\tilde\alpha_\mathrm{comm}^{(2)}\,|\bar t|^3/M^2$ for Strang splitting.  These errors affect both the OFT (Proposition 9) and the coherent term (Proposition 10), and are controlled by increasing $M$.

\medskip

**Parameter choices.**  To achieve a total error of $\varepsilon$, we allocate the error budget equally among the five sources (each target $\varepsilon/5$, or in practice simply ensuring each is $\mathcal{O}(\varepsilon)$).  The key parameter choices are:

- *Step size*: $\delta = \Theta(\varepsilon / (t_\mathrm{mix} \cdot M_\mathcal{A}))$, ensuring the weak-measurement and splitting errors are each $\mathcal{O}(\varepsilon)$.
- *Number of outer steps*: $L = \lceil t_\mathrm{mix}\log(2/\varepsilon)/\delta\rceil = \mathcal{O}(t_\mathrm{mix}^2\,M_\mathcal{A}\log(1/\varepsilon)/\varepsilon)$.
- *Frequency register qubits* ($r$ for the dissipative OFT): polylogarithmic in $1/\varepsilon$ for Gaussian or smooth Metropolis weights (see Table 5.1).
- *Time register qubits* ($r_-, r_+$ for the coherent term): polylogarithmic in $1/\varepsilon$ (equations (5.68), (5.78)).
- *Trotter steps* ($M, M_\pm$): chosen so that the Trotter error is $\mathcal{O}(\varepsilon)$; by Propositions 9 and 10, this requires $M = \mathcal{O}(\mathrm{poly}(\beta, \sigma^{-1}, \varepsilon^{-1/2}))$ for Strang splitting.
- *QSP degree*: $d = \mathcal{O}(\delta\alpha + \log(1/\varepsilon)/\log\log(1/\varepsilon))$.

\medskip

**Total step count and per-step cost.**  The total number of Lindbladian steps (outer loop iterations times jumps per iteration) is
$$
L \cdot M_\mathcal{A} = \mathcal{O}\!\left(\frac{t_\mathrm{mix}^2\,M_\mathcal{A}^2\,\log(1/\varepsilon)}{\varepsilon}\right).
$$
Each such step consists of one CoherentStep and one DissipativeStep.  The dominant cost within each step is the Hamiltonian simulation time:

- *DissipativeStep*: The controlled Trotterized OFT requires Hamiltonian simulation time $\mathcal{O}(t_0 N / 2) = \mathcal{O}(\sigma\sqrt{\log(1/\varepsilon)})$ (set by the time-domain truncation), distributed over $N = 2^r$ controlled time slices with $M$ Trotter steps each.
- *CoherentStep*: The block encoding requires outer Hamiltonian simulation time $\mathcal{O}(\beta\log(1/\varepsilon))$ (from the $b_-$ truncation window) and inner Hamiltonian simulation time $\mathcal{O}(\beta\sqrt{\log(1/\varepsilon)})$ for the Gaussian case or $\mathcal{O}(\beta\sqrt{\log(1/\varepsilon)}/(\sigma\sqrt{1+s}))$ for the smooth Metropolis.  The GQSP protocol repeats the block encoding $d = \mathcal{O}(\delta\alpha + \log(1/\varepsilon)/\log\log(1/\varepsilon))$ times.

\medskip

**Ancilla budget.**  The algorithm uses the following ancilla registers:

| Register | Size | Purpose |
|:---------|:-----|:--------|
| $\Omega$ | $r$ qubits | Frequency/time register for the dissipative OFT |
| $q_\gamma$ | 1 qubit | Boltzmann weight ancilla |
| $q_\delta$ | 1 qubit | Weak-measurement ancilla |
| $T_-$ | $r_-$ qubits | Outer time register for the coherent term |
| $T_+$ | $r_+$ qubits | Inner time register for the coherent term |
| $q_\mathrm{QSP}$ | 1 qubit | GQSP signal qubit |
| Trotter workspace | $\mathcal{O}(K)$ qubits | Ancillas for product-formula simulation of $H = \sum_k H_k$ |

The total ancilla count is $r + r_- + r_+ + 3 + \mathcal{O}(K)$, where $K$ is the number of terms in the Hamiltonian decomposition.  For the smooth Metropolis with $\sigma = 1/\beta$:
$$
r = \mathcal{O}\!\left(\log\!\left(|\mathcal{A}|(\|H\|+1)(\beta+1)^2\right) + \frac{1}{s}\log^2(1/\varepsilon)\right),
$$
$$
r_- = \mathcal{O}\!\left(\log\!\left(\frac{\beta^2\|H\|\,|\mathcal{A}|\,\|b_+\|_1\,\log(1/\varepsilon)}{\varepsilon}\right)\right), \qquad r_+ = \mathcal{O}\!\left(\log\!\left(\frac{\|H\|\,\|\sum_a A^{a\dagger}A^a\|\,|\mathcal{A}|\,\log(1/\varepsilon)}{\sigma(1+s)\varepsilon}\right)\right).
$$
All three register sizes scale polylogarithmically in $1/\varepsilon$, confirming that the ancilla overhead is modest.

\medskip

**Summary.**  Combining all the above, we state the complexity informally:

\medskip
\noindent\textit{Informal Proposition (CKG Gibbs sampler complexity).}  Given a system Hamiltonian $H = \sum_{k=1}^K H_k$ on $n$ qubits and a primitive CKG Lindbladian $\mathcal{L}$ with mixing time $t_\mathrm{mix}$, Algorithm 1 prepares a state $\varepsilon$-close to $\rho_\beta$ in trace distance using:
- $\mathcal{O}(t_\mathrm{mix}^2\,M_\mathcal{A}^2\,\log(1/\varepsilon)/\varepsilon)$ total Lindbladian step applications,
- $r + r_- + r_+ + 3 + \mathcal{O}(K)$ ancilla qubits (all polylogarithmic in $1/\varepsilon$),
- per-step Hamiltonian simulation time $\mathcal{O}(\beta\log(1/\varepsilon))$ for smooth Metropolis transitions,
- GQSP degree $d = \mathcal{O}(\delta\alpha + \log(1/\varepsilon)/\log\log(1/\varepsilon))$ calls to the block encoding per coherent step.

\noindent The $t_\mathrm{mix}^2/\varepsilon$ scaling arises from the combination of $L = \mathcal{O}(t_\mathrm{mix}/\delta)$ steps and $\delta = \mathcal{O}(\varepsilon/t_\mathrm{mix})$, which is inherent to the first-order weak-measurement scheme.

---

## Remarks

**Sequential vs.\ random jump selection.**  As noted above, the inner loop over jumps $a \in \mathcal{A}$ can be replaced by random sampling.  In the random-sampling variant, each outer step selects a single $a$ uniformly at random and applies the rescaled channel $e^{M_\mathcal{A}\delta\mathcal{L}_a}$.  This preserves KMS detailed balance on average (the convex combination of KMS DB channels is KMS DB), whereas the sequential sweep breaks it at $\mathcal{O}(\delta^2)$ through the Lie--Trotter commutator (equation (5.87)).  The spectral consequences are complementary: the sequential sweep's leading error is anti-KMS-adjoint and does not perturb the spectral gap (equation (5.89)), while the random-sampling error is symmetric and does.  Since both errors are $\mathcal{O}(\delta^2)$ per outer step, neither dominates asymptotically, and the choice is largely a matter of implementation preference.

\medskip

**Comparison with the CKBG compressed algorithm.**  Chen et al.\ \cite[Theorem III.2]{chen2023quantum} also describe a *compressed* variant of the weak-measurement scheme that achieves near-linear scaling in $t_\mathrm{mix}$ (as opposed to the quadratic scaling $t_\mathrm{mix}^2/\varepsilon$ of the simple scheme we use).  The compressed scheme uses quantum amplitude estimation to batch multiple weak-measurement steps, reducing the total number of Hamiltonian simulation calls.  We opt for the simple (uncompressed) scheme for several reasons: (i) it has smaller constant factors and is more amenable to near-term implementation; (ii) its error analysis composes transparently with the other error sources (Trotter, quadrature, generator splitting); (iii) for the system sizes accessible to current or near-term hardware, the quadratic-vs-linear distinction in $t_\mathrm{mix}$ is less important than the per-step cost.  If a near-linear scaling becomes essential, the compressed scheme can be straightforwardly substituted for the DissipativeStep, with the Trotter and quadrature analyses carrying over unchanged.

\medskip

**NISQ-friendly variant: Trotterizing the coherent evolution.**  An alternative to the QSP-based implementation of $e^{-i\delta B_a}$ is to directly Trotterize the coherent evolution.  Since $\mathcal{L}_{a,\mathrm{coh}}(\cdot) = -i[B_a, \cdot]$ generates a unitary channel, one could implement it via a product-formula simulation of $B_a$ itself, bypassing the block encoding, GQSP, and the associated ancilla registers $T_\pm$ and $q_\mathrm{QSP}$ entirely.  This would make the coherent and dissipative ancilla registers independent and simplify the circuit considerably.  The trade-off is that the product-formula simulation of $B_a$ -- which is itself defined through a double integral over Hamiltonian evolutions -- would incur Trotter errors that scale polynomially rather than logarithmically in $1/\varepsilon$, degrading the asymptotic $t_\mathrm{mix}$ scaling from the near-optimal $\mathcal{O}(t_\mathrm{mix}^2/\varepsilon)$ of the QSP-based scheme.  Nevertheless, for implementations prior to fault tolerance, where the constant factors and circuit depth matter more than asymptotic exponents, this Trotterized variant may be a practical alternative worth exploring.
