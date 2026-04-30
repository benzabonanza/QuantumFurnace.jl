# Review: Algorithm Subsection (Ch 5, pp. 85--87)

**Target**: `supplementary-informations/2_methods.tex`, lines 1363--1557 (three `algorithm` environments: CKG-GibbsSampler, DissipativeStep, CoherentStep).

**Reviewer**: Claude (draft-check agent)
**Date**: 2026-03-31


---

## 0. Executive Summary

The Algorithm subsection currently contains **only pseudocode** -- three well-structured algorithm floats with no accompanying prose whatsoever.  The pseudocode is largely correct relative to the preceding theory, but has several notation inconsistencies, one correctness issue with the Boltzmann rotation angle, a dangling `[REF to text]` placeholder, and -- most critically -- is missing roughly 3--5 pages of prose that any reader (or examiner) will expect.  This review catalogues (1) pseudocode correctness issues, (2) notation inconsistencies, (3) the missing prose/explanation, and (4) missing citations and cross-references.

---

## 1. Pseudocode Correctness

### 1.1 Boltzmann rotation angle (Algorithm 2, step 7) -- POTENTIAL ERROR

The rotation on $q_\gamma$ is written as:

$$R_Y\!\bigl(2\arcsin\!\sqrt{1 - \gamma(\bar\omega)}\bigr)$$

with the comment claiming the resulting state is $\sqrt{\gamma(\bar\omega)}|0\rangle + \sqrt{1-\gamma(\bar\omega)}|1\rangle$.

**Issue**: An $R_Y(\theta)$ rotation applied to $|0\rangle$ gives $\cos(\theta/2)|0\rangle + \sin(\theta/2)|1\rangle$.  With $\theta = 2\arcsin\sqrt{1-\gamma}$ we get $\sin(\theta/2) = \sqrt{1-\gamma}$ and $\cos(\theta/2) = \sqrt{\gamma}$.  So the result is indeed $\sqrt{\gamma}|0\rangle + \sqrt{1-\gamma}|1\rangle$ as stated in the comment -- the comment is correct.

However, this is confusing notation: the natural convention is $R_Y(2\arcsin\sqrt{\gamma(\bar\omega)})$ which maps $|0\rangle \to \sqrt{1-\gamma}|0\rangle + \sqrt{\gamma}|1\rangle$, and then you would condition on $q_\gamma = |1\rangle$ for acceptance.  The current convention (condition on $|0\rangle$) works but is the opposite of CKBG23's convention where they use $Y_{1-\gamma(\bar\omega)}$ (see CKBG23, Figure 4 and equation (3.10), where the Boltzmann weight qubit stores $\sqrt{\gamma(\bar\omega)}|0\rangle + \sqrt{1-\gamma(\bar\omega)}|1\rangle$ and acceptance is on the $|0\rangle$ outcome).

**Verdict**: Mathematically consistent internally, but the convention choice needs to be stated explicitly in accompanying prose.  The `[REF to text]` placeholder in the comment on line 1439 needs to be resolved -- it should point to a prose paragraph explaining the controlled rotation and its relation to the CKBG23 block-encoding (Figure 4 in that paper, or the equivalent passage in the thesis once written).

### 1.2 Weak-measurement uncomputation conditioned on wrong qubit? (Algorithm 2, steps 9--14)

The section header says "Reverse: uncompute $U$ (controlled on $q_\delta = |0\rangle$)" but steps 9--14 are written as **unconditional** operations (no explicit conditioning on $q_\delta$).  The comment text says "No-jump branch ($q_\delta = |0\rangle$): full uncomputation; jump branch ($q_\delta = |1\rangle$): $\Omega, q_\gamma, S$ remain entangled."

**Issue**: In the CKBG23 weak-measurement scheme (Figure 3, Theorem III.1), the reverse $U^\dagger$ is applied **controlled on the ancilla qubit being 0**.  In the thesis pseudocode, the uncomputation steps 9--14 appear to be unconditional.  If they are truly unconditional (applied to both branches), then:
- In the no-jump branch ($q_\delta = |0\rangle$): the full circuit $U^\dagger U = I$ gives correct uncomputation.
- In the jump branch ($q_\delta = |1\rangle$): applying $U^\dagger$ does NOT give a clean state; the ancillas remain entangled with the system.

In CKBG23, the $U^\dagger$ is controlled on the ancilla being $|0\rangle$ (the "no jump" outcome), which means the jump branch is untouched and the ancillas are traced out.  The thesis pseudocode should either:
(a) Explicitly state that steps 9--14 are controlled on $q_\delta = |0\rangle$, or
(b) Note that unconditional application also works because after tracing out the measurement outcomes, the effective channel is the same to $O(\delta^2)$ -- but this equivalence is non-obvious and should be proven or cited.

**Recommendation**: This is a significant point requiring careful prose explanation.  The current header comment is the only hint and it is insufficient.  The unconditional-vs-controlled distinction affects the circuit depth by a factor of ~2 (controlled operations cost more) and also changes the error analysis.  This needs to be resolved.

### 1.3 Total number of Lindbladian steps (Algorithm 1, line 1)

$$L \gets \lceil t_{\mathrm{mix}}(\mathcal{L}) \log(2/\varepsilon) / \delta \rceil$$

This formula counts the total number of **outer** Lindbladian steps.  The factor $\log(2/\varepsilon)$ comes from equation (3.22) in the preliminaries (mixing time definition: $\|e^{tL}(\rho_0) - \rho_\beta\|_1 \leq 2e^{-t/t_\mathrm{mix}}$ requires $t = t_\mathrm{mix}\log(2/\varepsilon)$).  However, each outer step sweeps over all $M_\mathcal{A}$ jumps, so the total number of **channel applications** is $L \times M_\mathcal{A}$.  This is correct as written (the outer loop runs $L$ times, the inner loop $M_\mathcal{A}$ times), but the prose should clarify that the effective Lindbladian time per outer step is $\delta \cdot M_\mathcal{A}$ (if using sequential sweep) or $\delta$ (if using random sampling with rescaling).

**Issue**: The sequential sweep vs. random sampling distinction discussed in the Generator splitting subsection (p. 81) is referenced in the comment on line 1387 but not resolved algorithmically.  In the sequential sweep case, $\delta$ is the step size per *individual jump*, and the effective time per outer step is $M_\mathcal{A} \delta$.  In the random sampling case, $\delta$ would be the full step size.  The formula $L = \lceil t_\mathrm{mix} \log(2/\varepsilon)/\delta \rceil$ is only correct if $\delta$ has been appropriately defined.  This should be clarified.

### 1.4 CoherentStep block encoding (Algorithm 3)

The block encoding of $B_a/\alpha$ is constructed correctly via the LCU-like structure with $T_-$ and $T_+$ registers.  The state preparation uses $\sqrt{t_0 |b_-(\bar t)|/\alpha_-} \cdot \mathrm{sgn}(b_-(\bar t))$ as amplitudes, which correctly encodes the signed kernel into the quantum state (the absolute value goes into the amplitude, the sign goes into the phase).  This matches equation (5.61).

**Minor issue**: The state preparation formula on line 1500 includes a factor $\sqrt{t_0 |b_-(\bar t)|/\alpha_-}$ times $\mathrm{sgn}(b_-(\bar t))$.  But $\mathrm{sgn}$ returns $\pm 1$, while quantum amplitudes can be complex.  For the $b_-$ kernel, which is real-valued (as stated on p. 60), this is fine.  But for $b_+^{(s,\eta)}$ from equation (5.43), the kernel is complex-valued (it involves $e^{-\sigma^2\beta^2 t(2t+i)(1+s)}$).  The $\mathrm{sgn}(b_+(\bar t'))$ notation on line 1510 should be replaced with something like $b_+(\bar t')/|b_+(\bar t')|$ (i.e., the phase factor) for a complex kernel, or the prose should note that for the Gaussian case $b_+$ is complex and the "sgn" notation means the unit-modulus phase.

### 1.5 QSP degree and Laurent polynomial (Algorithm 3, step 13)

The QSP step references Motlagh & Wiebe 2024 [MW24] correctly.  The degree bound
$$d = O(\delta\alpha + \log(1/\varepsilon_\mathrm{QSP})/\log\log(1/\varepsilon_\mathrm{QSP}))$$
is consistent with the Hamiltonian simulation cost from GQSP.  The target Laurent polynomial $P(e^{i\theta}) \approx e^{-i\delta\alpha\sin\theta}$ is correct: since the block encoding gives $B_a/\alpha$ as $\langle 0|U_{B_a}|0\rangle = B_a/\alpha$, and the eigenvalues of the block encoding unitary are $e^{\pm i\theta}$ where $\sin\theta = \lambda/\alpha$ for eigenvalue $\lambda$ of $B_a$, implementing $e^{-i\delta B_a}$ requires the polynomial $e^{-i\delta\alpha\sin\theta}$.

**No correctness issue here**, but this is exactly the kind of reasoning that MUST appear in the prose.


---

## 2. Notation Inconsistencies

### 2.1 `Ctrl-Trott` vs `Ctrl-Trot` (single t)
- Algorithm 2 uses `Ctrl-Trott` (lines 1419, 1421, 1425, 1427, 1460, 1464) -- double t.
- Algorithm 3 uses `Ctrl-Trot` (line 1503) -- single t.
- Algorithm 3 later uses `Ctrl-Trott` (lines 1519, 1523, 1527, 1534) -- double t.

This should be standardized to one spelling throughout.

### 2.2 Product formula notation: $S_p^{(M)}$ vs $S_2^{(M_\pm)}$
- Algorithm 2 uses $S_p^{(M)}(\bar t)$ -- generic order $p$, generic $M$.
- Algorithm 3 uses $S_2^{(M_-)}$ and $S_2^{(M_+)}$ -- specific order 2 (Strang), specific $M_\pm$.
- The Hamiltonian Trotterization subsection (p. 72--79) uses $S_p^{(M)}$ generically and specializes to $p=2$ for Strang.

**Issue**: Algorithm 2 leaves the product formula order generic ($p$) while Algorithm 3 specializes to $p=2$.  Since the Generator splitting subsection (p. 81--84) argues that first-order Lie-Trotter suffices for the coherent-dissipative split (the $O(\delta^2)$ error is already set by the weak-measurement), it is inconsistent for Algorithm 3 to use Strang ($p=2$) while Algorithm 2 uses generic $p$.  The choice should be made explicit and consistent.

### 2.3 Filter function notation
- Algorithm 1 inputs "filter function $f$ (5.3)" and "transition weight $\gamma$".
- Algorithm 2 inputs "filter function $f$, transition weight $\gamma$".
- But in the StatePrep step (Alg. 2, step 2), the amplitudes are $\bar f(\bar t)$ -- the discretized filter from equation (5.48), not $f$ itself.

The relationship $\bar f(\bar t) := \sqrt{t_0} f(\bar t)$ (absorbing the Riemann-sum prefactor) should be stated.  Currently the reader must know to look at equation (5.48) to understand the normalization.

### 2.4 Equation references use the `\eqref{}` macro for TeX equations
The pseudocode references (5.80), (5.3), (5.61), (5.40), (5.43), (5.35) which are all in the earlier subsections.  This is fine for the TeX source, but there is no forward reference to explain the notation in the Algorithm subsection itself -- reinforcing the need for introductory prose.

### 2.5 The `Trotter OFT` label in Algorithm 2
Step 2's section header says "Trotterized OFT (5.80)".  Equation (5.80) defines $\tilde A_a(\bar\omega)$ -- the Trotterized OFT.  This is correct but should be expanded in prose.

### 2.6 $M_\mathcal{A}$ vs $|\mathcal{A}|$
Algorithm 1 uses $M_\mathcal{A}$ for the number of jumps (line 1386).  The Generator splitting subsection (p. 82, eq. 5.82) defines $M_\mathcal{A} := |\mathcal{A}|$.  This notation is used consistently, but never defined within the Algorithm subsection itself.


---

## 3. Missing Prose -- Critical Gaps

This is the most important part of the review.  The Algorithm subsection currently has NO text -- only three algorithm floats.  The following prose elements are essential:

### 3.1 Introductory paragraph: what the algorithm does and how it relates to the theory

The subsection needs an opening paragraph that:
- Summarizes the full algorithm at a high level (continuous-time Markov chain simulation via Lie-Trotter splitting of the Lindbladian, with each step implemented via weak-measurement for the dissipative part and QSP-based Hamiltonian simulation for the coherent part).
- States the main result: the total cost to achieve $\varepsilon$-accuracy in trace distance.
- References the key prior work: CKBG23 for the weak-measurement scheme, CKG23 for the coherent term and KMS detailed balance, MW24 for GQSP.

### 3.2 StatePrep subroutine -- HOW to prepare the time registers

Both Algorithms 2 and 3 invoke `StatePrep` to prepare states like $|f\rangle = \sum \bar f(\bar t)|\bar t\rangle$ and $|b_\pm\rangle$.  The prose MUST address:

**(a) What is the canonical algorithm for this?**
- For the Gaussian filter $f(t) = (2\sigma^2/\pi)^{1/4} e^{-\sigma^2 t^2}$, the amplitudes $\bar f(\bar t) = \sqrt{t_0} f(\bar t)$ form a discretized Gaussian.  Gaussian states can be prepared efficiently.  CKBG23 (p. 24--25) mentions that Gaussian states are "relatively easy to prepare" and cites [MGB22] (McArdle, Gilyén, Berta 2022) and [BSG+22] for quantum-accessible window functions.
- For the $b_-$ kernel (equation 5.40), the amplitudes involve $\sin(-\beta\sigma t) e^{-2t^2} / \cosh(2\pi t/(\beta\sigma))$.  This is NOT a simple Gaussian.  How is this prepared?  Options include: (i) QSVT/QSP on a simpler state, (ii) a classically precomputed circuit synthesis.  The cost of this state preparation directly affects the overall complexity.
- For the $b_+$ kernel (equation 5.43), similarly complex amplitudes.

**(b) What is the cost?**
- CKBG23 states the cost of $\mathrm{Prep}_f$ is $O(\mathrm{poly}(\log N, \log(1/\varepsilon)))$ for Gaussian-like states (see footnote 33 and the discussion around equation (3.9)).
- For the $b_\pm$ kernels in the CoherentStep, this is less clear.  The $b_-$ amplitudes involve a product of functions that are not simple Gaussians.  A generic state-preparation circuit on $r$ qubits costs $O(2^r)$ gates, which would dominate the algorithm.  It is essential to either (i) cite a specific efficient preparation method, or (ii) argue that classical precomputation of an $O(r)$-depth circuit is possible because the amplitudes have special structure (e.g., they are smooth, rapidly decaying, etc.).

**(c) Reference to Chapter 4 (Quantum Circuits):**
Chapter 4 is currently a placeholder (all subsection titles, no content).  The "Basics" and "Weak-measurement techniques" subsections of Ch 4 are the natural home for the circuit-level details.  The Algorithm prose should at minimum include a forward reference: "For the circuit-level implementation of StatePrep, see Section 4.X."

### 3.3 Boltzmann ancilla rotation -- efficient implementation of controlled-$\gamma(\bar\omega)$

Algorithm 2, step 7 applies a rotation $R_Y(2\arcsin\sqrt{1-\gamma(\bar\omega)})$ controlled on the frequency register $|\bar\omega\rangle$.  This is a multi-controlled rotation parameterized by a classical function $\gamma$ evaluated at each grid point $\bar\omega$.  The prose MUST discuss:

**(a) How to implement this efficiently:**
- The naive approach is to decompose the controlled rotation into $N = 2^r$ individually controlled rotations, costing $O(N)$ gates -- exponential in $r$.
- For the Gaussian weight $\gamma_G$: CKBG23 (footnote 33, p. 25) suggests using controlled Hamiltonian simulation for $Y$ rotation: $\sum_{\bar\omega} e^{i(\bar\omega/2\|H\|)Y} \otimes |\bar\omega\rangle\langle\bar\omega|$ followed by QSVT to implement the polynomial $\arcsin\sqrt{1-\gamma(\bar\omega)}$.  This costs $O(\mathrm{poly}(\log N, \log(1/\varepsilon)))$.
- For Metropolis $\gamma_M$: the kink at $\omega = -\beta\sigma^2/2$ means the function is only Lipschitz, requiring higher-degree polynomial approximation.  The thesis itself notes this on p. 53--54 ("one is to use QSP, which has [...] but due to the kink we would need higher degree polynomial [...]").
- For smooth Metropolis $\gamma_M^{(s)}$: being Gevrey-1/2 (Proposition 7), it admits efficient polynomial approximation with QSP/QSVT.

**(b) The cost and its dependence on the choice of $\gamma$:**
This is a key practical bottleneck that directly connects to the Designing transition weights subsection.  The smooth Metropolis was introduced precisely to tame these costs.  The Algorithm prose should close this loop by stating the controlled-rotation cost for each $\gamma$ choice.

**(c) The `[REF to text]` placeholder on line 1439** must be resolved.

### 3.4 QSP for the coherent step -- reference and explanation

Algorithm 3, step 13 invokes "QSP [MW24]" with precomputed rotation angles.  The prose needs:

**(a) Why GQSP (Motlagh & Wiebe 2024) rather than standard QSP/QSVT:**
- The block encoding $U_{B_a}$ is constructed via LCU (Linear Combination of Unitaries -- the state-preparation/controlled-evolution/unpreparation sandwich).
- Standard QSP (Low & Chuang 2017, 2019) requires access to both $U$ and $U^\dagger$ with controlled applications, and the polynomial must satisfy parity constraints.
- GQSP removes these restrictions: it uses general SU(2) rotations as signal processing operators, the polynomial $P$ needs only $|P| \leq 1$ on the unit circle, and there are no parity constraints.  This is precisely why the Laurent polynomial $P(e^{i\theta}) \approx e^{-i\delta\alpha\sin\theta}$ (which has indefinite parity) can be directly implemented.
- Furthermore, GQSP only needs controlled-$U$ (not $U^\dagger$), which simplifies the circuit.

**(b) The degree bound:**
The stated bound $d = O(\delta\alpha + \log(1/\varepsilon_\mathrm{QSP})/\log\log(1/\varepsilon_\mathrm{QSP}))$ is from GQSP's Hamiltonian simulation application (Section V of MW24).  This should be cited more precisely.

**(c) What $\alpha = \|b_-\|_1 \|b_+\|_1$ means physically and numerically:**
- For Gaussian $\gamma_G$: $\|b_+\|_1 = O(1)$ (eq. 5.41), $\|b_-\|_1 = O(e^{\beta^2\sigma^2/8})$ (eq. 5.40), so $\alpha = O(e^{\beta^2\sigma^2/8})$.  With $\sigma = 1/\beta$: $\alpha = O(e^{1/8}) = O(1)$.
- For smooth Metropolis: $\|b_+^{(s,\eta)}\|_1 = O(\log(\beta\|H\|/\varepsilon))$ (eq. 5.43 and surrounding discussion), so $\alpha = O(\log(\beta\|H\|/\varepsilon))$.
- The QSP degree is $d = O(\delta \cdot \alpha)$ in the dominant regime, so this directly affects the per-step cost.  This information is scattered across the Time, time, time and Quadrature errors for $B$ subsections but never consolidated.

**(d) Reference to Chapter 4 (QSP subsection):**
Chapter 4's "Quantum Signal Processing" subsection is empty.  The Algorithm prose should state what QSP framework is used and reference both Ch 4 (once written) and MW24.

### 3.5 Weak-measurement technique -- how and why it works

Algorithm 2 implements the weak-measurement scheme from CKBG23, Theorem III.1 (p. 19, Figure 3).  The prose MUST explain:

**(a) The principle:**
The weak-measurement scheme exploits the quantum Zeno-like effect: for small step size $\delta$, the probability of a "jump" (measuring $q_\delta = |1\rangle$) is $O(\delta)$.  With probability $1 - O(\delta)$, no jump occurs and the system is approximately unperturbed.  The net channel on the system after tracing out the ancillas is $\rho \mapsto (I + \delta\mathcal{L})(\rho) + O(\delta^2)$, which approximates $e^{\delta\mathcal{L}}$ to first order.

**(b) Why the error is $O(\delta^2)$ per step and $O(t_\mathrm{mix} \cdot \delta)$ total:**
This is already discussed in the Generator splitting subsection (p. 81--84), but the Algorithm subsection should recapitulate: since we need $L = O(t_\mathrm{mix}/\delta)$ steps and each has $O(\delta^2)$ error, the total error is $O(t_\mathrm{mix} \cdot \delta)$.  Setting $\delta = \varepsilon / t_\mathrm{mix}$ gives $L = O(t_\mathrm{mix}^2/\varepsilon)$ steps.

**(c) Cost comparison with CKBG23's two methods:**
CKBG23 offers two algorithms: the simple weak-measurement scheme (Theorem III.1, $O(t^2/\varepsilon)$ scaling) and the compressed scheme (Theorem III.2, near-linear in $t$).  The thesis uses the simple one.  This choice should be justified -- likely because the simple scheme has lower constant factors and is more practically relevant.

### 3.6 Error budget -- combining all error sources at the algorithm level

This is perhaps the most important missing piece for a thesis.  The prior subsections individually analyze:
- Quadrature/discretization errors for $\mathcal{L}_\mathrm{diss}$: Table 5.1, equations (5.53), (5.60)
- Quadrature errors for $B$: equations (5.63)--(5.68), (5.75), (5.78)
- Trotterization errors: Propositions 9 and 10
- Generator splitting errors: equations (5.83)--(5.86) (jump-wise), (5.90)--(5.93) (coherent-dissipative)
- Weak-measurement error: $O(\delta^2)$ per step

The Algorithm subsection needs a **consolidated error budget** that:
1. Lists all error sources
2. Shows how to choose the algorithmic parameters ($\delta$, $r$, $M$, etc.) to achieve a target total error $\varepsilon$
3. States the resulting total complexity (gate count, ancilla count, Hamiltonian simulation time)

A natural format would be a **theorem or proposition** stating: "Given parameters [...], Algorithm 1 prepares a state $\varepsilon$-close to $\rho_\beta$ in trace distance using [total resources]."  This would be the capstone result of Chapter 5.

### 3.7 Ancilla budget

The three algorithms together use the following ancilla registers:
- $\Omega$: $r$ qubits (frequency/time register for the OFT in the dissipative step)
- $q_\gamma$: 1 qubit (Boltzmann)
- $q_\delta$: 1 qubit (weak measurement)
- $T_-$: $r_-$ qubits (outer time register for the coherent step)
- $T_+$: $r_+$ qubits (inner time register for the coherent step)
- $q_\mathrm{QSP}$: 1 qubit (QSP signal)
- Plus any ancillas needed for the Trotterized Hamiltonian simulation and state preparation

The prose should consolidate the total ancilla count: $r + r_- + r_+ + 3 + O(K)$ where $K$ is the number of terms in $H = \sum_k H_k$ (for Trotter).  The values of $r, r_\pm$ are determined by the quadrature error analysis and should be stated explicitly:
- $r = O(\log(|\mathcal{A}|(\|H\|+1)(\sigma+1/\sigma)(\beta+1) + \frac{1}{s}\log^2(1/\varepsilon)))$ for smooth Metropolis (eq. 5.60)
- $r_- = O(\log(\beta\sigma\|H\|\,|\mathcal{A}|\,\|b_+\|_1\log(1/\varepsilon)/\varepsilon))$ (eq. 5.68)
- $r_+ = O(\log(\|H\|\|\sum_a A^{a\dagger}A^a\|\,|\mathcal{A}|\log(1/\varepsilon)/(\sigma(1+s)\varepsilon)))$ (eq. 5.78)

### 3.8 The role of the coherent step in the algorithm

The CKG construction is distinguished from CKBG by the presence of the coherent term $B = -iG/2$ (where $G$ from equation 5.29 restores exact KMS detailed balance).  The Algorithm subsection should explicitly state:
- Without the CoherentStep (i.e., setting $B_a = 0$), the algorithm reduces to the CKBG Gibbs sampler with only approximate (GNS) detailed balance.
- The CoherentStep is what upgrades the algorithm from approximate to exact KMS detailed balance (modulo discretization/Trotter errors).
- This is the key conceptual contribution of CKG23 over CKBG23.


---

## 4. Missing Citations and Cross-References

### 4.1 `[REF to text]` placeholder (line 1439)
The Boltzmann rotation step in Algorithm 2 has `\Comment{see [REF to text]}`.  This needs to reference either:
- A prose paragraph in the Algorithm subsection explaining the controlled rotation
- The Designing transition weights subsection where $\gamma$ is defined
- CKBG23, Section III.B (explicit block-encodings), specifically the discussion of the controlled filter $W$ on p. 25

### 4.2 Missing citation for the weak-measurement scheme
Algorithm 2's caption says "weak-measurement simulation of $e^{\delta\mathcal{L}_{a,\mathrm{diss}}}$" but does not cite CKBG23, Theorem III.1.  This should be: "...weak-measurement simulation of $e^{\delta\mathcal{L}_{a,\mathrm{diss}}}$ \cite{chen2023quantum}".

### 4.3 Missing citation for block encoding via LCU
Algorithm 3 constructs a block encoding of $B_a/\alpha$ but does not cite the LCU framework.  The LCU technique is from Childs, Kothari, and Somma (2012) or Berry et al. (2015), and should be cited in the prose.  If Ch 4's "Linear Combination of Unitaries" subsection is written, it should be cross-referenced.

### 4.4 Missing cross-references to Chapter 4
Algorithm 2 uses: StatePrep, controlled Hamiltonian simulation, QFT, weak measurement.
Algorithm 3 uses: StatePrep, controlled Hamiltonian simulation, QSP.
All of these are listed as subsections of Chapter 4 (which is currently empty).  The prose should include forward references to Ch 4 for each of these primitives.

### 4.5 Missing cross-references to earlier Ch 5 subsections
The following should be explicitly cross-referenced in the prose:
- The OFT: equation (5.1), (5.47)--(5.48)
- The discretized Lindbladian: equation (5.47)
- The coherent term: equations (5.29), (5.39), (5.61)
- The generator splitting: equations (5.82)--(5.93)
- The quadrature error analysis: Table 5.1, equations (5.53), (5.60)
- The Trotter error analysis: Propositions 9 and 10

### 4.6 Missing citations for asymptotic Hamiltonian simulation costs
The total Hamiltonian simulation time is a key figure of merit.  The `src/simulation_time.jl` file in the codebase computes this numerically.  The prose should state the analytical scaling and cite the relevant complexity results (CKG23 Theorem I.2 for the overall scaling).

### 4.7 Ding, Li, Lin (DLL25) comparison
The Ding et al. papers (2024, 2025) propose alternative KMS detailed-balanced Lindbladians.  The thesis memory notes a planned "Ding & Chen comparison" chapter.  At minimum, the Algorithm subsection prose should note that the algorithmic structure (OFT + weak measurement + coherent correction) is common to both the CKG and DLL frameworks, with the differences being in the choice of filter functions and transition weights.


---

## 5. Structural Recommendations

### 5.1 Suggested prose structure

The Algorithm subsection should be organized as follows (suggested section headers in italics):

1. *Opening paragraph* (1/2 page): High-level summary of the full algorithm and its cost. State the main complexity result as an informal theorem.

2. *Dissipative step: weak-measurement simulation* (1 page): Explain the principle (Section 3.5 above), the connection to CKBG23, the error scaling, and the cost.  Discuss the StatePrep subroutine (Section 3.2) and the Boltzmann rotation (Section 3.3).

3. *Coherent step: block encoding and QSP* (1 page): Explain the LCU construction of $U_{B_a}$, why GQSP is used (Section 3.4), the subnormalization $\alpha$ and its impact on cost, and the StatePrep costs for $b_\pm$.

4. *Error budget and total complexity* (1--1.5 pages): Consolidate all error sources (Section 3.6), state the parameter choices, and give the total resource count (Section 3.7).

5. *Remarks on implementation variants* (1/2 page): Sequential vs random jump selection, Strang vs Lie-Trotter for coherent-dissipative splitting, comparison with CKBG23's compressed algorithm.

The three algorithm floats should be interspersed in this prose (Algorithm 1 after the opening, Algorithm 2 after the dissipative discussion, Algorithm 3 after the coherent discussion).

### 5.2 Consider adding a complexity table

A summary table listing the per-step and total costs for each choice of transition weight ($\gamma_G$, $\gamma_M$, $\gamma_M^{(s)}$) would be very valuable.  The columns could be: transition weight | estimating qubits $r$ | Hamiltonian sim time per step | QSP degree | total steps | total cost.

### 5.3 The comment on line 1360 (TeX source)

There is a commented-out remark at line 1360:
```
% Use Trotter splitting also for coherent and dissipative parts. t_mix
% scaling becomes poly with this compared to linear for Chen's, but
% doesn't need QSVT...
```
This is an important implementation consideration (a NISQ-friendly variant that avoids QSP entirely by Trotterizing the coherent evolution directly).  It should be either developed into a remark in the prose or removed.


---

## 6. Summary of Action Items

| Priority | Item | Section |
|----------|------|---------|
| **P0** | Write introductory prose paragraph | 3.1 |
| **P0** | Write error budget / total complexity result | 3.6 |
| **P0** | Explain weak-measurement scheme in prose | 3.5 |
| **P0** | Resolve `[REF to text]` placeholder (line 1439) | 4.1 |
| **P1** | Explain StatePrep subroutine and its cost | 3.2 |
| **P1** | Explain controlled-$\gamma(\bar\omega)$ rotation and its cost | 3.3 |
| **P1** | Explain QSP choice (GQSP) and degree bound | 3.4 |
| **P1** | Clarify unconditional-vs-controlled uncomputation (Alg 2, steps 9--14) | 1.2 |
| **P1** | Consolidate ancilla budget | 3.7 |
| **P1** | Explain role of coherent step (CKG vs CKBG) | 3.8 |
| **P2** | Standardize `Ctrl-Trott` vs `Ctrl-Trot` spelling | 2.1 |
| **P2** | Fix Trotter order inconsistency ($S_p^{(M)}$ vs $S_2^{(M)}$) | 2.2 |
| **P2** | Handle complex $\mathrm{sgn}$ for $b_+$ state prep | 1.4 |
| **P2** | Clarify $\delta$ definition (per-jump vs per-sweep) | 1.3 |
| **P2** | Add citation to CKBG23 in Alg 2 caption | 4.2 |
| **P2** | Add LCU citation in Alg 3 prose | 4.3 |
| **P2** | Add cross-references to Ch 4 subsections | 4.4 |
| **P3** | Add complexity summary table | 5.2 |
| **P3** | Develop or remove NISQ-variant comment (line 1360) | 5.3 |
| **P3** | Note DLL25 comparison where relevant | 4.7 |
