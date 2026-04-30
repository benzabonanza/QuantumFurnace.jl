# State Preparation Methods for the CKG Algorithm
<!-- Survey of quantum state preparation techniques for the time/frequency register states in the CKG Gibbs sampler -->

---

## The problem

Three state preparations appear in the CKG algorithm, each called many times per Lindbladian step:

| State | Function | Properties | Register size |
|-------|----------|------------|---------------|
| $\|f\rangle$ (dissipative OFT) | $\bar f(\bar t) = \sqrt{t_0}\,f(\bar t)$, Gaussian filter (5.3) | Real, **positive**, $C^\infty$, Gaussian decay | $r \approx 8$--$15$ qubits |
| $\|b_-\rangle$ (coherent outer) | $b_-(\bar t)$ from (5.40): $1/\cosh$ convolved with $\sin\cdot e^{-t^2}$ | Real, **signed** (changes sign), smooth, rapidly decaying, support $\mathcal{O}(\beta\sigma)$ | $r_- \approx 5$--$8$ qubits |
| $\|b_+\rangle$ (coherent inner) | $b_+(\bar t')$ from (5.41)/(5.43) | **Complex-valued** (Gaussian envelope $\times$ complex phase); Metropolis variant has regularised $1/t$ | $r_+ \approx 10$--$15$ qubits |

In the LCU framework, the prepared state encodes:
$$|b\rangle = \sum_j \sqrt{\frac{\Delta t\,|b(t_j)|}{\alpha}}\;\operatorname{sgn}(b(t_j))\,|j\rangle$$
where the absolute values go into the amplitudes (ensuring $\langle b|b\rangle = 1$) and the signs/phases go into the quantum state phases.

### Requirements
- **Efficiency**: state prep is called $d$ times per CoherentStep ($d$ = QSP degree, typically 1--2) and once per DissipativeStep, so it should be cheap relative to the Hamiltonian simulation that dominates the circuit.
- **Signed/complex amplitudes**: $b_-$ is real but changes sign; $b_+$ is complex. The preparation must encode these phases.
- **Precision**: the per-step channel error from quadrature is $\delta \cdot \|B_a - \bar{B}_a\|$, and we want this $\leq \delta^2/10$ (an order below the weak-measurement baseline $\delta^2$). This requires $\|B_a - \bar{B}_a\| \leq \delta/10$, which sets the quadrature target $\varepsilon = \delta/10$.
- **Ancilla budget**: ideally $\leq 4$--$6$ ancillae beyond the register itself.

### Per-jump register sizes

**Important**: the thesis quadrature bounds (5.63)--(5.78) are for the *full* coherent term $B = \sum_a B_a$, but the CoherentStep block-encodes a *single* $B_a$. Two factors reduce the required register sizes:
1. The sum $\sum_a$ disappears (no $|\mathcal{A}|$ factor).
2. The per-jump operator norm $\|A_a\|^2 = 1/M_\mathcal{A}$ (from $\|\sum A^\dagger A\| \leq 1$) gives an additional $1/M_\mathcal{A}$ suppression.

The per-jump formulas (replacing $|\mathcal{A}| \to 1/M_\mathcal{A}$ in eqs 5.68, 5.75):

$$r_-^{\text{per-jump}} = \Big\lceil\log_2\!\Big(\frac{\beta\sigma \cdot \|H\| \cdot \|b_+\|_1 \cdot \ln(1/\varepsilon)}{M_\mathcal{A} \cdot \varepsilon}\Big)\Big\rceil, \qquad r_+^{\text{per-jump},(G)} = \Big\lceil\log_2\!\Big(\frac{\beta\|H\| \cdot \ln(1/\varepsilon)}{M_\mathcal{A} \cdot \varepsilon}\Big)\Big\rceil$$

With $\beta\sigma = 1$, $\|H\| = 0.5$, $\|b_+\|_1 \approx 1/(2\pi)$, $M_\mathcal{A} = 33$ (single-site Paulis, $n = 11$), $\varepsilon = \delta/10$:

| | | $r_-$ | $r_+$ ($\beta{=}5$) | $r_+$ ($\beta{=}10$) | $r_+$ ($\beta{=}20$) |
|---|---|---|---|---|---|
| $\delta = 0.01$ | $\varepsilon = 10^{-3}$ | **5** | **10** | **11** | **11** |
| $\delta = 0.001$ | $\varepsilon = 10^{-4}$ | **8** | **13** | **14** | **15** |

These are **~10 qubits smaller** than the full-$B$ estimates due to the $M_\mathcal{A}^2$ reduction inside the log. Total qubit budget for CoherentStep ($\delta = 0.01$, $\beta = 10$): $n + r_- + r_+ + 1 = 11 + 5 + 11 + 1 = 28$ qubits.

---

## Methods

### 1. Grover--Rudolph (2002)
**Ref**: [arXiv:quant-ph/0208112](https://arxiv.org/abs/quant-ph/0208112)

- Prepares log-concave distributions (e.g. Gaussians) in $\mathcal{O}(r)$ controlled rotations.
- Each qubit of the register is prepared by a rotation conditioned on all previously prepared qubits, using the conditional probability $\Pr(x_k = 1 \mid x_1, \ldots, x_{k-1})$.
- **Handles only non-negative real amplitudes** --- no sign/phase encoding.
- **Best fit for**: the Gaussian filter $\bar f$ (positive, log-concave). Essentially optimal: $\mathcal{O}(r)$ gates.
- **Not suitable for**: $b_-$ (sign changes) or $b_+$ (complex).
- No fault-tolerant resource estimates in the original paper.

### 2. McArdle, Gily\'en, Berta (2022)
**Ref**: [arXiv:2210.14892](https://arxiv.org/abs/2210.14892) --- *"Quantum state preparation without coherent arithmetic"*

- Uses Quantum Eigenvalue Transformation (QET, a variant of QSP) to convert a block-encoding of $\sin(\pi x / 2^r)$ into the target function via polynomial approximation.
- **Gate count**: $\mathcal{O}(r \cdot d / F)$ where $d$ is the polynomial degree and $F$ is the $\ell^2$ filling fraction ($\|c\|_2 / (\|c\|_\infty \sqrt{N})$).
- **Ancillae**: 3--4 qubits.
- **Signed/complex amplitudes**: handled via LCU decomposition into even/odd polynomial parts (real case) or real/imaginary parts (complex case). Doubles the cost.
- **Smooth functions**: for $C^\infty$ or Gevrey functions, polynomial degree $d = \mathcal{O}(\log(1/\varepsilon))$ suffices (exponential convergence), giving total cost $\mathcal{O}(r \cdot \log(1/\varepsilon))$.
- **This is the method cited in CKBG23** [Chen et al. 2023, footnote] for preparing the Gaussian window states.
- **Best fit for**: all three of our functions. The smooth Gaussian cases ($\bar f$, $b_-$, Gaussian $b_+$) converge exponentially. The smooth Metropolis $b_+^{(s,\eta)}$ is Gevrey-$1/2$ (Proposition 7), so also converges exponentially with slightly larger constants.
- **Practical estimate**: for $r = 10$, $\varepsilon = 10^{-4}$: roughly $d \approx 15$--$25$ polynomial degree, $\sim 200$--$400$ T-gates per call.

### 3. O'Brien and Sunderhauf --- Piecewise QSVT (2025)
**Ref**: [arXiv:2409.07332](https://arxiv.org/abs/2409.07332) --- published in *Quantum* (2025)

- Approximates the target amplitude function by **piecewise polynomials** on $S$ segments, handling discontinuities and near-singularities gracefully.
- **Gate count**: $\mathcal{O}(d \cdot \max(r, \log(1/\varepsilon)))$ Toffoli gates for $S = \mathcal{O}(r)$ segments.
- **Ancillae**: $2r + \mathcal{O}(d \cdot \log(1/\varepsilon))$.
- Achieves $\sim 50\times$ fewer Toffolis than Kaiser-window methods for QPE applications.
- **Best fit for**: the smooth Metropolis $b_+^{(s,\eta)}$ (eq 5.43), which has a regularised $1/t$ singularity near $t' = 0$. A global polynomial would need high degree to approximate the singularity region, but a piecewise approximation can place segment boundaries at $t' = \pm\eta$ and use low-degree polynomials on each piece.
- Also a strong option for $b_-$ if the sign changes cause the global polynomial degree to be high.

### 4. Sanders, Low, Scherer, Berry (2019)
**Ref**: [arXiv:1807.03206](https://arxiv.org/abs/1807.03206) --- PRL

- Black-box method: treats $N$ amplitudes as a lookup table via QROM (quantum read-only memory).
- **Gate count**: $\mathcal{O}(N)$ T-gates. Does **not** exploit smoothness.
- For $r = 10$ ($N = 1024$): manageable. For $r = 15$ ($N = 32768$): may dominate the circuit.
- Handles arbitrary complex amplitudes natively.
- **Best fit for**: fallback option when the target function is irregular or when simplicity is preferred over optimality.

### 5. Bausch (2022)
**Ref**: [arXiv:2009.10709](https://arxiv.org/abs/2009.10709) --- published in *Quantum*

- Improves black-box state prep via structured amplitude amplification.
- Speedup depends on the "filling fraction" (ratio of $\ell^2$ to $\ell^\infty$ norm $\times \sqrt{N}$).
- For concentrated distributions (like narrow Gaussians on a wide grid), the filling fraction is small and the speedup is limited.
- **Less relevant** for our case since McArdle et al. and piecewise QSVT already exploit smoothness more directly.

### 6. Gosset, Kothari, Wu (2024)
**Ref**: [arXiv:2411.04790](https://arxiv.org/abs/2411.04790)

- Proves the **optimal T-count** for arbitrary (unstructured) state preparation: $\Theta(\sqrt{2^n \cdot \log(1/\varepsilon)} + \log(1/\varepsilon))$.
- For $n = 10$, $\varepsilon = 10^{-4}$: $\sim 100$ T-gates lower bound.
- Structured/smooth states can beat this via methods 2--3 above.

---

## Comparison for our functions

| Function | Best method | Exploit smoothness? | Handles signs/phases? | Gate cost | Notes |
|----------|------------|--------------------|-----------------------|-----------|-------|
| $\bar f$ (Gaussian filter) | Grover--Rudolph | Yes (log-concave) | N/A (positive) | $\mathcal{O}(r)$ | Simplest option; essentially free |
| $b_-$ (outer kernel, real signed) | McArdle et al. | Yes ($C^\infty$, Gevrey) | Yes (even/odd LCU) | $\mathcal{O}(r \cdot \log(1/\varepsilon))$ | Sign handled by decomposing into even + odd parts |
| $b_+^{(G)}$ (Gaussian inner) | McArdle et al. | Yes (Gaussian) | Yes (Re/Im LCU) | $\mathcal{O}(r \cdot \log(1/\varepsilon))$ | Complex phase handled by Re + iIm decomposition |
| $b_+^{(s,\eta)}$ (Metropolis inner) | Piecewise QSVT | Yes (piecewise smooth) | Yes | $\mathcal{O}(d \cdot r)$ | Segment boundaries at $\pm\eta$ handle the regularised singularity |
| Any (fallback) | Sanders et al. QROM | No | Yes | $\mathcal{O}(N)$ | Simple, generic, but doesn't scale with smoothness |

---

## Recommendations

### Qubit-efficient methods (zero or minimal ancillae)

For a partial fault-tolerance implementation where qubit count is the primary constraint, the QSP/QSVT-based methods (McArdle et al., piecewise QSVT) have a significant drawback: they each need their own signal/ancilla qubits (3--4 for McArdle et al., $2r + \mathcal{O}(d\log(1/\varepsilon))$ for piecewise QSVT).  These ancillae **cannot be shared** with the QSP ancilla $q_\mathrm{QSP}$ used for implementing $e^{-\mathrm{i}\delta B_a}$, since the state preparation occurs *inside* the block encoding circuit where those registers are in use.

The following zero-ancilla or low-ancilla approaches are better suited:

#### (a) Grover--Rudolph for positive log-concave functions (zero ancillae)

Applies directly to the Gaussian filter $\bar f$ and to the amplitude part $|b_\pm|$ when it is log-concave.  The algorithm uses $r$ layers of controlled $R_Y$ rotations, each conditioned on previously prepared qubits.  For a Gaussian, the conditional probabilities have closed-form expressions.  Total: $\mathcal{O}(r)$ multi-controlled rotations $\to$ $\mathcal{O}(r^2)$ elementary gates with **zero ancillae**.

#### (b) Uniformly controlled rotations (zero ancillae, generic)

The most general zero-ancilla approach.  Any state $|0\rangle^{\otimes r} \to \sum_j c_j |j\rangle$ can be decomposed as:
1. **Amplitude preparation**: a sequence of uniformly controlled $R_Y(\theta_j)$ rotations, producing $\sum_j \sqrt{|c_j|}\,|j\rangle$.
2. **Phase correction**: a sequence of uniformly controlled $R_Z(\phi_j)$ rotations, adding the phases $\arg(c_j)$.

Each uniformly controlled rotation block uses $\mathcal{O}(2^r)$ CNOTs + $\mathcal{O}(2^r)$ single-qubit rotations with **zero ancillae** (via Gray-code ordering, see Shende--Bullock--Markov 2006, Möttönen et al. 2004).  For $r = 10$: $\sim 2000$ elementary gates total.  This is the most qubit-efficient generic method and is entirely dominated by the Hamiltonian simulation cost.

#### (c) Exploiting function structure to reduce gate count

For our specific functions, significant savings are possible:

- **Gaussian filter $\bar f$**: positive, log-concave $\Rightarrow$ Grover--Rudolph.  $\mathcal{O}(r^2)$ gates, zero ancillae.

- **Gaussian $b_+$**: the modulus $|b_+(t)| \propto \exp(-4\beta\omega_\gamma t^2)$ is a Gaussian (log-concave!), preparable by Grover--Rudolph in $\mathcal{O}(r_+^2)$ gates.  The phase $\arg(b_+(t)) = -2\beta\omega_\gamma t$ is **linear in $t$**, which corresponds to a phase gradient $\sum_j e^{-\mathrm{i}\omega j} |j\rangle$ --- implementable with $\mathcal{O}(r_+)$ single-qubit $R_Z$ gates (one per qubit of the register, with angle $\omega \cdot 2^k$ for the $k$-th qubit).  **Total: $\mathcal{O}(r_+^2)$ gates, zero ancillae.**

- **Outer kernel $b_-$**: real-valued with a small number of sign changes ($K \sim 2$--$4$ on the grid).  Prepare $|b_-|$ as amplitudes (Grover--Rudolph if log-concave, or uniformly controlled $R_Y$ otherwise), then apply signs via a diagonal unitary --- a sequence of $\mathcal{O}(K)$ multi-controlled $Z$ gates, each flipping the sign in a contiguous region.  **Total: $\mathcal{O}(2^{r_-})$ gates worst case, $\mathcal{O}(r_-^2 + K \cdot r_-)$ if $|b_-|$ is log-concave, zero ancillae.**

- **Smooth Metropolis $b_+^{(s,\eta)}$**: the modulus is *not* log-concave (has a bump from the $\eta$-regularisation), so Grover--Rudolph doesn't apply directly.  Use the generic uniformly-controlled-rotation decomposition: $\mathcal{O}(2^{r_+})$ gates, zero ancillae.  For $r_+ = 10$: $\sim 2000$ gates, still small.

#### (d) Walsh--Hadamard truncation for smooth functions

The baseline $\mathcal{O}(2^r)$ gate count from (b) can be reduced for **smooth** functions via Walsh--Hadamard truncation.  A uniformly controlled $R_Y$ with angles $\theta_j$ ($j = 0,\ldots,2^r{-}1$) decomposes in the Walsh--Hadamard basis:
$$\theta_j = \sum_{S \subseteq \{1,\ldots,r\}} c_S \cdot (-1)^{\langle S, j\rangle}$$
Each non-zero Walsh coefficient $c_S$ of order $|S| = k$ requires $k$ CNOTs + 1 rotation in the circuit.  The full expansion has $2^r$ terms (all subsets $S$).  If the target function is smooth, high-order coefficients are suppressed and can be truncated.

**Decay rates by function regularity:**
- **$C^\infty$ / Gevrey** (Gaussian $\bar f$, Gaussian $b_+$, kernel $b_-$): Walsh coefficients decay super-polynomially with order.  Truncating to order $\leq p$ retains $\sum_{k=0}^{p}\binom{r}{k}$ terms --- e.g.\ for $p=3$: 93 terms at $r=10$ vs.\ 1024 full ($11\times$ reduction), 576 at $r=15$ vs.\ 32768 ($57\times$).
- **Piecewise smooth / discontinuous** (Metropolis $b_+^{(s,\eta)}$): the indicator $\mathbb{I}(|t'| \leq \eta)$ in (5.43) creates jump discontinuities at $t' = \pm\eta$, causing Walsh coefficients to decay only as $\mathcal{O}(1/|S|)$ --- same as Fourier coefficients of a step function.  Truncation gives **little benefit**; essentially the full $\mathcal{O}(2^r)$ decomposition is needed.

**Summary of Walsh truncation applicability:**

| Function | Regularity | Walsh decay | Truncation to order $p{=}3$ ($r{=}12$) | Benefit |
|----------|-----------|------------|--------------------------------------|---------|
| $\bar f$ (Gaussian) | $C^\infty$ | super-polynomial | 299 / 4096 terms | $\sim 14\times$ |
| $b_-$ (outer kernel) | $C^\infty$ | super-polynomial | 299 / 4096 terms | $\sim 14\times$ |
| Gaussian $b_+$ | $C^\infty$ | super-polynomial | 299 / 4096 terms | $\sim 14\times$ |
| Metropolis $b_+^{(s,\eta)}$ | piecewise, $C^0$ at $\pm\eta$ | $\mathcal{O}(1/|S|)$ | needs most terms | **negligible** |

**Piecewise workaround for Metropolis $b_+^{(s,\eta)}$:**  One could split the register into three regions ($[-T_+/2, -\eta]$, $[-\eta, \eta]$, $[\eta, T_+/2]$) using a comparator circuit ($\mathcal{O}(r)$ gates), and apply Walsh-truncated UCRs within each region where $b_+^{(s,\eta)}$ IS smooth.  This partially recovers the $C^\infty$ decay per piece, at the cost of the comparator overhead and $3\times$ the per-piece Walsh terms.  For $r_+ \leq 15$, this added complexity is probably not worth it --- the full $\mathcal{O}(2^{r_+})$ decomposition is simpler and still subdominant to Hamiltonian simulation.

### The crossover: zero-ancilla vs McArdle et al.

The choice between zero-ancilla methods and McArdle et al. depends on the register size $r$:

| $r$ | $N = 2^r$ | Zero-ancilla cost | McArdle et al. cost | Winner |
|---|---|---|---|---|
| $\leq 8$ | $\leq 256$ | $\leq 512$ gates | ~200 gates + 3--4 ancillae | **Zero-ancilla** (saves qubits, comparable gate count) |
| 9--12 | 512--4K | 1K--8K gates | ~300 gates + 3--4 ancillae | **Either** (zero-ancilla still tolerable, McArdle cheaper in gates) |
| 13--15 | 8K--32K | 16K--64K gates | ~400 gates + 3--4 ancillae | **McArdle et al.** (100$\times$ gate reduction worth the ancillae) |
| $\geq 16$ | $\geq 64$K | impractical | ~500 gates + 3--4 ancillae | **McArdle et al.** (only viable option) |

**Note on Grover--Rudolph**: the conceptual description gives $\mathcal{O}(r)$ "layers," but each layer is a uniformly controlled rotation with $2^{k-1}$ distinct angles at layer $k$, giving $\mathcal{O}(2^r)$ total gates in the generic case --- the same as the brute-force decomposition.  The $\mathcal{O}(r^2)$ scaling only holds when the rotation angles can be computed from a closed-form conditional probability (true for pure Gaussians, not necessarily for $b_-$ or Metropolis $b_+$).

### Practical recommendation for our parameter regime

From the per-jump register size estimates:
- $r_- \approx 5$--$8$: always in the zero-ancilla regime. Use Grover--Rudolph (if $|b_-|$ is log-concave) or generic uniformly controlled rotations. At most $\sim 256$ gates.
- $r_+ \approx 10$--$15$: **at the crossover**. For $\delta = 0.01$ ($r_+ \approx 10$--$11$): zero-ancilla is fine ($\sim 2000$ gates). For $\delta = 0.001$ ($r_+ \approx 13$--$15$): McArdle et al. becomes advantageous.

**For the thesis**: present both options. State that for $r \lesssim 12$ qubits, generic zero-ancilla decompositions (uniformly controlled rotations, $\mathcal{O}(2^r)$ gates) keep the state preparation cost subdominant to Hamiltonian simulation without requiring additional ancillae.  For $r \gtrsim 13$--$15$ qubits, the McArdle--Gily\'en--Berta method \cite{mcArdle2022quantum} achieves $\mathcal{O}(r \cdot \log(1/\varepsilon))$ gates at the cost of 3--4 ancilla qubits, and becomes the preferred choice as it avoids exponential gate scaling.

**For the code / numerical simulations**: the per-jump register sizes ($r_- \leq 8$, $r_+ \leq 15$) stay within the zero-ancilla regime for all practical $\delta$ values.  No additional state-prep ancillae are needed.

### Asymptotically optimal methods (for fault-tolerant regime)

If ancilla qubits are not constrained and T-gate count is the metric:

#### McArdle--Gily\'en--Berta (2022)

The natural choice, already cited in CKBG23.  Uses QET to convert a block-encoding of $\sin(\pi x/2^r)$ into the target function.  For smooth functions: $d = \mathcal{O}(\log(1/\varepsilon))$ polynomial degree, total cost $\mathcal{O}(r \cdot \log(1/\varepsilon))$ gates + 3--4 ancillae.  Handles complex amplitudes via Re+Im LCU at $2\times$ cost.  For $r = 10$, $\varepsilon = 10^{-4}$: $\sim 200$--$400$ T-gates per call.

#### Piecewise QSVT (O'Brien--Sunderhauf 2025)

Best for the smooth Metropolis $b_+^{(s,\eta)}$: segment boundaries at $\pm\eta$ handle the regularised singularity efficiently.  But ancilla cost ($2r + \mathcal{O}(d\log(1/\varepsilon))$) makes it impractical for qubit-constrained settings.

#### QROM (Sanders et al. 2019)

Fallback for generic/irregular functions.  $\mathcal{O}(N)$ T-gates, $\mathcal{O}(\log N)$ ancillae.  Acceptable for $r \leq 12$.

---

## References

1. Grover and Rudolph (2002). *Creating superpositions that correspond to efficiently integrable probability distributions.* [arXiv:quant-ph/0208112](https://arxiv.org/abs/quant-ph/0208112)
2. McArdle, Gily\'en, Berta (2022). *Quantum state preparation without coherent arithmetic.* [arXiv:2210.14892](https://arxiv.org/abs/2210.14892)
3. O'Brien and Sunderhauf (2025). *Piecewise QSVT state preparation.* [arXiv:2409.07332](https://arxiv.org/abs/2409.07332)
4. Sanders, Low, Scherer, Berry (2019). *Black-box quantum state preparation without arithmetic.* [arXiv:1807.03206](https://arxiv.org/abs/1807.03206), PRL.
5. Bausch (2022). *Fast black-box quantum state preparation.* [arXiv:2009.10709](https://arxiv.org/abs/2009.10709), Quantum.
6. Gosset, Kothari, Wu (2024). *Quantum state preparation with optimal T-count.* [arXiv:2411.04790](https://arxiv.org/abs/2411.04790)
7. Low, Kliuchnikov, Schaeffer (2018). *Trading T-gates for dirty qubits.* [arXiv:1812.00954](https://arxiv.org/abs/1812.00954)
