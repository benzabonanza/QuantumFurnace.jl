# Spectral Gap Estimation via Quantum Trajectory Simulation of CKG Lindbladians

## 1. Problem Statement

**Setup.** Consider the CKG Lindbladian $\mathcal{L}_\beta$ from Chen et al. (2023), Eq. (1.8):

$$\mathcal{L}_\beta[\rho] = -i[B, \rho] + \sum_{a \in \mathcal{A}} \int_{-\infty}^{\infty} \gamma(\omega) \left( \hat{A}^a(\omega) \rho \hat{A}^a(\omega)^\dagger - \frac{1}{2}\{\hat{A}^a(\omega)^\dagger \hat{A}^a(\omega), \rho\} \right) d\omega$$

with jump operators $A^a$ being single-site Paulis $\{X_i, Y_i, Z_i\}$ (self-adjoint, so $\{A^a\} = \{A^{a\dagger}\}$), the operator Fourier transform $\hat{A}^a(\omega)$ as in Eq. (1.9), a Gaussian or Metropolis filter $\gamma(\omega)$ satisfying KMS detailed balance $\gamma(\omega)/\gamma(-\omega) = e^{-\beta\omega}$, and the Hamiltonian being a Heisenberg chain $H = \sum_i \vec{S}_i \cdot \vec{S}_{i+1}$ (with or without external field).

For $n \leq 6$ qubits, you construct the full vectorized Lindbladian $\mathbb{L} \in \mathbb{C}^{4^n \times 4^n}$ and compute $\lambda_{\text{gap}} := -\text{Re}(\lambda_2(\mathbb{L}))$ exactly, where $\lambda_1 = 0$ is the steady-state eigenvalue.

**Goal.** For $n > 8$, estimate $\lambda_{\text{gap}}$ from quantum trajectory simulations by fitting the decay of an observable:

$$\langle O \rangle(t) - \langle O \rangle_{\text{ss}} \sim c \, e^{-\lambda_2 t}$$

where $\langle O \rangle_{\text{ss}} = \text{Tr}[O \rho_\beta]$ is the steady-state expectation value.

**Observed issue.** A persistent, systematic error in the estimated spectral gap compared to exact diagonalization results at small $n$.

---

## 2. The Simulation Chain: What Is Actually Being Simulated

Before cataloguing errors, it is essential to precisely identify what the trajectory simulator computes versus what the exact Lindbladian computes. There are **four distinct dynamical objects** in play:

**(i) The ideal CKG Lindbladian $\mathcal{L}_\beta$.** This is the continuous-time, continuous-frequency object from Chen et al. (2023), Eq. (1.8). Its spectral gap $\lambda_{\text{gap}}(\mathcal{L}_\beta)$ is the target quantity. In principle, $\rho_\beta$ is only an *approximate* fixed point of $\mathcal{L}_\beta$ (the deviation is bounded by the approximate detailed balance analysis, Proposition E.3 of Chen et al. (2023), which gives $\|\rho_{\text{fix}}(\mathcal{L}_\beta) - \rho_\beta\|_1 \leq 14\epsilon / \lambda_{\text{gap}}(H)$ where $\epsilon$ measures the violation of exact detailed balance).

**(ii) The discretized Lindbladian $\mathcal{L}_\beta^{\text{disc}}$.** In practice, you discretize the frequency integral and the time-domain filter $f(t)$. This introduces discretization error bounded by Chen et al. (2023), Appendix C, Lemma C.1. The vectorized $\mathbb{L}$ you construct for $n \leq 6$ is this object.

**(iii) The $\delta$-step weak measurement channel $\Phi_\delta$.** Each step of the trajectory simulator implements a CPTP map $\Phi_\delta$ that approximates $e^{\delta \mathcal{L}_\beta^{\text{disc}}}$. By Theorem III.1 of Chen et al. (2023), $\|(I + \delta \mathcal{L}) - \Phi_\delta\|_\diamond = O(\delta^2)$, and since $\|(I + \delta \mathcal{L}) - e^{\delta \mathcal{L}}\|_\diamond = O(\delta^2)$ as well, we get $\|e^{\delta \mathcal{L}} - \Phi_\delta\|_\diamond = O(\delta^2)$.

Moreover, by Corollary III.1 of Chen et al. (2023), each step randomly selects a single jump $A^a$ with probability $p(a)$ rather than applying all jumps simultaneously, giving the *randomized* channel $\tilde{\Phi}_\delta$ satisfying $\|e^{\delta \mathcal{L}} - \mathbb{E}_a[\Phi_\delta^{(a)}]\|_\diamond = O(\delta^2)$.

**(iv) The quantum trajectory (Monte Carlo wavefunction) unraveling.** Instead of propagating the density matrix under $\Phi_\delta$, each trajectory propagates a pure state $|\psi\rangle$ through the Kraus decomposition of $\Phi_\delta^{(a)}$. This is a stochastic process on pure states. The density matrix is recovered only in expectation: $\mathbb{E}[|\psi(t)\rangle\langle\psi(t)|] = \Phi_\delta^{r}[\rho_0]$ where $r = t/\delta$ and $\rho_0 = |\psi_0\rangle\langle\psi_0|$.

The spectral gap estimation task requires disentangling errors from each of these four levels.

---

## 3. Rigorous Error Catalogue

### Error 1: Trotterization Bias ($\delta$-step discretization)

**Statement.** The effective generator of the channel $\Phi_\delta$ applied $r = t/\delta$ times differs from $\mathcal{L}$ by a systematic $O(\delta)$ correction to the eigenvalues.

**Analysis.** Write $\Phi_\delta = e^{\delta \mathcal{L}} + \delta^2 E$ where $\|E\|_\diamond = O(1)$. The eigenvalues of $\Phi_\delta$ are:

$$\mu_k(\Phi_\delta) = e^{\delta \lambda_k} + O(\delta^2) = 1 + \delta \lambda_k + O(\delta^2)$$

After $r$ steps, the eigenvalues of $\Phi_\delta^r$ are $\mu_k^r$. The effective decay rate extracted from fitting is:

$$\tilde{\lambda}_k = \frac{1}{\delta} \ln|\mu_k(\Phi_\delta)| = \frac{1}{\delta} \ln|1 + \delta \lambda_k + O(\delta^2)| = \lambda_k + O(\delta)$$

**Consequence.** The estimated spectral gap has a **systematic bias of order $\delta$**. This is a *first-order* Trotter error in the gap estimate, even though the diamond-norm error per step is $O(\delta^2)$.

**Bound.** $|\tilde{\lambda}_{\text{gap}} - \lambda_{\text{gap}}| \leq C \delta$ where $C$ depends on $\|\mathcal{L}\|$ and the second-order terms in the Baker-Campbell-Hausdorff expansion.

**Diagnostic.** Run simulations at multiple values of $\delta$ and extrapolate $\delta \to 0$. If the estimated gap depends linearly on $\delta$, this is the dominant bias. **[This is a leading candidate for a persistent systematic error.]**


### Error 2: Randomized Jump Selection Variance

**Statement.** At each $\delta$-step, a single jump $A^a$ is selected uniformly at random from $\mathcal{A}$, rather than applying the full convex combination. By Corollary III.1 of Chen et al. (2023), $\|\mathbb{E}_a[\Phi_\delta^{(a)}] - e^{\delta \mathcal{L}}\|_\diamond = O(\delta^2)$, so the *expectation* is correct to the same order. However, this introduces additional variance.

**Analysis.** Consider the eigenvalue of the randomized channel. At each step, the effective channel is one of $\{\Phi_\delta^{(a)}\}_{a \in \mathcal{A}}$, each having potentially different spectra. The product $\Phi_\delta^{(a_r)} \circ \cdots \circ \Phi_\delta^{(a_1)}$ is a random product of CPTP maps. The spectral properties of random products of matrices differ from those of the expected product.

For the *density matrix evolution*, the expectation over jump selections gives the correct channel. But for *trajectory simulations*, you are running a single realization of the jump sequence, so the randomness of jump selection is already folded into the trajectory noise.

**Consequence.** No additional systematic bias to the spectral gap beyond Error 1, but increased variance in the estimator. The variance per trajectory scales as $O(1/|\mathcal{A}|)$ per step relative to the deterministic (all-jumps) version.


### Error 3: Trajectory Sampling Noise (Monte Carlo Variance)

**Statement.** Given $N_{\text{traj}}$ trajectories, the estimator for $\text{Tr}[O \rho(t)]$ is:

$$\hat{O}(t) = \frac{1}{N_{\text{traj}}} \sum_{i=1}^{N_{\text{traj}}} \langle \psi_i(t) | O | \psi_i(t) \rangle$$

This is an unbiased estimator with variance:

$$\text{Var}[\hat{O}(t)] = \frac{1}{N_{\text{traj}}} \text{Var}_{\text{traj}}[\langle \psi | O | \psi \rangle]$$

**Critical point about the variance.** The single-trajectory variance $\text{Var}_{\text{traj}}[\langle \psi(t) | O | \psi(t) \rangle]$ does **not** vanish as $t \to \infty$. Even after convergence to steady state, each trajectory remains a pure state $|\psi_i(\infty)\rangle$ that fluctuates around the thermal ensemble. The asymptotic variance is:

$$\text{Var}_\infty = \mathbb{E}[\langle \psi | O | \psi \rangle^2] - (\text{Tr}[O \rho_\beta])^2$$

For a generic observable on a thermal state, this is $O(1)$ (not small).

**Consequence for gap estimation.** The signal you are fitting is:

$$\hat{O}(t) = \underbrace{\text{Tr}[O \rho_\beta]}_{\text{thermal average}} + \underbrace{c \, e^{-\lambda_2 t}}_{\text{signal}} + \underbrace{\eta(t)}_{\text{noise, } O(1/\sqrt{N_{\text{traj}}})}$$

There exists a **signal-to-noise crossover time** $t^*$ where $|c| e^{-\lambda_2 t^*} \sim 1/\sqrt{N_{\text{traj}}}$, i.e.:

$$t^* \sim \frac{1}{\lambda_2} \ln(|c| \sqrt{N_{\text{traj}}})$$

For $t > t^*$, the exponential tail is buried in noise. Fitting data that includes points beyond $t^*$ **biases the estimated gap upward** (the fit sees a faster apparent decay as it tries to match the noise floor). **[This is a leading candidate for persistent overestimation of the gap.]**

**Diagnostic.** Compute the noise floor $\sigma = \sqrt{\text{Var}_\infty / N_{\text{traj}}}$, and restrict the fit to $t < t^*$. Alternatively, use a weighted least squares fit where the weight at time $t$ accounts for the noise level.


### Error 4: Observable Overlap with the Gap Mode

**Statement.** The spectral decomposition of $\mathcal{L}$ gives (assuming non-degenerate spectrum for simplicity):

$$\text{Tr}[O \, e^{t \mathcal{L}}[\rho_0]] - \text{Tr}[O \rho_\beta] = \sum_{k=2}^{4^n} c_k \, e^{\lambda_k t}$$

where $c_k = \text{Tr}[O \, R_k] \cdot \text{Tr}[L_k \, \rho_0]$ involves the right ($R_k$) and left ($L_k$) eigenvectors of $\mathcal{L}$.

The gap mode $k=2$ is only visible if $c_2 \neq 0$. If $|c_2| \ll |c_3|$, then the dominant observed decay rate at intermediate times is $|\text{Re}(\lambda_3)|$, not $\lambda_{\text{gap}} = |\text{Re}(\lambda_2)|$.

**Consequence.** A poor choice of observable or initial state leads to **systematic overestimation** of the spectral gap. You observe $\lambda_3$ (or a higher mode) instead of $\lambda_2$.

**Diagnostic.** This error manifests as a gap estimate that is too large and independent of $N_{\text{traj}}$ and $\delta$. To check: compare multiple observables. If they give different gap estimates, the one giving the *smallest* gap is closest to the truth (assuming sufficient data quality). The right observable should have non-trivial overlap with all symmetry sectors (see Error 6).

**Mitigation.** Use observables that break all symmetries of $H$. For the Heisenberg chain, good candidates include:
- Single-site $Z$ operator on a boundary site (breaks translation invariance)
- A two-point correlator $Z_1 Z_{n/2}$ 
- The energy density on a specific bond

Bad candidates include global conserved quantities like total $S_z = \sum_i Z_i$, which commute with the Hamiltonian and thus may have zero overlap with the gap mode.


### Error 5: Symmetry Sector Restriction

**Statement.** This is the most subtle and, in my assessment, the **most likely cause of a persistent systematic error** for the Heisenberg chain.

The Heisenberg Hamiltonian $H = \sum_i \vec{S}_i \cdot \vec{S}_{i+1}$ has $\mathrm{SU}(2)$ symmetry: $[H, S_\alpha^{\text{tot}}] = 0$ for $\alpha = x, y, z$. At minimum, $[H, S_z^{\text{tot}}] = 0$.

**Effect on the Lindbladian.** The Lindbladian $\mathcal{L}_\beta$ with single-site Pauli jumps $X_i, Y_i, Z_i$ does connect different $S_z^{\text{tot}}$ sectors (since $X_i$ and $Y_i$ change $S_z$ by $\pm 1$). So generically, $\mathcal{L}_\beta$ is **primitive** (unique fixed point). However, there can be **approximate symmetries** or **near-degeneracies**.

Specifically, the operator Fourier transform $\hat{A}^a(\omega)$ of a Pauli operator involves energy projections $\sum_{E_i - E_j = \omega} P_{E_i} A^a P_{E_j}$. For the Heisenberg chain, the energy eigenstates organize into $\mathrm{SU}(2)$ multiplets. The Fourier-transformed jumps $\hat{Z}_i(\omega)$ preserve $S_z^{\text{tot}}$ while $\hat{X}_i(\omega), \hat{Y}_i(\omega)$ change it. But due to the finite-width Gaussian filter $f(t)$, these selection rules are only approximate, and transitions between different $S_z^{\text{tot}}$ sectors may be **exponentially suppressed** by the filter overlap.

**Mechanism for persistent error.** If the inter-sector transitions are slow (small matrix elements), the Lindbladian effectively has a **block structure** with fast intra-sector equilibration and slow inter-sector mixing. The true spectral gap $\lambda_2$ corresponds to the slowest inter-sector process. But:

1. If your initial state $|\psi_0\rangle$ lies in a definite $S_z^{\text{tot}}$ sector, it must transit through the slow modes to reach the full thermal state. The trajectory may appear to have equilibrated within its sector long before reaching the global equilibrium.

2. If your observable $O$ commutes with $S_z^{\text{tot}}$ (e.g., $O = Z_1$), it may not distinguish between intra-sector and inter-sector equilibration. The observable could show fast decay (intra-sector gap) while the actual state is far from the true fixed point.

**For the Heisenberg chain without external field**: $\mathrm{SU}(2)$ symmetry means the Gibbs state is block-diagonal in total spin sectors. The spectral gap of $\mathcal{L}_\beta$ could be determined by the rate of transitions *between* different total spin sectors, which can be much slower than the rate of equilibration *within* a sector.

**With an external magnetic field $h \sum_i Z_i$**: The $\mathrm{SU}(2)$ symmetry is broken to $\mathrm{U}(1)$ ($S_z^{\text{tot}}$ conservation only). The Lindbladian's Pauli-$X$ and $Y$ jumps break this $\mathrm{U}(1)$, but the inter-sector matrix elements scale with the filter function overlap at the relevant Bohr frequencies. If $h$ is small, these frequencies are nearly degenerate within multiplets, and the Gaussian filter may not distinguish them efficiently.

**Diagnostic.** Compare the gap estimate with and without the external field. If adding a symmetry-breaking field changes the estimated gap significantly, then symmetry-sector restriction is the issue. Also: check whether the Lindbladian you construct for $n \leq 6$ has near-degenerate eigenvalues close to $\lambda_2$ that might correspond to different symmetry sectors.


### Error 6: Finite Evolution Time (Transient Contamination)

**Statement.** At short times, the observable decay is:

$$\text{Tr}[O \rho(t)] - \text{Tr}[O \rho_\beta] = c_2 e^{\lambda_2 t} + c_3 e^{\lambda_3 t} + \cdots$$

If $|c_3/c_2|$ is not negligibly small, the early-time signal is dominated by the faster-decaying modes. Fitting an exponential to data in the range $[0, T]$ yields an effective rate that is a weighted average of $|\lambda_2|, |\lambda_3|, \ldots$, biased toward $|\lambda_3|$ (upward bias on the gap).

**Remedy.** Fit only to the **late-time** data window $[t_{\min}, t_{\max}]$ where $t_{\min}$ is chosen large enough that all modes except $k=2$ have decayed:

$$t_{\min} \gg \frac{1}{|\text{Re}(\lambda_3)| - |\text{Re}(\lambda_2)|}$$

and $t_{\max} < t^*$ (the noise crossover from Error 3). The existence of a valid fitting window requires:

$$t_{\min} < t^* \quad \Longleftrightarrow \quad \frac{1}{\lambda_3 - \lambda_2} \ll \frac{1}{\lambda_2} \ln(|c_2| \sqrt{N_{\text{traj}}})$$

This places a **lower bound on $N_{\text{traj}}$** for reliable gap estimation.


### Error 7: Non-normality of the Lindbladian

**Statement.** The CKG Lindbladian satisfies only *approximate* detailed balance (Chen et al. (2023), Lemma II.2), meaning the decomposition $D(\rho, \mathcal{L}) = H + A$ has a non-zero anti-Hermitian part $A$ with $\|A\|_{2 \to 2}$ bounded but potentially non-negligible.

**Consequence.** For non-normal generators, the short-time and intermediate-time behavior can differ qualitatively from the asymptotic exponential decay. Specifically:
- Eigenvalues can be complex, meaning the decay is oscillatory: $c_k e^{(\text{Re}(\lambda_k) + i \,\text{Im}(\lambda_k))t}$.
- Pseudospectral effects can cause transient amplification before eventual decay.
- The left and right eigenvectors can be far from orthogonal, making the coefficients $c_k$ ill-conditioned.

By Proposition II.3 of Chen et al. (2023), if $\lambda_1(H)/\lambda_{\text{gap}}(H) \leq 1/100$ (which requires near-exact detailed balance), then $t_{\text{mix}} \leq 3 \ln(3\|\rho^{-1/2}\|) / \lambda_{\text{gap}}(H)$, and the Hermitian gap controls the decay. But if the anti-Hermitian part is large relative to the gap, fitting a pure real exponential to oscillatory data will produce biased estimates.

**Diagnostic.** Check whether the exact eigenvalues $\lambda_2, \lambda_3, \ldots$ for your $n \leq 6$ systems have significant imaginary parts. If $|\text{Im}(\lambda_2)| \gtrsim |\text{Re}(\lambda_2)|$, the signal is oscillatory and a real exponential fit is inadequate. Use a damped oscillation model $c \, e^{-\gamma t} \cos(\omega t + \phi)$ instead, extracting $\gamma = |\text{Re}(\lambda_2)|$.

---

## 4. Additional Systematic Effects

### Effect A: Filter Width and Secular Approximation

The Gaussian filter $f(t)$ of width $\sigma_t$ (or equivalently energy resolution $\sigma_E = 1/\sigma_t$) controls how well the Lindbladian approximates the Davies generator. From Chen et al. (2023), the fixed point deviation scales as $O(\sqrt{\beta/T} \cdot t_{\text{mix}})$ where $T$ is the time-domain truncation. If $\sigma_E$ is too large (wide filter), nearby Bohr frequencies are not resolved, and the Lindbladian fails to satisfy even approximate detailed balance for closely-spaced energy levels.

For the Heisenberg chain, the energy spectrum can have dense regions (especially at the center of the band), leading to poor resolution. This affects the fixed point accuracy and, indirectly, the spectral gap (via Proposition E.3 of Chen et al. (2023): $\|\rho_{\text{fix}} - \rho_\beta\|_1 \leq 14\epsilon/\lambda_{\text{gap}}$).


### Effect B: Trajectory-Level Stiffness

The quantum trajectory equation involves a stochastic process where, at each $\delta$-step, a jump occurs with probability $O(\delta)$ and no jump occurs with probability $1 - O(\delta)$. For small $\delta$, most steps are "no-jump" (non-Hermitian evolution by $I - \frac{\delta}{2} \sum_j L_j^\dagger L_j$, followed by renormalization). The rare jumps cause sudden changes in $|\psi\rangle$.

If $\delta$ is too large, the jump probability per step approaches $O(1)$, and the Trotterization error (Error 1) becomes significant. If $\delta$ is too small, you need many steps to reach a given evolution time, increasing computational cost without improving the gap estimate.

**Optimal $\delta$.** Balance Trotter bias ($O(\delta)$ on the gap) against computational cost ($O(1/\delta)$ steps). For gap estimation to precision $\epsilon_{\text{gap}}$, need $\delta \lesssim \epsilon_{\text{gap}} / C$ where $C \sim \|\mathcal{L}\|$.

---

## 5. Comprehensive Error Budget

Let $\lambda_{\text{gap}}^{\text{est}}$ denote the estimated spectral gap. The total error decomposes as:

$$\lambda_{\text{gap}}^{\text{est}} - \lambda_{\text{gap}} = \underbrace{O(\delta)}_{\text{Trotter bias}} + \underbrace{O(\delta_{\text{disc}})}_{\text{frequency discretization}} + \underbrace{\text{fit bias}}_{\text{Errors 4, 5, 6}} + \underbrace{\text{noise}}_{\text{Error 3, } O(1/\sqrt{N_{\text{traj}}})}$$

where:

| Error Source | Type | Direction | Scales With |
|:---|:---|:---|:---|
| $\delta$-step Trotter (Error 1) | Systematic | Can go either way | $\delta \|\mathcal{L}\|^2$ |
| Randomized jump (Error 2) | Variance | Neither | $1/(|\mathcal{A}| \cdot N_{\text{traj}})$ |
| Trajectory noise (Error 3) | Variance + fit bias | Overestimates gap | $1/\sqrt{N_{\text{traj}}}$ |
| Observable overlap (Error 4) | Systematic | Overestimates gap | Observable-dependent |
| Symmetry restriction (Error 5) | Systematic | **Overestimates gap** | Symmetry structure |
| Transient modes (Error 6) | Systematic | Overestimates gap | $1/(t_{\min}(\lambda_3 - \lambda_2))$ |
| Non-normality (Error 7) | Systematic | Either direction | $\|A\|_{2 \to 2} / \lambda_{\text{gap}}$ |

**Key observation.** Errors 3, 4, 5, and 6 all bias the gap **upward** (overestimate). If you are seeing a persistent overestimate, these are the prime suspects. If you are seeing a persistent *underestimate*, Error 1 or Error 7 is more likely.

---

## 6. Diagnosis Protocol for the Persistent Error

Given that you are seeing a persistent error, here is a systematic procedure to isolate it:

**Step 1: Characterize the error direction.** At $n \leq 6$, compare $\lambda_{\text{gap}}^{\text{est}}$ from trajectories to exact diagonalization. Is the trajectory estimate consistently *above* or *below* the true value?

**Step 2: $\delta$-convergence test.** Run trajectory simulations at $\delta, \delta/2, \delta/4$ with everything else fixed. Plot $\lambda_{\text{gap}}^{\text{est}}$ vs $\delta$. If linear dependence is observed, the Trotter bias (Error 1) is significant. Extrapolate to $\delta = 0$ via Richardson extrapolation.

**Step 3: Noise floor test.** Increase $N_{\text{traj}}$ by a factor of 4 (which halves the noise standard deviation). If the gap estimate changes, Error 3 is significant. In particular, check whether the fitting range extends past $t^*$.

**Step 4: Observable test.** Try several observables: $Z_1$, $X_1$, $Z_1 Z_2$, the energy $H$, a random traceless Hermitian. If different observables give different gap estimates, Error 4 or 5 is active. The smallest gap estimate is closest to the truth.

**Step 5: Symmetry test (critical for Heisenberg chain).** For the Heisenberg chain:
- Examine the exact vectorized Lindbladian $\mathbb{L}$ at $n \leq 6$. Check if $\lambda_2$ and $\lambda_3$ (or a cluster of eigenvalues near $\lambda_2$) correspond to different $S_z^{\text{tot}}$ sectors.
- Add a small symmetry-breaking field $h \sum_i Z_i$ with $h \ll J$. If the gap estimate changes dramatically, the true gap mode involves inter-sector transitions that your observable/initial state is not probing.
- Start trajectories from an initial state that is a superposition of different $S_z^{\text{tot}}$ sectors (e.g., $|\psi_0\rangle = |+\rangle^{\otimes n}$, which is a superposition over all $S_z$ values).

**Step 6: Fitting window test.** Vary $t_{\min}$ of the fitting window. If the gap estimate decreases as $t_{\min}$ increases (up to the noise limit), Error 6 (transient contamination) is active.

**Step 7: Non-normality check.** For $n \leq 6$, check $|\text{Im}(\lambda_2)|$ and $\|A\|$. If the eigenvalues are complex, switch to a damped-oscillation fitting model.

---

## 7. Recommended Estimator

Based on the above analysis, the robust procedure for spectral gap estimation is:

**(a) Multi-observable approach.** Compute the trajectory-averaged time series $\hat{O}_j(t)$ for several observables $O_j$. For each, subtract the long-time average: $\Delta_j(t) := \hat{O}_j(t) - \bar{O}_j$ where $\bar{O}_j = \frac{1}{T_{\text{avg}}} \int_{T - T_{\text{avg}}}^{T} \hat{O}_j(t') dt'$.

**(b) Windowed exponential fit.** For each observable, fit $|\Delta_j(t)|$ to $c_j e^{-\gamma_j t}$ on the window $[t_{\min}, t_{\max}]$. Choose $t_{\min}$ adaptively: start with $t_{\min} = 0$ and increase until $\gamma_j$ stabilizes. Choose $t_{\max}$ as the point where the signal-to-noise ratio drops below some threshold (e.g., 3).

**(c) Take the minimum.** The best estimate of the spectral gap is $\hat{\lambda}_{\text{gap}} = \min_j \gamma_j$, since observables can only *overestimate* the gap (by missing the true gap mode).

**(d) Richardson extrapolation in $\delta$.** Repeat at two or more values of $\delta$ and extrapolate to $\delta = 0$: $\hat{\lambda}_{\text{gap}}(0) \approx 2\hat{\lambda}_{\text{gap}}(\delta/2) - \hat{\lambda}_{\text{gap}}(\delta)$.

**(e) Bootstrap confidence intervals.** Resample the trajectories with replacement and repeat the fitting to get error bars.

---

## 8. Summary of Most Likely Causes of Persistent Error

In order of likelihood for the Heisenberg chain:

1. **Symmetry sector restriction (Error 5):** The Heisenberg chain's $\mathrm{SU}(2)$ symmetry means the true gap mode likely involves transitions between total-spin sectors that may be weakly coupled by the Fourier-transformed Pauli jumps. This would cause a persistent *overestimate* of the gap.

2. **Trotter bias (Error 1):** A finite $\delta$ shifts the effective eigenvalues by $O(\delta)$. This is easy to check via the convergence test and can bias in either direction.

3. **Trajectory noise floor masking the tail (Error 3 interacting with Error 6):** If $N_{\text{traj}}$ is insufficient, the exponential tail is lost in noise, and the fit captures faster-decaying modes. This causes *overestimation*.

4. **Observable overlap (Error 4):** If the observable is orthogonal to the gap eigenvector, you see a subleading mode. Multi-observable comparison diagnoses this immediately.

---

## References (from project knowledge)

- **Chen et al. (2023):** "Quantum Thermal State Preparation." Theorem III.1 (weak measurement), Corollary III.1 (randomized simulation), Proposition II.3 (mixing time from Hermitian gap), Proposition E.3 (fixed point accuracy), Lemma E.2 (Hermitian gap controls decay).
- **Ramkumar and Soleimanifar (2024):** "Mixing time of quantum Gibbs sampling for random sparse Hamiltonians." Theorem 2.1-2.3 (spectral gap bounds for graph-local and random jumps), matrix Bernstein inequality argument for concentration.
- **Chen et al. (2025):** "An efficient and exact noncommutative quantum Gibbs sampler." Exact KMS detailed balance construction, Metropolis coherent term.
- **Ding et al. (2025):** "End-to-End Efficient Quantum Thermal and Ground State Preparation Made Simple." Theorem 11 (time averaging), effective Lindblad dynamics from weak coupling, $O(\alpha^2)$ per-step effective evolution time.
- **Li and Wang:** "Simulating Markovian open quantum systems using higher-order series expansion." Theorem 11 (higher-order Lindbladian simulation), Kraus operator construction.
- **Lin (2025):** "Dissipative Preparation of Many-Body Quantum State." Eq. (6)-(7) for first-order Lindblad simulation error.
