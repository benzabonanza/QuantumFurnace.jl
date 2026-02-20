# Spectral Gap Estimation: Diagnosis and Improvement Pipeline

## Context for the AI Agent

You are working on a Julia codebase for quantum Gibbs sampling via the CKG Lindbladian (Chen, Kastoryano, Brandão, Gilyén, 2023). The system of interest is a Heisenberg spin chain $H = \sum_i \vec{S}_i \cdot \vec{S}_{i+1}$ (with or without an external field $h \sum_i Z_i$) as a first option, but later on Ising chain, and 2D Heisenberg will be added too. The Lindbladian uses single-site Pauli operators $\{X_i, Y_i, Z_i\}$ as jump operators, which are then operator-Fourier-transformed with a Gaussian filter to produce the actual Lindblad operators $\hat{A}^a(\omega)$.

There are two existing capabilities in the codebase:

1. **Exact Lindbladian construction**: For small systems ($n \leq 6$ qubits), the full vectorized Lindbladian $\mathbb{L} \in \mathbb{C}^{4^n \times 4^n}$ can be constructed and diagonalized to obtain the exact spectral gap $\lambda_{\mathrm{gap}} = -\mathrm{Re}(\lambda_2)$, where $\lambda_1 = 0$ is the steady-state eigenvalue. For $n = 6$, the observed spectral gap is approximately 0.11.

2. **Quantum trajectory simulator**: A Monte Carlo wavefunction simulator that evolves pure states $|\psi(t)\rangle$ through discrete $\delta$-steps. At each step, a random jump operator $A^a$ is selected, and the weak measurement scheme from Chen et al. (2023) Theorem III.1 / Corollary III.1 is applied. This involves two levels of stochastic sampling: (a) random selection of which jump $A^a$ to apply, and (b) within that $\delta$-step, whether a quantum jump occurs or the system evolves via the no-jump (non-Hermitian) branch. The trajectory-averaged expectation value $\hat{O}(t) = \frac{1}{N_{\mathrm{traj}}} \sum_i \langle \psi_i(t) | O | \psi_i(t) \rangle$ converges to $\mathrm{Tr}[O \, e^{t\mathcal{L}}[\rho_0]]$ as $N_{\mathrm{traj}} \to \infty$ and $\delta \to 0$.

**The problem**: Spectral gap estimates from trajectory simulations show a persistent systematic error compared to exact diagonalization at $n = 6$. The goal of this task is to (a) diagnose the source of this error and (b) implement improved estimation methods to obtain reliable spectral gap estimates that can be trusted at larger system sizes ($n > 8$) where exact diagonalization is infeasible.

---

## Overall Goal

Build a comprehensive diagnostic module and an improved spectral gap estimator, both validated against exact results at $n = 6$ (and ideally $n = 8$ via sparse diagonalization). The module should produce clear diagnostic plots that reveal which error source dominates, and the estimator should produce a spectral gap value with quantified uncertainty.

---

## Part 1: Diagnostic Infrastructure

### Task 1.1: Exact Reference Data at $n = 4, 6$

**What to compute and why.**

We need a complete picture of the exact Lindbladian spectrum to interpret all subsequent diagnostics. The vectorized Lindbladian $\mathbb{L}$ is a $4^n \times 4^n$ matrix acting on vectorized density matrices. Its eigenvalues $\{\lambda_k\}$ satisfy $\mathrm{Re}(\lambda_k) \leq 0$ for all $k$, with $\lambda_1 = 0$ corresponding to the steady state. The spectral gap is $\lambda_{\mathrm{gap}} = -\mathrm{Re}(\lambda_2)$, i.e., the real part of the eigenvalue second-closest to zero determines how fast the slowest-decaying mode relaxes.

Compute and store the following quantities from the exact vectorized Lindbladian at $n = 6$:

1. **The leading 20–30 eigenvalues** of $\mathbb{L}$, sorted by real part (closest to zero first). Use a sparse eigenvalue solver like `Arpack.jl` targeting eigenvalues nearest to zero (shift-invert mode). Store both real and imaginary parts. The imaginary parts tell us whether the gap mode is oscillatory, which would make a purely real exponential fit inappropriate.

2. **The corresponding right eigenvectors** (at minimum for eigenvalues $\lambda_2$ and $\lambda_3$). These are the "modes" of the Lindbladian. We need them to compute observable overlaps (Task 1.3).

3. **The Gibbs state** $\rho_\beta = e^{-\beta H} / \mathrm{Tr}(e^{-\beta H})$, vectorized. This is the steady state of the Lindbladian.

4. **The fixed point** $\rho_{\mathrm{fix}}$ of the Lindbladian, which is the right eigenvector of $\mathbb{L}$ corresponding to $\lambda_1 = 0$, reshaped back into a density matrix. Compute $\|\rho_{\mathrm{fix}} - \rho_\beta\|_1$ (trace distance). This quantifies the approximate detailed balance violation at the level of the fixed point. The user reports this is $\sim 10^{-6}$, but we should verify.

**Why this matters**: Every subsequent diagnostic references these exact values. Without them, we cannot distinguish estimation errors from true physics.


### Task 1.2: Anti-Hermitian Defect Diagnosis

**What to compute and why.**

The CKBG Lindbladian satisfies only *approximate* detailed balance when we have a TrotterDomain implementation of it with Trotter time evolutions in $B$ and the OFT's. This means that when we decompose the Lindbladian's action in the KMS inner product defined by the Gibbs state $\rho_\beta$, there is a non-zero anti-Hermitian component. Specifically, define the similarity-transformed generator (see Chen et al. 2023, Section II and Appendix E - note that this was written up for their approximate GNS DB case without the coherent term B, we assume - hopefully correctly - that the implemented approx. KMS DB algorithm also fulfils this):

$$D(\rho_\beta, \mathcal{L}) = \rho_\beta^{-1/4} \, \mathcal{L}[\rho_\beta^{1/4} \, (\cdot) \, \rho_\beta^{1/4}] \, \rho_\beta^{-1/4}$$

This is a superoperator that, when vectorized, becomes a $4^n \times 4^n$ matrix $\mathbb{D}$. Decompose it as $\mathbb{D} = \mathbb{H} + \mathbb{A}$ where $\mathbb{H} = \frac{1}{2}(\mathbb{D} + \mathbb{D}^\dagger)$ is the Hermitian part and $\mathbb{A} = \frac{1}{2}(\mathbb{D} - \mathbb{D}^\dagger)$ is the anti-Hermitian part. (Use the same vectorization methods that we already used for the run_lindbladian() simulations where we construct the full Lindbladian)

Compute and report:

1. $\|\mathbb{A}\|$ (operator norm of the anti-Hermitian part). This is the "size" of the detailed balance violation.
2. $\lambda_{\mathrm{gap}}(\mathbb{H})$ (spectral gap of the Hermitian part alone, i.e., its second-largest eigenvalue minus its largest, where "largest" means closest to zero since all eigenvalues are $\leq 0$).
3. The ratio $\|\mathbb{A}\| / \lambda_{\mathrm{gap}}(\mathbb{H})$. If this ratio is $\ll 0.01$, the Lindbladian is effectively normal and the Hermitian gap controls the transient decay (Proposition II.3 of Chen et al. 2023). If this ratio is $\gtrsim 0.1$, non-normality effects can cause oscillatory transients and the purely real exponential fitting model may be inadequate.
4. The imaginary parts of the leading eigenvalues of $\mathbb{L}$: $|\mathrm{Im}(\lambda_2)| / |\mathrm{Re}(\lambda_2)|$ and similarly for $\lambda_3$, $\lambda_4$. If these ratios are significant ($> 0.1$), the transient signal is oscillatory and fitting should use a damped-oscillation model $c \, e^{-\gamma t} \cos(\omega t + \phi)$ rather than a pure exponential.

**Implementation note**: The similarity transform involves $\rho_\beta^{\pm 1/4}$, which requires computing the matrix fourth root. Since $\rho_\beta$ is diagonal in the energy eigenbasis, this is straightforward: diagonalize $H$ to get $\rho_\beta$ in diagonal form, then raise the diagonal entries to the $\pm 1/4$ power. Now in TrotterDomain runs, we have everything in Trotter basis, the Gibbs state and the Linbdladian, so there come up with an efficient way to raise the Gibbs state to 1/4 powers.
The superoperator $D$ acts on vectorized operators $|X\rangle\!\rangle$ as $\mathbb{D} |X\rangle\!\rangle = |\rho_\beta^{-1/4} \mathcal{L}[\rho_\beta^{1/4} X \rho_\beta^{1/4}] \rho_\beta^{-1/4}\rangle\!\rangle$. In practice, form the matrix $\mathbb{D}$ by applying this map to each basis element of the vectorized operator space.

**Why this matters**: If the anti-Hermitian defect ratio is large, the observable decay signal has oscillatory components that a real exponential fit will systematically misestimate. This would explain a persistent error independent of $N_{\mathrm{traj}}$ or $\delta$. If the ratio is small (as we expect for the Gaussian-filtered CKG Lindbladian at moderate $\beta$ with $\sigma_E = \beta^{-1}$), we can safely use real exponential models and focus on other error sources.


### Task 1.3: Observable Overlap Analysis

**What to compute and why.**

The trajectory-averaged observable signal after subtracting the steady-state value is a sum of exponentials:

$$\hat{\Delta}(t) = \mathrm{Tr}[O \, e^{t\mathcal{L}}[\rho_0]] - \mathrm{Tr}[O \, \rho_\beta] = \sum_{k=2}^{4^n} c_k \, e^{\lambda_k t}$$

where $c_k = \mathrm{Tr}[O \, R_k] \cdot \mathrm{Tr}[L_k^\dagger \, (\rho_0 - \rho_\beta)]$ involves the right eigenvector $R_k$ and left eigenvector $L_k$ of $\mathcal{L}$. The coefficient $c_2$ is the "overlap" of the observable and initial state with the gap mode. If $|c_2|$ is small relative to $|c_3|$, the gap mode is invisible and any fit will return $|\mathrm{Re}(\lambda_3)|$ or worse.

For each candidate observable $O$ and initial state $\rho_0 = |\psi_0\rangle\langle\psi_0|$, compute the overlap coefficients $c_k$ for the leading 10–20 modes using the exact eigenvectors from Task 1.1. Specifically:

1. Use the right eigenvectors $R_k$ (reshaped into $2^n \times 2^n$ matrices from the vectorized form) and left eigenvectors $L_k$ of $\mathbb{L}$.
2. For each observable-initial-state pair, compute $c_k = \mathrm{Tr}[O \, R_k] \cdot \mathrm{Tr}[L_k^\dagger \, (\rho_0 - \rho_\beta)]$.
3. Report the ratio $|c_3/c_2|$, $|c_4/c_2|$, etc.

Test the following observables (all are Hermitian, traceless after subtracting the thermal average):

- $Z_1$ (Pauli-Z on the first site): simple, local, breaks translation symmetry.
- $X_1$ (Pauli-X on the first site): couples to different symmetry sectors than $Z_1$.
- $Z_1 Z_{\lfloor n/2 \rfloor}$ (two-point correlator): may have better overlap with long-wavelength modes.
- $H$ itself (the Hamiltonian): the energy is a natural slow observable.
- A random traceless Hermitian matrix (as a control).

You can use the build preset observables function already and can delete previous variants that are not in the above list. 

Test the following initial states:

- $|0\rangle^{\otimes n}$ (all spins up in Z basis): definite $S_z^{\mathrm{tot}} = n/2$ sector.
- $|+\rangle^{\otimes n}$ (all spins in X-plus state): superposition of all $S_z^{\mathrm{tot}}$ sectors.
- The maximally mixed state $I/2^n$: has support in all sectors by construction.

**Why this matters**: If the best observable-initial-state pair gives $|c_3/c_2| \gg 1$, then the fitting window must exclude a very long initial transient ($t_{\min} \sim \ln(|c_3/c_2|) / (\lambda_3 - \lambda_2)$), which may not leave enough signal before the noise floor. Conversely, if $|c_3/c_2| \lesssim 1$, the two-exponential fit should work immediately. This analysis also reveals whether symmetry sectors are causing the problem: if $c_2 = 0$ for a $Z$-basis initial state but $c_2 \neq 0$ for the $|+\rangle^{\otimes n}$ initial state, then the gap mode lives in a sector that the $Z$-basis state cannot access.


### Task 1.4: Lindbladian-Trotter Bias Diagnosis ($\delta$-convergence)

**What to compute and why.**

Each $\delta$-step of the trajectory simulator implements a CPTP map $\Phi_\delta$ that approximates $e^{\delta \mathcal{L}}$. The diamond-norm error per step is $O(\delta^2)$, but the *bias on the extracted spectral gap* is $O(\delta)$. This is because the eigenvalues of $\Phi_\delta$ are $\mu_k = 1 + \delta \lambda_k + O(\delta^2)$, and extracting the rate via $\tilde{\lambda}_k = \frac{1}{\delta} \ln |\mu_k|$ gives $\tilde{\lambda}_k = \lambda_k + O(\delta)$.

Run trajectory simulations at $n = 6$ with **three values of $\delta$**: a baseline value $\delta_0$ (whatever is currently used), $\delta_0/2$, and $\delta_0/4$. Keep everything else identical: same $N_{\mathrm{traj}}$, same total evolution time $t_{\mathrm{final}}$, same observable, same initial state, same random seed for reproducibility where possible.

For each $\delta$, extract the spectral gap estimate using the same method (we will improve the method later, but for this diagnostic use whatever is currently in the code). Plot $\hat{\lambda}_{\mathrm{gap}}$ vs $\delta$.

**Expected outcomes and interpretation**:

- If $\hat{\lambda}_{\mathrm{gap}}(\delta)$ shows a clear linear trend in $\delta$, the Trotter bias is significant and Richardson extrapolation to $\delta = 0$ is needed: $\hat{\lambda}_{\mathrm{gap}}(0) \approx 2\hat{\lambda}_{\mathrm{gap}}(\delta/2) - \hat{\lambda}_{\mathrm{gap}}(\delta)$.
- If the three values are indistinguishable within statistical error bars, $\delta$ is already small enough and the persistent error comes from elsewhere.

**Why this matters**: This is the easiest error source to check and to fix. If Trotter bias is the dominant error, the fix (Richardson extrapolation) is trivial and requires no changes to the fitting methodology.

---

## Part 2: The Effective Rate Plot ($\lambda_{\mathrm{eff}}(t)$)

### Task 2.1: Compute and Plot $\lambda_{\mathrm{eff}}(t)$

**What this is and why it is the single most useful diagnostic.**

The effective rate plot is a time-local, model-free estimate of the instantaneous decay rate of the observable signal. It requires no fitting whatsoever and immediately reveals the three time regimes: transient contamination (early times, rate too high), the golden window (intermediate times, plateau at the true gap), and the noise floor (late times, erratic behavior).

Given the trajectory-averaged observable signal $\hat{O}(t)$ from $N_{\mathrm{traj}}$ trajectories, define:

$$\hat{\Delta}(t) = \hat{O}(t) - \hat{O}_{\mathrm{ss}}$$

where $\hat{O}_{\mathrm{ss}}$ is the steady-state thermal expectation value. For the exact reference at $n = 6$, use $\hat{O}_{\mathrm{ss}} = \mathrm{Tr}[O \rho_\beta]$ from Task 1.1. For larger systems, estimate $\hat{O}_{\mathrm{ss}}$ as the time average of $\hat{O}(t)$ over the last $\sim 20\%$ of the simulation (after equilibration).

The effective rate is:

$$\lambda_{\mathrm{eff}}(t) = -\frac{1}{\tau} \ln \left| \frac{\hat{\Delta}(t + \tau)}{\hat{\Delta}(t)} \right|$$

where $\tau$ is the lag. Using $\tau = \delta t$ (the time spacing between data points) gives the highest time resolution but the noisiest estimate. Using $\tau = 3\text{--}5 \, \delta t$ smooths the noise at the cost of some time resolution. The choice depends on the noise level; start with $\tau = 3\delta t$ and adjust.

**Implementation details**:

```julia
function effective_rate(Δ::Vector{Float64}, dt::Float64; lag::Int=3)
    # Δ[i] is the trajectory-averaged signal minus steady-state value at time (i-1)*dt
    # lag is the number of time steps used for the finite-difference log-derivative
    # Returns (t_values, λ_eff_values) where NaN indicates invalid points
    
    n = length(Δ)
    t_vals = [(i-1)*dt for i in 1:(n-lag)]
    λ_eff = Vector{Float64}(undef, n-lag)
    
    τ = lag * dt
    for i in 1:(n-lag)
        d_now = Δ[i]
        d_later = Δ[i + lag]
        # Only compute if both values have the same sign and are nonzero
        # (sign change means the signal has crossed zero, which happens when 
        #  noise dominates or when oscillatory modes are present)
        if d_now != 0.0 && d_later != 0.0 && sign(d_now) == sign(d_later)
            λ_eff[i] = -log(abs(d_later / d_now)) / τ
        else
            λ_eff[i] = NaN
        end
    end
    
    return t_vals, λ_eff
end
```

**Plot specification**: Plot $\lambda_{\mathrm{eff}}(t)$ vs $t$ with the following overlays:

- A horizontal dashed line at the exact $\lambda_{\mathrm{gap}}$ from Task 1.1 (only for $n = 6$ where the exact value is known).
- A horizontal band showing the $\pm 1\sigma$ statistical uncertainty on $\hat{\Delta}(t)$, propagated through to $\lambda_{\mathrm{eff}}$. The uncertainty on $\lambda_{\mathrm{eff}}$ at time $t$ is approximately $\sigma_{\lambda}(t) \approx \frac{1}{\tau} \sqrt{\frac{\mathrm{Var}[\hat{\Delta}(t)]}{|\hat{\Delta}(t)|^2} + \frac{\mathrm{Var}[\hat{\Delta}(t+\tau)]}{|\hat{\Delta}(t+\tau)|^2}}$, where $\mathrm{Var}[\hat{\Delta}(t)] = \mathrm{Var}_{\mathrm{traj}}[\langle\psi(t)|O|\psi(t)\rangle] / N_{\mathrm{traj}}$.
- Vertical dashed lines marking the estimated $t_{\min}$ (where $\lambda_{\mathrm{eff}}$ first enters the plateau) and $t_{\max}$ (where the error bars on $\lambda_{\mathrm{eff}}$ exceed some threshold, e.g., $30\%$ of the plateau value). These define the golden fitting window.

**What to look for in the plot**:

- A clear plateau at the exact $\lambda_{\mathrm{gap}}$ confirms that the method works and identifies the correct fitting window. The previous persistent error was likely caused by fitting data from outside this window.
- A plateau at a value *above* the exact $\lambda_{\mathrm{gap}}$ suggests the gap mode ($\lambda_2$) has poor observable overlap and a higher mode is being measured. Cross-reference with Task 1.3.
- No plateau at all (monotonically decreasing $\lambda_{\mathrm{eff}}$ until the noise floor) means either $N_{\mathrm{traj}}$ is insufficient or the gap between $\lambda_2$ and $\lambda_3$ is too small for the modes to separate within the available signal-to-noise.
- A plateau at a value *below* the exact $\lambda_{\mathrm{gap}}$ would be surprising and would indicate a bug or a problem with the steady-state subtraction.
- Oscillatory behavior in $\lambda_{\mathrm{eff}}$ (beyond statistical noise) indicates complex eigenvalues, confirming that the anti-Hermitian defect (Task 1.2) is relevant and a damped-oscillation fit model is needed.

**Why this matters**: This single plot is the Rosetta Stone for the entire problem. It diagnoses transient contamination, noise floor limitations, and fitting window selection simultaneously, with zero model assumptions. Every subsequent fitting method in this pipeline should be validated against what this plot shows.


### Task 2.2: Error Bar Estimation for $\lambda_{\mathrm{eff}}(t)$ via Bootstrap

To get reliable error bars on $\lambda_{\mathrm{eff}}(t)$, use bootstrap resampling over trajectories. The procedure is:

1. You have $N_{\mathrm{traj}}$ individual trajectory time series $\{O_i(t)\}_{i=1}^{N_{\mathrm{traj}}}$.
2. For each bootstrap sample $b = 1, \ldots, N_{\mathrm{boot}}$ (use $N_{\mathrm{boot}} = 200$), draw $N_{\mathrm{traj}}$ trajectories with replacement, compute the bootstrap-averaged signal $\hat{O}^{(b)}(t)$, subtract the steady-state value, and compute $\lambda_{\mathrm{eff}}^{(b)}(t)$.
3. At each time $t$, the error bar on $\lambda_{\mathrm{eff}}(t)$ is the standard deviation across bootstrap samples.

This correctly propagates the trajectory-level variance into the effective rate estimate, accounting for correlations between different time points in the same trajectory.

---

## Part 3: Improved Spectral Gap Estimator (Transient Regime)

### Task 3.1: Two-or-multi-Exponential Fit

**What this is and why it solves the contamination problem.**

The trajectory-averaged signal is a sum of exponentials $\hat{\Delta}(t) = \sum_k c_k e^{\lambda_k t}$. A single-exponential fit $\hat{\Delta}(t) \approx c_1 e^{-\gamma_1 t}$ attempts to capture this entire sum with one term. At any finite time, the faster modes ($k \geq 3$) contribute and pull $\gamma_1$ upward (toward $|\mathrm{Re}(\lambda_3)|$ or a weighted average of the faster rates). This is the "excited-state contamination" problem well-known in lattice QCD.

The two-exponential fit models the signal as:

$$\hat{\Delta}(t) \approx c_1 \, e^{-\gamma_1 t} + c_2 \, e^{-\gamma_2 t}$$

with the constraint $0 < \gamma_1 < \gamma_2$. Here $\gamma_1$ is our estimate of $\lambda_{\mathrm{gap}}$ and the second term is a nuisance exponential that absorbs contamination from all faster modes. We do not need $\gamma_2$ to accurately represent any single eigenvalue; it just needs to soak up the fast transient contribution so that $\gamma_1$ can settle on the true gap.

**Implementation details**:

Use a nonlinear least-squares fit (e.g., `LsqFit.jl` in Julia). The fit should be performed on a *selected time window* $[t_{\min}, t_{\max}]$ informed by the effective rate plot from Task 2.1:

- $t_{\min}$: Can be set to $0$ for the two-exponential model (since the second exponential handles early-time contamination), or to a small positive value to exclude the very first few steps where $O(\delta^2)$-per-step errors accumulate nonlinearly.
- $t_{\max}$: Set to the point where the effective rate plot shows the error bars on $\lambda_{\mathrm{eff}}$ becoming larger than $\sim 30\%$ of the plateau value, or equivalently, where $|\hat{\Delta}(t)| < 3 \sigma_{\hat{\Delta}}(t)$ (signal-to-noise ratio drops below 3).

**Initialization strategy (critical for convergence)**:

Two-exponential fits are notoriously sensitive to initial conditions. Use the following robust initialization that avoids local minima:

1. From the effective rate plot, read off the plateau value as $\gamma_1^{(0)}$ and the early-time value (at $t \approx 0$) as $\gamma_2^{(0)}$.
2. For the amplitudes, use: $c_1^{(0)} = \hat{\Delta}(t_{\mathrm{mid}}) / e^{-\gamma_1^{(0)} t_{\mathrm{mid}}}$ where $t_{\mathrm{mid}}$ is in the middle of the plateau region, and $c_2^{(0)} = \hat{\Delta}(0) - c_1^{(0)}$.
3. Pass bounded constraints to the optimizer: $\gamma_1 \in [0.01, \gamma_1^{(0)} \cdot 3]$, $\gamma_2 \in [\gamma_1, 10 \cdot \gamma_2^{(0)}]$, $c_1$ and $c_2$ unconstrained in sign but bounded in magnitude by $10 \cdot |\hat{\Delta}(0)|$.

Alternatively, for a more robust (but slightly more complex) initialization, use the **Prony two-point method**: pick two time points $t_a$ and $t_b = t_a + \tau$ in the signal (not too early, not too late). Form $r_1 = \hat{\Delta}(t_a + \tau/2) / \hat{\Delta}(t_a)$ and $r_2 = \hat{\Delta}(t_b) / \hat{\Delta}(t_a + \tau/2)$. The two decay constants $z_1 = e^{-\gamma_1 \tau/2}$ and $z_2 = e^{-\gamma_2 \tau/2}$ satisfy $z_1 + z_2 = r_1 + r_2$ and $z_1 z_2 = r_1 r_2$ (approximately, in the noiseless two-exponential case). Solve this quadratic to get initial guesses for $\gamma_1, \gamma_2$.

**Error bars via bootstrap**: Use the same bootstrap resampling from Task 2.2. For each bootstrap sample, recompute the trajectory average, redo the two-exponential fit, and record $\gamma_1^{(b)}$. The standard deviation of $\{\gamma_1^{(b)}\}$ gives the statistical uncertainty on the spectral gap estimate. Also report the mean and median of the bootstrap distribution, since the mean can be biased by occasional failed fits (where the optimizer lands in a local minimum). The median is more robust.

**Validation**: At $n = 4, 6$, compare the extracted $\gamma_1$ (with bootstrap error bars) to the exact $\lambda_{\mathrm{gap}}$. The result should agree within $2\sigma$. If it doesn't, check whether the discrepancy has the same sign and magnitude as the Trotter bias (Task 1.4), in which case apply Richardson extrapolation (Task 3.2) and check again.


### Task 3.2: Richardson Extrapolation in $\delta$

**What this is and why it eliminates the leading Trotter bias.**

The $\delta$-step channel $\Phi_\delta$ approximates $e^{\delta \mathcal{L}}$ with a $O(\delta^2)$ diamond-norm error per step. Over $r = t/\delta$ steps, the total error in the extracted decay rate is $O(\delta)$ (the per-step $O(\delta^2)$ error accumulates over $O(1/\delta)$ steps). This means:

$$\hat{\lambda}_{\mathrm{gap}}(\delta) = \lambda_{\mathrm{gap}} + a \cdot \delta + O(\delta^2)$$

for some constant $a$ that depends on the Lindbladian's higher-order structure. Richardson extrapolation eliminates the leading-order term by combining estimates at two values of $\delta$:

$$\hat{\lambda}_{\mathrm{gap}}^{\mathrm{Rich}} = 2 \hat{\lambda}_{\mathrm{gap}}(\delta/2) - \hat{\lambda}_{\mathrm{gap}}(\delta)$$

This cancels the $O(\delta)$ term and gives a result with $O(\delta^2)$ residual bias. If three $\delta$ values are available ($\delta$, $\delta/2$, $\delta/4$), one can verify the linear scaling by checking that the three points lie on a line when plotted against $\delta$, and use the Richardson-extrapolated value from the two smallest $\delta$ values as the final estimate.

**Implementation**: This requires running the trajectory simulation at multiple $\delta$ values (Task 1.4 already generates this data). For each $\delta$, apply the two-exponential fit from Task 3.1 to extract $\hat{\lambda}_{\mathrm{gap}}(\delta)$. Then combine via the formula above. Propagate uncertainties assuming the estimates at different $\delta$ are independent (they use different trajectory data):

$$\sigma_{\mathrm{Rich}} = \sqrt{4\sigma^2(\delta/2) + \sigma^2(\delta)}$$

where $\sigma(\delta)$ is the bootstrap uncertainty from Task 3.1 at step size $\delta$.


### Task 3.3: Fitting Window Optimization

**What this is and why automatic window selection matters.**

The quality of the two-exponential fit depends on the time window $[t_{\min}, t_{\max}]$. If $t_{\max}$ extends into the noise floor, the fit tries to accommodate noise as signal and biases the result. If $t_{\min}$ is too large, we throw away useful data. If $t_{\min}$ is too small and the signal has contributions from more than two exponentials at early times, the two-exponential model is also insufficient.

Implement the following automatic window selection:

1. **$t_{\max}$ selection**: Compute the signal-to-noise ratio $\mathrm{SNR}(t) = |\hat{\Delta}(t)| / \sigma_{\hat{\Delta}}(t)$ where $\sigma_{\hat{\Delta}}(t)$ is the standard error of $\hat{\Delta}(t)$ across trajectories. Set $t_{\max}$ as the last time point where $\mathrm{SNR}(t) > 3$.

2. **$t_{\min}$ selection via stability test**: Perform the two-exponential fit over $[t_{\min}, t_{\max}]$ for a sequence of values $t_{\min} = 0, \delta t, 2\delta t, \ldots$. Plot the extracted $\gamma_1$ vs $t_{\min}$. A robust estimate produces a $\gamma_1$ that is *stable* (does not change beyond statistical error) as $t_{\min}$ varies. If $\gamma_1$ decreases as $t_{\min}$ increases (before the error bars blow up), the two-exponential model is insufficient and there is still contamination from a third mode. Choose $t_{\min}$ as the smallest value where $\gamma_1$ has stabilized.

The $t_{\min}$ stability plot is another diagnostic that belongs in the thesis. It demonstrates that the result is robust to the choice of fitting window. Should definitely figure out if we always have a good window, or if somehow the window size becomes constraining and stops us from making a good spectral gap estimate.

---

## Part 4: Symmetry Sector Analysis

### Task 4.1: Symmetry Labels on Lindbladian Eigenvalues

**What to compute and why.**

The isotropic Heisenberg chain $H = \sum_i \vec{S}_i \cdot \vec{S}_{i+1}$ commutes with total $S_z = \sum_i Z_i$ and total $\vec{S}^2$. The Lindbladian, acting on the space of density matrices, has a block structure related to these symmetries. Specifically, a density matrix element $|E_i\rangle\langle E_j|$ carries quantum number $\Delta S_z = S_z(E_i) - S_z(E_j)$.

The single-site $Z$ jump operators preserve $S_z^{\mathrm{tot}}$ (and thus $\Delta S_z$), while $X$ and $Y$ jump operators change $S_z^{\mathrm{tot}}$ by $\pm 1$ (and thus change $\Delta S_z$). After Fourier transformation, the Lindblad operators $\hat{A}^a(\omega)$ still respect these selection rules (up to the smearing introduced by the Gaussian filter).

For the exact Lindbladian at $n = 6$:

1. Compute $S_z^{\mathrm{tot}}$ for each energy eigenstate of $H$.
2. Label each Lindbladian eigenvector $R_k$ by the $\Delta S_z$ quantum number it carries. Formally, $R_k$ as a density matrix has support on matrix elements $|E_i\rangle\langle E_j|$ with $S_z(E_i) - S_z(E_j) = \Delta S_z(k)$. For the vectorized Lindbladian, this means examining which blocks of the eigenvector have support.
3. Report the $\Delta S_z$ labels for $\lambda_2, \lambda_3, \lambda_4, \ldots$. If $\lambda_2$ carries $\Delta S_z = 0$ (same as the steady state), it describes decay of populations. If it carries $\Delta S_z \neq 0$, it describes decay of coherences.

**Why this matters**: If the gap mode ($\lambda_2$) lives in a $\Delta S_z \neq 0$ sector, then observables that are diagonal in the energy basis (like $Z_1$ or $H$) have exactly zero overlap with it ($c_2 = 0$), and you are measuring a higher mode. This is particularly relevant for the Heisenberg chain where the spectrum organizes into $\mathrm{SU}(2)$ multiplets. An external field $h \sum_i Z_i$ breaks $\mathrm{SU}(2)$ to $\mathrm{U}(1)$ and can rearrange which mode is slowest.

### Task 4.2: Test with External Field

Run the complete diagnostic (Tasks 1.1 through 2.1) at $n = 4, 6$ both without and with a small external field $h = 0.1 J$ (where $J$ is the Heisenberg coupling strength). Compare the exact spectral gaps, the effective rate plots, and the two-exponential fit results. If the field significantly changes the estimated gap from trajectory simulations (while the exact gap changes only mildly), the symmetry sector issue is the primary error source.

---

## Part 5: Comprehensive Validation and Output

### Task 5.1: Produce a Summary Dashboard

After all diagnostics are computed, produce a single figure (or a small set of figures) that summarizes the findings. The dashboard should include:

1. **Panel A**: Exact Lindbladian spectrum — real vs imaginary parts of the leading 20 eigenvalues, with $\Delta S_z$ labels as colors/markers.
2. **Panel B**: Anti-Hermitian defect metrics — $\|\mathbb{A}\| / \lambda_{\mathrm{gap}}(\mathbb{H})$ and $|\mathrm{Im}(\lambda_2)/\mathrm{Re}(\lambda_2)|$.
3. **Panel C**: Observable overlap coefficients $|c_k|$ for the best observable-initial-state pair, showing how much signal comes from each mode.
4. **Panel D**: The effective rate plot $\lambda_{\mathrm{eff}}(t)$ with the exact gap overlaid, golden window marked, and bootstrap error bands.
5. **Panel E**: $\delta$-convergence plot — $\hat{\lambda}_{\mathrm{gap}}$ vs $\delta$ with Richardson extrapolation.
6. **Panel F**: Two-exponential fit results — the fitted curve overlaid on the data, with $\gamma_1$ and its bootstrap confidence interval compared to the exact value.
7. **Panel G**: $t_{\min}$ stability plot — $\gamma_1$ vs $t_{\min}$ showing the plateau.

### Task 5.2: Final Spectral Gap Estimate at $n = 4, 6$

Report the best estimate of $\lambda_{\mathrm{gap}}$ from trajectories, obtained by combining the two-exponential fit with Richardson extrapolation, along with:

- The bootstrap $1\sigma$ confidence interval.
- The exact value from diagonalization.
- The discrepancy in units of $\sigma$.

If the discrepancy is $< 2\sigma$, the method is validated. If not, check which of the diagnostics (anti-Hermitian defect, observable overlap, symmetry sector) suggests the remaining error source, and report which it is.

### DEFERRED TO A LATER MILESTONE: Task 5.3: Extend to $n = 8$ (Sparse Diagonalization)

We will do this at a later point. Would be too much to pack this into the same milestone as the rest.

For $n = 8$, the vectorized Lindbladian is $65536 \times 65536$. This is too large for dense diagonalization but amenable to sparse eigenvalue computation using KrylovKit.jl methods. For this we need to write up a fast apply_lindbladian function but possibly also an apply_algorithmic_lindbladian that is the effective weak measurement Lindbladian evolution faitfhul to Chen's quantum algorithm. Compute the leading 10 eigenvalues of the exact Lindbladian at $n = 8$ and repeat the validation from Task 5.2. This gives a second independent reference point before extrapolating to larger systems.



---

## Summary of Error Sources and Expected Impact

For reference, here is a table of all error sources, how each diagnostic detects them, and what the fix is. This should guide interpretation of the results.

**Trotter bias ($O(\delta)$)**: Detected by Task 1.4 ($\delta$-convergence). Fixed by Richardson extrapolation (Task 3.2). Direction: can bias either way. Typically the easiest to fix.

**Transient contamination from faster modes**: Detected by the effective rate plot (Task 2.1) showing $\lambda_{\mathrm{eff}}$ decreasing toward a plateau. Fixed by the two-exponential fit (Task 3.1). Direction: overestimates the gap. This is the most likely dominant error source.

**Trajectory sampling noise / noise floor**: Detected by the effective rate plot showing erratic behavior at late times. Mitigated by increasing $N_{\mathrm{traj}}$ or restricting $t_{\max}$ (Task 3.3). Direction: biases toward overestimation if the fit extends into the noise floor.

**Poor observable overlap with gap mode**: Detected by Task 1.3 (overlap coefficients). Fixed by choosing a better observable. Direction: overestimates the gap by measuring a faster mode.

**Symmetry sector restriction**: Detected by Task 4.1 (symmetry labels) and Task 4.2 (field comparison). Fixed by choosing observables/initial states that couple to the correct sector, or by adding a small symmetry-breaking field. Direction: overestimates the gap.

**Anti-Hermitian defect (non-normality)**: Detected by Task 1.2 (the ratio $\|\mathbb{A}\|/\lambda_{\mathrm{gap}}(\mathbb{H})$). If significant, requires switching to a damped-oscillation fit model. Direction: causes oscillatory bias in real exponential fits.

---

## Code Organization Suggestion

Create a module (or a set of scripts) organized as follows:

```
src/diagnostics/
    exact_reference.jl      # Tasks 1.1, 1.2, 1.3, 4.1: exact Lindbladian analysis
    trotter_convergence.jl   # Task 1.4: δ-convergence runs
    effective_rate.jl        # Tasks 2.1, 2.2: λ_eff computation and bootstrap
    two_exp_fit.jl           # Tasks 3.1, 3.2, 3.3: fitting and Richardson
    symmetry_analysis.jl     # Tasks 4.1, 4.2: symmetry sector labels
    dashboard.jl             # Task 5.1: summary plots
```

Each script should be independently runnable and should save its results (numerical values and plots) to a designated output directory. The dashboard script reads the saved results and composes the summary figure.
