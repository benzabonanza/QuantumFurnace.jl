# Domain Pitfalls: v2.2 Hamiltonian Simulation Time Counting

**Domain:** Adding quantum algorithm cost accounting (Hamiltonian simulation time counting) to an existing classical quantum Gibbs sampler simulator (QuantumFurnace.jl). The new code computes how expensive the *quantum algorithm* would be, NOT how expensive the classical simulation is.
**Researched:** 2026-03-04
**Confidence:** HIGH (grounded in direct codebase analysis + Chen 2023/2025 paper formulas + existing OFT/B/transition code inspection)

**Relationship to prior research:** This document covers pitfalls specific to the v2.2 milestone. The v2.0 PITFALLS.md covered codebase restructure risks and the v2.1 PITFALLS.md covered threading and mixing time estimation. Those pitfalls remain relevant but are not repeated here. This document covers NEW pitfalls introduced by adding Hamiltonian simulation time counting to the existing system.

**Key distinction driving all pitfalls:** The simulator already computes OFT, B terms, transition weights, and QPE grids for *classical simulation*. The cost counter must compute the *quantum algorithm's* cost from these same parameters. The simulator's truncations, optimizations, and NUFFT shortcuts are valid for simulation but do NOT reflect the quantum computer's actual work. Confusing "what our simulator does" with "what the quantum algorithm does" is the central risk of this milestone.

---

## Critical Pitfalls

Mistakes that produce silently wrong cost estimates, potentially by orders of magnitude.

---

### CRIT-01: Truncated Energy Labels vs Full QPE Grid -- The Classical Simulator Does NOT Use the Full Grid

**What goes wrong:**
The existing code creates energy labels via `_create_energy_labels(num_energy_bits, w0)` which produces a full grid of `2^num_energy_bits` points: `w0 * [-N/2 : N/2-1]`. But immediately after, `_truncate_energy_labels(energy_labels, config)` (in `furnace_utensils.jl` lines 152-201) aggressively truncates this grid to only the range where the integrand `gamma(w) * G(w,nu1) * G(w,nu2)` exceeds `1e-12`. For Gaussian transitions with narrow sigma, this can reduce thousands of grid points to a few hundred.

The classical simulator is correct to truncate: the integrand truly IS negligible outside this range, and skipping zero-contribution points saves computation. But the *quantum algorithm* running on a quantum computer uses the **full** `2^r` QPE grid. Every grid point requires a controlled-Hamiltonian simulation of duration `|t_k| = |k * t0|`. The quantum computer cannot skip grid points because QPE is a unitary operation that coherently processes all time points simultaneously.

**The pitfall:** If the cost counter sums `|t_k|` only over the truncated `energy_labels` (as returned by `_precompute_labels`), it will massively undercount the quantum algorithm's Hamiltonian simulation time. The truncated set might have 200 points while the full set has 4096. Since the sum `sum(|t_k|)` over the full grid is `t0 * N^2/4`, underestimating N by a factor of ~4-5 underestimates the total cost by a factor of ~16-25.

**Why it happens:**
The existing `_precompute_data` functions return the truncated `energy_labels` and never store the full grid separately. A developer implementing cost counting will naturally reach for `precomputed_data.energy_labels` since it is the array used everywhere in the simulator. It feels correct -- it is the energy grid the code uses. But the cost counter needs the grid the *quantum computer* uses.

**Concrete example from the codebase:**
For `num_energy_bits=12, w0=0.05, sigma=0.08, beta=1.0`:
- Full grid: 4096 points from -102.4 to +102.35
- Truncated grid (from `_truncate_energy_labels`): typically ~200-400 points in [-1, 1]
- `sum(|k|)` for full grid: `t0 * 4096^2 / 4` ~ `t0 * 4.2e6`
- `sum(|k|)` for truncated grid: `t0 * ~4e4` (100x smaller)

**However, there is a subtlety about the energy sum itself:**
The OFT cost is `sum_omega [sum_k |t_k|]` where the outer sum is over energy labels omega. The question is: does the *quantum algorithm* iterate over the full energy grid or the truncated one? On a quantum computer, the algorithm evaluates the channel at each energy omega in sequence (classically controlled), and CAN skip omegas where gamma(omega) is negligible. So the OUTER sum (over omega) uses the truncated energy grid, but the INNER sum (over time/QPE) uses the full QPE grid. Conflating these two grids is the core mistake.

**Prevention:**
- The cost counter MUST use `_create_energy_labels(num_energy_bits, w0)` (the FULL grid) for summing QPE time costs per energy evaluation.
- The cost counter SHOULD use the truncated energy grid for counting how many energy evaluations are performed.
- Create a clear separation: `simulation_energy_labels` (truncated, for counting outer sum) vs `qpe_time_points` (full, for inner QPE time sum). Never mix them.
- Add a docstring on the cost counting function: "This uses the FULL QPE grid for time sums, not the truncated simulation grid."
- Write a test: for the inner sum, verify against the closed form `sum_{k=-N/2}^{N/2-1} |k| * t0 = t0 * N * (N-2) / 4 + t0 * N/2`.

**Detection:**
- Cost estimates that seem suspiciously low -- orders of magnitude below what the paper quotes.
- Cost estimates where the inner time sum changes with `sigma` (it should NOT -- sigma controls the truncation/filter, not the QPE grid).
- Sanity check: `inner_time_sum ≈ t0 * N^2 / 4` for `N = 2^num_energy_bits`.

**Which phase:** Must be addressed in the FIRST phase of cost counting. This is the foundational grid computation.

---

### CRIT-02: OFT Cost Must Sum |t_k| Without Filter Function Weighting

**What goes wrong:**
The quantum algorithm implements the OFT (Operator Fourier Transform) as:

```
A(omega) = sum_k f(t_k) * e^{iHt_k} * A * e^{-iHt_k}
```

where `f(t_k)` is the filter function (Gaussian envelope: `exp(-sigma^2 * t_k^2)`), and each term requires Hamiltonian simulation for time `|t_k|`. The *total Hamiltonian simulation time* for one OFT evaluation at frequency `omega` is:

```
T_OFT(omega) = 2 * sum_k |t_k|
```

The factor of 2 comes from the `e^{+iHt}` and `e^{-iHt}` evolutions (forward and backward). The filter function `f(t_k)` determines the *amplitude* of each term's contribution, NOT whether the time evolution actually runs. On a quantum computer, QPE implements ALL grid points coherently. The filter weights appear in the success probability, not in the gate count.

**Three common mistakes:**

1. **Forgetting the factor of 2:** The `e^{+iHt}` requires one Hamiltonian simulation forward, and `e^{-iHt}` requires one backward. Each costs `|t_k|` in Hamiltonian simulation time. Total is `2 * |t_k|` per term. Looking at the existing `time_oft!` code (ofts.jl lines 22-72), it computes `caches.U.diag .= exp.(1im .* hamiltonian.eigvals .* t)` -- this is the eigenvalue-domain shortcut for `e^{iHt}`, but the quantum computer implements both directions explicitly.

2. **Weighting by the filter function:** Multiplying each `|t_k|` by `exp(-sigma^2 * t_k^2)` before summing. This is what the NUFFT prefactors do (`base_weights = exp.(-sigma^2 .* time_labels .^ 2)` in nufft.jl line 44). Using these weights for cost counting gives the *classical simulation's* effective contribution, NOT the quantum computer's gate count. The quantum computer does NOT skip terms with small filter values.

3. **Using the truncated time grid:** The simulator uses `_truncate_time_labels_for_oft` (in `time_domain.jl`) to drop time points where `exp(-t^2 * sigma^2)` falls below `1e-12`. This is valid for the classical Riemann sum approximation. The quantum computer's QPE uses ALL `2^r` time points.

**Prevention:**
- Define the OFT cost as a pure function of the grid parameters:
  ```julia
  function oft_ham_sim_time(num_energy_bits::Int, t0::Float64)
      N = 2^num_energy_bits
      # Full QPE grid: k = -N/2, ..., N/2-1
      # |t_k| = |k| * t0
      # sum_{k=-N/2}^{N/2-1} |k| = N^2/4 - N/2 + N/2 = N^2/4 (for even N)
      # More precisely: 2 * sum_{k=1}^{N/2-1} k + N/2 = N*(N-2)/4 + N/2
      total_abs_time = t0 * (N * (N - 2) / 4 + N / 2)
      return 2.0 * total_abs_time  # factor of 2 for forward + backward evolution
  end
  ```
- This function should have NO dependency on `sigma`, `transition`, or any truncation.
- Test against a brute-force loop: `2 * t0 * sum(abs(k) for k in -N/2:N/2-1)`.
- Document the factor-of-2 convention prominently.

**Detection:**
- Cost estimate that changes with `sigma` is WRONG (sigma affects the classical simulation's accuracy, not the quantum algorithm's gate count).
- Cost estimate that is orders of magnitude smaller than the exact formula `T ~ t0 * N^2 / 2` is suspicious.

**Which phase:** OFT cost counting phase. Must include factor-of-2 from the start.

---

### CRIT-03: B Coherent Term Cost Is a Double Sum Over Time, Not a Single Sum

**What goes wrong:**
The coherent B term (Chen 2025, KMS construction) involves two nested time sums:

```
B = t0^2 * sum_t b_minus(t) * sum_s b_plus(s) * [time evolution terms involving s and t]
```

Looking at the existing `B_time` function (coherent.jl lines 105-148), the structure is:

1. **Inner loop over `s`** (keys of `b_plus`): For each `s`, compute `U(s*beta) * A' * U(-2*s*beta) * A * U(s*beta)`. The Hamiltonian simulation time for this inner step involves:
   - `U(s*beta)` = evolution for time `|s * beta|`
   - `U(-2*s*beta)` = evolution for time `|2*s * beta|`
   - `U(s*beta)` again = evolution for time `|s * beta|`
   - Total inner evolution per (s): `4 * |s * beta|`

2. **Outer loop over `t`** (keys of `b_minus`): For each `t`, compute `U(t/sigma)' * [inner result] * U(t/sigma)`. This requires Hamiltonian simulation for time `2 * |t / sigma|` per outer step.

The total Hamiltonian simulation time for B is:

```
T_B = sum_t sum_s [4 * |s * beta| + 2 * |t / sigma|]
    = |b_minus_terms| * sum_s 4*|s*beta| + |b_plus_terms| * sum_t 2*|t/sigma|
```

This is O(N_outer * N_inner) in the number of loop iterations but the *time* decomposes into two independent sums scaled by the opposite count.

**Common mistakes:**

1. **Counting only the outer sum:** Computing `sum_t |t/sigma|` and ignoring the inner sum over `s`. The inner sum can have hundreds of terms, each requiring Hamiltonian simulation.

2. **Double-counting the t0^2 factor:** The existing code returns `B .* t0^2` (coherent.jl line 148). This `t0^2` is the Riemann sum weight (dx_inner * dx_outer), NOT a Hamiltonian simulation time. The cost counter should account for the actual time evolution durations, not the Riemann sum weights.

3. **Confusing truncated dictionaries with full sums:** The existing code computes `b_minus` and `b_plus` as truncated dictionaries (via `_compute_truncated_func` in coherent.jl lines 228-232) that drop terms below `1e-12`. For cost counting, you need to know how many terms the *quantum algorithm* would need. The truncation is an approximation -- the quantum computer must implement enough terms for the desired accuracy. However, unlike the OFT case, the B term's truncation IS algorithmically valid: the quantum computer also truncates terms with negligible b-function values, because B is implemented as a finite sum (not a QPE circuit).

**Connection to existing code:**
The `B_trotter` function (coherent.jl lines 152-198) computes `num_t0_steps = Int(round(s * beta / trotter.t0))`, giving the discrete time in Trotter steps. The cost counter should use these actual step counts, not the continuous time values, for Trotter-specific cost reporting.

**Prevention:**
- Write the B cost as an explicit double loop over the truncated b_minus and b_plus dictionaries:
  ```julia
  function b_term_ham_sim_time(b_minus_dict, b_plus_dict, beta, sigma)
      n_outer = length(b_minus_dict)
      n_inner = length(b_plus_dict)
      inner_time_sum = sum(4 * abs(s * beta) for s in keys(b_plus_dict))
      outer_time_sum = sum(2 * abs(t / sigma) for t in keys(b_minus_dict))
      return n_outer * inner_time_sum + n_inner * outer_time_sum
  end
  ```
- Actually compute the truncated dicts using `_compute_truncated_func` (this is cheap, O(N) where N = 2^r).
- Verify against the known complexity from the paper: the B term cost should scale roughly as `O(beta^2 / sigma)` times the number of grid points.

**Detection:**
- B cost that is smaller than OFT cost is possible but should be checked -- for KMS with coherent term at large beta, B can be significant.
- B cost that does not depend on `beta` is WRONG.
- B cost that has only single-sum scaling O(N) instead of a product-of-sums indicates a missed nesting level.

**Which phase:** B coherent term cost phase. Should come after OFT cost is validated.

---

### CRIT-04: Weak Measurement Channel Factor of 2 -- U and Controlled-U^dagger

**What goes wrong:**
The quantum Gibbs sampler implements weak measurements via the CPTP channel (Chen Eq. 3.2):

```
E(rho) = K0 * rho * K0^dagger + delta * sum_w L_w * rho * L_w^dagger + U_res * rho * U_res^dagger
```

On the quantum computer, implementing this channel requires:
1. **Computing `A(w)` for each `w`:** This is the OFT, which costs `2 * sum_k |t_k|` per omega value.
2. **Implementing the channel:** The weak measurement protocol uses `U_channel` (implementing the dissipative part) and `controlled-U_channel^dagger` (for the non-demolition measurement). This requires TWO applications:
   - One forward: `U` (implementing the channel)
   - One controlled backward: `controlled-U^dagger` (for measurement reversal)

This means each step of the quantum algorithm costs **2x** the Hamiltonian simulation time of a single OFT computation for the dissipative part.

**The easy-to-miss factor:**
Looking at the existing `_build_cptp_channel` (furnace_utensils.jl lines 183-201), the code computes `alpha = 1 - sqrt(1 - delta)` and constructs `K0 = I - alpha * R`. This is the *mathematical* form of the channel. The *quantum implementation* doubles the OFT cost.

**Connection to the coherent unitary:**
The KMS construction also applies `U_B = exp(-i * delta * B)` (see `_precompute_coherent_unitary` in coherent.jl lines 57-99). The coherent unitary is a single unitary application (NOT controlled), so it contributes 1x its Hamiltonian simulation time, not 2x. Applying the 2x factor uniformly to both the dissipative part AND the coherent part is wrong -- only the dissipative channel gets the 2x.

**Prevention:**
- Clearly separate the cost into three components per step:
  1. OFT dissipative cost: `2 * T_OFT` (U + controlled-U^dagger)
  2. B coherent cost: `1 * T_B` (single unitary)
  3. Total per step: `2 * T_OFT + T_B`
- Document which factor applies to which component.
- Write a test: `per_step_cost == 2 * oft_cost + b_cost`, NOT `2 * (oft_cost + b_cost)`.

**Detection:**
- If the total cost is reported as `2 * (T_OFT + T_B)` instead of `2 * T_OFT + T_B`, the B term is over-counted. For large beta where B is significant, this error can be substantial.
- Check: does the cost increase when switching from GNS (no B) to KMS (with B)? It should increase by exactly `T_B`, not `2 * T_B`.

**Which phase:** Must be addressed when combining OFT and B costs into per-step totals.

---

### CRIT-05: Confusing Classical Simulation Mixing Time with Quantum Algorithm Mixing Time

**What goes wrong:**
The `config.mixing_time` parameter controls how long the *classical simulation* runs to evolve the density matrix. The v2.1 milestone added `estimate_mixing_time` which extrapolates the mixing time from simulation data. These are properties of the *classical simulation's* approximation to the Lindbladian dynamics.

The quantum algorithm's mixing time is the mixing time of the *ideal* Lindbladian on the quantum computer. For an exact implementation, the quantum and classical mixing times should agree (up to approximation errors from Trotter, QPE truncation, etc.). But the classical simulation has additional sources of error:

1. **Energy grid truncation** (CRIT-01) -- makes the classical Lindbladian slightly different from the ideal one.
2. **Riemann sum approximation** -- the discrete sum over omega approximates a continuous integral. Coarser grids give a different Lindbladian.
3. **Floor effect** (documented in MEMORY.md) -- the trace distance plateaus at a delta-dependent floor, not at zero. The floor scales as `k_energy * delta + floor_Trotter`. This floor affects the classical mixing time estimate but is irrelevant to the quantum algorithm's ideal mixing time.
4. **Trotter error** -- the TrotterDomain simulator approximates time evolution. The quantum computer uses the same Trotter circuit, so this error IS shared.
5. **Bohr frequency collisions** (documented in MEMORY.md) -- even-n chains with Z-only disorder have exact spectral symmetry causing collisions. This affects both classical and quantum mixing times.

**The pitfall:** Using the classically estimated mixing time (from `estimate_mixing_time`) as the quantum algorithm's mixing time for cost counting. If the classical estimate is 150 time units due to a floor effect, but the ideal Lindbladian converges in 100 time units, the cost estimate is 50% too high. Conversely, if the classical Lindbladian has a larger spectral gap than the ideal one (due to truncation effects), the classical estimate could be too optimistic.

**Prevention:**
- The cost counter should accept `mixing_time` as an explicit input parameter, not derive it from the simulation. The user provides their best estimate of the quantum algorithm's mixing time.
- Provide utilities to compute "what-if" costs for a range of mixing times: `compute_cost(config, ham; mixing_times=10:10:200)`.
- Document clearly: "The mixing time used for cost estimation is the quantum algorithm's mixing time, which may differ from the classical simulation's observed mixing time due to approximation errors."
- If using the classical estimate, provide a warning when the estimated floor is a significant fraction of the target epsilon.
- Consider providing a separate field like `mixing_time_source::Symbol` (`:user_specified` vs `:classical_estimate`) in the cost result struct.

**Detection:**
- Cost estimate that changes dramatically when `delta` is varied (delta affects the classical floor but not the quantum mixing time).
- Cost estimate that gives different answers for BohrDomain vs TrotterDomain for the same Hamiltonian (they approximate the same quantum algorithm differently).
- Compare the classical estimate from `estimate_mixing_time` against the Krylov spectral gap (`run_krylov_spectrum`) -- significant discrepancies indicate the classical estimate is unreliable for cost counting.

**Which phase:** Must be addressed from the start. This is a design-level decision about the cost counting API.

---

### CRIT-06: Transition Weight Function Differences Affect the Effective Energy Grid Size

**What goes wrong:**
The codebase supports three transition weight families:

1. **Gaussian** (`with_linear_combination=false`): `gamma(w) = exp(-(w + w_gamma)^2 / (2 * sigma_gamma^2))`. Infinite support, Gaussian decay.
2. **Metropolis** (`a=0, with_linear_combination=true`): `gamma(w) = exp(-beta * max(w + beta*sigma^2/2, 0))`. One-sided exponential decay with a kink.
3. **Smooth Metropolis** (`a>0, b>=0`): Intermediate smoothness. Tails determined by `a` and `b` parameters.

**The critical difference for cost counting:**
The number of *effective* energy grid points (those with non-negligible `gamma(w)`) varies dramatically between transition types. `_truncate_energy_labels` uses `cutoff=1e-12` to determine this. For the same `num_energy_bits`:
- Gaussian with narrow sigma_gamma: few effective points (Gaussian decays fast)
- Metropolis: asymmetric, more effective points on the negative-w side
- Smooth Metropolis: intermediate

The total OFT cost = `n_effective_energies * T_OFT_per_energy`. Getting `n_effective_energies` wrong scales the total cost linearly.

**The subtlety about symmetrization:**
The existing `_truncate_energy_labels` SYMMETRIZES the truncated range (line 197: `sym_limit = max(abs(start), abs(end))`). This means Metropolis transitions keep energy points where `gamma(w) = 0` because the symmetric range extends to cover both tails. This is the classical simulation's choice (for Hermitian jump operator symmetry exploitation). The quantum algorithm could potentially use an asymmetric range, but for correctness must cover at least the range where `gamma(w)` is non-negligible.

**For cost counting:**
The cost counter should use the same truncated energy grid that the quantum algorithm would use. The symmetrized grid from `_truncate_energy_labels` is a reasonable choice (conservative -- includes some zero-contribution points for Metropolis), but the cost counter should document this.

**Prevention:**
- Use `_truncate_energy_labels` (or equivalent) for counting the number of energy evaluations. This is what the quantum algorithm also needs.
- For the per-energy QPE time cost, use the full `2^r` QPE grid (CRIT-01). Do NOT confuse the energy grid truncation with the time grid truncation.
- Add test cases for all three transition types, verifying the effective energy count changes appropriately.
- Report `n_effective_energies` as a diagnostic in the cost struct.

**Detection:**
- n_effective_energies that does not change between Gaussian and Metropolis transitions (it should change significantly).
- n_effective_energies that equals the full grid size `2^r` (truncation should reduce it for reasonable parameters).

**Which phase:** OFT cost implementation phase. Affects the outer sum count.

---

## Moderate Pitfalls

Issues that cause confusing or inaccurate results, but not orders-of-magnitude errors.

---

### MOD-01: Delta Step Size vs Total Cost -- Off-by-One and Delta-Dependence

**What goes wrong:**
The total quantum algorithm cost is:

```
T_total = ceil(mixing_time / delta) * T_per_step
```

Two common issues:

1. **Off-by-one in step count:** Is it `floor(mixing_time / delta)` or `ceil(mixing_time / delta)`? The existing code uses `num_steps = Int(config.mixing_time / config.delta)` (implicit truncation for exact division). For cost counting, `ceil` is more conservative (guarantees convergence).

2. **Per-step cost independence from delta:** The OFT cost per step does NOT depend on delta. The B coherent term cost per step also does NOT depend on delta. Only the *number* of steps depends on delta. This means the cost counter cleanly factors as `n_steps * per_step_cost`. If delta sneaks into `per_step_cost`, something is wrong.

3. **Total cost cancels delta in the limit:** For small delta, `n_steps ~ mixing_time / delta` and `per_step_cost` is delta-independent. So `T_total ~ mixing_time * per_step_cost / delta`. Wait -- this grows without bound as delta shrinks. This is physically correct: smaller delta means more steps, each with the same QPE cost. The key question for users: what is the "right" delta? The answer depends on the desired approximation accuracy of the weak measurement. The cost counter should accept delta as input, not choose it.

**Prevention:**
- Use `ceil(Int, mixing_time / delta)` for step count.
- Report both `cost_per_step` and `total_cost` separately so users can adjust mixing_time or delta independently.
- Add a note: "Per-step cost is independent of delta. Total cost = ceil(mixing_time/delta) * per_step_cost."
- Provide a utility for parameter sweeps: varying delta shows total cost scaling as 1/delta.

**Detection:**
- Total cost that does not scale linearly with `1/delta` indicates a bug.
- Per-step cost that changes with `delta` is WRONG.

**Which phase:** Final assembly phase.

---

### MOD-02: The Rescaling Factor Converts Between Physical and Rescaled Hamiltonians

**What goes wrong:**
The HamHam constructor rescales the Hamiltonian so its spectrum fits in `[0, 0.45]` (hamiltonian.jl lines 394-405). The `rescaling_factor` field stores this conversion: `H_physical = rescaling_factor * (H_rescaled - shift * I)`.

ALL Bohr frequencies, eigvals, and derived quantities in the HamHam are in rescaled units. The parameters `t0`, `w0`, `sigma` in Config are also in rescaled units. The cost counter's output is naturally in rescaled units.

**The trap:** If the cost counter reports in rescaled units without stating this, and the user compares costs across Hamiltonians with different rescaling factors, the comparison is meaningless. A 3-qubit Hamiltonian (rescaling ~ 20) and a 6-qubit Hamiltonian (rescaling ~ 200) would appear to have similar costs in rescaled units when the physical costs differ by 10x.

**Prevention:**
- Report in rescaled units (matching all other internal time quantities) AND store `rescaling_factor` in the result struct.
- Add a helper: `physical_time(cost::HamSimCost) = cost.total_ham_sim_time * cost.rescaling_factor`.
- For cross-system comparisons, always convert to physical units first.
- The struct docstring must state: "All times are in rescaled Hamiltonian units (spectrum normalized to [0, 0.45]). Multiply by rescaling_factor for physical Hamiltonian time."

**Detection:**
- Cross-system comparison where a larger system appears cheaper (likely a rescaling factor issue).
- Cost that does not change when the overall coupling strength of the Hamiltonian is scaled (it should scale linearly via `rescaling_factor`).

**Which phase:** API and struct design phase. Must be decided early.

---

### MOD-03: Per-Jump vs Summed-Over-Jumps Cost Reporting Ambiguity

**What goes wrong:**
The quantum algorithm applies one randomly chosen jump operator per step (uniform over `n_jumps` operators). Looking at `run_thermalize` (furnace.jl): `idx = rand(rng, 1:length(jumps))` -- ONE jump is selected per step.

The per-step OFT cost is for ONE jump at all (truncated) energies, not for all jumps. The QPE time sum (`sum_k |t_k|`) is the same regardless of which jump is chosen -- it depends only on the grid, not the jump operator.

**But there is a subtlety with the B coherent term:** The existing code computes per-jump B_k for Thermalize mode (`_precompute_coherent_unitary`, coherent.jl lines 57-99) but total B = sum_k B_k for Lindbladian mode. The quantum algorithm applies `exp(-i * delta * B_k)` for the chosen jump k, so the per-step B cost is for ONE jump's B_k.

Since all jumps have the same OFT cost (same grid, same QPE circuit) and the same B_k cost structure (same grid, same b_minus/b_plus), the per-step cost is jump-independent for symmetric jump sets (Pauli operators on periodic chains). For asymmetric jump sets, report worst-case.

**Prevention:**
- Report per-step cost as the cost for ONE jump operator's channel.
- Total cost: `ceil(mixing_time / delta) * T_per_step_one_jump`.
- Add a breakdown: `per_step_oft_cost`, `per_step_b_cost`, `per_step_total`, `num_steps`, `total_cost`.

**Detection:**
- Total cost that scales with `n_jumps` is suspicious (per-step cost should be jump-independent for Pauli jumps).

**Which phase:** API design phase.

---

### MOD-04: GNS vs KMS Construction Changes Which Cost Components Exist

**What goes wrong:**
The GNS construction has no coherent B term (`with_coherent(GNS()) = false`), while KMS has it (`with_coherent(KMS()) = true`). The cost counter must handle both:

- GNS: `T_per_step = 2 * T_OFT`
- KMS: `T_per_step = 2 * T_OFT + T_B`

**Additionally**, GNS uses a shifted transition weight `gamma(w + beta * sigma^2 / 2)` (documented in `_pick_transition_gns` docstring). This shift changes the effective energy support, affecting `n_effective_energies`. The cost counter must dispatch on the construction type for both:
1. Whether B cost is included (structural difference)
2. Which transition function determines the energy truncation (parametric difference)

**Prevention:**
- Dispatch on `config.construction`:
  ```julia
  per_step_cost(config::Config{<:Any, <:Any, GNS}, ...) = 2 * oft_cost
  per_step_cost(config::Config{<:Any, <:Any, KMS}, ...) = 2 * oft_cost + b_cost
  ```
- Use `Union{Nothing, Float64}` for `b_cost` in the result struct (nothing for GNS).
- For energy truncation, use `pick_transition(config)` which already dispatches correctly on GNS/KMS.

**Detection:**
- GNS result with non-nothing B cost is WRONG.
- KMS result with nothing B cost is WRONG.

**Which phase:** API design phase.

---

### MOD-05: Domain Independence -- Same Quantum Algorithm, Different Classical Computation

**What goes wrong:**
The four domain types (BohrDomain, EnergyDomain, TimeDomain, TrotterDomain) are different levels of approximation for the *classical simulation*. They ALL approximate the SAME quantum algorithm. The cost counter should produce the same result regardless of which domain the classical simulator uses.

**The complication:** BohrDomain and EnergyDomain do not have `t0` (it is `nothing` in their Config). The cost counter needs `t0` to compute QPE time sums. For these domains, `t0` must be computed from the Fourier dual: `t0 = 2*pi / (2^num_energy_bits * w0)`.

Additionally, `num_energy_bits` is `nothing` for BohrDomain (which operates at the Bohr frequency level, not a discretized grid). BohrDomain users must specify `num_energy_bits` explicitly for cost counting.

**Prevention:**
- The cost counter should accept `num_energy_bits`, `w0` (or `t0`) as explicit parameters, independent of the domain.
- If extracting from Config, gracefully handle `nothing` values:
  ```julia
  r = something(config.num_energy_bits, error("num_energy_bits required for cost counting"))
  t0_val = something(config.t0, 2*pi / (2^r * config.w0))
  ```
- For TrotterDomain, additionally report `num_trotter_steps` as a separate metric.
- Document: "The cost counter computes the quantum algorithm's cost, independent of classical simulation domain."

**Detection:**
- Running cost computation with different domains gives different results for the same underlying parameters: WRONG.
- BohrDomain user getting an error because `config.t0 === nothing`: handle gracefully.

**Which phase:** API design phase.

---

### MOD-06: HamHam Does Not Store beta -- Must Be Passed Separately

**What goes wrong:**
`beta` (inverse temperature) is not a field on `HamHam{T}`. It is used during construction (`HamHam(raw, beta)`) to compute the Gibbs state, but the scalar value is not retained. The cost counter needs `beta` for:
- B term cost (b_plus depends on `beta`)
- Transition weight computation (gamma depends on `beta`)
- Energy truncation (truncation range depends on `beta`)

**Prevention:**
- Require `beta` as a mandatory argument to the cost function.
- Alternatively, accept a `Config` object which already contains `beta`.
- The cleanest API: `compute_ham_sim_cost(config::Config, hamiltonian::HamHam; mixing_time::Float64)`.

**Detection:**
- MethodError when trying to access `hamiltonian.beta` (field does not exist).

**Which phase:** API design phase.

---

### MOD-07: w0 and t0 Are Fourier Duals -- Not Free Parameters

**What goes wrong:**
The codebase enforces `t0 * w0 = 2*pi / 2^num_energy_bits` (see `validate_config!` in misc_tools.jl lines 192-195). If the cost counter accepts `r`, `w0`, AND `t0` as independent inputs, inconsistent values could be passed.

**Prevention:**
- Accept `r` and `w0`, compute `t0 = 2*pi / (2^r * w0)`.
- Or accept a Config (which already validates this relation).
- Add an assertion: `@assert isapprox(t0 * w0, 2*pi / 2^r)`.

**Detection:**
- Cost values that don't match hand calculations -- check the `t0 * w0` relation first.

**Which phase:** First implementation phase.

---

## Minor Pitfalls

Issues that cause inconvenience or confusion but not wrong numbers.

---

### MIN-01: Cost Struct Should Clearly Separate Per-Step from Total

**What goes wrong:**
Users need costs at multiple levels:
- Per OFT evaluation (one omega value, one direction)
- Per step (one application of the channel, all omegas, both directions)
- Per query (full mixing time evolution)

Without clear naming, users compare apples to oranges.

**Prevention:**
- Use a hierarchical result struct with self-documenting field names:
  ```julia
  struct HamSimCost
      # Per-step breakdown
      oft_time_per_energy::Float64     # T_QPE = 2 * sum_k |t_k| for one omega
      n_effective_energies::Int        # number of truncated energy grid points
      oft_time_per_step::Float64       # = oft_time_per_energy * n_effective_energies
      b_time_per_step::Union{Nothing, Float64}
      total_time_per_step::Float64     # = 2*oft_per_step + b_per_step
      # Total
      num_steps::Int
      total_ham_sim_time::Float64
      # Metadata
      mixing_time_used::Float64
      delta_used::Float64
      num_energy_bits::Int
      t0::Float64
      w0::Float64
      rescaling_factor::Float64
      construction::Symbol             # :GNS or :KMS
  end
  ```
- Provide a `Base.show` method for human-readable output.

**Which phase:** Struct definition phase.

---

### MIN-02: Trotter Domain Adds Trotter Step Counting

**What goes wrong:**
For TrotterDomain, each `e^{iHt_k}` requires `num_trotter_steps_per_t0 * |k|` Trotter steps. Users on TrotterDomain expect circuit-level cost, not just Hamiltonian simulation time.

**Prevention:**
- For TrotterDomain, report both `total_ham_sim_time` and `total_trotter_steps` as separate fields.
- `total_trotter_steps = num_trotter_steps_per_t0 * sum_k |k|` for the QPE grid.

**Which phase:** Can be deferred if TrotterDomain is not the initial focus.

---

### MIN-03: Cost Counter Must Be Deterministic

**What goes wrong:**
The cost counter is an analytical formula. It should be fully deterministic with no RNG dependency. If accidentally coupled to `run_thermalize` (e.g., to get mixing time), it inherits stochastic behavior.

**Prevention:**
- Pure function of config parameters and Hamiltonian properties.
- Input: Config (or explicit grid params), HamHam (for rescaling_factor), mixing_time.
- Output: HamSimCost struct.
- Test: two identical calls give bitwise identical results.

**Which phase:** All phases. Design principle.

---

### MIN-04: B Term Computation Requires _compute_b_plus Variant Dispatch

**What goes wrong:**
The B cost depends on which `b_plus` function is used (Gaussian vs Metropolis vs Smooth Metropolis). The dispatch logic in `_select_b_plus_calculator` (furnace_utensils.jl lines 143-156) checks `with_linear_combination`, then `a` and `b` values. The cost counter must replicate this dispatch exactly to get the correct truncated dictionary size.

**Prevention:**
- Reuse `_select_b_plus_calculator` directly (it returns the function + args).
- Then call `_compute_truncated_func` to get the truncated dictionary.
- This ensures the cost counter uses the same truncation as the simulator.

**Detection:**
- B cost that does not change when switching transition weight types is suspicious.

**Which phase:** B cost computation phase.

---

## Phase-Specific Warning Summary

| Phase | Pitfall | Severity | Key Mitigation |
|-------|---------|----------|----------------|
| API/struct design | CRIT-05: Classical vs quantum mixing time | Critical | Accept mixing_time as explicit input |
| API/struct design | MOD-02: Rescaling factor | Moderate | Store in struct; document units |
| API/struct design | MOD-03: Per-jump vs summed | Moderate | Per-step = one jump's cost |
| API/struct design | MOD-04: GNS vs KMS dispatch | Moderate | Dispatch on construction; Nothing for absent B |
| API/struct design | MOD-05: Domain independence | Moderate | Accept grid params explicitly |
| API/struct design | MOD-06: beta not on HamHam | Moderate | Require beta as argument |
| OFT cost | CRIT-01: Full vs truncated grid | Critical | Full QPE grid for inner time sum; truncated for outer energy count |
| OFT cost | CRIT-02: No filter weighting, factor of 2 | Critical | Pure grid formula; no sigma dependency in time sum |
| OFT cost | CRIT-06: Transition type energy range | Critical | Use _truncate_energy_labels for outer count |
| B cost | CRIT-03: Double sum structure | Critical | Explicit double loop; reuse truncated b-dicts |
| B cost | MIN-04: b_plus variant dispatch | Minor | Reuse _select_b_plus_calculator |
| Assembly | CRIT-04: Factor of 2 for U+cU^dagger | Critical | 2*OFT + 1*B, NOT 2*(OFT+B) |
| Assembly | MOD-01: Delta step count | Moderate | ceil(mixing_time/delta); linear 1/delta scaling |
| Assembly | MOD-07: w0/t0 consistency | Moderate | Compute one from other; validate |
| Reporting | MIN-01: Clear struct | Minor | Separate per-step, per-query, metadata |
| Reporting | MIN-02: Trotter steps | Minor | Report separately for TrotterDomain |
| Reporting | MIN-03: Determinism | Minor | Pure function, no RNG |

---

## Recommended Phase Ordering Based on Pitfall Dependencies

1. **API and struct design** first: Define `HamSimCost` struct, decide units convention (MOD-02), accept explicit mixing time (CRIT-05), handle GNS/KMS dispatch (MOD-04), domain-independent API (MOD-05). This is the foundation that all computation phases build on.

2. **OFT cost computation** second: Implement the full-grid `sum(|t_k|)` formula (CRIT-01, CRIT-02). Count effective energies via `_truncate_energy_labels` (CRIT-06). Include the factor-of-2 for U+controlled-U^dagger (CRIT-04). Validate against closed-form `t0 * N^2 / 2`. This is the simplest computation and typically the largest cost component.

3. **B coherent term cost** third: Implement the double-sum formula using actual truncated b-dicts (CRIT-03). Handle the b_plus variant dispatch (MIN-04). Only applies to KMS. Validate by comparing truncated dict sizes against expectations.

4. **Assembly and total cost** fourth: Combine OFT + B costs with step count (MOD-01). Apply correct factors (2*OFT + B, CRIT-04). Multiply by `ceil(mixing_time / delta)`.

5. **Trotter-specific and domain-specific reporting** last (optional): Add Trotter step counting for TrotterDomain users (MIN-02). This is an enhancement.

---

## Sources

- Direct codebase analysis of QuantumFurnace.jl (HIGH confidence)
  - `src/furnace_utensils.jl`: `_precompute_labels`, `_truncate_energy_labels`, `_create_energy_labels`, `_select_b_plus_calculator`, `_build_cptp_channel`
  - `src/nufft.jl`: `_prepare_oft_nufft_prefactors` (filter weights, base_weights)
  - `src/coherent.jl`: `B_time`, `B_trotter`, `_compute_b_minus`, `_compute_b_plus`, `_compute_b_plus_metro`, `_compute_b_plus_smooth`
  - `src/ofts.jl`: `oft!`, `time_oft!`, `trotter_oft!`
  - `src/energy_domain.jl`: `pick_transition`, `_pick_transition_kms`, `_pick_transition_gns`
  - `src/time_domain.jl`: `_truncate_time_labels_for_oft`
  - `src/structs.jl`: `Config`, `Workspace`, domain and construction types
  - `src/hamiltonian.jl`: `HamHam`, `_rescaling_and_shift_factors`
  - `src/jump_workers.jl`: `_jump_contribution!`
  - `src/furnace.jl`: `run_thermalize` (jump selection), `construct_lindbladian`
  - `src/misc_tools.jl`: `validate_config!` (Fourier relation enforcement)
- [Chen et al. (2023), "An efficient and exact noncommutative quantum Gibbs sampler", arXiv:2311.09207](https://arxiv.org/abs/2311.09207) (HIGH confidence -- foundational paper for GNS construction)
- [Chen, Kastoryano, Gilyen (2025), "Efficient quantum Gibbs samplers with KMS detailed balance condition", Comm. Math. Phys.](https://link.springer.com/article/10.1007/s00220-025-05235-3) (HIGH confidence -- KMS construction paper with B term)
- Project MEMORY.md: floor analysis (`floor = k_energy * delta + floor_Trotter`), Bohr frequency collision analysis, bi-exponential fitting diagnostics (HIGH confidence -- direct prior work)
