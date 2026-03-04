# Feature Landscape: v2.2 Hamiltonian Simulation Time Counter

**Domain:** Analytical cost counting for quantum Gibbs sampler algorithm (Chen 2023 / Chen-Kastoryano-Gilyen 2024)
**Researched:** 2026-03-04
**Confidence:** HIGH (codebase analysis verified against paper abstracts and existing implementation)

## Physics Background: What We Are Counting

The quantum algorithm (Chen et al. arXiv:2311.09207, arXiv:2404.05998) prepares Gibbs states by simulating a Lindbladian on a quantum computer. The dominant cost is **Hamiltonian simulation queries** -- calls to `e^{+/- iHt}` for various times `t`. This milestone counts those queries analytically, without running any quantum or classical simulation.

The quantum algorithm has this per-step structure (for one delta step):
1. Pick a jump operator `A^a` (uniformly at random from `n_jumps` operators)
2. Apply `U` containing the OFT -- requires Hamiltonian simulation queries for `e^{iHt_k}` at each time grid point
3. Apply `controlled-U^dagger` -- same cost as step 2
4. If KMS construction: apply coherent unitary `e^{-i*delta*B}` -- B itself is built from a double time integral requiring its own Hamiltonian simulation queries

**Per-step Ham sim time = 2 * OFT_time + B_time** (for KMS)
**Per-step Ham sim time = 2 * OFT_time** (for GNS, which has no B)
**Total cost = (mixing_time / delta) * per_step_cost * n_jumps** (since each step picks one jump, and we run mixing_time/delta steps total)

### QPE Time Grid (the `r` Estimating Qubits)

The existing codebase already models this grid via `num_energy_bits` (called `r` in the papers):
- `N = 2^r` equally spaced points
- Energy labels: `w_k = w0 * k` for `k in -N/2 : N/2-1`
- Time labels: `t_k = t0 * k` for the same range (via `time_labels = energy_labels .* (t0 / w0)`)
- Fourier relation: `t0 * w0 = 2*pi / N`

The codebase function `_create_energy_labels(num_energy_bits, w0)` already generates this grid. The time grid is derived from `energy_labels .* (t0 / w0)`.

### OFT Hamiltonian Simulation Time

The OFT (Operator Fourier Transform) computes:
```
A(omega) = sum_k f(t_k) * e^{iH*t_k} * A * e^{-iH*t_k}
```
where `f(t_k) = exp(-t_k^2 * sigma^2) * exp(-i*omega*t_k)` (Gaussian-windowed Fourier prefactor).

On the quantum computer, this is implemented via LCU (linear combination of unitaries). Each term in the sum requires two Hamiltonian simulation queries (`e^{+iHt_k}` and `e^{-iHt_k}`), each for time `|t_k|`. But the LCU protocol uses these once (not per-omega), so the total OFT simulation time is:

```
OFT_time = sum_{k} |c_k| * |t_k|
```

where `|c_k|` is the weight (the LCU 1-norm contribution) at time point `t_k`. In the Gaussian case, `|c_k| = |f(t_k)| = exp(-t_k^2 * sigma^2)`.

However, the classical simulation truncates the time labels (via `_truncate_time_labels_for_oft`), whereas the quantum computer uses the FULL grid. The quantum algorithm's cost should use ALL `N = 2^r` time points before truncation.

**Key insight from codebase**: The function `_truncate_time_labels_for_oft(time_labels, sigma)` computes a cutoff at `sqrt(log(prefactor/tolerance)) / sigma`. The *classical* simulation truncates, but the *quantum* algorithm pays for ALL time points (the LCU protocol processes them all). For cost counting, we sum over the FULL grid, though in practice the Gaussian weights `exp(-t_k^2 * sigma^2)` make distant time points contribute negligibly.

### B Coherent Term Hamiltonian Simulation Time

The B operator (KMS construction only) has a double-time structure (from `coherent.jl`):

```
B = t0^2 * sum_t b_minus(t) * U_t^dag * [sum_s b_plus(s) * U_s * A^dag * U_s^{-2} * A * U_s] * U_t
```

where `U_t = diag(exp(i * eigvals * t))`. Each inner-sum `s` point requires Hamiltonian simulation for time `|s*beta|`, and each outer-sum `t` point requires time `|t/sigma|`.

**B_time = sum_t |b_minus(t)| * |t/sigma| * [sum_s |b_plus(s)| * |s*beta|]**

But the `b_minus` and `b_plus` are already truncated in the codebase via `_compute_truncated_func`, so only non-negligible terms contribute. The quantum algorithm must implement this double sum via Hamiltonian simulation queries.

### Total Algorithm Cost

```
steps = ceil(mixing_time / delta)
per_step_time = 2 * OFT_time + B_time    # KMS
per_step_time = 2 * OFT_time              # GNS
total_ham_sim_time = steps * per_step_time * n_jumps
```

Note: the `n_jumps` factor appears because each step randomly picks one jump, but over `steps` iterations all jumps are used. The expectation is that each jump is used `steps` times (not `steps / n_jumps` -- the random pick means each step pays for one jump's OFT, but the paper's total cost accounts for running all steps).

**Correction from milestone context**: The user specifies "per_step_cost = B_time + 2 * OFT_time" and "total = steps * per_step_cost". The `n_jumps` factor may or may not be included depending on whether we count per-jump or per-step. Both should be available.

## Table Stakes

Features users expect. Missing = the time counter is incomplete.

| Feature | Why Expected | Complexity | Dependencies | Notes |
|---------|--------------|------------|--------------|-------|
| `SimTimeBudget` result struct | Users need structured access to all cost components | Low | None | Immutable struct, follows `MixingTimeEstimate` pattern |
| OFT Ham sim time computation | This is the PRIMARY cost component, the whole point of the milestone | Med | QPE grid, transition weights, sigma | Sum over full QPE grid with Gaussian weights |
| B coherent term Ham sim time | Completes the full cost picture for KMS | Med | b_minus/b_plus functions, beta, sigma | Double-sum structure; reuse existing `_compute_b_minus`, `_compute_b_plus` formulas |
| Total step count | `ceil(mixing_time / delta)` | Low | mixing_time, delta | Trivial arithmetic |
| Per-step cost | `2 * OFT_time + B_time` (KMS) or `2 * OFT_time` (GNS) | Low | OFT_time, B_time | Simple formula |
| Total Ham sim time | `steps * per_step_cost` | Low | step count, per-step cost | The headline number |
| All 3 transition weight functions | Gaussian, Metropolis (a=0), smooth Metropolis (a>0, b>0) | Low | Already implemented in codebase | Affects B_time via `_compute_b_plus` variant selection |
| Accept `HamHam` for parameter extraction | Users should not have to manually extract `rescaling_factor`, `eigvals`, `nu_min` | Low | HamHam struct | Read `ham.rescaling_factor`, `ham.nu_min` directly |
| QPE grid parameter reporting | Users need `r`, `N`, `w0`, `t0`, energy range for paper tables | Low | Config or raw parameters | Pure arithmetic |
| KMS vs GNS construction toggle | GNS has no B term; cost structure differs | Low | `with_coherent` trait | Already exists as `with_coherent(::KMS) = true`, `with_coherent(::GNS) = false` |

## Differentiators

Features that set the time counter apart from naive cost estimates.

| Feature | Value Proposition | Complexity | Dependencies | Notes |
|---------|-------------------|------------|--------------|-------|
| Exact numerical OFT time (not asymptotic) | Papers give O-tilde bounds; we compute the EXACT sum for specific parameters | Med | Full grid enumeration | This is the key value-add over reading papers |
| Exact numerical B time (not asymptotic) | Same: exact double-sum for specific (n, beta, sigma, a, b) | Med | b_minus/b_plus evaluation | |
| `estimate_simulation_time(::MixingTimeEstimate, ...)` overload | Chains naturally with existing workflow: estimate mixing time, then compute cost | Low | `MixingTimeEstimate` struct | One-line wrapper that extracts `mixing_time` from estimate |
| Per-component breakdown in struct | Users see `oft_time`, `b_time`, `per_step_time`, `total_time`, `n_steps` separately | Low | Already separate fields | Enables analysis of which component dominates |
| Rescaled vs unrescaled time reporting | The Hamiltonian is rescaled by `rescaling_factor`; report both rescaled (algorithm-internal) and physical times | Low | `ham.rescaling_factor` | Physical time = rescaled_time * rescaling_factor |
| Sensitivity to `r` (num_energy_bits) | Show how cost scales with QPE resolution -- useful for resource estimation papers | Low | Sweep over `r` values | User can call in a loop; function is fast |
| Cost comparison across constructions | Given same (n, beta, sigma), compare KMS vs GNS cost | Low | Two calls | User composes, but struct makes comparison easy |
| `n_jumps` factoring | Report per-jump and per-step costs separately, so user can decide whether to include n_jumps in total | Low | n_jumps parameter | Per-jump OFT cost does not depend on n_jumps |

## Anti-Features

Features to explicitly NOT build.

| Anti-Feature | Why Avoid | What to Do Instead |
|--------------|-----------|-------------------|
| Config-based API | Config is a simulation struct with `mixing_time`, `delta`, `sim`, `domain` -- none relevant for cost counting; creates coupling to simulation infrastructure | Accept raw parameters: `r`, `w0`, `t0`, `sigma`, `beta`, `a`, `b`, `delta`, `mixing_time`, `n_jumps` |
| Domain-dispatched counting | The quantum algorithm cost is domain-independent. EnergyDomain vs TrotterDomain is a classical simulation distinction, not a quantum cost distinction | Single function with no domain dispatch |
| Running actual simulations | Cost counting must be purely analytical -- no matrix construction, no eigendecomposition, no NUFFT | Pure arithmetic on scalar parameters + grid sums |
| Modifying HamHam or Config structs | Adding fields to immutable structs breaks BSON serialization | New standalone struct |
| Plotting or visualization | Presentation is a separate concern | Return data struct; let scripts/notebooks plot |
| Asymptotic complexity analysis | The counter gives exact numerical costs, not big-O | Compute actual numerical values |
| Gate-level compilation | The "time" is Hamiltonian simulation time (sum of `|t_k|`), not gate count. Gate count depends on the specific Hamiltonian simulation method chosen (Trotter, QSP, etc.) | Report Ham sim time only; gate compilation is a separate concern |
| Trotter step counting for quantum algorithm | Trotter step count for the quantum algorithm's internal Hamiltonian simulation is a deeper concern. QuantumFurnace's `TrotterDomain` is about the *classical* simulation, not the quantum one | Only report Ham sim time; user converts to gates if needed |

## Feature Dependencies

```
SimTimeBudget struct                 (no deps -- pure data container)
    |
    +--- _qpe_grid(r, w0)           (no deps -- generates N=2^r energy/time labels)
    |
    +--- _oft_ham_sim_time(          (depends on: grid, sigma)
    |        time_labels, sigma)
    |        = sum_k exp(-t_k^2 * sigma^2) * |t_k|   [Gaussian weighting]
    |
    +--- _b_ham_sim_time(            (depends on: grid, beta, sigma, a, b, transition_type)
    |        time_labels, beta,
    |        sigma, a, b, ...)
    |        = double sum over b_minus, b_plus with time factors
    |
    +--- _transition_weight_b_plus(  (depends on: a, b, sigma, beta -- selects Metro/smooth/Gauss)
    |        ..., a, b, sigma, beta)
    |
    v
compute_simulation_time(             (depends on: all above + mixing_time, delta, n_jumps)
    r, w0, sigma, beta, a, b,
    delta, mixing_time, n_jumps;
    construction=KMS())
    |
    +--- HamHam convenience overload (depends on: compute_simulation_time + HamHam fields)
    |        Extracts: rescaling_factor, num_qubits
    |
    +--- MixingTimeEstimate overload (depends on: compute_simulation_time + MixingTimeEstimate fields)
             Extracts: mixing_time from estimate
```

## Detailed Formula Inventory

### 1. OFT Hamiltonian Simulation Time

The codebase's `time_oft!` in `src/ofts.jl` shows the prefactors:
```julia
@. caches.prefactors = exp(-time_labels^2 * sigma^2 - 1im * energy * time_labels)
```

The magnitude (LCU 1-norm contribution per time point) is:
```
|c_k| = exp(-t_k^2 * sigma^2)
```

This is energy-independent (the `exp(-i*omega*t_k)` has unit magnitude). The OFT cost per energy label is:
```
OFT_time_per_omega = sum_{k=0}^{N-1} exp(-t_k^2 * sigma^2) * |t_k|
```

But the quantum algorithm processes ALL omega values via QPE simultaneously, so:
```
OFT_time = sum_{k=0}^{N-1} exp(-t_k^2 * sigma^2) * |t_k|
```

This is the same regardless of omega -- the LCU processes the sum once.

### 2. B Coherent Term Hamiltonian Simulation Time

From `coherent.jl`, the B term has structure:
```
B = t0^2 * sum_t b_minus(t) * [conj(U_t) .* inner_sum .* U_t_row]
```
where `U_t = exp(i * eigvals * t/sigma)` and the inner sum is:
```
inner_sum = sum_s b_plus(s) * [U_s .* (A' * diag(U_s^{-2}) * A) .* U_s_row]
```
where `U_s = exp(i * eigvals * s*beta)`.

The quantum implementation requires Hamiltonian simulation for each time evolution:
- Outer loop: time `|t/sigma|` per outer-sum point with weight `|b_minus(t)|`
- Inner loop: time `|s*beta|` per inner-sum point with weight `|b_plus(s)|`

```
B_inner_time = sum_s |b_plus(s)| * |s * beta|     [for each A^a]
B_outer_time = sum_t |b_minus(t)| * |t / sigma|
B_time = B_outer_time + B_inner_time               [total per jump]
```

Note: the inner and outer sums involve DIFFERENT time evolutions (inner uses `s*beta`, outer uses `t/sigma`), so their costs add rather than multiply.

### 3. Transition Weight Impact

The transition weight function (gamma) affects the B term through the choice of `b_plus`:
- **Gaussian**: `_compute_b_plus(t, beta, w_gamma, sigma_gamma)` -- simple exponential decay
- **Metropolis** (a=0): `_compute_b_plus_metro(t, beta, sigma, eta)` -- has `1/t` singularity regulated by eta
- **Smooth Metropolis** (a>0, b>0): `_compute_b_plus_smooth(t, beta, sigma, a, b)` -- Glauber-like

The OFT cost is independent of transition weight (the Gaussian window `exp(-t^2*sigma^2)` is the same). Only the B term cost varies with weight function choice.

### 4. Total Cost Assembly

```
n_steps = ceil(Int, mixing_time / delta)
per_step = 2 * OFT_time + B_time            # KMS
per_step = 2 * OFT_time                     # GNS (no B)
total_ham_sim_time = n_steps * per_step
```

## MVP Recommendation

Prioritize (in implementation order):
1. **`SimTimeBudget` struct** -- the output container; needed by everything else
2. **`_qpe_grid` helper** -- generates the full N=2^r time/energy grid
3. **OFT time counting** -- `sum_k exp(-t_k^2 * sigma^2) * |t_k|`; this is the core computation
4. **B time counting** -- evaluates `b_minus`, `b_plus` on the grid and computes the double sum
5. **`compute_simulation_time` main API** -- assembles OFT + B + steps into total
6. **HamHam convenience overload** -- extracts parameters from HamHam for easy use
7. **MixingTimeEstimate chaining overload** -- syntactic sugar for the common workflow

Defer: Nothing. All features are small enough (~200-300 lines total including tests) to ship in one milestone. The above ordering is the dependency order for implementation.

## API Surface

### Primary Entry Point
```julia
compute_simulation_time(;
    r::Int,                    # QPE estimating qubits (num_energy_bits)
    w0::Float64,               # Energy grid spacing
    sigma::Float64,            # Gaussian width
    beta::Float64,             # Inverse temperature
    delta::Float64,            # Step size
    mixing_time::Float64,      # Total mixing time
    n_jumps::Int,              # Number of jump operators
    construction::AbstractConstruction = KMS(),  # KMS or GNS
    a::Float64 = 0.0,         # Transition weight parameter
    b::Float64 = 0.0,         # Transition weight parameter
    eta::Float64 = 0.0,       # Metropolis regularization
    # Gaussian-specific (only if a==0 && b==0 && not Metro)
    w_gamma::Union{Nothing,Float64} = nothing,
    sigma_gamma::Union{Nothing,Float64} = nothing,
) -> SimTimeBudget
```

### Convenience Overloads
```julia
# Extract params from HamHam + Config-like parameters
compute_simulation_time(ham::HamHam; r, sigma, beta, delta, mixing_time, n_jumps, ...)

# Chain from MixingTimeEstimate (extracts mixing_time)
compute_simulation_time(est::MixingTimeEstimate; r, w0, sigma, beta, delta, n_jumps, ...)
```

### Result Struct
```julia
struct SimTimeBudget
    # Per-step components
    oft_time::Float64           # sum_k exp(-t_k^2*sigma^2) * |t_k|
    b_time::Float64             # B coherent term time (0.0 for GNS)
    per_step_time::Float64      # 2 * oft_time + b_time

    # Total cost
    n_steps::Int                # ceil(mixing_time / delta)
    total_time::Float64         # n_steps * per_step_time

    # Grid info
    r::Int                      # QPE estimating qubits
    grid_points::Int            # N = 2^r
    w0::Float64                 # Energy spacing
    t0::Float64                 # Time spacing
    t_max::Float64              # Maximum |t_k| in grid

    # Parameters used
    sigma::Float64
    beta::Float64
    delta::Float64
    mixing_time::Float64
    n_jumps::Int
    construction::Symbol        # :KMS or :GNS
end
```

## Sources

- **Codebase analysis** (HIGH confidence):
  - `src/energy_domain.jl`: `_create_energy_labels(num_energy_bits, w0)` -- QPE grid generation
  - `src/ofts.jl`: `time_oft!` -- OFT prefactor formula `exp(-t^2*sigma^2 - i*omega*t)`
  - `src/coherent.jl`: `B_time`, `_compute_b_minus`, `_compute_b_plus` variants -- B term structure
  - `src/furnace_utensils.jl`: `_precompute_labels`, `oft_domain_prefactor` -- grid and normalization
  - `src/misc_tools.jl`: `validate_config!` -- Fourier relation `t0 * w0 = 2*pi / 2^N`
  - `src/mixing.jl`: `MixingTimeEstimate` -- pattern for result struct
  - `src/structs.jl`: `Config`, domain types, `with_coherent` trait
- **Paper abstracts** (MEDIUM confidence -- abstracts only, not full text):
  - [Chen et al. 2023 (arXiv:2311.09207)](https://arxiv.org/abs/2311.09207): "An efficient and exact noncommutative quantum Gibbs sampler" -- cost proportional to mixing_time * beta, polylog in precision
  - [Chen-Kastoryano-Gilyen 2024 (arXiv:2404.05998)](https://arxiv.org/abs/2404.05998): "Efficient quantum Gibbs samplers with KMS detailed balance condition" -- total Hamiltonian simulation time metric, truncation time T = O-tilde(beta * log(t_mix/epsilon))
  - [Chen et al. 2023 (arXiv:2303.18224)](https://arxiv.org/abs/2303.18224): "Quantum Thermal State Preparation" -- original GNS construction
- **Milestone context** (HIGH confidence): User-provided specification of v2.2 scope
