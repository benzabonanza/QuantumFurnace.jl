# Phase 29: Eigensolver Integration - Research

**Researched:** 2026-02-24
**Domain:** KrylovKit eigsolve wrapping for matrix-free Lindbladian spectral gap computation
**Confidence:** HIGH

## Summary

Phase 29 wraps KrylovKit's `eigsolve` around the matrix-free matvec infrastructure built in Phases 27-28 to provide a single-call `krylov_spectral_gap()` API for computing Lindbladian spectral gaps without constructing the full dense Liouvillian. The codebase already has all matvec functions (`apply_lindbladian!`, `apply_adjoint_lindbladian!`) for all four domains (Energy, Time, Trotter, Bohr) and the `KrylovWorkspace` construction. What remains is: (1) adding KrylovKit as a dependency, (2) implementing `apply_delta_channel!` for the CPTP channel path, (3) wrapping `eigsolve` with config-type dispatch, retry logic, and memory guards, (4) defining the `KrylovGapResult` struct, and (5) implementing the channel-to-Lindbladian eigenvalue conversion.

KrylovKit v0.10.2 provides `eigsolve(f, x0, howmany, which; kwargs...)` where `f` is any callable representing the linear map. This is precisely what we need: `f` will be the closure `v -> vec(apply_lindbladian!(ws, reshape(v, dim, dim), config, ham))`. The Arnoldi algorithm handles the general (non-Hermitian) Lindbladian, and `:LR` targeting finds eigenvalues with largest real part (steady state at Re(lambda)~0). For the CPTP channel path, the channel E(rho) = rho + delta*L(rho) has eigenvalues mu = 1 + delta*lambda_L, so `:LM` targeting (largest magnitude, near 1) finds the steady state and gap mode. The Lindbladian gap is recovered as `lambda_L = (mu - 1) / delta`.

**Primary recommendation:** Use KrylovKit Arnoldi with `:LR` for the Lindbladian path and `:LM` for the channel path. Implement a retry strategy that increases krylovdim by 50% per retry (30 -> 45 -> 68) for up to 3 retries. The channel-to-Lindbladian conversion is the exact linear formula `lambda_L = (mu_E - 1) / delta` (no log approximation needed since E = I + delta*L is exact by construction, not an exponential).

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

#### API surface design
- Single function `krylov_spectral_gap(config, ham)` -- config type dispatches the path
- `LindbladConfig` -> Lindbladian eigsolve with `:LR` targeting (largest real part)
- `ThermalizeConfig` -> CPTP channel eigsolve with `:LM` targeting (largest magnitude), using delta from config
- All four domains supported from day one (matvec already exists from Phase 28)
- Workspace allocated internally -- no user-facing workspace kwarg
- Conservative defaults: `krylovdim=30`, `howmany=4`, `tol=1e-10` -- user tunes down for speed
- KrylovKit kwargs passed through as keyword arguments for full control

#### Result struct contents
- `KrylovGapResult` stores top 2 eigenvalues and eigenvectors (fixed-point + gap mode)
- Eigenvectors reshaped to dim x dim matrices (not raw vectors)
- Unified `gap` field -- always represents Lindbladian spectral gap regardless of path
- Channel path auto-converts channel eigenvalues back to Lindbladian gap internally
- Convergence info and matvec count included as metadata
- **RESEARCH FLAG:** Exact channel->Lindbladian gap conversion formula needs careful verification during planning (log vs linear approximation, complex eigenvalue handling)

#### CPTP channel strategy
- `apply_delta_channel!` computes E(rho) = (I + delta*L)(rho) using delta from ThermalizeConfig
- Delta used as-is from config -- no validation or auto-selection
- Channel matvec lives in same Krylov module alongside Lindbladian matvec
- `:LM` targeting for channel (eigenvalues cluster near 1)
- `:LR` targeting for Lindbladian (steady state at 0, gap from lambda_2)

#### Failure and guard behavior
- Pre-flight memory estimate: `krylovdim * 4^n * 16 * 1.5` bytes
- Threshold: 80% of `Sys.free_memory()` -- adapts to the machine
- Memory guard issues `@warn` and proceeds (not a hard error)
- Convergence retry strategy: Claude's discretion (research KrylovKit behavior to pick between doubling vs incremental krylovdim increase)
- Final failure: simple error message -- "KrylovKit failed to converge: n_converged/howmany eigenvalues after N retries"

### Claude's Discretion
- Convergence retry strategy (krylovdim increase amount and number of retries)
- Internal workspace allocation details
- KrylovKit algorithm selection (Arnoldi vs other iterative methods)
- Exact channel->Lindbladian gap conversion formula (after research)

### Deferred Ideas (OUT OF SCOPE)
None -- discussion stayed within phase scope
</user_constraints>

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| KrylovKit.jl | 0.10.2 | Matrix-free iterative eigensolver | De facto Julia standard for Krylov methods; supports function-based linear maps, Arnoldi/Lanczos, thick restarts, and AD |
| LinearAlgebra (stdlib) | -- | BLAS operations, matrix reshaping | Already used throughout project |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| Base.Sys | stdlib | `Sys.free_memory()` for memory estimation | Pre-flight memory guard |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| KrylovKit | Arpack.jl | Already in Project.toml but requires matrix/LinearMap input, not bare functions. KrylovKit's function-based API is cleaner for the matvec pattern. |
| KrylovKit | ArnoldiMethod.jl | KrylovKit is more actively maintained (v0.10.2 Oct 2025) and has richer API (`:LR`, `:LM` targeting built in) |

**Installation:**
```julia
] add KrylovKit
```

Then add to `Project.toml` `[deps]` and `[compat]` sections. The `using KrylovKit` goes in `QuantumFurnace.jl`.

## Architecture Patterns

### Recommended File Structure
```
src/
├── krylov_workspace.jl      # [EXISTS] KrylovWorkspace struct + constructor
├── krylov_matvec.jl         # [EXISTS] apply_lindbladian!, apply_adjoint_lindbladian!
├── krylov_eigsolve.jl       # [NEW] krylov_spectral_gap(), apply_delta_channel!, KrylovGapResult
└── QuantumFurnace.jl        # Add: using KrylovKit, include krylov_eigsolve.jl, exports
```

### Pattern 1: Config-Type Dispatch for Eigsolve Path Selection
**What:** Julia multiple dispatch selects Lindbladian vs channel path based on config type.
**When to use:** When the same function name should have different behavior for different config types.
**Example:**
```julia
# Lindbladian path: config isa AbstractLiouvConfig
function krylov_spectral_gap(
    config::AbstractLiouvConfig,
    hamiltonian::HamHam,
    jumps::Vector{JumpOp};
    trotter=nothing,
    krylovdim=30, howmany=4, tol=1e-10,
    krylov_kwargs...
)
    ws = KrylovWorkspace(config, hamiltonian, jumps; trotter=trotter)
    dim = size(hamiltonian.data, 1)

    # Build matvec closure: Vector{ComplexF64} -> Vector{ComplexF64}
    function lindbladian_matvec(v::AbstractVector)
        rho = reshape(v, dim, dim)
        L_rho = apply_lindbladian!(ws, rho, config, hamiltonian)
        return vec(L_rho)
    end

    x0 = vec(Matrix{ComplexF64}(I(dim) / dim))  # maximally mixed initial guess
    vals, vecs, info = eigsolve(lindbladian_matvec, x0, howmany, :LR,
        Arnoldi(; krylovdim=krylovdim, tol=tol, verbosity=0); krylov_kwargs...)
    # ... post-process
end

# Channel path: config isa AbstractThermalizeConfig
function krylov_spectral_gap(
    config::AbstractThermalizeConfig,
    hamiltonian::HamHam,
    jumps::Vector{JumpOp};
    trotter=nothing,
    krylovdim=30, howmany=4, tol=1e-10,
    krylov_kwargs...
)
    # ... build channel_matvec using apply_delta_channel!
    # ... eigsolve with :LM targeting
    # ... convert channel eigenvalues to Lindbladian gap
end
```

### Pattern 2: KrylovKit Function-Based Linear Map
**What:** Pass a closure directly to KrylovKit eigsolve instead of building a matrix or LinearMap.
**When to use:** Matrix-free eigensolve where the linear map is defined by an in-place function.
**Example:**
```julia
# Source: KrylovKit.jl documentation (https://jutho.github.io/KrylovKit.jl/stable/man/eig/)
# eigsolve(f, x0, howmany, which; kwargs...)
vals, vecs, info = eigsolve(
    v -> vec(apply_lindbladian!(ws, reshape(v, dim, dim), config, ham)),
    rand(ComplexF64, dim^2),  # random initial vector
    4,    # howmany eigenvalues
    :LR,  # largest real part
    Arnoldi(krylovdim=30, tol=1e-10, maxiter=100, verbosity=0)
)
```
**Critical note:** The closure captures `ws` which is mutated in-place. KrylovKit calls the function sequentially (not in parallel), so this is safe. However, the input vector `v` and the output must be *different* objects -- KrylovKit does NOT support in-place mutation of the input. Since `apply_lindbladian!` writes to `ws.rho_out` (a separate buffer), and we `vec()` the result (creating a new vector), this is automatically satisfied.

### Pattern 3: Convergence Retry with krylovdim Increase
**What:** When `info.converged < howmany`, retry with larger Krylov subspace.
**When to use:** Always -- Krylov methods can fail to converge for difficult spectra.
**Example:**
```julia
function _eigsolve_with_retry(f, x0, howmany, which;
    krylovdim=30, tol=1e-10, maxiter=100, max_retries=3)

    current_krylovdim = krylovdim
    for attempt in 1:(max_retries + 1)
        vals, vecs, info = eigsolve(f, x0, howmany, which,
            Arnoldi(; krylovdim=current_krylovdim, tol=tol, maxiter=maxiter, verbosity=0))

        if info.converged >= howmany
            return vals, vecs, info
        end

        if attempt <= max_retries
            new_krylovdim = ceil(Int, current_krylovdim * 1.5)
            @warn "KrylovKit: $(info.converged)/$(howmany) converged. " *
                  "Retrying with krylovdim=$new_krylovdim (attempt $(attempt+1)/$(max_retries+1))"
            current_krylovdim = new_krylovdim
        end
    end

    error("KrylovKit failed to converge: $(info.converged)/$(howmany) eigenvalues " *
          "after $(max_retries + 1) attempts (final krylovdim=$(current_krylovdim))")
end
```

### Pattern 4: apply_delta_channel! Implementation
**What:** Compute E(rho) = rho + delta * L(rho) by reusing the existing Lindbladian matvec.
**When to use:** CPTP channel eigsolve path (ThermalizeConfig).
**Example:**
```julia
function apply_delta_channel!(
    ws::KrylovWorkspace{T},
    rho::Matrix{T},
    delta::Real,
    config_liouv::AbstractLiouvConfig,
    hamiltonian::HamHam,
) where {T<:Complex}
    # L(rho) -> ws.rho_out
    apply_lindbladian!(ws, rho, config_liouv, hamiltonian)
    # E(rho) = rho + delta * L(rho) -> ws.rho_out
    @. ws.rho_out = rho + delta * ws.rho_out
    return ws.rho_out
end
```
**Key insight:** `apply_delta_channel!` needs a `LiouvConfig` (not `ThermalizeConfig`) because the Lindbladian matvec dispatches on `AbstractLiouvConfig`. The `krylov_spectral_gap` method for `ThermalizeConfig` must extract a `LiouvConfig` from the `ThermalizeConfig` fields (they share the same Lindbladian parameters) or build one internally.

### Anti-Patterns to Avoid
- **Allocating new vectors inside the KrylovKit closure:** KrylovKit calls the function O(krylovdim * maxiter) times. Allocating a new `dim x dim` matrix each call would destroy performance. The existing `ws.rho_out` buffer handles this -- just `vec()` at the end (which is a view, not a copy, for contiguous arrays).
- **Using `vec(copy(ws.rho_out))` instead of `vec(ws.rho_out)`:** Actually, `vec()` on a `Matrix` returns a view (reshape). KrylovKit stores the result before calling the function again, so the view is safe. But verify this: if KrylovKit reuses the vector before storing, we'd need `copy`. **Correction:** KrylovKit does NOT reuse the input -- it stores each Krylov vector. But `vec()` returns a *reshaped view* of `ws.rho_out`, and `ws.rho_out` gets overwritten on the next call. So we MUST return `copy(vec(ws.rho_out))` or equivalently `vec(copy(ws.rho_out))` to avoid aliasing. **This is critical -- KrylovKit builds an orthonormal basis from the returned vectors, and if they all alias the same memory, the algorithm breaks completely.**
- **Forgetting to handle the jumps argument:** The CONTEXT.md specifies `krylov_spectral_gap(config, ham)` but the workspace needs `jumps`. The function must accept `jumps` as a required argument (likely `krylov_spectral_gap(config, hamiltonian, jumps; ...)`).

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Iterative eigensolver | Custom Arnoldi/Lanczos | KrylovKit.eigsolve | Thick restarts, convergence tracking, numerical stability, orthogonalization -- all handled |
| Memory estimation | Platform-specific `/proc/meminfo` parsing | `Sys.free_memory()` | Julia stdlib, cross-platform |
| Vector orthogonalization | Manual Gram-Schmidt in Krylov loop | KrylovKit's `ModifiedGramSchmidtIR` (default) | Numerically stable, with iterative refinement |

**Key insight:** KrylovKit handles the entire Krylov subspace management internally. Our code only needs to provide the linear map (matvec closure) and interpret the results. The retry logic wraps KrylovKit, not replaces it.

## Common Pitfalls

### Pitfall 1: Vector Aliasing in KrylovKit Closure
**What goes wrong:** KrylovKit stores each Krylov basis vector. If the closure returns a view of `ws.rho_out` (which gets overwritten each call), all stored vectors alias the same memory, producing garbage.
**Why it happens:** `vec(ws.rho_out)` returns a reshaped view, not a copy. The workspace reuses `rho_out` on every call.
**How to avoid:** Always return `copy(vec(ws.rho_out))` from the closure. This allocates one `dim^2`-element vector per KrylovKit call, which is unavoidable and expected.
**Warning signs:** All eigenvalues are identical, eigenvectors are all the same, `info.normres` is enormous.

### Pitfall 2: Config Type Mismatch for Channel Path
**What goes wrong:** `apply_lindbladian!` dispatches on `AbstractLiouvConfig`, but the channel path receives an `AbstractThermalizeConfig`. Calling `apply_lindbladian!(ws, rho, thermalize_config, ham)` will fail with a MethodError.
**Why it happens:** `ThermalizeConfig` is not a subtype of `AbstractLiouvConfig`. The type hierarchy is: `AbstractConfig > AbstractLiouvConfig > LiouvConfig/LiouvConfigGNS` and `AbstractConfig > AbstractThermalizeConfig > ThermalizeConfig/ThermalizeConfigGNS`.
**How to avoid:** The channel path must construct a corresponding `LiouvConfig` from the `ThermalizeConfig` fields. Both share the same Lindbladian parameters (num_qubits, beta, sigma, domain, etc.). Build a `LiouvConfig` (or `LiouvConfigGNS` for GNS configs) internally, stripping the `mixing_time` and `delta` fields. Then construct the `KrylovWorkspace` from this `LiouvConfig`.
**Warning signs:** `MethodError: no method matching apply_lindbladian!(::KrylovWorkspace, ::Matrix, ::ThermalizeConfig, ::HamHam)`.

### Pitfall 3: Wrong Eigenvalue Sorting for Lindbladian
**What goes wrong:** Using `:LM` (largest magnitude) for the Lindbladian misses the spectral gap. The steady state eigenvalue is near 0 (smallest magnitude), and `:LM` would find the most negative eigenvalues instead.
**Why it happens:** Lindbladian eigenvalues have Re(lambda) <= 0. The steady state is at lambda ~ 0, and the gap is from lambda_2 (next closest to 0). `:LR` (largest real part) correctly finds eigenvalues closest to 0 on the real axis.
**How to avoid:** Use `:LR` for Lindbladian, `:LM` for channel. This is already specified in the locked decisions.
**Warning signs:** Returned "spectral gap" is absurdly large (it's actually the spectral radius, not the gap).

### Pitfall 4: Channel Eigenvalue Conversion for Complex Eigenvalues
**What goes wrong:** The conversion formula `lambda_L = (mu_E - 1) / delta` produces complex Lindbladian eigenvalues, and taking `abs(real(lambda_L_2))` for the gap is correct but the imaginary part should not be ignored in the result struct.
**Why it happens:** Lindbladians can have complex eigenvalues (oscillatory modes). The channel E = I + delta*L inherits this: mu_E = 1 + delta*lambda_L is also complex.
**How to avoid:** Store the full complex eigenvalue in the result. The `gap` field should be `abs(real(lambda_L_2))` (following the convention in `diagnostics.jl:extract_leading_eigendata`). The full complex eigenvalue is available in the eigenvalues vector.
**Warning signs:** None obvious -- just ensure consistency with the existing dense diagnostics convention.

### Pitfall 5: KrylovWorkspace Input Mutation
**What goes wrong:** `apply_lindbladian!(ws, rho, config, ham)` reads from `rho` and writes to `ws.rho_out`. But KrylovKit may pass the *same* vector as input that it previously received as output (if using in-place operations). Since we reshape the input vector, the `rho` matrix could alias `ws.rho_out` if KrylovKit reuses memory.
**Why it happens:** KrylovKit manages its own vector storage. The input to our closure is a vector from KrylovKit's Krylov basis, not our workspace.
**How to avoid:** Since KrylovKit passes fresh vectors (from its orthonormal basis), and we reshape them (creating a view into KrylovKit's memory, not `ws.rho_out`), there is no aliasing issue on the *input* side. The aliasing concern is only on the *output* side (Pitfall 1). Still, verify this during implementation by checking that `pointer(rho) != pointer(ws.rho_out)`.
**Warning signs:** NaN/Inf in eigenvalues after a few iterations.

### Pitfall 6: krylovdim Must Exceed howmany
**What goes wrong:** If `krylovdim <= howmany`, KrylovKit cannot build a large enough subspace to extract the requested eigenvalues.
**Why it happens:** The Krylov subspace of dimension `krylovdim` is used to approximate `howmany` eigenvalues. The subspace must be strictly larger.
**How to avoid:** The defaults (krylovdim=30, howmany=4) satisfy this. Add a guard: `krylovdim > howmany || error("krylovdim must be > howmany")`.
**Warning signs:** KrylovKit error or immediate convergence failure.

## Code Examples

### KrylovKit eigsolve with Function Input
```julia
# Source: https://jutho.github.io/KrylovKit.jl/stable/man/eig/
using KrylovKit

# Define linear map as a function
f(v) = A * v  # or any callable

# Compute 4 eigenvalues with largest real part
vals, vecs, info = eigsolve(f, rand(ComplexF64, n), 4, :LR,
    Arnoldi(krylovdim=30, tol=1e-10, maxiter=100, verbosity=0))

# Check convergence
@assert info.converged >= 4 "Only $(info.converged) eigenvalues converged"

# info fields: converged, residual, normres, numops, numiter
```

### ConvergenceInfo Fields
```julia
# Source: https://jutho.github.io/KrylovKit.jl/stable/man/eig/
# info.converged::Int      -- number of converged eigenvalues
# info.residual::Vector    -- residual vectors
# info.normres::Vector     -- norms of residuals
# info.numops::Int         -- number of linear map applications (matvec count)
# info.numiter::Int        -- number of Krylov subspace rebuilds (restarts)
```

### Eigenvalue Sorting Values
```julia
# Source: https://jutho.github.io/KrylovKit.jl/stable/man/eig/
# :LM  -- largest magnitude (default)
# :LR  -- largest real part
# :SR  -- smallest real part
# :LI  -- largest imaginary part
# :SI  -- smallest imaginary part
# EigSorter(by; rev=false) -- custom sorting
```

### Pre-flight Memory Estimation
```julia
# Memory for Krylov basis vectors: krylovdim vectors of length dim^2
# Each ComplexF64 is 16 bytes
# Factor 1.5 for internal KrylovKit overhead (Hessenberg matrix, temporaries)
n_qubits = config.num_qubits
estimated_bytes = krylovdim * (4^n_qubits) * 16 * 1.5
available = Sys.free_memory()
if estimated_bytes > 0.8 * available
    @warn "Krylov memory estimate $(round(estimated_bytes / 1e9; digits=2)) GB " *
          "exceeds 80% of free memory $(round(available / 1e9; digits=2)) GB. " *
          "Consider reducing krylovdim or num_qubits."
end
```

### ThermalizeConfig to LiouvConfig Conversion
```julia
# ThermalizeConfig and LiouvConfig share all Lindbladian parameters.
# Build a LiouvConfig (or LiouvConfigGNS) from a ThermalizeConfig:
function _thermalize_to_liouv_config(tc::ThermalizeConfig)
    LiouvConfig(
        num_qubits = tc.num_qubits,
        with_coherent = tc.with_coherent,
        with_linear_combination = tc.with_linear_combination,
        domain = tc.domain,
        beta = tc.beta,
        sigma = tc.sigma,
        gaussian_parameters = tc.gaussian_parameters,
        a = tc.a, b = tc.b,
        num_energy_bits = tc.num_energy_bits,
        t0 = tc.t0, w0 = tc.w0,
        eta = tc.eta,
        num_trotter_steps_per_t0 = tc.num_trotter_steps_per_t0,
    )
end

function _thermalize_to_liouv_config(tc::ThermalizeConfigGNS)
    LiouvConfigGNS(
        num_qubits = tc.num_qubits,
        with_coherent = false,
        with_linear_combination = tc.with_linear_combination,
        domain = tc.domain,
        beta = tc.beta,
        sigma = tc.sigma,
        gaussian_parameters = tc.gaussian_parameters,
        a = tc.a, b = tc.b,
        num_energy_bits = tc.num_energy_bits,
        t0 = tc.t0, w0 = tc.w0,
        eta = tc.eta,
        num_trotter_steps_per_t0 = tc.num_trotter_steps_per_t0,
    )
end
```

## Channel-to-Lindbladian Gap Conversion (Research Flag Resolution)

### Mathematical Analysis

The CPTP channel is defined as `E(rho) = rho + delta * L(rho)`, i.e., `E = I + delta * L` as a superoperator.

**Eigenvalue relationship:** If `L * v = lambda_L * v`, then `E * v = (I + delta * L) * v = v + delta * lambda_L * v = (1 + delta * lambda_L) * v`. So:

```
mu_E = 1 + delta * lambda_L
lambda_L = (mu_E - 1) / delta
```

This is an **exact linear relationship** (not an approximation). There is no log involved because we define E = I + delta*L directly, NOT as E = exp(delta*L). The exponential relationship `E = exp(delta*L)` would give `lambda_L = log(mu_E) / delta`, but that is NOT what this codebase uses.

**Verification against existing code:** In `jump_workers.jl`, the thermalization channel `_jump_contribution!` for `ThermalizeConfig` builds the CPTP map using `K0 = I - alpha*R` and Kraus operators with `delta`-scaled jump terms. The structure is `E(rho) = K0*rho*K0' + delta * sum(L_a*rho*L_a') + U_res*rho*U_res'`, which to first order in delta gives `E(rho) ~ rho + delta*L(rho) + O(delta^2)`. So the `apply_delta_channel!` definition `E(rho) = rho + delta*L(rho)` is the **first-order Euler approximation** of the physical channel, and the eigenvalue conversion `lambda_L = (mu_E - 1) / delta` is exact for this linear approximation.

**Complex eigenvalue handling:** Both `lambda_L` and `mu_E` are generally complex. The conversion preserves the full complex structure. The spectral gap is `abs(real(lambda_L_2))`, following the convention in `diagnostics.jl:186`.

**Confidence:** HIGH -- the formula is a direct algebraic identity from the definition E = I + delta*L. No numerical approximation involved.

### Recommendation

Use the exact linear formula: `lambda_L = (mu_E - 1) / delta`. Sort channel eigenvalues by magnitude (`:LM` targeting finds largest |mu_E|, which corresponds to mu ~ 1 for the steady state). After conversion, sort by `abs(real(lambda_L))` to identify the gap mode (second smallest |Re|).

## Convergence Retry Strategy (Discretion Area Resolution)

### Analysis

KrylovKit's Arnoldi uses thick restarts (Krylov-Schur). The key parameters affecting convergence are:

1. **`krylovdim`** (default 30): Maximum Krylov subspace dimension before restart. Larger = more accurate approximation per cycle, but more memory.
2. **`maxiter`** (default 100): Number of restart cycles. More cycles = more chances to converge, but more compute.
3. **`tol`** (default 1e-12): Convergence tolerance on residual norm.

When convergence fails (`info.converged < howmany`), the most effective remedy is increasing `krylovdim` because it directly improves the quality of the eigenvalue approximation from the projected Hessenberg matrix. Increasing `maxiter` only helps if the algorithm was making progress but ran out of iterations.

### Recommendation: 50% krylovdim Increase, 3 Retries

**Strategy:** On each retry, increase `krylovdim` by 50% (multiply by 1.5, round up).

**Sequence:** 30 -> 45 -> 68 -> 102

**Rationale:**
- 50% increase is substantial enough to make a difference (unlike +5 or +10 which may not help)
- Not as aggressive as doubling (which wastes memory if a small increase would suffice)
- 3 retries gives krylovdim=102, which is 3.4x the original -- enough for most practical cases
- Memory grows linearly with krylovdim, so 3.4x is manageable

**maxiter:** Keep at 100 (default) for all attempts. If 100 restarts with krylovdim=102 does not converge, the problem is genuinely hard and the user should know.

### Alternative Considered: Fixed Step Increase
Adding a fixed amount (+20 per retry: 30 -> 50 -> 70 -> 90) gives a more predictable sequence but does not scale well -- for small problems krylovdim=30 is already large relative to dim^2, while for large problems +20 may be insufficient. The 50% multiplicative strategy scales naturally.

## Algorithm Selection (Discretion Area Resolution)

### Recommendation: Arnoldi

The Lindbladian is a general (non-Hermitian, non-normal) linear operator. Lanczos requires Hermiticity, which the Lindbladian does not satisfy (except in the KMS-similarity-transformed frame, which is not what we compute). Therefore, **Arnoldi is the only valid choice**.

KrylovKit's Arnoldi implementation uses the Krylov-Schur variant with thick restarts, which is numerically stable and efficient for non-Hermitian problems.

**Do NOT use Lanczos** even though the KMS-transformed Lindbladian is approximately Hermitian (see `diagnostics.jl:compute_anti_hermitian_defect`). The matvec computes `L(rho)`, not `D^{-1}*L*D(rho)`, so Arnoldi is required.

## KrylovGapResult Struct Design

### Recommended Fields

```julia
struct KrylovGapResult{T<:AbstractFloat}
    # Core spectral data
    eigenvalues::Vector{Complex{T}}     # Top eigenvalues (howmany), sorted by |Re(lambda)|
    spectral_gap::T                     # abs(real(eigenvalues[2]))
    fixed_point::Matrix{Complex{T}}     # Steady-state density matrix (from eigvec 1)
    gap_mode::Matrix{Complex{T}}        # Gap mode operator (from eigvec 2)

    # Metadata
    converged::Int                      # Number of converged eigenvalues
    matvec_count::Int                   # info.numops
    num_restarts::Int                   # info.numiter
    normres::Vector{T}                  # Residual norms for each eigenvalue

    # Channel-path info (nothing for Lindbladian path)
    channel_eigenvalues::Union{Nothing, Vector{Complex{T}}}  # Raw channel eigenvalues before conversion
    delta_used::Union{Nothing, T}       # Delta from ThermalizeConfig (nothing for Lindbladian path)
end
```

**Notes:**
- `eigenvalues` stores Lindbladian eigenvalues regardless of path (channel eigenvalues are converted internally)
- `fixed_point` is the eigenvector corresponding to eigenvalues[1], reshaped to dim x dim, Hermitianized, and trace-normalized (mirroring `furnace.jl:25-28` and `diagnostics.jl:209-218`)
- `gap_mode` is the eigenvector corresponding to eigenvalues[2], reshaped to dim x dim (not normalized as a density matrix)
- `channel_eigenvalues` preserves the raw channel eigenvalues for diagnostics (allows user to verify the conversion)

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Dense `eigen(L)` for all sizes | Dense for n<=6, Krylov for n>=7 | Phase 29 (this) | Enables spectral gap computation for 7+ qubits |
| Arpack `eigs(L, nev=2)` on dense matrix | KrylovKit `eigsolve(f, x0, 4, :LR)` on matvec closure | Phase 29 (this) | Matrix-free: O(dim^2) memory instead of O(dim^4) |
| No CPTP channel eigsolve | Channel eigsolve with `:LM` targeting | Phase 29 (this) | Direct comparison of Lindbladian and channel spectral properties |

**Deprecated/outdated:**
- `run_lindbladian` in `furnace.jl` uses `Arpack.eigs` on the full dense Liouvillian. This remains for backward compatibility but Phase 29's `krylov_spectral_gap` is the matrix-free replacement.
- The commented-out `linearmaps_liouv.jl` was an earlier attempt at matrix-free eigsolve using `LinearMaps.jl` + `Arpack`. It was abandoned as "sadly, slow" (line 16). KrylovKit's direct function-based API avoids the LinearMap wrapper overhead.

## Open Questions

1. **vec() Copy Requirement Verification**
   - What we know: `vec(Matrix)` returns a reshaped view (not a copy) in Julia. KrylovKit stores each Krylov basis vector. If the view aliases `ws.rho_out`, subsequent calls overwrite previously stored vectors.
   - What's unclear: Whether Julia's `vec()` on a Matrix always returns a view (it should for contiguous column-major storage, which `Matrix` guarantees).
   - Recommendation: Use `copy(vec(ws.rho_out))` in the closure to be safe. The allocation (one `dim^2` vector per KrylovKit call) is unavoidable and dwarfed by the matvec compute cost for n >= 4. Verify with `@allocated` that only the expected allocation occurs.

2. **Initial Vector Choice**
   - What we know: KrylovKit accepts any starting vector `x0`. A random vector is the default, but a better guess can improve convergence.
   - What's unclear: Whether `I/dim` (maximally mixed) or a random density matrix provides faster convergence for Lindbladian eigsolve.
   - Recommendation: Use `vec(Matrix{ComplexF64}(I(dim) / dim))` as default (maximally mixed state has overlap with the steady state). For the channel path, the same choice works since the channel's dominant eigenvalue (mu~1) corresponds to the steady state.

3. **jumps Argument Placement**
   - What we know: The CONTEXT.md says `krylov_spectral_gap(config, ham)` but `KrylovWorkspace` requires `jumps`. The user always has jumps available at the call site.
   - What's unclear: Whether the user intent is `krylov_spectral_gap(config, ham)` with jumps somehow inferred, or `krylov_spectral_gap(config, ham, jumps)`.
   - Recommendation: Use `krylov_spectral_gap(config, hamiltonian, jumps; trotter=nothing, kwargs...)` -- jumps is essential and cannot be inferred. The CONTEXT.md shorthand likely omitted it for brevity.

## Sources

### Primary (HIGH confidence)
- [KrylovKit.jl eigsolve docs](https://jutho.github.io/KrylovKit.jl/stable/man/eig/) -- full API reference for eigsolve, ConvergenceInfo, eigenvalue targeting
- [KrylovKit.jl algorithms docs](https://jutho.github.io/KrylovKit.jl/stable/man/algorithms/) -- Arnoldi/Lanczos constructor parameters, default values
- [KrylovKit.jl GitHub](https://github.com/Jutho/KrylovKit.jl) -- version 0.10.2, October 2025
- Codebase: `src/krylov_matvec.jl` -- existing apply_lindbladian! for all 4 domains
- Codebase: `src/krylov_workspace.jl` -- KrylovWorkspace struct and constructor
- Codebase: `src/structs.jl` -- Config type hierarchy (AbstractLiouvConfig vs AbstractThermalizeConfig)
- Codebase: `src/furnace.jl` -- existing run_lindbladian with Arpack eigs (reference implementation)
- Codebase: `src/diagnostics.jl` -- LindbladianResult, spectral gap convention (abs(real(eigenvalues[2])))

### Secondary (MEDIUM confidence)
- [Julia Discourse: Sys.free_memory](https://discourse.julialang.org/t/how-to-check-how-free-ram-is-available/15288) -- community confirmation of Sys.free_memory() availability

### Tertiary (LOW confidence)
- None

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- KrylovKit is the de facto Julia Krylov library, well-documented, v0.10.2 verified
- Architecture: HIGH -- patterns directly mirror existing codebase conventions (config dispatch, workspace pattern, matvec closures)
- Channel conversion formula: HIGH -- exact algebraic identity from E = I + delta*L
- Convergence retry strategy: MEDIUM -- 50% increase heuristic is reasonable but not validated on this specific problem class; may need tuning
- Pitfalls: HIGH -- aliasing issue (Pitfall 1) and config type mismatch (Pitfall 2) are verified from code inspection

**Research date:** 2026-02-24
**Valid until:** 2026-03-24 (KrylovKit API is stable; internal codebase patterns unlikely to change)
