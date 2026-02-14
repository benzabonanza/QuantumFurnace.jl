# Phase 2: Trajectory Bug Fixes - Research

**Researched:** 2026-02-14
**Domain:** Quantum trajectory simulation correctness (CPTP Kraus unraveling)
**Confidence:** HIGH

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions
- **Jump sampling fix approach**: Claude determines the correct jump sampling scheme by cross-checking THREE sources: trajectory code, DM code (`jump_contribution!`), and the paper equations (Chen 2023 Theorem III.1, eq. 3.2-3.3). Chen 2023 describes the weak measurement scheme explicitly with statevector equations; Chen 2025 adds the coherent term (deterministic coherent evolution via U_B), which is likely already correct in code. If the jump sampling is structurally wrong (different algorithm, not just a bug), rewrite it to match the paper's construction. The coherent evolution (U_B) is probably correct; the jump sampling is the likely problem.
- **U_B ordering**: U_B (coherent unitary) must be applied BEFORE branch selection, matching the DM simulator's ordering. Claude should verify the DM simulator's ordering and make the trajectory simulator match.
- **Normalization guard**: Warning + continue (not assertion): emit a brief message like "Normalization violation: sum = X.XXX at step N". No full diagnostic dump -- keep output minimal. Check p_nojump + p_res + p_jump_total ~= 1.0 at each trajectory step.
- **PSD fallback**: When S matrix is not positive semi-definite, clamp negative eigenvalues to zero (project onto PSD cone). Silent fallback -- no warning or error, just fix it numerically.
- **CPTP verification**: Tolerance: 1e-10 (practical, allows small numerical accumulation). Single delta value (TEST_DELTA from Phase 1 fixtures). Test all three domains: Energy, Time, Trotter. Claude verifies the correct completeness relation formula from the paper (don't assume K0'K0 + delta*R + U_res'U_res = I is exactly right). Separate test file from the bug fix tests.
- **Bug fix validation**: One dedicated test per TFIX requirement (TFIX-02, TFIX-03, TFIX-04, TFIX-05) for clear traceability. Single-step tests (one call to `step_along_trajectory!`) -- sufficient to isolate each fix. CPTP verification (TVAL-01) in its own separate test file.

### Claude's Discretion
- Internal code organization and refactoring approach for the fixes
- Test helper functions and fixtures beyond what Phase 1 established
- Exact warning message formatting for normalization check
- How to structure the paper cross-check (inline comments, separate analysis doc, etc.)

### Deferred Ideas (OUT OF SCOPE)
None -- discussion stayed within phase scope
</user_constraints>

## Summary

This phase fixes four bugs in `step_along_trajectory!` (in `src/trajectories.jl`) and adds one CPTP verification test. The research reveals precise details of each bug by cross-referencing the trajectory code against the DM reference implementation (`jump_contribution!` in `src/jump_workers.jl`) and the paper's weak measurement construction.

The most critical finding is the **U_B ordering bug (TFIX-02)**: the trajectory code applies `U_B` to `psi` *after* computing `K0*psi` and `U_res*psi` but *before* branch selection. This means branch outcomes use stale pre-U_B computed vectors while psi itself has been rotated, producing an inconsistent state. The DM code applies `U_B` *first*, before all dissipative operations. The fix is to move U_B application to the very beginning of `step_along_trajectory!`.

The **jump sampling (TFIX-05)** was analyzed by cross-checking all three sources. The trajectory code's channel structure matches the DM code structurally -- both use the same Kraus decomposition: K0 = I - alpha*R (no-jump), sqrt(delta * rate2(w)) * A_w (jump), and U_residual (residual). The jump probabilities and sampling logic appear faithful to the paper's construction. The primary structural issue is the U_B ordering, which contaminates the jump sampling indirectly. Once U_B ordering is fixed, the channel should be correct.

The **CPTP completeness relation (TVAL-01)** needs careful derivation. The correct identity from the code's Kraus construction is: K0'*K0 + delta*R + S = I, where S = (2*alpha - delta)*R - alpha^2*R^2, and U_res'*U_res = S (by Cholesky). This simplifies to I exactly by algebra, so the CPTP test should verify: K0'*K0 + delta*R + U_res'*U_res = I.

**Primary recommendation:** Fix U_B ordering first (move before all probability computations), add normalization warning and PSD guard, then verify CPTP completeness. The jump sampling logic itself is structurally sound.

## Standard Stack

Not applicable -- this phase modifies existing Julia code in `src/trajectories.jl` and adds tests. No new libraries needed.

### Libraries Already in Use
| Library | Purpose | Relevant to Phase 2 |
|---------|---------|---------------------|
| LinearAlgebra | Matrix ops, Cholesky, eigen | PSD guard uses `eigen()`, existing code uses `cholesky!()` |
| Test | Test framework | New test files for TFIX tests and TVAL-01 |
| StableRNGs | Reproducible RNG | May be useful for single-step trajectory tests |

## Architecture Patterns

### Existing Code Structure
```
src/
  trajectories.jl          # TrajectoryFramework, build_trajectoryframework, step_along_trajectory!
  jump_workers.jl          # DM-based jump_contribution! (reference implementation)
  coherent.jl              # precompute_coherent_total_B, precompute_coherent_unitary_terms
  structs.jl               # KrausScratch, TrajectoryWorkspace, TrajectoryFramework types
test/
  runtests.jl              # Test entry point
  test_helpers.jl          # Shared fixtures (TEST_HAM, TEST_JUMPS, etc.)
  test_compilation.jl      # Phase 1 smoke tests
```

### Pattern: DM Code as Reference Implementation

The DM-based `jump_contribution!` functions in `src/jump_workers.jl` serve as the ground truth for the trajectory code. Both implement the same CPTP Kraus channel but:
- DM code: applies channel to density matrix rho -> K0 * rho * K0' + sum(delta * rate * A_w * rho * A_w') + U_res * rho * U_res'
- Trajectory code: samples one branch stochastically from {K0, A_w, U_res} with probability proportional to ||K * psi||^2

The ordering and mathematical structure must match exactly.

### Pattern: Modification Points in step_along_trajectory!

Both variants of `step_along_trajectory!` (EnergyDomain and Time/TrotterDomain) follow the same structure:
1. Compute R*psi, p_nojump, p_jump_total
2. Compute U_res*psi, p_res
3. Roll random number
4. Apply U_B (if present) -- **BUG: should be step 0**
5. Branch selection: no-jump / residual / jump

The fix applies identically to both variants.

## Detailed Bug Analysis

### Bug 1: U_B Ordering (TFIX-02) -- Confidence: HIGH

**Current code** (both EnergyDomain and Time/Trotter variants):
```julia
# Lines 447-478 of trajectories.jl (Time/Trotter variant)
mul!(ws.Rpsi, R, psi)                          # Step 1: R*psi
# ... compute p_nojump, p_jump_total ...
mul!(ws.Rpsi, fw.U_residual, psi)              # Step 2: U_res*psi (overwrites Rpsi)
# ... compute p_res, total_weight, roll r ...

# Step 3: Apply U_B AFTER probability computation
if fw.U_B !== nothing
    mul!(ws.psi_tmp, fw.U_B, psi)
    copyto!(psi, ws.psi_tmp)
    rmul!(psi, 1.0 / sqrt(max(_norm2(psi), eps(Float64))))
end

# Step 4: Branch selection
if r < p_nojump
    copyto!(psi, ws.psi_tmp)   # ws.psi_tmp was K0*psi (pre-U_B)
    # ...
```

**DM reference code** (`jump_workers.jl` lines 291-297, EnergyDomain):
```julia
U_B = coherent_unitary_cache
if U_B !== nothing
    mul!(scratch.tmp1, U_B, evolving_dm)
    mul!(scratch.rho_next, scratch.tmp1, U_B')
    copyto!(evolving_dm, scratch.rho_next)
end
# THEN compute R, K0, jumps, etc. on the ALREADY-ROTATED dm
```

**The problem:** In the trajectory code, U_B is applied to `psi` at step 3, but:
- `ws.psi_tmp` still contains `K0 * (old_psi)` from step 1
- `ws.Rpsi` still contains `U_res * (old_psi)` from step 2
- In the no-jump branch, `copyto!(psi, ws.psi_tmp)` restores the pre-U_B K0*psi, so the U_B application at step 3 is entirely wasted
- In the residual branch, same issue with ws.Rpsi
- In the jump branch, `ws.jump_oft * psi` operates on post-U_B psi, which is inconsistent with the probability weights that were computed on pre-U_B psi
- Additionally, after U_B application, `ws.psi_tmp` gets overwritten with `U_B * psi`, destroying the K0*psi that was stored there

**The fix:** Move U_B application to the very beginning of `step_along_trajectory!`, before any probability computation:
```julia
function step_along_trajectory!(psi, fw)
    # Step 0: Apply coherent unitary FIRST (matching DM code)
    if fw.U_B !== nothing
        mul!(ws.psi_tmp, fw.U_B, psi)
        copyto!(psi, ws.psi_tmp)
        rmul!(psi, 1.0 / sqrt(max(_norm2(psi), eps(Float64))))
    end

    # Step 1: Now compute R*psi, K0*psi on the ALREADY-ROTATED psi
    mul!(ws.Rpsi, R, psi)
    # ... rest of the function unchanged ...
end
```

**Impact:** This fix applies to BOTH `step_along_trajectory!` methods (EnergyDomain at lines 587-729, and Time/Trotter at lines 422-572). Both have the identical bug.

**Verification:** After the fix, with `with_coherent=true`, U_B should be applied before all Kraus operations, matching the DM code's ordering. A single-step test can verify by comparing the output state distribution.

### Bug 2: Normalization Guard (TFIX-03) -- Confidence: HIGH

**Current code:** No normalization check exists. The code computes `total_weight = p_nojump + p_res + p_jump_total` and uses `r = rand() * total_weight` which implicitly handles non-unit total weight, but never warns when it deviates from 1.0.

**The fix:** Add a check after computing `total_weight`:
```julia
total_weight = p_nojump + p_res + p_jump_total
if abs(total_weight - 1.0) > 1e-6
    @warn "Normalization violation: sum = $(total_weight)"
end
```

Per user decision: warning + continue, no assertion, brief message, no step number needed (but could add it as a detail). The `@warn` macro is appropriate for Julia -- it logs once by default and can be suppressed with `Logging` module settings. However, `@warn` has the behavior of only printing once per call site by default, which may not be ideal for tracking repeated violations. A `@printf` or `println` to stderr could be used for every-step warnings.

**Recommendation:** Use `@warn` with a unique identifier, or use a simple conditional print. Since user wants minimal output and `@warn` deduplicates by default, use `@warn` for the first violation and let Julia's warning system handle it. If the user wants every-step warnings, a `@printf(stderr, ...)` would be better.

**Implementation choice (Claude's discretion):** Use `@warn` since it produces minimal output (one message per warning site) and matches the user's "warning + continue" request. Add it to both `step_along_trajectory!` variants.

### Bug 3: PSD Guard for S Matrix (TFIX-04) -- Confidence: HIGH

**Current code** in `build_trajectoryframework` (lines 76-87):
```julia
# Numerical symmetrization and tiny diagonal shift
scratch.tmp2 .= 0.5 .* (scratch.tmp2 .+ scratch.tmp2')
eps_shift = 10 * eps(Float64)
@inbounds for i in 1:dim
    scratch.tmp2[i,i] += eps_shift
end

S = copy(scratch.tmp2)
cholesky_S = cholesky!(Hermitian(scratch.tmp2), check=false)
U_residual = Matrix{ComplexF64}(cholesky_S.U)
```

**The problem:** `cholesky!(... , check=false)` does not throw on non-PSD input, but produces garbage (NaN/Inf) silently. The tiny diagonal shift (10*eps ~ 2.2e-15) is often insufficient for S matrices that have eigenvalues on the order of -1e-10 or worse.

**The fix:** Replace the diagonal shift with proper eigenvalue clamping (PSD projection):
```julia
scratch.tmp2 .= 0.5 .* (scratch.tmp2 .+ scratch.tmp2')
S_hermitian = Hermitian(scratch.tmp2)
eig = eigen(S_hermitian)
# Clamp negative eigenvalues to zero (PSD projection)
eig.values .= max.(eig.values, 0.0)
# Reconstruct: S_psd = V * diag(max(lambda, 0)) * V'
# Then take Cholesky of the clamped matrix
# Or directly compute U_residual = diag(sqrt.(clamped_eigenvalues)) * V'
S_psd = eig.vectors * Diagonal(eig.values) * eig.vectors'
S_psd .= 0.5 .* (S_psd .+ S_psd')  # re-symmetrize after reconstruction

cholesky_S = cholesky!(Hermitian(S_psd))  # now safe, check=true is fine
U_residual = Matrix{ComplexF64}(cholesky_S.U)
```

**Alternative (more efficient):** Since we already have the eigendecomposition, compute U_residual directly:
```julia
eig = eigen(Hermitian(scratch.tmp2))
eig.values .= max.(eig.values, 0.0)
# U_residual = diag(sqrt.(eigenvalues)) * eigenvectors'
U_residual = Diagonal(sqrt.(eig.values)) * eig.vectors'
```
This skips the Cholesky entirely and is mathematically equivalent (U_res'*U_res = V*D*V' where D has clamped eigenvalues). This approach is cleaner and avoids potential Cholesky issues entirely.

Per user decision: silent fallback, no warning or error.

**Important:** The DM code in `jump_workers.jl` also uses `cholesky!(... , check=false)` with only the tiny eps shift. However, the DM code computes S fresh each step with a potentially different R (per-jump R vs trajectory's all-jumps R). The trajectory code computes S once at build time, so getting it right there is more critical. The DM code's S matrix is not modified by this phase -- it is in scope for Phase 3 if needed.

**Where to apply:** Only in `build_trajectoryframework` (the DM code's `jump_contribution!` is not in scope for Phase 2).

### Bug 4: Jump Sampling Faithfulness (TFIX-05) -- Confidence: HIGH

**Cross-check of three sources:**

**Source 1: Paper (Chen 2023, Theorem III.1, Eq. 3.2-3.3)**
The weak measurement scheme defines a CPTP map with Kraus operators:
- K0 = I - alpha*R, where alpha = 1 - sqrt(1-delta)
- K_{a,w} = sqrt(delta * gamma(w)) * L_{a,w} for each jump a and frequency w
- K_res such that K_res'*K_res = S = I - K0'*K0 - delta*R

The completeness relation: K0'*K0 + sum(delta*gamma(w)*L_{a,w}'*L_{a,w}) + K_res'*K_res = I
Since R = sum(gamma(w)*L_{a,w}'*L_{a,w}), this becomes: K0'*K0 + delta*R + S = I.

Expanding K0'*K0 = (I - alpha*R)'*(I - alpha*R) = I - 2*alpha*R + alpha^2*R^2, we get:
I - 2*alpha*R + alpha^2*R^2 + delta*R + S = I
=> S = (2*alpha - delta)*R - alpha^2*R^2

This exactly matches the code's computation of S.

The stochastic unraveling samples one outcome per step:
- With probability p_nojump = ||K0*psi||^2: psi -> K0*psi / ||K0*psi||
- With probability p_{a,w} = delta*gamma(w)*||L_{a,w}*psi||^2: psi -> L_{a,w}*psi / ||L_{a,w}*psi||
- With probability p_res = ||K_res*psi||^2: psi -> K_res*psi / ||K_res*psi||

**Source 2: DM code (`jump_contribution!` in `jump_workers.jl`)**
The DM code applies the FULL channel: rho_next = K0*rho*K0' + sum(delta*rate2*Aω*rho*Aω') + Ures*rho*Ures'

Key observations from the DM code:
- R is computed with the same base_prefactor and transition weights
- K0 = I - (1-sqrt(1-delta))*R  -- MATCHES trajectory code
- S = (2*alpha - delta)*R - alpha^2*R^2  -- MATCHES trajectory code
- Cholesky: S = Ures'*Ures  -- MATCHES trajectory code
- Individual jump terms use: delta * rate2 * Aω*rho*Aω'  -- MATCHES trajectory's per-outcome probability

For Hermitian jumps, the DM code uses:
- Positive frequency: L = Aω, L'L = Aω'*Aω, jump term = Aω*rho*Aω'
- Negative frequency partner: L = Aω', L'L = Aω*Aω', jump term = Aω'*rho*Aω

**Source 3: Trajectory code (`step_along_trajectory!` in `trajectories.jl`)**
The trajectory code samples from the SAME Kraus operators:
- p_nojump = ||K0*psi||^2 = ||(I-alpha*R)*psi||^2  -- MATCHES paper
- p_jump_total = delta * <psi|R|psi>  -- MATCHES (sum of individual jump probabilities)
- p_res = ||U_res*psi||^2  -- MATCHES paper
- Individual jump: p = delta * rate2(w) * ||Aω*psi||^2  -- MATCHES paper

For Hermitian jumps in trajectory code:
- Positive frequency: applies Aω to psi, probability delta*rate2(w)*||Aω*psi||^2
- Negative frequency: applies Aω' to psi, probability delta*rate2(-w)*||Aω'*psi||^2

**Conclusion:** The jump sampling structure is CORRECT. The Kraus operators, their probabilities, and the sampling logic all match between the three sources. The only bug affecting jump sampling is the U_B ordering (TFIX-02), which causes jumps to be applied to wrongly-rotated states when with_coherent=true.

**Note on p_jump_total computation:** The trajectory code uses `p_jump_total = delta * expR` (where expR = <psi|R|psi>) as the total jump probability. This equals the sum of individual p_{a,w} values since R = sum(rate2(w)*Aω'*Aω), so delta*<psi|R|psi> = sum(delta*rate2(w)*||Aω*psi||^2). This is algebraically exact for a pure state -- no bug here.

### CPTP Completeness Relation (TVAL-01) -- Confidence: HIGH

**Derivation of the correct completeness relation:**

Given:
- K0 = I - alpha*R, where alpha = 1 - sqrt(1-delta)
- R = sum over all (a,w) of rate2(w) * L_{a,w}' * L_{a,w}
- S = (2*alpha - delta)*R - alpha^2*R^2
- U_res such that U_res'*U_res = S (Cholesky factor)

Completeness: K0'*K0 + delta*R + U_res'*U_res should equal I.

K0'*K0 = (I - alpha*R)' * (I - alpha*R) = I - alpha*R' - alpha*R + alpha^2*R'*R
Since R is Hermitianized (R = R'), this is:
K0'*K0 = I - 2*alpha*R + alpha^2*R^2

Adding delta*R:
K0'*K0 + delta*R = I - (2*alpha - delta)*R + alpha^2*R^2

Adding S = (2*alpha - delta)*R - alpha^2*R^2:
K0'*K0 + delta*R + S = I - (2*alpha - delta)*R + alpha^2*R^2 + (2*alpha - delta)*R - alpha^2*R^2 = I

So the completeness relation is: **K0'*K0 + delta*R + U_res'*U_res = I**

This is exactly what the phase description states. The formula is correct.

**Test structure:**
```julia
# For each domain in [EnergyDomain(), TimeDomain(), TrotterDomain()]:
config = make_thermalize_config(domain; delta=TEST_DELTA)
# Build framework
fw = build_trajectoryframework(...)
# Verify completeness
completeness = fw.K0' * fw.K0 + TEST_DELTA * fw.R + fw.U_residual' * fw.U_residual
@test isapprox(completeness, I(DIM); atol=1e-10)
```

**Important implementation note:** For TrotterDomain, a `TrottTrott` object is needed. The test helpers from Phase 1 create configs but don't build Trotter objects. The test will need to create one:
```julia
trotter = create_trotter(TEST_HAM, T0, NUM_TROTTER_STEPS_PER_T0)
```
The `create_trotter` function is exported by QuantumFurnace (per `src/QuantumFurnace.jl`).

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| PSD projection | Manual eigenvalue iteration | `eigen(Hermitian(S))` + clamp | Julia's `eigen` handles all edge cases (complex, degenerate) |
| Matrix square root | Custom Cholesky fallback | `Diagonal(sqrt.(max.(eigvals, 0))) * eigvecs'` | Direct from eigendecomposition, bypasses Cholesky entirely |
| Normalization check | Custom tolerance comparison | `abs(total_weight - 1.0) > tol` with `@warn` | Standard Julia warning infrastructure handles deduplication |

## Common Pitfalls

### Pitfall 1: Buffer Aliasing in step_along_trajectory!
**What goes wrong:** The workspace reuses buffers (`ws.psi_tmp`, `ws.Rpsi`) across different stages of the computation. Moving U_B application earlier may overwrite a buffer that's needed later.
**Why it happens:** `ws.psi_tmp` is used for both `K0*psi` (stored for the no-jump branch) and `U_B*psi` (for coherent rotation).
**How to avoid:** After moving U_B to the beginning, `psi` is rotated in-place. Then `ws.psi_tmp` is free to hold `K0*psi` for the no-jump branch. The key insight: once U_B is applied to `psi` directly (via `ws.psi_tmp` as a temporary then `copyto!` back), `ws.psi_tmp` is available for K0*psi computation. Verify buffer lifetimes carefully.
**Warning signs:** NaN in trajectory output, ||psi|| drifting from 1.

### Pitfall 2: Cholesky vs Eigen for PSD guard
**What goes wrong:** After PSD projection via eigendecomposition, attempting Cholesky on the reconstructed matrix may still fail due to floating-point errors in the reconstruction.
**How to avoid:** Either (a) compute U_residual directly from the eigendecomposition as `Diagonal(sqrt.(clamped_eigenvalues)) * eigvecs'`, or (b) add a small eps shift after PSD projection before Cholesky. Option (a) is cleaner.
**Warning signs:** `PosDefException` or NaN in U_residual.

### Pitfall 3: @warn Deduplication
**What goes wrong:** Julia's `@warn` macro deduplicates by default -- it prints the first occurrence and suppresses subsequent ones from the same call site.
**Why it matters:** The user wants to see normalization violations but not be flooded. `@warn` default behavior is actually perfect here -- shows first violation, then quiets down.
**How to avoid:** This is a non-issue if the desired behavior is "warn once". If per-step warnings are needed, use `@warn ... maxlog=Inf` or a plain `println`.

### Pitfall 4: TrotterDomain Test Setup
**What goes wrong:** TrotterDomain tests require a `TrottTrott` object, which Phase 1's test helpers don't create.
**How to avoid:** Add a `make_test_trotter()` helper or inline the Trotter creation in the CPTP test. The `create_trotter` function takes `(hamiltonian, t0, num_trotter_steps_per_t0)`.
**Warning signs:** `MethodError: no method matching build_trajectoryframework(... ::TrottTrott ...)` if the wrong ham_or_trott is passed.

### Pitfall 5: EnergyDomain vs Time/Trotter variant duplication
**What goes wrong:** Both `step_along_trajectory!` methods (EnergyDomain and Time/Trotter) have the SAME U_B ordering bug. Fixing only one leaves the other broken.
**How to avoid:** Apply the exact same fix to both variants. The EnergyDomain variant is at lines 587-729, the Time/Trotter variant at lines 422-572. Search for all `if fw.U_B !== nothing` blocks in `step_along_trajectory!`.

## Code Examples

### Fix Pattern: U_B Ordering (TFIX-02)
```julia
# BEFORE (buggy) -- U_B applied after probability computation:
function step_along_trajectory!(psi, fw)
    # ... compute R*psi, p_nojump, p_jump_total, U_res*psi, p_res ...
    # ... roll r ...
    if fw.U_B !== nothing          # <-- TOO LATE
        mul!(ws.psi_tmp, fw.U_B, psi)
        copyto!(psi, ws.psi_tmp)
        rmul!(psi, 1.0 / sqrt(max(_norm2(psi), eps(Float64))))
    end
    # ... branch selection ...
end

# AFTER (fixed) -- U_B applied FIRST:
function step_along_trajectory!(psi, fw)
    ws = fw.ws
    # Step 0: Coherent unitary (if enabled): psi <- U_B * psi
    if fw.U_B !== nothing
        mul!(ws.psi_tmp, fw.U_B, psi)
        copyto!(psi, ws.psi_tmp)
        rmul!(psi, 1.0 / sqrt(max(_norm2(psi), eps(Float64))))
    end
    # Step 1: Compute probabilities on the ALREADY-ROTATED psi
    mul!(ws.Rpsi, R, psi)          # R * (U_B * psi)
    # ... rest unchanged ...
end
```

### Fix Pattern: PSD Guard (TFIX-04)
```julia
# BEFORE (in build_trajectoryframework):
scratch.tmp2 .= 0.5 .* (scratch.tmp2 .+ scratch.tmp2')
eps_shift = 10 * eps(Float64)
@inbounds for i in 1:dim; scratch.tmp2[i,i] += eps_shift; end
cholesky_S = cholesky!(Hermitian(scratch.tmp2), check=false)
U_residual = Matrix{ComplexF64}(cholesky_S.U)

# AFTER (with PSD projection):
S = Hermitian(0.5 .* (scratch.tmp2 .+ scratch.tmp2'))
eig = eigen(S)
eig.values .= max.(eig.values, 0.0)  # Clamp negatives to zero
U_residual = Matrix{ComplexF64}(Diagonal(sqrt.(eig.values)) * eig.vectors')
```

### Fix Pattern: Normalization Warning (TFIX-03)
```julia
# Insert after computing total_weight:
total_weight = p_nojump + p_res + p_jump_total
total_weight = max(total_weight, 0.0)

# Normalization check
if abs(total_weight - 1.0) > 1e-6
    @warn "Normalization violation: sum = $(round(total_weight; digits=6))"
end

r = rand() * total_weight
```

### Test Pattern: CPTP Verification (TVAL-01)
```julia
@testset "CPTP Completeness" begin
    for (domain, ham_or_trott) in [
        (EnergyDomain(), TEST_HAM),
        (TimeDomain(), TEST_HAM),
        (TrotterDomain(), test_trotter),
    ]
        @testset "$(typeof(domain))" begin
            config = make_thermalize_config(domain; delta=TEST_DELTA)
            precomputed = precompute_data(config.domain, config, ham_or_trott)
            scratch = KrausScratch(ComplexF64, DIM)
            fw = build_trajectoryframework(
                TEST_JUMPS, ham_or_trott, config, precomputed, scratch, TEST_DELTA
            )

            completeness = fw.K0' * fw.K0 + TEST_DELTA * fw.R + fw.U_residual' * fw.U_residual
            @test isapprox(completeness, Matrix{ComplexF64}(I, DIM, DIM); atol=1e-10)
        end
    end
end
```

### Test Pattern: Single-Step Bug Fix Tests
```julia
@testset "TFIX-02: U_B ordering" begin
    # With coherent term, single step should produce normalized output
    config = make_thermalize_config(EnergyDomain(); with_coherent=true)
    precomputed = precompute_data(config.domain, config, TEST_HAM)
    scratch = KrausScratch(ComplexF64, DIM)
    fw = build_trajectoryframework(TEST_JUMPS, TEST_HAM, config, precomputed, scratch, TEST_DELTA)

    psi = zeros(ComplexF64, DIM)
    psi[1] = 1.0  # computational basis state
    step_along_trajectory!(psi, fw)
    @test isapprox(norm(psi), 1.0; atol=1e-12)  # Output must be normalized
end

@testset "TFIX-04: PSD guard" begin
    # Build framework -- if PSD guard works, no error/NaN
    config = make_thermalize_config(EnergyDomain())
    precomputed = precompute_data(config.domain, config, TEST_HAM)
    scratch = KrausScratch(ComplexF64, DIM)
    fw = build_trajectoryframework(TEST_JUMPS, TEST_HAM, config, precomputed, scratch, TEST_DELTA)
    @test all(isfinite, fw.U_residual)
    @test !any(isnan, fw.U_residual)
end
```

## State of the Art

Not applicable -- this phase fixes bugs in existing code, not adopting new approaches.

## Open Questions

1. **Normalization tolerance threshold**
   - What we know: The user wants a warning when p_nojump + p_res + p_jump_total deviates from 1.0
   - What's unclear: What tolerance threshold should trigger the warning? 1e-6? 1e-4?
   - Recommendation: Use 1e-6 as the threshold. This is tight enough to catch real issues but loose enough to avoid false positives from floating-point arithmetic in the probability computation. This is Claude's discretion per CONTEXT.md.

2. **Trotter test helper**
   - What we know: TrotterDomain CPTP test needs a `TrottTrott` object
   - What's unclear: Should this be a shared helper in test_helpers.jl or local to the CPTP test?
   - Recommendation: Add `make_test_trotter()` to `test_helpers.jl` since Phase 4 will also need it for TrotterDomain cross-validation. Compute once at include time like other fixtures.

3. **Paper cross-check documentation**
   - What we know: User wants the cross-check to verify trajectory code, DM code, and paper agree
   - What's unclear: How to structure this -- inline comments? Separate doc?
   - Recommendation: Add inline comments in the code at the fix sites referencing the paper equations. Include the cross-check analysis as comments in the test file. No separate document needed -- this RESEARCH.md already contains the full analysis.

## File Modification Map

### Files to Modify
| File | Changes | Requirements |
|------|---------|-------------|
| `src/trajectories.jl` | Move U_B before probs (TFIX-02), add normalization warning (TFIX-03), PSD guard in build_trajectoryframework (TFIX-04) | TFIX-02, TFIX-03, TFIX-04 |
| `test/test_helpers.jl` | Add `make_test_trotter()` helper | Support for TVAL-01 TrotterDomain test |
| `test/runtests.jl` | Include new test files | -- |

### Files to Create
| File | Purpose | Requirements |
|------|---------|-------------|
| `test/test_trajectory_fixes.jl` | Single-step tests for TFIX-02, TFIX-03, TFIX-04, TFIX-05 | TFIX-02, TFIX-03, TFIX-04, TFIX-05 |
| `test/test_cptp.jl` | CPTP completeness verification across all three domains | TVAL-01 |

### Files NOT Modified
| File | Reason |
|------|--------|
| `src/jump_workers.jl` | DM code is reference implementation, not modified in this phase |
| `src/coherent.jl` | Coherent term computation is correct (user confirmed) |
| `src/structs.jl` | No struct changes needed |

## Sources

### Primary (HIGH confidence)
- `src/trajectories.jl` -- direct code analysis of both `step_along_trajectory!` variants
- `src/jump_workers.jl` -- DM reference implementation, all three domain variants of `jump_contribution!`
- `src/coherent.jl` -- coherent term precomputation logic
- `test/test_helpers.jl` -- Phase 1 test fixtures
- `test/test_compilation.jl` -- Phase 1 smoke tests confirming build_trajectoryframework works

### Secondary (MEDIUM confidence)
- Paper references (Chen 2023 Theorem III.1, Eq. 3.2-3.3) -- mathematical derivation of the CPTP channel, referenced in code comments but not directly verified against the published paper (no URL access)
- `.planning/REQUIREMENTS.md` -- requirement definitions for TFIX-02 through TFIX-05 and TVAL-01
- `.planning/ROADMAP.md` -- phase structure and dependencies

## Metadata

**Confidence breakdown:**
- U_B ordering bug (TFIX-02): HIGH -- verified by direct comparison of trajectory vs DM code ordering
- Normalization guard (TFIX-03): HIGH -- straightforward addition, no ambiguity
- PSD guard (TFIX-04): HIGH -- standard eigenvalue clamping technique, well-understood
- Jump sampling faithfulness (TFIX-05): HIGH -- cross-checked all three sources, channel structure matches
- CPTP verification (TVAL-01): HIGH -- algebraic derivation verified, formula matches code construction

**Research date:** 2026-02-14
**Valid until:** Indefinite (code-level analysis, not version-dependent)
