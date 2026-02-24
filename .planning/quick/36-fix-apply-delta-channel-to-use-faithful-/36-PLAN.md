---
phase: quick-36
plan: 01
type: execute
wave: 1
depends_on: []
files_modified:
  - src/krylov_workspace.jl
  - src/krylov_eigsolve.jl
  - test/test_krylov_eigsolve.jl
autonomous: true

must_haves:
  truths:
    - "apply_delta_channel! implements Chen Eq. 3.2 with K0, jump sandwich, and U_residual terms"
    - "Channel path spectral gap matches dense Lindbladian reference to O(delta^2) accuracy"
    - "Precomputed channel matrices (K0, U_residual, U_coherent) are computed once at workspace construction, not per matvec"
    - "Existing channel eigsolve tests pass with tighter or equal tolerance"
  artifacts:
    - path: "src/krylov_workspace.jl"
      provides: "KrylovWorkspace with channel-specific precomputed fields and ThermalizeConfig constructor"
    - path: "src/krylov_eigsolve.jl"
      provides: "Faithful Chen channel apply_delta_channel! and updated channel path in krylov_spectral_gap"
    - path: "test/test_krylov_eigsolve.jl"
      provides: "Updated tests verifying faithful channel vs dense Lindbladian and CPTP property"
  key_links:
    - from: "src/krylov_workspace.jl"
      to: "src/krylov_eigsolve.jl"
      via: "KrylovWorkspace channel fields (channel_K0, channel_U_residual, etc.) used in apply_delta_channel!"
    - from: "src/krylov_eigsolve.jl"
      to: "src/jump_workers.jl"
      via: "Mirrors _finalize_kraus_step! logic for K0, S, U_residual computation"
---

<objective>
Replace the Euler-approximate `apply_delta_channel!` (rho + delta*L(rho)) with Chen's faithful
CPTP quantum channel (Eq. 3.2), matching the thermalization code in `_finalize_kraus_step!`.
Precompute rho-independent channel matrices (R_total, K0, U_residual, U_coherent) at workspace
construction time so the per-matvec cost is only the rho-dependent sandwich terms.

Purpose: The Euler approximation is O(delta) accurate and not CPTP. The Chen channel is O(delta^2)
accurate and exactly trace-preserving, matching what `run_thermalization` actually implements.
This eliminates a systematic error source in the Krylov channel eigsolve path.

Output: Updated krylov_workspace.jl, krylov_eigsolve.jl, test_krylov_eigsolve.jl
</objective>

<execution_context>
@./.claude/get-shit-done/workflows/execute-plan.md
@./.claude/get-shit-done/templates/summary.md
</execution_context>

<context>
@src/krylov_workspace.jl
@src/krylov_eigsolve.jl
@src/krylov_matvec.jl
@src/jump_workers.jl (lines 160-212: _finalize_kraus_step!, lines 312-389: EnergyDomain _jump_contribution!)
@src/coherent.jl (lines 1-100: _precompute_coherent_total_B, _precompute_coherent_unitary_terms)
@src/kraus.jl (lines 1-15: KrausScratch struct definition)
@test/test_krylov_eigsolve.jl
@test/test_helpers.jl
</context>

<tasks>

<task type="auto">
  <name>Task 1: Add channel fields to KrylovWorkspace and ThermalizeConfig constructor</name>
  <files>src/krylov_workspace.jl</files>
  <action>
Add new optional fields to the `KrylovWorkspace` struct for the precomputed channel matrices.
These are `Union{Nothing, Matrix{T}}` fields, defaulting to `nothing` for the Lindbladian path:

```
channel_K0::Union{Nothing, Matrix{T}}           # I - alpha * R_total
channel_U_residual::Union{Nothing, Matrix{T}}    # sqrt_psd(S) from Chen Eq. 3.2
channel_U_coherent::Union{Nothing, Matrix{T}}    # exp(-i*delta*B_total) for exponentiated coherent
channel_rho_jump::Union{Nothing, Matrix{T}}      # scratch for jump sandwich accumulation
```

Add these fields AFTER `rho_out` in the struct definition. Update the existing `KrylovWorkspace`
constructor for `AbstractLiouvConfig` to pass `nothing` for all four new fields.

Then add a NEW constructor method dispatching on `AbstractThermalizeConfig`:

```julia
function KrylovWorkspace(
    config::AbstractThermalizeConfig,
    hamiltonian::HamHam,
    jumps::Vector{JumpOp};
    trotter::Union{TrottTrott, Nothing}=nothing,
)
```

This constructor should:

1. Call the SAME precomputation steps as the LiouvConfig constructor: `_precompute_data`,
   `_precompute_coherent_total_B`, extract `jump_eigenbases`, `jump_hermitian`, allocate
   `jump_oft`, `tmp1`, `tmp2`, `LdagL`, `rho_out`. Use `_thermalize_to_liouv_config(config)`
   to get a LiouvConfig for `_precompute_data` and `_precompute_coherent_total_B` calls
   (these dispatch on `AbstractLiouvConfig`/`AbstractConfig`).

   NOTE: `_thermalize_to_liouv_config` is defined in `krylov_eigsolve.jl`. To avoid circular
   dependency, INLINE the config conversion logic directly in this constructor, OR move
   `_thermalize_to_liouv_config` to krylov_workspace.jl. Prefer moving it -- it's a pure
   data transformation with no eigsolve dependency.

2. Compute R_total (physics convention, NOT kron convention) by looping over all jumps and
   energy labels, accumulating `rate^2 * (L_omega' * L_omega)` for each frequency. This mirrors
   the R accumulation in `_jump_contribution!` for EnergyDomain (jump_workers.jl:339-381)
   but without any rho-dependent terms. The domain dispatch logic:

   **EnergyDomain:**
   ```
   prefactor = (w0 / (sigma * sqrt(2pi))) * gamma_norm_factor
   for each jump eigenbasis:
     for w in energy_labels (half-grid for hermitian):
       L_omega = oft(eigenbasis, w, bohr_freqs, sigma)  # Gaussian filter
       rate2 = prefactor * transition(w)
       R += rate2 * (L_omega' * L_omega)
       if hermitian && w > 0:
         rate2_neg = prefactor * transition(-w)
         R += rate2_neg * (L_omega * L_omega')  # negative freq: L†L -> LL†
   ```

   **TimeDomain/TrotterDomain:**
   Same structure but OFT uses NUFFT prefactors instead of Gaussian:
   ```
   prefactor = w0 * t0^2 * (sigma * sqrt(2/pi)) / (2pi) * gamma_norm_factor
   L_omega = eigenbasis .* nufft_prefactor_matrix
   ```

   **BohrDomain:**
   ```
   for each jump eigenbasis, for each nu_2 in bohr_dict keys:
     alpha_A = alpha(bohr_freqs, nu_2) .* eigenbasis
     build A_nu2_dag by scatter: A_nu2_dag[j,i] = conj(eigenbasis[i,j]) for (i,j) in bucket
     R += gamma_norm_factor * (A_nu2_dag * alpha_A)
   ```

   Use `ws.tmp1`, `ws.LdagL` as scratch during this one-time computation. After the loop,
   call `hermitianize!(R)`.

3. Compute K0, S, U_residual following `_finalize_kraus_step!` logic (jump_workers.jl:176-198):
   ```julia
   delta = config.delta
   alpha_chen = 1 - sqrt(1 - delta)  # NOT the Bohr alpha function!
   K0 = Matrix{CT}(I, dim, dim) - alpha_chen * R_total
   # S = (2*alpha_chen - delta)*R - alpha_chen^2 * R^2
   R2 = R_total * R_total
   S = (2*alpha_chen - delta) .* R_total .- (alpha_chen^2) .* R2
   hermitianize!(S)
   eig = eigen(Hermitian(S))
   eig.values .= max.(eig.values, 0.0)
   U_residual = Matrix{CT}(Diagonal(sqrt.(eig.values)) * eig.vectors')
   ```

4. Compute U_coherent = `exp(-1im * delta * Hermitian(B_total))` if B_total is not nothing,
   else nothing. This is a single matrix exponential at construction time. Note: Use
   `Hermitian(B_total)` for numerically stable exponentiation (B is Hermitian by construction).

5. Allocate `channel_rho_jump = zeros(CT, dim, dim)` scratch.

6. Return `KrylovWorkspace{CT, typeof(precomputed_data)}(...)` with all fields including the
   new channel fields.

IMPORTANT: The variable name `alpha` is already used in BohrDomain precomputed_data as a
function `alpha(bohr_freqs, nu_2)`. Use `alpha_chen` for the Chen channel parameter
`1 - sqrt(1 - delta)` to avoid shadowing.
  </action>
  <verify>
Julia REPL or test: construct a KrylovWorkspace from a ThermalizeConfig and verify the channel
fields are populated (not nothing). Verify the LiouvConfig constructor still produces nothing
for channel fields. Run `julia --project -e "using QuantumFurnace"` to verify no syntax errors.
  </verify>
  <done>
KrylovWorkspace struct has 4 new Union{Nothing,Matrix} fields. LiouvConfig constructor passes
nothing for them. ThermalizeConfig constructor precomputes R_total, K0, U_residual, U_coherent
and populates them. `_thermalize_to_liouv_config` moved to krylov_workspace.jl.
  </done>
</task>

<task type="auto">
  <name>Task 2: Rewrite apply_delta_channel! and update channel eigsolve path</name>
  <files>src/krylov_eigsolve.jl, test/test_krylov_eigsolve.jl</files>
  <action>
**Part A: Rewrite `apply_delta_channel!` (krylov_eigsolve.jl lines 214-226)**

Replace the Euler formula with the faithful Chen channel. The new function uses precomputed
channel matrices from the workspace. The rho-dependent computation per matvec call is:

```julia
function apply_delta_channel!(ws::KrylovWorkspace{T}, rho::Matrix{T}) where {T<:Complex}
    delta = ...  # no longer needed as parameter -- everything is precomputed

    # 1. Apply coherent unitary (if present)
    #    rho1 = U_coherent * rho * U_coherent'
    #    Store in ws.rho_out as intermediate
    U_coh = ws.channel_U_coherent
    if U_coh !== nothing
        mul!(ws.tmp1, U_coh, rho)          # tmp1 = U_coh * rho
        mul!(ws.rho_out, ws.tmp1, U_coh')   # rho_out = U_coh * rho * U_coh'
        rho_eff = ws.rho_out  # alias -- rho after coherent rotation
    else
        rho_eff = rho
    end

    # 2. Accumulate rho_jump = delta * sum_l L_l * rho_eff * L_l'
    #    This reuses the EXISTING dissipator loop from apply_lindbladian! but
    #    only the SANDWICH terms (no anticommutator).
    #    However, instead of re-looping, we use a simpler approach:
    #    Call apply_lindbladian! to get L(rho_eff), then extract the sandwich part.
    #
    #    Actually, it's cleaner to compute rho_jump directly. The sandwich terms
    #    in physics convention are: delta * sum rate^2 * L * rho * L'
    #    This matches the accumulation in _jump_contribution! for rho_jump.
    #
    #    For simplicity and correctness, compute rho_jump using:
    #    rho_jump = apply_lindbladian!(ws, rho_eff, config_liouv, ham) * delta + delta*{R, rho_eff}/2 + rho_eff - rho_eff
    #    NO -- that's circular.
    #
    #    CORRECT APPROACH: The sandwich term rho_jump = delta * sum L_l rho L_l'
    #    can be extracted from the Lindbladian as:
    #      L(rho) = sum D_l(rho) + coherent = sum (L rho L' - 0.5{L'L, rho}) + coherent
    #    So: sum L rho L' = L(rho) - coherent + 0.5{R, rho}
    #    where R = sum L'L.
    #    And: rho_jump = delta * (L(rho) - coherent + 0.5{R, rho})
    #
    #    BUT we already have R precomputed and L(rho) is computable via apply_lindbladian!.
    #    HOWEVER, apply_lindbladian! uses the KRON convention (conj(L)*rho*L^T, (L'L)^T),
    #    not the physics convention (L*rho*L', L'L). The R we precomputed is in physics
    #    convention. So we can't simply add back 0.5{R, rho} to the kron-convention L(rho).
    #
    #    SIMPLEST CORRECT APPROACH: Compute the sandwich terms DIRECTLY in physics convention
    #    by looping over jumps and frequencies, similar to the _jump_contribution! code
    #    but only the sandwich part. This avoids convention mixing.

    # Compute rho_jump = delta * sum_l rate^2 * L_omega * rho_eff * L_omega'
    # (physics convention sandwich, matches jump_workers.jl accumulation)
    fill!(ws.channel_rho_jump, 0)
    _accumulate_channel_jump_sandwich!(ws, rho_eff)

    # 3. Assemble: rho_out = K0 * rho_eff * K0' + rho_jump + U_res * rho_eff * U_res'
    K0 = ws.channel_K0
    U_res = ws.channel_U_residual

    # If we applied coherent and rho_eff aliases rho_out, we need to copy rho_eff first
    # to avoid overwriting it when we write to rho_out.
    # Solution: if U_coh !== nothing, copy rho_out (which IS rho_eff) to channel_rho_jump's
    # place? No -- channel_rho_jump is in use.
    # Better: use tmp2 to hold rho_eff if it aliases rho_out.
    if U_coh !== nothing
        copyto!(ws.tmp2, ws.rho_out)  # tmp2 = rho_eff (safe copy)
        rho_eff_safe = ws.tmp2
    else
        rho_eff_safe = rho  # original rho, not aliased
    end

    mul!(ws.tmp1, K0, rho_eff_safe)         # tmp1 = K0 * rho_eff
    mul!(ws.rho_out, ws.tmp1, K0')          # rho_out = K0 * rho_eff * K0'

    ws.rho_out .+= ws.channel_rho_jump       # += rho_jump

    mul!(ws.tmp1, U_res, rho_eff_safe)       # tmp1 = U_res * rho_eff
    mul!(ws.rho_out, ws.tmp1, U_res', 1.0, 1.0)  # += U_res * rho_eff * U_res'

    return ws.rho_out
end
```

WAIT -- the above pseudocode has a problem. We need a separate function to accumulate the
sandwich terms. Create a new helper function `_accumulate_channel_jump_sandwich!` that fills
`ws.channel_rho_jump` with `delta * sum rate^2 * L * rho * L'`. This needs to be domain-aware.

ACTUALLY, re-think the approach. The channel-path `krylov_spectral_gap` for
`AbstractThermalizeConfig` currently converts to LiouvConfig and uses `apply_lindbladian!`.
Instead, the new approach should:

1. The workspace stores the domain config (via precomputed_data) so we know which domain to use.
2. `apply_delta_channel!` needs to know the domain to loop over frequencies correctly.
3. Pass `config` and `hamiltonian` to `apply_delta_channel!` so it can do the sandwich loop.

NEW SIGNATURE:
```julia
function apply_delta_channel!(
    ws::KrylovWorkspace{T},
    rho::Matrix{T},
    config::AbstractLiouvConfig,  # keep LiouvConfig for dissipator dispatch
    hamiltonian::HamHam,
) where {T<:Complex}
```

The function needs access to the domain-specific dissipator loop for sandwich-only terms.
Create domain-dispatched `_accumulate_jump_sandwich!` helpers (one per domain) that mirror
`apply_lindbladian!` but ONLY accumulate `delta * L * rho * L'` (no anticommutator, no coherent).

For EnergyDomain:
```julia
function _accumulate_jump_sandwich!(
    out::Matrix{T},       # accumulator for rho_jump
    ws::KrylovWorkspace{T},
    rho::Matrix{T},       # rho_eff after coherent rotation
    delta::Real,
    config::AbstractLiouvConfig{EnergyDomain},
    hamiltonian::HamHam,
) where {T<:Complex}
    (; transition, gamma_norm_factor, energy_labels) = ws.precomputed_data
    bohr_freqs = hamiltonian.bohr_freqs
    inv_4sigma2 = 1.0 / (4 * config.sigma^2)
    prefactor = (config.w0 / (config.sigma * sqrt(2*pi))) * gamma_norm_factor

    for (k, eigenbasis) in enumerate(ws.jump_eigenbases)
        is_herm = ws.jump_hermitian[k]
        if is_herm
            for w_raw in energy_labels
                w_raw > 1e-12 && continue
                w = abs(w_raw)
                _krylov_oft!(ws.jump_oft, eigenbasis, bohr_freqs, w, inv_4sigma2)
                rate2 = prefactor * transition(w)
                # Sandwich: delta * rate2 * L * rho * L'  (physics convention)
                # L = ws.jump_oft
                mul!(ws.tmp1, rho, ws.jump_oft')          # tmp1 = rho * L'
                mul!(out, ws.jump_oft, ws.tmp1, delta*rate2, 1.0)  # out += delta*rate2 * L * rho * L'
                if w > 1e-12
                    rate2_neg = prefactor * transition(-w)
                    # Neg freq: L_neg = L' (adjoint), so sandwich = L' * rho * L
                    mul!(ws.tmp1, rho, ws.jump_oft)        # tmp1 = rho * L
                    mul!(out, ws.jump_oft', ws.tmp1, delta*rate2_neg, 1.0)  # out += delta*rate2_neg * L' * rho * L
                end
            end
        else
            for w in energy_labels
                _krylov_oft!(ws.jump_oft, eigenbasis, bohr_freqs, w, inv_4sigma2)
                rate2 = prefactor * transition(w)
                mul!(ws.tmp1, rho, ws.jump_oft')
                mul!(out, ws.jump_oft, ws.tmp1, delta*rate2, 1.0)
            end
        end
    end
    return nothing
end
```

Create analogous methods for TimeDomain/TrotterDomain (same but NUFFT prefactor OFT) and
BohrDomain (Bohr bucket iteration with A_nu2 and alpha_A sandwich).

For BohrDomain, the sandwich in physics convention is:
`delta * gamma_norm_factor * alpha_A * rho * A_nu2_dag'`
Wait -- need to think about what the physics-convention sandwich is for BohrDomain.

In the thermalization code (jump_workers.jl:277):
```
rho_jump += scaled_delta * B_nu2 * (rho * A_nu2_dag)  -- where B_nu2 = alpha_A
```
This is: `delta * gamma_norm_factor * alpha_A * rho * A_nu2_dag`

where `A_nu2_dag[j,i] = conj(eigenbasis[i,j])` for (i,j) in bucket, so
`A_nu2_dag` is NOT the transpose of `A_nu2`. It's a scattered sparse construction.
And the product `rho * A_nu2_dag` is done entrywise in the thermalization code.

So for BohrDomain:
```julia
function _accumulate_jump_sandwich!(
    out::Matrix{T},
    ws::KrylovWorkspace{T},
    rho::Matrix{T},
    delta::Real,
    config::AbstractLiouvConfig{BohrDomain},
    hamiltonian::HamHam,
) where {T<:Complex}
    (; alpha, gamma_norm_factor) = ws.precomputed_data
    dim = size(rho, 1)
    A_nu2_dag = zeros(T, dim, dim)  # one alloc per call, acceptable

    for (k, eigenbasis) in enumerate(ws.jump_eigenbases)
        for nu_2 in keys(hamiltonian.bohr_dict)
            @. ws.jump_oft = alpha(hamiltonian.bohr_freqs, nu_2) * eigenbasis  # alpha_A = B_nu2

            fill!(A_nu2_dag, 0)
            indices = hamiltonian.bohr_dict[nu_2]
            @inbounds for idx in indices
                i = idx[1]; j = idx[2]
                A_nu2_dag[j, i] = conj(eigenbasis[i, j])
            end

            # rho_jump += delta * gamma_norm_factor * alpha_A * rho * A_nu2_dag
            mul!(ws.tmp1, rho, A_nu2_dag)  # tmp1 = rho * A_nu2_dag
            mul!(out, ws.jump_oft, ws.tmp1, delta * gamma_norm_factor, 1.0)
        end
    end
    return nothing
end
```

For TimeDomain/TrotterDomain:
Same as EnergyDomain but with NUFFT prefactors:
```julia
function _accumulate_jump_sandwich!(
    out::Matrix{T},
    ws::KrylovWorkspace{T},
    rho::Matrix{T},
    delta::Real,
    config::AbstractLiouvConfig{D},
    hamiltonian::HamHam,
) where {T<:Complex, D<:Union{TimeDomain, TrotterDomain}}
    (; transition, gamma_norm_factor, energy_labels, oft_nufft_prefactors) = ws.precomputed_data
    prefactor = config.w0 * config.t0^2 * (config.sigma * sqrt(2/pi)) / (2*pi) * gamma_norm_factor

    for (k, eigenbasis) in enumerate(ws.jump_eigenbases)
        is_herm = ws.jump_hermitian[k]
        if is_herm
            for w_raw in energy_labels
                w_raw > 1e-12 && continue
                w = abs(w_raw)
                nufft_pf = _prefactor_view(oft_nufft_prefactors, w)
                @. ws.jump_oft = eigenbasis * nufft_pf
                rate2 = prefactor * transition(w)
                mul!(ws.tmp1, rho, ws.jump_oft')
                mul!(out, ws.jump_oft, ws.tmp1, delta*rate2, 1.0)
                if w > 1e-12
                    rate2_neg = prefactor * transition(-w)
                    mul!(ws.tmp1, rho, ws.jump_oft)
                    mul!(out, ws.jump_oft', ws.tmp1, delta*rate2_neg, 1.0)
                end
            end
        else
            for w in energy_labels
                nufft_pf = _prefactor_view(oft_nufft_prefactors, w)
                @. ws.jump_oft = eigenbasis * nufft_pf
                rate2 = prefactor * transition(w)
                mul!(ws.tmp1, rho, ws.jump_oft')
                mul!(out, ws.jump_oft, ws.tmp1, delta*rate2, 1.0)
            end
        end
    end
    return nothing
end
```

Now `apply_delta_channel!` becomes:

```julia
function apply_delta_channel!(
    ws::KrylovWorkspace{T},
    rho::Matrix{T},
    config::AbstractLiouvConfig,
    hamiltonian::HamHam,
) where {T<:Complex}
    K0 = ws.channel_K0
    U_res = ws.channel_U_residual
    U_coh = ws.channel_U_coherent
    delta = ...  # need delta for sandwich scaling

    # 1. Coherent rotation
    if U_coh !== nothing
        mul!(ws.tmp1, U_coh, rho)
        mul!(ws.LdagL, ws.tmp1, U_coh')  # LdagL as scratch for rho_eff
        rho_eff = ws.LdagL
    else
        rho_eff = rho
    end

    # 2. Jump sandwich
    fill!(ws.channel_rho_jump, 0)
    _accumulate_jump_sandwich!(ws.channel_rho_jump, ws, rho_eff, delta, config, hamiltonian)

    # 3. Assemble K0 * rho_eff * K0' + rho_jump + U_res * rho_eff * U_res'
    mul!(ws.tmp1, K0, rho_eff)
    mul!(ws.rho_out, ws.tmp1, K0')
    ws.rho_out .+= ws.channel_rho_jump
    mul!(ws.tmp1, U_res, rho_eff)
    mul!(ws.rho_out, ws.tmp1, U_res', 1.0, 1.0)

    return ws.rho_out
end
```

The delta is needed for the sandwich scaling. Store it in the workspace or pass it. Since it's
already in the ThermalizeConfig and the workspace is config-specific, add a `channel_delta`
field to KrylovWorkspace (Union{Nothing, Float64}).

**Part B: Update `krylov_spectral_gap` for `AbstractThermalizeConfig`**

In the channel path (krylov_eigsolve.jl:356-438), change:
1. Construct workspace using ThermalizeConfig constructor: `ws = KrylovWorkspace(config, hamiltonian, jumps; trotter=trotter)` -- this now dispatches to the new constructor.
2. Still need a LiouvConfig for the sandwich loop dispatch: `config_liouv = _thermalize_to_liouv_config(config)` (now in krylov_workspace.jl).
3. Update channel_matvec closure to call the new `apply_delta_channel!(ws, rho, config_liouv, hamiltonian)`.
4. Everything else (eigsolve, sorting, conversion) stays the same.

**Part C: Update `test_krylov_eigsolve.jl`**

1. **Update Testset 1** ("apply_delta_channel! round-trip vs dense"): The old test compared
   `apply_delta_channel!` output to `(I + delta*L)*vec(rho)`. With the faithful channel, the
   output is no longer `rho + delta*L(rho)`. Instead, test that:
   - The channel is trace-preserving: `tr(E(rho)) approx tr(rho)` to high precision
   - The channel is positive: eigenvalues of E(rho) are non-negative (for valid density matrix input)
   - The channel output differs from the Euler approximation by O(delta^2):
     `norm(E_chen(rho) - E_euler(rho)) < C * delta^2` for some reasonable C

   Replace the test body:
   ```julia
   @testset "apply_delta_channel! faithful Chen channel" begin
       config_therm = make_thermalize_config(EnergyDomain(); with_coherent=true, delta=0.01)
       config_liouv = make_liouv_config(EnergyDomain(); with_coherent=true)
       ws = KrylovWorkspace(config_therm, TEST_HAM, TEST_JUMPS)
       delta = config_therm.delta

       # Dense Lindbladian for Euler comparison
       L_dense = construct_lindbladian(TEST_JUMPS, config_liouv, TEST_HAM)
       I_d2 = Matrix{ComplexF64}(LinearAlgebra.I(DIM^2))

       for _ in 1:5
           rho = Matrix(random_density_matrix(NUM_QUBITS))
           # Faithful Chen channel
           apply_delta_channel!(ws, rho, config_liouv, TEST_HAM)
           rho_chen = copy(ws.rho_out)

           # Trace preservation
           @test isapprox(real(tr(rho_chen)), real(tr(rho)); atol=1e-10)

           # Positivity (eigenvalues >= -eps)
           eigs = eigvals(Hermitian(rho_chen))
           @test all(eigs .> -1e-10)

           # O(delta^2) close to Euler
           v_euler = (I_d2 + delta * L_dense) * vec(rho)
           @test norm(vec(rho_chen) - v_euler) < 50 * delta^2  # C=50 generous bound
       end
   end
   ```

2. **Testset 5** ("Channel eigsolve accuracy"): Should pass with potentially tighter tolerance
   since the channel is now faithful. Keep rtol=1e-3 (the conversion formula lambda=(mu-1)/delta
   still introduces O(delta) error in the eigenvalue mapping, even though the channel itself is
   more accurate). The test should still pass.

3. **Add new testset** verifying that the faithful channel's steady state matches Gibbs better
   than the Euler channel did (optional, only if easy to verify).
  </action>
  <verify>
Run `julia --project -e "using Pkg; Pkg.test()"` and verify all tests pass. Specifically:
- The new "apply_delta_channel! faithful Chen channel" testset passes
- Testset 5 "Channel eigsolve accuracy" still passes with rtol=1e-3
- All other existing testsets unaffected
  </verify>
  <done>
apply_delta_channel! implements Chen Eq. 3.2 with precomputed K0, U_residual, U_coherent.
The channel path in krylov_spectral_gap uses the ThermalizeConfig workspace constructor.
Tests verify trace preservation, positivity, O(delta^2) agreement with Euler, and channel
eigsolve accuracy against dense reference. All tests pass.
  </done>
</task>

</tasks>

<verification>
1. `julia --project -e "using Pkg; Pkg.test()"` -- all tests pass
2. Channel eigsolve spectral gap matches dense reference to rtol=1e-3
3. apply_delta_channel! output is trace-preserving (|tr(E(rho)) - tr(rho)| < 1e-10)
4. apply_delta_channel! output has non-negative eigenvalues for valid density matrix input
5. Faithful channel differs from Euler by O(delta^2): norm(E_chen - E_euler) < C*delta^2
</verification>

<success_criteria>
- KrylovWorkspace has ThermalizeConfig constructor that precomputes K0, U_residual, U_coherent
- apply_delta_channel! uses Chen Eq. 3.2 (K0*rho*K0' + rho_jump + U_res*rho*U_res')
- Channel path krylov_spectral_gap uses new workspace constructor
- All existing tests pass; new tests verify CPTP property and O(delta^2) accuracy
</success_criteria>

<output>
After completion, create `.planning/quick/36-fix-apply-delta-channel-to-use-faithful-/36-SUMMARY.md`
</output>
