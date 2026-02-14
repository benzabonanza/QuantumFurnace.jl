---
phase: quick-7
plan: 01
type: execute
wave: 1
depends_on: []
files_modified:
  - src/trajectories.jl
  - test/test_cptp.jl
autonomous: true
must_haves:
  truths:
    - "Each trajectory step randomly selects ONE operator a=rand(1:N_jumps) and applies only that operator's channel"
    - "Per-operator Kraus channel is CPTP: K0_a'*K0_a + delta_eff*R_a + U_res_a'*U_res_a = I for each a"
    - "Rates are rescaled by 1/p_jump (=N_jumps) so physical mixing time semantics match DM run_thermalization"
    - "All three domains (Energy, Time, Trotter) use per-operator stepping"
    - "Existing test suite passes (compilation, CPTP, DM tests)"
  artifacts:
    - path: "src/trajectories.jl"
      provides: "PerOperatorKraus struct, refactored TrajectoryFramework with per_operator vector, per-operator step_along_trajectory!"
      contains: "struct PerOperatorKraus"
    - path: "test/test_cptp.jl"
      provides: "Updated CPTP test verifying per-operator channel completeness"
      contains: "per_operator"
  key_links:
    - from: "src/trajectories.jl (build_trajectoryframework)"
      to: "src/trajectories.jl (precompute_R)"
      via: "Loop over individual jumps, calling precompute_R for single-jump vector [jump]"
      pattern: "precompute_R.*\\[jump\\]"
    - from: "src/trajectories.jl (step_along_trajectory!)"
      to: "src/trajectories.jl (PerOperatorKraus)"
      via: "Random operator selection then per_op field access"
      pattern: "rand.*1:.*per_operator"
    - from: "src/trajectories.jl (build_trajectoryframework)"
      to: "src/coherent.jl (precompute_coherent_unitary_terms)"
      via: "Per-operator B_a and U_B_a computation with 1/p_jump rescaling"
      pattern: "precompute_coherent_unitary_terms"
---

<objective>
Refactor the trajectory simulator to use per-operator Lie-Trotter splitting instead of a single combined channel.

Purpose: The current `step_along_trajectory!` precomputes a SINGLE total R = sum_a R^a, single K0, single U_residual, and single U_B from ALL operators simultaneously. This implements the full Lindblad channel in one step, NOT the per-operator Lie-Trotter splitting that Chen's quantum algorithm actually uses. The DM simulator (`run_thermalization` in `src/furnace.jl`) already does the correct thing: randomly picks operator `a`, applies only that operator's channel with rates rescaled by `1/p_jump`. This refactor makes the trajectory simulator match that structure.

Output: Refactored `src/trajectories.jl` with per-operator Kraus data and stepping, updated `test/test_cptp.jl`.
</objective>

<execution_context>
@./.claude/get-shit-done/workflows/execute-plan.md
@./.claude/get-shit-done/templates/summary.md
</execution_context>

<context>
@src/trajectories.jl
@src/furnace.jl (run_thermalization -- reference for per-operator random selection with rescale_by_inv_prob)
@src/jump_workers.jl (jump_contribution! -- reference for per-operator channel structure)
@src/coherent.jl (precompute_coherent_unitary_terms -- already computes per-operator B_a/U_B_a)
@test/test_cptp.jl
@test/test_helpers.jl
</context>

<tasks>

<task type="auto">
  <name>Task 1: Add PerOperatorKraus struct and refactor TrajectoryFramework + build_trajectoryframework</name>
  <files>src/trajectories.jl</files>
  <action>
1. Define a new struct `PerOperatorKraus{T}` near the top of the file (after `TrajectoryWorkspace`, before `TrajectoryFramework`):
```julia
struct PerOperatorKraus{T}
    R::Matrix{T}            # R^a for this operator
    K0::Matrix{T}           # I - alpha * R^a
    U_residual::Matrix{T}   # sqrt(S^a) via eigendecomposition
    U_B::Union{Nothing, Matrix{T}}  # exp(-i * delta_eff * B^a), or nothing
end
```

2. Modify the `TrajectoryFramework` struct:
   - REMOVE fields: `B`, `U_B`, `R`, `K0`, `U_residual`
   - ADD field: `per_operator::Vector{PerOperatorKraus{T}}`
   - ADD field: `n_jumps::Int` (convenience, = length(per_operator))
   - Keep all other fields: `domain`, `jumps`, `ham_or_trott`, `config`, `precomputed_data`, `delta`, `alpha`, `ws`

3. Refactor `build_trajectoryframework`:
   - Compute `p_jump = 1.0 / length(jumps)` and `delta_eff = delta / p_jump` (= delta * N_jumps) for rate rescaling.
   - Use `alpha = 1 - sqrt(1 - delta_eff)` instead of `alpha = 1 - sqrt(1 - delta)`. This is the key rescaling: each per-operator channel uses an effective delta that is N_jumps times larger, so when you randomly pick one of N_jumps operators each step, the net effect per unit time matches.
   - NOTE: `delta_eff` may be > 1 if `delta * N_jumps > 1`. Guard: assert `delta_eff < 1` or warn. For our test system (delta=0.01, N_jumps=12), delta_eff=0.12 which is fine.
   - Loop over each jump operator `a` in `1:length(jumps)`:
     a. Compute per-operator `R_a = precompute_R(domain, [jumps[a]], ham_or_trott, config, precomputed_data, scratch)` -- pass a single-element vector `[jumps[a]]`. Copy the result from scratch.R immediately since it gets overwritten: `R_a = copy(scratch.R)`.
     b. Rescale: `R_a .*= (1.0 / p_jump)` -- this matches `jump_weight_scaling = gamma_norm_factor / jump_prob` in `jump_contribution!` when `rescale_by_inv_prob=true`, except `precompute_R` already includes `gamma_norm_factor` so we only multiply by `1/p_jump`.

     WAIT -- re-examine: `precompute_R` already multiplies by `gamma_norm_factor` (it uses `base_prefactor` which includes it). And the DM code's `jump_contribution!` uses `jump_weight_scaling = gamma_norm_factor / jump_prob` when rescaling, replacing `gamma_norm_factor` with `gamma_norm_factor / p_jump`. So the trajectory's per-operator R^a should also be scaled by `1/p_jump` on top of what `precompute_R` gives. This is correct.

     c. Build `K0_a = I - alpha * R_a` (where alpha uses delta_eff as computed above).
     d. Build `S_a = (2*alpha - delta_eff)*R_a - alpha^2 * R_a^2`, PSD-guard via eigendecomposition (same as current code), get `U_residual_a`.
     e. For the coherent part: use `precompute_coherent_unitary_terms(jumps, hamiltonian_for_coherent, config, precomputed_data; trotter=trotter_if_needed, delta_scale=1.0/p_jump)` to get the per-jump U_B vector. This function already exists and supports `delta_scale`. Call it ONCE before the loop, get back `coherent_unitaries::Union{Nothing, Vector{Matrix{ComplexF64}}}`. Then `U_B_a = coherent_unitaries === nothing ? nothing : coherent_unitaries[a]`.

     IMPORTANT: For the coherent call we need the original `hamiltonian` (HamHam), not `ham_or_trott` which could be TrottTrott. Looking at `precompute_coherent_unitary_terms` signature, it takes `hamiltonian::HamHam` as positional arg and `trotter` as keyword. But in `build_trajectoryframework` we receive `ham_or_trott::Union{HamHam, TrottTrott}`. The current code uses `precompute_coherent_total_B(jumps, ham_or_trott, ...)` which accepts both. But `precompute_coherent_unitary_terms` needs a `HamHam`. Resolution: for TrotterDomain, `ham_or_trott` is a TrottTrott, and `precompute_coherent_unitary_terms` needs both. We need to pass the hamiltonian separately. BUT `build_trajectoryframework` doesn't receive the hamiltonian when TrotterDomain is used. Looking at `run_trajectories`, it has both `hamiltonian` and `trotter` -- but only passes `ham_or_trott` to `build_trajectoryframework`.

     SIMPLER APPROACH for coherent: Instead of calling `precompute_coherent_unitary_terms` (which needs HamHam), compute per-operator B and U_B inline using the pattern already in `build_trajectoryframework`:
     - If `config.with_coherent`:
       - For each jump `a`, compute `B_a` using the SAME approach as `precompute_coherent_total_B` but passing `[jumps[a]]` instead of `jumps`. Then `U_B_a = exp(-1im * delta_eff * Hermitian(0.5 * (B_a + B_a')))`.
       - Actually, `precompute_coherent_total_B(jumps, ...)` sums over all jumps internally. Passing `[jumps[a]]` gives the single-jump B^a. Then scale by `1/p_jump` (= N_jumps) because the DM code scales delta by `1/p_jump`. But wait, `precompute_coherent_total_B` returns `B * gamma_norm_factor`. The DM code's `coherent_unitaries` use `delta_scale = 1/p_jump`. So `U_B_a = exp(-1im * (delta/p_jump) * Hermitian(B_a))`.
       - Implementation: call `precompute_coherent_total_B([jumps[a]], ham_or_trott, config, precomputed_data)` to get B_a (already scaled by gamma_norm_factor). Then `U_B_a = exp(-1im * delta_eff * Hermitian(0.5 * (B_a + B_a')))`.

     f. Construct `PerOperatorKraus(R_a, K0_a, U_residual_a, U_B_a)`, push to vector.

   - Construct `TrajectoryFramework` with the `per_operator` vector and `n_jumps`.
   - The `delta` stored in framework is still the ORIGINAL delta (not delta_eff), since it controls how many steps to take for a given total_time. The `alpha` stored should be the per-operator alpha (from delta_eff) since it's used in step_along_trajectory!. Actually, store `delta` as the original and `alpha` from delta_eff. The `delta` in the framework determines num_steps in `_evolve_along_trajectory!`. The per-operator channel uses `delta_eff` for the Kraus weights, but that's already baked into K0_a, U_residual_a. The only place `delta` appears in `step_along_trajectory!` is in `p_jump_total = delta * expR`. This should become `delta_eff * expR` for the per-operator version. So also store `delta_eff::Float64` in the framework.

   REVISED TrajectoryFramework fields:
   ```julia
   struct TrajectoryFramework{T,C,H,PD,D<:AbstractDomain}
       domain::D
       jumps::Vector{JumpOp}
       ham_or_trott::H
       config::C
       precomputed_data::PD
       per_operator::Vector{PerOperatorKraus{T}}
       n_jumps::Int
       delta::Float64          # original delta (for time stepping)
       delta_eff::Float64      # delta / p_jump = delta * n_jumps (for Kraus probabilities)
       alpha::Float64          # 1 - sqrt(1 - delta_eff)
       ws::TrajectoryWorkspace{T}
   end
   ```
  </action>
  <verify>
  `cd /Users/bence/code/QuantumFurnace.jl && julia --project -e "using QuantumFurnace; println(\"Module loads OK\")"` compiles without error. The struct definitions and build function are syntactically valid.
  </verify>
  <done>
  `PerOperatorKraus` struct defined. `TrajectoryFramework` holds `per_operator::Vector{PerOperatorKraus}` instead of single R/K0/U_residual/B/U_B. `build_trajectoryframework` loops over jumps computing per-operator Kraus data with 1/p_jump rescaling. Module compiles.
  </done>
</task>

<task type="auto">
  <name>Task 2: Refactor step_along_trajectory! for per-operator branching (both variants)</name>
  <files>src/trajectories.jl</files>
  <action>
Refactor BOTH `step_along_trajectory!` methods (EnergyDomain variant and TimeDomain/TrotterDomain variant) to use per-operator random selection. The two variants differ only in how they build A_{a,omega} (oft! vs NUFFT prefactors).

For both variants, the new logic is:

1. **Random operator selection**: `a = rand(1:fw.n_jumps)`, then `per_op = fw.per_operator[a]`, `jump = fw.jumps[a]`.

2. **Coherent rotation** (if present): Apply `per_op.U_B` (per-operator, not total):
   ```julia
   if per_op.U_B !== nothing
       mul!(ws.psi_tmp, per_op.U_B, psi)
       copyto!(psi, ws.psi_tmp)
       rmul!(psi, 1.0 / sqrt(max(_norm2(psi), eps(Float64))))
   end
   ```

3. **Probabilities from per-operator R_a, K0_a, U_residual_a**:
   ```julia
   mul!(ws.Rpsi, per_op.R, psi)                    # R_a * psi
   expR = max(real(dot(psi, ws.Rpsi)), 0.0)         # <psi|R_a|psi>

   # K0_a * psi = psi - alpha * R_a * psi (but K0_a is precomputed, use it directly)
   mul!(ws.psi_tmp, per_op.K0, psi)                 # K0_a * psi
   p_nojump = _norm2(ws.psi_tmp)

   p_jump_total = fw.delta_eff * expR               # delta_eff, NOT delta

   mul!(ws.Rpsi, per_op.U_residual, psi)            # U_res_a * psi (reuse buffer)
   p_res = _norm2(ws.Rpsi)
   ```
   Note: use `mul!(ws.psi_tmp, per_op.K0, psi)` directly instead of building K0*psi manually, since K0_a is precomputed as a dense matrix. This is cleaner and avoids needing separate alpha*Rpsi computation.

4. **Normalization check and random draw**: Same as current code (check total_weight ~ 1, warn if off, draw r).

5. **Branch: no-jump / residual / dissipative jump**: Same 3-way structure as current, but:
   - No-jump: `psi <- ws.psi_tmp / ||...||` (K0_a * psi, already computed)
   - Residual: `psi <- ws.Rpsi / ||...||` (U_res_a * psi, already computed)
   - Jump: iterate only over frequencies omega for the SINGLE selected jump `a` (not all jumps). The inner loop iterates `energy_labels` for `jump` only (no outer loop over `fw.jumps`).

For the **EnergyDomain jump branch**: use `oft!(ws.jump_oft, jump, w, ham, cfg.sigma)` -- same as current but only for the selected `jump`.

For the **TimeDomain/TrotterDomain jump branch**: use `prefactor_view(oft_prefactors, w)` and `@. ws.jump_oft = jump.in_eigenbasis * pref` -- same as current but only for the selected `jump`.

In the jump probability computation, the rate weight is `fw.delta_eff * base_prefactor * transition(w) / p_jump` -- NO WAIT. The `base_prefactor` already includes `gamma_norm_factor` from `precomputed_data`. The per-operator R_a was rescaled by `1/p_jump`. For consistency, the jump probability weights in the dissipative branch should also use `1/p_jump` scaling. The current formula is `p = delta * rate2(w) * ||A_w * psi||^2`. With rescaling: `p = delta_eff * (base_prefactor / p_jump) * transition(w) * ||A_w * psi||^2`. But `base_prefactor` includes `gamma_norm_factor`. And `per_op.R` was built as `precompute_R([jump]) * (1/p_jump)`. The `precompute_R` uses `base_prefactor * transition(w)` to weight each term. So `R_a = (1/p_jump) * sum_w base_prefactor * transition(w) * A_w' * A_w`. The jump probabilities should be `delta_eff * (base_prefactor / p_jump) * transition(w) * ||A_w psi||^2` to be consistent with `p_jump_total = delta_eff * <psi|R_a|psi>`.

Simplification: define `scaled_prefactor = base_prefactor / p_jump` (= `base_prefactor * n_jumps`). Then each individual jump outcome probability is `delta_eff * scaled_prefactor * transition(w) * n2`.

6. **Keep the same fallback logic** if cumulative scan doesn't reach target due to rounding.

7. Keep `_evolve_along_trajectory!` unchanged -- it uses `fw.delta` for num_steps which is still the original delta. Each step randomly picks an operator, which is the Lie-Trotter splitting.

8. Keep `run_trajectories` unchanged (it just calls `build_trajectoryframework` and `_evolve_along_trajectory!`).
  </action>
  <verify>
  `cd /Users/bence/code/QuantumFurnace.jl && julia --project -e "using QuantumFurnace; println(\"Module loads OK\")"` compiles without error. Then run a quick smoke test:
  ```julia
  julia --project -e "
    using QuantumFurnace, LinearAlgebra, BSON
    ham = BSON.load(\"hamiltonians/heis_disordered_periodic_n3.bson\")[:hamiltonian]
    ham = finalize_hamham(ham, 10.0)
    jumps = [JumpOp(Matrix{ComplexF64}(I,8,8), Matrix{ComplexF64}(I,8,8), true, true)]
    cfg = ThermalizeConfig(num_qubits=3, with_coherent=false, with_linear_combination=true, domain=EnergyDomain(), beta=10.0, sigma=0.1, a=1/3, b=0.4, num_energy_bits=12, w0=0.05, t0=2pi/(4096*0.05), num_trotter_steps_per_t0=10, mixing_time=1.0, delta=0.01)
    pd = precompute_data(cfg.domain, cfg, ham)
    sc = KrausScratch(ComplexF64, 8)
    fw = build_trajectoryframework(jumps, ham, cfg, pd, sc, 0.01)
    psi = zeros(ComplexF64, 8); psi[1] = 1.0
    step_along_trajectory!(psi, fw)
    println(\"Step OK, norm = \", norm(psi))
  "
  ```
  </verify>
  <done>
  Both `step_along_trajectory!` methods (EnergyDomain and TimeDomain/TrotterDomain) randomly select one operator per step, use that operator's precomputed R_a/K0_a/U_residual_a/U_B_a, and iterate only that operator's frequencies for jump sampling. Rate rescaling by 1/p_jump is consistent throughout.
  </done>
</task>

<task type="auto">
  <name>Task 3: Update CPTP test for per-operator channel verification</name>
  <files>test/test_cptp.jl</files>
  <action>
Update `test/test_cptp.jl` to verify per-operator CPTP completeness. The old test checked `K0'*K0 + delta*R + U_res'*U_res = I` for the single combined channel. Now each per-operator channel must satisfy: `K0_a'*K0_a + delta_eff*R_a + U_res_a'*U_res_a = I`.

1. For each domain (EnergyDomain, TimeDomain, TrotterDomain):
   - Build the framework as before (same `build_trajectoryframework` call).
   - Loop over `fw.per_operator`: for each `per_op`:
     ```julia
     completeness = per_op.K0' * per_op.K0 + fw.delta_eff * per_op.R + per_op.U_residual' * per_op.U_residual
     @test isapprox(completeness, identity; atol=1e-10)
     ```
   - Test that `fw.n_jumps == length(TEST_JUMPS)` (= 12 for 4-qubit system).

2. Keep the same tolerance (1e-10) per user decision from Phase 2.

3. Update the testset name to reflect per-operator structure, e.g., "CPTP Per-Operator Completeness (TVAL-01)".
  </action>
  <verify>
  Run the CPTP test: `cd /Users/bence/code/QuantumFurnace.jl && julia --project -e "using Pkg; Pkg.test()" 2>&1 | tail -40`. All tests pass, including the updated CPTP test.
  </verify>
  <done>
  CPTP test verifies per-operator completeness (`K0_a'*K0_a + delta_eff*R_a + U_res_a'*U_res_a = I`) for every operator in all three domains. Full test suite passes.
  </done>
</task>

</tasks>

<verification>
1. Module compiles: `julia --project -e "using QuantumFurnace"`
2. All existing tests pass: `julia --project -e "using Pkg; Pkg.test()"`
3. Per-operator CPTP completeness verified for all three domains
4. `TrajectoryFramework` has `per_operator::Vector{PerOperatorKraus}` field
5. `step_along_trajectory!` uses `rand(1:fw.n_jumps)` for operator selection
6. Rates are rescaled by `1/p_jump` (= N_jumps) consistently in R_a, K0_a, U_residual_a, U_B_a, and jump probabilities
</verification>

<success_criteria>
- PerOperatorKraus struct exists with R, K0, U_residual, U_B fields
- TrajectoryFramework stores per_operator vector instead of single combined operators
- build_trajectoryframework computes per-operator Kraus data with 1/p_jump rescaling
- step_along_trajectory! (both variants) randomly selects one operator per step
- CPTP test passes for each per-operator channel individually
- Full test suite passes (compilation, trajectory fixes, DM tests, scaling)
</success_criteria>

<output>
After completion, create `.planning/quick/7-refactor-step-along-trajectory-to-per-op/7-SUMMARY.md`
</output>
