---
phase: quick
plan: 22
type: execute
autonomous: true

must_haves:
  truths:
    - "The trajectory per-step CPTP channel uses bare delta (NOT delta_eff = delta * n_jumps) for alpha, K0, S, and jump probabilities — matching the DM run_thermalization approach"
    - "R_a scaling by 1/p_jump = n_jumps is kept (this is the correct single compensation)"
    - "Coherent U_B still uses delta * n_jumps (matching DM coherent_unitaries scaling)"
    - "CPTP completeness tests pass with the updated delta"
    - "All existing tests pass"
---

<objective>
Fix double-counting bug in trajectory simulator: the per-operator CPTP channel scales both R by n_jumps AND uses delta_eff = delta * n_jumps, causing the trajectory to evolve n_jumps times faster than the Lindbladian. The DM simulator correctly uses only R scaling with bare delta. This fix makes the trajectory match.

Root cause: In TrajectoryFramework constructor (trajectories.jl), `alpha = 1 - sqrt(1 - delta_eff)` uses `delta_eff = delta * n_jumps` while R_a is already scaled by n_jumps. The DM's `_finalize_kraus_step!` uses `alpha = 1 - sqrt(1 - config.delta)` with the same R scaling. The trajectory should match.

Verified numerically: with the fix, the per-step map spectral gap matches the Lindbladian (ratio 1.0012 vs 12.17 with the bug).
</objective>

<context>
@src/trajectories.jl - TrajectoryFramework constructor and step functions
@src/jump_workers.jl - DM _finalize_kraus_step! (reference for correct approach)
@test/test_cptp.jl - CPTP completeness tests (uses fw.delta_eff, needs update)
@test/test_trajectory_fixes.jl - Trajectory fix tests (uses fw.delta_eff)
@test/test_gns_trajectory.jl - GNS trajectory tests (uses fw.delta_eff)
</context>

<tasks>

<task type="auto">
  <name>Task 1: Fix TrajectoryFramework constructor to use bare delta for CPTP channel</name>
  <files>src/trajectories.jl</files>
  <action>
In the `TrajectoryFramework` constructor function (starts around line 108), make these changes:

1. **Line 92 (struct comment)**: Change the `delta_eff` field comment from:
   ```
   delta_eff::Float64       # delta / p_jump = delta * n_jumps (for per-operator Kraus probabilities)
   ```
   to:
   ```
   delta_eff::Float64       # CPTP channel time parameter (equals delta; R scaling handles 1/p_jump compensation)
   ```

2. **Line 93 (struct comment)**: Change alpha comment from:
   ```
   alpha::Float64           # α = 1 - sqrt(1-δ_eff)
   ```
   to:
   ```
   alpha::Float64           # α = 1 - sqrt(1-δ)
   ```

3. **Lines 119-124**: Change the delta_eff/alpha computation. Replace:
   ```julia
   p_jump = 1.0 / n_jumps
   delta_eff = delta / p_jump   # = delta * n_jumps

   @assert delta_eff < 1.0 "delta_eff = $(delta_eff) >= 1.0: delta=$(delta) * n_jumps=$(n_jumps) is too large for per-operator splitting"

   alpha = 1 - sqrt(1 - delta_eff)
   ```
   with:
   ```julia
   p_jump = 1.0 / n_jumps

   # The per-operator CPTP channel uses bare delta (NOT delta*n_jumps).
   # R_a is already scaled by 1/p_jump = n_jumps to compensate for random
   # operator selection. Using delta*n_jumps would double-count.
   # (The DM _finalize_kraus_step! uses the same approach: scaled R, bare delta.)
   @assert delta < 1.0 "delta = $(delta) >= 1.0: too large for CPTP channel"

   alpha = 1 - sqrt(1 - delta)
   ```

4. **Line 140 (coherent U_B)**: Keep the n_jumps scaling for coherent term but use a local variable. Replace:
   ```julia
   per_op_U_B[a] = exp(-1im * delta_eff * Hermitian(B_a))
   ```
   with:
   ```julia
   # Coherent uses delta/p_jump = delta*n_jumps (matching DM coherent_unitaries scaling)
   per_op_U_B[a] = exp(-1im * (delta / p_jump) * Hermitian(B_a))
   ```

5. **Line 162-164 (S computation)**: Change delta_eff to delta. Replace:
   ```julia
   # S_a = (2*alpha - delta_eff)*R_a - alpha^2 * R_a^2
   mul!(scratch.tmp1, R_a, R_a)   # tmp1 := R_a^2
   s1 = 2 * alpha - delta_eff
   ```
   with:
   ```julia
   # S_a = (2*alpha - delta)*R_a - alpha^2 * R_a^2
   mul!(scratch.tmp1, R_a, R_a)   # tmp1 := R_a^2
   s1 = 2 * alpha - delta
   ```

6. **Line 207 (struct construction)**: Change `delta_eff` value stored in struct from `delta_eff` to `delta`:
   Replace:
   ```julia
   Float64(delta_eff),
   ```
   with:
   ```julia
   Float64(delta),       # delta_eff field = delta (R scaling handles 1/p_jump)
   ```

7. **Lines 863-868 (step docstring)**: Update the docstring comment. Replace:
   ```
   # Per-operator channel structure (Chen 2023, adapted for Lie-Trotter splitting):
   #   Pick a ∈ {1,...,N_jumps} uniformly at random
   #   K0_a = I - alpha*R_a, where alpha = 1 - sqrt(1-delta_eff), delta_eff = delta*N_jumps
   #   K_{a,w} = sqrt(delta_eff * scaled_rate(w)) * L_{a,w}  (jump operators for operator a)
   #   U_res_a: U_res_a'*U_res_a = S_a  (residual for operator a)
   # Rates rescaled by 1/p_jump so net effect per unit time matches DM run_thermalization.
   ```
   with:
   ```
   # Per-operator channel structure (Chen 2023, adapted for Lie-Trotter splitting):
   #   Pick a ∈ {1,...,N_jumps} uniformly at random
   #   K0_a = I - alpha*R_a, where alpha = 1 - sqrt(1-delta), R_a scaled by n_jumps
   #   K_{a,w} = sqrt(delta * scaled_rate(w)) * L_{a,w}  (jump operators for operator a)
   #   U_res_a: U_res_a'*U_res_a = S_a  (residual for operator a)
   # R_a rates rescaled by 1/p_jump; CPTP channel uses bare delta (matching DM).
   ```

8. **Line 878 (TimeDomain step function)**: Change `delta_eff` variable to use `fw.delta` instead. Replace:
   ```julia
   delta_eff = fw.delta_eff
   ```
   with:
   ```julia
   delta = fw.delta
   ```
   Then find-and-replace ALL occurrences of `delta_eff` within this function body (up to the `end` at line 999) with `delta`. These are at lines:
   - 907: `p_jump_total = delta_eff * expR` → `p_jump_total = delta * expR`
   - 951: `p = delta_eff * (scaled_prefactor * transition(w)) * n2` → `p = delta * (scaled_prefactor * transition(w)) * n2`
   - 965: `p = delta_eff * (scaled_prefactor * transition(-w)) * n2` → `p = delta * ...`
   - 983: `p = delta_eff * (scaled_prefactor * transition(w)) * n2` → `p = delta * ...`

9. **Line 1016 (EnergyDomain step function)**: Same change. Replace:
   ```julia
   delta_eff = fw.delta_eff
   ```
   with:
   ```julia
   delta = fw.delta
   ```
   Then replace all `delta_eff` in this function body (up to `end` at line 1136) with `delta`. These are at lines:
   - 1045: `p_jump_total = delta_eff * expR` → `delta * expR`
   - 1089: `p = delta_eff * ...` → `p = delta * ...`
   - 1102: `p = delta_eff * ...` → `p = delta * ...`
   - 1119: `p = delta_eff * ...` → `p = delta * ...`
  </action>
  <verify>
Run `julia --project=. -e "using QuantumFurnace; println(\"OK\")"` to verify compilation.
  </verify>
  <done>
- TrajectoryFramework constructor uses `alpha = 1 - sqrt(1 - delta)` (bare delta)
- R_a scaling by n_jumps is preserved
- Coherent U_B still uses delta * n_jumps
- Both step functions use bare delta for jump probabilities
  </done>
</task>

<task type="auto">
  <name>Task 2: Update CPTP completeness tests and other tests referencing delta_eff</name>
  <files>test/test_cptp.jl, test/test_trajectory_fixes.jl, test/test_gns_trajectory.jl</files>
  <action>
The CPTP completeness identity `K0†K0 + δ·R + S = I` uses whatever δ was used to construct K0 and S. Since we now use bare delta, the tests should use `fw.delta` (which equals `fw.delta_eff` after the fix, but using `fw.delta` is semantically correct).

1. **test/test_cptp.jl**: Update comments (lines 5-13) to reference delta instead of delta_eff. Change all three completeness checks (lines 30, 46, 62) from:
   ```julia
   completeness = per_op.K0' * per_op.K0 + fw.delta_eff * per_op.R + per_op.U_residual' * per_op.U_residual
   ```
   to:
   ```julia
   completeness = per_op.K0' * per_op.K0 + fw.delta * per_op.R + per_op.U_residual' * per_op.U_residual
   ```

2. **test/test_trajectory_fixes.jl**: Update lines 59, 61, and 113 similarly:
   - Line 61: `fw.delta_eff` → `fw.delta`
   - Line 113: `fw.delta_eff` → `fw.delta`

3. **test/test_gns_trajectory.jl**: Update line 59:
   - `fw.delta_eff` → `fw.delta`
  </action>
  <verify>
Run the CPTP tests: `julia --project=. -e 'using Test; include("test/runtests.jl")' 2>&1 | head -50` or run individual test files to verify they pass.
  </verify>
  <done>
- All CPTP completeness tests use fw.delta
- All tests pass
  </done>
</task>

<task type="auto">
  <name>Task 3: Run full test suite and fix any failures</name>
  <files>test/</files>
  <action>
Run the full test suite. The trajectory now evolves at the correct (slower) rate, matching the Lindbladian. Tests that check convergence behavior may need their parameters adjusted (more steps or relaxed tolerance) since the trajectory was previously evolving n_jumps times faster.

Run: `julia --project=. test/runtests.jl`

Expected: Most tests should pass. Tests that relied on the (buggy) fast convergence rate might need:
- Increased ntraj or total_time
- Relaxed tolerance on convergence checks
- Updated expected values

Fix any failures by adjusting test parameters, NOT by reverting the fix.
  </action>
  <verify>
Full test suite passes: `julia --project=. test/runtests.jl`
  </verify>
  <done>
- All tests pass with the corrected trajectory dynamics
  </done>
</task>

<task type="auto">
  <name>Task 4: Update validation script to remove n_jumps normalization (no longer needed)</name>
  <files>experiments/validate_gap_estimation.jl</files>
  <action>
With the trajectory fix, the fitted gap should now directly approximate the Lindbladian gap (no n_jumps factor). The n_jumps normalization and residual factor analysis from plan 24-03 are no longer needed.

Update validate_gap_estimation.jl:
1. Remove the normalized analysis section (n_jumps division, residual_factor computation)
2. Revert the pass criterion to check the raw cross_validate_gap result:
   ```julia
   passed = cv.within_ci || cv.relative_error < 0.3
   ```
3. Update the return type back to just CrossValidationResult
4. Update the summary section to show the raw results
5. Update header comments

Run: `julia --project=. experiments/validate_gap_estimation.jl`

The fitted gap should now be close to the exact gap (within CI or relative_error < 0.3).
  </action>
  <verify>
Run `julia --project=. experiments/validate_gap_estimation.jl` and verify:
- Both n=4 and n=6 output PASS
- Fitted gap is close to exact gap (not 20x off)
- OVERALL: PASS
  </verify>
  <done>
- Validation script works with direct comparison (no n_jumps normalization)
- Both systems pass cross-validation
  </done>
</task>

</tasks>

<verification>
1. `julia --project=. -e "using QuantumFurnace"` compiles
2. `julia --project=. test/runtests.jl` all tests pass
3. `julia --project=. experiments/validate_gap_estimation.jl` both systems PASS with OVERALL: PASS
4. No changes to DM simulator (src/jump_workers.jl, src/furnace.jl)
</verification>

<success_criteria>
- Trajectory per-step CPTP channel uses bare delta (matching DM)
- R_a × n_jumps scaling preserved (single compensation)
- All tests pass
- Validation script shows fitted gap ≈ exact gap (no n_jumps factor)
</success_criteria>
