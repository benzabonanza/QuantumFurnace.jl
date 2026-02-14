---
phase: quick-8
plan: 01
type: execute
wave: 1
depends_on: []
files_modified:
  - src/furnace.jl
  - src/trajectories.jl
  - test/trajectory_validation/run_trajectory_validation.jl
autonomous: true
must_haves:
  truths:
    - "TrotterDomain Liouvillian fixed point is within ~1e-6 of Gibbs (not ~0.004)"
    - "TrotterDomain trajectory convergence threshold can be tightened from 0.02 to match reduced domain approx error"
    - "All existing tests still pass (no regressions)"
  artifacts:
    - path: "src/furnace.jl"
      provides: "Trotter-basis-transformed jumps for construct_lindbladian and run_thermalization"
      contains: "trafo_from_eigen_to_trotter"
    - path: "src/trajectories.jl"
      provides: "Trotter-basis-transformed jumps for precompute_R and step_along_trajectory"
      contains: "trafo_from_eigen_to_trotter"
  key_links:
    - from: "src/furnace.jl"
      to: "src/jump_workers.jl"
      via: "jump_contribution! receives jump.in_eigenbasis already in correct basis"
      pattern: "jumps_for_diss|jump_trott"
    - from: "src/trajectories.jl"
      to: "src/jump_workers.jl"
      via: "precompute_R receives transformed jumps"
      pattern: "jumps_for_diss|jump_trott"
---

<objective>
Fix TrotterDomain Gibbs fixed point distance from ~0.004 to ~1e-6 by correcting
a basis mismatch in the dissipative Liouvillian construction and trajectory code.

Purpose: The OFT identity A(w)_{ij} = A_{ij} * P_{ij}(w) requires A to be in the
same eigenbasis as the frequencies used to compute P. For TrotterDomain, P uses
Trotter quasi-Bohr frequencies, so A must be in the Trotter eigenbasis. Currently
jump.in_eigenbasis (Hamiltonian eigenbasis) is used directly, causing a basis mismatch.

Output: Corrected furnace.jl and trajectories.jl with Trotter-transformed jumps;
tightened test thresholds confirming the fix.
</objective>

<execution_context>
@./.claude/get-shit-done/workflows/execute-plan.md
@./.claude/get-shit-done/templates/summary.md
</execution_context>

<context>
@.planning/PROJECT.md
@.planning/ROADMAP.md
@.planning/STATE.md
@src/furnace.jl
@src/trajectories.jl
@src/jump_workers.jl
@src/coherent.jl (lines 220-260 for B_trotter basis transform pattern)
@src/trotter_domain.jl (TrottTrott struct with trafo_from_eigen_to_trotter field)
@src/structs.jl (JumpOp struct definition)
@test/trajectory_validation/run_trajectory_validation.jl
@test/test_dm_detailed_balance.jl
</context>

<tasks>

<task type="auto">
  <name>Task 1: Transform jump operators to Trotter eigenbasis in all entry points</name>
  <files>src/furnace.jl, src/trajectories.jl</files>
  <action>
The bug: `jump.in_eigenbasis` is in Hamiltonian eigenbasis, but for TrotterDomain the NUFFT
prefactors use Trotter quasi-Bohr frequencies. The element-wise product `A .* P` only holds
when A is in the same basis as the frequencies defining P.

The fix pattern (already used in coherent.jl:221 for B_trotter):
```julia
U = trotter.trafo_from_eigen_to_trotter
jump_in_trotter = U * jump.in_eigenbasis * U'
```

Apply this transformation at three entry points, creating a `jumps_for_diss` vector
that replaces `jumps` in downstream code for TrotterDomain:

**1. `construct_lindbladian` in furnace.jl (line ~72):**
Before the jump loop (line 72), when `config.domain isa TrotterDomain`, create:
```julia
jumps_for_diss = if config.domain isa TrotterDomain
    U = trotter.trafo_from_eigen_to_trotter
    [JumpOp(j.data, U * j.in_eigenbasis * U', j.orthogonal, j.hermitian) for j in jumps]
else
    jumps
end
```
Then use `jumps_for_diss` in the for-loop on line 72 (`for (k, jump) in pairs(jumps_for_diss)`)
and also pass `jumps_for_diss` to `precompute_coherent_total_B` on line 66. Note: the coherent
B computation in `precompute_coherent_total_B` already handles its own basis transform internally
for TrotterDomain (via B_trotter in coherent.jl:221), so passing transformed jumps would
double-transform. Therefore ONLY change the dissipative loop (line 72), NOT line 66.

**2. `run_thermalization` in furnace.jl (line ~107-132):**
After the `ham_or_trott` assignment block (line ~101), create `jumps_for_diss`:
```julia
jumps_for_diss = if config.domain isa TrotterDomain
    U = trotter.trafo_from_eigen_to_trotter
    [JumpOp(j.data, U * j.in_eigenbasis * U', j.orthogonal, j.hermitian) for j in jumps]
else
    jumps
end
```
Then use `jumps_for_diss` in `jump_contribution!` call (line 122, the `jump = jumps[idx]` on
line 120 should become `jump = jumps_for_diss[idx]`). Keep using original `jumps` for
`precompute_coherent_unitary_terms` (line 108) -- that function handles its own Trotter basis
transform internally.

**3. `build_trajectoryframework` in trajectories.jl (line 49-125):**
After line 57 (`dim = size(jumps[1].data, 1)`), determine if domain is TrotterDomain and
create transformed jumps:
```julia
jumps_for_diss = if config.domain isa TrotterDomain && ham_or_trott isa TrottTrott
    U = ham_or_trott.trafo_from_eigen_to_trotter
    [JumpOp(j.data, U * j.in_eigenbasis * U', j.orthogonal, j.hermitian) for j in jumps]
else
    collect(jumps)
end
```
Use `jumps_for_diss` in:
- `precompute_R` call on line 84: `precompute_R(config.domain, [jumps_for_diss[a]], ...)`
- `precompute_coherent_total_B` on line 71 should stay with original `jumps[a]` (handles own transform)
- The `collect(jumps)` stored in the framework at line 114 should store `jumps_for_diss`
  so that `step_along_trajectory!` uses the transformed jumps (since it accesses `fw.jumps[a].in_eigenbasis`)

**4. `precompute_R` for Union{TimeDomain,TrotterDomain} in trajectories.jl (line 197-249):**
No changes needed here -- it receives jumps from callers who will now pass transformed jumps.

**5. `step_along_trajectory!` for Union{TimeDomain,TrotterDomain} in trajectories.jl (line 445-584):**
No changes needed here -- it reads from `fw.jumps[a]` which will now contain transformed jumps
since `build_trajectoryframework` stores `jumps_for_diss`.

IMPORTANT: Do NOT touch `precompute_coherent_total_B` or `precompute_coherent_unitary_terms`
calls -- these already handle TrotterDomain basis transform internally via `B_trotter()` in
coherent.jl. Transforming jumps before passing to these would double-transform.
  </action>
  <verify>
Run the full test suite:
```bash
cd /Users/bence/code/QuantumFurnace.jl && julia --project -e 'using Pkg; Pkg.test()'
```
All tests must pass including DMTST-02 (domain hierarchy). Check the DMTST-02 log output
to confirm TrotterDomain Gibbs distance has dropped from ~0.004 to ~1e-6 or smaller.
  </verify>
  <done>
TrotterDomain Liouvillian fixed point trace distance to Gibbs is ~1e-6 or better (was ~0.004).
All DMTST tests pass. run_thermalization and trajectory code use Trotter-basis-transformed
jump operators for the dissipative contribution.
  </done>
</task>

<task type="auto">
  <name>Task 2: Tighten TrotterDomain test thresholds</name>
  <files>test/trajectory_validation/run_trajectory_validation.jl</files>
  <action>
After the fix in Task 1, TrotterDomain's domain approximation error drops from ~0.004 to ~1e-6.
Update test thresholds that were artificially inflated to accommodate the basis mismatch:

1. **Fixed point -> Gibbs check (line 197):**
   Change `@test fp_gibbs_dist < 0.01` to `@test fp_gibbs_dist < 1e-4`
   (allowing generous margin above expected ~1e-6)

2. **Trajectory -> Gibbs check (line 239):**
   The comment on line 236 says "domain approximation offset (~0.005)". Update the comment
   to reflect that domain approximation is now ~1e-6, so the dominant error source is purely
   statistical noise from finite trajectory count (~0.01 for 10k trajectories).
   Change `@test dist_traj_gibbs < 0.02` to `@test dist_traj_gibbs < 0.015`
   (statistical noise ~0.01 is now the dominant term; 0.015 provides margin)

3. **Update comment on line 182:**
   Change "offset (~0.005)" to "offset (~1e-6)" to reflect corrected domain approximation error.

4. **Update comment on line 236-237:**
   Reflect that domain approx is now negligible, total expected distance is ~0.01 (statistical only).
  </action>
  <verify>
Run the trajectory validation tests:
```bash
cd /Users/bence/code/QuantumFurnace.jl && julia --project -e '
    include("test/test_helpers.jl")
    include("test/trajectory_validation/run_trajectory_validation.jl")
'
```
All trajectory validation tests pass with tightened thresholds.
  </verify>
  <done>
TrotterDomain test thresholds reflect the corrected domain approximation error (~1e-6 instead
of ~0.004). Fixed point test uses 1e-4 threshold. Trajectory test uses 0.015 threshold.
All comments updated to reflect actual error budget.
  </done>
</task>

</tasks>

<verification>
1. Run full test suite: `julia --project -e 'using Pkg; Pkg.test()'` -- all pass
2. Check DMTST-02 log output: TrotterDomain distance to Gibbs should be ~1e-6 (was ~0.004)
3. Check trajectory validation log: fixed point -> Gibbs distance should be ~1e-6
4. Verify no regressions in EnergyDomain or TimeDomain results (no basis transform applied for those)
</verification>

<success_criteria>
- TrotterDomain Liouvillian fixed point trace distance to Gibbs: < 1e-4 (expected ~1e-6)
- All existing test suites pass without regression
- Test thresholds tightened to reflect corrected error budget
- Basis transform applied consistently in construct_lindbladian, run_thermalization, and trajectory framework
</success_criteria>

<output>
After completion, create `.planning/quick/8-fix-trotterdomain-gibbs-fixed-point-dist/8-SUMMARY.md`
</output>
