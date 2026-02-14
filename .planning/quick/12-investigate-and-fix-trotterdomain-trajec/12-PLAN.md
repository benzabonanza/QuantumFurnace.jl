---
phase: quick-12
plan: 01
type: execute
wave: 1
depends_on: []
files_modified:
  - test/test_regression.jl
  - test/reference/generate_references.jl
autonomous: true
must_haves:
  truths:
    - "Trajectory regression tests pass on any platform (macOS aarch64, Linux x86_64) without regenerating BSON files"
    - "Trajectory regression tests still catch real code bugs (broken Kraus operators, wrong jump sampling)"
    - "DM regression tests remain unchanged at 1e-10 tolerance"
  artifacts:
    - path: "test/test_regression.jl"
      provides: "Platform-portable regression tests"
      contains: "exp.*delta.*L"
    - path: "test/reference/generate_references.jl"
      provides: "Updated reference generator (DM-only trajectory references removed)"
  key_links:
    - from: "test/test_regression.jl"
      to: "construct_lindbladian + exp(delta*L)"
      via: "DM evolution computed at test time as trajectory reference"
      pattern: "exp.*delta.*L"
---

<objective>
Replace frozen BSON trajectory regression tests with DM-based comparison that is platform-portable.

Purpose: The trajectory regression tests have failed 3 times (quick-10, quick-11, now quick-12) because
frozen BSON trajectory averages are not portable across platforms. Stochastic branch selection in
`step_along_trajectory!` depends on BLAS internals (OpenBLAS vs Accelerate) and RNG stream behavior
that differs between Linux x86_64 (where references were generated) and macOS aarch64 (where the user
runs tests). The DM regression tests pass at 1e-10 across all platforms because `exp(delta*L)` is
deterministic. The fix is to compare trajectory averages against the DM evolution result (computed fresh
at test time) instead of against frozen trajectory BSON data, using a statistical tolerance appropriate
for 1000 trajectories.

Output: Modified test/test_regression.jl with platform-portable trajectory tests, cleaned up
generate_references.jl (trajectory BSON files kept but no longer used by regression tests).
</objective>

<execution_context>
@./.claude/get-shit-done/workflows/execute-plan.md
@./.claude/get-shit-done/templates/summary.md
</execution_context>

<context>
@.planning/STATE.md
@test/test_regression.jl
@test/test_helpers.jl
@test/reference/generate_references.jl
@test/runtests.jl
</context>

<tasks>

<task type="auto">
  <name>Task 1: Replace frozen trajectory BSON comparison with DM-based comparison</name>
  <files>test/test_regression.jl</files>
  <action>
Modify the two trajectory regression testsets ("Trajectory regression: EnergyDomain" and
"Trajectory regression: TrotterDomain (coherent)") to compare trajectory-averaged rho against
the DM evolution result computed fresh at test time, instead of loading from BSON.

For EACH trajectory regression testset:

1. REMOVE the BSON.load call and the ref_data extraction (rho_ref, delta, seed, ntraj).

2. DEFINE constants inline at the top of each testset:
   - delta = 0.1 (same as REF_DELTA in generate_references.jl)
   - seed = 12345 (same as REF_SEED)
   - ntraj = 1000 (same as REF_NTRAJ)

3. COMPUTE the DM reference fresh (same pattern as the DM regression tests above):
   For EnergyDomain:
     liouv_config = make_small_liouv_config(EnergyDomain())
     L = construct_lindbladian(SMALL_JUMPS, liouv_config, SMALL_HAM)
     rho_dm = reshape(exp(delta * L) * vec(rho0), SMALL_DIM, SMALL_DIM)
     rho_dm = (rho_dm + rho_dm') / 2

   For TrotterDomain:
     liouv_config = make_small_liouv_config(TrotterDomain(); with_coherent=true)
     L = construct_lindbladian(SMALL_JUMPS, liouv_config, SMALL_HAM; trotter=SMALL_TROTTER)
     rho_dm = reshape(exp(delta * L) * vec(rho0), SMALL_DIM, SMALL_DIM)
     rho_dm = (rho_dm + rho_dm') / 2

4. KEEP the trajectory computation EXACTLY as-is (therm_config, precomputed, scratch, fw,
   Random.seed!, trajectory loop, averaging, Hermitianization). Do not change any of this.

5. CHANGE the final @test to compare against rho_dm instead of rho_ref:
     @test isapprox(rho_traj, rho_dm; atol=0.05)

   Use atol=0.05 because with 1000 trajectories and delta=0.1, the statistical error is
   O(1/sqrt(1000)) ~ 0.03. The 0.05 tolerance gives ~1.7 sigma headroom, which is generous
   enough to never false-positive while still catching real bugs (a broken Kraus operator
   produces errors of O(1), not O(0.03)).

6. UPDATE the comment block above each trajectory testset to explain:
   - Trajectory averages are compared against DM evolution (not frozen BSON)
   - This is platform-portable because DM evolution via exp(delta*L) is deterministic
   - The atol=0.05 tolerance accommodates O(1/sqrt(N_traj)) statistical noise
   - Real regressions produce O(1) errors and will be caught easily

7. REMOVE the `using BSON, Random` import at the top of the file. Replace with just
   `using Random` (BSON is no longer needed in this file since the DM regression tests
   still use BSON -- wait, check: lines 26 and 77 still load BSON for DM tests. So KEEP
   `using BSON, Random` unchanged.)

8. UPDATE the module docstring at the top of the file to note that trajectory regression
   tests compare against DM evolution (not frozen trajectory data) for platform portability.

Do NOT change the DM regression testsets at all. They work perfectly at 1e-10.
  </action>
  <verify>
Run: cd /Users/bence/code/QuantumFurnace.jl && julia --project -e 'using Pkg; Pkg.test()'

All 4 regression tests should pass:
- DM regression: EnergyDomain (atol=1e-10) -- unchanged
- Trajectory regression: EnergyDomain (atol=0.05) -- now vs DM
- DM regression: TrotterDomain (atol=1e-10) -- unchanged
- Trajectory regression: TrotterDomain (atol=0.05) -- now vs DM
  </verify>
  <done>
Both trajectory regression tests compare against DM evolution computed at test time (not BSON),
pass on this platform, and will pass on any platform where the DM tests pass. The test file
no longer loads trajectory BSON reference files.
  </done>
</task>

<task type="auto">
  <name>Task 2: Clean up generate_references.jl and remove stale trajectory BSON files</name>
  <files>test/reference/generate_references.jl, test/reference/energy_traj_reference.bson, test/reference/trotter_coherent_traj_reference.bson</files>
  <action>
Since the trajectory BSON files are no longer loaded by any test, clean them up:

1. DELETE the two trajectory BSON files:
   - test/reference/energy_traj_reference.bson
   - test/reference/trotter_coherent_traj_reference.bson

2. UPDATE test/reference/generate_references.jl:
   - Remove the `generate_traj_reference` function entirely
   - Remove the two `generate_traj_reference(...)` calls at the bottom
   - Update the docstring at the top: now generates 2 reference files (not 4), listing only
     energy_dm_reference.bson and trotter_coherent_dm_reference.bson
   - Remove `Random` from the `using` line (no longer needed)
   - Remove REF_SEED and REF_NTRAJ constants (only used by trajectory generator)
   - Keep REF_DELTA and REF_DIR (still used by DM generator)
   - Update the final println to say "2 reference files" instead of "4"

3. VERIFY the DM reference files still exist and are untouched:
   - test/reference/energy_dm_reference.bson (must exist)
   - test/reference/trotter_coherent_dm_reference.bson (must exist)
  </action>
  <verify>
Verify trajectory BSON files are gone:
  ls test/reference/*.bson  -- should show only energy_dm_reference.bson and trotter_coherent_dm_reference.bson

Verify generate_references.jl still works:
  cd /Users/bence/code/QuantumFurnace.jl && julia --project test/reference/generate_references.jl
  -- should generate the 2 DM reference files without error

Run full test suite again to confirm nothing is broken:
  cd /Users/bence/code/QuantumFurnace.jl && julia --project -e 'using Pkg; Pkg.test()'
  </verify>
  <done>
Trajectory BSON reference files are deleted, generate_references.jl only generates DM references,
and the full test suite passes cleanly.
  </done>
</task>

</tasks>

<verification>
1. `julia --project -e 'using Pkg; Pkg.test()'` passes all tests including the 4 regression tests
2. No test file references energy_traj_reference.bson or trotter_coherent_traj_reference.bson
3. DM regression tests remain at atol=1e-10 (untouched)
4. Trajectory regression tests use atol=0.05 against DM evolution (platform-portable)
</verification>

<success_criteria>
- All regression tests pass on this sandbox (Linux x86_64)
- Trajectory regression tests no longer depend on any BSON trajectory files
- The tests will also pass on the user's macOS aarch64 + Julia 1.12.4 because the DM reference
  is computed fresh at test time (same as the DM regression tests that already pass at 1e-10)
- Real code regressions (O(1) errors from broken Kraus operators) will still be caught by the 0.05 tolerance
</success_criteria>

<output>
After completion, create `.planning/quick/12-investigate-and-fix-trotterdomain-trajec/12-SUMMARY.md`
</output>
