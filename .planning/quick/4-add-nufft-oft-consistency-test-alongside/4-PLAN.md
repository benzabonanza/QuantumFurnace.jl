---
phase: quick-4
plan: 01
type: execute
wave: 1
depends_on: []
files_modified: [test/test_dm_scaling.jl]
autonomous: true

must_haves:
  truths:
    - "NUFFT time OFT matches time_oft! to near machine precision (<1e-10)"
    - "NUFFT trotter OFT matches trotter_oft! to near machine precision (<1e-10)"
    - "NUFFT time OFT matches analytical oft! within TOL_QUADRATURE"
    - "NUFFT trotter OFT matches analytical oft! within 0.1 (Trotter error)"
  artifacts:
    - path: "test/test_dm_scaling.jl"
      provides: "DMTST-06b NUFFT OFT consistency testset"
      contains: "DMTST-06b"
  key_links:
    - from: "test/test_dm_scaling.jl"
      to: "src/nufft.jl"
      via: "QuantumFurnace.prepare_oft_nufft_prefactors, QuantumFurnace.prefactor_view"
      pattern: "prefactor_view"
    - from: "test/test_dm_scaling.jl"
      to: "src/furnace_utensils.jl"
      via: "precompute_data(TimeDomain(), ...) and precompute_data(TrotterDomain(), ...)"
      pattern: "precompute_data"
---

<objective>
Add a new testset "DMTST-06b: NUFFT OFT consistency" to test/test_dm_scaling.jl that verifies the NUFFT-based OFT method gives results consistent with the analytical OFT and the existing time_oft!/trotter_oft! methods.

Purpose: Ensure the NUFFT acceleration of OFT computation is mathematically faithful to the direct summation methods it replaces.
Output: New passing testset in test/test_dm_scaling.jl
</objective>

<execution_context>
@./.claude/get-shit-done/workflows/execute-plan.md
@./.claude/get-shit-done/templates/summary.md
</execution_context>

<context>
@test/test_dm_scaling.jl
@test/test_helpers.jl
@src/nufft.jl
@src/furnace_utensils.jl (lines 82-133: precompute_data for TimeDomain/TrotterDomain)
@src/ofts.jl (time_oft!, trotter_oft!, oft!)
</context>

<tasks>

<task type="auto">
  <name>Task 1: Add DMTST-06b NUFFT OFT consistency testset</name>
  <files>test/test_dm_scaling.jl</files>
  <action>
Append a new `@testset "DMTST-06b: NUFFT OFT consistency"` block after the existing DMTST-06 testset (after line 185) in test/test_dm_scaling.jl. The testset should:

1. **Setup (shared with existing DMTST-06 pattern):**
   - `jump = TEST_JUMPS[1]`
   - `w = -3 * W0` (same test energy as DMTST-06)
   - Compute `energy_oft_prefactor = 1 / sqrt(SIGMA * sqrt(2 * pi))`
   - Compute analytical `A_energy` via `oft!(A_energy, jump, w, TEST_HAM, SIGMA)` scaled by `energy_oft_prefactor`
   - Reconstruct `energy_labels`, `time_labels_full`, `oft_time_labels` (same as DMTST-06 lines 151-153)
   - Compute `time_oft_prefactor = T0 * sqrt(SIGMA * sqrt(2 / pi) / (2 * pi))`
   - Create `caches = OFTCaches(DIM)`

2. **Time NUFFT OFT:**
   - Get NUFFT prefactors via `precompute_data(TimeDomain(), config_time, TEST_HAM)` where `config_time = make_liouv_config(TimeDomain())`
   - Extract the NUFFT prefactor matrix: `nufft_pf = QuantumFurnace.prefactor_view(precomputed_time.oft_nufft_prefactors, w)`
   - Compute `A_nufft_time = jump.in_eigenbasis .* nufft_pf`
   - Scale: `A_nufft_time .*= time_oft_prefactor`
   - Also compute `A_time` via `QuantumFurnace.time_oft!(A_time, caches, jump, w, TEST_HAM, oft_time_labels, SIGMA)` scaled by `time_oft_prefactor` (same as DMTST-06) for direct comparison

3. **Trotter NUFFT OFT:**
   - Get NUFFT prefactors via `precompute_data(TrotterDomain(), config_trott, TEST_TROTTER)` where `config_trott = make_liouv_config(TrotterDomain())`
   - Transform jump to Trotter eigenbasis: `U = TEST_TROTTER.trafo_from_eigen_to_trotter; jump_trott = JumpOp(jump.data, U * jump.in_eigenbasis * U', jump.orthogonal, jump.hermitian)`
   - Extract NUFFT prefactor matrix: `nufft_pf_trott = QuantumFurnace.prefactor_view(precomputed_trott.oft_nufft_prefactors, w)`
   - Compute `A_nufft_trott = jump_trott.in_eigenbasis .* nufft_pf_trott` (in Trotter eigenbasis)
   - Scale: `A_nufft_trott .*= time_oft_prefactor`
   - Transform back to H-eigenbasis: `A_nufft_trott_in_eigen = U' * A_nufft_trott * U`
   - Also compute `A_trott` via `QuantumFurnace.trotter_oft!(A_trott, caches, jump_trott, w, TEST_TROTTER, oft_time_labels, SIGMA)` scaled by `time_oft_prefactor`, transformed back: `A_trott_in_eigen = U' * A_trott * U` (same as DMTST-06)

4. **Diagnostics (println):**
   - `dist_nufft_time_vs_time = norm(A_nufft_time - A_time)` -- should be ~0 (< 1e-10)
   - `dist_nufft_trott_vs_trott = norm(A_nufft_trott_in_eigen - A_trott_in_eigen)` -- should be ~0 (< 1e-10)
   - `dist_nufft_time_vs_energy = norm(A_nufft_time - A_energy)` -- should match TOL_QUADRATURE
   - `dist_nufft_trott_vs_energy = norm(A_nufft_trott_in_eigen - A_energy)` -- should be < 0.1

5. **Assertions:**
   ```julia
   # NUFFT time OFT matches direct time_oft! (both compute same sum, NUFFT uses FFT with eps=1e-12)
   @test norm(A_nufft_time - A_time) < 1e-10

   # NUFFT trotter OFT matches direct trotter_oft!
   @test norm(A_nufft_trott_in_eigen - A_trott_in_eigen) < 1e-10

   # NUFFT time OFT matches analytical (same tolerance as DMTST-06 time vs energy)
   @test norm(A_nufft_time - A_energy) < TOL_QUADRATURE

   # NUFFT trotter OFT matches analytical within Trotter error bound
   @test norm(A_nufft_trott_in_eigen - A_energy) < 0.1
   ```

**Important details:**
- The `precompute_data` for TimeDomain/TrotterDomain returns `energy_labels` which are the TRUNCATED energy labels. Verify that `w = -3 * W0 = -0.15` exists in the NUFFT prefactors' energy_to_index dict. If not, pick a `w` value that is in the dict. Add an explicit `@test haskey(precomputed_time.oft_nufft_prefactors.energy_to_index, w)` as a sanity check before proceeding.
- For the Trotter case, `jump_trott.in_eigenbasis` is the jump operator expressed in the Trotter eigenbasis (not H-eigenbasis), matching what `trotter_oft!` expects.
- The NUFFT prefactors for TrotterDomain use `TEST_TROTTER.bohr_freqs` (Trotter Bohr frequencies), consistent with `trotter_oft!`.
  </action>
  <verify>
Run the test suite:
```bash
cd /Users/bence/code/QuantumFurnace.jl && julia --project -e 'using Pkg; Pkg.test()'
```
Or run just the scaling test file:
```bash
cd /Users/bence/code/QuantumFurnace.jl && julia --project -e '
  using Test
  include("test/test_helpers.jl")
  include("test/test_dm_scaling.jl")
'
```
All tests pass, including the new DMTST-06b testset. Verify diagnostic output shows:
- NUFFT vs direct distances are < 1e-10
- NUFFT time vs analytical distance is < TOL_QUADRATURE (1e-6)
- NUFFT trotter vs analytical distance is < 0.1
  </verify>
  <done>
DMTST-06b testset exists in test/test_dm_scaling.jl and passes. Four assertions verify:
(1) NUFFT time OFT matches time_oft! < 1e-10,
(2) NUFFT trotter OFT matches trotter_oft! < 1e-10,
(3) NUFFT time OFT matches analytical OFT < TOL_QUADRATURE,
(4) NUFFT trotter OFT matches analytical OFT < 0.1.
  </done>
</task>

</tasks>

<verification>
- `julia --project -e 'using Pkg; Pkg.test()'` passes with no failures
- DMTST-06b testset appears in test output
- Diagnostic printlns confirm expected magnitude of distances
</verification>

<success_criteria>
- New DMTST-06b testset added after DMTST-06 in test/test_dm_scaling.jl
- All 4 NUFFT consistency assertions pass
- Existing DMTST-03 through DMTST-06 tests unaffected
</success_criteria>

<output>
After completion, create `.planning/quick/4-add-nufft-oft-consistency-test-alongside/4-SUMMARY.md`
</output>
