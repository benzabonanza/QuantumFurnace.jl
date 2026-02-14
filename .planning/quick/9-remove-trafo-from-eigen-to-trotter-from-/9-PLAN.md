---
phase: quick-9
plan: 01
type: execute
wave: 1
depends_on: []
files_modified:
  - src/trotter_domain.jl
  - src/furnace.jl
  - src/coherent.jl
  - src/trajectories.jl
  - src/ofts.jl
  - test/test_dm_scaling.jl
autonomous: true

must_haves:
  truths:
    - "TrottTrott struct has no trafo_from_eigen_to_trotter field"
    - "All TrotterDomain simulations produce the same numerical results as before"
    - "All existing tests pass unchanged (except for the field removal adaptation in test code)"
  artifacts:
    - path: "src/trotter_domain.jl"
      provides: "TrottTrott struct without trafo_from_eigen_to_trotter"
      contains: "mutable struct TrottTrott"
    - path: "src/furnace.jl"
      provides: "Jump basis transforms using trotter.eigvecs directly"
    - path: "src/coherent.jl"
      provides: "B_trotter functions using trotter.eigvecs directly"
    - path: "src/trajectories.jl"
      provides: "TrajectoryFramework jump transforms using trotter.eigvecs directly"
  key_links:
    - from: "src/furnace.jl"
      to: "trotter.eigvecs"
      via: "trotter.eigvecs' * j.data * trotter.eigvecs"
      pattern: "trotter\\.eigvecs.*j\\.data.*trotter\\.eigvecs"
    - from: "src/coherent.jl"
      to: "trotter.eigvecs"
      via: "trotter.eigvecs' * jump.data * trotter.eigvecs"
      pattern: "trotter\\.eigvecs.*jump\\.data.*trotter\\.eigvecs"
---

<objective>
Remove the `trafo_from_eigen_to_trotter` field from the TrottTrott struct and eliminate all
redundant basis transformations that use it.

Purpose: This field stores `trottU_eigvecs' * hamiltonian.eigvecs`, which is a derived quantity
from `trotter.eigvecs` and `hamiltonian.eigvecs` that are already available at all call sites.
The transformation `U * jump.in_eigenbasis * U'` (H-eigen to Trotter-eigen) is mathematically
equivalent to `trotter.eigvecs' * jump.data * trotter.eigvecs` (computational basis to Trotter-eigen
directly). Removing this redundant field simplifies the struct and makes the code clearer.

Output: Cleaned TrottTrott struct, all call sites updated, all tests passing.
</objective>

<execution_context>
@./.claude/get-shit-done/workflows/execute-plan.md
@./.claude/get-shit-done/templates/summary.md
</execution_context>

<context>
@src/trotter_domain.jl (TrottTrott struct definition and constructor)
@src/furnace.jl (construct_lindbladian lines 74-77, run_thermalization lines 118-121)
@src/coherent.jl (B_trotter single-jump lines 221-222, B_trotter multi-jump lines 254+266)
@src/trajectories.jl (build_trajectoryframework lines 69-74)
@src/ofts.jl (comment on line 136)
@test/test_dm_scaling.jl (DMTST-05 line 116, DMTST-06 lines 162-163, DMTST-06b lines 230-231)
</context>

<tasks>

<task type="auto">
  <name>Task 1: Remove field from struct/constructor and update all source call sites</name>
  <files>
    src/trotter_domain.jl
    src/furnace.jl
    src/coherent.jl
    src/trajectories.jl
    src/ofts.jl
  </files>
  <action>
  **1. `src/trotter_domain.jl`:**
  - Remove the `trafo_from_eigen_to_trotter::Matrix{ComplexF64}` field (line 17) from the TrottTrott struct.
  - Remove the docstring line referencing it (line 10).
  - In the constructor `TrottTrott(hamiltonian::HamHam, ...)`:
    - Remove the line `unitary_from_eigen_to_trotter = trottU_eigvecs' * hamiltonian.eigvecs` (line 25).
    - Remove `unitary_from_eigen_to_trotter` from the return TrottTrott(...) call (line 32).
  - The struct should now have 4 fields: `t0`, `num_trotter_steps_per_t0`, `eigvals_t0`, `eigvecs`, `bohr_freqs`.

  **2. `src/furnace.jl` -- `construct_lindbladian` (around line 74):**
  Replace the block:
  ```julia
  jumps_for_diss = if config.domain isa TrotterDomain
      U = trotter.trafo_from_eigen_to_trotter
      JumpOp[JumpOp(j.data, U * j.in_eigenbasis * U', j.orthogonal, j.hermitian) for j in jumps]
  else
      jumps
  end
  ```
  With:
  ```julia
  jumps_for_diss = if config.domain isa TrotterDomain
      JumpOp[JumpOp(j.data, trotter.eigvecs' * j.data * trotter.eigvecs, j.orthogonal, j.hermitian) for j in jumps]
  else
      jumps
  end
  ```

  **3. `src/furnace.jl` -- `run_thermalization` (around line 118):**
  Apply the same replacement pattern:
  ```julia
  jumps_for_diss = if config.domain isa TrotterDomain
      JumpOp[JumpOp(j.data, trotter.eigvecs' * j.data * trotter.eigvecs, j.orthogonal, j.hermitian) for j in jumps]
  else
      jumps
  end
  ```

  **4. `src/coherent.jl` -- `B_trotter` single-jump variant (around line 221):**
  Replace:
  ```julia
  U = trotter.trafo_from_eigen_to_trotter
  jump_in_trotter = U * jump.in_eigenbasis * U'
  ```
  With:
  ```julia
  jump_in_trotter = trotter.eigvecs' * jump.data * trotter.eigvecs
  ```

  **5. `src/coherent.jl` -- `B_trotter` multi-jump variant (around line 254):**
  Remove:
  ```julia
  U = trotter.trafo_from_eigen_to_trotter
  ```
  And replace:
  ```julia
  jump_a_trotter = U * jump_a.in_eigenbasis * U'
  ```
  With:
  ```julia
  jump_a_trotter = trotter.eigvecs' * jump_a.data * trotter.eigvecs
  ```

  **6. `src/trajectories.jl` -- `build_trajectoryframework` (around line 69):**
  Replace:
  ```julia
  jumps_for_diss = if config.domain isa TrotterDomain && ham_or_trott isa TrottTrott
      U = ham_or_trott.trafo_from_eigen_to_trotter
      JumpOp[JumpOp(j.data, U * j.in_eigenbasis * U', j.orthogonal, j.hermitian) for j in jumps]
  else
      JumpOp[j for j in jumps]
  end
  ```
  With:
  ```julia
  jumps_for_diss = if config.domain isa TrotterDomain && ham_or_trott isa TrottTrott
      JumpOp[JumpOp(j.data, ham_or_trott.eigvecs' * j.data * ham_or_trott.eigvecs, j.orthogonal, j.hermitian) for j in jumps]
  else
      JumpOp[j for j in jumps]
  end
  ```

  **7. `src/ofts.jl` -- comment on line 136:**
  Remove or update the comment that references `trotter.trafo_from_eigen_to_trotter`. Replace with a note like:
  ```julia
  # To compare trotter OFT with energy OFT: transform via trotter.eigvecs and hamiltonian.eigvecs
  ```
  </action>
  <verify>
  ```bash
  cd /Users/bence/code/QuantumFurnace.jl && grep -rn "trafo_from_eigen_to_trotter" src/
  ```
  Should return zero matches. The struct should have 5 fields (t0, num_trotter_steps_per_t0, eigvals_t0, eigvecs, bohr_freqs).
  </verify>
  <done>No source files reference trafo_from_eigen_to_trotter. TrottTrott struct has 5 fields.</done>
</task>

<task type="auto">
  <name>Task 2: Update test code and verify all tests pass</name>
  <files>
    test/test_dm_scaling.jl
  </files>
  <action>
  **1. `test/test_dm_scaling.jl` -- DMTST-05 (line 116):**
  This line transforms `B_trott` (in Trotter eigenbasis) back to H-eigenbasis for comparison:
  ```julia
  B_trott_in_eigen = TEST_TROTTER.trafo_from_eigen_to_trotter' * B_trott * TEST_TROTTER.trafo_from_eigen_to_trotter
  ```
  Replace with the equivalent using eigvecs directly. Note the math:
  `trafo' = (trotter.eigvecs' * ham.eigvecs)' = ham.eigvecs' * trotter.eigvecs`
  So `trafo' * X * trafo = (ham.eigvecs' * trotter.eigvecs) * X * (trotter.eigvecs' * ham.eigvecs)`.
  Replace with:
  ```julia
  U_t2e = TEST_TROTTER.eigvecs' * TEST_HAM.eigvecs  # Trotter eigenbasis <- H eigenbasis
  B_trott_in_eigen = U_t2e' * B_trott * U_t2e
  ```

  **2. `test/test_dm_scaling.jl` -- DMTST-06 (lines 162-163):**
  Replace:
  ```julia
  U = TEST_TROTTER.trafo_from_eigen_to_trotter  # H-eigen -> Trotter eigen
  jump_trott = JumpOp(jump.data, U * jump.in_eigenbasis * U', jump.orthogonal, jump.hermitian)
  ```
  With:
  ```julia
  jump_trott = JumpOp(jump.data, TEST_TROTTER.eigvecs' * jump.data * TEST_TROTTER.eigvecs, jump.orthogonal, jump.hermitian)
  ```
  Also update lines 167-168 where `U'` and `U` are used for transforming OFT result back:
  ```julia
  A_trott_in_eigen = U' * A_trott * U
  ```
  Replace with:
  ```julia
  U_t2e = TEST_TROTTER.eigvecs' * TEST_HAM.eigvecs
  A_trott_in_eigen = U_t2e' * A_trott * U_t2e
  ```

  **3. `test/test_dm_scaling.jl` -- DMTST-06b (lines 230-231):**
  Same pattern. Replace:
  ```julia
  U = TEST_TROTTER.trafo_from_eigen_to_trotter
  jump_trott = JumpOp(jump.data, U * jump.in_eigenbasis * U', jump.orthogonal, jump.hermitian)
  ```
  With:
  ```julia
  jump_trott = JumpOp(jump.data, TEST_TROTTER.eigvecs' * jump.data * TEST_TROTTER.eigvecs, jump.orthogonal, jump.hermitian)
  ```
  And update the `U'`/`U` transforms of OFT results (around line 236):
  ```julia
  A_nufft_trott_in_eigen = U' * A_nufft_trott * U
  ```
  and (around line 243):
  ```julia
  A_trott_in_eigen = U' * A_trott * U
  ```
  Replace both with:
  ```julia
  U_t2e = TEST_TROTTER.eigvecs' * TEST_HAM.eigvecs
  A_nufft_trott_in_eigen = U_t2e' * A_nufft_trott * U_t2e
  ...
  A_trott_in_eigen = U_t2e' * A_trott * U_t2e
  ```
  (Define `U_t2e` once near the top of the testset for DMTST-06b.)

  **4. Run the full test suite:**
  ```bash
  cd /Users/bence/code/QuantumFurnace.jl && julia --project -e 'using Pkg; Pkg.test()'
  ```
  All tests must pass. The key tests to watch:
  - DMTST-05 (coherent term B consistency) -- verifies B_trotter still matches B_bohr
  - DMTST-06 (OFT consistency) -- verifies Trotter OFT still matches energy OFT
  - DMTST-06b (NUFFT OFT consistency) -- same for NUFFT variant
  - DMTST-01/02 (Gibbs fixed point) -- verifies TrotterDomain Lindbladian still correct
  - All trajectory tests
  </action>
  <verify>
  ```bash
  cd /Users/bence/code/QuantumFurnace.jl && grep -rn "trafo_from_eigen_to_trotter" test/ src/
  ```
  Should return zero matches across the entire codebase (excluding .planning/ docs).

  ```bash
  cd /Users/bence/code/QuantumFurnace.jl && julia --project -e 'using Pkg; Pkg.test()'
  ```
  All tests pass.
  </verify>
  <done>Zero references to trafo_from_eigen_to_trotter in source or test code. All tests pass with identical numerical results.</done>
</task>

</tasks>

<verification>
1. `grep -rn "trafo_from_eigen_to_trotter" src/ test/` returns no matches.
2. TrottTrott struct has exactly 5 fields: t0, num_trotter_steps_per_t0, eigvals_t0, eigvecs, bohr_freqs.
3. `julia --project -e 'using Pkg; Pkg.test()'` -- all tests pass.
4. Key numerical tests (DMTST-01, DMTST-02, DMTST-05, DMTST-06, DMTST-06b) show identical thresholds.
</verification>

<success_criteria>
- TrottTrott struct no longer has trafo_from_eigen_to_trotter field
- All basis transformations in source code use trotter.eigvecs directly (via `trotter.eigvecs' * jump.data * trotter.eigvecs`)
- All basis transformations in test code use the explicit eigvec products
- Full test suite passes with no regressions
</success_criteria>

<output>
After completion, create `.planning/quick/9-remove-trafo-from-eigen-to-trotter-from-/9-SUMMARY.md`
</output>
