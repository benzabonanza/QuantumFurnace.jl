# Architecture Patterns: QuantumFurnace.jl Restructure

**Domain:** Julia scientific computing codebase restructure (8,312 LOC src, 5,071 LOC test)
**Researched:** 2026-02-25
**Overall confidence:** HIGH (based on complete codebase read, every source file examined)

## Current Architecture

### Module Structure

Single module (`QuantumFurnace`), flat `src/` directory, 28 source files loaded via sequential `include()` in `QuantumFurnace.jl`. The include order is load-order-dependent (e.g., `structs.jl` must precede `trajectories.jl` because `TrajectoryWorkspace` references `ConvergenceData`).

### Logical Layers (Current)

```
Layer 0: Foundation
  constants.jl (4), hamiltonian.jl (331), trotter_domain.jl (206)
  structs.jl (358), qi_tools.jl (220), misc_tools.jl (326)

Layer 1: Domain Physics
  time_domain.jl (19), nufft.jl (96), ofts.jl (110)
  errors.jl (1), kraus.jl (14)
  energy_domain.jl (165), bohr_domain.jl (184)
  coherent.jl (394)

Layer 2: Simulation Engines
  jump_workers.jl (461) -- dense Liouvillian construction + DM thermalization inner loop
  trajectories.jl (1139) -- trajectory engine (biggest file)
  furnace_utensils.jl (135) -- _precompute_data() dispatchers
  furnace.jl (163) -- run_lindbladian(), run_thermalization(), construct_lindbladian()

Layer 3: Krylov
  krylov_workspace.jl (498) -- KrylovWorkspace struct + constructors
  krylov_matvec.jl (586) -- apply_lindbladian!(), apply_adjoint_lindbladian!()
  krylov_eigsolve.jl (568) -- krylov_spectral_gap(), apply_delta_channel!()

Layer 4: Analysis
  log_sobolev.jl (206), convergence.jl (395), fitting.jl (217)
  gap_estimation.jl (361), diagnostics.jl (573)

Layer 5: Persistence
  results.jl (453) -- ExperimentResult, BSON serialization
```

### Four Simulation Paths

| Path | Entry Point | Config Types | Workspace | Inner Loop |
|------|------------|--------------|-----------|------------|
| Dense Lindbladian | `run_lindbladian()` | AbstractLiouvConfig | LindbladianWorkspace | `_jump_contribution!` (vectorized) |
| DM Thermalization | `run_thermalization()` | AbstractThermalizeConfig | KrausScratch | `_jump_contribution!` (Kraus) |
| Krylov Spectrum | `krylov_spectral_gap()` | AbstractLiouvConfig or AbstractThermalizeConfig | KrylovWorkspace | `apply_lindbladian!` / `apply_delta_channel!` |
| Trajectories | `run_trajectories()` | AbstractThermalizeConfig | TrajectoryWorkspace + TrajectoryFramework | `step_along_trajectory!` |

### Config Hierarchy (Current -- 4 concrete types, massive field duplication)

```
AbstractConfig{D,T}
  +-- AbstractLiouvConfig{D,T}
  |     +-- LiouvConfig{D,T}       (14 fields)
  |     +-- LiouvConfigGNS{D,T}    (14 identical fields, with_coherent forced false)
  +-- AbstractThermalizeConfig{D,T}
        +-- ThermalizeConfig{D,T}      (16 fields = LiouvConfig + mixing_time + delta)
        +-- ThermalizeConfigGNS{D,T}   (16 identical fields, with_coherent forced false)
```

**Problem:** 4 structs with 14-16 fields each, ~60 lines of field definitions are copy-pasted across them. Every field change requires 4 simultaneous edits. The KMS/GNS distinction is only `with_coherent=false` plus different alpha function selection (dispatched at runtime via `_pick_alpha`).

### Workspace Types (Current -- 4 types with overlapping fields)

| Type | Fields | Used By |
|------|--------|---------|
| `LindbladianWorkspace{T}` | Id, jump_tmp, jump_conj, jump_dag_jump, jump2_jump1 | `construct_lindbladian` |
| `KrausScratch{T}` | jump_oft, LdagL, R, rho_jump, K0, tmp1, tmp2, rho_next | `run_thermalization` |
| `KrylovWorkspace{T,PD}` | jump_oft, tmp1, tmp2, LdagL, rho_out + channel fields + G matrices | `apply_lindbladian!`, `apply_delta_channel!` |
| `TrajectoryWorkspace{T}` | jump_oft, psi_tmp, Rpsi, rho_acc | `step_along_trajectory!` |

**Overlap:** `jump_oft`, `tmp1`/`tmp2`, `LdagL` appear in 3 of 4 workspace types.

### Prefactor Duplication (The Worst Offender)

The domain-dependent scalar prefactor is computed identically in **10+ locations**:

1. `jump_workers.jl:65` -- dense Liouvillian EnergyDomain
2. `jump_workers.jl:109` -- dense Liouvillian Time/TrotterDomain
3. `jump_workers.jl:334` -- DM thermalization EnergyDomain
4. `jump_workers.jl:410` -- DM thermalization Time/TrotterDomain
5. `krylov_matvec.jl:170` -- Krylov EnergyDomain
6. `krylov_matvec.jl:483` -- Krylov Time/TrotterDomain
7. `krylov_workspace.jl:258` -- R_total accumulation EnergyDomain
8. `krylov_workspace.jl:303` -- R_total accumulation Time/TrotterDomain
9. `trajectories.jl:186-189` -- TrajectoryFramework builder
10. `krylov_eigsolve.jl:231,277` -- Channel sandwich accumulation

There are two formulas:
- **EnergyDomain:** `config.w0 / (config.sigma * sqrt(2 * pi)) * gamma_norm_factor`
- **Time/TrotterDomain:** `config.w0 * config.t0^2 * (config.sigma * sqrt(2 / pi)) / (2 * pi) * gamma_norm_factor`

These should be a single `_dissipator_prefactor(config, gamma_norm_factor)` function dispatched on domain type.

---

## Recommended Restructure Architecture

### Target Directory Layout

```
src/
  QuantumFurnace.jl           # Module entry point, exports, include()s
  types/
    domains.jl                # AbstractDomain, BohrDomain, EnergyDomain, TimeDomain, TrotterDomain
    configs.jl                # AbstractConfig hierarchy (unified)
    hamiltonian.jl            # HamHam, TrottTrott (merged from hamiltonian.jl + trotter_domain.jl)
    jump_operator.jl          # JumpOp struct
    results.jl                # All result structs consolidated
    workspaces.jl             # All workspace structs
  physics/
    transition_functions.jl   # pick_transition, create_alpha, create_f (from energy/bohr/misc)
    coherent.jl               # B operators (B_time, B_trotter, B_bohr)
    oft.jl                    # oft!(), NUFFT prefactors, _prefactor_view
    precompute.jl             # _precompute_data, _precompute_labels (from furnace_utensils.jl)
    dissipator.jl             # Shared dissipator logic: prefactors, half-grid, sandwiches
  simulation/
    lindbladian.jl            # construct_lindbladian, run_lindbladian (vectorized path)
    thermalization.jl         # run_thermalization (DM Kraus path)
    krylov.jl                 # KrylovWorkspace construction, apply_lindbladian!, krylov_spectral_gap
    trajectories.jl           # TrajectoryFramework, step_along_trajectory!, run_trajectories
  analysis/
    convergence.jl
    fitting.jl
    gap_estimation.jl
    diagnostics.jl
  util/
    qi_tools.jl
    pauli.jl                  # X, Y, Z, pad_term, pauli_string_to_matrix (from misc_tools.jl)
    validation.jl             # validate_config!, error helpers
    persistence.jl            # ExperimentResult serialization (save/load)
    constants.jl
```

### Why This Layout

1. **types/ first:** All struct definitions in one place. Include order becomes trivial -- types before anything else.
2. **physics/ is the shared kernel:** `dissipator.jl` extracts the duplicated prefactor/half-grid/sandwich logic into shared functions that all 4 simulation paths call.
3. **simulation/ per path:** Each simulation engine gets its own file, but they share the physics/ layer.
4. **analysis/ is leaf-level:** These files only depend on simulation outputs, not on simulation internals.

---

## Component Boundaries and Dependencies

### Dependency Graph (What Depends on What)

```
types/domains.jl         -> (nothing)
types/configs.jl         -> types/domains.jl
types/hamiltonian.jl     -> types/domains.jl
types/jump_operator.jl   -> (nothing)
types/results.jl         -> types/configs.jl
types/workspaces.jl      -> types/jump_operator.jl

physics/transition_functions.jl  -> types/configs.jl
physics/oft.jl                   -> types/jump_operator.jl, types/hamiltonian.jl
physics/coherent.jl              -> types/jump_operator.jl, types/hamiltonian.jl, physics/oft.jl
physics/precompute.jl            -> types/configs.jl, physics/transition_functions.jl, physics/oft.jl
physics/dissipator.jl            -> types/configs.jl, physics/precompute.jl, physics/oft.jl

simulation/lindbladian.jl    -> physics/dissipator.jl, types/workspaces.jl
simulation/thermalization.jl -> physics/dissipator.jl, types/workspaces.jl
simulation/krylov.jl         -> physics/dissipator.jl, types/workspaces.jl
simulation/trajectories.jl   -> physics/dissipator.jl, types/workspaces.jl

analysis/*                   -> simulation/*, types/results.jl
util/persistence.jl          -> types/results.jl, types/configs.jl
```

### Critical Coupling Points (Verified by Code Read)

| Coupling | Specific Files/Functions Affected | Severity |
|----------|----------------------------------|----------|
| Config field access | Every `config.sigma`, `config.w0`, `config.t0`, `config.beta` -- appears in 20+ functions across 8 files | HIGH |
| `_precompute_data()` NamedTuple shape | furnace_utensils.jl returns domain-specific NamedTuples; consumed by jump_workers.jl, trajectories.jl, krylov_workspace.jl, coherent.jl, krylov_matvec.jl, krylov_eigsolve.jl | HIGH |
| `JumpOp.in_eigenbasis` convention | Every simulation path accesses this field; jump_workers.jl, trajectories.jl, krylov_matvec.jl all iterate over it | MEDIUM |
| Prefactor formulas | 10+ locations compute identical `base_prefactor` from config fields | HIGH (duplication) |
| `_thermalize_to_liouv_config()` | krylov_workspace.jl:190-232 manually copies all 14 config fields | HIGH (breaks if config fields change) |
| `_config_to_dict()` / `_reconstruct_config()` | results.jl serialization hard-codes all config field names | HIGH (breaks if config fields change) |
| Test factory functions | test_helpers.jl:228-410 has 8 factory functions constructing configs with all 14-16 fields | MEDIUM |
| `include()` order | QuantumFurnace.jl:102-128 must be in exact dependency order | LOW (just bookkeeping) |

---

## Restructure Sequencing: Build Order That Minimizes Breakage

### Guiding Principles

1. **Leaves first, roots last.** Change what nothing depends on first. Change what everything depends on last.
2. **One axis of change per phase.** Never change struct layout and function signatures simultaneously.
3. **Tests green after every phase.** Each phase has a defined test-passing checkpoint.
4. **Backward-compatible intermediate states.** Use re-exports, type aliases, and forwarding methods during transition.

### Phase Dependency Graph

```
Phase 1: File Reorganization (include-path-only)
   |
   v
Phase 2: Function Deduplication (extract shared helpers)
   |
   v
Phase 3: Workspace Consolidation (unify naming, reduce redundancy)
   |
   v
Phase 4: Config Hierarchy Redesign (reduce 4 config types to 2)
   |
   v
Phase 5: Result Struct Cleanup (depends on new Config types)
   |
   v
Phase 6: Test Cleanup (depends on all above being stable)
```

**Parallelization:** Phases 2 and 3 are partially parallelizable (they touch different files), but both must complete before Phase 4. Phases 5 and 6 are sequential.

---

### Phase 1: File Reorganization (Include-Path-Only Changes)

**Goal:** Move files into subdirectories, update `include()` calls. Zero logic changes.

**Risk:** LOW. Julia's module system does not care about file paths -- `include()` just evaluates the file content. As long as include order is preserved, this is a pure rename operation.

**Files affected:**
- `src/QuantumFurnace.jl` -- rewrite all 27 `include()` paths
- All 28 source files -- rename/move

**Specific steps:**
1. Create `src/types/`, `src/physics/`, `src/simulation/`, `src/analysis/`, `src/util/`
2. Move files (use `git mv` to preserve history):
   - `structs.jl` -> split into `types/domains.jl`, `types/configs.jl`, `types/jump_operator.jl`, `types/results.jl`, `types/workspaces.jl`
   - `hamiltonian.jl` + `trotter_domain.jl` -> `types/hamiltonian.jl`
   - `energy_domain.jl` + `bohr_domain.jl` + part of `misc_tools.jl` -> `physics/transition_functions.jl`
   - `coherent.jl` -> `physics/coherent.jl`
   - `ofts.jl` + `nufft.jl` -> `physics/oft.jl`
   - `furnace_utensils.jl` -> `physics/precompute.jl`
   - `jump_workers.jl` -> split: vectorized part to `simulation/lindbladian.jl`, Kraus part to `simulation/thermalization.jl`
   - `furnace.jl` -> split between `simulation/lindbladian.jl` and `simulation/thermalization.jl`
   - `krylov_workspace.jl` + `krylov_matvec.jl` + `krylov_eigsolve.jl` -> `simulation/krylov.jl`
   - `trajectories.jl` -> `simulation/trajectories.jl`
   - etc.
3. Update `include()` in module entry point
4. Run full test suite

**Test strategy:** Run `] test` after each file move. The test suite does NOT reference source file paths -- it only uses `using QuantumFurnace`. So renaming source files cannot break tests as long as `include()` order is correct.

**Checkpoint:** All 600+ tests pass with new file layout.

---

### Phase 2: Function Deduplication (Extract Shared Helpers)

**Goal:** Extract duplicated code into shared functions without changing any public API or struct definitions.

**Risk:** MEDIUM. Changing internal function signatures could break zero-allocation guarantees if function boundaries cause unexpected allocations. Must verify with allocation tests.

**Specific extractions:**

#### 2A: Dissipator Prefactor Helper
**Current:** 10+ locations compute `base_prefactor` from config fields.
**Target:** Single function `_dissipator_prefactor(config, gamma_norm_factor)` dispatched on domain.

```julia
function _dissipator_prefactor(config::AbstractConfig{EnergyDomain}, gamma_norm_factor)
    return config.w0 / (config.sigma * sqrt(2 * pi)) * gamma_norm_factor
end

function _dissipator_prefactor(config::AbstractConfig{D}, gamma_norm_factor) where {D<:Union{TimeDomain,TrotterDomain}}
    return config.w0 * config.t0^2 * (config.sigma * sqrt(2 / pi)) / (2 * pi) * gamma_norm_factor
end
```

**Files touched:** jump_workers.jl, krylov_matvec.jl, krylov_workspace.jl, krylov_eigsolve.jl, trajectories.jl
**LOC saved:** ~30 duplicated lines removed, ~8 lines added

#### 2B: Half-Grid Iteration Pattern
**Current:** The "iterate w_raw, skip positive, abs(w_raw), mirror for hermitian" pattern is duplicated ~12 times across all simulation paths.
**Target:** Consider a callback-based helper or macro.

**CRITICAL CAUTION:** This MUST be benchmarked before adoption. If the callback causes allocations or prevents inlining, keep the pattern duplicated and just document it. Performance trumps DRY in hot paths. Closures capturing mutable state may allocate.

**Recommendation:** Start by just extracting the prefactor (2A) which is pure computation, no closure risk. Defer half-grid extraction until profiling confirms it is allocation-free.

#### 2C: OFT Computation Unification
**Current:** `oft!()` in ofts.jl:1-5 and `_krylov_oft!()` in krylov_matvec.jl:11-19 do the same computation.
**Target:** Single `_oft_energy!()` function used by all paths.

```julia
@inline function _oft_energy!(out, eigenbasis, bohr_freqs, energy, inv_4sigma2)
    @. out = eigenbasis * exp(-(energy - bohr_freqs)^2 * inv_4sigma2)
end
```

**Files touched:** ofts.jl, krylov_matvec.jl, jump_workers.jl (EnergyDomain paths that call `oft!`)

#### 2D: Identify Krylov/DM Shared Logic (Assessment Only)
**Current:** `apply_delta_channel!()` in krylov_eigsolve.jl and `_jump_contribution!` for AbstractThermalizeConfig in jump_workers.jl both implement the Chen CPTP channel. The user's design notes explicitly call out this overlap.
**Target:** Identify exactly which functions can be shared. The key difference: DM thermalization computes R^a, K0^a, U_residual^a **per jump per step** while Krylov uses summed R_total, K0_total, U_residual_total **precomputed once**.

**Assessment:** The inner sandwich accumulation (`rate * L * rho * L'`) is sharable. The outer channel structure (per-operator vs summed) is NOT. Extract sandwich helpers; keep channel structure separate.

**Test strategy:** After each extraction:
1. Run `test_allocation.jl` to verify zero-allocation invariants
2. Run `test_regression.jl` for numerical correctness
3. Run `test_krylov_crossvalidation.jl` for Krylov vs dense agreement

**Checkpoint:** All tests pass, `test_allocation.jl` allocations unchanged.

---

### Phase 3: Workspace Consolidation

**Goal:** Unify workspace field naming. Reduce redundancy where safe.

**Risk:** MEDIUM. Workspaces are used in hot paths. Changing field names requires updating every reference.

**Strategy:** Do NOT merge all 4 workspace types into one. They serve genuinely different purposes:
- `LindbladianWorkspace` handles `dim^2 x dim^2` vectorized representation -- fundamentally different scale
- `KrausScratch` and `KrylovWorkspace` share `jump_oft, LdagL, tmp1, tmp2` but KrylovWorkspace has additional channel and G-matrix fields
- `TrajectoryWorkspace` uses vectors (psi_tmp, Rpsi) not matrices -- fundamentally different data shapes

**Concrete plan:**

1. **Keep 4 workspace types** but unify field names:
   - All matrix scratch: `jump_oft`, `tmp1`, `tmp2`, `LdagL`
   - Rename `LindbladianWorkspace.jump_tmp` -> `jump_oft`
   - Rename `LindbladianWorkspace.jump_conj` -> `tmp1`
   - Rename `LindbladianWorkspace.jump_dag_jump` -> `LdagL`
   - Rename `LindbladianWorkspace.jump2_jump1` -> `tmp2`

2. **Consider removing LindbladianWorkspace.Id:** The identity matrix is only used once in `_vectorize_liouvillian_coherent!`. Can be replaced with `Matrix{CT}(I, dim, dim)` at the call site or passed as an argument.

3. **Document workspace aliasing:** Add comments documenting which workspace fields are aliased during multi-step operations (e.g., `rho_eff = ws.LdagL` in `apply_delta_channel!`).

**Files touched:** structs.jl (workspace definitions), jump_workers.jl (vectorized path field references), all functions that access LindbladianWorkspace fields.

**Test strategy:** Run `test_allocation.jl` after each rename to catch regressions.

**Checkpoint:** All tests pass. Workspace field names are consistent across types.

---

### Phase 4: Config Hierarchy Redesign

**Goal:** Eliminate the 4-config explosion. Reduce to 2 concrete types with KMS/GNS as a runtime field.

**Risk:** HIGH. This is the most invasive change. Config types are referenced in:
- Every simulation entry point (furnace.jl, trajectories.jl)
- Every `_precompute_data` dispatcher (furnace_utensils.jl)
- Every `_jump_contribution!` method (jump_workers.jl)
- KrylovWorkspace constructors (krylov_workspace.jl)
- `_thermalize_to_liouv_config` (krylov_workspace.jl)
- All serialization code (results.jl)
- 8 test factory functions (test_helpers.jl)
- 14 direct config constructor calls in test files

**Strategy:** Phased internal migration with backward compatibility at each step.

**Step 4A: Add db_type flag to existing configs**
- Add `db_type::Symbol = :kms` to LiouvConfig and ThermalizeConfig
- Add `db_type::Symbol = :gns` to LiouvConfigGNS and ThermalizeConfigGNS
- Replace `config isa Union{LiouvConfigGNS, ThermalizeConfigGNS}` checks with `config.db_type == :gns` (7 occurrences in results.jl, diagnostics, companion text)
- **Files touched:** structs.jl only for field addition; results.jl for isa-check replacement
- Run tests: must pass unchanged (new field has a default)

**Step 4B: Migrate dispatch to use db_type**
- `_pick_alpha(config)` already dispatches on config type for KMS vs GNS. Change to dispatch on `config.db_type`.
- `_select_b_plus_calculator(config)` only handles KMS types. Guard with `config.db_type == :kms`.
- **Files touched:** energy_domain.jl (or wherever `_pick_alpha` lives), furnace_utensils.jl, coherent.jl
- Run tests: must pass

**Step 4C: Deprecate GNS config types**
- Add constructors that create `LiouvConfig(; kwargs..., db_type=:gns, with_coherent=false)` when `LiouvConfigGNS(; kwargs...)` is called
- Same for `ThermalizeConfigGNS`
- **Files touched:** structs.jl (add forwarding constructors)
- Run tests: must pass via forwarding

**Step 4D: Remove GNS config types**
- Delete `LiouvConfigGNS` and `ThermalizeConfigGNS` struct definitions
- Delete their outer constructors
- Simplify `_thermalize_to_liouv_config` to one method (no GNS variant needed)
- Update test factories to use `db_type=:gns` keyword
- Update all 14 test constructor calls
- **Files touched:** structs.jl, krylov_workspace.jl, test_helpers.jl, test_gns_trajectory.jl, test_krylov_crossvalidation.jl, test_dm_detailed_balance.jl

**Step 4E: Simplify _thermalize_to_liouv_config**
After removing GNS types, this function becomes:
```julia
function _thermalize_to_liouv_config(tc::ThermalizeConfig)
    LiouvConfig(; (f => getfield(tc, f) for f in fieldnames(LiouvConfig))...)
end
```
This is now a single method that works for both KMS and GNS.

**Recommendation on Liouv/Thermalize unification:** Keep 2 concrete types (`LiouvConfig` and `ThermalizeConfig`). The dispatch separation (AbstractLiouvConfig vs AbstractThermalizeConfig) is load-bearing in furnace.jl, trajectories.jl, krylov_eigsolve.jl, and krylov_workspace.jl. Merging them would require runtime checks everywhere instead of compile-time dispatch.

**Test strategy:** Run full suite after each sub-step (4A through 4E). The sub-steps are designed so tests pass after each one.

**Checkpoint:** 2 config types (LiouvConfig, ThermalizeConfig), `db_type` field for KMS/GNS, all tests pass.

---

### Phase 5: Result Struct Cleanup

**Goal:** Consolidate result types in one location, simplify serialization.

**Risk:** LOW-MEDIUM. Results are leaf-level (only constructed, not deeply coupled).

**Specific changes:**

1. **Consolidate all result struct definitions into `types/results.jl`:**
   - `LindbladianResult` (currently structs.jl)
   - `DMSimulationResult` (currently structs.jl)
   - `ConvergenceData` (currently structs.jl)
   - `TrajectoryResult` (currently trajectories.jl:23-30)
   - `ObservableTrajectoryResult` (currently trajectories.jl:37-66)
   - `PerOperatorKraus` (currently trajectories.jl:72-77)
   - `KrylovGapResult` (currently krylov_eigsolve.jl:40-51)
   - `SpectralGapResult` (currently gap_estimation.jl:37-50)
   - `FitResult` (currently fitting.jl)
   - Diagnostics result structs (currently diagnostics.jl:33-90)
   - `ExperimentResult` (currently results.jl:18-23)

2. **Simplify BSON serialization:** After Config redesign (Phase 4), `_config_to_dict` and `_reconstruct_config` simplify because there are only 2 config types + `db_type` field. The `config_type` tag becomes just `d[:db_type] = config.db_type`.

3. **Clean up ConvergenceData backward compat:** The 6-arg outer constructor exists for BSON backward compatibility from Phase 16->17 transition. Decide if old data format support is still needed.

**Files touched:** results.jl (major simplification), structs.jl (struct definitions move out), trajectories.jl (struct definitions move out), krylov_eigsolve.jl (struct definition moves out), gap_estimation.jl (struct definition moves out), diagnostics.jl (struct definitions move out), fitting.jl (struct definition moves out), test_results.jl

**Test strategy:** Run `test_results.jl` (362 LOC) after changes. Test BSON round-trip serialization.

**Checkpoint:** All result types in one location, serialization simplified, all tests pass.

---

### Phase 6: Test Cleanup

**Goal:** Reduce test helper duplication, consolidate factory functions, slim test infrastructure.

**Risk:** LOW. Tests don't affect production code.

**Specific changes:**

1. **Consolidate factory functions:** After Phase 4 removes GNS types, reduce from 8 to 4 factories:
   - `make_liouv_config(domain; db_type=:kms, with_coherent=true)`
   - `make_thermalize_config(domain; db_type=:kms, with_coherent=true, delta=..., mixing_time=...)`
   - `make_small_liouv_config(domain; db_type=:kms, with_coherent=false)`
   - `make_small_thermalize_config(domain; db_type=:kms, with_coherent=false, delta=..., mixing_time=...)`

2. **Parameterize test systems:** Currently `make_test_system()` and `make_small_test_system()` are nearly identical (only num_qubits differs). Replace with:
   ```julia
   make_test_system(; num_qubits=4, trotter=nothing)
   ```

3. **Remove old_tests/**: 7 files in `test/old_tests/` are unused by `runtests.jl`. Archive or delete.

4. **Add test info output:** Add `@info` at the start of each `@testset` block for visibility during test runs.

5. **Review thresholds:** The user flagged that some testing thresholds have been "off" and "more deceiving than useful." Systematically review:
   - `test_convergence.jl` (763 LOC, largest test file) -- check convergence thresholds
   - `test_krylov_crossvalidation.jl` (503 LOC) -- check tolerance values
   - `test_dm_scaling.jl` (268 LOC) -- check scaling tolerances

**Files touched:** test_helpers.jl, all test files (factory function signature changes)

**Checkpoint:** Test LOC reduced, factory functions consolidated, old_tests archived, all tests pass.

---

## Integration Points Between Phases

### Phase 1 -> Phase 2

Phase 1 moves files to subdirectories. Phase 2 adds new shared functions (e.g., `_dissipator_prefactor`) that need to go in `physics/dissipator.jl`. The file layout from Phase 1 determines WHERE these new functions live.

**Integration action:** Phase 2 creates `physics/dissipator.jl` and adds it to the `include()` list from Phase 1.

### Phase 2 -> Phase 3

Phase 2 extracts shared functions. Phase 3 renames workspace fields. Some extracted functions will reference workspace field names (e.g., `ws.jump_oft`).

**Integration action:** Complete Phase 2 first using current field names, then Phase 3 renames uniformly. This avoids coordinating two simultaneous changes.

### Phase 3 -> Phase 4

Phase 3 touches workspace constructors. Phase 4 changes config types. `KrylovWorkspace` constructors dispatch on config type (`AbstractLiouvConfig` / `AbstractThermalizeConfig`).

**Integration action:** Phase 4 updates `KrylovWorkspace` constructor signatures. This is straightforward because the constructors already use abstract types, not concrete ones.

### Phase 4 -> Phase 5

Phase 4 changes config struct definitions. Phase 5 updates serialization code that hard-codes config field names.

**Integration action:** Phase 5 MUST follow Phase 4. The `_config_to_dict` / `_reconstruct_config` functions reference specific config types and field names.

### Phase 4 -> Phase 6

Phase 4 changes config constructors. Phase 6 updates test factory functions that construct configs.

**Integration action:** Phase 6 MUST follow Phase 4. Test factories construct configs by name.

---

## Strategy for Keeping Tests Green Throughout

### The Test Safety Net

| Test File | What It Guards | Earliest Phase Affected |
|-----------|---------------|------------------------|
| test_compilation.jl (56 LOC) | Module loads, exports exist | Phase 1 (include paths) |
| test_aqua.jl (9 LOC) | No ambiguities, no unbound args | Phase 4 (config changes) |
| test_allocation.jl (171 LOC) | Zero-allocation hot paths | Phase 2, 3 |
| test_regression.jl (157 LOC) | Numerical correctness vs reference | Phase 2 |
| test_cptp.jl (67 LOC) | CPTP channel preservation | Phase 2, 3 |
| test_dm_detailed_balance.jl (83 LOC) | DM thermalization physics | Phase 2, 3, 4 |
| test_dm_scaling.jl (268 LOC) | Multi-domain DM correctness | Phase 2, 3, 4 |
| test_krylov_matvec.jl (401 LOC) | Krylov matvec correctness | Phase 2, 3 |
| test_krylov_eigsolve.jl (210 LOC) | Krylov eigsolve | Phase 2, 3, 4 |
| test_krylov_crossvalidation.jl (503 LOC) | Krylov vs dense agreement | Phase 2, 3 |
| test_trajectory_fixes.jl (128 LOC) | Trajectory CPTP and normalization | Phase 2, 3 |
| test_gns_trajectory.jl (125 LOC) | GNS trajectory path | Phase 4 (GNS type removal) |
| test_results.jl (362 LOC) | BSON serialization round-trip | Phase 4, 5 |
| test_convergence.jl (763 LOC) | Convergence tracking | Phase 2 |
| test_fitting.jl (155 LOC) | Exponential fitting | Phase 5 |
| test_observable_trajectories.jl (85 LOC) | Observable trajectory path | Phase 2 |
| test_gap_estimation.jl (279 LOC) | Gap estimation pipeline | Phase 2 |
| test_diagnostics.jl (524 LOC) | Exact diagnostics | Phase 2 |
| test_threading.jl (195 LOC) | Multi-threaded correctness | Phase 2, 3 |
| test_workspace_independence.jl (91 LOC) | Workspace isolation | Phase 3 |

### Testing Protocol Per Phase

**Before each phase:**
```bash
] test   # Baseline: all pass
```

**Phase-specific fast feedback loops:**

| Phase | Fast Test Set (run after each file change) |
|-------|-------------------------------------------|
| 1 | test_compilation.jl (verifies include paths and exports) |
| 2 | test_allocation.jl + test_regression.jl + test_krylov_crossvalidation.jl |
| 3 | test_allocation.jl + test_workspace_independence.jl |
| 4 | test_aqua.jl + test_gns_trajectory.jl + test_results.jl + test_dm_scaling.jl |
| 5 | test_results.jl |
| 6 | Full suite (test cleanup is the final validation) |

### Rollback Strategy

Each phase should be a single Git branch with atomic commits per sub-step. If a phase breaks tests:

1. `git stash` or `git checkout` to last green commit
2. Investigate which specific change broke the test
3. Fix forward (preferred) or revert the specific commit

**Critical:** Never proceed to the next phase with failing tests. Each phase boundary IS a green test suite.

---

## Patterns to Follow

### Pattern 1: Domain Dispatch via Type Parameters (Already Used, Keep)
**What:** Julia's multiple dispatch handles domain-specific logic without if-else chains.
**Example:** `_precompute_data(config::AbstractConfig{EnergyDomain}, ...)` vs `_precompute_data(config::AbstractConfig{TimeDomain}, ...)`
**Preserve this pattern.** Do not consolidate into runtime `if domain isa ...` branches. The parametric dispatch on `{D<:AbstractDomain}` gives Julia's compiler type-stable, specialized code paths.

### Pattern 2: Precompute Once, Iterate Many (Already Used, Keep)
**What:** `_precompute_data()` runs once at setup; its NamedTuple is passed to every hot-path call.
**Preserve this pattern.** The precomputed_data NamedTuple is the right abstraction for domain-varying cached data.

### Pattern 3: Workspace Pre-allocation (Already Used, Strengthen)
**What:** All scratch matrices allocated once, reused via in-place operations.
**Strengthen by:** Adding comments documenting which workspace fields are aliased during multi-step operations (e.g., `rho_eff = ws.LdagL` in `apply_delta_channel!`).

### Pattern 4: Config Validation at Entry Point Only
**What:** `validate_config!` is called once in `run_lindbladian`, `run_thermalization`, etc.
**Preserve:** Do not add validation inside hot paths.

## Anti-Patterns to Avoid

### Anti-Pattern 1: Deep Abstraction for Hot Paths
**What:** Creating abstract helper functions that wrap 2-3 BLAS calls behind a function boundary.
**Why bad:** Each function boundary can prevent inlining and cause allocations (especially closures that capture mutable state).
**Instead:** Inline short sequences. Use `@inline` annotation. Verify with `@allocated`. The half-grid iteration pattern is a good candidate for leaving duplicated if extraction causes allocations.

### Anti-Pattern 2: Merging Fundamentally Different Workspace Types
**What:** Creating a single "UnifiedWorkspace" with all possible fields from all 4 simulation paths.
**Why bad:** Adds unused fields to every hot-path struct; cache pollution; unclear field ownership; makes it impossible to reason about which fields are in use at any given point.
**Instead:** Keep separate types, unify naming conventions.

### Anti-Pattern 3: Runtime Dispatch Where Compile-Time Dispatch Exists
**What:** Replacing `f(config::AbstractLiouvConfig{EnergyDomain})` with `if config.domain isa EnergyDomain`.
**Why bad:** Loses Julia's specialization; hot paths become dynamically dispatched; type instability propagates.
**Instead:** Keep parametric dispatch. The domain type parameter `D` is explicitly designed for this.

### Anti-Pattern 4: Changing Public API During Restructure
**What:** Renaming `run_lindbladian` to `run_liouvillian` or changing function signatures while also moving files.
**Why bad:** Compounds changes; makes it impossible to verify that restructure is behavior-preserving.
**Instead:** Pure structural changes first (Phases 1-3). API changes (if desired) as a separate, later milestone.

---

## Scalability Considerations

| Concern | Current (n<=6 qubits) | At n=8 (dim=256) | At n=10+ |
|---------|----------------------|-------------------|----------|
| Struct field access | Negligible | Negligible | Negligible |
| Config construction | ~microseconds | Same | Same |
| Workspace allocation | ~1ms (dim=64) | ~100ms (dim^2 matrices) | Krylov-only path |
| File reorganization | No runtime impact | No runtime impact | No runtime impact |
| Test runtime | ~2-5 min full suite | Not affected by restructure | Not affected |

The restructure has zero impact on computational scalability. All changes are organizational/structural, not algorithmic.

---

## Sources

- All findings derived from direct codebase read of 28 source files and 22 test files in `/Users/bence/code/QuantumFurnace.jl/`
- `supplementary-informations/quantumfurnace-structure.md` -- author's design notes and future plans for DLL config, Hamiltonian accumulator, qiskit circuit estimation
- Julia documentation on module `include()` semantics (training data, HIGH confidence)
- Julia performance tips on struct field access, type stability, and `@inline` (training data, HIGH confidence)
