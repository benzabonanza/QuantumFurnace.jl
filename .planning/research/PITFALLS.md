# Domain Pitfalls: Major Codebase Restructure of QuantumFurnace.jl

**Domain:** Refactoring a Julia numerical simulation package (8,312 LOC src + 5,071 LOC test, 600+ tests)
**Researched:** 2026-02-25
**Confidence:** HIGH (grounded in direct codebase analysis + Julia documentation + known Julia compiler issues)

**Relationship to prior research:** This document supersedes the v1.5 PITFALLS.md (2026-02-20) which covered Krylov eigensolving pitfalls. This document covers pitfalls specific to the major codebase restructure: type hierarchy redesign, workspace consolidation, code deduplication, test restructuring, and file reorganization.

---

## Critical Pitfalls

Mistakes that cause silent numerical errors, allocation regressions, or test suite invalidation.

---

### CRIT-01: Abstract Type Field Boxing During Workspace Consolidation

**What goes wrong:**
When consolidating `KrausScratch`, `KrylovWorkspace`, `TrajectoryWorkspace`, and `LindbladianWorkspace` into fewer structs, it is tempting to use abstract types for fields that vary across use cases. For example, replacing the concrete `jump_eigenbases::Vector{Matrix{T}}` in `KrylovWorkspace` with a field typed `::Vector{Matrix{<:Complex}}` (the same abstract type that `JumpOp.in_eigenbasis` uses). In Julia, a struct field with an abstract type parameter causes the compiler to heap-allocate (box) every access to that field. The existing codebase already discovered this -- see `krylov_workspace.jl` lines 49-51:

```julia
# Concrete-typed jump data for zero-allocation hot path
# (JumpOp.in_eigenbasis is Matrix{<:Complex} -- abstract element type causes boxing)
jump_eigenbases::Vector{Matrix{T}}
```

Consolidating workspaces risks reintroducing exactly this pattern if the unified struct tries to hold matrices of different element types via abstract field types.

**How to detect:**
- `@allocated` tests in `test_allocation.jl` and `test_krylov_matvec.jl` will catch this: `@test allocs == 0` for `apply_lindbladian!` and `step_along_trajectory!`.
- Run `@code_warntype` on hot-path functions after consolidation. Any field access showing `Any` or `Matrix{<:Complex}` (red in REPL) instead of `Matrix{ComplexF64}` indicates boxing.
- AllocCheck.jl (`@check_allocs`) provides static analysis: annotate `apply_lindbladian!` and `step_along_trajectory!` to catch boxing at compile time, not just at runtime.

**How to prevent:**
- Every workspace struct that participates in a hot path must be parameterized on its concrete element type: `struct UnifiedWorkspace{T<:Complex} ... fields::Matrix{T} ... end`.
- When extracting data from abstractly-typed sources (`JumpOp.in_eigenbasis::Matrix{<:Complex}`), convert to the concrete type at construction time (as `KrylovWorkspace` already does), never in the hot loop.
- Use `@assert isconcretetype(T)` in workspace inner constructors during development to enforce this invariant.
- Write a `_verify_concrete_types(ws)` test helper that walks struct fields with `fieldtype()` and asserts concreteness.

**Which phase:** Workspace consolidation phase. Must be the gate criterion for merging workspace changes.

**Detection priority:** Run `@code_warntype apply_lindbladian!(ws, rho, config, ham)` after every workspace struct modification.

---

### CRIT-02: Union{T, Nothing} Fields Defeating Constant Propagation in Hot Paths

**What goes wrong:**
The existing config structs (`LiouvConfig`, `ThermalizeConfig`) have many `Union{T, Nothing}` fields (e.g., `a::Union{T, Nothing}`, `t0::Union{T, Nothing}`). Julia's compiler treats `Union{T, Nothing}` as a small union and handles it efficiently with "union splitting" -- but only when the access pattern is simple. In a hot path, accessing a `Union{T, Nothing}` field forces a branch on every access even if the value is always `T` at runtime. The `TrajectoryFramework` struct already works around this by extracting values into concrete-typed copies:

```julia
# Hot-path fields with concrete types (avoid accessing abstract-typed config/precomputed_data in step loop)
scaled_prefactor::Float64
sigma::Float64
transition::F
energy_labels::Vector{Float64}
```

When redesigning the config hierarchy (e.g., introducing `Config{S,D,DB,T}`), if these extraction patterns are not preserved, every `config.t0` access in a hot loop becomes a union-split branch.

**How to detect:**
- `@code_warntype` on `step_along_trajectory!` -- look for `Union{Float64, Nothing}` in inferred types.
- Allocation tests: `@test allocs == 0` for `step_along_trajectory!` will fail if union splitting causes boxing (happens when the union has more than 4 variants or when the compiler gives up).
- Performance regression: benchmark `step_along_trajectory!` before and after config redesign.

**How to prevent:**
- Maintain the "framework extraction" pattern: hot-path structs (`TrajectoryFramework`, `KrylovWorkspace`) copy values from config into concrete-typed fields at construction time.
- In the new config hierarchy, either (a) eliminate `Union{T, Nothing}` entirely by using required fields per domain variant, or (b) ensure no hot-path function ever reads a `Union{T, Nothing}` field directly.
- Consider using domain-specific config subtypes that only have the fields relevant to that domain (no Nothing fields needed).

**Which phase:** Type hierarchy redesign phase. This is a design decision, not a bug fix -- the new hierarchy must encode this rule.

---

### CRIT-03: Breaking Zero-Allocation Guarantees via Lazy Wrappers in mul!/BLAS Calls

**What goes wrong:**
The hot paths use `BLAS.gemm!('N', 'N', ...)` and `mul!` for zero-allocation matrix multiply. When deduplicating code (e.g., merging similar sandwich computations from `krylov_matvec.jl`), it is tempting to write:

```julia
# WRONG: creates lazy Adjoint wrapper, may allocate in BLAS path
mul!(out, A', B)  # A' creates Adjoint(A) wrapper
```

While `mul!(C, A', B)` *usually* dispatches to `gemm!('C', 'N', ...)` without allocation in recent Julia versions, the behavior depends on the exact types involved. The existing code carefully avoids this by using explicit BLAS.gemm! with character flags:

```julia
# CORRECT: zero allocation, explicit transpose flag
BLAS.gemm!('T', 'N', CT, L_op, rho, ZT, ws.tmp1)
```

During deduplication, replacing the 4 near-identical `_accumulate_sandwich*` functions (`_accumulate_sandwich!`, `_accumulate_sandwich_adj_L!`, `_accumulate_adjoint_sandwich!`, `_accumulate_adjoint_sandwich_adj_L!`) with a single parametric function risks introducing adjoint wrappers in the generic path.

**How to detect:**
- `test_krylov_matvec.jl` has `@test allocs == 0` for `apply_lindbladian!` across all domains (EnergyDomain, TimeDomain, TrotterDomain). These tests are the primary firewall.
- `test_allocation.jl` has `@test allocs == 0` for `step_along_trajectory!`.
- `@allocated` inside a function barrier (the codebase already uses this pattern -- see `test_allocation.jl` line 152).

**How to prevent:**
- When deduplicating sandwich functions, parameterize on the BLAS transpose character ('N', 'T', 'C') rather than using Julia `adjoint()` or `transpose()` wrappers:
  ```julia
  function _accumulate_sandwich!(out, L, rho, scalar, ws; transA='N', transB='T')
      BLAS.gemm!(transA, 'N', CT, L, rho, ZT, ws.tmp1)
      BLAS.gemm!('N', transB, CT, ws.tmp1, L, ZT, ws.tmp2)
      BLAS.axpy!(T(scalar), ws.tmp2, out)
  end
  ```
- Run allocation tests after every deduplication step, not just at the end.
- Use `@code_llvm` on the BLAS call to verify it maps to a single `ccall` without intermediate allocation.

**Which phase:** Deduplication phase. Must run allocation regression tests as a gate after each function merge.

---

### CRIT-04: BSON Deserialization Breakage from Struct Field Reordering

**What goes wrong:**
BSON.jl serializes Julia structs by field order. The existing codebase already handles one instance of this: `test_helpers.jl` has a custom `_load_test_hamiltonian` that manually parses legacy HamHam BSON files because the field layout changed. The `results.jl` module converts to Dicts for BSON safety specifically to avoid this problem.

During the refactor, if you:
1. Reorder fields in `LiouvConfig`, `ThermalizeConfig`, `HamHam`, `ConvergenceData`, or `ExperimentResult`
2. Add/remove fields from any BSON-serialized struct
3. Change type parameters (e.g., `LiouvConfig{D,T}` to `Config{S,D,DB,T}`)

...any existing `.bson` files (reference data in `test/reference/`, hamiltonians in `hamiltonians/`, saved experiment results) become unloadable.

**How to detect:**
- `test_results.jl` tests save/load round-trips.
- `test_regression.jl` loads frozen BSON reference data (`energy_dm_reference.bson`, `trotter_coherent_dm_reference.bson`).
- `test_helpers.jl` loads hamiltonian BSON files via `_load_test_hamiltonian`.
- All three will fail if struct layouts change without migration.

**How to prevent:**
- The Dict-based serialization in `results.jl` (already present via `_config_to_dict` / `_dict_to_experiment`) is the correct pattern. Extend it: never BSON-serialize parametric structs directly.
- For `HamHam`: the `_load_test_hamiltonian` function in `test_helpers.jl` already uses `BSON.parse` (raw bytes) + manual reconstruction. This must be updated when fields change. Keep the migration code versioned.
- For the new config hierarchy: write `_config_to_dict` for the new types FIRST, then migrate old configs by loading via the old dict path and saving via the new one.
- If renaming `LiouvConfig{D,T}` to `Config{S,D,DB,T}`, add a compatibility `_reconstruct_config` path that handles both old and new dict formats. The existing code in `results.jl` already uses `get()` with defaults for forward compatibility -- extend this pattern.

**Which phase:** Type hierarchy redesign phase. Must update serialization BEFORE changing struct definitions.

---

### CRIT-05: Closure Capture Causing Heap Allocation in Transition Functions

**What goes wrong:**
The `transition` function stored in `precomputed_data` is a closure that captures config parameters (beta, sigma, a, b). This closure is called inside the innermost loop of every hot path (`step_along_trajectory!`, `apply_lindbladian!`). Julia closures that capture mutable variables or variables whose types are not inferrable at compile time will allocate on each call.

The existing code handles this well -- `TrajectoryFramework` stores `transition::F` as a parametric type parameter, so the compiler knows the exact closure type and can inline. But during deduplication, if you extract a common "iterate over frequencies and apply weights" function that takes `transition` as a plain `Function` argument:

```julia
# WRONG: Function is abstract, kills inlining and causes allocation
function _iterate_frequencies(transition::Function, ...)

# CORRECT: parametric, compiler knows exact closure type
function _iterate_frequencies(transition::F, ...) where {F}
```

**How to detect:**
- `@code_warntype` on the extracted function -- `transition` should show as a concrete closure type, not `Function` or `Any`.
- `@allocated` tests catch the resulting allocations.
- Performance benchmark: closure boxing in the inner loop causes 10-100x slowdown because every call triggers dynamic dispatch.

**How to prevent:**
- Every function that receives a callable in a hot path must use `where {F}` parameterization: `function foo(f::F, ...) where {F}`.
- Never type-annotate callbacks as `::Function` in performance-critical code.
- The deduplication extracting shared iteration logic must preserve the parametric typing that `TrajectoryFramework{T,D,F,P}` already uses for the `transition` field.

**Which phase:** Deduplication phase. Applies to every shared function extracted from hot paths.

---

### CRIT-06: Thread Safety Regression from Shared Workspace Mutation

**What goes wrong:**
The multi-threaded trajectory execution (see `trajectories.jl` lines 538-567) creates per-thread workspaces (`ws_per_task`) and uses `BLAS.set_num_threads(1)` to prevent OpenBLAS internal threading from conflicting with Julia tasks. The `TrajectoryFramework` is immutable and shared across threads, while `TrajectoryWorkspace` is per-thread mutable.

During workspace consolidation, if the "read-only shared data" and "per-thread mutable scratch" are merged into a single struct, or if a "consolidate workspaces" step accidentally makes `TrajectoryFramework` mutable (e.g., adding a scratch buffer to it for "convenience"), the per-thread isolation breaks. This causes:
1. Data races (silent numerical corruption, non-deterministic results)
2. BLAS threading conflicts (segfaults or wrong results from concurrent gemm! calls on shared buffers)

**How to detect:**
- `test_threading.jl` tests bitwise determinism: `@test result1.rho_mean == result2.rho_mean` (not `isapprox` -- actual bitwise equality with `==`).
- `test_threading.jl` tests serial-threaded agreement: `@test isapprox(result_threaded.rho_mean, rho_ref; atol=1e-13)`.
- `test_workspace_independence.jl` tests workspace isolation.
- Run with `julia --check-bounds=yes -t 4` and verify determinism.

**How to prevent:**
- Maintain strict separation: framework/config structs are immutable (`struct`, not `mutable struct`) and shared. Workspace structs are mutable and per-thread.
- After consolidation, grep for any `mutable struct` that is shared across `Threads.@spawn` boundaries.
- The `BLAS.set_num_threads(1)` / restore pattern in `trajectories.jl` must survive any refactoring. Wrap it in a `_with_serial_blas(f)` helper that is used consistently.
- Add a threading stress test: 10,000 trajectories with 4+ threads, assert bitwise determinism.

**Which phase:** Workspace consolidation phase. Thread safety tests must be the merge gate.

---

## Moderate Pitfalls

Issues that cause test failures or performance regression but not silent data corruption.

---

### MOD-01: @allocated Tests Breaking from Testset-Lambda Interaction (Julia Bug #50796)

**What goes wrong:**
Julia has a known bug ([GitHub issue #50796](https://github.com/JuliaLang/julia/issues/50796)) where `@allocated` inside a `@testset` reports spurious allocations if any lambda function is defined elsewhere in the same testset scope. The existing codebase already works around this with function barriers:

```julia
# From test_allocation.jl line 150-164:
# Julia's @allocated in global/testset scope can show spurious allocations
# from boxing local variables; a function barrier ensures proper optimization.
function _measure_step_allocs(fw, ws, psi0)
    psi = copy(psi0)
    rng = Xoshiro(999)
    for _ in 1:100  # warmup
        step_along_trajectory!(psi, fw, ws, rng)
    end
    copyto!(psi, psi0)
    rng2 = Xoshiro(999)
    return @allocated step_along_trajectory!(psi, fw, ws, rng2)
end
allocs = _measure_step_allocs(fw, ws, psi0)
@test allocs == 0
```

During test restructuring, if allocation tests are moved into testsets that also contain lambda expressions, or if the function barrier pattern is lost, previously passing `@test allocs == 0` tests will spuriously fail.

**How to detect:**
- Tests that intermittently fail with small non-zero allocations (typically 16-80 bytes) in `@allocated` checks.
- Failures only when run as part of the full test suite, not in isolation.

**How to prevent:**
- Preserve the function barrier pattern for ALL `@allocated` tests: wrap the warmup+measure in a plain function, call it from testset.
- Never put `@allocated` in the same function scope as anonymous functions, `map(x -> ..., ...)`, or generator expressions.
- Consider migrating critical allocation tests to use [AllocCheck.jl](https://github.com/JuliaLang/AllocCheck.jl) (`@check_allocs`) for static analysis. This eliminates the flakiness entirely because it analyzes LLVM IR, not runtime behavior.
- Document this rule: "All allocation tests must use function barriers."

**Which phase:** Test restructure phase. Add as a review checklist item.

---

### MOD-02: Tolerance Drift When Sharing Numerical Fixtures Across Restructured Test Files

**What goes wrong:**
The current test suite uses shared fixtures from `test_helpers.jl` computed once at include time (line 146-149):

```julia
const TEST_SYSTEM = make_test_system()
const TEST_HAM = TEST_SYSTEM.hamiltonian
const TEST_JUMPS = TEST_SYSTEM.jumps
const TEST_GIBBS = TEST_SYSTEM.gibbs
```

Tolerance constants are also shared (line 80-82):

```julia
const TOL_EXACT = 1e-12          # machine precision identities
const TOL_QUADRATURE = 1e-6      # quadrature / discretization errors
TOL_DELTA(delta) = 5.0 * delta   # unraveling error, C = 5.0
```

When restructuring tests:
1. If fixtures are recomputed per-file instead of shared, floating-point non-associativity means the "same" Gibbs state may differ by ULP. Tests comparing against `TEST_GIBBS` with `atol=1e-12` will fail.
2. If tolerance constants are copied into individual test files instead of centralized, they will drift as developers adjust thresholds locally.
3. If `NUM_ENERGY_BITS`, `W0`, `T0` are changed in some test files but not others, cross-validation tests become meaningless.

**How to detect:**
- Tests that pass individually but fail when run together (fixture recomputation order matters).
- Tests that pass on one platform but fail on another (floating-point non-determinism from recomputation).
- Grep for hardcoded tolerance values: `grep -rn "atol=1e-" test/ | grep -v "TOL_"` should return zero hits (or only justified exceptions).

**How to prevent:**
- Keep ONE `test_helpers.jl` (or rename to `test_fixtures.jl`) that is included exactly once from `runtests.jl`. Never re-include it.
- All tolerance constants must be defined in the shared fixture file and referenced by name.
- Pin random seeds in fixture generation (already done: fixtures are deterministic from BSON-loaded Hamiltonians).
- The 146 instances of `isapprox|atol|rtol|TOL` across 16 test files should all reference named constants.

**Which phase:** Test restructure phase. Run full test suite after every restructuring step, not just the modified file.

---

### MOD-03: Include Order Breakage from File Reorganization

**What goes wrong:**
Julia's `include()` is textual inclusion -- files are evaluated in order in `QuantumFurnace.jl` (lines 102-129). The current include chain has implicit dependencies:

```
constants.jl        -> (no deps)
hamiltonian.jl      -> constants.jl (HamHam uses X, Y, Z from constants)
trotter_domain.jl   -> hamiltonian.jl (TrottTrott uses HamHam)
structs.jl          -> (defines AbstractDomain, configs, workspaces)
qi_tools.jl         -> (standalone utilities)
...
krylov_workspace.jl -> jump_workers.jl, coherent.jl (uses _precompute_data, _precompute_coherent_total_B)
krylov_matvec.jl    -> krylov_workspace.jl (uses KrylovWorkspace)
krylov_eigsolve.jl  -> krylov_matvec.jl (uses apply_lindbladian!)
...
results.jl          -> (uses everything, must be last)
```

Renaming or splitting files (e.g., extracting a shared `workspace.jl` from `krylov_workspace.jl` + `trajectories.jl`) will fail if the new file is included before its dependencies. Julia will not error on an undefined function reference at include time -- it only errors when the function is *called*. This means `using QuantumFurnace` succeeds but tests crash at runtime.

**How to detect:**
- `test_compilation.jl` catches some load-order failures (exists in the test suite).
- Run `using QuantumFurnace` in a fresh Julia session after every file rename/move.
- Add a smoke test that calls `methods(apply_lindbladian!)` to verify method tables are populated.

**How to prevent:**
- Draw the dependency graph before reorganizing. Map every function/type to the file that defines it and every file to the functions/types it references.
- Move files in topological order: first rename files that have no dependents, then rename files whose dependents have already been updated.
- After every file move, run `using QuantumFurnace` in a fresh Julia REPL.
- Consider adding an "all modules load" smoke test:
  ```julia
  @testset "All exports defined" begin
      for name in names(QuantumFurnace)
          @test isdefined(QuantumFurnace, name)
      end
  end
  ```

**Which phase:** File reorganization phase. Do in small, testable increments -- one file rename per commit.

---

### MOD-04: precomputed_data::Any Defeating Type Inference When Consolidating Workspaces

**What goes wrong:**
The `TrajectoryFramework` has `precomputed_data::Any` (line 84 of `trajectories.jl`). This was an intentional choice because the `NamedTuple` type varies by domain (different fields for Energy vs Time vs Trotter). The hot path works around this by extracting needed values into concrete-typed fields at construction time.

During consolidation, if a unified workspace struct tries to store precomputed data for multiple domains in the same field, the `::Any` typing pattern spreads. If hot-path code starts accessing `ws.precomputed_data.transition` instead of `ws.transition`, the compiler cannot infer the return type, causing dynamic dispatch on every frequency iteration.

The `KrylovWorkspace{T, PD}` already has the correct pattern -- parameterized on the precomputed data type `PD<:NamedTuple`. The `PD` type parameter lets the compiler know the exact NamedTuple layout.

**How to detect:**
- `@code_warntype` on `step_along_trajectory!` or `apply_lindbladian!` -- red `Any` annotations on precomputed_data access.
- Allocation regression in the hot path.

**How to prevent:**
- When consolidating, either (a) keep the `PD` parametrization: `struct UnifiedWorkspace{T, PD<:NamedTuple}`, or (b) extract all hot-path values into concrete fields at construction time (the `TrajectoryFramework` pattern).
- Never access `precomputed_data` fields inside a loop that runs per-frequency or per-trajectory-step.
- The `KrylovWorkspace{T,PD}` pattern is the model to follow.

**Which phase:** Workspace consolidation phase. Design decision at the start, verified by `@code_warntype` after.

---

### MOD-05: Over-Parameterization Making the Config Hierarchy Unusable

**What goes wrong:**
The proposed new config type `Config{S,D,DB,T}` has 4 type parameters: Simulation, Domain, DetailedBalance variant, and Float type. This is already borderline for Julia ergonomics. If additional parameters are added (e.g., `Config{S,D,DB,T,LC}` for linear combination variant), the type becomes unwieldy:
1. Every function signature must carry `where {S,D,DB,T}` constraints.
2. Method specialization explodes: 4 simulations x 4 domains x 2 DB variants x 1 float type = 32 potential specializations per method.
3. Error messages become incomprehensible: `MethodError: no method matching foo(::Config{ThermalizeSimulation, TimeDomain, KMS, Float64})`.
4. Construction becomes verbose.

**How to detect:**
- User frustration at the API level.
- Method ambiguity errors when new methods are added.
- Compile time increases from excessive specialization.

**How to prevent:**
- Limit to 2-3 type parameters maximum. Use the domain singleton for dispatch (already a type parameter). Use DB as a type tag only if GNS vs KMS affects hot-path method selection; otherwise, make it a runtime field.
- The current pattern of 4 separate structs (`LiouvConfig`, `LiouvConfigGNS`, `ThermalizeConfig`, `ThermalizeConfigGNS`) is verbose but has no parametric complexity. The new hierarchy should reduce struct count without adding more than one new type parameter.
- Consider trait-based dispatch (e.g., `is_gns(config)::Bool` or `detailed_balance_type(config)::Symbol`) instead of type-parameterization for the DB variant. Reserve type parameters for things that affect compilation (element type, domain).
- Prototype the API before implementing -- write the function signatures first, then design the types.

**Which phase:** Type hierarchy redesign phase.

---

### MOD-06: Losing @inline Annotations During Function Extraction

**What goes wrong:**
Several hot-path functions have explicit `@inline` annotations:
- `_krylov_oft!` (krylov_matvec.jl:11)
- `_norm2` (trajectories.jl:344)
- `_accumulate_density_matrix!` (trajectories.jl:346)
- `_prefactor_view` (nufft.jl:93)
- `_apply_coherent_unitary!` (jump_workers.jl:145)

When extracting these into shared utility files or deduplicating them, if the `@inline` annotation is lost, the compiler may choose not to inline them. For small functions called millions of times in the inner loop (like `_norm2`), this adds function call overhead and prevents SIMD and constant propagation across the call boundary.

**How to detect:**
- `@code_native` on the hot-path caller -- if the small function is not inlined, you will see a `call` instruction to it.
- Benchmark regression.

**How to prevent:**
- When moving a function, copy the `@inline` annotation with it.
- Add a comment `# @inline: required for hot path` to every `@inline`-annotated function so reviewers know it is intentional.
- After deduplication, verify inlining with `@code_native` on the callers.

**Which phase:** Deduplication phase. Trivial to prevent with attention; hard to debug if missed.

---

## Minor Pitfalls

Issues that cause inconvenience or technical debt but not correctness problems.

---

### MIN-01: Export List Staleness After File Reorganization

**What goes wrong:**
`QuantumFurnace.jl` has a large export list (lines 29-100) organized by functional area. After reorganizing files, exported symbols may be defined in different files than expected, or new symbols from consolidated code may not be exported. Downstream code using `using QuantumFurnace` will silently lose access to symbols that were previously exported.

**How to prevent:**
- After every file reorganization, run `test_aqua.jl` (which checks for unbound exports).
- Maintain the export list grouped by functional area with comments indicating which file defines each group.

**Which phase:** File reorganization phase.

---

### MIN-02: Git History Loss from File Renames

**What goes wrong:**
If files are renamed AND modified in the same commit, `git log --follow` cannot track the rename. This makes it impossible to trace the history of performance-critical code back to when allocation optimizations were introduced.

**How to prevent:**
- Rename files in pure-rename commits (no content changes). Then modify content in subsequent commits.
- Use `git mv old.jl new.jl` explicitly.
- Verify with `git diff --stat` that the rename commit shows `rename` not `delete+add`.

**Which phase:** File reorganization phase. One commit per rename.

---

### MIN-03: Test Runtime Blowup from Fixture Duplication

**What goes wrong:**
The test suite computes expensive fixtures once at `include("test_helpers.jl")` time: `make_test_system()` builds full 4-qubit and 3-qubit Hamiltonians, computes eigenstates, creates jump operators, builds Trotter objects. This takes several seconds. If test restructuring causes fixtures to be recomputed per-file or per-testset, the test suite runtime could increase dramatically.

**How to prevent:**
- Keep the "compute once, share everywhere" pattern via module-level `const` values.
- If splitting into sub-test-modules, pass fixtures as arguments rather than recomputing.
- Time the full test suite before and after restructuring; flag changes that increase runtime by >20%.

**Which phase:** Test restructure phase.

---

## Phase-Specific Warning Summary

| Phase | Likely Pitfall | Severity | Key Mitigation |
|-------|---------------|----------|----------------|
| Type hierarchy redesign | CRIT-04: BSON breakage from struct changes | Critical | Update serialization layer FIRST, before changing any struct |
| Type hierarchy redesign | CRIT-02: Union{T,Nothing} in hot paths | Critical | Maintain framework extraction pattern for hot-path values |
| Type hierarchy redesign | MOD-05: Over-parameterization | Moderate | Cap at 2-3 type params; prototype signatures first |
| Workspace consolidation | CRIT-01: Abstract field boxing | Critical | Parameterize on concrete element type; test with @allocated |
| Workspace consolidation | CRIT-06: Thread safety regression | Critical | Immutable shared data, mutable per-thread scratch |
| Workspace consolidation | MOD-04: precomputed_data::Any spreading | Moderate | Parameterize on NamedTuple type OR extract to concrete fields |
| Deduplication | CRIT-03: Lazy wrapper allocations in BLAS | Critical | Use BLAS.gemm! with char flags, not adjoint wrappers |
| Deduplication | CRIT-05: Closure capture allocation | Critical | Always use `where {F}` for callable parameters |
| Deduplication | MOD-06: Lost @inline annotations | Moderate | Copy annotations; verify with @code_native |
| Test restructure | MOD-01: @allocated + testset interaction | Moderate | Function barrier pattern; consider AllocCheck.jl |
| Test restructure | MOD-02: Tolerance drift from split fixtures | Moderate | Single fixture file; named constants; never hardcode atol |
| Test restructure | MIN-03: Runtime blowup from re-fixture | Minor | Compute once, share everywhere |
| File reorganization | MOD-03: Include order breakage | Moderate | Draw dependency graph; one rename per commit; smoke test |
| File reorganization | MIN-01: Stale export list | Minor | Run Aqua.jl after every reorganization |
| File reorganization | MIN-02: Git history loss | Minor | Pure rename commits; use git mv |

---

## Recommended Verification Gates Per Phase

### After Type Hierarchy Redesign
1. `test_results.jl` passes (BSON round-trip for both old and new format)
2. `test_regression.jl` passes (frozen reference data loads)
3. `_load_test_hamiltonian` in `test_helpers.jl` still loads hamiltonian BSON files
4. `@code_warntype` on `step_along_trajectory!` shows no red annotations
5. All 600+ existing tests pass with no tolerance changes

### After Workspace Consolidation
6. `test_allocation.jl` passes (all allocation bounds)
7. `test_krylov_matvec.jl` passes (`@test allocs == 0` for all domains)
8. `test_threading.jl` passes (bitwise determinism with `==`)
9. `test_workspace_independence.jl` passes
10. No `mutable struct` shared across `Threads.@spawn` boundaries

### After Deduplication
11. `test_allocation.jl` passes (unchanged thresholds)
12. `test_krylov_matvec.jl` passes (unchanged `allocs == 0`)
13. Benchmark: no >5% regression in `step_along_trajectory!` or `apply_lindbladian!`
14. `@code_warntype` on all extracted functions shows concrete types
15. All `@inline` annotations preserved (grep count matches pre-dedup)

### After File Reorganization
16. `using QuantumFurnace` succeeds in fresh Julia REPL
17. `test_aqua.jl` passes (no unbound exports)
18. All 600+ tests pass
19. `git log --follow` works for renamed files (verify for 3+ key files)

### After Test Restructure
20. All 600+ tests pass
21. Test suite runtime within 20% of pre-refactor baseline
22. No hardcoded tolerances outside `test_helpers.jl` (grep verification)
23. All `@allocated` tests wrapped in function barriers (grep verification)

---

## Sources

- [Julia Performance Tips -- Abstract Type Fields and Parametric Types](https://docs.julialang.org/en/v1/manual/performance-tips/) (HIGH confidence)
- [Julia Issue #50796 -- @allocated Spurious Allocations in @testset](https://github.com/JuliaLang/julia/issues/50796) (HIGH confidence)
- [AllocCheck.jl -- Static Allocation Analysis](https://github.com/JuliaLang/AllocCheck.jl) (HIGH confidence)
- [Julia Discourse -- Abstract Type Field Boxing Allocations](https://discourse.julialang.org/t/approach-to-avoid-memory-allocation-with-struct-with-abstract-type-fields/24822) (MEDIUM confidence)
- [Julia Discourse -- Allocations for Abstract Field of Struct](https://discourse.julialang.org/t/allocations-for-abstract-field-of-struct/96088) (MEDIUM confidence)
- [BSON.jl -- Serialization by Field Order](https://github.com/JuliaIO/BSON.jl) (MEDIUM confidence)
- [Julia Discourse -- Circular Dependencies in include()](https://discourse.julialang.org/t/circular-dependency/35571) (MEDIUM confidence)
- Direct codebase analysis of QuantumFurnace.jl src/ and test/ directories (HIGH confidence)
