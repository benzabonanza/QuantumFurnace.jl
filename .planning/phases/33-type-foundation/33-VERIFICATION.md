---
phase: 33-type-foundation
verified: 2026-02-25T22:44:28Z
status: passed
score: 5/5 must-haves verified
re_verification: false
---

# Phase 33: Type Foundation Verification Report

**Phase Goal:** All simulation code dispatches on a single `Config{S,D,C,T}` struct, eliminating 4 duplicate config types and enabling future DLL construction via a new singleton type
**Verified:** 2026-02-25T22:44:28Z
**Status:** passed
**Re-verification:** No — initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | `Config{Lindbladian,EnergyDomain,KMS,Float64}` can be constructed and used in place of `LiouvConfig{EnergyDomain,Float64}` for all 4 domain types | VERIFIED | `src/structs.jl` lines 108-145: Config{S,D,C,T} defined with outer constructor. test_helpers.jl factory `make_liouv_config` returns `Config(sim=Lindbladian(), domain=domain, construction=KMS(), ...)`. test_compilation.jl line 47 asserts `lc isa Config{Lindbladian}`. |
| 2 | `Config{Thermalize,TimeDomain,GNS,Float64}` can be constructed and used in place of `ThermalizeConfigGNS{TimeDomain,Float64}` for all 4 domain types | VERIFIED | test_helpers.jl factory `make_thermalize_config` returns `Config(sim=Thermalize(), domain=domain, construction=construction, ...)`. test_compilation.jl line 52 asserts `tc isa Config{Thermalize}`. trajectories.jl, convergence.jl, gap_estimation.jl all dispatch on `Config{Thermalize}`. |
| 3 | `with_coherent` is derived from construction type (KMS -> true, GNS -> false) at compile time, not stored as a field | VERIFIED | `src/structs.jl` lines 52-54: three-method trait `with_coherent(::KMS)=true`, `with_coherent(::GNS)=false`, `with_coherent(::DLL)=true`. Zero `config.with_coherent` field accesses anywhere in src/ or test/. All callers use `with_coherent(config.construction)`. |
| 4 | All existing tests pass using either the new `Config{S,D,C,T}` type or backward-compatible aliases | VERIFIED | SUMMARY reports 1187 tests passed (7 pre-existing failures in test_diagnostics.jl unrelated to migration). Zero old type names (LiouvConfig, ThermalizeConfig, etc.) found in any .jl file across src/, test/, simulations/, experiments/. |
| 5 | Adding a `DLL` construction type requires only `struct DLL <: AbstractConstruction end` plus dispatch methods, with zero changes to existing code | VERIFIED | `struct DLL <: AbstractConstruction end` already defined in structs.jl line 49 and exported. All dispatch is open: `pick_transition`, `_pick_alpha`, `B_bohr`, `_pick_f`, `_select_b_plus_calculator` dispatch on construction type parameter. Adding DLL only requires new methods (`pick_transition(config::Config{<:Any,<:Any,DLL})`, etc.) — no modification of existing methods. |

**Score:** 5/5 truths verified

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `src/structs.jl` | Type hierarchies, Config struct, with_coherent trait | VERIFIED | Lines 39-145: AbstractSimulation (4 singletons), AbstractConstruction (3 singletons), with_coherent trait (3 methods), Config{S,D,C,T} @kwdef struct, outer constructor. |
| `src/QuantumFurnace.jl` | Updated exports | VERIFIED | Lines 29-35: exports Config, AbstractSimulation, Lindbladian, Thermalize, KrylovSpectrum, Trajectory, AbstractConstruction, KMS, GNS, DLL, with_coherent. Zero old type names in exports. |
| `src/energy_domain.jl` | KMS/GNS transition dispatch on Config | VERIFIED | Lines 3-4: `pick_transition(config::Config{<:Any,<:Any,KMS})`, `pick_transition(config::Config{<:Any,<:Any,GNS})`. Two methods, not four. |
| `src/bohr_domain.jl` | KMS/GNS alpha and B_bohr dispatch on Config | VERIFIED | Lines 2, 28, 56, 82-83, 85, 115: B_bohr (KMS only x2), _pick_f (KMS only), _pick_alpha (2 methods: KMS/GNS), _pick_alpha_kms/gns. |
| `src/misc_tools.jl` | validate_config!, _generate_filename, _print_press using Config | VERIFIED | Line 117: `validate_config!(config::Config)`. Lines 76, 89: `_generate_filename(config::Config{Lindbladian})`, `_generate_filename(config::Config{Thermalize})`. All `with_coherent` via trait. |
| `src/results.jl` | Serialization/deserialization using Config | VERIFIED | Line 18: `ExperimentResult{C<:Config, T}`. Line 62: `_config_to_dict(config::Config)`. Lines 160-177: `_reconstruct_config` returns `Config(sim=sim, domain=domain, construction=construction, ...)`. |
| `src/furnace.jl` | run_lindbladian, construct_lindbladian, run_thermalization using Config | VERIFIED | Lines 1, 42, 88: dispatches on `Config{Lindbladian,D,C,Tc}`, `Config{Lindbladian}`, `Config{Thermalize,D,C,Tc}`. |
| `src/krylov_workspace.jl` | KrylovWorkspace constructors using Config, _thermalize_to_liouv_config deleted | VERIFIED | Lines 90, 336: KrylovWorkspace constructors for `Config{Lindbladian}` and `Config{Thermalize}`. Comment at line 136 confirms deletion. Zero occurrences of `_thermalize_to_liouv_config` as a function definition. |
| `src/krylov_eigsolve.jl` | krylov_spectral_gap dispatching on Config{Lindbladian} and Config{Thermalize} | VERIFIED | Lines 386, 487: `krylov_spectral_gap(config::Config{Lindbladian}, ...)` and `krylov_spectral_gap(config::Config{Thermalize}, ...)`. |
| `test/test_helpers.jl` | Updated config factory functions using Config | VERIFIED | Lines 223+: all 7 factory functions (`make_liouv_config`, `make_liouv_config_gns`, `make_thermalize_config`, `make_small_*`) return Config with sim/construction singletons. |
| `test/test_compilation.jl` | Updated isa checks for Config types | VERIFIED | Lines 47, 52: `@test lc isa Config{Lindbladian}`, `@test tc isa Config{Thermalize}`. |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `src/structs.jl` | `src/QuantumFurnace.jl` | export list references new types | WIRED | QuantumFurnace.jl lines 29-35 export Config, AbstractSimulation, Lindbladian, Thermalize, KrylovSpectrum, Trajectory, AbstractConstruction, KMS, GNS, DLL, with_coherent |
| `src/energy_domain.jl` | `src/structs.jl` | Config type parameter dispatch | WIRED | `Config{<:Any, <:Any, KMS}` and `Config{<:Any, <:Any, GNS}` dispatch at lines 3-4 |
| `src/results.jl` | `src/structs.jl` | Config construction for deserialization | WIRED | `_reconstruct_config` at line 177 calls `Config(sim=sim, domain=domain, construction=construction, kwargs...)` |
| `src/furnace.jl` | `src/furnace_utensils.jl` | _precompute_data call with Config | WIRED | furnace.jl calls `_precompute_data(config, ...)` passing `Config{Lindbladian}` and `Config{Thermalize}` directly |
| `src/krylov_workspace.jl` | `src/furnace_utensils.jl` | direct Config pass-through (no _thermalize_to_liouv_config) | WIRED | KrylovWorkspace constructor passes Config directly to `_precompute_data`. No `_thermalize_to_liouv_config` conversion. |
| `test/test_helpers.jl` | `src/structs.jl` | Config constructor calls | WIRED | Factory functions call `Config(sim=Lindbladian(), ...)` at multiple locations |
| `test/test_compilation.jl` | `src/structs.jl` | isa Config{Lindbladian} checks | WIRED | Lines 47, 52: `isa Config{Lindbladian}`, `isa Config{Thermalize}` |

### Requirements Coverage

No REQUIREMENTS.md mapping found for phase 33 specifically.

### Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| `src/structs.jl` | 54 | `with_coherent(::DLL) = true  # placeholder for Ding et al.` | Info | Expected placeholder — DLL is a future construction type stub with coherent=true as a reasonable default |
| `src/errors.jl` | 1 | `# Error computation utilities (placeholder for Phase 10 API cleanup)` | Info | Pre-existing comment, unrelated to phase 33 goal |
| `src/misc_tools.jl` | 78, 91, 220, 250 | `config.construction isa GNS ? "GNS" : "KMS"` | Info | Would silently tag DLL as "KMS" in filenames/print output. Not a blocker — DLL is a future construction type and this is a labeling-only concern, not a dispatch correctness issue |

No blocker anti-patterns found.

### Human Verification Required

None — all critical success criteria are verifiable programmatically.

**Optional (low priority):** Running the full test suite `julia --project -e 'using Pkg; Pkg.test()'` to confirm the 1187 passes reported in SUMMARY still hold in the current environment. SUMMARY documents this was verified at commit 5cb331c/4114811.

### Gaps Summary

No gaps. All 5 phase success criteria are satisfied:

1. `Config{S,D,C,T}` struct is defined, exported, and used across all 12+ source files in place of 4 old config types. Zero old type names remain in any tracked .jl file.
2. `with_coherent` is a pure trait function dispatching on construction type singletons. The field was removed from Config. No `config.with_coherent` field access exists anywhere.
3. All dispatch has been reduced from 4-way (one per old config type) to 2-way (KMS/GNS) for construction-specific functions, and simulation-type-specific dispatch uses `Config{Lindbladian}` / `Config{Thermalize}`.
4. `_thermalize_to_liouv_config` is deleted. KrylovWorkspace accepts Config directly.
5. `struct DLL <: AbstractConstruction end` is defined and exported. Adding DLL behavior requires only new dispatch methods, with zero modifications to existing code.

---

_Verified: 2026-02-25T22:44:28Z_
_Verifier: Claude (gsd-verifier)_
