# Phase 15: Data Architecture - Research

**Researched:** 2026-02-16
**Domain:** Experiment result serialization via BSON in Julia
**Confidence:** HIGH

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions
- Single flat struct `ExperimentResult{C,T}` parameterized on config type `C<:AbstractConfig` and element type `T`
- Reuse existing result types (TrajectoryResult, LindbladianResult, etc.) where they already carry the right data -- embed them rather than duplicating fields
- Full config object embedded in the result -- everything needed to re-run the experiment is stored
- Config type parameter `C` preserves KMS vs GNS distinction at the type level (dispatch-friendly)
- Descriptive file names: e.g., `kms_n4_beta10_trotter_20260216.bson`
- Default output directory `results/` in project root, with optional path override
- Subdirectories by construction type: `results/kms/`, `results/approx_gns/` (and later `results/ding/`)
- Overwriting existing files is allowed silently -- re-running replaces the old result
- Companion `.txt` file saved alongside each `.bson` with key parameters for browsing without Julia
- Metadata scope: seed, thread count, Julia version, timestamp, git commit hash (auto-captured), Hamiltonian parameters (n, J, h) plus coefficient matrices and term vectors, wall-clock time, total trajectory count
- No package version field -- git hash is sufficient
- No schema version number -- keep it simple
- Missing fields on load should fill with `nothing`/defaults (best-effort partial load)
- Experiments are cheap enough to re-run if schema diverges significantly

### Claude's Discretion
- Whether to use dedicated save_experiment/load_experiment wrapper functions vs raw BSON.bson/BSON.load
- Whether ExperimentResult wraps TrajectoryResult as a field or flattens its contents
- Exact companion .txt format and content
- Auto-capture mechanism for git hash (LibGit2 vs shell command)

### Deferred Ideas (OUT OF SCOPE)
- `results/ding/` directory for Ding et al. (2024) construction -- future milestone
- Index/manifest file listing all experiments in a sweep -- could be useful for Phase 18 but not needed for the data layer itself
</user_constraints>

## Summary

Phase 15 builds the data persistence layer for QuantumFurnace experiment results. The core deliverable is an `ExperimentResult{C,T}` struct that bundles a config, trajectory result, metadata, and Hamiltonian provenance into a single serializable object, plus `save_experiment`/`load_experiment` wrapper functions that handle BSON serialization, directory creation, companion text files, and metadata auto-capture.

The primary technical challenge is BSON.jl 0.3's behavior with parametric types, abstract-typed fields, and LinearAlgebra wrapper types like `Hermitian`. The codebase already has extensive experience with BSON pitfalls (the legacy HamHam loader in `misc_tools.jl` and `test_helpers.jl` uses `BSON.parse` + manual field reconstruction to work around deserialization failures). The approach for ExperimentResult should avoid these problems by design: store concrete types (unwrap `Hermitian` to `Matrix` before saving), avoid abstract-typed fields in the BSON payload (convert domain singletons to strings), and use `BSON.bson`/`BSON.load` with explicit Dict wrapping rather than relying on automatic struct serialization.

**Primary recommendation:** Use dedicated `save_experiment`/`load_experiment` wrapper functions that convert ExperimentResult to/from a plain Dict for BSON serialization, avoiding BSON.jl's known issues with parametric structs and abstract type fields.

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| BSON.jl | 0.3.9 | Binary serialization of experiment results | Already a project dependency; used for Hamiltonians and regression data |
| LibGit2 | stdlib | Auto-capture git commit hash at save time | Julia stdlib, no extra dependency; provides `head_oid` for programmatic commit hash |
| Dates | stdlib | Timestamp generation for metadata | Julia stdlib; `Dates.now()` and `Dates.format` |
| Pkg | stdlib | Project root path detection | Already imported in module; `Pkg.project().path` for finding `results/` |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| Printf | stdlib | Companion `.txt` file formatting | Already used throughout codebase for `@printf` |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| BSON.jl | JLD2.jl | JLD2 handles parametric types better, but BSON is already a dep and user decided against JLD2 |
| LibGit2 | `read(\`git rev-parse HEAD\`, String)` | Shell command is simpler but fails if git not on PATH; LibGit2 is a stdlib with no external dependency |

**Installation:**
No new dependencies needed. BSON, Pkg, LibGit2, Dates, Printf are all already available.

## Architecture Patterns

### Recommended Project Structure
```
src/
    results.jl              # ExperimentResult struct + save/load functions
    QuantumFurnace.jl       # Add include("results.jl") + exports
results/
    kms/                    # KMS experiment results (created on first save)
    approx_gns/             # GNS experiment results (created on first save)
```

### Pattern 1: Dict-Based BSON Serialization (Recommended)
**What:** Convert ExperimentResult to a plain Dict{Symbol,Any} before calling `BSON.bson`, and reconstruct from Dict on load.
**When to use:** Always -- avoids all BSON.jl struct serialization pitfalls.
**Why:** The codebase already uses this pattern for regression data (`generate_references.jl` saves `Dict(:rho => ..., :delta => ..., :domain => string(typeof(domain)))`) and it loads cleanly with `BSON.load`. The HamHam legacy loader (`misc_tools.jl`) demonstrates what happens when you rely on automatic struct serialization -- you need `BSON.parse` + `raise_recursive` + manual field extraction.
**Example:**
```julia
# Save: convert to Dict
function save_experiment(result::ExperimentResult, path::String)
    d = _experiment_to_dict(result)
    mkpath(dirname(path))
    BSON.bson(path, d)
    _write_companion_txt(result, replace(path, ".bson" => ".txt"))
end

# Load: reconstruct from Dict
function load_experiment(path::String)
    d = BSON.load(path)
    return _dict_to_experiment(d)
end
```

### Pattern 2: Hermitian Unwrapping for Serialization
**What:** Store `Hermitian{T, Matrix{T}}` fields as plain `Matrix{T}` in the BSON Dict, re-wrap on load.
**When to use:** Any field that holds a `Hermitian` (like `TrajectoryResult.rho_mean` or Gibbs states).
**Why:** BSON.jl can struggle with LinearAlgebra wrapper types. The existing regression test data stores density matrices as `Matrix(rho_dm)` (see `generate_references.jl` line 47). Unwrapping to `Matrix` and re-wrapping to `Hermitian` on load is trivial and avoids deserialization issues.
**Example:**
```julia
# In _experiment_to_dict:
:rho_mean => Matrix(result.trajectory_result.rho_mean)

# In _dict_to_experiment:
rho = Hermitian(d[:rho_mean])  # or just Matrix if not Hermitian-typed
```

### Pattern 3: Domain/Config Type Serialization via String Tags
**What:** Store domain type and config type as string tags, reconstruct on load via a lookup table.
**When to use:** For the `domain` field (a singleton like `TrotterDomain()`) and config type discrimination (KMS vs GNS).
**Why:** BSON cannot reliably round-trip zero-field singleton structs or parametric type parameters. String tags are robust.
**Example:**
```julia
# Save
:domain => string(typeof(config.domain))  # "TrotterDomain"
:config_type => config isa LiouvConfigGNS ? "GNS" : "KMS"

# Load -- reconstruct domain
domain = _string_to_domain(d[:domain])  # Dict("TrotterDomain" => TrotterDomain(), ...)
```

### Pattern 4: Embed TrajectoryResult as a Field (Not Flattened)
**What:** Store TrajectoryResult as a nested sub-dict within the ExperimentResult Dict.
**When to use:** For the trajectory result data (rho_mean, n_trajectories, seed, times, measurements_mean).
**Why:** TrajectoryResult has only 5 fields, all concrete types (Matrix, Int, Union{Nothing,Vector}, Union{Nothing,Matrix}). Embedding it as a named sub-dict preserves the logical grouping and makes reconstruction straightforward. Flattening would create name collisions with top-level metadata fields (both have `seed`).
**Example:**
```julia
# In the BSON Dict
:trajectory => Dict(
    :rho_mean => Matrix(result.trajectory_result.rho_mean),
    :n_trajectories => result.trajectory_result.n_trajectories,
    :seed => result.trajectory_result.seed,
    :times => result.trajectory_result.times,
    :measurements_mean => result.trajectory_result.measurements_mean,
)
```

### Anti-Patterns to Avoid
- **Direct struct serialization with `@save`:** BSON.jl's `@save` macro stores the full type tag including parametric parameters. If the struct definition changes (field added/removed), loading fails silently or errors. Always go through Dict.
- **Storing closures or functions:** The `precomputed_data` NamedTuple contains closures (the `transition` function). Never include precomputed_data or TrajectoryFramework in saved results -- they cannot round-trip.
- **Storing the full HamHam object:** HamHam has 14 fields including a `Dict{T, Vector{CartesianIndex{2}}}` (bohr_dict) and `Hermitian` gibbs -- both problematic for BSON. Store only the reproduction-relevant subset: num_qubits, base_coeffs, disordering_coeffs, base_terms, periodic.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Binary serialization | Custom binary format | BSON.jl via Dict | Already a dependency, well-tested for Julia arrays/matrices |
| Git commit hash | Shell `git rev-parse` | `LibGit2.head_oid(repo)` | Works without git on PATH, pure Julia stdlib |
| Directory creation | Manual checks | `mkpath(dirname(path))` | Julia stdlib handles nested creation and existing dirs |
| Timestamp formatting | `string(now())` | `Dates.format(now(), dateformat"yyyy-mm-dd_HH:MM:SS")` | Reproducible, sortable format |

**Key insight:** The Dict-based serialization pattern is already proven in this codebase (regression references, Hamiltonian storage). Building a custom serializer or relying on BSON's automatic struct handling would be strictly worse.

## Common Pitfalls

### Pitfall 1: BSON Cannot Round-Trip Abstract-Typed Fields
**What goes wrong:** A struct field typed as `config::AbstractConfig{D,T}` will not deserialize correctly because BSON stores the abstract type tag but cannot reconstruct the concrete subtype.
**Why it happens:** BSON.jl 0.3 serializes Julia types by their full type path. When loading, it tries to reconstruct the type, but abstract type parameters and Union fields cause ambiguity.
**How to avoid:** Convert configs to Dict representation before saving. On load, reconstruct using the concrete config constructor (LiouvConfig, ThermalizeConfig, etc.) from the Dict fields. The `config_type` string tag ("KMS" or "GNS") plus the `domain` string tag suffice to pick the right constructor.
**Warning signs:** `MethodError` or `TypeError` on `BSON.load` when the file contains a parametric struct with abstract fields.

### Pitfall 2: Hermitian/Diagonal Wrapper Types in BSON
**What goes wrong:** `Hermitian{ComplexF64, Matrix{ComplexF64}}` or `Diagonal{...}` types may not deserialize to the same wrapper type. BSON may return a plain Matrix or fail with a type mismatch.
**Why it happens:** These are thin wrappers around `Matrix` that BSON doesn't always handle transparently.
**How to avoid:** Always unwrap to `Matrix()` before saving. Re-wrap to `Hermitian()` on load if needed. The existing `generate_references.jl` already does `Matrix(rho_dm)`.
**Warning signs:** Type assertion failures when loading, or fields that are `Matrix` when you expected `Hermitian`.

### Pitfall 3: Domain Singleton Struct Serialization
**What goes wrong:** `TrotterDomain()` is a zero-field singleton struct. BSON can serialize it but the round-trip may fail if the module/namespace is wrong.
**Why it happens:** BSON stores the type as `QuantumFurnace.TrotterDomain` -- on load, it needs the module in scope to resolve this.
**How to avoid:** Store as a string `"TrotterDomain"` and reconstruct from a lookup Dict on load. This is what `generate_references.jl` already does: `string(typeof(domain))`.
**Warning signs:** `UndefVarError` or namespace resolution errors on load.

### Pitfall 4: Config Fields with `Union{T, Nothing}` and `Tuple{T,T} | Tuple{Nothing,Nothing}`
**What goes wrong:** The many `Union{T, Nothing}` fields in LiouvConfig/ThermalizeConfig can cause BSON to lose type information, loading `nothing` as the wrong type or losing the numeric precision.
**Why it happens:** BSON's internal representation doesn't distinguish between typed `nothing` and untyped `nothing`.
**How to avoid:** When converting config to Dict, store each field as-is (BSON handles scalars and nothing fine in Dict context). On reconstruction, the `@kwdef` constructor with default values handles `nothing` fields naturally.
**Warning signs:** Type instability in loaded configs; `InexactError` when converting loaded values.

### Pitfall 5: `BSON.load` Namespace Resolution
**What goes wrong:** `BSON.load(path)` in a different module context fails to resolve custom types.
**Why it happens:** Default namespace is `Main`. If ExperimentResult is defined in QuantumFurnace, BSON needs to know that.
**How to avoid:** Since we use Dict-based serialization (no custom types in the BSON), this is automatically avoided. The Dict contains only standard Julia types (String, Float64, Int, Matrix, Vector, Nothing).
**Warning signs:** `KeyError` or type resolution errors when loading BSON files from scripts outside the QuantumFurnace module.

### Pitfall 6: Forgetting `results/` is Gitignored
**What goes wrong:** Developer saves results, commits, but the BSON files aren't tracked.
**Why it happens:** `results/` is in `.gitignore` (confirmed: line 13 of `.gitignore`).
**How to avoid:** This is correct behavior -- experiment results should NOT be in git. Document it clearly. The companion `.txt` files also go in `results/` and are also gitignored.
**Warning signs:** None -- this is by design.

## Code Examples

### ExperimentResult Struct Definition
```julia
# Source: designed for this phase based on codebase analysis

"""
    ExperimentResult{C<:AbstractConfig, T<:AbstractFloat}

Complete experiment result with config, trajectory data, and metadata.
Parameterized on config type `C` (preserves KMS/GNS distinction) and
element type `T` (typically Float64).
"""
struct ExperimentResult{C<:AbstractConfig, T<:AbstractFloat}
    # Core results
    config::C
    trajectory_result::TrajectoryResult{Complex{T}}

    # Hamiltonian provenance (not the full HamHam -- just what's needed to reconstruct)
    hamiltonian_params::Dict{Symbol, Any}

    # Metadata
    metadata::Dict{Symbol, Any}
end
```

### Save Function with Metadata Auto-Capture
```julia
using LibGit2, Dates

function save_experiment(result::ExperimentResult, path::String)
    d = _experiment_to_dict(result)
    mkpath(dirname(path))
    BSON.bson(path, d)
    _write_companion_txt(result, replace(path, ".bson" => ".txt"))
    return path
end

function _capture_git_hash()
    try
        project_root = Pkg.project().path |> dirname
        repo = LibGit2.GitRepo(project_root)
        hash = string(LibGit2.head_oid(repo))
        close(repo)
        return hash
    catch
        return "unknown"
    end
end
```

### Load Function with Best-Effort Partial Loading
```julia
function load_experiment(path::String)
    d = BSON.load(path)
    return _dict_to_experiment(d)
end

function _dict_to_experiment(d::Dict)
    # Reconstruct config from Dict fields
    config = _reconstruct_config(d[:config])

    # Reconstruct TrajectoryResult
    traj_d = d[:trajectory]
    traj = TrajectoryResult(
        traj_d[:rho_mean],
        traj_d[:n_trajectories],
        traj_d[:seed],
        get(traj_d, :times, nothing),
        get(traj_d, :measurements_mean, nothing),
    )

    # Metadata -- use get() for forward compatibility (missing fields = nothing)
    metadata = get(d, :metadata, Dict{Symbol,Any}())
    ham_params = get(d, :hamiltonian_params, Dict{Symbol,Any}())

    T = eltype(real(traj.rho_mean))
    C = typeof(config)
    return ExperimentResult{C,T}(config, traj, ham_params, metadata)
end
```

### Companion Text File
```julia
function _write_companion_txt(result::ExperimentResult, path::String)
    open(path, "w") do io
        cfg = result.config
        meta = result.metadata

        println(io, "=== QuantumFurnace Experiment Result ===")
        println(io, "")
        println(io, "Date:       ", get(meta, :timestamp, "unknown"))
        println(io, "Git:        ", get(meta, :git_hash, "unknown"))
        println(io, "Julia:      ", get(meta, :julia_version, "unknown"))
        println(io, "")
        println(io, "--- Config ---")
        println(io, "Type:       ", cfg isa Union{LiouvConfigGNS, ThermalizeConfigGNS} ? "GNS" : "KMS")
        println(io, "Domain:     ", typeof(cfg.domain))
        println(io, "n_qubits:   ", cfg.num_qubits)
        println(io, "beta:       ", cfg.beta)
        println(io, "sigma:      ", cfg.sigma)
        if hasproperty(cfg, :mixing_time)
            println(io, "mix_time:   ", cfg.mixing_time)
            println(io, "delta:      ", cfg.delta)
        end
        println(io, "")
        println(io, "--- Results ---")
        traj = result.trajectory_result
        println(io, "N_traj:     ", traj.n_trajectories)
        println(io, "Seed:       ", traj.seed)
        println(io, "Threads:    ", get(meta, :n_threads, "unknown"))
        println(io, "Wall time:  ", get(meta, :wall_time_seconds, "unknown"), " s")
        println(io, "rho dim:    ", size(traj.rho_mean, 1), "x", size(traj.rho_mean, 2))
    end
end
```

### Filename Generation
```julia
function _generate_experiment_filename(config::AbstractConfig)
    db_str = (config isa Union{LiouvConfigGNS, ThermalizeConfigGNS}) ? "gns" : "kms"
    domain_str = lowercase(replace(string(typeof(config.domain)), "Domain" => ""))
    n_str = "n$(config.num_qubits)"
    beta_str = "beta$(Int(config.beta))"
    date_str = Dates.format(Dates.now(), dateformat"yyyymmdd")
    return "$(db_str)_$(n_str)_$(beta_str)_$(domain_str)_$(date_str).bson"
end
```

### Config Reconstruction from Dict
```julia
function _reconstruct_config(d::Dict)
    domain = _string_to_domain(d[:domain])
    config_type = d[:config_type]  # "KMS" or "GNS"

    # Determine if this is a Liouvillian or Thermalize config
    has_mixing = haskey(d, :mixing_time) && d[:mixing_time] !== nothing

    if config_type == "GNS" && has_mixing
        ThermalizeConfigGNS(; _dict_to_kwargs(d, domain)...)
    elseif config_type == "GNS"
        LiouvConfigGNS(; _dict_to_kwargs(d, domain)...)
    elseif has_mixing
        ThermalizeConfig(; _dict_to_kwargs(d, domain)...)
    else
        LiouvConfig(; _dict_to_kwargs(d, domain)...)
    end
end

const DOMAIN_LOOKUP = Dict(
    "TrotterDomain" => TrotterDomain(),
    "TimeDomain" => TimeDomain(),
    "EnergyDomain" => EnergyDomain(),
    "BohrDomain" => BohrDomain(),
)

_string_to_domain(s::AbstractString) = DOMAIN_LOOKUP[s]
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `BSON.@save` with raw structs | Dict-based `BSON.bson(path, Dict(...))` | Learned from HamHam legacy loader | Avoids struct versioning and type resolution issues |
| `BSON.load` assuming struct types | `BSON.load` returning Dict of primitives | Proven in regression test framework (Phase 5) | Clean round-trip without module namespace issues |

**Deprecated/outdated:**
- The commented-out save code in `simulations/main_liouv.jl` and `main_thermalize.jl` (using `BSON.bson(path, Dict("results" => result))`) saves raw result structs -- this is fragile and will break if struct definitions change. Phase 15 replaces this with proper Dict conversion.

## Discretion Recommendations

### 1. Use Dedicated save_experiment/load_experiment Wrapper Functions
**Recommendation: YES, use wrappers.**
**Rationale:** Raw `BSON.bson`/`BSON.load` works for simple data but doesn't handle: (a) directory creation, (b) companion text file, (c) metadata auto-capture (git hash, timestamp, Julia version), (d) Dict conversion for safe serialization, (e) type reconstruction on load. Wrappers centralize all of this. The existing `_generate_filename` pattern in `misc_tools.jl` already suggests the codebase wants this abstraction.

### 2. Embed TrajectoryResult as a Field (Don't Flatten)
**Recommendation: Embed as a field.**
**Rationale:** TrajectoryResult has 5 fields with clear semantics. Embedding preserves the logical grouping and avoids name collisions (e.g., `seed` appears in both TrajectoryResult and metadata). In the BSON Dict, it becomes a nested sub-Dict under the `:trajectory` key. On load, it's reconstructed directly from the sub-Dict fields.

### 3. Companion .txt Format
**Recommendation: Simple key-value format, human-readable.**
**Rationale:** The txt file's purpose is quick browsing without Julia. A minimal key-value format (see code example above) covers the essential information: date, config type, domain, num_qubits, beta, N_traj, seed, wall time. No need for JSON or structured format -- this is a convenience file, not a data interchange format.

### 4. Git Hash via LibGit2 (Not Shell Command)
**Recommendation: Use LibGit2.**
**Rationale:** LibGit2 is a Julia stdlib with zero external dependencies. The API is simple: `LibGit2.head_oid(LibGit2.GitRepo(path))` returns a `GitHash` that converts to string. This works even if `git` is not on PATH (e.g., on HPC clusters). The try/catch wrapper handles the case where the repo is not a git repo (returns `"unknown"`).

## Open Questions

1. **LindbladianResult / HotSpectralResults embedding**
   - What we know: The CONTEXT.md mentions "reuse existing result types (TrajectoryResult, LindbladianResult, etc.)" but there is no `LindbladianResult` struct in the codebase -- there is `HotSpectralResults`. TrajectoryResult exists and is well-defined.
   - What's unclear: Whether ExperimentResult should also support embedding HotSpectralResults (Liouvillian spectral analysis) or only TrajectoryResult.
   - Recommendation: For this phase, focus on TrajectoryResult only (matching the DATA-01/DATA-02 requirements about trajectory experiments). HotSpectralResults embedding can be added later if needed. The `ExperimentResult{C,T}` struct can be extended with a Union field or a second result type parameter in a future phase.

2. **Hamiltonian parameter storage: "J, h" naming**
   - What we know: The user wants "n, J, h" plus coefficient matrices. The HamHam struct stores `base_terms` (Pauli matrices), `base_coeffs`, `disordering_term`, and `disordering_coeffs`.
   - What's unclear: "J" and "h" are physics labels (coupling and field strengths), but the code uses `base_coeffs` (typically `[1.0, 1.0, 1.0]` for Heisenberg) and `disordering_coeffs` (per-site random values). There's no single "J" and "h" field.
   - Recommendation: Store `base_coeffs`, `disordering_coeffs`, `base_terms` (the Pauli matrix lists), `num_qubits`, and `periodic` in the `hamiltonian_params` Dict. This captures everything needed to reconstruct via `HamHam(terms, coeffs, ...)` without storing the full eigendecomposition. Label them with their code names, not physics names, for consistency.

## Sources

### Primary (HIGH confidence)
- Codebase: `src/structs.jl` -- all config structs, result types, type hierarchy
- Codebase: `src/trajectories.jl` -- TrajectoryResult definition and run_trajectories return type
- Codebase: `src/misc_tools.jl` -- existing `_load_hamiltonian_bson` (BSON.parse + manual reconstruction pattern), `_generate_filename`
- Codebase: `test/reference/generate_references.jl` -- Dict-based BSON save pattern (proven working)
- Codebase: `test/test_regression.jl` -- BSON.load with Dict-based access pattern (proven working)
- Codebase: `test/test_helpers.jl` -- `_load_test_hamiltonian` (manual BSON reconstruction)
- Codebase: `Project.toml` -- BSON 0.3, dependencies list
- Codebase: `Manifest.toml` -- BSON 0.3.9 confirmed
- Julia docs: [LibGit2 stdlib](https://docs.julialang.org/en/v1/stdlib/LibGit2/) -- `head_oid`, `GitRepo` API

### Secondary (MEDIUM confidence)
- [BSON.jl GitHub](https://github.com/JuliaIO/BSON.jl) -- README, namespace resolution with `@__MODULE__`
- [BSON.jl Issue #48](https://github.com/JuliaIO/BSON.jl/issues/48) -- Array type not preserved on load (abstract type parameter issue)
- [BSON.jl Issue #69](https://github.com/JuliaIO/BSON.jl/issues/69) -- World age issue with closures
- [Julia Discourse: LibGit2 commit hash](https://discourse.julialang.org/t/how-to-get-the-commit-of-head-of-a-julia-package/10877) -- Programmatic HEAD commit retrieval

### Tertiary (LOW confidence)
- None -- all findings verified against codebase or official documentation.

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- BSON.jl already in project, all tools are Julia stdlibs
- Architecture: HIGH -- Dict-based serialization pattern proven in existing codebase (regression tests, Hamiltonian loader)
- Pitfalls: HIGH -- All pitfalls derived from actual issues encountered in this codebase (HamHam legacy loader) or verified via BSON.jl GitHub issues

**Research date:** 2026-02-16
**Valid until:** 2026-03-16 (stable domain -- BSON.jl 0.3 is mature and unlikely to change)
