#!/usr/bin/env julia
# migrate_bson_beta_phys.jl  (qf-6vr.7 / Task 7)
#
# Legacy BSON archaeology utility: annotate pre-qf-6vr sidecars with the
# β_phys / β_alg / rescaling_factor triple so downstream readers
# (`fit_scaling(::Vector{<:NamedTuple})`, plot scripts) can consume them
# under the qf-6vr contract without touching the original BSON.
#
# Operation:
#   1. Walk `scripts/output/sweep_S*/` and `drafts/figures/numerics/sweep_cache/`
#      (configurable below). For each `*.bson` whose name carries `_beta<β_alg>`
#      (the legacy filename grammar), parse it.
#   2. Load the `:result` Dict. If `:beta_phys` is already present, skip
#      (already migrated). Otherwise, read `:n` and `:beta` (= β_alg) plus
#      `:family` (if recorded; defaults to `:disordered`).
#   3. Look up the corresponding Hamiltonian fixture's `rescaling_factor`
#      under `hamiltonians/heis_<family>_periodic_n<n>.bson`. Compute
#      `β_phys = β_alg / rescaling_factor`.
#   4. Write a NEW sidecar adjacent to the original with `.betaphys.bson`
#      appended to the basename (or, for channel sweeps, `_betaphys.bson`
#      to fit the existing dot-conventions). The original is preserved.
#   5. Skip cells where the fixture file is missing — emit a `@warn`.
#
# Idempotent: re-running this script over a directory that's already been
# migrated is a no-op (skips files with `:beta_phys` already present).
#
# Usage:
#   julia --project scripts/migrate_bson_beta_phys.jl
#   julia --project scripts/migrate_bson_beta_phys.jl --dry-run
#   julia --project scripts/migrate_bson_beta_phys.jl --dir scripts/output/sweep_S3_channel

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using QuantumFurnace
using BSON
using Printf

# --- Configurable defaults -----------------------------------------------

const DEFAULT_DIRS = String[
    joinpath(@__DIR__, "output"),
    joinpath(@__DIR__, "..", "drafts", "figures", "numerics", "sweep_cache"),
]

const HAM_DIR = joinpath(@__DIR__, "..", "hamiltonians")

# Map a record's `:family` (Symbol or String) to the on-disk fixture pattern.
# When the record carries no `:family`, fall back to the current 1D family
# `heis_xxx_zzdisordered_periodic_n*` (qf-2kd: find_typical + [[Z],[Z,Z]]).
function _ham_path_for(record::AbstractDict)
    n_raw = get(record, :n, nothing)
    n_raw === nothing && return (nothing, "missing :n field")
    n = Int(n_raw)
    family_raw = get(record, :family, nothing)
    family = if family_raw === nothing
        "xxx_zzdisordered"
    elseif family_raw isa Symbol
        String(family_raw)
    else
        String(family_raw)
    end
    fname = "heis_$(family)_periodic_n$(n).bson"
    path = joinpath(HAM_DIR, fname)
    return (path, "n=$n, family=$family")
end

# Pull the legacy β_alg from a result Dict, handling both common key shapes.
function _legacy_beta(record::AbstractDict)
    for key in (:beta_alg, :beta)
        v = get(record, key, nothing)
        v !== nothing && v isa Real && isfinite(v) && return (Float64(v), key)
    end
    return (nothing, nothing)
end

# --- Argument parsing ----------------------------------------------------

function _parse_args()
    dry_run = false
    target_dirs = String[]
    args = copy(ARGS)
    while !isempty(args)
        a = popfirst!(args)
        if a == "--dry-run"
            dry_run = true
        elseif a == "--dir"
            push!(target_dirs, popfirst!(args))
        elseif startswith(a, "--dir=")
            push!(target_dirs, a[length("--dir=")+1:end])
        else
            error("Unknown argument: $a (supported: --dry-run, --dir <path>)")
        end
    end
    if isempty(target_dirs)
        append!(target_dirs, DEFAULT_DIRS)
    end
    return (; dry_run = dry_run, target_dirs = target_dirs)
end

# --- Migration core ------------------------------------------------------

mutable struct MigrationCounters
    scanned::Int
    skipped_no_beta::Int
    skipped_already_migrated::Int
    skipped_no_fixture::Int
    skipped_already_emitted::Int
    written::Int
    errored::Int
end
MigrationCounters() = MigrationCounters(0, 0, 0, 0, 0, 0, 0)

function _migrated_sidecar_path(original::String)
    if endswith(original, ".bson")
        return original[1:end-5] * ".betaphys.bson"
    else
        return original * ".betaphys"
    end
end

function _migrate_one!(c::MigrationCounters, path::String; dry_run::Bool)
    c.scanned += 1
    local loaded
    try
        loaded = BSON.load(path, QuantumFurnace)
    catch err
        @warn "Failed to load BSON; skipping" path err
        c.errored += 1
        return
    end
    if !haskey(loaded, :result) || !(loaded[:result] isa AbstractDict)
        @warn "BSON missing :result Dict; skipping" path
        c.errored += 1
        return
    end
    record = loaded[:result]
    if haskey(record, :beta_phys)
        c.skipped_already_migrated += 1
        return
    end
    β_alg, _ = _legacy_beta(record)
    if β_alg === nothing
        c.skipped_no_beta += 1
        return
    end
    ham_path, ham_descr = _ham_path_for(record)
    if ham_path === nothing || !isfile(ham_path)
        @warn "Hamiltonian fixture missing; cannot derive β_phys" path ham_descr ham_path
        c.skipped_no_fixture += 1
        return
    end

    out_path = _migrated_sidecar_path(path)
    if isfile(out_path)
        c.skipped_already_emitted += 1
        return
    end

    raw = try
        QuantumFurnace._parse_hamiltonian_bson(ham_path)
    catch err
        @warn "Failed to parse Hamiltonian fixture; cannot derive β_phys" ham_path err
        c.errored += 1
        return
    end
    rescale = Float64(raw.rescaling_factor)
    β_phys  = β_alg / rescale

    new_record = Dict{Symbol, Any}(record)  # copy
    new_record[:beta_alg]         = β_alg
    new_record[:beta_phys]        = β_phys
    new_record[:rescaling_factor] = rescale

    if dry_run
        @info "[dry-run] would write annotated BSON" path out_path β_alg β_phys rescale
    else
        try
            BSON.bson(out_path, Dict(:result => new_record))
            c.written += 1
        catch err
            @warn "Failed to write annotated BSON" out_path err
            c.errored += 1
            return
        end
    end
    return
end

function _walk!(c::MigrationCounters, root::String; dry_run::Bool)
    isdir(root) || return
    for (dirpath, _, files) in walkdir(root)
        for f in files
            endswith(f, ".bson") || continue
            endswith(f, ".betaphys.bson") && continue  # don't migrate our own output
            _migrate_one!(c, joinpath(dirpath, f); dry_run = dry_run)
        end
    end
end

# --- Main ----------------------------------------------------------------

function main()
    cfg = _parse_args()
    println("=== qf-6vr.7 BSON archaeology migration ===")
    @printf("dry-run         : %s\n", cfg.dry_run)
    @printf("target dirs     :\n")
    for d in cfg.target_dirs
        @printf("  • %s\n", d)
    end

    counters = MigrationCounters()
    for dir in cfg.target_dirs
        _walk!(counters, dir; dry_run = cfg.dry_run)
    end

    println("\n=== Summary ===")
    @printf("Scanned                 : %d\n", counters.scanned)
    @printf("Written (annotated)     : %d\n", counters.written)
    @printf("Skipped (no β)          : %d\n", counters.skipped_no_beta)
    @printf("Skipped (already β_phys): %d\n", counters.skipped_already_migrated)
    @printf("Skipped (fixture miss)  : %d\n", counters.skipped_no_fixture)
    @printf("Skipped (output exists) : %d\n", counters.skipped_already_emitted)
    @printf("Errored                 : %d\n", counters.errored)
end

main()
