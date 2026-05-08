#=
scripts/scratch_fair_comparison_dirichlet_sweep.jl

Headline KMS-Dirichlet fair-comparison sweep for beads issue qf-mto.4.

Runs three samplers
   • CKG smooth-Metropolis (a = 0, s = 0.25; thesis defaults; γ_norm baked in by production)
   • DLL Metropolis (S = 2)
   • DLL Gaussian
across n ∈ {3, 4, 5} and β ∈ {1, 5, 10, 20} — 36 cells total.

Per cell we materialise the dense Schrödinger Lindbladian L_S from the
production matrix-free `apply_lindbladian!` (so γ_norm and Lamb-shift are
both included for CKG; matched DLL-native normalisation for DLL) and
compute, via the freshly added `src/kms_geometry.jl` utilities,

   λ              = `spectral_gap_kms(L_S, σ_β).gap`             — KMS-Poincaré gap
   Λ_max          = `max_dirichlet_rate_kms(L_S, σ_β)`           — largest Dirichlet rate
   ρ_intrinsic    = `intrinsic_mixing_ratio(L_S, σ_β) = λ/Λ_max` — scale-invariant ratio
   tr_alpha       = `dissipator_trace_alpha(α_per_coupling)`     — Σ α(ν,ν), per coupling
   one_to_one     = `dissipator_one_to_one_norm_bound(L_a_list)` — DLL only (CKG: NaN)
   hs_norm        = `hs_operator_norm(L_S)`                      — opnorm(vec(L_S))
   τ_mix^pred     = (1/λ) · log(1 / (p_min · ε)),  ε = 1e-3      — gap-to-mixing bound
   τ_mix^meas     = lookup from drafts/figures/numerics/ckg_vs_dll_taumix.bson (NaN if absent)

The per-coupling Kossakowski α used for `tr_alpha`:
  • DLL:  rank-1 outer product `freq_kernel(filter, ν)·conj(freq_kernel(filter, ν'))`
          via `dll_kossakowski_bohr(filter, bohr_freqs)`. This is the *same*
          α that production code uses inside `dll_lindblad_op_bohr` per coupling.
  • CKG smooth-Metro: K×K matrix `α[k, l] = create_alpha(ν_k, ν_l, β, σ, a, s)`
          via `bohr_domain.jl::create_alpha`, evaluated on the unique Bohr
          frequencies (Eq. 4.6 of Chen et al. 2025 with the smooth-Metropolis
          γ used in production). γ_norm is *not* baked into α — it lives in
          B (`coherent.jl:38`); we report bare α here, consistent with how
          DLL reports bare α.

Sanity assertions per cell:
  • λ > 0  and  Λ_max ≥ λ
  • ρ_intrinsic ∈ (0, 1]
  • DLL Gauss ρ_intrinsic at (n=3, β=20) is the smallest of the three at that cell
    (collapse from prior findings; matches CKG≈DLL Metro headline tie).

Output:
  • `scripts/output/fair_comparison_kms_dirichlet.bson`
      Vector of NamedTuple, one per cell with all of the above fields.
  • Stdout: 3 sampler tables (n, β, λ, Λ_max, ρ_intrinsic, tr_alpha,
    one_to_one_bound, hs_norm, τ_mix^pred, τ_mix^meas).

Run:   julia --project scripts/scratch_fair_comparison_dirichlet_sweep.jl

Cost estimate: n=5 dense L_S is 1024×1024; eigendecomp of 1024×1024 ~5s,
build_dense_superoperator (apply_lindbladian! × 1024 columns) ~ a few minutes
per (n, β, sampler) cell at n=5. Total wall time ≈ 30 min on a laptop.
=#

using QuantumFurnace
using LinearAlgebra
using Printf
using BSON
using Dates

# ============================================================================
# Configuration & paths
# ============================================================================

const N_VALUES         = (3, 4, 5)
const BETA_VALUES      = (1.0, 5.0, 10.0, 20.0)
# PHYSICS CHECK: smooth-Metropolis defaults a=0, s=0.25 are the locked thesis
# convention (see MEMORY.md and src/lindblad_action.jl::sweep_mixing_times).
const CKG_A            = 0.0
const CKG_S            = 0.25
# PHYSICS CHECK: DLL Metropolis cutoff S=2.0 is the value used throughout
# the existing DLL Metropolis tests / `make_dll_n3_system` family.
const DLL_METRO_S      = 2.0
# PHYSICS CHECK: ε = 1e-3 matches the target_epsilon used in
# drafts/figures/numerics/ckg_vs_dll_taumix.bson, so τ_mix^pred and τ_mix^meas
# share the same precision target. log(1/(p_min·ε)) = log(1/(p_min·1e-3)).
const EPSILON          = 1e-3
const NUM_ENERGY_BITS  = 12
const W0               = 0.05
const T0               = 2π / (2^NUM_ENERGY_BITS * W0)
const NUM_TROTTER_STEPS_PER_T0 = 10

const PROJECT_ROOT     = dirname(@__DIR__)
const OUTPUT_DIR       = joinpath(PROJECT_ROOT, "scripts", "output")
const OUTPUT_BSON      = joinpath(OUTPUT_DIR, "fair_comparison_kms_dirichlet.bson")
const TAUMIX_BSON      = joinpath(PROJECT_ROOT, "drafts", "figures", "numerics", "ckg_vs_dll_taumix.bson")


# ============================================================================
# System builders (n-parametric variants of the helpers in
# scripts/scratch_kms_geometry.jl, used here to make this script self-contained
# and to extend the n=3 helpers to n ∈ {3,4,5})
# ============================================================================

"""
    _build_jump_set(ham::HamHam, n::Integer) -> Vector{JumpOp}

Standard 3n single-site Pauli jump set (X/Y/Z on each site), normalised by
sqrt(3n). Mirrors `make_dll_n3_system` and `_build_jump_set` in
`src/lindblad_action.jl`. Construction-agnostic — same jumps used for CKG
and DLL.
"""
function _build_jump_set(ham::HamHam, n::Integer)
    paulis = ([X], [Y], [Z])
    num_jumps = length(paulis) * n
    jump_norm = sqrt(num_jumps)
    jumps = JumpOp[]
    for pauli in paulis
        for site in 1:n
            op = Matrix(pad_term(pauli, n, site)) ./ jump_norm
            op_eb = ham.eigvecs' * op * ham.eigvecs
            push!(jumps, JumpOp(op, op_eb, op == transpose(op), op == op'))
        end
    end
    return jumps
end

"""
    ckg_smooth_metro_config(n, β) -> Config{Lindbladian, EnergyDomain, KMS}

CKG smooth-Metropolis EnergyDomain config at (n, β) with locked thesis
defaults (a=0, s=0.25). EnergyDomain is the production fast path for CKG
(see qf-lkb.11). Production γ_norm is applied automatically by
`_precompute_data + apply_lindbladian!`.
"""
function ckg_smooth_metro_config(n::Integer, β::Real)
    return Config(;
        sim = Lindbladian(),
        domain = EnergyDomain(),
        construction = KMS(),
        num_qubits = n,
        with_linear_combination = true,
        beta = float(β),
        sigma = 1.0 / float(β),
        a = CKG_A,
        s = CKG_S,
        num_energy_bits = NUM_ENERGY_BITS,
        w0 = W0,
        t0 = T0,
        num_trotter_steps_per_t0 = NUM_TROTTER_STEPS_PER_T0,
    )
end

"""
    dll_config(n, β, filter) -> Config{Lindbladian, BohrDomain, DLL}

DLL BohrDomain config at (n, β) with the supplied filter (Metropolis or
Gaussian). BohrDomain is the production fast path for DLL (see qf-hur).
"""
function dll_config(n::Integer, β::Real, filter::AbstractFilter)
    return Config(;
        sim = Lindbladian(),
        domain = BohrDomain(),
        construction = DLL(),
        num_qubits = n,
        with_linear_combination = true,
        beta = float(β),
        sigma = 1.0 / float(β),
        a = float(β) / 30.0,
        s = 0.4,
        num_energy_bits = NUM_ENERGY_BITS,
        w0 = W0,
        t0 = T0,
        num_trotter_steps_per_t0 = NUM_TROTTER_STEPS_PER_T0,
        filter = filter,
    )
end

"""
    build_dense_lindbladian(config, ham, jumps) -> (L_super::Matrix{ComplexF64}, dim::Int)

Materialise the Schrödinger Lindbladian as a dense d²×d² matrix from the
production matrix-free `apply_lindbladian!` (γ_norm + Lamb-shift included).
"""
function build_dense_lindbladian(config::Config{Lindbladian}, ham::HamHam,
                                  jumps::Vector{JumpOp})
    dim = size(ham.data, 1)
    ws  = Workspace(config, ham, jumps)
    L_apply! = function(out, X)
        apply_lindbladian!(ws, X, config, ham)
        copyto!(out, ws.scratch.rho_out)
        return out
    end
    L_super = build_dense_superoperator(L_apply!, dim)
    return L_super, dim
end


# ============================================================================
# Per-coupling Kossakowski matrices
# ============================================================================

"""
    ckg_kossakowski_smooth_metro(β, σ, a, s, bohr_freqs) -> Matrix{ComplexF64}

CKG smooth-Metropolis Kossakowski α^{CKG}_{ν, ν'} = `create_alpha(ν, ν', β,
σ, a, s)` (Eq. 4.6 of Chen et al. 2025). Returns a K×K complex matrix on
the unique-Bohr-frequency grid. Reported "bare" — γ_norm is *not* applied
here because (i) DLL has no analogous normalisation and (ii) γ_norm in
production lives in B, not in α.
"""
function ckg_kossakowski_smooth_metro(β::Real, σ::Real, a::Real, s::Real,
                                       bohr_freqs::AbstractVector{<:Real})
    K = length(bohr_freqs)
    α = zeros(ComplexF64, K, K)
    @inbounds for q in 1:K, p in 1:K
        α[p, q] = create_alpha(bohr_freqs[p], bohr_freqs[q], β, σ, a, s)
    end
    return α
end

"""
    sampler_kossakowski(name, ham, β, filter_or_nothing) -> Matrix{ComplexF64}

Per-coupling Kossakowski α for a given sampler. CKG uses the smooth-Metro
analytic form on the unique Bohr-frequency grid; DLL uses the production
rank-1 outer product `dll_kossakowski_bohr(filter, ham)`.
"""
function sampler_kossakowski(sampler_name::Symbol, ham::HamHam, β::Real,
                              filter::Union{Nothing, AbstractFilter})
    if sampler_name === :ckg_smooth_metro
        bohr_freqs = sort!(collect(keys(ham.bohr_dict)))
        σ = 1.0 / float(β)
        return ckg_kossakowski_smooth_metro(β, σ, CKG_A, CKG_S, bohr_freqs)
    elseif sampler_name in (:dll_metro, :dll_gauss)
        @assert filter !== nothing  "DLL sampler requires a filter"
        α, _ = dll_kossakowski_bohr(filter, ham)
        return Matrix{ComplexF64}(α)
    else
        throw(ArgumentError("Unknown sampler: $sampler_name"))
    end
end


# ============================================================================
# DLL Lindblad-operator list (for the 1→1 norm bound)
# ============================================================================

"""
    dll_lindblad_ops(jumps, ham, filter) -> Vector{Matrix{ComplexF64}}

List of per-coupling DLL Lindblad operators in the eigenbasis (one per
jump), used by `dissipator_one_to_one_norm_bound`.
"""
function dll_lindblad_ops(jumps::Vector{JumpOp}, ham::HamHam, filter::AbstractFilter)
    return [Matrix{ComplexF64}(dll_lindblad_op_bohr(jump, ham, filter)) for jump in jumps]
end


# ============================================================================
# τ_mix^measured lookup from existing taumix BSON (β=1, 2, 5, 10 grid)
# ============================================================================

"""
    load_measured_taumix(path) -> Dict{Tuple{Int,Float64,Symbol}, Float64}

Read the existing CKG-vs-DLL τ_mix sweep BSON and return a flat
`(n, β, sampler) → τ_mix` dictionary. Sampler keys are
`:ckg_smooth_metro`, `:dll_metro`, `:dll_gauss` to match this script's
naming. β values not present (e.g. β=20) silently absent — caller falls
back to NaN.
"""
function load_measured_taumix(path::AbstractString)
    out = Dict{Tuple{Int, Float64, Symbol}, Float64}()
    isfile(path) || return out
    data = BSON.load(path)
    pairs = (
        (:results_ckg,        :ckg_smooth_metro),
        (:results_dll_metro,  :dll_metro),
        (:results_dll_gauss,  :dll_gauss),
    )
    for (key, sampler) in pairs
        haskey(data, key) || continue
        for r in data[key]
            out[(Int(r.n), Float64(r.beta), sampler)] = Float64(r.mixing_time)
        end
    end
    return out
end


# ============================================================================
# Per-cell driver
# ============================================================================

"""
    run_cell(n, β, sampler_name, ham, taumix_lut) -> NamedTuple

Build the dense Lindbladian for one (n, β, sampler) cell, run all
KMS-geometry diagnostics, and return a NamedTuple ready to push to the
result vector.
"""
function run_cell(n::Integer, β::Real, sampler_name::Symbol, ham::HamHam,
                  taumix_lut::Dict{Tuple{Int, Float64, Symbol}, Float64})
    jumps = _build_jump_set(ham, n)

    # Sampler-specific config + (for DLL) filter object.
    if sampler_name === :ckg_smooth_metro
        cfg = ckg_smooth_metro_config(n, β)
        filter_obj = nothing
        sampler_label = "CKG smooth-Metro"
    elseif sampler_name === :dll_metro
        filter_obj = DLLMetropolisFilter(float(β); S = DLL_METRO_S)
        cfg = dll_config(n, β, filter_obj)
        sampler_label = "DLL Metro"
    elseif sampler_name === :dll_gauss
        filter_obj = DLLGaussianFilter(float(β))
        cfg = dll_config(n, β, filter_obj)
        sampler_label = "DLL Gauss"
    else
        throw(ArgumentError("Unknown sampler: $sampler_name"))
    end

    # Dense Lindbladian + Gibbs state (already diagonal in eigenbasis).
    L_super, dim = build_dense_lindbladian(cfg, ham, jumps)
    gibbs        = Matrix{ComplexF64}(ham.gibbs)

    # KMS-geometry diagnostics.
    gap_info = spectral_gap_kms(L_super, gibbs)
    λ        = gap_info.gap
    Λ_max    = max_dirichlet_rate_kms(L_super, gibbs)
    ρ_intr   = intrinsic_mixing_ratio(L_super, gibbs)
    hs_norm  = hs_operator_norm(L_super)

    # Per-coupling Kossakowski → Tr(α).
    α_matrix = sampler_kossakowski(sampler_name, ham, β, filter_obj)
    tr_alpha = dissipator_trace_alpha(α_matrix)

    # 1→1 bound: DLL only (extracting CKG L_a's requires α eigendecomp; deferred).
    one_to_one = if filter_obj === nothing
        NaN
    else
        dissipator_one_to_one_norm_bound(dll_lindblad_ops(jumps, ham, filter_obj))
    end

    # τ_mix^predicted = (1/λ) · log(1 / (p_min · ε)) with ε = 1e-3.
    # PHYSICS CHECK: this is the standard gap-to-mixing upper bound; p_min is
    # the smallest Gibbs eigenvalue (≈ exp(-β·E_max)/Z), hardest-to-reach mass.
    p_min       = real(minimum(eigvals(Hermitian(gibbs))))
    p_min_safe  = max(p_min, eps(Float64))    # guard against numerical zero
    τ_mix_pred  = (1 / λ) * log(1 / (p_min_safe * EPSILON))

    # τ_mix^measured from existing sweep, NaN if absent.
    τ_mix_meas = get(taumix_lut, (Int(n), Float64(β), sampler_name), NaN)

    # ----- Sanity assertions (cell-local) -----
    @assert λ > 0      "λ must be > 0 (KMS-DBC ⇒ nondegenerate gap); got λ = $λ"
    @assert Λ_max ≥ λ - 1e-12  "Λ_max < λ at (n=$n, β=$β, $sampler_label)"
    @assert 0 < ρ_intr ≤ 1 + 1e-12  "ρ_intrinsic out of (0, 1] at (n=$n, β=$β, $sampler_label): ρ=$ρ_intr"

    return (;
        n = Int(n),
        β = Float64(β),
        sampler_name = String(sampler_label),
        sampler_key  = sampler_name,
        λ            = λ,
        Λ_max        = Λ_max,
        ρ_intrinsic  = ρ_intr,
        tr_alpha     = tr_alpha,
        one_to_one_bound = one_to_one,
        hs_norm      = hs_norm,
        tau_mix_predicted = τ_mix_pred,
        tau_mix_measured  = τ_mix_meas,
        p_min        = p_min,
    )
end


# ============================================================================
# Top-level sweep + reporting
# ============================================================================

function print_sampler_table(label::AbstractString, rows::Vector{<:NamedTuple})
    println("\n", "="^140)
    println("Sampler: ", label)
    println("="^140)
    @printf("%-3s  %-6s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s\n",
        "n", "β", "λ", "Λ_max", "ρ_intr", "Tr(α)", "1→1 bound", "‖L‖_HS",
        "τ_mix^pred", "τ_mix^meas")
    println("-"^140)
    for r in rows
        one_str  = isnan(r.one_to_one_bound) ? "      n/a   " : @sprintf("%12.4e", r.one_to_one_bound)
        meas_str = isnan(r.tau_mix_measured) ? "      n/a   " : @sprintf("%12.4f", r.tau_mix_measured)
        @printf("%-3d  %-6.1f  %12.4e  %12.4e  %12.4e  %12.4e  %s  %12.4e  %12.4f  %s\n",
            r.n, r.β, r.λ, r.Λ_max, r.ρ_intrinsic, r.tr_alpha,
            one_str, r.hs_norm, r.tau_mix_predicted, meas_str)
    end
    println("="^140)
end

function print_headline_comparison(rows::Vector{<:NamedTuple})
    println("\n", "="^110)
    println("HEADLINE: ρ_intrinsic head-to-head — CKG smooth-Metro vs DLL Metro vs DLL Gauss")
    println("="^110)
    @printf("%-3s  %-6s  %16s  %16s  %16s  %14s  %14s\n",
        "n", "β", "ρ(CKG sM)", "ρ(DLL Metro)", "ρ(DLL Gauss)",
        "ρ(CKG)/ρ(DLL_M)", "ρ(DLL_M)/ρ(DLL_G)")
    println("-"^110)
    for n in N_VALUES, β in BETA_VALUES
        # Find the three rows for this (n, β) cell.
        cell_rows = filter(r -> r.n == n && r.β == β, rows)
        ck = first(filter(r -> r.sampler_key === :ckg_smooth_metro, cell_rows)).ρ_intrinsic
        dm = first(filter(r -> r.sampler_key === :dll_metro,        cell_rows)).ρ_intrinsic
        dg = first(filter(r -> r.sampler_key === :dll_gauss,        cell_rows)).ρ_intrinsic
        ratio_cm_dm = ck / dm
        ratio_dm_dg = dm / dg
        @printf("%-3d  %-6.1f  %16.6e  %16.6e  %16.6e  %14.3f  %14.3f\n",
            n, β, ck, dm, dg, ratio_cm_dm, ratio_dm_dg)
    end
    println("="^110)
end

function run_sweep()
    println("="^72)
    println("KMS-Dirichlet fair-comparison sweep (qf-mto.4)")
    println("="^72)
    println("Started:        ", Dates.now())
    println("n ∈ ",  N_VALUES)
    println("β ∈ ",  BETA_VALUES)
    println("samplers: CKG smooth-Metro (a=$CKG_A, s=$CKG_S), DLL Metro (S=$DLL_METRO_S), DLL Gauss")
    println("ε (mixing target) = ", EPSILON)
    println()

    # Pre-load measured τ_mix once (same lookup used for every cell).
    taumix_lut = load_measured_taumix(TAUMIX_BSON)
    if isempty(taumix_lut)
        println("[note] No measured-τ_mix BSON found at $TAUMIX_BSON — τ_mix^meas will be NaN everywhere.")
    else
        println("[note] Measured τ_mix lookup loaded: ", length(taumix_lut), " entries from $TAUMIX_BSON")
    end
    println()

    samplers = (:ckg_smooth_metro, :dll_metro, :dll_gauss)
    results = NamedTuple[]

    for n in N_VALUES
        # Load the n-qubit Hamiltonian once per n; reuse across β and samplers
        # by passing β to load_hamiltonian — but Gibbs state depends on β, so
        # we actually need to reload per (n, β). Inexpensive (n ≤ 5).
        for β in BETA_VALUES
            ham = load_hamiltonian("heis", n; beta = float(β))
            for sampler in samplers
                t_start = time_ns()
                row = run_cell(n, β, sampler, ham, taumix_lut)
                elapsed = (time_ns() - t_start) / 1e9
                push!(results, row)
                @printf("  done: n=%d β=%-5.1f %-22s  λ=%.4e  Λ_max=%.4e  ρ=%.4e  (%.1fs)\n",
                    row.n, row.β, row.sampler_name, row.λ, row.Λ_max, row.ρ_intrinsic, elapsed)
            end
        end
    end

    # ----- Per-sampler tables -----
    rows_ckg   = [r for r in results if r.sampler_key === :ckg_smooth_metro]
    rows_metro = [r for r in results if r.sampler_key === :dll_metro]
    rows_gauss = [r for r in results if r.sampler_key === :dll_gauss]
    print_sampler_table("CKG smooth-Metro (a=$CKG_A, s=$CKG_S)", rows_ckg)
    print_sampler_table("DLL Metropolis (S=$DLL_METRO_S)",       rows_metro)
    print_sampler_table("DLL Gaussian",                          rows_gauss)
    print_headline_comparison(results)

    # ----- Cross-cell sanity assertions -----
    # (a) For each sampler at fixed n, λ should generally decrease as β increases.
    #     Use a soft check (warn only) since for some samplers β=1 may be lower
    #     than β=5 due to weak-coupling regime peculiarities.
    println("\n[soft sanity] λ vs β (expected to weakly decrease with β at fixed n):")
    for sampler in samplers
        for n in N_VALUES
            cell_rows = sort(filter(r -> r.n == n && r.sampler_key === sampler, results); by = r -> r.β)
            λs = [r.λ for r in cell_rows]
            βs = [r.β for r in cell_rows]
            decreasing = all(λs[i] >= λs[i+1] - 1e-3 * λs[i] for i in 1:length(λs)-1)
            tag = decreasing ? "OK " : "INV"
            @printf("  [%s] %-22s n=%d  λ(β=%s) = %s\n",
                tag, String(first(cell_rows).sampler_name), n,
                join((@sprintf("%.1f", β) for β in βs), ", "),
                join((@sprintf("%.3e", λ) for λ in λs), ", "))
        end
    end

    # (b) Headline collapse: at (n=3, β=20), ρ_intrinsic for DLL Gauss should
    #     be the smallest of the three.
    cell_n3_β20 = filter(r -> r.n == 3 && r.β == 20.0, results)
    if length(cell_n3_β20) == 3
        ρ_ck = first(filter(r -> r.sampler_key === :ckg_smooth_metro, cell_n3_β20)).ρ_intrinsic
        ρ_dm = first(filter(r -> r.sampler_key === :dll_metro,        cell_n3_β20)).ρ_intrinsic
        ρ_dg = first(filter(r -> r.sampler_key === :dll_gauss,        cell_n3_β20)).ρ_intrinsic
        @printf("\n[headline] (n=3, β=20):  ρ(CKG sM) = %.4e   ρ(DLL Metro) = %.4e   ρ(DLL Gauss) = %.4e\n",
                ρ_ck, ρ_dm, ρ_dg)
        if ρ_dg ≤ ρ_ck && ρ_dg ≤ ρ_dm
            println("[headline] DLL Gauss has the smallest ρ_intrinsic at (n=3, β=20) — confirms collapse from prior findings.")
        else
            println("[headline] ! WARNING: DLL Gauss ρ_intrinsic NOT smallest at (n=3, β=20). Expected (ρ_DG ≤ ρ_CK and ρ_DG ≤ ρ_DM).")
        end
    end

    # ----- Persist results -----
    isdir(OUTPUT_DIR) || mkpath(OUTPUT_DIR)
    BSON.bson(OUTPUT_BSON, Dict(
        :results        => results,
        :n_values       => collect(N_VALUES),
        :beta_values    => collect(BETA_VALUES),
        :samplers       => collect(samplers),
        :ckg_a          => CKG_A,
        :ckg_s          => CKG_S,
        :dll_metro_S    => DLL_METRO_S,
        :epsilon        => EPSILON,
        :timestamp      => string(Dates.now()),
        :script_path    => @__FILE__,
    ))
    println("\nWrote $(length(results)) cells → $OUTPUT_BSON")
    println("Finished:       ", Dates.now())
    return results
end


# ============================================================================
# Main
# ============================================================================

if abspath(PROGRAM_FILE) == @__FILE__
    run_sweep()
end
