#!/usr/bin/env julia
#
# α(ν, ν) diagonal-of-Kossakowski analysis for fair-comparison hypothesis H2
# (epic qf-mto.5).
#
# Hypothesis H2: the per-Bohr-frequency transition rate α(ν, ν) (the
# Kossakowski diagonal) is what controls the spectral gap; off-diagonal /
# rank structure is secondary. If H2 is correct, CKG smooth-Metropolis and
# DLL Metropolis should have *similar* α_diag(ν) curves at every (n, β),
# because they share the Metropolis acceptance shape, while DLL Gaussian
# should diverge — especially at low T (large β).
#
# What this script does:
#   1. For each (n, β, sampler) on n ∈ {3, 4, 5}, β ∈ {1, 5, 10, 20}:
#      a. Build the disordered Heisenberg fixture, get Bohr frequencies B_H.
#      b. Compute α(ν, ν) for each ν ∈ B_H, summed over the 3n Pauli
#         couplings (X, Y, Z on each site). Since the Kossakowski formula
#         is coupling-independent for both CKG and DLL, this aggregate is
#         (3n) · α^{single}(ν, ν).
#   2. Compute the L²(ν) distance between α^{CKG smooth-Metro} and
#      each DLL variant. Normalise by ‖α^CKG‖₂ so the result is scale-free.
#      Uniform weight on the ν grid (no transition-probability weighting —
#      H2 is a per-frequency claim, so no extra weights).
#   3. Report tables; persist BSON for plotting.
#
# Sanity checks:
#   - Σ_ν α_diag(ν) ≈ dissipator_trace_alpha(α_total)         (consistency)
#   - α_diag(ν) > 0 for all ν                                  (positivity)
#   - DLL Gaussian α_diag drops to ~0 outside |βν| ≲ 2         (collapse signature)
#
# Filter conventions (PHYSICS CHECK):
#   - CKG smooth-Metro: σ = 1/β, a = 0, s = 0.25 — locked thesis-numerics
#     defaults (see MEMORY.md and src/lindblad_action.jl::sweep_mixing_times).
#   - DLL Metropolis: S = 2.0 — flat-top region [-1, 1] strictly contains the
#     fixture's Bohr range [-0.45, 0.45], so the Hörmander bump is invisible.
#   - DLL Gaussian: no S parameter (compact support is implicit).
#
# Run: julia --project scripts/scratch_alpha_diagonal_analysis.jl

using Printf
using LinearAlgebra
using BSON
using QuantumFurnace

# ---------------------------------------------------------------------------
# Sweep parameters
# ---------------------------------------------------------------------------

const NS    = (3, 4, 5)
const BETAS = (1.0, 5.0, 10.0, 20.0)

# CKG smooth-Metropolis parameters (locked thesis-numerics convention).
const CKG_A = 0.0
const CKG_S = 0.25

# DLL Metropolis bump radius (flat top [-S/2, S/2] = [-1, 1] ⊃ B_H).
const DLL_S = 2.0

const SAMPLERS = ("CKG smooth-Metro", "DLL Metropolis", "DLL Gaussian")

# ---------------------------------------------------------------------------
# Hamiltonian + Bohr-frequency setup
# ---------------------------------------------------------------------------

"""
    load_fixture(n, beta) -> (; ham, bohr_freqs, num_jumps)

Load the disordered-Heisenberg n-qubit fixture at given β. Returns the
HamHam (with bohr_dict + Gibbs populated), the sorted unique Bohr-frequency
vector, and the number of single-site Pauli couplings (3n).
"""
function load_fixture(n::Int, beta::Real)
    source_root = dirname(@__DIR__)
    ham_path = joinpath(source_root, "hamiltonians", "heis_xxx_zzdisordered_periodic_n$(n).bson")
    ham = QuantumFurnace._load_hamiltonian_bson(ham_path, Float64(beta))
    bohr_freqs = sort!(collect(keys(ham.bohr_dict)))
    num_jumps = 3 * n  # X, Y, Z on each site
    return (; ham, bohr_freqs, num_jumps)
end

# ---------------------------------------------------------------------------
# Kossakowski-diagonal builders
# ---------------------------------------------------------------------------
#
# Both CKG and DLL Kossakowski formulas are *coupling-independent* in this
# setup (single-site Paulis, no anisotropy), so summing over the 3n
# couplings amounts to multiplying the per-coupling α by 3n.
#
# CKG: α(ν, ν') = create_alpha(ν, ν', β, σ, a, s) with σ = 1/β.
# DLL: α(ν, ν') = freq_kernel(filter, ν) · conj(freq_kernel(filter, ν'))
#      → diagonal: α(ν, ν) = |freq_kernel(filter, ν)|²
#                          = |q(ν)|² · e^{-βν/2}  (since freq_kernel = q · e^{-βν/4})
#
# We return the per-Bohr-frequency "summed-over-couplings" rate
# α_diag(ν) = Σ_a α^a(ν, ν), which is the H2-relevant total per-frequency
# dissipation strength.

"""
    alpha_diag_ckg(beta, bohr_freqs, num_jumps) -> Vector{Float64}

CKG smooth-Metropolis Kossakowski diagonal aggregated over 3n couplings.
"""
function alpha_diag_ckg(beta::Real, bohr_freqs::AbstractVector{<:Real}, num_jumps::Int)
    sigma = 1.0 / Float64(beta)
    out = Vector{Float64}(undef, length(bohr_freqs))
    @inbounds for (i, ν) in enumerate(bohr_freqs)
        out[i] = num_jumps * create_alpha(Float64(ν), Float64(ν), Float64(beta),
                                          sigma, CKG_A, CKG_S)
    end
    return out
end

"""
    alpha_diag_dll(filter, bohr_freqs, num_jumps) -> Vector{Float64}

DLL Kossakowski diagonal aggregated over 3n couplings:
    α_diag(ν) = num_jumps · |freq_kernel(filter, ν)|²
"""
function alpha_diag_dll(filter::AbstractFilter, bohr_freqs::AbstractVector{<:Real},
                        num_jumps::Int)
    out = Vector{Float64}(undef, length(bohr_freqs))
    @inbounds for (i, ν) in enumerate(bohr_freqs)
        f = Float64(freq_kernel(filter, Float64(ν)))
        out[i] = num_jumps * f * f
    end
    return out
end

# ---------------------------------------------------------------------------
# L²(ν) distance — uniform weight, normalised
# ---------------------------------------------------------------------------
#
# L²(ν) distance with uniform weights on the unique Bohr grid:
#   d²(α, α') = Σ_ν |α(ν) - α'(ν)|²
# Normalisation: dist_norm = √(d²) / ‖α‖₂  → scale-invariant.
#
# Choice rationale (PHYSICS CHECK): no transition-probability weighting.
# H2 claims that the per-frequency rate is the driver, so we treat all
# Bohr frequencies equally — a uniformly-distributed mismatch is a
# uniform contradiction of H2. Multiplicities of degenerate Bohr
# transitions are not folded in either; we work on the *unique* Bohr
# grid (consistent with how the rest of the analysis presents B_H).
"""
    l2_distance_normalised(alpha_ref, alpha_other) -> Float64

‖alpha_ref - alpha_other‖₂ / ‖alpha_ref‖₂ on the (uniform) Bohr grid.
Returns Inf if alpha_ref ≡ 0 (numerically guarded by max-relative).
"""
function l2_distance_normalised(alpha_ref::AbstractVector{<:Real},
                                alpha_other::AbstractVector{<:Real})
    @assert length(alpha_ref) == length(alpha_other)
    nref = norm(alpha_ref)
    if nref == 0
        return Inf
    end
    return norm(alpha_ref .- alpha_other) / nref
end

# ---------------------------------------------------------------------------
# Per-cell driver
# ---------------------------------------------------------------------------

"""
    run_cell(n, beta) -> NamedTuple

For a given (n, β), build the three α_diag vectors and the L² distances
(DLL variants vs CKG smooth-Metro reference). Returns a NamedTuple.
"""
function run_cell(n::Int, beta::Float64)
    fix = load_fixture(n, beta)
    bohr_freqs = fix.bohr_freqs
    num_jumps = fix.num_jumps

    α_ckg  = alpha_diag_ckg(beta, bohr_freqs, num_jumps)
    α_dllm = alpha_diag_dll(DLLMetropolisFilter(beta; S=DLL_S), bohr_freqs, num_jumps)
    α_dllg = alpha_diag_dll(DLLGaussianFilter(beta), bohr_freqs, num_jumps)

    # Sanity: positivity
    @assert all(>(-1e-15), α_ckg)  "CKG α_diag has negative entries"
    @assert all(>(-1e-15), α_dllm) "DLL Metro α_diag has negative entries"
    @assert all(>(-1e-15), α_dllg) "DLL Gauss α_diag has negative entries"

    # Sanity: cross-check Σ_ν α_diag(ν) vs dissipator_trace_alpha applied
    # to the per-coupling Kossakowski matrix (DLL only — for CKG we don't
    # have a packaged matrix builder; sum-equals-trace is then trivial).
    α_dll_metro_mat = dll_kossakowski_bohr(DLLMetropolisFilter(beta; S=DLL_S), bohr_freqs)
    α_dll_gauss_mat = dll_kossakowski_bohr(DLLGaussianFilter(beta), bohr_freqs)
    sum_metro_check = num_jumps * dissipator_trace_alpha(α_dll_metro_mat)
    sum_gauss_check = num_jumps * dissipator_trace_alpha(α_dll_gauss_mat)
    rel_metro = abs(sum(α_dllm) - sum_metro_check) / max(abs(sum_metro_check), eps())
    rel_gauss = abs(sum(α_dllg) - sum_gauss_check) / max(abs(sum_gauss_check), eps())
    @assert rel_metro < 1e-12 "DLL Metro Σα_diag mismatch with Tr(α): rel=$(rel_metro)"
    @assert rel_gauss < 1e-12 "DLL Gauss Σα_diag mismatch with Tr(α): rel=$(rel_gauss)"

    dist_metro = l2_distance_normalised(α_ckg, α_dllm)
    dist_gauss = l2_distance_normalised(α_ckg, α_dllg)

    return (;
        n,
        beta,
        bohr_freqs,
        α_ckg,
        α_dll_metro = α_dllm,
        α_dll_gauss = α_dllg,
        sum_α_ckg   = sum(α_ckg),
        sum_α_metro = sum(α_dllm),
        sum_α_gauss = sum(α_dllg),
        dist_metro,
        dist_gauss,
    )
end

# ---------------------------------------------------------------------------
# Run sweep
# ---------------------------------------------------------------------------

println("="^96)
println("α(ν,ν) diagonal-of-Kossakowski analysis — H2 hypothesis test (qf-mto.5)")
println("="^96)
@printf("Sweep: n ∈ %s, β ∈ %s\n", NS, BETAS)
@printf("Samplers: %s\n", SAMPLERS)
@printf("DLL Metropolis bump radius S = %.1f (flat-top [-%.1f, %.1f] ⊃ B_H ⊂ [-0.45, 0.45])\n\n",
        DLL_S, DLL_S/2, DLL_S/2)

cells = NamedTuple[]
for n in NS
    for β in BETAS
        cell = run_cell(n, β)
        push!(cells, cell)
    end
end

# ---------------------------------------------------------------------------
# Tables — L² distance (DLL_x vs CKG smooth-Metro) per (n, β)
# ---------------------------------------------------------------------------

println("="^96)
println("L² distance table: ‖α_DLL_Metro - α_CKG‖₂ / ‖α_CKG‖₂  (smaller = H2 supported)")
println("="^96)
@printf("%-6s | ", "n \\ β")
for β in BETAS
    @printf("%-12s ", @sprintf("β=%.0f", β))
end
println()
println("-"^96)
for n in NS
    @printf("%-6d | ", n)
    for β in BETAS
        cell = cells[findfirst(c -> c.n == n && c.beta == β, cells)]
        @printf("%-12.4f ", cell.dist_metro)
    end
    println()
end

println()
println("="^96)
println("L² distance table: ‖α_DLL_Gauss - α_CKG‖₂ / ‖α_CKG‖₂  (expected larger, esp. high β)")
println("="^96)
@printf("%-6s | ", "n \\ β")
for β in BETAS
    @printf("%-12s ", @sprintf("β=%.0f", β))
end
println()
println("-"^96)
for n in NS
    @printf("%-6d | ", n)
    for β in BETAS
        cell = cells[findfirst(c -> c.n == n && c.beta == β, cells)]
        @printf("%-12.4f ", cell.dist_gauss)
    end
    println()
end

# ---------------------------------------------------------------------------
# Per-frequency inspection at the (n=3, β=20) extreme cell — for sanity
# ---------------------------------------------------------------------------

println("\n", "="^96)
println("Sample α_diag(ν) at n=3, β=20 (max-stress cell): a few representative ν")
println("="^96)
cell_extreme = cells[findfirst(c -> c.n == 3 && c.beta == 20.0, cells)]
@printf("%-12s | %-12s | %-12s | %-12s\n", "ν", "CKG sM", "DLL Metro", "DLL Gauss")
println("-"^60)
# Pick a sparse selection: most-negative, near-zero, most-positive, plus a few intermediate
νs = cell_extreme.bohr_freqs
K = length(νs)
sample_indices = unique(round.(Int, range(1, K; length=11)))
for i in sample_indices
    @printf("%-12.6f | %-12.4e | %-12.4e | %-12.4e\n",
            νs[i],
            cell_extreme.α_ckg[i],
            cell_extreme.α_dll_metro[i],
            cell_extreme.α_dll_gauss[i])
end

# ---------------------------------------------------------------------------
# Headline H2 verdict
# ---------------------------------------------------------------------------

dist_metro_max = maximum(c.dist_metro for c in cells)
dist_metro_avg = sum(c.dist_metro for c in cells) / length(cells)
dist_gauss_max = maximum(c.dist_gauss for c in cells)
dist_gauss_avg = sum(c.dist_gauss for c in cells) / length(cells)

println("\n", "="^96)
println("HEADLINE — H2 verdict")
println("="^96)
@printf("DLL Metropolis vs CKG smooth-Metro:  max dist_norm = %.4f, avg = %.4f\n",
        dist_metro_max, dist_metro_avg)
@printf("DLL Gaussian   vs CKG smooth-Metro:  max dist_norm = %.4f, avg = %.4f\n",
        dist_gauss_max, dist_gauss_avg)

const H2_THRESHOLD = 0.30  # "small" if normalised L² distance below 30%
if dist_metro_max < H2_THRESHOLD
    println("\n→ H2 SUPPORTED: DLL Metro α_diag stays within $(round(Int, 100*H2_THRESHOLD))% of CKG smooth-Metro across all $(length(cells)) cells.")
    println("  The Metropolis-shaped diagonal is the structural commonality.")
else
    println("\n→ H2 NOT cleanly confirmed at the $(round(Int, 100*H2_THRESHOLD))% level.")
    println("  Worst cell: dist_metro = $(round(dist_metro_max; digits=3)).")
    println("  Inspect per-cell breakdown above.")
end

# Cells where the prediction failed (DLL Metro genuinely diverges from CKG)
fail_cells = [(c.n, c.beta, c.dist_metro) for c in cells if c.dist_metro > H2_THRESHOLD]
if !isempty(fail_cells)
    println("\nCells where DLL Metro vs CKG dist > $(H2_THRESHOLD):")
    for (n, β, d) in fail_cells
        @printf("  n=%d, β=%.1f: dist=%.4f\n", n, β, d)
    end
end

# ---------------------------------------------------------------------------
# Persist BSON for plotting
# ---------------------------------------------------------------------------

out_dir = joinpath(@__DIR__, "output")
mkpath(out_dir)
out_path = joinpath(out_dir, "alpha_diagonal_analysis.bson")

# Flatten into per (n, β, sampler) records for the plotting script.
records = NamedTuple[]
for c in cells
    push!(records, (; c.n, beta=c.beta, sampler="CKG smooth-Metro",
                       ν_grid=c.bohr_freqs, α_diag=c.α_ckg, dist_to_ckg=0.0))
    push!(records, (; c.n, beta=c.beta, sampler="DLL Metropolis",
                       ν_grid=c.bohr_freqs, α_diag=c.α_dll_metro,
                       dist_to_ckg=c.dist_metro))
    push!(records, (; c.n, beta=c.beta, sampler="DLL Gaussian",
                       ν_grid=c.bohr_freqs, α_diag=c.α_dll_gauss,
                       dist_to_ckg=c.dist_gauss))
end

BSON.bson(out_path, Dict(
    :records => records,
    :ns      => collect(NS),
    :betas   => collect(BETAS),
    :samplers => collect(SAMPLERS),
    :ckg_a    => CKG_A,
    :ckg_s    => CKG_S,
    :dll_S    => DLL_S,
))
@info "Saved alpha-diagonal analysis BSON" out_path

println("\nDone.")
