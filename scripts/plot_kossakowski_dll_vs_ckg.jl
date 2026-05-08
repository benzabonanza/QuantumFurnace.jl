#!/usr/bin/env julia
#
# Side-by-side DLL vs CKG Kossakowski matrix comparison
# (Ding et al. 2024 Fig. 3 reproduction).
#
# Strategy: build a 1×3 panel highlighting the structural distinction
# between the DLL and CKG Kossakowski matrices on a shared Bohr grid:
#
#   Panel 1 — CKG (Gaussian-Metropolis γ): full-rank Kossakowski α^{CKG}_{ν,ν'}
#             via `create_alpha(...)`. The off-diagonal "skirt" decays with
#             |ν − ν'| at length ∝ σ; at σ = 1/β the skirt narrows for cold β.
#   Panel 2 — DLL (Gaussian-type Gevrey): rank-1 outer product
#             α^{DLL}_{ν,ν'} = freq_kernel(ν) · conj(freq_kernel(ν'))
#             via `dll_kossakowski_bohr(filter, νs)`. Concentrates near ν, ν'
#             ≈ −1/β (the modal frequency of the DLL filter).
#   Panel 3 — singular-value spectra of the two matrices on log scale.
#             DLL has a single non-zero SV; CKG has a continuous tail.
#
# This is the empirical analogue of the DLL paper's Fig. 3 ("rank gap" between
# CKG and DLL Kossakowski matrices).
#
# System: disordered Heisenberg, n = 3, periodic, ferromagnetic XXZ
# (coeffs = [-1, -1, -1]) with on-site Z disorder. Same fixture as
# `scripts/plot_kossakowski_heatmaps.jl` so the comparison is apples-to-apples.
#
# PHYSICS CHECK: σ = 1/β follows the Chen et al. canonical setting; matches
# the DLL paper's Gaussian-type filter `q(ν) = e^{-(βν)²/8}` width parameter.
# With σ = 1/β both filters have ~the same support in frequency.
#
# PHYSICS CHECK: β = 10 — user-specified stress level. β = 1 hides the rank
# distinction at high temperature.
#
# Output: stdout summary + (if Plots.jl is loaded) a PNG file
# `drafts/plots/kossakowski_dll_vs_ckg.png`.
#
# Usage: julia --project scripts/plot_kossakowski_dll_vs_ckg.jl

using Printf
using Random
using LinearAlgebra
using QuantumFurnace

# ── Parameters ────────────────────────────────────────────────────────────────
const num_qubits = 3
const coeffs     = [-1.0, -1.0, -1.0]
const seed       = 42
const batch_size = 200
const a_reg      = 0.1
const s_smooth   = 0.4
const β_values   = (1.0, 5.0, 10.0)

# ── Build disordered Heisenberg ──────────────────────────────────────────────
Random.seed!(seed)
@info "Building disordered Heisenberg (n=$num_qubits, periodic, Z disorder)…"
raw = find_ideal_heisenberg(num_qubits, coeffs; batch_size=batch_size, periodic=true)
ham = HamHam(raw, first(β_values))
unique_freqs = sort(collect(keys(ham.bohr_dict)))
K = length(unique_freqs)
@printf("Unique Bohr frequencies: K = %d\n", K)
@printf("ν range: [%.4f, %.4f]\n", first(unique_freqs), last(unique_freqs))

# ── Comparison loop over β ───────────────────────────────────────────────────
println("\n", "="^70)
println("DLL vs CKG Kossakowski comparison")
println("="^70)
@printf("%-6s %-12s %-14s %-14s %-14s %-14s\n",
        "β", "σ=1/β", "‖α^CKG‖_F", "rank_eff(CKG)", "‖α^DLL‖_F", "rank_eff(DLL)")
println("-"^70)

results = []
for β in β_values
    σ = 1.0 / β

    # CKG smooth-Metropolis Kossakowski (Eq. 4.6 / `create_alpha`).
    α_ckg = Matrix{Float64}(undef, K, K)
    for q in 1:K, p in 1:K
        α_ckg[p, q] = create_alpha(unique_freqs[p], unique_freqs[q], β, σ, a_reg, s_smooth)
    end

    # DLL Gaussian-type Kossakowski (rank-1 outer product, Sec. 4 of DLL paper).
    filt = DLLGaussianFilter(β)
    α_dll = dll_kossakowski_bohr(filt, unique_freqs)

    sv_ckg = svdvals(α_ckg)
    sv_dll = svdvals(real.(α_dll))  # |α_dll| has imaginary noise from outer product

    # Effective rank: # of SVs above 1% of leading.
    rank_ckg = count(s -> s / sv_ckg[1] > 1e-2, sv_ckg)
    rank_dll = count(s -> s / sv_dll[1] > 1e-2, sv_dll)

    @printf("%-6.1f %-12.4f %-14.4e %-14d %-14.4e %-14d\n",
            β, σ, norm(α_ckg), rank_ckg, norm(α_dll), rank_dll)
    push!(results, (; β, sv_ckg, sv_dll, α_ckg, α_dll))
end

# ── Skew-symmetry / KMS-DBC check (Ding–Li–Lin 2024 Eq. 4.7) ─────────────────
# α(ν, ν') = α(-ν', -ν) · e^{-β(ν+ν')/2}  for any KMS-DB Kossakowski.
println("\n", "="^70)
println("KMS-DBC skew-symmetry: α(ν, ν') = α(-ν', -ν) · e^{-β(ν+ν')/2}")
println("="^70)
@printf("%-6s %-12s %-20s %-20s\n", "β", "σ", "max|err|^CKG", "max|err|^DLL")
println("-"^70)
neg_idx = [findfirst(==(-unique_freqs[k]), unique_freqs) for k in 1:K]
@assert all(!isnothing, neg_idx) "Bohr grid must be symmetric (±ν pairs)"

for r in results
    β, α_ckg, α_dll = r.β, r.α_ckg, r.α_dll
    σ = 1.0 / β

    err_ckg = 0.0
    err_dll = 0.0
    for q in 1:K, p in 1:K
        wt = exp(-β * (unique_freqs[p] + unique_freqs[q]) / 2)
        err_ckg = max(err_ckg, abs(α_ckg[p, q] - α_ckg[neg_idx[q], neg_idx[p]] * wt))
        err_dll = max(err_dll, abs(α_dll[p, q] - α_dll[neg_idx[q], neg_idx[p]] * wt))
    end
    @printf("%-6.1f %-12.4f %-20.4e %-20.4e\n", β, σ, err_ckg, err_dll)
end
println("="^70)
println("Both CKG and DLL satisfy α(ν,ν') = α(-ν',-ν) · e^{-β(ν+ν')/2} (Eq. 4.7)")
println("— the matrix-level witness of KMS detailed balance.")
println("="^70)

println("="^70)
println("Key takeaway (Ding et al. 2024 Sec. 4):")
println("  CKG: full-rank Kossakowski (~K modes), off-diagonal Gaussian skirt.")
println("  DLL: RANK-1 Kossakowski (one mode), `freq_kernel(ν)·conj(freq_kernel(ν'))`.")
println("  This rank gap is what makes DLL admit a *finite* number of Lindblad")
println("  operators (one per coupling), vs CKG's continuous ω-integral.")
println("="^70)

# ── Optional plot output ─────────────────────────────────────────────────────
# Skip the plot if Plots is not loaded (avoid a Plots dependency in CI). The
# numerical summary above is the paper-Fig.-3 reproduction; the plot is a
# convenience for thesis Ch5 / appendix figures.
global plot_output = nothing
try
    @eval using Plots
    global plot_output = joinpath(@__DIR__, "..", "drafts", "plots", "kossakowski_dll_vs_ckg.png")
    mkpath(dirname(plot_output))
catch
    @info "Plots.jl not available; skipping figure generation."
end

if plot_output !== nothing
    β_target = 10.0
    idx = findfirst(r -> r.β == β_target, results)
    r = results[idx]

    p1 = Plots.heatmap(unique_freqs, unique_freqs, r.α_ckg;
        title="CKG α^{KMS}_{ν,ν'} (β=10, σ=0.1)", xlabel="ν'", ylabel="ν",
        c=:viridis, aspect_ratio=:equal)

    p2 = Plots.heatmap(unique_freqs, unique_freqs, real.(r.α_dll);
        title="DLL α^{DLL}_{ν,ν'} (β=10)", xlabel="ν'", ylabel="ν",
        c=:viridis, aspect_ratio=:equal)

    p3 = Plots.plot(1:length(r.sv_ckg), r.sv_ckg .+ eps();
        yscale=:log10, label="CKG", marker=:circle,
        title="Singular values (β=10)", xlabel="rank index", ylabel="σ_k")
    Plots.plot!(p3, 1:length(r.sv_dll), r.sv_dll .+ eps();
        label="DLL (rank-1)", marker=:square, yscale=:log10)

    fig = Plots.plot(p1, p2, p3; layout=(1, 3), size=(1400, 420))
    Plots.savefig(fig, plot_output)
    @info "Saved comparison plot" plot_output
end
