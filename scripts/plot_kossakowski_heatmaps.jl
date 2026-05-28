#!/usr/bin/env julia
#
# Plot 3: Kossakowski matrix α_{ν1,ν2} heat maps for the disordered
# Heisenberg Hamiltonian.
#
# Strategy: build a 1×3 panel that pedagogically isolates the
# off-diagonal contribution Chen et al. add over Davies.
#
#   Panel 1 — Davies (GNS DB) at β_phys = 0.5:
#             α^GNS_{ν,ν} on the diagonal, zero elsewhere. Davies only
#             admits GNS detailed balance; per Corollary cor:kossakowski-shift
#             of 2_methods.tex (line 413), α^GNS_{ν1,ν2} = α^KMS_{ν1−β_alg·σ²/2,
#             ν2−β_alg·σ²/2}, so the diagonal entries are evaluated as
#             create_alpha(ν − β_alg·σ²/2, ν − β_alg·σ²/2, β_alg, σ, 0, s).
#             The diagonal-only structure does not change qualitatively
#             with β, so a single panel suffices.
#   Panel 2 — KMS Chen at β_phys = 0.5: full α^KMS_{ν1,ν2} = create_alpha(...).
#   Panel 3 — KMS Chen at β_phys = 1.0: full α^KMS_{ν1,ν2} (colder).
#
# qf-6vr (β_phys / β_alg split): all formulas above take β = β_alg, where
# β_alg = β_phys · raw.rescaling_factor. The user-facing temperature is
# β_phys (legacy values 10, 20 were actually β_alg — see CLAUDE.md);
# under the canonical β_phys grid {0.25, 0.5, 1.0} they map to
# β_phys ∈ {0.5, 1.0} for the warm/cold pair.
#
# All three panels use the smooth Metropolis filter (a = 0, s = 0.25),
# matching scripts/plot_transition_weights.jl and plot_b_plus_b_minus.jl.
# Colour scale is shared across all three panels so the visual gap
# between Davies (just diagonal) and KMS Chen (diagonal + skirt) is
# directly the off-diagonal contribution.
#
# System: n = 3 disordered Heisenberg, find_ideal disorder realization
# (live seed scan: pick the seed that MAXIMIZES the smallest positive
# Bohr gap ν_min, giving the most uniform spectrum and hence the
# smoothest visual rendering of the α heatmap). This reconstructs the
# behavior of the legacy `find_ideal_heisenberg` selector (removed in
# qf-2kd) inline; for n=3 the 1000-seed scan runs in a few seconds.
# At BATCH_SIZE=1000 this lands on seed=519 with ν_min ≈ 0.033 (rescaled),
# ≈ 1.9× the find_typical fixture's ν_min ≈ 0.017 and producing zero
# cluster gaps above σ — the figure renders as a continuous Gaussian
# skirt with no visible block boundaries.
#
# This is a REPRESENTATION-ONLY figure: numerical thesis claims use the
# canonical find_typical fixture on disk; find_ideal is fine here because
# the figure illustrates the qualitative shape of α(ν₁, ν₂), not a
# specific numerical comparison.
#
# PHYSICS CHECK: σ = 1/β_alg follows the Chen et al. canonical setting on
# the rescaled spectrum and matches plots 1 and 2 in this epic. With
# σ = 1/β_alg the KMS DB Kossakowski matrix is skew-symmetric
# (eq. line 160 of 2_methods.tex).
#
# PHYSICS CHECK: s = 0.25 matches plots 1, 2 of this epic.
#
# PHYSICS CHECK: β_phys values 0.5, 1.0 (warm vs cold) on the canonical
# β_phys grid. With n=3 rescaling_factor ≈ 21.86, β_alg ≈ {10.9, 21.9}
# — very close to the legacy β_alg pair {10, 20}, so the physical
# regime is essentially preserved. Bohr spectrum (rescaled) spans
# [-0.45, 0.45], σ = 1/β_alg. The Gaussian factor exp(-(ν1-ν2)²/(8σ²))
# controls off-diagonal decay; the cold panel shows a tighter
# off-diagonal skirt than the warm one.
#
# Output: drafts/plots/kossakowski_heatmap_ckg.{pdf,png}
#   (the `_ckg` suffix mirrors the DLL companion `plot_kossakowski_dll_heatmaps.jl`
#    so both are addressable as `kossakowski_heatmap_{ckg,dll}*` in the thesis).
#
# Colour map: "ocean inverted" gradient (deep navy → mint), per
# memory: reference_gradient_palettes. Same dark-base / bright-peak feel as
# `:inferno` but anchored at navy (#003147) instead of pure black, which reads
# more elegantly against the white thesis page.
#
# Usage:  julia --project scripts/plot_kossakowski_heatmaps.jl

using Printf
using LinearAlgebra
using Plots
using QuantumFurnace

# ── Parameters ────────────────────────────────────────────────────────────────
const num_qubits     = 3
const a_reg          = 0.0
const s_smooth       = 0.25
# qf-6vr canonical β_phys grid: legacy β_alg {5, 10, 20} -> β_phys {0.25, 0.5, 1.0}.
# Warm/cold pair used here mirrors the legacy {10, 20} entries.
const β_phys_warm    = 0.5
const β_phys_cold    = 1.0

# ── find_ideal_heisenberg reconstruction: max-ν_min over seed scan ───────────
const FIND_IDEAL_BATCH = 1000
const coeffs            = [-1.0, -1.0, -1.0]
const disordering_terms = Vector{Matrix{ComplexF64}}[[Z], [Z, Z]]

function find_ideal_seed(n::Int, batch::Int)
    best_seed = 0
    best_nu   = -Inf
    for s in 1:batch
        raw = build_heis_1d(n, coeffs;
                            seed              = s,
                            periodic          = true,
                            disordering_terms = disordering_terms,
                            disorder_strength = 1.0)
        if raw.nu_min > best_nu
            best_nu   = raw.nu_min
            best_seed = s
        end
    end
    return best_seed, best_nu
end

@info "find_ideal scan: n=$num_qubits, batch=$FIND_IDEAL_BATCH, maximizing ν_min(rescaled)…"
ideal_seed, ideal_nu = find_ideal_seed(num_qubits, FIND_IDEAL_BATCH)
@printf("Selected seed = %d   ν_min(rescaled) = %.6f\n", ideal_seed, ideal_nu)

raw = build_heis_1d(num_qubits, coeffs;
                    seed              = ideal_seed,
                    periodic          = true,
                    disordering_terms = disordering_terms,
                    disorder_strength = 1.0)
# β passed to HamHam here only initialises the (unused) Gibbs state; α takes β_alg
# explicitly below.
ham = HamHam(raw, 1.0)
const rescale  = ham.rescaling_factor
const β_alg_warm = β_phys_warm * rescale
const β_alg_cold = β_phys_cold * rescale

unique_freqs = sort(collect(keys(ham.bohr_dict)))
Nfreq = length(unique_freqs)
@printf("Fixture: find_ideal n=%d seed=%d (live scan of %d seeds)\n",
        num_qubits, ideal_seed, FIND_IDEAL_BATCH)
@printf("rescaling_factor = %.4f\n", rescale)
@printf("β_phys = %.3f -> β_alg = %.3f   (warm)\n", β_phys_warm, β_alg_warm)
@printf("β_phys = %.3f -> β_alg = %.3f   (cold)\n", β_phys_cold, β_alg_cold)
@printf("Unique Bohr frequencies: %d   (expected %d for n=%d)\n",
        Nfreq, 2^num_qubits * (2^num_qubits - 1) + 1, num_qubits)
@printf("ν_min (smallest positive Bohr): %.6f\n", ham.nu_min)
@printf("ν range (rescaled): [%.4f, %.4f]\n", first(unique_freqs), last(unique_freqs))

# ── Compute Kossakowski matrices α_{ν1,ν2} ───────────────────────────────────
# β below is always β_alg (algorithm-space inverse temperature, σ = 1/β_alg
# on the rescaled Bohr grid).
function compute_alpha_matrix(unique_freqs::Vector{Float64}, β_alg::Real;
                              σ::Real=1/β_alg, a::Real=a_reg, s::Real=s_smooth)
    N = length(unique_freqs)
    A = zeros(Float64, N, N)
    for j in 1:N, i in 1:N
        A[i, j] = create_alpha(unique_freqs[i], unique_freqs[j], β_alg, σ, a, s)
    end
    return A
end

@info "Computing α at β_phys=$(β_phys_warm) (β_alg=$(round(β_alg_warm, digits=3))), σ=$(round(1/β_alg_warm, digits=4))…"
α_warm = compute_alpha_matrix(unique_freqs, β_alg_warm)
@info "Computing α at β_phys=$(β_phys_cold) (β_alg=$(round(β_alg_cold, digits=3))), σ=$(round(1/β_alg_cold, digits=4))…"
α_cold = compute_alpha_matrix(unique_freqs, β_alg_cold)

# Davies (GNS DB) panel: diagonal of α^GNS at the warm β, off-diagonals zeroed.
# Per Corollary cor:kossakowski-shift (2_methods.tex line 413):
#   α^GNS_{ν1,ν2} = α^KMS_{ν1 − β_alg·σ²/2, ν2 − β_alg·σ²/2},
# so on the diagonal α^GNS(ν, ν) = α^KMS(ν − β_alg·σ²/2, ν − β_alg·σ²/2).
# At σ = 1/β_alg the shift β_alg·σ²/2 = 1/(2β_alg) is small compared to
# the rescaled ν range, but it correctly identifies the Davies generator
# with GNS DB — the only kind of detailed balance Davies satisfies.
σ_warm = 1 / β_alg_warm
shift_warm = β_alg_warm * σ_warm^2 / 2   # = 1/(2β_alg) for σ = 1/β_alg
α_davies = zeros(size(α_warm))
for i in 1:Nfreq
    ν = unique_freqs[i]
    α_davies[i, i] = create_alpha(ν - shift_warm, ν - shift_warm,
                                   β_alg_warm, σ_warm, a_reg, s_smooth)
end
@printf("Davies shift β_alg·σ²/2 at β_phys=%.3f: %.4f   (ν range [%.2f, %.2f])\n",
        β_phys_warm, shift_warm, first(unique_freqs), last(unique_freqs))

# Sanity: skew symmetry of α at σ = 1/β_alg (Eq. line 160 of 2_methods.tex):
#   α(ν1, ν2) = α(-ν2, -ν1) · exp(-β_alg (ν1 + ν2) / 2)
function skew_residual(α::Matrix, freqs::Vector, β_alg::Real)
    res = 0.0
    for j in 1:length(freqs), i in 1:length(freqs)
        ν1, ν2 = freqs[i], freqs[j]
        i_mirror = findfirst(≈(-ν2), freqs)
        j_mirror = findfirst(≈(-ν1), freqs)
        i_mirror === nothing && continue
        j_mirror === nothing && continue
        rhs = α[i_mirror, j_mirror] * exp(-β_alg * (ν1 + ν2) / 2)
        res = max(res, abs(α[i, j] - rhs))
    end
    return res
end
@printf("Skew-symmetry residual α_warm: %.2e   (should be ≪ 1)\n",
        skew_residual(α_warm, unique_freqs, β_alg_warm))
@printf("Skew-symmetry residual α_cold: %.2e   (should be ≪ 1)\n",
        skew_residual(α_cold, unique_freqs, β_alg_cold))

# Shared colour scale across all three panels
α_max = max(maximum(α_warm), maximum(α_cold))
α_min = 0.0
@printf("α range across panels: [%.3e, %.3e]\n", α_min, α_max)

# Off-diagonal magnitude diagnostics
function offdiag_max(α::Matrix)
    m = 0.0
    for j in 1:size(α, 2), i in 1:size(α, 1)
        i == j && continue
        m = max(m, α[i, j])
    end
    return m
end
@printf("max diag(α_warm) = %.3e   max off-diag(α_warm) = %.3e   ratio = %.2f\n",
        maximum(diag(α_warm)), offdiag_max(α_warm),
        maximum(diag(α_warm)) / offdiag_max(α_warm))
@printf("max diag(α_cold) = %.3e   max off-diag(α_cold) = %.3e   ratio = %.2f\n",
        maximum(diag(α_cold)), offdiag_max(α_cold),
        maximum(diag(α_cold)) / offdiag_max(α_cold))

# ── Plotting defaults (match other thesis-plot scripts) ──────────────────────
default(
    fontfamily              = "Computer Modern",
    titlefontsize           = 11,
    guidefontsize           = 11,
    tickfontsize            = 8,
    legendfontsize          = 8,
    framestyle              = :axes,
    grid                    = false,
    margin                  = 1Plots.mm,
    # Hide axis lines and tick marks while keeping tick *labels* visible:
    # ν values render in the default text colour, but the L-shape of
    # bottom + left axis lines is painted white (i.e. invisible on a
    # white canvas).
    foreground_color_axis   = :white,
    foreground_color_border = :white,
    tick_direction          = :none,
)

# Index-based axes with ~7 sparse ν labels — uniform cells are easier to read
# than the non-uniformly-spaced raw Bohr frequencies.
tick_idx = round.(Int, range(1, Nfreq; length=7))
tick_lbl = [@sprintf("%.2f", unique_freqs[i]) for i in tick_idx]

# ── Colour map ────────────────────────────────────────────────────────────────
# "Ocean inverted": deep navy (#003147) zero-anchor → mint (#B7E6A5) peaks.
# Same dark-base / bright-peak feel as `:inferno`, but with a non-black base
# so the figure stays elegant on white paper.
const CMAP = cgrad(["#003147", "#045275", "#00718B", "#089099",
                    "#46AEA0", "#7CCBA2", "#B7E6A5"])

# All three panels: no individual colorbar, identical bottom margin so they
# render at the same physical size; ν₂ is labelled only on the middle panel
# (visually centred, serves as the shared x-axis label for all three).
function make_heatmap(α::Matrix, title::String;
                      show_y::Bool=true, show_xlabel::Bool=false,
                      left_margin=0Plots.mm, right_margin=0Plots.mm)
    return heatmap(
        α;
        c             = CMAP,
        clims         = (α_min, α_max),
        xticks        = (tick_idx, tick_lbl),
        yticks        = show_y ? (tick_idx, tick_lbl) : nothing,
        xlabel        = show_xlabel ? "\$\\nu_2\$" : "",
        ylabel        = show_y ? "\$\\nu_1\$" : "",
        title         = title,
        colorbar      = false,
        aspect_ratio  = :equal,
        xrotation     = 45,
        framestyle    = :axes,
        left_margin   = left_margin,
        right_margin  = right_margin,
        bottom_margin = 7Plots.mm,
        top_margin    = 1Plots.mm,
    )
end

# Margins tightened: only p1 needs left-margin for the "ν₁" ylabel; inter-panel
# gaps go to 0 since p2/p3 have show_y=false (no y-tick labels to collide with
# p1/p2's right edges) and the heatmaps anyway sit inside white-axis frames
# that do not visually clash when butted together.
p1 = make_heatmap(α_davies, "Davies (GNS DB), \$\\beta=$(β_phys_warm)\$";
                  show_y=true,  show_xlabel=false,
                  left_margin=7Plots.mm, right_margin=0Plots.mm)
p2 = make_heatmap(α_warm,   "CKG (KMS DB), \$\\beta=$(β_phys_warm)\$";
                  show_y=false, show_xlabel=true,
                  left_margin=0Plots.mm, right_margin=0Plots.mm)
p3 = make_heatmap(α_cold,   "CKG (KMS DB), \$\\beta=$(β_phys_cold)\$";
                  show_y=false, show_xlabel=false,
                  left_margin=0Plots.mm, right_margin=0Plots.mm)

# Shared vertical colorbar — thin gradient column as a 4th subplot.
n_cb     = 256
cb_y     = collect(range(α_min, α_max; length=n_cb))
cb_data  = repeat(reshape(cb_y, n_cb, 1), 1, 2)
cb_ticks = round.(range(α_min, α_max; length=5); digits=2)
cbar_panel = heatmap(
    [0.0, 1.0], cb_y, cb_data;
    c              = CMAP,
    clims          = (α_min, α_max),
    xticks         = nothing,
    yticks         = cb_ticks,
    ymirror        = true,
    xlabel         = "\$\\alpha_{\\nu_1,\\nu_2}\$",
    title          = "",
    colorbar       = false,
    framestyle     = :axes,
    grid           = false,
    aspect_ratio   = :auto,
    top_margin     = 4Plots.mm,
    bottom_margin  = 8Plots.mm,
    left_margin    = 0Plots.mm,
    right_margin   = 0Plots.mm,
)

l = @layout [grid(1,3) b{0.025w}]
plt_all = plot(p1, p2, p3, cbar_panel;
               layout = l,
               size   = (1300, 460),
               dpi    = 200,
               # Leading "\n" pushes the visible title text down by one line
               # height, clearing the GR backend's top-clip of ascenders on
               # letters like "K" and "M". The blank first line takes the
               # canvas-top clip; the visible title sits one line below it,
               # safely inside the canvas. `plot_titlefontvalign = :bottom`
               # was tried first but is ignored by GR for supertitles.
               # vspan tightened to 0.12 (was 0.20): the 2-line title only
               # needs ~10 % of the canvas; the freed vertical space goes
               # into the heatmap squares (aspect_ratio=:equal is height-
               # limited so this is the dominant lever for heatmap size).
               plot_title = "\nKossakowski matrix \$\\alpha_{\\nu_1,\\nu_2}\$ (smooth Metropolis, 1D Heisenberg)",
               plot_titlefontsize = 15,
               plot_titlevspan   = 0.16)

out_pdf = joinpath(@__DIR__, "..", "drafts", "plots", "kossakowski_heatmap_ckg.pdf")
out_png = joinpath(@__DIR__, "..", "drafts", "plots", "kossakowski_heatmap_ckg.png")
mkpath(dirname(out_pdf))
savefig(plt_all, out_pdf)
savefig(plt_all, out_png)
@printf("Saved %s\n", relpath(out_pdf, joinpath(@__DIR__, "..")))
@printf("Saved %s\n", relpath(out_png, joinpath(@__DIR__, "..")))
