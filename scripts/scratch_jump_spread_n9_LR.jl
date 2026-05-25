#!/usr/bin/env julia
# scratch_jump_spread_n9_LR.jl  (qf-vr4)
#
# Appendix figure: Heisenberg-picture operator spreading of a single-site
# Pauli A = σ^x_5 on the 1D AFM Heisenberg chain (n = 9, seed = 46, Z+ZZ
# disorder, periodic) at the two TimeDomain-relevant times t_max and
# t_median pulled from `_precompute_labels` for each canonical
# β_phys ∈ {0.25, 0.5, 1.0}, overlaid with the Lieb–Robinson lightcone for
# the 1D XXX chain.
#
# This is a PURE OPERATOR-SPREAD measurement. We do NOT apply any
# Lindbladian L to any state. We do NOT apply any channel Φ_δ to any state.
# We do NOT compute jump_oft, jump_contribution, or any per-step dissipator
# matvec. We do NOT need a ρ_0 — Heisenberg-picture operator spreading is
# state-independent. The TimeDomain config enters ONLY as the source of the
# physically-meaningful t-values for each (n, β_phys) cell, via the
# truncated Gaussian time grid that the simulator's OFT actually visits.
#
# Per-site weight definition:
#   w_k(t) := (1/d) [‖A(t)‖_HS² − (1/2) ‖Tr_k[A(t)]‖_HS²], with d = 2^n.
# By the Pauli decomposition A(t) = Σ_P c_P P this equals
#   w_k(t) = Σ_{P : P_k ≠ I_k} |c_P|²,
# i.e. the sum of squared Pauli coefficients over strings with a
# non-identity Pauli at site k. w_k = 1 iff A(t) has full Pauli support on
# site k; w_k = 0 iff A(t) is identity on site k. Σ_k w_k ∈ [1, n] is the
# "average Pauli weight" of A(t). At t = 0 with A = σ^x_5: w_5 = 1, others 0.
#
# Lieb–Robinson overlay: tight 1D XXX velocity v_LR = 2π (spinon group
# velocity in our Pauli convention h_{i,i+1} = X⊗X + Y⊗Y + Z⊗Z = 4 S·S),
# see drafts/numerics-extras/lieb-robinson-notes.md.
#
# Outputs:
#   drafts/figures/numerics/jump_spread_LR_n9.{png,pdf}   3-panel figure
#   drafts/figures/numerics/jump_spread_LR_n9.bson         sidecar
#
# Usage (from repo root):
#   julia --project scripts/scratch_jump_spread_n9_LR.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using LinearAlgebra
using SparseArrays
using Statistics
using Random
using Printf
using BSON
using QuantumFurnace
using QuantumFurnace: X, Y, Z, HamHam, Config, Lindbladian, TimeDomain, KMS,
                      beta_alg
# Internal API consumed for the t-grid; documented in furnace_utensils.jl.
# We touch only label construction — never the simulator hot paths.
using QuantumFurnace: _precompute_labels, _truncate_time_labels_for_oft
import QuantumFurnace
ENV["GKSwstype"] = "100"
using Plots

# --- Hard-coded constants ------------------------------------------------

const N            = 9                        # chain length
const SEED         = 46                       # canonical single-seed (qf-e4z.34)
const SITE_X0      = 5                        # middle site (1-indexed)
const PAULI_LABEL  = "X"                      # σ^x is the canonical single-site jump
const BETA_PHYS_ALL = (0.25, 0.5, 1.0)        # canonical β_phys grid
const R_D          = 8                        # canonical TimeDomain dissipative register
const TAIL_C       = 8.0                      # ω-range = 2(‖H‖ + TAIL_C · σ_alg)  (matches qf-e4z.35)
const SMOOTH_S     = 0.25                     # smooth-Metro s
const SMOOTH_A     = 0.0                      # smooth-Metro a
const WEIGHT_THRESH = 1e-3                    # threshold on Gaussian time weight (vs max)
const OFT_TOL      = 1e-12                    # tolerance for OFT time-label truncation

# Lieb–Robinson velocities in our Pauli convention.
# v_LR_tight   = 2π          (spinon group velocity for 1D XXX, our units)
# v_LR_BH      = 12 · e      (BH 2006 rigorous bound; not drawn, only quoted)
const V_LR_TIGHT   = 2.0 * pi
const V_LR_BH      = 12.0 * exp(1.0)
const LR_DECAY_RATE = 2.0                     # site-weight envelope: exp(-2(|k-x0| - v_LR t))

# --- File paths ----------------------------------------------------------

const REPO_ROOT    = joinpath(@__DIR__, "..")
const FIXTURE_DIR  = joinpath(REPO_ROOT, "scripts", "output", "multiseed_fixtures")
mkpath(FIXTURE_DIR)
const FIXTURE_PATH = joinpath(FIXTURE_DIR,
    "heis_xxx_disordered_periodic_n$(N)_seed$(SEED).bson")
const FIG_DIR      = joinpath(REPO_ROOT, "drafts", "figures", "numerics")
mkpath(FIG_DIR)
const FIG_BASE     = joinpath(FIG_DIR, "jump_spread_LR_n9")

# --- Fixture: generate n=9, seed=46 if missing ---------------------------

# Same recipe as scripts/scratch_multiseed_disordered_fixtures.jl
# (`run_1d_heisenberg`): XXX with [Z]+[Z,Z] disorder, ε = 0.1.
function _ensure_fixture()
    if isfile(FIXTURE_PATH)
        @printf("[fixture] %s already exists; loading.\n", basename(FIXTURE_PATH))
        return
    end
    @printf("[fixture] Generating %s (n=%d, seed=%d, Z+ZZ ε=0.1)\n",
            basename(FIXTURE_PATH), N, SEED)
    coeffs = [1.0, 1.0, 1.0]
    raw = QuantumFurnace.build_heis_1d(N, coeffs; seed=SEED, periodic=true,
                                        disorder_strength=0.1)
    BSON.bson(FIXTURE_PATH, hamiltonian=raw)
    @printf("  saved (R=%.4f, nu_min=%.4e)\n", raw.rescaling_factor, raw.nu_min)
end

function _load_ham(beta_phys::Real)
    raw = QuantumFurnace._parse_hamiltonian_bson(FIXTURE_PATH)
    return HamHam(raw; beta_phys = float(beta_phys))
end

# --- TimeDomain config for the (n, β_phys) cell --------------------------

function _build_cfg(ham::HamHam, beta_phys::Real)
    β_alg = beta_alg(ham, float(beta_phys))
    σ_alg = 1.0 / β_alg                                  # canonical CKG σ at this cell
    H_norm = maximum(abs, ham.eigvals)
    ω_range = 2.0 * (H_norm + TAIL_C * σ_alg)
    w0_D = ω_range / 2.0^R_D
    t0_D = 2π / (2.0^R_D * w0_D)
    cfg = Config(
        sim = Lindbladian(),
        domain = TimeDomain(),
        construction = KMS(),
        num_qubits = N,
        with_linear_combination = true,
        beta = β_alg,
        beta_phys = float(beta_phys),
        sigma = σ_alg,
        a = SMOOTH_A, s = SMOOTH_S,
        gaussian_parameters = (nothing, nothing),
        num_energy_bits_D = R_D, w0_D = w0_D, t0_D = t0_D,
        num_trotter_steps_per_t0 = 1,
        filter = nothing,                                # CKG smooth-Metro Gaussian
    )
    return cfg, β_alg, σ_alg, w0_D, t0_D
end

# --- Pick t_max and t_median from the weight-thresholded TimeDomain grid -

function _pick_t_values(cfg::Config, σ_alg::Real)
    # The simulator's truncated OFT time grid for this cell.
    _, time_labels_full = _precompute_labels(cfg)
    oft_labels = _truncate_time_labels_for_oft(time_labels_full, σ_alg;
                                               tolerance=OFT_TOL)
    weights = @. exp(-σ_alg^2 * oft_labels^2)
    w_max = maximum(weights)
    survives = weights .> (WEIGHT_THRESH * w_max)
    surv_t = oft_labels[survives]
    surv_w = weights[survives]
    abs_t = abs.(surv_t)
    # Take t_max = max |t| in the surviving window (the simulator's
    # "furthest physical time"), and t_median = median |t| (typical time).
    t_max_alg = maximum(abs_t)
    t_median_alg = median(abs_t)
    return (
        oft_labels = oft_labels, weights = weights, survives = survives,
        t_max_alg = t_max_alg, t_median_alg = t_median_alg,
    )
end

# --- Build A = σ^x on site x0 of n-qubit chain (computational basis) -----

function _single_site_pauli(P::AbstractMatrix{<:Complex}, x0::Integer, n::Integer)
    @assert 1 <= x0 <= n
    tensor_list = SparseMatrixCSC{ComplexF64}[]
    for k in 1:n
        push!(tensor_list, k == x0 ? sparse(ComplexF64.(P)) : sparse(I, 2, 2))
    end
    return Matrix(kron(tensor_list...))
end

# --- Partial trace over site k (1-indexed) -------------------------------
#
# For an n-qubit operator A in the computational basis with site k traced
# out, the result is a (2^(n-1))×(2^(n-1)) matrix on the remaining n−1 sites.
# We use the standard reshape-and-sum identity: bit k in the computational
# basis is the (n-k)-th bit from the right (lattice site k in our kron
# convention pads as I⊗…⊗P_k⊗…⊗I starting from the left).
function _partial_trace_site(A::AbstractMatrix{<:Complex}, k::Integer, n::Integer)
    d = 2^n
    @assert size(A) == (d, d)
    # Lattice site k corresponds to the (n-k+1)-th tensor factor counting from
    # the right end; bit position in the binary index is `bit_pos = n - k`.
    bit_pos = n - k
    bit_mask = 1 << bit_pos
    d_left  = 2^(n - bit_pos - 1)        # number of basis states for sites < k (in kron order: those on the LEFT)
    d_right = 2^bit_pos                  # number of basis states for sites > k
    @assert d_left * 2 * d_right == d
    out = zeros(eltype(A), d ÷ 2, d ÷ 2)
    # For each high-bit value (sites < k), low-bit value (sites > k):
    @inbounds for i_lo in 0:(d_right - 1)
        for j_lo in 0:(d_right - 1)
            for i_hi in 0:(d_left - 1)
                for j_hi in 0:(d_left - 1)
                    s = zero(eltype(A))
                    for b in 0:1
                        # Build full index: (i_hi bits) (b) (i_lo bits)
                        i_full = i_hi * (2 * d_right) + b * d_right + i_lo + 1
                        j_full = j_hi * (2 * d_right) + b * d_right + j_lo + 1
                        s += A[i_full, j_full]
                    end
                    # Output index without bit k:
                    i_red = i_hi * d_right + i_lo + 1
                    j_red = j_hi * d_right + j_lo + 1
                    out[i_red, j_red] += s
                end
            end
        end
    end
    return out
end

# --- Site weights w_k from A(t) ------------------------------------------

function _site_weights(A::AbstractMatrix{<:Complex}, n::Integer)
    d = 2^n
    A_hs2 = real(tr(A' * A))                 # ‖A‖_HS² = tr(A† A)
    w = zeros(Float64, n)
    for k in 1:n
        Tk = _partial_trace_site(A, k, n)
        Tk_hs2 = real(tr(Tk' * Tk))
        w[k] = (A_hs2 - 0.5 * Tk_hs2) / d
    end
    return w, A_hs2
end

# --- Heisenberg-picture A(t) under H_alg --------------------------------

function _evolve_pauli(A0::AbstractMatrix{<:Complex}, H_alg::AbstractMatrix{<:Complex}, t_alg::Real)
    # U = exp(-i H_alg t_alg). Dense; n=9 → 512×512, microseconds.
    U = exp(-1im * t_alg * Hermitian(H_alg))
    return U' * A0 * U
end

# --- LR envelope on the figure (site-weight bound) -----------------------
#
# Inside lightcone |k - x0| ≤ v_LR · |t_phys|: w_k ≲ 1 (drawn as a flat
# ceiling). Outside lightcone: w_k ≲ exp(-LR_DECAY_RATE · (|k - x0| -
# v_LR · |t_phys|)). This is schematic and is meant only to visually
# anchor whether the data respects causality, not to give a sharp bound.

function _lr_envelope(sites::AbstractVector{<:Integer}, x0::Integer,
                       v_LR::Real, t_phys::Real)
    return [min(1.0, exp(-LR_DECAY_RATE * (abs(k - x0) - v_LR * t_phys)))
            for k in sites]
end

# --- Per-cell computation ------------------------------------------------

function _run_cell(beta_phys::Real)
    ham = _load_ham(beta_phys)
    cfg, β_alg, σ_alg, w0_D, t0_D = _build_cfg(ham, beta_phys)
    R = ham.rescaling_factor

    info = _pick_t_values(cfg, σ_alg)
    t_max_alg, t_median_alg = info.t_max_alg, info.t_median_alg
    t_max_phys = t_max_alg / R
    t_median_phys = t_median_alg / R

    # Operator A = σ^x on site x0. Heisenberg-picture evolution under H_alg.
    A0 = _single_site_pauli(X, SITE_X0, N)
    H_alg = ham.data

    A_at_tmax    = _evolve_pauli(A0, H_alg, t_max_alg)
    A_at_tmedian = _evolve_pauli(A0, H_alg, t_median_alg)
    A_at_zero    = A0

    w_tmax, _    = _site_weights(A_at_tmax, N)
    w_tmedian, _ = _site_weights(A_at_tmedian, N)
    w_zero, _    = _site_weights(A_at_zero, N)

    return (;
        beta_phys = float(beta_phys), beta_alg = β_alg, sigma_alg = σ_alg,
        R = R, w0_D = w0_D, t0_D = t0_D, r_D = R_D,
        t_max_alg = t_max_alg, t_max_phys = t_max_phys,
        t_median_alg = t_median_alg, t_median_phys = t_median_phys,
        oft_labels = info.oft_labels, oft_weights = info.weights,
        oft_survives = info.survives,
        w_zero = w_zero, w_tmedian = w_tmedian, w_tmax = w_tmax,
    )
end

# --- Plotting ------------------------------------------------------------

function _plot_panel(cell, sites; show_legend::Bool=false, show_ylabel::Bool=false)
    yticks_arg = show_ylabel ? (0.0:0.25:1.0) : ([0.0, 0.25, 0.5, 0.75, 1.0], ["", "", "", "", ""])
    p = plot(; xlabel="site k",
              ylabel = show_ylabel ? "\$w_k\$" : "",
              left_margin = show_ylabel ? 8Plots.mm : 0Plots.mm,
              right_margin = 0Plots.mm,
              bottom_margin = 8Plots.mm,
              top_margin = 4Plots.mm,
              framestyle=:box, dpi=150,
              xticks = collect(sites),
              yticks = yticks_arg,
              ylims = (0.0, 1.05),
              legend = (show_legend ? :topright : false),
              titlefontsize=11, legendfontsize=8, tickfontsize=8, guidefontsize=12)
    env_max = _lr_envelope(sites, SITE_X0, V_LR_TIGHT, cell.t_max_phys)
    plot!(p, sites, env_max; lc=:grey, ls=:dot, lw=1.5,
          label = show_legend ? "LR envelope (\$v = 2\\pi\$, \$t_{\\max}\$)" : "")

    plot!(p, sites, cell.w_zero;     lc=:black, lw=1.0, ls=:dash,
          label = show_legend ? "\$t = 0\$ (\$\\sigma^x_5\$)" : "")
    plot!(p, sites, cell.w_tmedian;  lc=:blue,  lw=2.0, marker=:circle, ms=3,
          label = show_legend ? @sprintf("\$t_{\\mathrm{mid}} = %.3f\$",
                                          cell.t_median_phys) : "")
    plot!(p, sites, cell.w_tmax;     lc=:red,   lw=2.0, marker=:diamond, ms=3,
          label = show_legend ? @sprintf("\$t_{\\max} = %.3f\$",
                                          cell.t_max_phys) : "")
    title!(p, @sprintf("\$\\beta = %.2g\$", cell.beta_phys))
    return p
end

function _save_figure(cells)
    sites = 1:N
    panels = [_plot_panel(cells[i], sites;
                          show_legend = (i == length(cells)),
                          show_ylabel = (i == 1))
              for i in 1:length(cells)]
    fig = plot(panels...; layout = (1, length(cells)),
                size = (1300, 440), dpi = 150,
                plot_title = "\$A(t) = \\sigma^x_5(t)\$ on 1D Heisenberg",
                plot_titlefontsize = 13,
                plot_titlevspan = 0.10)
    savefig(fig, FIG_BASE * ".png")
    savefig(fig, FIG_BASE * ".pdf")
    @printf("[figure] wrote %s.{png,pdf}\n", FIG_BASE)
end

# --- Main ----------------------------------------------------------------

function main()
    @printf("[start] qf-vr4 jump-spread vs LR on n=%d, seed=%d, x0=%d, A=σ^%s\n",
            N, SEED, SITE_X0, PAULI_LABEL)
    _ensure_fixture()

    cells = NamedTuple[]
    for bp in BETA_PHYS_ALL
        @printf("\n== β_phys = %.2g ==\n", bp)
        cell = _run_cell(bp)
        push!(cells, cell)
        @printf("  R = %.4f, β_alg = %.4f, σ_alg = %.4g\n",
                cell.R, cell.beta_alg, cell.sigma_alg)
        @printf("  r_D = %d, w0_D = %.4g, t0_D = %.4g\n",
                cell.r_D, cell.w0_D, cell.t0_D)
        @printf("  t_max_alg    = %9.4f → t_max_phys    = %.4f\n",
                cell.t_max_alg, cell.t_max_phys)
        @printf("  t_median_alg = %9.4f → t_median_phys = %.4f\n",
                cell.t_median_alg, cell.t_median_phys)
        @printf("  v_LR_tight · t_max_phys = %.3f  (chain length n = %d)\n",
                V_LR_TIGHT * cell.t_max_phys, N)
        @printf("  v_LR_BH    · t_max_phys = %.3f  (rigorous BH bound)\n",
                V_LR_BH    * cell.t_max_phys)
        Σw = sum(cell.w_tmax)
        @printf("  Σ_k w_k(t_max) = %.4f  (avg Pauli weight; t=0 baseline = 1)\n", Σw)
        # Interpretation flag for the appendix:
        if V_LR_TIGHT * cell.t_max_phys > (N - 1)
            @printf("  ▶ HYPOTHESIS CONFIRMED: lightcone (v_LR·t_max) EXCEEDS chain length.\n")
        else
            @printf("  ▶ HYPOTHESIS REFUTED at this cell: lightcone is %.2f sites < n-1 = %d.\n",
                     V_LR_TIGHT * cell.t_max_phys, N - 1)
        end
    end

    _save_figure(cells)

    # Sidecar BSON.
    bson_path = FIG_BASE * ".bson"
    BSON.bson(bson_path,
        n = N, seed = SEED, pauli = PAULI_LABEL, x0 = SITE_X0,
        beta_phys_grid = collect(BETA_PHYS_ALL),
        v_LR_tight = V_LR_TIGHT, v_LR_BH = V_LR_BH,
        r_D = R_D, weight_threshold = WEIGHT_THRESH,
        cells = cells,
    )
    @printf("\n[sidecar] wrote %s\n", bson_path)

    @printf("\n[done] qf-vr4. See %s.{png,pdf} and the lieb-robinson-notes.md.\n", FIG_BASE)
end

main()
