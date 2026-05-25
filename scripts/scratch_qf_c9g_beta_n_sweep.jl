#!/usr/bin/env julia
# scratch_qf_c9g_beta_n_sweep.jl  (qf-c9g — 2D TFIM ordered-phase verification)
#
# Question (raised by colleague 2026-05-25): is the 2D TFIM ORDERED Lindbladian
# gap closing across n=4,6,8 (9.96e-2 → 3.06e-3 → 1.57e-4 at β_phys=2, h=1)
# genuine ordered-phase physics (Z₂ free-energy / surface-tension barrier) or a
# trivial "Gibbs ≈ GS doublet → exponentially small tunnel splitting" artefact
# that would happen in any system with a quasi-degenerate doublet?
#
# Literature (drafts/qf-c9g-lit-survey.md):
#   - Gamarnik–Kiani–Zlokapa 2024 (arXiv:2411.04300): T_mix ≥ exp[n^(1/2-o(1))]
#     for 2D TFIM at constant β ≥ β*, h ≤ h*, exact CKG KMS-DB sampler class.
#   - Classical analogue: Martinelli-Olivieri 1994 — gap ≤ exp(-τ(β)·L).
#   - Diagnostic checklist: (i) closing must persist over β ∈ (β_c, ∞), not just
#     deep cold; (ii) Binder U_4 → 2/3, ⟨m²⟩ → plateau; (iii) R_2 lives in
#     doublet × bulk block (already verified at n=4 β=2 in qf-biz Follow-up A);
#     (iv) Ω(1) gap in paramagnetic control β_phys=0.1.
#
# This script delivers (i) and (iv) via a uniform 15-cell n×β grid; (ii) is read
# from Subtask-B' diagnostics; (iii) was verified by qf-biz and is checked again
# at multiple β in a separate dense-L diagnostic script.
#
# SWEEP A — Krylov spectral gap (matrix-free) and Subtask-B' Gibbs diagnostics
#   (n, β_phys) ∈ {4, 6, 8} × {0.10, 0.25, 0.50, 1.0, 2.0} = 15 cells
#   Canonical CKG smooth-Metro (s=0.25, a=0, r_D=7, kdim=60, howmany=4).
#
# SWEEP B — Hamiltonian doublet structure (β-independent, dense)
#   At each n: dense `eigen(H_phys)`, report Δ_k = E_{k+1} - E_1 for k=1..5.
#
# Comparison datums to verify:
#   - β=0.25, n=4,6,8 gap_phys should ≈ qf-1jj DISORDERED: 1.910, 1.575, 1.741.
#   - β=2.0,  n=4,6,8 gap_phys should ≈ qf-1jj ORDERED:    9.96e-2, 3.06e-3, 1.57e-4.
#   - β=0.10, n=6 gap_phys should ≈ qf-biz Check 1:        4.00.
#   - β=0.50, n=6 gap_phys should ≈ qf-biz Check 1:        0.286.
#   - β=1.0,  n=6 gap_phys should ≈ qf-biz Check 1:        0.0124.
#
# Threading: JULIA_NUM_THREADS=max, BLAS=1 (per `.claude/rules/julia-code.md`).

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using QuantumFurnace
using QuantumFurnace: _construct_2d_heisenberg_base, _construct_disordering_terms,
                     _rescaling_and_shift_factors, pad_term, _build_jump_set,
                     pick_transition
using LinearAlgebra
using BSON
using Printf
using Dates

BLAS.set_num_threads(1)
println("[init] $(Dates.now())  Julia threads = $(Threads.nthreads())  BLAS threads = $(BLAS.get_num_threads())")
flush(stdout)

const OUTPUT_DIR = joinpath(@__DIR__, "output", "qf_c9g_ordered_gap_mechanism")
mkpath(OUTPUT_DIR)

const J_COUPLING = 1.0
const H_FIELD    = 1.0
const KRYLOVDIM  = 60          # safety margin above qf-1jj kdim=40 (qf-65e showed kdim=40 is
                                # tight at n=8 β=2 for the *trajectory* predictor; krylov_spectral_gap
                                # with the qf-8fr GUE-seed worked at kdim=40 in qf-1jj — but +20 margin
                                # is free at d=256).
const HOWMANY    = 4
const TAIL_C     = 8.0
const R_D        = 7            # canonical thesis r_D; matches qf-1jj for direct comparison
const T_C_AT_H1  = 2.07         # Hesselmann–Wessel 2016 (HW Fig. 4)

const N_LIST     = [4, 6, 8]
const GEOM       = Dict(4 => (2, 2), 6 => (2, 3), 8 => (2, 4))   # Lx × Ly
const BETA_GRID  = [0.10, 0.25, 0.50, 1.0, 2.0]
const N_HSPEC_K  = 5            # how many Δ_k to report from H_phys spectrum

# ============================================================================
# TFIM raw fixture (clean, no disorder; qf-8fr 1e-10 GUE Krylov seed handles SSB)
# ============================================================================

function build_clean_tfim_raw(Lx::Integer, Ly::Integer;
                              J::Float64 = 1.0, h::Float64 = 1.0)
    n = Lx * Ly
    H_bond  = _construct_2d_heisenberg_base(Lx, Ly, [[Z, Z]], [-J];
                                             periodic_x = true, periodic_y = true)
    H_field = _construct_disordering_terms([[X]], [fill(-h, n)], n)
    H_phys  = Hermitian(Matrix(H_bond) + Matrix(H_field))
    rescaling_factor, shift = _rescaling_and_shift_factors(H_phys)
    d = 2^n
    rescaled = Matrix(H_phys) ./ rescaling_factor .+ shift * I(d)
    rescaled_eigvals, rescaled_eigvecs = eigen(Hermitian(rescaled))
    nu_min = minimum(diff(rescaled_eigvals))
    return (
        matrix              = rescaled,
        terms               = [[Z, Z], [X]],
        base_coeffs         = [-J / rescaling_factor, -h / rescaling_factor],
        disordering_terms   = nothing,
        disordering_coeffs  = nothing,
        eigvals             = rescaled_eigvals,
        eigvecs             = rescaled_eigvecs,
        nu_min              = nu_min,
        shift               = shift,
        rescaling_factor    = rescaling_factor,
        periodic            = true,
        H_phys_dense        = Matrix{ComplexF64}(H_phys),
    )
end

function build_ckg_energy_cfg(n::Integer, beta_phys_val::Real, ham; r_D::Integer = R_D)
    β_alg_val = beta_alg(ham, float(beta_phys_val))
    σ = 1.0 / β_alg_val
    H_norm = maximum(abs, ham.eigvals)
    omega_range = 2.0 * (H_norm + TAIL_C * σ)
    w0_D = omega_range / 2.0^r_D
    t0_D = 2π / (2.0^r_D * w0_D)
    return Config(
        sim = Lindbladian(),
        domain = EnergyDomain(),
        construction = KMS(),
        num_qubits = n,
        with_linear_combination = true,
        beta = β_alg_val,
        beta_phys = float(beta_phys_val),
        sigma = σ,
        a = 0.0, s = 0.25,
        gaussian_parameters = (nothing, nothing),
        num_energy_bits_D = r_D, w0_D = w0_D, t0_D = t0_D,
        num_trotter_steps_per_t0 = 10,
        filter = nothing,
    )
end

# ============================================================================
# Subtask B' Gibbs diagnostics (cheap dense; reused from qf-biz)
# ============================================================================

function z2_projector(n::Integer)
    P = Matrix{ComplexF64}(pad_term([X], n, 1))
    for i in 2:n
        P = P * Matrix(pad_term([X], n, i))
    end
    return P
end

function gibbs_alg(eigvals_alg, eigvecs_alg, beta_alg_val::Real)
    weights = exp.(-beta_alg_val .* (eigvals_alg .- minimum(eigvals_alg)))
    weights ./= sum(weights)
    return eigvecs_alg * Diagonal(weights) * eigvecs_alg'
end

function vn_entropy(rho::AbstractMatrix)
    F = eigen(Hermitian((rho + rho') / 2))
    p = max.(real(F.values), 0.0)
    p ./= sum(p)
    s = 0.0
    for pi in p
        pi > 1e-18 && (s -= pi * log(pi))
    end
    return s
end

function subtask_b_prime(raw, n::Integer; beta_phys_val::Real)
    d = 2^n
    P = z2_projector(n)
    Z_ops = [Matrix{ComplexF64}(pad_term([Z], n, i)) for i in 1:n]

    β_alg_val = beta_phys_val * raw.rescaling_factor
    σβ = gibbs_alg(raw.eigvals, raw.eigvecs, β_alg_val)

    Pp = (I(d) + P) / 2
    Pm = (I(d) - P) / 2
    w_plus  = real(tr(Pp * σβ))
    w_minus = real(tr(Pm * σβ))

    M  = sum(Z_ops) ./ n
    M2 = M * M
    M4 = M2 * M2
    m1 = real(tr(M  * σβ))
    m2 = real(tr(M2 * σβ))
    m4 = real(tr(M4 * σβ))
    U4 = 1.0 - m4 / (3.0 * m2^2)

    S  = vn_entropy(σβ)
    Sn = S / log(d)
    eff_rank = exp(S)

    Δ1_phys = (raw.eigvals[2] - raw.eigvals[1]) * raw.rescaling_factor
    Δ3_phys = length(raw.eigvals) >= 4 ?
        (raw.eigvals[4] - raw.eigvals[1]) * raw.rescaling_factor : NaN

    return (;
        n,
        beta_phys = beta_phys_val,
        beta_alg = β_alg_val,
        T_over_Tc = (1.0 / beta_phys_val) / T_C_AT_H1,
        w_plus, w_minus,
        m1, m2, m4, U4,
        S, Sn, eff_rank,
        Δ1_phys, Δ3_phys,
    )
end

# ============================================================================
# Krylov L gap (matrix-free, canonical CKG smooth-Metro)
# ============================================================================

function compute_l_gap(raw, n::Integer, beta_phys_val::Real; kdim::Integer = KRYLOVDIM)
    ham = HamHam(raw; beta_phys=float(beta_phys_val))
    jumps = _build_jump_set(ham, n)
    cfg = build_ckg_energy_cfg(n, beta_phys_val, ham)
    t0 = time()
    res = krylov_spectral_gap(cfg, ham, jumps;
        krylovdim = kdim, howmany = HOWMANY, tol = 1e-10)
    wall = time() - t0
    return (;
        gap_alg = res.spectral_gap,
        gap_phys = res.spectral_gap * ham.rescaling_factor,
        eigenvalues = res.eigenvalues,
        converged = res.converged,
        matvec = res.matvec_count,
        wall = wall,
        rescaling_factor = ham.rescaling_factor,
        beta_alg = cfg.beta,
        sigma = cfg.sigma,
    )
end

# ============================================================================
# SWEEP B — Hamiltonian doublet structure (β-independent, cheap)
# ============================================================================

function hamiltonian_spectrum_row(raw, n::Integer; k_max::Integer = N_HSPEC_K)
    F = eigen(Hermitian(raw.H_phys_dense))
    eigvals_phys = real(F.values)
    # already sorted ascending after Hermitian eigen
    Δ_phys = Float64[]
    for k in 1:k_max
        if length(eigvals_phys) >= k + 1
            push!(Δ_phys, eigvals_phys[k+1] - eigvals_phys[1])
        else
            push!(Δ_phys, NaN)
        end
    end
    # bulk gap = E_4 - E_1 (above the doublet); doublet splitting = E_2 - E_1
    bulk_gap = length(eigvals_phys) >= 4 ?
        (eigvals_phys[4] - eigvals_phys[1]) : NaN
    return (;
        n,
        eigvals_low = eigvals_phys[1:min(end, k_max+1)],
        Δ_phys,                   # Δ_k for k=1..k_max
        doublet_split = eigvals_phys[2] - eigvals_phys[1],
        bulk_gap,                  # gap to first non-doublet excitation
        ratio_bulk_over_doublet = bulk_gap / (eigvals_phys[2] - eigvals_phys[1]),
    )
end

# ============================================================================
# SWEEP A — driver
# ============================================================================

function run_sweep_a()
    println("\n" * "="^120)
    println("SWEEP A — Krylov L gap on n × β grid (canonical CKG smooth-Metro)")
    println("="^120)

    # Cache raw fixtures across the n loop
    raws = Dict{Int,NamedTuple}()
    for n in N_LIST
        (Lx, Ly) = GEOM[n]
        raws[n] = build_clean_tfim_raw(Lx, Ly; J=J_COUPLING, h=H_FIELD)
        @printf("[fixture] n=%d (%dx%d)  R = %.4f  shift = %+.4f  nu_min = %.3g\n",
            n, Lx, Ly, raws[n].rescaling_factor, raws[n].shift, raws[n].nu_min)
        flush(stdout)
    end

    @printf("\n%-3s %-7s %-7s %-8s %-9s %-8s %-9s %-9s %-11s %-11s %-9s %-9s %-7s\n",
        "n", "β_phys", "T/T_c", "β_alg", "<M_z²>", "U_4", "S/log d", "eff_rank",
        "Δ_1(phys)", "λ_L(phys)", "λ_L/ΔE_1", "|λ3/λ2|", "wall")
    println("-"^140)

    rows = NamedTuple[]
    for β_phys_val in BETA_GRID
        for n in N_LIST
            raw = raws[n]
            diag = subtask_b_prime(raw, n; beta_phys_val=β_phys_val)
            gap = compute_l_gap(raw, n, β_phys_val)
            # ratio of 3rd-smallest |Re λ| to 2nd-smallest (the "doublet vs bulk" signature)
            re_eigs = sort([abs(real(λ)) for λ in gap.eigenvalues])
            λ3_over_λ2 = length(re_eigs) >= 3 ? re_eigs[3] / max(re_eigs[2], eps()) : NaN
            ratio_gap_dE1 = gap.gap_phys / diag.Δ1_phys

            @printf("%-3d %-7.3g %-7.3g %-8.3g %-9.4f %-8.4f %-9.4f %-9.3g %-11.4e %-11.4e %-9.3g %-9.3g %-7.1f\n",
                n, β_phys_val, diag.T_over_Tc, diag.beta_alg, diag.m2, diag.U4, diag.Sn,
                diag.eff_rank, diag.Δ1_phys, gap.gap_phys, ratio_gap_dE1,
                λ3_over_λ2, gap.wall)
            flush(stdout)

            push!(rows, merge(diag, (
                n = n,
                Lx = GEOM[n][1], Ly = GEOM[n][2],
                gap_alg = gap.gap_alg,
                gap_phys = gap.gap_phys,
                gap_eigenvalues = gap.eigenvalues,
                gap_converged = gap.converged,
                matvec = gap.matvec,
                wall_gap = gap.wall,
                rescaling_factor = gap.rescaling_factor,
                sigma_alg = gap.sigma,
                λ3_over_λ2 = λ3_over_λ2,
                ratio_gap_dE1 = ratio_gap_dE1,
            )))

            # Resilience: re-save BSON after every cell so a script kill doesn't lose data.
            BSON.bson(joinpath(OUTPUT_DIR, "qf_c9g_beta_n_sweep_PARTIAL.bson"), Dict(
                :sweep_a_rows => [Dict(pairs(r)...) for r in rows],
                :cells_done => length(rows),
                :total_cells => length(BETA_GRID) * length(N_LIST),
            ))
        end
        println("-"^140)   # divider between β blocks
        flush(stdout)
    end

    return rows
end

# ============================================================================
# SWEEP B — driver
# ============================================================================

function run_sweep_b()
    println("\n" * "="^120)
    println("SWEEP B — Hamiltonian spectrum (β-independent)")
    println("="^120)
    @printf("%-3s  %-12s %-12s %-12s %-12s %-12s  %-9s\n",
        "n", "Δ_1(phys)", "Δ_2(phys)", "Δ_3(phys)", "Δ_4(phys)", "Δ_5(phys)",
        "bulk/doublet")
    println("-"^100)

    rows = NamedTuple[]
    for n in N_LIST
        (Lx, Ly) = GEOM[n]
        raw = build_clean_tfim_raw(Lx, Ly; J=J_COUPLING, h=H_FIELD)
        r = hamiltonian_spectrum_row(raw, n)
        @printf("%-3d  %-12.4e %-12.4e %-12.4e %-12.4e %-12.4e  %-9.3g\n",
            r.n, r.Δ_phys[1], r.Δ_phys[2], r.Δ_phys[3], r.Δ_phys[4], r.Δ_phys[5],
            r.ratio_bulk_over_doublet)
        push!(rows, r)
    end
    return rows
end

# ============================================================================
# Main
# ============================================================================

function main()
    t0 = time()
    sweep_b_rows = run_sweep_b()   # cheap, do first
    sweep_a_rows = run_sweep_a()
    wall = time() - t0

    bson_path = joinpath(OUTPUT_DIR, "qf_c9g_beta_n_sweep.bson")
    BSON.bson(bson_path, Dict(
        :sweep_a_rows => [Dict(pairs(r)...) for r in sweep_a_rows],
        :sweep_b_rows => [Dict(pairs(r)...) for r in sweep_b_rows],
        :n_list => N_LIST,
        :geom => GEOM,
        :beta_grid => BETA_GRID,
        :T_c_at_h1 => T_C_AT_H1,
        :h => H_FIELD,
        :J => J_COUPLING,
        :r_D => R_D,
        :krylovdim => KRYLOVDIM,
        :howmany => HOWMANY,
        :wall_total => wall,
    ))
    println("\n[done] wall = $(round(wall; digits=1)) s")
    println("[done] sidecar: $bson_path")
end

main()
