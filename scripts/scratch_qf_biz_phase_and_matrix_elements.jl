#!/usr/bin/env julia
# scratch_qf_biz_phase_and_matrix_elements.jl  (qf-biz, qf-1jj follow-up)
#
# Two physics checks for the 2D TFIM ORDERED gap closing (qf-1jj draft
# `drafts/2d-tfim-ordered-vs-disordered.md`):
#
#   CHECK 1 — β sweep at n=6 (2×3 ladder), h=1.
#     For β_phys ∈ {0.1, 0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0}:
#       - Subtask B' dense-Gibbs diagnostics (cheap)
#       - L gap λ_L^phys via `krylov_spectral_gap` (canonical CKG smooth-Metro,
#         r_D = 7, kdim = 40, howmany = 4).
#     Goal: confirm β_phys = 2 is deep in the ordered phase; see how the L gap
#     varies as β crosses T_c(h=1) ≈ 2.07.
#
#   CHECK 2 — matrix elements at n=4, n=6 (β_phys = 2 ORDERED only).
#     Lowest-two eigenvectors |ψ_1⟩, |ψ_2⟩ of H_phys, then for each canonical
#     1/√(3n)-normalised single-site Pauli A_a ∈ {X_i, Y_i, Z_i}_{i=1..n}:
#       me_a = ⟨ψ_2|A_a|ψ_1⟩
#       M² = Σ_a |me_a|²  (the candidate "matrix-element bottleneck" factor)
#     Compare n=4 vs n=6 to decide if M² closes with n.
#     Also: γ_KMS(±ΔE_1^alg) so we can compare λ_L^alg / (γ × M²) across n.
#
# Skipped deliberately (qf-biz scope; user asked for the minimum that answers
# the physics questions):
#   - BohrDomain dense cross-check (canonical r_D=7 already validated at n=6 in qf-1jj)
#   - Trajectory τ_mix prediction (orthogonal — gap is the focus)
#   - Multiple Krylov dim / multiple seeds (clean Hamiltonian, deterministic)
#   - n=8 matrix elements (heavy, user said n=4,6 only)
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
println("[init] Julia threads = ", Threads.nthreads(), ", BLAS threads = ", BLAS.get_num_threads())

const OUTPUT_DIR = joinpath(@__DIR__, "output", "qf_biz_phase_and_matrix_elements")
mkpath(OUTPUT_DIR)

const J_COUPLING = 1.0
const H_FIELD    = 1.0
const KRYLOVDIM  = 40
const HOWMANY    = 4
const TAIL_C     = 8.0
const R_D        = 7
const T_C_AT_H1  = 2.07          # Hesselmann–Wessel 2016 (HW Fig. 4)

# CHECK 1: β sweep at fixed n=6 (2×3 ladder)
const N_CHECK1     = 6
const LX_CHECK1    = 2
const LY_CHECK1    = 3
# Spans deep-disordered (T/T_c ≈ 4) through T_c (β ≈ 0.48) to deep-ordered (T/T_c ≈ 0.16)
const BETA_GRID    = [0.1, 0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]

# CHECK 2: matrix elements at n=4 and n=6, ORDERED β_phys = 2 only
const N_CHECK2          = [(n=4, Lx=2, Ly=2), (n=6, Lx=2, Ly=3)]
const BETA_PHYS_CHECK2  = 2.0

# ---------------------------------------------------------------------------
# Build clean 2D TFIM (same recipe as qf-1jj phase-diagnostics / sweep, inline
# here so this script is self-contained).
# ---------------------------------------------------------------------------

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
        # extras useful for Check 2
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
        a = 0.0, s = 0.25,                  # canonical thesis convention (CLAUDE.md 2026-05-12)
        gaussian_parameters = (nothing, nothing),
        num_energy_bits_D = r_D, w0_D = w0_D, t0_D = t0_D,
        num_trotter_steps_per_t0 = 10,
        filter = nothing,
    )
end

# ---------------------------------------------------------------------------
# Subtask B' diagnostics (dense, cheap).
# ---------------------------------------------------------------------------

function z2_projector(n::Integer)
    P = Matrix{ComplexF64}(pad_term([X], n, 1))
    for i in 2:n
        P = P * Matrix(pad_term([X], n, i))
    end
    return P
end

"""Gibbs σ_β = e^{-β_alg H_alg}/Z built from cached eigendecomp of H_alg."""
function gibbs_alg(eigvals_alg::AbstractVector, eigvecs_alg::AbstractMatrix,
                   beta_alg_val::Real)
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

"""
    subtask_b_prime(raw, n, beta_phys_val) -> NamedTuple

PHYSICS CHECK: Gibbs-state phase membership diagnostics at one (β_phys, n)
cell on the clean 2D TFIM raw. Mirrors `scripts/scratch_2d_tfim_phase_diagnostics.jl`.
"""
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

    M = sum(Z_ops) ./ n
    M2 = M * M
    M4 = M2 * M2
    m1 = real(tr(M  * σβ))
    m2 = real(tr(M2 * σβ))
    m4 = real(tr(M4 * σβ))
    U4 = 1.0 - m4 / (3.0 * m2^2)

    S  = vn_entropy(σβ)
    Sn = S / log(d)
    eff_rank = exp(S)

    # ΔE_1^phys and Δ_4^phys: rescaling-frame-invariant up to multiplication by R.
    # raw.eigvals are H_alg eigenvalues (sorted ascending after `eigen`).
    Δ1_phys = (raw.eigvals[2] - raw.eigvals[1]) * raw.rescaling_factor
    Δ4_phys = length(raw.eigvals) >= 5 ?
        (raw.eigvals[5] - raw.eigvals[1]) * raw.rescaling_factor : NaN

    return (;
        n,
        beta_phys = beta_phys_val,
        beta_alg = β_alg_val,
        T_over_Tc = (1.0 / beta_phys_val) / T_C_AT_H1,
        w_plus, w_minus,
        m1, m2, m4, U4,
        S, Sn, eff_rank,
        Δ1_phys, Δ4_phys,
    )
end

# ---------------------------------------------------------------------------
# L gap via Krylov (matrix-free, canonical CKG smooth-Metro).
# ---------------------------------------------------------------------------

function compute_l_gap(raw, n::Integer, beta_phys_val::Real)
    ham = HamHam(raw; beta_phys=float(beta_phys_val))
    jumps = _build_jump_set(ham, n)
    cfg = build_ckg_energy_cfg(n, beta_phys_val, ham)
    t0 = time()
    res = krylov_spectral_gap(cfg, ham, jumps;
        krylovdim = KRYLOVDIM, howmany = HOWMANY, tol = 1e-10)
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
        cfg = cfg,
    )
end

# ---------------------------------------------------------------------------
# CHECK 1 — β sweep at n=6.
# ---------------------------------------------------------------------------

function run_check1()
    println("\n" * "="^100)
    println("CHECK 1 — β sweep at n=$(N_CHECK1) ($(LX_CHECK1)×$(LY_CHECK1) ladder), h=$(H_FIELD)")
    println("        T_c(h=1) ≈ $(T_C_AT_H1) (HW 2016)  ⇒  β_c ≈ $(round(1/T_C_AT_H1; sigdigits=3))")
    println("="^100)

    raw = build_clean_tfim_raw(LX_CHECK1, LY_CHECK1; J=J_COUPLING, h=H_FIELD)
    println("[check1] R = $(round(raw.rescaling_factor; digits=4)), shift = $(round(raw.shift; digits=4)), nu_min = $(round(raw.nu_min; sigdigits=4))")

    @printf("\n%-8s %-7s %-8s %-9s %-8s %-9s %-9s %-9s %-9s %-11s %-11s %-9s %-7s\n",
        "β_phys", "T/T_c", "β_alg", "<M_z²>", "U_4", "S/log d", "eff_rank",
        "w_+", "w_-", "Δ_1(phys)", "λ_L(phys)", "λ_L/ΔE_1", "wall")
    println("-"^140)

    rows = NamedTuple[]
    for β_phys_val in BETA_GRID
        diag = subtask_b_prime(raw, N_CHECK1; beta_phys_val=β_phys_val)
        gap = compute_l_gap(raw, N_CHECK1, β_phys_val)
        ratio = gap.gap_phys / diag.Δ1_phys

        @printf("%-8.3g %-7.3g %-8.3g %-9.4f %-8.4f %-9.4f %-9.3g %-9.4f %-9.4f %-11.4e %-11.4e %-9.3g %-7.1f\n",
            β_phys_val, diag.T_over_Tc, diag.beta_alg, diag.m2, diag.U4, diag.Sn,
            diag.eff_rank, diag.w_plus, diag.w_minus,
            diag.Δ1_phys, gap.gap_phys, ratio, gap.wall)

        push!(rows, merge(diag, (
            gap_alg = gap.gap_alg,
            gap_phys = gap.gap_phys,
            gap_eigenvalues = gap.eigenvalues,
            gap_converged = gap.converged,
            matvec = gap.matvec,
            wall_gap = gap.wall,
            ratio_gap_dE1 = ratio,
        )))
    end

    bson_path = joinpath(OUTPUT_DIR, "check1_beta_sweep_n6.bson")
    BSON.bson(bson_path, Dict(
        :rows => [Dict(pairs(r)...) for r in rows],
        :n => N_CHECK1, :Lx => LX_CHECK1, :Ly => LY_CHECK1,
        :h => H_FIELD, :J => J_COUPLING,
        :beta_grid => BETA_GRID, :T_c_at_h1 => T_C_AT_H1,
    ))
    println("\n[check1] sidecar: $bson_path")
    return rows
end

# ---------------------------------------------------------------------------
# CHECK 2 — matrix elements between the Z₂ doublet at n=4 and n=6.
# ---------------------------------------------------------------------------

"""
    matrix_element_diagnostic(raw, n, beta_phys_val) -> NamedTuple

Compute |⟨ψ_2|A_a|ψ_1⟩|² for every canonical jump operator A_a in the
1/√(3n)-normalised 3n-Pauli set. Also report bare-Pauli matrix elements
(without the normalisation factor), per-Pauli-type breakdown, Z₂ parities
of |ψ_1⟩, |ψ_2⟩, and γ_KMS(±ΔE_1^alg) at the canonical CKG config.
"""
function matrix_element_diagnostic(raw, n::Integer, beta_phys_val::Real)
    d = 2^n

    # Dense H_phys eigendecomposition (used for matrix elements at PHYS frame).
    F = eigen(Hermitian(raw.H_phys_dense))
    perm = sortperm(real(F.values))
    eigvals_phys_sorted = real(F.values)[perm]
    eigvecs_phys_sorted = F.vectors[:, perm]

    psi_1 = eigvecs_phys_sorted[:, 1]
    psi_2 = eigvecs_phys_sorted[:, 2]
    ΔE_1_phys = eigvals_phys_sorted[2] - eigvals_phys_sorted[1]

    # PHYSICS CHECK — Z₂ parity sanity:
    # |ψ_1⟩ is the symmetric (P=+1) ground state, |ψ_2⟩ is the
    # antisymmetric (P=-1) tunnelling partner. Verify ⟨ψ_i|P|ψ_i⟩ = ±1.
    P = z2_projector(n)
    psi_1_z2 = real(psi_1' * P * psi_1)
    psi_2_z2 = real(psi_2' * P * psi_2)
    @printf("[check2 n=%d] Z₂ parities: ⟨ψ_1|P|ψ_1⟩ = %+.6f,  ⟨ψ_2|P|ψ_2⟩ = %+.6f\n",
        n, psi_1_z2, psi_2_z2)

    # Canonical jump set (matches _build_jump_set in src/lindblad_action.jl:1025).
    # 3n single-site Paulis, each normalised by 1/√(3n).
    paulis = ([X], [Y], [Z])
    pauli_labels = ["X", "Y", "Z"]
    num_jumps = length(paulis) * n
    jump_norm_sq = float(num_jumps)   # |1/√(3n) · σ_α|² scales |me|² by 1/(3n)

    me2_norm = Float64[]              # |⟨2|A_a|1⟩|² WITH 1/(3n) normalisation
    me2_bare = Float64[]              # |⟨2|σ_α^(i)|1⟩|² WITHOUT normalisation
    me_labels = String[]
    me2_by_type = Dict("X"=>0.0, "Y"=>0.0, "Z"=>0.0)
    for (pi, pauli) in enumerate(paulis)
        for site in 1:n
            A_bare = Matrix{ComplexF64}(pad_term(pauli, n, site))
            me = psi_2' * A_bare * psi_1
            bare2 = abs2(me)
            push!(me2_bare, bare2)
            push!(me2_norm, bare2 / jump_norm_sq)
            push!(me_labels, "$(pauli_labels[pi])_$site")
            me2_by_type[pauli_labels[pi]] += bare2
        end
    end

    M2_norm = sum(me2_norm)            # the "matrix-element bottleneck" candidate
    M2_bare = sum(me2_bare)

    # γ at this Bohr frequency (alg frame, canonical CKG smooth-Metro)
    ΔE_1_alg = ΔE_1_phys / raw.rescaling_factor
    ham = HamHam(raw; beta_phys=float(beta_phys_val))
    cfg = build_ckg_energy_cfg(n, beta_phys_val, ham)
    γ_plus  = pick_transition(cfg, ΔE_1_alg)
    γ_minus = pick_transition(cfg, -ΔE_1_alg)
    γ_zero  = pick_transition(cfg, 0.0)

    return (;
        n,
        beta_phys = float(beta_phys_val),
        eigvals_phys_low5 = eigvals_phys_sorted[1:min(5, end)],
        ΔE_1_phys, ΔE_1_alg,
        psi_1_z2, psi_2_z2,
        M2_norm, M2_bare,
        me2_norm, me2_bare, me_labels, me2_by_type,
        γ_plus, γ_minus, γ_zero,
        rescaling_factor = raw.rescaling_factor,
        beta_alg = cfg.beta, sigma = cfg.sigma,
    )
end

function run_check2()
    println("\n" * "="^100)
    println("CHECK 2 — matrix elements at n ∈ {4, 6}, ORDERED β_phys = $(BETA_PHYS_CHECK2)")
    println("="^100)

    diags = NamedTuple[]
    gaps  = Dict{Int, NamedTuple}()

    for cell in N_CHECK2
        println("\n--- n = $(cell.n)  ($(cell.Lx)×$(cell.Ly) ladder) ---")
        raw  = build_clean_tfim_raw(cell.Lx, cell.Ly; J=J_COUPLING, h=H_FIELD)
        diag = matrix_element_diagnostic(raw, cell.n, BETA_PHYS_CHECK2)
        gap  = compute_l_gap(raw, cell.n, BETA_PHYS_CHECK2)
        gaps[cell.n] = gap

        @printf("ΔE_1^phys = %.4e    ΔE_1^alg = %.4e    β_alg = %.4g    σ_alg = %.4g\n",
            diag.ΔE_1_phys, diag.ΔE_1_alg, diag.beta_alg, diag.sigma)
        @printf("Lowest 5 E_phys = [%.6e, %.6e, %.6e, %.6e, %.6e]\n",
            diag.eigvals_phys_low5...)
        @printf("M²_norm (sum, 1/(3n)-normalised) = %.4e    M²_bare (sum, raw Paulis) = %.4e\n",
            diag.M2_norm, diag.M2_bare)
        @printf("By Pauli type (bare):  |X|²_sum = %.4e   |Y|²_sum = %.4e   |Z|²_sum = %.4e\n",
            diag.me2_by_type["X"], diag.me2_by_type["Y"], diag.me2_by_type["Z"])
        @printf("γ(+ΔE_1) = %.4f    γ(-ΔE_1) = %.4f    γ(0) = %.4f\n",
            diag.γ_plus, diag.γ_minus, diag.γ_zero)
        @printf("λ_L^alg  = %.4e    λ_L^phys = %.4e    matvecs = %d    converged = %s\n",
            gap.gap_alg, gap.gap_phys, gap.matvec, gap.converged)

        # Per-Pauli matrix-element table (top 6 contributors by |me|²_bare)
        idx_sorted = sortperm(diag.me2_bare; rev=true)[1:min(6, end)]
        println("Top contributors to M²_bare:")
        for i in idx_sorted
            @printf("  %-6s : |me|² = %.4e\n", diag.me_labels[i], diag.me2_bare[i])
        end

        push!(diags, diag)
    end

    # --- Comparison n=4 vs n=6 -------------------------------------------
    n4 = diags[1]; n6 = diags[2]
    println("\n" * "="^100)
    println("MATRIX-ELEMENT COMPARISON — n=4 → n=6 at β_phys = $(BETA_PHYS_CHECK2)")
    println("="^100)
    @printf("ΔE_1^phys:       n=4: %.4e    n=6: %.4e    ratio n6/n4: %.6f\n",
        n4.ΔE_1_phys, n6.ΔE_1_phys, n6.ΔE_1_phys / n4.ΔE_1_phys)
    @printf("M²_norm:         n=4: %.4e    n=6: %.4e    ratio n6/n4: %.6f\n",
        n4.M2_norm, n6.M2_norm, n6.M2_norm / n4.M2_norm)
    @printf("M²_bare:         n=4: %.4e    n=6: %.4e    ratio n6/n4: %.6f\n",
        n4.M2_bare, n6.M2_bare, n6.M2_bare / n4.M2_bare)
    @printf("|Y|² sum:        n=4: %.4e    n=6: %.4e    ratio n6/n4: %.6f\n",
        n4.me2_by_type["Y"], n6.me2_by_type["Y"],
        n6.me2_by_type["Y"] / max(n4.me2_by_type["Y"], 1e-30))
    @printf("|Z|² sum:        n=4: %.4e    n=6: %.4e    ratio n6/n4: %.6f\n",
        n4.me2_by_type["Z"], n6.me2_by_type["Z"],
        n6.me2_by_type["Z"] / max(n4.me2_by_type["Z"], 1e-30))
    @printf("λ_L^phys:        n=4: %.4e    n=6: %.4e    ratio n6/n4: %.6f\n",
        gaps[4].gap_phys, gaps[6].gap_phys, gaps[6].gap_phys / gaps[4].gap_phys)
    @printf("λ_L^alg:         n=4: %.4e    n=6: %.4e    ratio n6/n4: %.6f\n",
        gaps[4].gap_alg, gaps[6].gap_alg, gaps[6].gap_alg / gaps[4].gap_alg)

    # Attribution: if λ_L ≈ M² × γ × const, then
    #   λ_L(n6) / λ_L(n4) ≈ (M²(n6)/M²(n4)) · (γ(n6)/γ(n4))
    # If M² ratio is O(1) but λ_L ratio is exponentially small, the matrix-element
    # bottleneck story FAILS and the L gap closing has a different mechanism.
    γ_n4 = max(n4.γ_plus, n4.γ_minus)
    γ_n6 = max(n6.γ_plus, n6.γ_minus)
    γ_ratio = γ_n6 / γ_n4
    M2_ratio = n6.M2_norm / n4.M2_norm
    λ_ratio = gaps[6].gap_alg / gaps[4].gap_alg
    println("\nAttribution of L gap ratio (alg frame, n=4 → n=6):")
    @printf("  γ_KMS(±ΔE_1) ratio:    %.4e\n", γ_ratio)
    @printf("  M²_norm ratio:         %.4e\n", M2_ratio)
    @printf("  λ_L^alg ratio:         %.4e\n", λ_ratio)
    @printf("  Predicted (γ × M²):    %.4e\n", γ_ratio * M2_ratio)
    @printf("  λ_actual / Predicted:  %.4f  (should be ~1 if M²×γ explains the closing)\n",
        λ_ratio / (γ_ratio * M2_ratio))

    bson_path = joinpath(OUTPUT_DIR, "check2_matrix_elements.bson")
    BSON.bson(bson_path, Dict(
        :diags => [Dict(pairs(d)...) for d in diags],
        :gaps  => Dict(k => Dict(:gap_alg=>v.gap_alg, :gap_phys=>v.gap_phys,
                                  :eigenvalues=>v.eigenvalues,
                                  :converged=>v.converged, :matvec=>v.matvec,
                                  :rescaling_factor=>v.rescaling_factor,
                                  :beta_alg=>v.beta_alg, :sigma=>v.sigma)
                       for (k, v) in gaps),
        :beta_phys => BETA_PHYS_CHECK2,
    ))
    println("\n[check2] sidecar: $bson_path")
    return diags, gaps
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

function main()
    println("[main] start  $(now())")
    println("[main] J = $(J_COUPLING), h = $(H_FIELD)  (clean 2D TFIM, no disorder)")
    println("[main] CKG canonical: smooth-Metro s=0.25, a=0, r_D=$(R_D), kdim=$(KRYLOVDIM)")
    t0 = time()
    rows1 = run_check1()
    diags2, gaps2 = run_check2()
    wall = time() - t0
    println("\n[main] done   $(now())   total wall = $(round(wall; digits=1)) s")
end

main()
