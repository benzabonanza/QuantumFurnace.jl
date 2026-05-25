#!/usr/bin/env julia
# scratch_qf_c9g_slow_mode_n4.jl  (qf-c9g — slow-mode operator-space character)
#
# Question: at n=4 ORDERED β_phys=2, qf-biz Follow-up A established that the
# slow Lindbladian mode R_2 lives 89% in the doublet × bulk off-diagonal block
# (the surface-tension signature). Does this PERSIST across moderately-cold β
# where Gibbs is NOT just the doublet (β=0.5, eff_rank ≈ 5 at n=4)? If yes,
# the mechanism is the same surface-tension story all the way down to β_c;
# if not, the picture must be revised.
#
# Method: dense `construct_lindbladian` at n=4 (d=16, d²=256 — trivial), then
# `eigen(L_dense)`, sort by |Re λ|, and decompose R_2 in the energy eigenbasis
# (block fractions: doublet × doublet, doublet × bulk, bulk × bulk + Frobenius
# overlaps with M_z and various doublet-coherence templates).
#
# 5 cells: β_phys ∈ {0.1, 0.25, 0.5, 1.0, 2.0} at n=4 (Lx=Ly=2). β=2 case
# replicates qf-biz Follow-up A exactly as a sanity check.
#
# Wall budget: ~5 s per cell × 5 = ~30 s.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using QuantumFurnace
using QuantumFurnace: _construct_2d_heisenberg_base, _construct_disordering_terms,
                     _rescaling_and_shift_factors, pad_term, _build_jump_set
using LinearAlgebra
using BSON
using Printf
using Dates

BLAS.set_num_threads(1)
println("[init] $(Dates.now())  Julia threads = $(Threads.nthreads())  BLAS threads = $(BLAS.get_num_threads())")

const OUTPUT_DIR = joinpath(@__DIR__, "output", "qf_c9g_ordered_gap_mechanism")
mkpath(OUTPUT_DIR)

const J_COUPLING = 1.0
const H_FIELD    = 1.0
const TAIL_C     = 8.0
const R_D        = 7
const N          = 4
const LX, LY     = 2, 2
const BETA_GRID  = [0.10, 0.25, 0.50, 1.0, 2.0]
const T_C_AT_H1  = 2.07

# ============================================================================
# Same TFIM raw + config builders as the main sweep script
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
    return (
        matrix              = rescaled,
        terms               = [[Z, Z], [X]],
        base_coeffs         = [-J / rescaling_factor, -h / rescaling_factor],
        disordering_terms   = nothing,
        disordering_coeffs  = nothing,
        eigvals             = rescaled_eigvals,
        eigvecs             = rescaled_eigvecs,
        nu_min              = minimum(diff(rescaled_eigvals)),
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

function z2_projector(n::Integer)
    P = Matrix{ComplexF64}(pad_term([X], n, 1))
    for i in 2:n
        P = P * Matrix(pad_term([X], n, i))
    end
    return P
end

# ============================================================================
# Slow-mode decomposition
# ============================================================================

"""
    decompose_slow_mode(L_dense, raw, n) -> NamedTuple

Given a dense Lindbladian, find the slowest-decaying non-zero eigenmode R_2
(reshape vec → d×d), and decompose it in the energy eigenbasis. Identifies
how much of |R_2|² lives in doublet × doublet, doublet × bulk, bulk × bulk
blocks, plus Z₂ parity character and order-parameter overlap.

The 'doublet' is defined as the two lowest energy eigenstates {ψ_1, ψ_2}.
"""
function decompose_slow_mode(L_dense::AbstractMatrix, raw, n::Integer;
                              beta_phys_val::Real)
    d = 2^n
    F = eigen(L_dense)
    # sort eigenvalues by |Re λ| ascending (so λ_1 ≈ 0 fixed point, λ_2 = slow mode)
    perm = sortperm([abs(real(λ)) for λ in F.values])
    eigvals_sorted = F.values[perm]
    eigvecs_sorted = F.vectors[:, perm]

    λ_1 = eigvals_sorted[1]
    λ_2 = eigvals_sorted[2]

    # Reshape R_2 from vec form (d² entries) back to d × d
    R2_vec = eigvecs_sorted[:, 2]
    R_2 = reshape(R2_vec, (d, d))
    R_2_norm = sqrt(real(tr(R_2' * R_2)))

    # Energy-eigenbasis representation: R_2 in {ψ_k}-basis
    U = raw.eigvecs    # columns are ψ_k (algorithm-frame H eigenvectors)
    R_2_ebasis = U' * R_2 * U
    abs2_block = abs2.(R_2_ebasis)

    # block fractions in energy eigenbasis
    mass_doublet = sum(abs2_block[1:2, 1:2])
    mass_db_bulk = sum(abs2_block[1:2, 3:end]) + sum(abs2_block[3:end, 1:2])
    mass_bulk    = sum(abs2_block[3:end, 3:end])
    total_mass   = sum(abs2_block)

    # Z₂ parity character of R_2:  tr(R_2† P R_2 P) / ‖R_2‖²
    P = z2_projector(n)
    z2_parity = real(tr(R_2' * P * R_2 * P)) / max(real(tr(R_2' * R_2)), eps())

    # M_z overlap (Frobenius): ⟨R_2, M_z⟩ / (‖R_2‖·‖M_z‖)
    Z_ops = [Matrix{ComplexF64}(pad_term([Z], n, i)) for i in 1:n]
    M_z = sum(Z_ops)
    Mz_norm = sqrt(real(tr(M_z' * M_z)))
    overlap_Mz = abs(tr(R_2' * M_z)) / max(R_2_norm * Mz_norm, eps())

    # Within-doublet coherence overlap: |⟨R_2, |ψ_1⟩⟨ψ_2| + h.c.⟩| (symm)
    psi_1 = U[:, 1]; psi_2 = U[:, 2]
    coh_sym = (psi_1 * psi_2' + psi_2 * psi_1') ./ sqrt(2)
    coh_sym_norm = sqrt(real(tr(coh_sym' * coh_sym)))
    overlap_doublet_coh_sym = abs(tr(R_2' * coh_sym)) / max(R_2_norm * coh_sym_norm, eps())

    coh_anti = (im) .* (psi_1 * psi_2' - psi_2 * psi_1') ./ sqrt(2)
    coh_anti_norm = sqrt(real(tr(coh_anti' * coh_anti)))
    overlap_doublet_coh_anti = abs(tr(R_2' * coh_anti)) / max(R_2_norm * coh_anti_norm, eps())

    # Doublet imbalance (population difference): |ψ_1⟩⟨ψ_1| − |ψ_2⟩⟨ψ_2|
    imbalance = (psi_1 * psi_1' - psi_2 * psi_2') ./ sqrt(2)
    imbalance_norm = sqrt(real(tr(imbalance' * imbalance)))
    overlap_imbalance = abs(tr(R_2' * imbalance)) / max(R_2_norm * imbalance_norm, eps())

    # Top off-diagonal entries (largest |R_2_ebasis[i,j]| with i ≠ j)
    top_entries = NamedTuple[]
    flat_idx = sortperm(vec(abs2_block); rev = true)
    seen = 0
    for k in flat_idx
        i, j = ind2sub_ij(d, k)
        if i != j
            push!(top_entries, (
                i = i, j = j,
                E_i = raw.eigvals[i] * raw.rescaling_factor,
                E_j = raw.eigvals[j] * raw.rescaling_factor,
                ω_phys = (raw.eigvals[i] - raw.eigvals[j]) * raw.rescaling_factor,
                abs2 = abs2_block[i, j],
                amp = R_2_ebasis[i, j],
            ))
            seen += 1
            seen >= 8 && break
        end
    end

    return (;
        n,
        beta_phys = beta_phys_val,
        T_over_Tc = (1.0 / beta_phys_val) / T_C_AT_H1,
        d = d,
        λ_1_alg = λ_1,
        λ_2_alg = λ_2,
        λ_2_phys = λ_2 * raw.rescaling_factor,
        R_2_frob_norm = R_2_norm,
        mass_doublet_doublet = mass_doublet / total_mass,
        mass_doublet_bulk    = mass_db_bulk / total_mass,
        mass_bulk_bulk       = mass_bulk    / total_mass,
        z2_parity,
        overlap_Mz,
        overlap_doublet_coh_sym,
        overlap_doublet_coh_anti,
        overlap_imbalance,
        top_entries,
    )
end

# linear index k → (i,j) for d × d Fortran-order reshape (Julia column-major)
ind2sub_ij(d::Integer, k::Integer) = (mod1(k, d), div(k - 1, d) + 1)

# ============================================================================
# Main
# ============================================================================

function main()
    println("\n" * "="^120)
    println("SWEEP C — slow-mode operator-space character at n=$N (h=$H_FIELD, J=$J_COUPLING)")
    println("="^120)

    raw = build_clean_tfim_raw(LX, LY; J=J_COUPLING, h=H_FIELD)
    @printf("[fixture] n=%d  R = %.4f  shift = %+.4f  nu_min = %.3g\n",
        N, raw.rescaling_factor, raw.shift, raw.nu_min)

    @printf("\n%-7s %-7s %-7s %-10s %-12s %-14s %-12s %-9s %-9s %-10s\n",
        "β_phys", "T/T_c", "phase", "|λ_2|^phys",
        "M(d×d)/|R₂|²", "M(d×bulk)/|R₂|²", "M(b×b)/|R₂|²", "Z₂par", "⟨R₂,M_z⟩", "⟨R₂,σ⟩(d-coh)")
    println("-"^140)

    rows = NamedTuple[]
    for β_phys_val in BETA_GRID
        ham = HamHam(raw; beta_phys=float(β_phys_val))
        jumps = _build_jump_set(ham, N)
        cfg = build_ckg_energy_cfg(N, β_phys_val, ham)
        t0 = time()
        L_dense = construct_lindbladian(jumps, cfg, ham; include_coherent=true)
        decomp = decompose_slow_mode(L_dense, raw, N; beta_phys_val=β_phys_val)
        wall = time() - t0

        phase_str = decomp.T_over_Tc < 1.0 ? "ORD" : "DIS"
        @printf("%-7.3g %-7.3g %-7s %-10.4e %-12.4f %-14.4f %-12.4f %-9.3g %-9.3g %-10.3g\n",
            decomp.beta_phys, decomp.T_over_Tc, phase_str,
            abs(decomp.λ_2_phys),
            decomp.mass_doublet_doublet, decomp.mass_doublet_bulk, decomp.mass_bulk_bulk,
            decomp.z2_parity, decomp.overlap_Mz, decomp.overlap_doublet_coh_sym)

        push!(rows, merge(decomp, (wall = wall,)))
    end

    println("\nTop-three off-diagonal entries per cell (where R_2 lives in energy eigenbasis):")
    for r in rows
        @printf("  β=%-5.3g:  top 3:  ", r.beta_phys)
        for k in 1:min(3, length(r.top_entries))
            e = r.top_entries[k]
            @printf("(%d,%d)[ω=%+.3f, |amp|²=%.3f]  ", e.i, e.j, e.ω_phys, e.abs2 / r.R_2_frob_norm^2)
        end
        println()
    end

    bson_path = joinpath(OUTPUT_DIR, "qf_c9g_slow_mode_n4.bson")
    BSON.bson(bson_path, Dict(
        :rows => [Dict(pairs(r)...) for r in rows],
        :n => N, :Lx => LX, :Ly => LY,
        :beta_grid => BETA_GRID,
        :T_c_at_h1 => T_C_AT_H1,
        :J => J_COUPLING, :h => H_FIELD,
    ))
    println("\n[done] sidecar: $bson_path")
end

main()
