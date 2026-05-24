#!/usr/bin/env julia
# scratch_qf_biz_followup.jl   (qf-biz, follow-ups recommended by physics-checker)
#
# Two cheap follow-ups to `scratch_qf_biz_phase_and_matrix_elements.jl`:
#
#   FOLLOW-UP A — slow-mode eigenmode at n=4 ORDERED β=2 (dense L, d²=256).
#     Build dense `construct_lindbladian`, eigendecompose, identify the operator
#     structure of the slowest non-zero mode R_2. Compute:
#       - Z₂ parity ⟨R_2 | P R_2 P⟩ / ⟨R_2 | R_2⟩
#       - overlap with normalised "doublet coherence" (|ψ_1⟩⟨ψ_2| + h.c.)/√2
#       - overlap with normalised "doublet imbalance" (|ψ_1⟩⟨ψ_1| - |ψ_2⟩⟨ψ_2|)/√2
#       - overlap with magnetisation M_z (Frob-normalised, intensive)
#       - matrix R_2 in the energy eigenbasis (top-few significant entries)
#     Goal: figure out what the slow mode actually IS.
#
#   FOLLOW-UP B — kinky Metropolis (s=0, a=0) at n=6 β_phys=2 (single cell).
#     If the L gap under kinky-Metro is monotonically slower at colder β (vs
#     the smooth-Metro non-monotone behaviour seen in Check 1), the
#     non-monotonicity is a smooth-Metro filter artefact. Otherwise it's
#     physical.

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
const TAIL_C     = 8.0
const R_D        = 7
const KRYLOVDIM  = 40
const HOWMANY    = 4

function build_clean_tfim_raw(Lx, Ly; J=1.0, h=1.0)
    n = Lx * Ly
    H_bond  = _construct_2d_heisenberg_base(Lx, Ly, [[Z, Z]], [-J];
                                             periodic_x=true, periodic_y=true)
    H_field = _construct_disordering_terms([[X]], [fill(-h, n)], n)
    H_phys = Hermitian(Matrix(H_bond) + Matrix(H_field))
    R, s = _rescaling_and_shift_factors(H_phys)
    d = 2^n
    rescaled = Matrix(H_phys) ./ R .+ s * I(d)
    rev = rescaled
    rescaled_eigvals, rescaled_eigvecs = eigen(Hermitian(rev))
    nu_min = minimum(diff(rescaled_eigvals))
    return (
        matrix=rescaled, terms=[[Z,Z],[X]],
        base_coeffs=[-J/R, -h/R],
        disordering_terms=nothing, disordering_coeffs=nothing,
        eigvals=rescaled_eigvals, eigvecs=rescaled_eigvecs,
        nu_min=nu_min, shift=s, rescaling_factor=R, periodic=true,
        H_phys_dense=Matrix{ComplexF64}(H_phys),
    )
end

function z2_projector(n)
    P = Matrix{ComplexF64}(pad_term([X], n, 1))
    for i in 2:n
        P = P * Matrix(pad_term([X], n, i))
    end
    return P
end

function build_cfg(n, beta_phys_val, ham; s_val=0.25, a_val=0.0, r_D_val=R_D)
    β_alg_val = beta_alg(ham, float(beta_phys_val))
    σ = 1.0 / β_alg_val
    H_norm = maximum(abs, ham.eigvals)
    omega_range = 2.0 * (H_norm + TAIL_C * σ)
    w0_D = omega_range / 2.0^r_D_val
    t0_D = 2π / (2.0^r_D_val * w0_D)
    return Config(
        sim=Lindbladian(), domain=EnergyDomain(), construction=KMS(),
        num_qubits=n, with_linear_combination=true,
        beta=β_alg_val, beta_phys=float(beta_phys_val), sigma=σ,
        a=a_val, s=s_val,
        gaussian_parameters=(nothing, nothing),
        num_energy_bits_D=r_D_val, w0_D=w0_D, t0_D=t0_D,
        num_trotter_steps_per_t0=10, filter=nothing,
    )
end

# ---------------------------------------------------------------------------
# FOLLOW-UP A — dense L eigenmode at n=4 ORDERED β=2.
# ---------------------------------------------------------------------------

function followup_a()
    println("\n" * "="^100)
    println("FOLLOW-UP A — dense L eigenmode at n=4 (2×2 ladder), β_phys=2 ORDERED")
    println("="^100)

    n = 4
    raw = build_clean_tfim_raw(2, 2; J=J_COUPLING, h=H_FIELD)
    ham = HamHam(raw; beta_phys=2.0)
    jumps = _build_jump_set(ham, n)
    cfg = build_cfg(n, 2.0, ham)

    d = 2^n
    println("[A] d = $d, d² = $(d*d). Building dense L...")
    t0 = time()
    L_dense = construct_lindbladian(jumps, cfg, ham; include_coherent=true)
    println("[A] L_dense built in $(round(time()-t0; digits=3)) s, size = $(size(L_dense))")

    # Diagonalise. Sort by |Re|.
    t0 = time()
    F = eigen(L_dense)
    perm = sortperm(F.values; by = v -> abs(real(v)))
    vals = F.values[perm]
    vecs = F.vectors[:, perm]
    println("[A] eigen in $(round(time()-t0; digits=3)) s")

    @printf("[A] Lowest 6 |Re λ|: %s\n",
        join([@sprintf("%.4e", abs(real(v))) for v in vals[1:6]], "  "))

    # Slowest non-zero mode: λ_2 (skip λ_1 ≈ 0 which is the steady state)
    λ_2 = vals[2]
    R_2_vec = vecs[:, 2]
    R_2 = reshape(R_2_vec, d, d)
    # Frobenius-normalise
    R_2 ./= norm(R_2)

    @printf("[A] λ_2 = %.4e %+.4e·i    |Re λ_2| = %.4e (= L gap = %.4e from qf-1jj)\n",
        real(λ_2), imag(λ_2), abs(real(λ_2)), 2.6243e-3)

    # Z₂ parity of R_2: project under P R_2 P; ratio = ±1 if pure parity
    P = z2_projector(n)
    PR2P = P * R_2 * P
    z2_inner = real(tr(R_2' * PR2P))   # tr(R_2† P R_2 P)
    z2_self  = real(tr(R_2' * R_2))    # ‖R_2‖² = 1
    @printf("[A] Z₂ parity tr(R_2† P R_2 P) / ‖R_2‖² = %+.6f   (should be ±1)\n",
        z2_inner / z2_self)

    # H_phys eigenvectors for projections
    Fp = eigen(Hermitian(raw.H_phys_dense))
    perm_p = sortperm(real(Fp.values))
    psi = Fp.vectors[:, perm_p]
    psi_1 = psi[:, 1]
    psi_2 = psi[:, 2]

    # Build operator-basis projections (Frobenius-normalised d×d operators)
    A_coh_sym  = (psi_1 * psi_2' + psi_2 * psi_1') / sqrt(2.0)
    A_coh_anti = (1im * (psi_1 * psi_2' - psi_2 * psi_1')) / sqrt(2.0)
    A_imbalance = (psi_1 * psi_1' - psi_2 * psi_2') / sqrt(2.0)
    # Magnetisation M = (1/n) Σ_i σ_z^(i), Frob-normalised
    M = sum(Matrix{ComplexF64}(pad_term([Z], n, i)) for i in 1:n) ./ n
    M ./= norm(M)

    # Inner products with R_2 (Frobenius-Hilbert-Schmidt)
    function fro_inner(A, B)
        return tr(A' * B)
    end

    o_coh_sym  = fro_inner(A_coh_sym,  R_2)
    o_coh_anti = fro_inner(A_coh_anti, R_2)
    o_imb      = fro_inner(A_imbalance, R_2)
    o_M        = fro_inner(M,           R_2)

    @printf("[A] |⟨R_2 | (|ψ_1⟩⟨ψ_2|+h.c.)/√2⟩| = %.4f\n", abs(o_coh_sym))
    @printf("[A] |⟨R_2 | i(|ψ_1⟩⟨ψ_2|-h.c.)/√2⟩|= %.4f\n", abs(o_coh_anti))
    @printf("[A] |⟨R_2 | (|ψ_1⟩⟨ψ_1|-|ψ_2⟩⟨ψ_2|)/√2⟩|= %.4f\n", abs(o_imb))
    @printf("[A] |⟨R_2 | M_z (Frob-norm)⟩|         = %.4f\n", abs(o_M))

    # R_2 in the energy eigenbasis: show top entries
    R_2_eb = psi' * R_2 * psi
    @printf("[A] R_2 in energy eigenbasis — top 8 |entries|:\n")
    entries = [(abs(R_2_eb[i, j]), i, j) for i in 1:d, j in 1:d if i != j || true]
    entries_sorted = sort(vec(entries); by = x -> x[1], rev=true)
    for k in 1:min(8, length(entries_sorted))
        m, i, j = entries_sorted[k]
        @printf("    (%d,%d): %.4e\n", i, j, m)
    end

    # Projection: fraction of ‖R_2‖² that lives in
    # (a) doublet ⊗ doublet block (i,j ∈ {1,2})
    # (b) doublet ⊗ bulk     (i ∈ {1,2}, j > 2 OR i > 2, j ∈ {1,2})
    # (c) bulk ⊗ bulk        (i,j > 2)
    block_dd = 0.0; block_db = 0.0; block_bb = 0.0
    for i in 1:d, j in 1:d
        w = abs2(R_2_eb[i, j])
        i_in_dbl = (i <= 2)
        j_in_dbl = (j <= 2)
        if i_in_dbl && j_in_dbl
            block_dd += w
        elseif i_in_dbl ⊻ j_in_dbl
            block_db += w
        else
            block_bb += w
        end
    end
    @printf("[A] ‖R_2‖² mass by block: doublet×doublet = %.4f, doublet×bulk = %.4f, bulk×bulk = %.4f\n",
        block_dd, block_db, block_bb)

    return (;
        λ_2, R_2_eb,
        z2_parity = z2_inner / z2_self,
        o_coh_sym, o_coh_anti, o_imb, o_M,
        block_dd, block_db, block_bb,
    )
end

# ---------------------------------------------------------------------------
# FOLLOW-UP B — kinky Metropolis check at n=6 β_phys=2.
# ---------------------------------------------------------------------------

function followup_b()
    println("\n" * "="^100)
    println("FOLLOW-UP B — kinky Metro (s=0, a=0) at n=6 β_phys=2 vs smooth-Metro reference")
    println("="^100)

    n = 6
    raw = build_clean_tfim_raw(2, 3; J=J_COUPLING, h=H_FIELD)
    ham = HamHam(raw; beta_phys=2.0)
    jumps = _build_jump_set(ham, n)

    results = NamedTuple[]

    for β_val in [1.0, 1.5, 2.0, 2.5, 3.0]
        ham_b = HamHam(raw; beta_phys=β_val)
        jumps_b = _build_jump_set(ham_b, n)

        # Kinky (s=0, a=0)
        cfg_k = build_cfg(n, β_val, ham_b; s_val=0.0, a_val=0.0)
        t0 = time()
        res_k = krylov_spectral_gap(cfg_k, ham_b, jumps_b;
            krylovdim=KRYLOVDIM, howmany=HOWMANY, tol=1e-10)
        w_k = time() - t0

        @printf("β_phys=%.2g  kinky   λ_phys = %.4e   matvec=%d   wall=%.1fs\n",
            β_val, res_k.spectral_gap * ham_b.rescaling_factor, res_k.matvec_count, w_k)
        push!(results, (β_phys=β_val, kind=:kinky,
            gap_phys=res_k.spectral_gap * ham_b.rescaling_factor,
            gap_alg=res_k.spectral_gap, matvec=res_k.matvec_count, wall=w_k))
    end

    println("\n[B] Comparison kinky vs smooth (from Check 1):")
    @printf("%-9s  %-13s  %-13s  %-7s\n", "β_phys", "λ_phys kinky", "λ_phys smooth", "ratio")
    smooth = Dict(1.0=>1.2371e-2, 1.5=>2.6454e-3, 2.0=>3.0648e-3, 2.5=>4.1815e-3, 3.0=>5.5741e-3)
    for r in results
        sm = smooth[r.β_phys]
        @printf("%-9.2g  %-13.4e  %-13.4e  %-7.3g\n", r.β_phys, r.gap_phys, sm, r.gap_phys/sm)
    end

    return results
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

function main()
    println("[main] start  $(now())")
    t0 = time()
    res_a = followup_a()
    res_b = followup_b()
    wall = time() - t0
    println("\n[main] done   $(now())   total wall = $(round(wall; digits=1)) s")

    bson_path = joinpath(OUTPUT_DIR, "followup_a_b.bson")
    BSON.bson(bson_path, Dict(
        :followup_a => Dict(
            :lambda_2 => res_a.λ_2,
            :z2_parity => res_a.z2_parity,
            :o_coh_sym => abs(res_a.o_coh_sym),
            :o_coh_anti => abs(res_a.o_coh_anti),
            :o_imb => abs(res_a.o_imb),
            :o_M => abs(res_a.o_M),
            :block_dd => res_a.block_dd,
            :block_db => res_a.block_db,
            :block_bb => res_a.block_bb,
        ),
        :followup_b => [Dict(pairs(r)...) for r in res_b],
    ))
    println("\n[main] sidecar: $bson_path")
end

main()
