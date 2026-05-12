#!/usr/bin/env julia
#
# Operator-norm of the coherent term — does DLL G overstep the LCU/block-
# encoding norm-1 budget?
#
# Motivation
# ----------
# Ding et al. 2024 Remark 23 (paper p. 24) and Eq. 3.31 state that DLL's
# coherent term G has a block-encoding normalisation constant that grows
# linearly in β. Their Table 1 carries this through to a β² overall total
# cost vs CKG's β. We test the underlying claim numerically:
#
#   - Build the DLL Gaussian G and DLL Metropolis G via `dll_coherent_op_bohr`
#     across β ∈ {1, 2, 5, 10, 15, 20, 30}, n ∈ {3, 4, 5}.
#   - For comparison, build CKG smooth-Metropolis B via `B_bohr`.
#   - Report `opnorm(G)` for each. An LCU implementation block-encodes
#     `G / ‖G‖`, so any ‖G‖ > 1 means a normalisation factor must be
#     introduced (and contributes linearly to the simulation cost
#     normalisation in Lindblad-LCU schemes).
#
# Construction note
# -----------------
# B_bohr returns the CKG coherent term scaled by `gamma_norm_factor` in
# the production path (`_precompute_coherent_B` does the `rmul!`). Post-
# qf-etx.2 `gamma_norm_factor = 1 / pick_gamma_sup(config) = 1.0` for every
# standard family, so the scaled and unscaled forms now coincide.
#
# PHYSICS CHECK: jump set is the standard single-site Paulis (X, Y, Z) on
# each site with the conventional `1/√(3n)` normalisation, matching
# `sweep_mixing_times`. The coherent term G scales with the sum over
# couplings, so this normalisation matters for comparison consistency.
#
# PHYSICS CHECK: rescaled Hamiltonians have `‖H‖ ≈ 0.45` and Bohr freqs
# in [-0.45, 0.45], so β·max|ν| ranges from 0.45 (β=1) to 13.5 (β=30).
# The transition between linear-in-β and saturated regimes is expected
# around β ~ 1/max|ν| ~ 2.

using Printf
using LinearAlgebra
using QuantumFurnace
include(joinpath(@__DIR__, "..", "test", "test_helpers.jl"))

# ── Sweep grid ────────────────────────────────────────────────────────────────
const NS    = [3, 4, 5]
const BETAS = [1.0, 2.0, 5.0, 10.0, 15.0, 20.0, 30.0]
const SIGMA_CKG_FACTOR = 1.0   # σ = SIGMA_CKG_FACTOR / β   (sweep convention)
const A_CKG = 0.0
const S_CKG = 0.25
const S_DLL = 5.0  # bump support — must be ≥ 2·max|ν_BH| (~0.9)

# ── Helpers (small, inline, no library additions) ─────────────────────────────

# Recreate the production jump set from `_build_jump_set` (private in
# lindblad_action.jl) without importing it.
function build_paulis_jumps(ham::HamHam, n::Int)
    paulis = ([X], [Y], [Z])
    num_jumps = length(paulis) * n
    jump_norm = sqrt(num_jumps)
    jumps = JumpOp[]
    for pauli in paulis, site in 1:n
        op = Matrix(pad_term(pauli, n, site)) ./ jump_norm
        op_eb = ham.eigvecs' * op * ham.eigvecs
        push!(jumps, JumpOp(op, op_eb,
                            op == transpose(op), op == op'))
    end
    return jumps
end

# Build a minimal BohrDomain Config for B_bohr at given (n, β).
function ckg_bohr_config(n::Int, β::Float64)
    return Config(
        sim = Lindbladian(),
        domain = BohrDomain(),
        construction = KMS(),
        num_qubits = n,
        with_linear_combination = true,
        beta = β,
        sigma = SIGMA_CKG_FACTOR / β,   # matches sweep_mixing_times convention
        a = A_CKG,
        s = S_CKG,
        filter = nothing,
    )
end

# ── Main sweep ────────────────────────────────────────────────────────────────
println("="^96)
println("Coherent-term op-norm check — does DLL G break the ‖·‖ ≤ 1 LCU budget?")
println("="^96)
@printf("%-3s %-5s | %-10s %-12s %-12s | %-10s %-10s\n",
        "n", "β", "‖B_CKG‖", "‖G_DLLgauss‖", "‖G_DLLmetro‖",
        "‖H‖", "max|β·ν|")
println("-"^96)

# Cache results for the slope check at the end.
results = Dict{Tuple{Int,Float64}, NamedTuple}()

for n in NS
    ham_path = joinpath(@__DIR__, "..", "hamiltonians",
                         "heis_xxx_zzdisordered_periodic_n$(n).bson")
    isfile(ham_path) || (@warn "missing fixture, skipping" n; continue)

    # β-independent: ham eigvals + jumps. (Loaded once per n.)
    ham_ref = _load_test_hamiltonian(ham_path, first(BETAS))
    jumps   = build_paulis_jumps(ham_ref, n)
    H_norm  = maximum(abs, ham_ref.eigvals)              # ‖H‖ for rescaled fixture
    ν_max   = maximum(abs, keys(ham_ref.bohr_dict))

    for β in BETAS
        # Re-load to get the right Gibbs state for the β (Bohr dict / eigvals
        # are β-independent but `_load_test_hamiltonian` also stamps `gibbs`).
        ham = _load_test_hamiltonian(ham_path, β)

        config_ckg = ckg_bohr_config(n, β)
        B_ckg      = B_bohr(ham, jumps, config_ckg)
        G_dll_g    = dll_coherent_op_bohr(jumps, ham, DLLGaussianFilter(β),         β)
        G_dll_m    = dll_coherent_op_bohr(jumps, ham, DLLMetropolisFilter(β; S=S_DLL), β)

        norm_B   = opnorm(Matrix(B_ckg))
        norm_G_g = opnorm(Matrix(G_dll_g))
        norm_G_m = opnorm(Matrix(G_dll_m))

        @printf("%-3d %-5.1f | %-10.4e %-12.4e %-12.4e | %-10.4f %-10.2f\n",
                n, β, norm_B, norm_G_g, norm_G_m, H_norm, β * ν_max)

        results[(n, β)] = (; norm_B, norm_G_g, norm_G_m, H_norm, ν_max)
    end
    println("-"^96)
end

# ── Slope diagnostic: does ‖G_DLL‖ scale linearly with β at large β? ──────────
# Linear regression of log10(‖G‖) on log10(β) over the high-β tail.
function slope_loglog(βs::Vector{Float64}, ys::Vector{Float64})
    x = log10.(βs); y = log10.(ys)
    n = length(x); x̄ = sum(x)/n; ȳ = sum(y)/n
    return sum((x .- x̄) .* (y .- ȳ)) / sum((x .- x̄).^2)
end

println()
println("="^96)
println("β-scaling slope (log10–log10 fit) over high-β tail β ∈ {5, 10, 15, 20, 30}")
println("Expected per Ding et al. 2024 Remark 23: slope ≈ 1.0 for DLL G")
println("="^96)
high_β = [5.0, 10.0, 15.0, 20.0, 30.0]
@printf("%-3s | %-12s %-12s %-12s\n", "n", "slope_B_CKG", "slope_G_dllG", "slope_G_dllM")
println("-"^60)
for n in NS
    haskey(results, (n, first(high_β))) || continue
    Bs = [results[(n, β)].norm_B    for β in high_β]
    Ggs = [results[(n, β)].norm_G_g for β in high_β]
    Gms = [results[(n, β)].norm_G_m for β in high_β]
    @printf("%-3d | %-12.4f %-12.4f %-12.4f\n",
            n, slope_loglog(high_β, Bs),
               slope_loglog(high_β, Ggs),
               slope_loglog(high_β, Gms))
end

# ── Headline ──────────────────────────────────────────────────────────────────
println()
println("="^96)
println("Headline: does the DLL coherent term G overstep ‖G‖ ≤ 1?")
println("="^96)
for n in NS, β in BETAS
    r = get(results, (n, β), nothing)
    r === nothing && continue
    over_g = r.norm_G_g > 1
    over_m = r.norm_G_m > 1
    if over_g || over_m
        @printf("n=%d β=%-4.1f  ‖G_DLL_Gauss‖=%.3f%s  ‖G_DLL_Metro‖=%.3f%s\n",
                n, β,
                r.norm_G_g, over_g ? "  ← exceeds 1" : "",
                r.norm_G_m, over_m ? "  ← exceeds 1" : "")
    end
end

println("\nDone.")
