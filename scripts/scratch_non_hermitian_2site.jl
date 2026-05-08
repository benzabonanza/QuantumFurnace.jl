#!/usr/bin/env julia
# qf-bm1 follow-up: 2-site non-Hermitian fixture test.
#
# Per physics-checker review (.claude-memory/qf_bm1_physics_review.md A.3),
# the vanishing condition for the Q2 contact term is
#   [H, Σ_a A^a† A^a] = 0
# (summed over jumps). For paired single-site (σ⁺_j, σ⁻_j) the per-pair sum
# is exactly identity, so cancellation is automatic on any Hamiltonian.
# For paired 2-site jumps like A = σ⁺_1 σ⁻_2 paired with A† = σ⁻_1 σ⁺_2:
#   A†A + AA† = (I - Z₁Z₂)/2
# which does NOT trivially commute with H_XXX (e.g., [Z₁Z₂, X₂X₃] ≠ 0).
# So if the proof's contact-term claim is correct, B_bohr and B_time should
# differ measurably for this fixture, and the rel_norm of the time-domain
# Lindbladian should exceed the Bohr-domain rel_norm.
#
# This script measures that gap.

using QuantumFurnace
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "scratch_non_hermitian_kms_db.jl"))

const SIGMA_PLUS_2x2  = ComplexF64[0 1; 0 0]
const SIGMA_MINUS_2x2 = ComplexF64[0 0; 1 0]

"""
    build_2site_nh_pair(ham; site1=1, site2=2, basis=ham.eigvecs) -> Vector{JumpOp}

Paired (σ⁺_{site1} σ⁻_{site2}, σ⁻_{site1} σ⁺_{site2}). Both jumps marked
`hermitian=false`. With `site1=1, site2=2` the per-pair sum is
(I - Z_1 Z_2)/2 which does NOT commute with the n=3 H_XXX boundary terms
involving site 2 ↔ site 3. This is the multi-site stress fixture.
"""
function build_2site_nh_pair(ham::HamHam; site1::Int=1, site2::Int=2,
                              basis::AbstractMatrix=ham.eigvecs)
    n = Int(log2(size(ham.data, 1)))
    @assert 1 <= site1 <= n && 1 <= site2 <= n && site1 != site2
    norm_fac = 1.0 / sqrt(2)
    ops = JumpOp[]
    for (op_a, op_b) in (
            (SIGMA_PLUS_2x2,  SIGMA_MINUS_2x2),
            (SIGMA_MINUS_2x2, SIGMA_PLUS_2x2),
        )
        op = Matrix(pad_term([op_a], n, site1)) *
             Matrix(pad_term([op_b], n, site2)) .* norm_fac
        op_eb = basis' * op * basis
        push!(ops, JumpOp(op, op_eb, false, false))
    end
    return Vector{JumpOp}(ops)
end

function run_2site_diagnostic(beta::Real)
    ham = _load_n3_ham(beta)
    println("=" ^ 78)
    @printf("β = %g — n=3 paired 2-site (σ⁺_1 σ⁻_2, σ⁻_1 σ⁺_2)\n", beta)
    println("=" ^ 78)

    domains = [BohrDomain(), EnergyDomain(), TimeDomain(), TrotterDomain()]
    for dom in domains
        try
            cfg = nh_config(dom; beta=beta, num_qubits=3)
            if dom isa TrotterDomain
                trotter = make_trotter_for_config(ham, cfg)
                jumps = build_2site_nh_pair(ham; basis=trotter.eigvecs)
                L = Matrix{ComplexF64}(construct_lindbladian(jumps, cfg, ham; trotter=trotter))
                gibbs_T = trotter_gibbs_target(ham, trotter)
                res = verify_detailed_balance(L, gibbs_T)
            else
                jumps = build_2site_nh_pair(ham)
                L = Matrix{ComplexF64}(construct_lindbladian(jumps, cfg, ham))
                res = verify_detailed_balance(L, ham.gibbs)
            end
            @printf("  %-12s  rel_norm = %.3e   fp_resid = %.3e   D_norm = %.3e\n",
                    string(typeof(dom).name.name),
                    res.relative_norm, res.fixed_point_residual, res.discriminant_norm)
        catch e
            @printf("  %-12s  ERROR: %s\n", string(typeof(dom).name.name),
                    first(sprint(showerror, e), 200))
        end
    end
    println()
end

if abspath(PROGRAM_FILE) == @__FILE__
    for β in (5.0, 10.0)
        run_2site_diagnostic(β)
    end
end
