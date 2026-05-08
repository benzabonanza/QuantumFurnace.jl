#!/usr/bin/env julia
#
# n=5 spot check: Krylov-SVD vs dense opnorm crossover.
#
# Targeted subset: β=10 only, one representative r_D per filter.
# Goal: determine if n=5 (d=32, d²=1024, superop=16MB) is the crossover.
#
# Launch: JULIA_NUM_THREADS=8 julia --project scripts/scratch_krylov_vs_dense_n5_spot.jl

using Printf
using LinearAlgebra
using Random
using QuantumFurnace
using QuantumFurnace: construct_lindbladian, _load_hamiltonian_bson,
    hs_operator_norm_krylov

Random.seed!(20260506)

const TAIL_C = 8.0
const ETA    = 1e-3
const N      = 5
const BETA   = 10.0
const SIGMA  = 1.0 / BETA

# Spot checks: one interesting r per filter + kinky at two r values
const CASES = [
    (:gaussian, 5),
    (:smooth,   5),
    (:kinky,    8),
    (:kinky,   10),
    (:kinky,   12),
]

ham_path = joinpath(@__DIR__, "..", "hamiltonians",
    "heis_xxx_zzdisordered_periodic_n$(N).bson")
ham = _load_hamiltonian_bson(ham_path, BETA)

jump_norm = sqrt(3 * N)
jumps = JumpOp[]
for pauli in ([X], [Y], [Z])
    for site in 1:N
        op = Matrix(pad_term(pauli, N, site)) ./ jump_norm
        op_eb = ham.eigvecs' * op * ham.eigvecs
        push!(jumps, JumpOp(op, op_eb, op == transpose(op), op == op'))
    end
end

H_norm = opnorm(ham.data)
omega_range = 2 * (H_norm + TAIL_C * SIGMA)
d = size(ham.data, 1)

function build_cfg(domain, filter::Symbol, r::Int, w0::Real)
    t0 = 2π / (2^r * w0)
    common = (
        sim = Lindbladian(), domain = domain, construction = KMS(),
        num_qubits = N, beta = BETA, sigma = SIGMA,
        num_energy_bits = r, w0 = w0, t0 = t0,
    )
    if filter === :gaussian
        sigma_gamma = SIGMA
        w_gamma = BETA * (SIGMA^2 + sigma_gamma^2) / 2
        return Config(; common...,
            with_linear_combination = false,
            gaussian_parameters = (w_gamma, sigma_gamma),
        )
    elseif filter === :smooth
        return Config(; common...,
            with_linear_combination = true, a = 0.0, s = 0.25, eta = ETA,
        )
    elseif filter === :kinky
        return Config(; common...,
            with_linear_combination = true, a = 0.0, s = 0.0, eta = ETA,
        )
    end
end

println("n=$N  d=$d  d²=$(d^2)  β=$BETA  ‖H‖=$(round(H_norm, digits=4))  ω-range=$(round(omega_range, digits=4))")
println("Julia threads: $(Threads.nthreads()), BLAS threads: $(BLAS.get_num_threads())")
println()
@printf "%-9s %4s | %14s %10s | %14s %10s | %7s  %s\n" "filter" "r_D" "‖ΔL‖_dense" "t_dense" "‖ΔL‖_krylov" "t_krylov" "ratio" "agreement"
println("-" ^ 95)

for (filt, r) in CASES
    w0 = omega_range / 2^r
    cfg_ref  = build_cfg(BohrDomain(),   filt, r, w0)
    cfg_test = build_cfg(EnergyDomain(), filt, r, w0)

    # --- Dense path ---
    # Post-qf-etx.2: gnf is grid-independent so plain L_test - L_ref is the error.
    GC.gc()
    t_dense = @elapsed begin
        L_ref  = construct_lindbladian(jumps, cfg_ref, ham)
        L_test = construct_lindbladian(jumps, cfg_test, ham)
        err_dense = opnorm(L_test .- L_ref)
    end

    # --- Krylov path ---
    GC.gc()
    t_krylov = @elapsed begin
        ws_ref  = Workspace(cfg_ref, ham, jumps)
        ws_test = Workspace(cfg_test, ham, jumps)
        buf = zeros(ComplexF64, d, d)
        diff!(out, X) = begin
            apply_lindbladian!(ws_ref,  X, cfg_ref,  ham)
            copyto!(buf, ws_ref.scratch.rho_out)
            apply_lindbladian!(ws_test, X, cfg_test, ham)
            out .= ws_test.scratch.rho_out .- buf; nothing
        end
        adj!(out, X) = begin
            apply_adjoint_lindbladian!(ws_ref,  X, cfg_ref,  ham)
            copyto!(buf, ws_ref.scratch.rho_out)
            apply_adjoint_lindbladian!(ws_test, X, cfg_test, ham)
            out .= ws_test.scratch.rho_out .- buf; nothing
        end
        err_krylov = hs_operator_norm_krylov(diff!, adj!, d)
    end

    ratio = t_dense / t_krylov
    rel_diff = err_dense > 1e-15 ? abs(err_dense - err_krylov) / err_dense : abs(err_dense - err_krylov)
    ok = rel_diff < 1e-4 ? "✓" : "✗ $(round(rel_diff, sigdigits=2))"
    @printf "%-9s %4d | %14.6e %9.3fs | %14.6e %9.3fs | %6.2fx  %s\n" string(filt) r err_dense t_dense err_krylov t_krylov ratio ok
end
