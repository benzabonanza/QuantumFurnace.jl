#!/usr/bin/env julia
#
# Dissipative-only Krylov quadrature error sweep for smooth Metropolis at n=5,6,7.
#
# Measures ‖L_diss(Time, r_D) - L_diss(Energy, r_ref)‖_op via Krylov-SVD,
# with include_coherent=false so only the dissipator part is compared.
#
# Reference: EnergyDomain at r_ref=8 (smooth Metro saturates by r=6).
# Test: TimeDomain at varying r_D — includes NUFFT approximation error.
#
# Launch: JULIA_NUM_THREADS=8 julia --project scripts/scratch_dissipative_krylov_quadrature.jl

using Printf
using LinearAlgebra
using Random
using QuantumFurnace
using QuantumFurnace: construct_lindbladian, _load_hamiltonian_bson,
    hs_operator_norm_krylov

Random.seed!(20260506)

const TAIL_C = 8.0
const ETA    = 1e-3
const BETA   = 10.0
const SIGMA  = 1.0 / BETA

const R_REF    = 8
const R_GRID   = [3, 4, 5, 6]
const NS       = [5, 6, 7]

function load_fixture(n::Int, β::Real)
    ham_path = joinpath(@__DIR__, "..", "hamiltonians",
        "heis_xxx_zzdisordered_periodic_n$(n).bson")
    ham = _load_hamiltonian_bson(ham_path, β)
    jump_norm = sqrt(3 * n)
    jumps = JumpOp[]
    for pauli in ([X], [Y], [Z])
        for site in 1:n
            op = Matrix(pad_term(pauli, n, site)) ./ jump_norm
            op_eb = ham.eigvecs' * op * ham.eigvecs
            push!(jumps, JumpOp(op, op_eb, op == transpose(op), op == op'))
        end
    end
    H_norm = opnorm(ham.data)
    return ham, jumps, H_norm
end

function build_cfg(n::Int, r::Int, w0::Real; domain=EnergyDomain())
    t0 = 2π / (2^r * w0)
    Config(;
        sim = Lindbladian(), domain = domain, construction = KMS(),
        num_qubits = n, beta = BETA, sigma = SIGMA,
        num_energy_bits = r, w0 = w0, t0 = t0,
        with_linear_combination = true, a = 0.0, s = 0.25, eta = ETA,
    )
end

function krylov_diss_error(cfg_ref, cfg_test, ham, jumps)
    d = size(ham.data, 1)
    ws_ref  = Workspace(cfg_ref, ham, jumps)
    ws_test = Workspace(cfg_test, ham, jumps)

    # Post-qf-etx.2: gamma_norm_factor is grid-independent (= 1.0), so the
    # difference closure is plain L_ref - L_test (no /gnf rescaling needed).
    diff!(out, X) = begin
        apply_lindbladian!(ws_ref, X, cfg_ref, ham; include_coherent=false)
        copyto!(out, ws_ref.scratch.rho_out)
        apply_lindbladian!(ws_test, X, cfg_test, ham; include_coherent=false)
        out .-= ws_test.scratch.rho_out
        nothing
    end
    adj!(out, X) = begin
        apply_adjoint_lindbladian!(ws_ref, X, cfg_ref, ham; include_coherent=false)
        copyto!(out, ws_ref.scratch.rho_out)
        apply_adjoint_lindbladian!(ws_test, X, cfg_test, ham; include_coherent=false)
        out .-= ws_test.scratch.rho_out
        nothing
    end
    return hs_operator_norm_krylov(diff!, adj!, d)
end

println("Dissipative-only Krylov quadrature error: TimeDomain(r_D) vs EnergyDomain(r_ref)")
println("Filter: smooth Metro (s=0.25), β=$BETA  σ=$SIGMA  r_ref=$R_REF  TAIL_C=$TAIL_C")
println("Julia threads: $(Threads.nthreads()), BLAS threads: $(BLAS.get_num_threads())")
println()

for n in NS
    ham, jumps, H_norm = load_fixture(n, BETA)
    d = size(ham.data, 1)
    omega_range = 2 * (H_norm + TAIL_C * SIGMA)
    w0_ref = omega_range / 2^R_REF
    cfg_ref = build_cfg(n, R_REF, w0_ref)

    println("=" ^ 70)
    @printf "n=%d  d=%d  ‖H‖=%.4f  ω-range=%.4f  ref=Energy(r=%d)  test=Time(r_D)\n" n d H_norm omega_range R_REF
    println("=" ^ 70)
    @printf "%4s | %14s %10s\n" "r_D" "‖ΔL_diss‖" "wall"
    println("-" ^ 35)
    flush(stdout)

    for r in R_GRID
        w0 = omega_range / 2^r
        cfg_test = build_cfg(n, r, w0; domain=TimeDomain())

        GC.gc()
        t = @elapsed err = krylov_diss_error(cfg_ref, cfg_test, ham, jumps)
        @printf "%4d | %14.6e %9.3fs\n" r err t
        flush(stdout)
    end
    println()
    flush(stdout)
end
