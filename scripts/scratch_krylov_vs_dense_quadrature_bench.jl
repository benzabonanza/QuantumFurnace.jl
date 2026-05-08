#!/usr/bin/env julia
#
# Krylov-SVD vs dense opnorm cross-check for dissipative quadrature error.
#
# Purpose: verify that the matrix-free Krylov path (apply_lindbladian! closures +
# hs_operator_norm_krylov) matches the dense path (construct_lindbladian + opnorm)
# at n=3,4, and measure wall times for both.
#
# Launch with: JULIA_NUM_THREADS=8 julia --project scripts/scratch_krylov_vs_dense_quadrature_bench.jl
#
# Sweep: n ∈ {3,4}, β ∈ {5,10,20}, filter ∈ {gaussian, smooth, kinky},
# r_D values from the canonical recipe (qf-7xt).

using Printf
using LinearAlgebra
using Random
using QuantumFurnace
using QuantumFurnace: construct_lindbladian, _load_hamiltonian_bson,
    hs_operator_norm_krylov

Random.seed!(20260506)

const TAIL_C = 8.0
const ETA    = 1e-3
const BETAS  = [5.0, 10.0, 20.0]
const NS     = [3, 4]

# Per-filter r_D sweep: representative points from the canonical recipe.
# Gaussian/smooth saturate fast; kinky needs the full range.
const R_GRID = Dict(
    :gaussian => [3, 4, 5, 6],
    :smooth   => [3, 4, 5, 6],
    :kinky    => [6, 8, 10, 12],
)

# ── Fixture builders ────────────────────────────────────────────────────────

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

function build_cfg(domain, filter::Symbol, n::Int, β::Real, r::Int, w0::Real)
    σ  = 1.0 / β
    t0 = 2π / (2^r * w0)
    common = (
        sim = Lindbladian(), domain = domain, construction = KMS(),
        num_qubits = n, beta = β, sigma = σ,
        num_energy_bits = r, w0 = w0, t0 = t0,
    )
    if filter === :gaussian
        sigma_gamma = σ
        w_gamma = β * (σ^2 + sigma_gamma^2) / 2
        return Config(; common...,
            with_linear_combination = false,
            gaussian_parameters = (w_gamma, sigma_gamma),
        )
    elseif filter === :smooth
        return Config(; common...,
            with_linear_combination = true,
            a = 0.0, s = 0.25, eta = ETA,
        )
    elseif filter === :kinky
        return Config(; common...,
            with_linear_combination = true,
            a = 0.0, s = 0.0, eta = ETA,
        )
    end
end

# ── Dense path: construct_lindbladian + opnorm ──────────────────────────────
# Post-qf-etx.2: gnf is grid-independent so plain L_test - L_ref is the error.

function dense_error(cfg_ref, cfg_test, ham, jumps)
    L_ref  = construct_lindbladian(jumps, cfg_ref, ham)
    L_test = construct_lindbladian(jumps, cfg_test, ham)
    return opnorm(L_test .- L_ref)
end

# ── Krylov path: apply_lindbladian! closures + hs_operator_norm_krylov ──────

function krylov_error(cfg_ref, cfg_test, ham, jumps)
    d = size(ham.data, 1)
    ws_ref  = Workspace(cfg_ref, ham, jumps)
    ws_test = Workspace(cfg_test, ham, jumps)
    buf = zeros(ComplexF64, d, d)

    diff!(out, X) = begin
        apply_lindbladian!(ws_ref,  X, cfg_ref,  ham)
        copyto!(buf, ws_ref.scratch.rho_out)
        apply_lindbladian!(ws_test, X, cfg_test, ham)
        out .= ws_test.scratch.rho_out .- buf
        nothing
    end
    adj!(out, X) = begin
        apply_adjoint_lindbladian!(ws_ref,  X, cfg_ref,  ham)
        copyto!(buf, ws_ref.scratch.rho_out)
        apply_adjoint_lindbladian!(ws_test, X, cfg_test, ham)
        out .= ws_test.scratch.rho_out .- buf
        nothing
    end
    return hs_operator_norm_krylov(diff!, adj!, d)
end

# ── Main sweep ──────────────────────────────────────────────────────────────

println("Julia threads: $(Threads.nthreads()), BLAS threads: $(BLAS.get_num_threads())")
println()

for n in NS
    for β in BETAS
        ham, jumps, H_norm = load_fixture(n, β)
        σ = 1.0 / β
        omega_range = 2 * (H_norm + TAIL_C * σ)
        d = size(ham.data, 1)

        println("=" ^ 90)
        @printf "n=%d  d=%d  β=%.0f  σ=%.3f  ‖H‖=%.4f  ω-range=%.4f\n" n d β σ H_norm omega_range
        println("=" ^ 90)
        @printf "%-9s %4s %12s | %14s %10s | %14s %10s | %10s\n" "filter" "r_D" "w0" "‖ΔL‖_dense" "t_dense" "‖ΔL‖_krylov" "t_krylov" "agreement"
        println("-" ^ 90)

        for filt in (:gaussian, :smooth, :kinky)
            for r in R_GRID[filt]
                w0 = omega_range / 2^r
                cfg_ref  = build_cfg(BohrDomain(),   filt, n, β, r, w0)
                cfg_test = build_cfg(EnergyDomain(), filt, n, β, r, w0)

                GC.gc()
                t_dense = @elapsed err_dense = dense_error(cfg_ref, cfg_test, ham, jumps)

                # For the Krylov path, reference = BohrDomain, test = EnergyDomain
                GC.gc()
                t_krylov = @elapsed err_krylov = krylov_error(cfg_ref, cfg_test, ham, jumps)

                rel_diff = err_dense > 1e-15 ? abs(err_dense - err_krylov) / err_dense : abs(err_dense - err_krylov)
                ok = rel_diff < 1e-4 ? "✓" : "✗ $(round(rel_diff, sigdigits=2))"
                @printf "%-9s %4d %12.4e | %14.6e %9.3fs | %14.6e %9.3fs | %s\n" string(filt) r w0 err_dense t_dense err_krylov t_krylov ok
            end
        end
        println()
    end
end
