#!/usr/bin/env julia
#
# Probe: does bumping BLAS threads from 4 → 8 (or higher) speed up the dense
# quadrature pipeline?  At n=5, r=10, kinky, β=10 the build is dominated by
# `kron!` on d×d=32×32 GEMMs — small enough that BLAS parallel scaling may
# stall.  Measure to find out.

using Printf, LinearAlgebra, Random
using QuantumFurnace
using QuantumFurnace: construct_lindbladian, _load_hamiltonian_bson

Random.seed!(20260506)

const TAIL_C = 8.0
const ETA    = 1e-3

function build_jumps(ham::HamHam, n::Int)
    jumps = JumpOp[]
    jump_norm = sqrt(3 * n)
    for pauli in ([X], [Y], [Z]), site in 1:n
        op = Matrix(pad_term(pauli, n, site)) ./ jump_norm
        op_eb = ham.eigvecs' * op * ham.eigvecs
        push!(jumps, JumpOp(op, op_eb, op == transpose(op), op == op'))
    end
    jumps
end

function kinky_cfg(domain, n, β, r, w0)
    σ  = 1/β; t0 = 2π / (2^r * w0)
    Config(; sim=Lindbladian(), domain=domain, construction=KMS(),
        num_qubits=n, beta=β, sigma=σ, num_energy_bits=r, w0=w0, t0=t0,
        with_linear_combination=true, a=0.0, s=0.0, eta=ETA)
end

n, β, r = 5, 10.0, 10
σ = 1/β
ham = _load_hamiltonian_bson(
    "/Users/bence/code/QuantumFurnace.jl/hamiltonians/heis_xxx_zzdisordered_periodic_n$(n).bson", β)
jumps = build_jumps(ham, n)
H_norm = opnorm(ham.data)
omega_range = 2 * (H_norm + TAIL_C * σ)
w0 = omega_range / 2^r
cfg_b = kinky_cfg(BohrDomain(), n, β, r, w0)
cfg_t = kinky_cfg(TimeDomain(), n, β, r, w0)

# Warmup once
_ = construct_lindbladian(jumps, cfg_b, ham; include_coherent=false)
_ = construct_lindbladian(jumps, cfg_t, ham; include_coherent=false)

println("\n=== n=5 r=10 kinky β=10 BLAS thread sweep ===")
@printf "%-7s %12s %12s %12s\n" "BLAS" "build_b (s)" "build_t (s)" "opnorm (s)"
for nthr in (1, 2, 4, 8)
    BLAS.set_num_threads(nthr)
    GC.gc()
    t_b = @elapsed L_b = construct_lindbladian(jumps, cfg_b, ham; include_coherent=false)
    GC.gc()
    t_t = @elapsed L_t = construct_lindbladian(jumps, cfg_t, ham; include_coherent=false)
    GC.gc()
    # Post-qf-etx.2: gnf is grid-independent so plain L_t - L_b suffices.
    diff = L_t .- L_b
    GC.gc()
    t_op = @elapsed _ = opnorm(diff)
    @printf "%-7d %12.3f %12.3f %12.3f\n" nthr t_b t_t t_op
end
