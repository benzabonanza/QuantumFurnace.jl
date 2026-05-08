#!/usr/bin/env julia
#
# Single-point timing probe at n=6, r=8 (cheapest meaningful) — extrapolates
# the dense build cost at higher r.  Runs ONLY the kinky filter, BohrDomain
# reference, since EnergyDomain ref doesn't probe quadrature error and only
# adds wall-time.

using Printf
using LinearAlgebra
using Random
using QuantumFurnace
using QuantumFurnace: construct_lindbladian, _load_hamiltonian_bson

Random.seed!(20260506)

const TAIL_C = 8.0
const ETA    = 1e-3

function build_jumps(ham::HamHam, n::Int)
    jumps = JumpOp[]
    jump_norm = sqrt(3 * n)
    for pauli in ([X], [Y], [Z])
        for site in 1:n
            op = Matrix(pad_term(pauli, n, site)) ./ jump_norm
            op_eb = ham.eigvecs' * op * ham.eigvecs
            push!(jumps, JumpOp(op, op_eb, op == transpose(op), op == op'))
        end
    end
    return jumps
end

function kinky_cfg(domain, n::Int, β::Real, r::Int, w0::Real)
    σ  = 1.0 / β
    t0 = 2π / (2^r * w0)
    return Config(;
        sim = Lindbladian(), domain = domain, construction = KMS(),
        num_qubits = n, beta = β, sigma = σ,
        num_energy_bits = r, w0 = w0, t0 = t0,
        with_linear_combination = true,
        a = 0.0, s = 0.0, eta = ETA,
    )
end

readmem() = begin
    txt = read("/proc/meminfo", String)
    avail = match(r"MemAvailable:\s+(\d+)", txt).captures[1] |> x -> parse(Int, x)
    return avail / 1024
end

n = 6; β = 10.0; r = 8
σ = 1/β
ham = _load_hamiltonian_bson(
    "/Users/bence/code/QuantumFurnace.jl/hamiltonians/heis_xxx_zzdisordered_periodic_n$(n).bson", β)
jumps = build_jumps(ham, n)
H_norm = opnorm(ham.data)
omega_range = 2 * (H_norm + TAIL_C * σ)
w0 = omega_range / 2^r
cfg_b = kinky_cfg(BohrDomain(), n, β, r, w0)
cfg_t = kinky_cfg(TimeDomain(), n, β, r, w0)

println("\n=== n=$n r=$r kinky β=$β single-point timing ===")
@printf "%-30s %12s\n" "stage" "wall (s)"
@printf "(start) MemAvail = %.1f MiB\n" readmem()

GC.gc()
t_b = @elapsed L_b = construct_lindbladian(jumps, cfg_b, ham; include_coherent=false)
@printf "%-30s %12.3f\n" "build L_b (Bohr)" t_b
@printf "MemAvail = %.1f MiB\n" readmem()

GC.gc()
t_t = @elapsed L_t = construct_lindbladian(jumps, cfg_t, ham; include_coherent=false)
@printf "%-30s %12.3f\n" "build L_t (Time)" t_t
@printf "MemAvail = %.1f MiB\n" readmem()

GC.gc()
# Post-qf-etx.2: gnf is grid-independent so plain L_t - L_b is the error.
t_diff = @elapsed err = opnorm(L_t .- L_b)
@printf "%-30s %12.3f\n" "opnorm" t_diff
@printf "‖ΔL‖ = %.3e\n" err
@printf "%-30s %12.3f\n" "TOTAL" (t_b + t_t + t_diff)
