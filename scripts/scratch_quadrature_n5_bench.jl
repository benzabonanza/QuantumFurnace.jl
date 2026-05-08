#!/usr/bin/env julia
#
# Cost probe: dense quadrature-error pipeline at n=4 → predicts n=5 feasibility.
#
# Strategy
# --------
# 1) At n=4 (cheap), measure wall-time of the SAME pipeline used in
#    `scratch_{energy,time}_ref_convergence.jl`:
#        build L_b   (BohrDomain)
#        build L_t   (TimeDomain   or EnergyDomain)
#        opnorm(L_t − L_b)        (post-qf-etx.2: gnf is grid-independent)
#    for the kinky filter at a few r values (r ∈ {6, 9, 12}).
# 2) Then run a SINGLE point at n=5 (β=10) for r=10 (kinky) — the cheapest
#    interesting r — to fit the scaling and decide dense vs Krylov, and Bohr-
#    vs-Energy reference, before committing to a sweep.
#
# We only print wall times, not convergence values — this is purely a cost
# probe.  All quadrature numbers come from the existing convergence scripts.
#
# Why kinky only?  Per the canonical recipe (qf-7xt) only kinky needs r_D ≥ 10
# at ε=1e-6 — Gaussian and smooth saturate at r_D ≤ 5 and so were never
# bottlenecks.  The user is asking specifically about the kinky case.

using Printf
using LinearAlgebra
using Random
using QuantumFurnace
using QuantumFurnace: construct_lindbladian, _load_hamiltonian_bson

Random.seed!(20260506)

const TAIL_C = 8.0
const ETA    = 1e-3

# ── Setup builders ──────────────────────────────────────────────────────────
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

function build_kinky_cfg(domain, n::Int, β::Real, r::Int, w0::Real)
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

# ── n=4 calibration sweep across (ref_domain, r) for kinky ──────────────────
function bench_one(n::Int, β::Real, r::Int; ref::Symbol)
    σ  = 1.0 / β
    ham_path = "/Users/bence/code/QuantumFurnace.jl/hamiltonians/heis_xxx_zzdisordered_periodic_n$(n).bson"
    ham = _load_hamiltonian_bson(ham_path, β)
    jumps = build_jumps(ham, n)

    H_norm = opnorm(ham.data)
    omega_range = 2 * (H_norm + TAIL_C * σ)
    w0 = omega_range / 2^r

    if ref === :bohr
        cfg_ref = build_kinky_cfg(BohrDomain(), n, β, r, w0)
    elseif ref === :energy
        cfg_ref = build_kinky_cfg(EnergyDomain(), n, β, r, w0)
    end
    cfg_t = build_kinky_cfg(TimeDomain(), n, β, r, w0)

    GC.gc()
    t_ref = @elapsed L_ref = construct_lindbladian(jumps, cfg_ref, ham; include_coherent = false)
    GC.gc()
    t_t   = @elapsed L_t   = construct_lindbladian(jumps, cfg_t, ham; include_coherent = false)
    GC.gc()
    # Post-qf-etx.2: gnf is grid-independent so plain L_t - L_ref is the error.
    t_diff = @elapsed err = opnorm(L_t .- L_ref)
    return (; t_ref, t_t, t_diff, t_total = t_ref + t_t + t_diff, err)
end

println("\n=== n=4 calibration (kinky, β=10), all 3 stages timed ===")
@printf "%-8s %4s %10s %10s %10s %10s %14s\n" "ref" "r" "t_ref(s)" "t_time(s)" "t_diff(s)" "total(s)" "‖ΔL‖"
for ref in (:bohr, :energy)
    for r in (6, 9, 12)
        res = bench_one(4, 10.0, r; ref)
        @printf "%-8s %4d %10.3f %10.3f %10.3f %10.3f %14.6e\n" string(ref) r res.t_ref res.t_t res.t_diff res.t_total res.err
    end
    println()
end

# ── n=5 single point at r=10 (cheapest interesting) ────────────────────────
println("\n=== n=5 single point (kinky, β=10, r=10), all 3 stages timed ===")
@printf "%-8s %4s %10s %10s %10s %10s %14s\n" "ref" "r" "t_ref(s)" "t_time(s)" "t_diff(s)" "total(s)" "‖ΔL‖"
for ref in (:bohr, :energy)
    res = bench_one(5, 10.0, 10; ref)
    @printf "%-8s %4d %10.3f %10.3f %10.3f %10.3f %14.6e\n" string(ref) 10 res.t_ref res.t_t res.t_diff res.t_total res.err
end
