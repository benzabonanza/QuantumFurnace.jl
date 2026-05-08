#!/usr/bin/env julia
# Sanity check: does the GQSP path achieve comparable thermalization to the
# matrix-exp baseline? Both should reach a similar trace distance to the Gibbs
# state, with the gap between them bounded by the polynomial truncation error
# (O((δα)^{d+1}) per step × n_steps).

using Pkg
using QuantumFurnace
using LinearAlgebra
using Printf
using Random
using BSON

# Mirror the test_helpers locked parameters
const NUM_QUBITS = 3
const BETA = 10.0
const SIGMA = 1.0 / BETA
const NUM_ENERGY_BITS = 12
const W0 = 0.05
const T0 = 2pi / (2^NUM_ENERGY_BITS * W0)
const NUM_TROTTER_STEPS_PER_T0 = 10

function build_cfg(domain; with_gqsp::Bool, gqsp_degree::Int=1, delta::Float64, mixing_time::Float64)
    Config(;
        sim=Thermalize(), domain=domain, construction=KMS(),
        num_qubits=NUM_QUBITS, with_linear_combination=true,
        beta=BETA, sigma=SIGMA, a=BETA / 30.0, s=0.4,
        num_energy_bits=NUM_ENERGY_BITS, w0=W0, t0=T0,
        num_trotter_steps_per_t0=NUM_TROTTER_STEPS_PER_T0,
        mixing_time=mixing_time, delta=delta,
        with_gqsp=with_gqsp, gqsp_degree=gqsp_degree,
    )
end

function load_n3_system()
    project_root = Pkg.project().path |> dirname
    ham_path = joinpath(project_root, "hamiltonians", "heis_disordered_periodic_n$(NUM_QUBITS).bson")
    raw = open(ham_path) do io
        BSON.parse(io)
    end
    fields = raw[:hamiltonian][:data]
    cache = IdDict()
    init = QuantumFurnace
    raw_nt = (
        matrix = BSON.raise_recursive(fields[1], cache, init)::Matrix{ComplexF64},
        terms = Vector{Vector{Matrix{ComplexF64}}}(BSON.raise_recursive(fields[4], cache, init)),
        base_coeffs = BSON.raise_recursive(fields[5], cache, init)::Vector{Float64},
        disordering_term = let dt = BSON.raise_recursive(fields[6], cache, init)
            dt === nothing ? nothing : Vector{Matrix{ComplexF64}}(dt)
        end,
        disordering_coeffs = let dc = BSON.raise_recursive(fields[7], cache, init)
            dc === nothing ? nothing : Vector{Float64}(dc)
        end,
        eigvals = BSON.raise_recursive(fields[8], cache, init)::Vector{Float64},
        eigvecs = BSON.raise_recursive(fields[9], cache, init)::Matrix{ComplexF64},
        nu_min = Float64(fields[10]),
        shift = Float64(fields[11]),
        rescaling_factor = Float64(fields[12]),
        periodic = Bool(fields[13]),
    )
    return HamHam(raw_nt, BETA)
end

ham = load_n3_system()

# Build jumps as in test_helpers
function build_jumps(ham, num_qubits)
    jump_paulis = [[X], [Y], [Z]]
    n = length(jump_paulis) * num_qubits
    jump_norm = sqrt(n)
    jumps = JumpOp[]
    for pauli in jump_paulis
        for site in 1:num_qubits
            jop = Matrix(pad_term(pauli, num_qubits, site)) ./ jump_norm
            jeig = ham.eigvecs' * jop * ham.eigvecs
            herm = jop == jop'
            push!(jumps, JumpOp(jop, jeig, jop == transpose(jop), herm))
        end
    end
    return jumps
end

jumps = build_jumps(ham, NUM_QUBITS)

println("== n=3 Heisenberg, KMS, TimeDomain. Compare GQSP d=1 vs matrix-exp ==\n")

@printf("%6s  %8s  %5s  %15s  %15s  %12s\n",
        "δ", "T_mix", "steps", "td(exp,Gibbs)", "td(gqsp,Gibbs)", "‖Δ ρ‖")

for (δ, T) in [(1e-2, 0.5), (1e-2, 2.0), (5e-3, 2.0), (5e-3, 5.0), (1e-2, 5.0)]
    cfg_e = build_cfg(TimeDomain(); with_gqsp=false, delta=δ, mixing_time=T)
    cfg_g = build_cfg(TimeDomain(); with_gqsp=true, gqsp_degree=1, delta=δ, mixing_time=T)
    seed = 1234
    r_e = run_thermalize(jumps, cfg_e, ham; rng=Xoshiro(seed))
    r_g = run_thermalize(jumps, cfg_g, ham; rng=Xoshiro(seed))
    n_steps = Int(ceil(T / δ))
    td_e = r_e.trace_distances[end]
    td_g = r_g.trace_distances[end]
    diff = opnorm(r_e.final_dm - r_g.final_dm)
    @printf("%6.0e  %8.2f  %5d  %15.6e  %15.6e  %12.3e\n",
            δ, T, n_steps, td_e, td_g, diff)
end
