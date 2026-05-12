#!/usr/bin/env julia
# Audit the GQSP coherent-step trace drift end-to-end:
#   1. Compute α_a (block-encoding norm) and ‖B_a‖_op for each jump in the n=3 fixture
#   2. Compute ε_a = ‖U_a^† U_a − I‖_op for U_a = f_1(B_a/α_a) (the per-step non-unitarity)
#   3. Compare to the analytical bound ε_a ≲ (δ_eff α_a)² (slope-2 from Bessel tail)
#   4. Predict cumulative trace drift after n_steps and compare to observed run_thermalize drift
#
# This nails down whether the trace defect we observe (~7e-7 over 50 steps) is the expected
# polynomial-truncation effect or a hidden bug.

using Pkg
using QuantumFurnace
using LinearAlgebra
using Printf
using Random
using BSON

const NUM_QUBITS = 3
const BETA = 10.0
const SIGMA = 1.0 / BETA
const NUM_ENERGY_BITS = 12
const W0 = 0.05
const T0 = 2pi / (2^NUM_ENERGY_BITS * W0)
const NUM_TROTTER_STEPS_PER_T0 = 10

function load_n3_system()
    project_root = Pkg.project().path |> dirname
    ham_path = joinpath(project_root, "hamiltonians", "heis_xxx_zzdisordered_periodic_n$(NUM_QUBITS).bson")
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

# ----------------------------------------------------------------------
# Build config exactly matching test_gqsp_thermalize.jl smoke test
# (delta=1e-3, mixing_time=0.05; in run_thermalize the per-jump effective
#  step is delta_eff = delta * n_jumps = 9e-3 because rescale_by_inv_prob=true)
# ----------------------------------------------------------------------
ham = load_n3_system()
jumps = build_jumps(ham, NUM_QUBITS)
n_jumps = length(jumps)

cfg = Config(;
    sim=Thermalize(), domain=TimeDomain(), construction=KMS(),
    num_qubits=NUM_QUBITS, with_linear_combination=true,
    beta=BETA, sigma=SIGMA, a=BETA / 30.0, s=0.4,
    num_energy_bits=NUM_ENERGY_BITS, w0=W0, t0=T0,
    num_trotter_steps_per_t0=NUM_TROTTER_STEPS_PER_T0,
    mixing_time=0.05, delta=1e-3,
    with_gqsp=true, gqsp_degree=1,
)

precomputed = QuantumFurnace._precompute_data(cfg, ham)
γ_nf = precomputed.gamma_norm_factor
b_minus = precomputed.b_minus
b_plus = precomputed.b_plus
@printf("γ_nf = %.6f\n", γ_nf)
@printf("‖b_-‖_ℓ¹ = %.6f, ‖b_+‖_ℓ¹ = %.6f\n", sum(abs, values(b_minus)), sum(abs, values(b_plus)))
@printf("t0 = %.6f, t0² = %.6e\n\n", cfg.t0, cfg.t0^2)

# ----------------------------------------------------------------------
# For each jump: build B_a, α_a, U_a; measure non-unitarity ε_a = ‖U_a†U_a − I‖
# ----------------------------------------------------------------------
delta_eff = cfg.delta * n_jumps   # what _precompute_coherent_unitary actually uses

@printf("delta = %.4f, delta_eff = delta × n_jumps = %.4f\n\n", cfg.delta, delta_eff)
@printf("%4s  %8s  %8s  %8s  %14s  %14s  %14s\n",
        "a", "‖A‖²", "‖B_a‖", "α_a", "‖B_a‖/α_a", "δ_eff·α_a", "ε=‖U†U−I‖")

global ε_total = 0.0
global αs = Float64[]
for a in 1:n_jumps
    global ε_total
    global αs
    jump = jumps[a]
    B = QuantumFurnace.B_time([jump], ham, b_minus, b_plus, cfg.t0, cfg.beta, cfg.sigma)
    rmul!(B, γ_nf)
    α_a = QuantumFurnace._gqsp_block_encoding_alpha(jump, b_minus, b_plus, cfg.t0, γ_nf)
    push!(αs, α_a)
    Bnorm = opnorm(Hermitian((B + B') / 2))
    A_norm_sq = opnorm(jump.data)^2
    U_a = QuantumFurnace._gqsp_apply_polynomial(B, α_a, delta_eff, 1)
    ε_a = opnorm(U_a' * U_a - I(size(U_a, 1)))
    ε_total += ε_a
    @printf("%4d  %8.4f  %8.4f  %8.4f  %14.6e  %14.6e  %14.6e\n",
            a, A_norm_sq, Bnorm, α_a, Bnorm / α_a, delta_eff * α_a, ε_a)
end
ε_avg = ε_total / n_jumps
@printf("\n  avg ε per jump = %.6e\n", ε_avg)
α_max = maximum(αs)
@printf("  max α_a = %.4f, max δ_eff α_a = %.4f\n", α_max, delta_eff * α_max)

# ----------------------------------------------------------------------
# Predicted vs observed cumulative trace drift
# ----------------------------------------------------------------------
n_steps = Int(ceil(cfg.mixing_time / cfg.delta))
@printf("\nn_steps = %d\n", n_steps)
@printf("Predicted total trace drift (linear accumulation): n_steps × ε_avg ≲ %.3e\n",
        n_steps * ε_avg)
@printf("Per-step Bessel-tail prediction (δ_eff·α_avg)²/2: %.3e\n",
        (delta_eff * (sum(αs) / n_jumps))^2 / 2)

# ----------------------------------------------------------------------
# Now actually run run_thermalize and compare
# ----------------------------------------------------------------------
seed = 1234
result = run_thermalize(jumps, cfg, ham; rng=Xoshiro(seed))
trace_defect = abs(tr(result.final_dm) - 1)
@printf("\nObserved trace defect after run_thermalize: %.6e\n", trace_defect)

# Same test with d=2 should reduce drift by another factor of (δα)
cfg_d2 = Config(;
    sim=Thermalize(), domain=TimeDomain(), construction=KMS(),
    num_qubits=NUM_QUBITS, with_linear_combination=true,
    beta=BETA, sigma=SIGMA, a=BETA / 30.0, s=0.4,
    num_energy_bits=NUM_ENERGY_BITS, w0=W0, t0=T0,
    num_trotter_steps_per_t0=NUM_TROTTER_STEPS_PER_T0,
    mixing_time=0.05, delta=1e-3,
    with_gqsp=true, gqsp_degree=2,
)
result_d2 = run_thermalize(jumps, cfg_d2, ham; rng=Xoshiro(seed))
@printf("Observed trace defect with d=2:                %.6e\n", abs(tr(result_d2.final_dm) - 1))

# Matrix-exp baseline (no drift expected)
cfg_exp = Config(;
    sim=Thermalize(), domain=TimeDomain(), construction=KMS(),
    num_qubits=NUM_QUBITS, with_linear_combination=true,
    beta=BETA, sigma=SIGMA, a=BETA / 30.0, s=0.4,
    num_energy_bits=NUM_ENERGY_BITS, w0=W0, t0=T0,
    num_trotter_steps_per_t0=NUM_TROTTER_STEPS_PER_T0,
    mixing_time=0.05, delta=1e-3,
    with_gqsp=false, gqsp_degree=1,
)
result_exp = run_thermalize(jumps, cfg_exp, ham; rng=Xoshiro(seed))
@printf("Observed trace defect with matrix-exp:         %.6e\n", abs(tr(result_exp.final_dm) - 1))
