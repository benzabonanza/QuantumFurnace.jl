#!/usr/bin/env julia
#
# Bi-exponential mixing time verification:
#   Phase 1: Run a truncated simulation, fit bi-exponential, extrapolate steps for ε = 1e-4
#   Phase 2: Run Thermalize() with exactly that many steps and verify convergence
#
# With coherent unitary, floor ≈ 0.068 × δ = 6.8e-5, safely below target 1e-4.
#
# Usage:  OPENBLAS_NUM_THREADS=4 julia -t4 --project scripts/biexp_mixing_verify.jl

using QuantumFurnace
using LinearAlgebra
using BSON
using Random

BLAS.set_num_threads(4)
println("BLAS threads: $(BLAS.get_num_threads()), Julia threads: $(Threads.nthreads())")

# ── 1. Load Hamiltonian ───────────────────────────────────────────────────────
num_qubits = 3
beta       = 10.0

ham_path = joinpath(@__DIR__, "..", "hamiltonians",
                    "heis_disordered_periodic_n$(num_qubits).bson")

raw = open(ham_path) do io; BSON.parse(io) end
fields = raw[:hamiltonian][:data]

cache = IdDict()
data_matrix      = BSON.raise_recursive(fields[1],  cache, QuantumFurnace)::Matrix{ComplexF64}
base_terms       = Vector{Vector{Matrix{ComplexF64}}}(BSON.raise_recursive(fields[4], cache, QuantumFurnace))
base_coeffs      = BSON.raise_recursive(fields[5],  cache, QuantumFurnace)::Vector{Float64}
disordering_term = let dt = BSON.raise_recursive(fields[6], cache, QuantumFurnace)
    dt === nothing ? nothing : Vector{Matrix{ComplexF64}}(dt)
end
disordering_coeffs = let dc = BSON.raise_recursive(fields[7], cache, QuantumFurnace)
    dc === nothing ? nothing : Vector{Float64}(dc)
end
eigvals_vec      = BSON.raise_recursive(fields[8],  cache, QuantumFurnace)::Vector{Float64}
eigvecs_mat      = BSON.raise_recursive(fields[9],  cache, QuantumFurnace)::Matrix{ComplexF64}

raw_nt = (
    matrix             = data_matrix,
    terms              = base_terms,
    base_coeffs        = base_coeffs,
    disordering_term   = disordering_term,
    disordering_coeffs = disordering_coeffs,
    eigvals            = eigvals_vec,
    eigvecs            = eigvecs_mat,
    nu_min             = Float64(fields[10]),
    shift              = Float64(fields[11]),
    rescaling_factor   = Float64(fields[12]),
    periodic           = Bool(fields[13]),
)

hamiltonian = HamHam(raw_nt, beta)
println("Hamiltonian: $(num_qubits)-qubit disordered Heisenberg (periodic)")
println("  dim = $(2^num_qubits), β = $beta")

# ── 2. Build jump operators ───────────────────────────────────────────────────
jump_paulis = [[X], [Y], [Z]]
n_jumps     = length(jump_paulis) * num_qubits
norm_factor = sqrt(n_jumps)

jumps = JumpOp[]
for pauli in jump_paulis
    for site in 1:num_qubits
        op       = Matrix(pad_term(pauli, num_qubits, site)) ./ norm_factor
        op_eigen = hamiltonian.eigvecs' * op * hamiltonian.eigvecs
        push!(jumps, JumpOp(op, op_eigen, op == transpose(op), op == op'))
    end
end
println("  $(length(jumps)) jump operators")

# ── 3. Parameters ────────────────────────────────────────────────────────────
delta  = 0.001   # floor ≈ 0.068×δ = 6.8e-5 < 1e-4 target
target = 1e-4
sigma  = 1.0 / beta

# ══════════════════════════════════════════════════════════════════════════════
# PHASE 1: TRUNCATED RUN + BI-EXPONENTIAL EXTRAPOLATION
# ══════════════════════════════════════════════════════════════════════════════
println("\n" * "="^70)
println("PHASE 1: TRUNCATED EXTRAPOLATION RUN  (δ = $delta)")
println("="^70)

# Run for T=45 (~14 time constants at gap≈0.32), enough to clearly see both modes
extrap_mixing_time = 45.0
extrap_save_every  = 100
extrap_num_steps   = round(Int, extrap_mixing_time / delta)

config_extrap = Config(;
    sim                       = Thermalize(),
    domain                    = EnergyDomain(),
    construction              = KMS(),
    num_qubits                = num_qubits,
    with_linear_combination   = true,
    beta                      = beta,
    sigma                     = sigma,
    a                         = beta / 30.0,
    b                         = 0.4,
    num_energy_bits           = 12,
    w0                        = 0.05,
    t0                        = 2π / (2^12 * 0.05),
    num_trotter_steps_per_t0  = 10,
    mixing_time               = extrap_mixing_time,
    delta                     = delta,
)

println("  T = $extrap_mixing_time → $extrap_num_steps steps, save_every = $extrap_save_every")
println("  Running...")
flush(stdout)

result_extrap = run_thermalize(jumps, config_extrap, hamiltonian;
    rng        = MersenneTwister(42),
    save_every = extrap_save_every,
)

wall1 = result_extrap.metadata[:wall_time_seconds]
println("  Done in $(round(wall1; digits=1))s")
println("  Final trace distance: $(round(result_extrap.trace_distances[end]; sigdigits=4))")
println("  Data points: $(length(result_extrap.trace_distances))")

# ── Bi-exponential fit ────────────────────────────────────────────────────────
est_bi = estimate_mixing_time(result_extrap;
    model          = :biexp,
    skip_initial   = 0.2,
    target_epsilon = target,
    extrapolate    = true,
)

bifit = est_bi.biexp_fit_result
println("\n--- Bi-exponential fit ---")
println("  d(t) = $(round(bifit.amplitude_fast; sigdigits=4)) × exp(-$(round(bifit.gap_fast; sigdigits=4)) t)")
println("       + $(round(bifit.amplitude; sigdigits=4)) × exp(-$(round(bifit.gap; sigdigits=4)) t)")
println("       + $(round(bifit.offset; sigdigits=4))")
println("  R²         = $(round(est_bi.r_squared; sigdigits=6))")
println("  Slow gap   = $(round(bifit.gap; sigdigits=4))  (spectral gap)")
println("  Fast gap   = $(round(bifit.gap_fast; sigdigits=4))")
println("  Offset C   = $(round(bifit.offset; sigdigits=4))  (vs target ε = $target)")
println("  Converged  = $(est_bi.converged)")

if est_bi.mixing_time_extrapolated === nothing || isnan(est_bi.mixing_time_extrapolated)
    println("\nERROR: Bi-exponential extrapolation failed!")
    println("  Offset C = $(round(bifit.offset; sigdigits=4)) ≥ target ε = $target ?")
    exit(1)
end

predicted_time  = est_bi.mixing_time_extrapolated
predicted_steps = ceil(Int, predicted_time / delta)

println("\n  Extrapolated crossing: t = $(round(predicted_time; sigdigits=5))")
println("  Predicted delta steps: $predicted_steps")

# Also show what single-exp would have predicted, for comparison
est_single = estimate_mixing_time(result_extrap;
    model          = :single,
    skip_initial   = 0.2,
    target_epsilon = target,
    extrapolate    = true,
)
if est_single.mixing_time_extrapolated !== nothing && !isnan(est_single.mixing_time_extrapolated)
    single_steps = ceil(Int, est_single.mixing_time_extrapolated / delta)
    println("\n  (Single-exp would predict: $single_steps steps, offset C = $(round(est_single.offset; sigdigits=4)))")
end

# ══════════════════════════════════════════════════════════════════════════════
# PHASE 2: VERIFICATION — run exactly the predicted number of steps
# ══════════════════════════════════════════════════════════════════════════════
println("\n" * "="^70)
println("PHASE 2: VERIFICATION RUN  ($predicted_steps steps, δ = $delta)")
println("="^70)

verify_mixing_time = predicted_steps * delta
verify_save_every  = max(1, div(predicted_steps, 1000))

config_verify = Config(;
    sim                       = Thermalize(),
    domain                    = EnergyDomain(),
    construction              = KMS(),
    num_qubits                = num_qubits,
    with_linear_combination   = true,
    beta                      = beta,
    sigma                     = sigma,
    a                         = beta / 30.0,
    b                         = 0.4,
    num_energy_bits           = 12,
    w0                        = 0.05,
    t0                        = 2π / (2^12 * 0.05),
    num_trotter_steps_per_t0  = 10,
    mixing_time               = verify_mixing_time,
    delta                     = delta,
)

println("  T = $(round(verify_mixing_time; sigdigits=5)) → $predicted_steps steps")
println("  save_every = $verify_save_every")
println("  Running...")
flush(stdout)

result_verify = run_thermalize(jumps, config_verify, hamiltonian;
    rng        = MersenneTwister(42),
    save_every = verify_save_every,
)

wall2 = result_verify.metadata[:wall_time_seconds]
final_dist = result_verify.trace_distances[end]
min_dist   = minimum(result_verify.trace_distances)

println("  Done in $(round(wall2; digits=1))s")
println("  Data points: $(length(result_verify.trace_distances))")

# ── Results ───────────────────────────────────────────────────────────────────
crossing_idx = findfirst(d -> d ≤ target, result_verify.trace_distances)

println("\n" * "="^70)
println("RESULTS")
println("="^70)
println()
println("  Target ε          = $target")
println("  Predicted steps   = $predicted_steps")
println("  Final trace dist  = $(round(final_dist; sigdigits=4))")
println("  Min trace dist    = $(round(min_dist; sigdigits=4))")
println()

if crossing_idx !== nothing
    crossing_time = result_verify.time_steps[crossing_idx]
    crossing_step = round(Int, crossing_time / delta)
    pct_err = abs(crossing_step - predicted_steps) / crossing_step * 100

    println("  ✓ Target REACHED at step $crossing_step  (t = $(round(crossing_time; sigdigits=5)))")
    println("    Bi-exp predicted $predicted_steps steps  →  error: $(round(pct_err; digits=2))%")

    if final_dist ≤ target
        println("    Final distance $(round(final_dist; sigdigits=4)) ≤ $target  ✓  PASS")
    else
        println("    Final distance $(round(final_dist; sigdigits=4)) > $target")
        println("    (target was reached at step $crossing_step but distance rose above it by end)")
    end
else
    println("  ✗ Target NOT reached after $predicted_steps steps")
    println("    Final distance: $(round(final_dist; sigdigits=4))")
    println("    Gap to target:  $(round(final_dist - target; sigdigits=2))")
end

println()
println("  Total wall time: $(round(wall1 + wall2; digits=1))s")
