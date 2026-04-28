#!/usr/bin/env julia
#
# Two-phase mixing time demo with δ = 0.001:
#   Phase 1: Run a moderate-length simulation, extrapolate steps needed for ε = 1e-4
#   Phase 2: Run the full simulation for that many steps and verify convergence
#
# With coherent unitary, floor ≈ 0.068 × δ = 6.8e-5, safely below target 1e-4.
# Need T ≈ 30 time units (gap ≈ 0.32) → 30,000 steps.
#
# Usage:  OPENBLAS_NUM_THREADS=4 julia -t4 --project scripts/mixing_time_extrapolate_verify.jl

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

println("\n" * "="^70)
println("PHASE 1: EXTRAPOLATION RUN  (δ = $delta)")
println("="^70)

# Run for T=15 (~5 time constants at gap≈0.32), enough to clearly fit the decay
# 15 / 0.001 = 15,000 steps; save_every=10 → 1500 data points
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
    s                         = 0.4,
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

# Fit and extrapolate — single-exponential (original)
est = estimate_mixing_time(result_extrap;
    model          = :single,
    skip_initial   = 0.2,
    target_epsilon = target,
    extrapolate    = true,
)

println("\n--- Single-exponential fit ---")
println("Fit: d(t) = $(round(est.amplitude; sigdigits=4)) × exp(-$(round(est.fitted_gap; sigdigits=4)) t) + $(round(est.offset; sigdigits=4))")
println("  R² = $(round(est.r_squared; sigdigits=5))")
println("  Gap CI (95%): $(round.(est.gap_ci; sigdigits=4))")
println("  Converged: $(est.converged)")
println("  Offset vs target: C = $(round(est.offset; sigdigits=4)) vs ε = $target")

if est.mixing_time_actual !== nothing
    println("\n  Actual crossing at t = $(round(est.mixing_time_actual; sigdigits=5))")
    println("    → $(round(Int, est.mixing_time_actual / delta)) delta steps")
end

if est.mixing_time_extrapolated !== nothing && !isnan(est.mixing_time_extrapolated)
    println("  Extrapolated crossing at t = $(round(est.mixing_time_extrapolated; sigdigits=5))")
    println("    → $(ceil(Int, est.mixing_time_extrapolated / delta)) delta steps")
else
    println("\n  WARNING: Single-exp extrapolation failed!")
    println("  Fitted offset C=$(round(est.offset; sigdigits=4)) vs target ε=$(target)")
end

single_steps = if est.mixing_time_extrapolated !== nothing && !isnan(est.mixing_time_extrapolated)
    ceil(Int, est.mixing_time_extrapolated / delta)
else
    nothing
end

# Fit and extrapolate — bi-exponential (Phase 43)
est_bi = estimate_mixing_time(result_extrap;
    model          = :biexp,
    skip_initial   = 0.2,
    target_epsilon = target,
    extrapolate    = true,
)

println("\n--- Bi-exponential fit ---")
bifit = est_bi.biexp_fit_result
println("Fit: d(t) = $(round(bifit.amplitude_fast; sigdigits=4)) × exp(-$(round(bifit.gap_fast; sigdigits=4)) t) + $(round(bifit.amplitude; sigdigits=4)) × exp(-$(round(bifit.gap; sigdigits=4)) t) + $(round(bifit.offset; sigdigits=4))")
println("  R² = $(round(est_bi.r_squared; sigdigits=5))")
println("  Slow gap (spectral gap): $(round(est_bi.fitted_gap; sigdigits=4))")
println("  Fast gap: $(round(bifit.gap_fast; sigdigits=4))")
println("  Converged: $(est_bi.converged)")
println("  Offset vs target: C = $(round(est_bi.offset; sigdigits=4)) vs ε = $target")

if est_bi.mixing_time_extrapolated !== nothing && !isnan(est_bi.mixing_time_extrapolated)
    println("  Extrapolated crossing at t = $(round(est_bi.mixing_time_extrapolated; sigdigits=5))")
    println("    → $(ceil(Int, est_bi.mixing_time_extrapolated / delta)) delta steps")
else
    println("\n  WARNING: Bi-exp extrapolation failed!")
end

# Use bi-exponential prediction for the verification run (more accurate)
biexp_steps = if est_bi.mixing_time_extrapolated !== nothing && !isnan(est_bi.mixing_time_extrapolated)
    ceil(Int, est_bi.mixing_time_extrapolated / delta)
else
    nothing
end

# Choose the best prediction available
needed_steps = if biexp_steps !== nothing
    biexp_steps
elseif single_steps !== nothing
    single_steps
else
    println("\nERROR: Both extrapolations failed!")
    exit(1)
end
needed_time = needed_steps * delta

println("\n--- Comparison ---")
if single_steps !== nothing
    println("  Single-exp prediction: $single_steps steps")
end
if biexp_steps !== nothing
    println("  Bi-exp prediction:     $biexp_steps steps")
end
println("  Using bi-exp prediction: $needed_steps steps (physical time = $(round(needed_time; sigdigits=5)))")

# ── 4. PHASE 2: VERIFICATION RUN ─────────────────────────────────────────────
println("\n" * "="^70)
println("PHASE 2: VERIFICATION RUN  ($needed_steps steps, δ = $delta)")
println("="^70)

# Use the extrapolated mixing time + 20% margin (to see the floor clearly)
verify_mixing_time = needed_time * 1.2
verify_num_steps   = round(Int, verify_mixing_time / delta)
# ~1000 data points
verify_save_every  = max(1, div(verify_num_steps, 1000))

config_verify = Config(;
    sim                       = Thermalize(),
    domain                    = EnergyDomain(),
    construction              = KMS(),
    num_qubits                = num_qubits,
    with_linear_combination   = true,
    beta                      = beta,
    sigma                     = sigma,
    a                         = beta / 30.0,
    s                         = 0.4,
    num_energy_bits           = 12,
    w0                        = 0.05,
    t0                        = 2π / (2^12 * 0.05),
    num_trotter_steps_per_t0  = 10,
    mixing_time               = verify_mixing_time,
    delta                     = delta,
)

println("  T = $(round(verify_mixing_time; sigdigits=5)) → $verify_num_steps steps")
println("  save_every = $verify_save_every")
println("  Running...")
flush(stdout)

result_verify = run_thermalize(jumps, config_verify, hamiltonian;
    rng        = MersenneTwister(42),
    save_every = verify_save_every,
)

wall2 = result_verify.metadata[:wall_time_seconds]
final_dist = result_verify.trace_distances[end]
println("  Done in $(round(wall2; digits=1))s")
println("  Data points: $(length(result_verify.trace_distances))")

# Find when we actually cross the target
crossing_idx = findfirst(d -> d ≤ target, result_verify.trace_distances)

println("\n" * "="^70)
println("VERIFICATION RESULTS")
println("="^70)
println()
println("  Target ε = $target")
println("  Final trace distance = $(round(final_dist; sigdigits=4))")
println()

if crossing_idx !== nothing
    crossing_time = result_verify.time_steps[crossing_idx]
    crossing_step = round(Int, crossing_time / delta)
    println("  Target REACHED at t = $(round(crossing_time; sigdigits=5)) (step $crossing_step)")

    if single_steps !== nothing
        single_pct_err = abs(crossing_step - single_steps) / crossing_step * 100
        println("    Single-exp predicted: $single_steps steps  (error: $(round(single_pct_err; digits=1))%)")
    end
    if biexp_steps !== nothing
        biexp_pct_err = abs(crossing_step - biexp_steps) / crossing_step * 100
        println("    Bi-exp predicted:     $biexp_steps steps  (error: $(round(biexp_pct_err; digits=1))%)")
    end
    println("    Actual:               $crossing_step steps")
else
    println("  Target NOT reached after $verify_num_steps steps")
    println("    Final distance: $(round(final_dist; sigdigits=4))")
    println("    Predicted $needed_steps steps would suffice -- but it didn't")
end

println()
println("  Minimum trace distance seen: $(round(minimum(result_verify.trace_distances); sigdigits=4))")
println("  Total wall time: $(round(wall1 + wall2; digits=1))s")
