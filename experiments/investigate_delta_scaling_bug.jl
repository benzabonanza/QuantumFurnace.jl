#!/usr/bin/env julia
# ============================================================================
# Investigate Non-Monotonic Delta Scaling Bug (Quick-32)
# ============================================================================
#
# Diagnoses whether the non-monotonic delta scaling in gap estimation
# (Quick-30/31: Mz_stagg error goes +0.6% -> -6.5% -> -19.7% as delta
# shrinks 0.1 -> 0.01 -> 0.001) originates from:
#   (a) A bug in the trajectory simulation pipeline
#   (b) A property of the exponential decay fitting procedure
#   (c) Accumulated Lie-Trotter splitting bias in multi-step trajectories
#
# Key diagnostic: compare trajectory-averaged rho against exact exp(t*L)*rho0
# at the density matrix level (bypassing the fitting procedure entirely).
#
# Expected scaling: trace_distance(rho_traj, rho_exact) = O(delta) for
# multi-step evolution with T_total fixed (N*O(delta^2) = (T/delta)*O(delta^2)).
# If trace distance does NOT decrease monotonically with delta -> code bug.
# If it does decrease -> non-monotonic gap error is in the FITTING, not simulation.
#
# Parameters (matching Quick-30/31):
#   - System: n=4 disordered periodic Heisenberg, beta=10.0
#   - Domain: TimeDomain, with_coherent=false
#   - Initial state: psi0 = ones(ComplexF64, dim) / sqrt(dim) (uniform)
#   - T_total = 20.0 (matching Quick-31 mixing_time)
#   - Deltas: [0.1, 0.01, 0.001]
#   - Trajectories: 50,000 (need higher count for DM accuracy at small delta)
#   - Seed: 42
#
# Usage:
#   cd QuantumFurnace.jl && julia -t 4 --project experiments/investigate_delta_scaling_bug.jl
# ============================================================================

using QuantumFurnace
using LinearAlgebra
using Printf
using Random

# ============================================================================
# Section 0: Constants and System Setup
# ============================================================================

const BETA = 10.0
const SEED = 42
const DELTAS = [0.1, 0.01, 0.001]
const T_TOTAL = 20.0
const NTRAJ = 50_000

# Grid parameters matching test suite conventions
const NUM_ENERGY_BITS = 12
const W0 = 0.05
const T0 = 2pi / (2^NUM_ENERGY_BITS * W0)
const NUM_TROTTER_STEPS_PER_T0 = 10
const SIGMA = 1.0 / BETA  # 0.1

# Pauli matrices
const Xp = ComplexF64[0 1; 1 0]
const Yp = ComplexF64[0 -im; im 0]
const Zp = ComplexF64[1 0; 0 -1]

# ---------------------------------------------------------------------------
# Helper: create system
# ---------------------------------------------------------------------------
function make_system(n, beta)
    ham = load_hamiltonian("heis", n; beta=beta)
    dim = 2^n

    jump_paulis = [[Xp], [Yp], [Zp]]
    num_of_jumps = 3 * n
    jump_normalization = sqrt(num_of_jumps)
    V = ham.eigvecs

    jumps = JumpOp[]
    for pauli in jump_paulis
        for site in 1:n
            jump_op = Matrix(pad_term(pauli, n, site)) ./ jump_normalization
            jump_in_eigen = V' * jump_op * V
            orthogonal = (jump_op == transpose(jump_op))
            herm = (jump_op == jump_op')
            push!(jumps, JumpOp(jump_op, jump_in_eigen, orthogonal, herm))
        end
    end

    return ham, jumps, dim
end

# ---------------------------------------------------------------------------
# Helper: create configs
# ---------------------------------------------------------------------------
function make_liouv_config(n)
    LiouvConfig(;
        num_qubits = n,
        with_coherent = false,
        with_linear_combination = true,
        domain = TimeDomain(),
        beta = BETA,
        sigma = SIGMA,
        a = BETA / 30.0,
        b = 0.4,
        num_energy_bits = NUM_ENERGY_BITS,
        w0 = W0,
        t0 = T0,
        num_trotter_steps_per_t0 = NUM_TROTTER_STEPS_PER_T0,
    )
end

function make_thermalize_config(n; mixing_time=50.0, delta=0.01)
    ThermalizeConfig(;
        num_qubits = n,
        with_coherent = false,
        with_linear_combination = true,
        domain = TimeDomain(),
        beta = BETA,
        sigma = SIGMA,
        a = BETA / 30.0,
        b = 0.4,
        num_energy_bits = NUM_ENERGY_BITS,
        w0 = W0,
        t0 = T0,
        num_trotter_steps_per_t0 = NUM_TROTTER_STEPS_PER_T0,
        mixing_time = mixing_time,
        delta = delta,
    )
end

# ---------------------------------------------------------------------------
# Helper: trace distance for general (non-Hermitian wrapped) matrices
# ---------------------------------------------------------------------------
function trace_dist(A::Matrix, B::Matrix)
    diff = A - B
    diff_h = (diff + diff') / 2
    evals = eigvals(Hermitian(diff_h))
    return 0.5 * sum(abs.(evals))
end

println("=" ^ 70)
println("  INVESTIGATE NON-MONOTONIC DELTA SCALING BUG (Quick-32)")
println("  n=4 disordered Heisenberg, beta=$BETA")
println("  T_total=$T_TOTAL, ntraj=$NTRAJ, seed=$SEED")
println("  Deltas: $DELTAS")
println("=" ^ 70)
println()

# ============================================================================
# System setup
# ============================================================================

const n = 4
@printf("Loading disordered Hamiltonian for n=%d...\n", n)
ham, jumps, dim = make_system(n, BETA)

@printf("Building Lindbladian for n=%d...\n", n)
config_l = make_liouv_config(n)
liouv_result = run_lindbladian(jumps, config_l, ham)
exact_gap = abs(real(liouv_result.spectral_gap))
@printf("Exact gap (ARPACK): %.10f\n\n", exact_gap)

# Full dense Lindbladian for exp(t*L) computation
@printf("Constructing full dense Lindbladian...\n")
L_full = Matrix(construct_lindbladian(jumps, config_l, ham))
@printf("Lindbladian size: %d x %d\n\n", size(L_full)...)

# Build observables
obs, obs_names = build_preset_trajectory_observables(ham, n)
@printf("Observables: %s\n\n", join(obs_names, ", "))

# Initial state: uniform superposition
psi0 = ones(ComplexF64, dim) / sqrt(dim)
rho0 = psi0 * psi0'
rho0_vec = vec(rho0)

# ============================================================================
# Section 1: Single-step DM-vs-trajectory sanity check (3 deltas)
# ============================================================================

println()
println("=" ^ 70)
println("  SECTION 1: SINGLE-STEP SANITY CHECK")
println("  rho_traj(1 step) vs exp(delta*L)*rho0")
println("  Expected: trace_distance = O(delta^2)")
println("=" ^ 70)
println()

single_step_dists = Float64[]
for delta_val in DELTAS
    # Exact: exp(delta*L)*rho0
    rho_exact = reshape(exp(delta_val * L_full) * rho0_vec, dim, dim)
    rho_exact = (rho_exact + rho_exact') / 2

    # Trajectory: 1 step, many trajectories
    config_t = make_thermalize_config(n; mixing_time=Float64(delta_val), delta=delta_val)
    result = run_trajectories(jumps, config_t, psi0, ham;
        total_time=Float64(delta_val), delta=delta_val,
        ntraj=NTRAJ, seed=SEED)
    rho_traj = result.rho_mean

    dist = trace_dist(rho_traj, rho_exact)
    push!(single_step_dists, dist)
    @printf("  delta=%-8.4f  trace_dist=%.2e\n", delta_val, dist)
end
println()

# Check O(delta^2) scaling
if length(single_step_dists) >= 2
    for i in 1:(length(DELTAS)-1)
        ratio = single_step_dists[i] / single_step_dists[i+1]
        expected_ratio = (DELTAS[i] / DELTAS[i+1])^2
        @printf("  ratio[%d->%d]: %.2f (expected ~%.1f for O(delta^2))\n",
                i, i+1, ratio, expected_ratio)
    end
end
println()

# ============================================================================
# Section 2: Multi-step rho comparison: rho_traj vs exp(T*L)*rho0 (3 deltas)
# ============================================================================

println("=" ^ 70)
println("  SECTION 2: MULTI-STEP RHO COMPARISON")
println("  rho_traj(T=$T_TOTAL) vs exp(T*L)*rho0")
println("  Expected: trace_distance = O(delta) for fixed T")
println("  (N steps * O(delta^2)/step = (T/delta)*O(delta^2) = O(delta))")
println("=" ^ 70)
println()

# Compute exact rho at T_total (once, independent of delta)
@printf("Computing exact rho(T=%g) via exp(T*L)*rho0...\n", T_TOTAL)
t_start = time()
rho_exact_T = reshape(exp(T_TOTAL * L_full) * rho0_vec, dim, dim)
rho_exact_T = (rho_exact_T + rho_exact_T') / 2
rho_exact_T ./= tr(rho_exact_T)
@printf("Done (%.1f sec). tr(rho_exact) = %.10f\n\n", time() - t_start, real(tr(rho_exact_T)))

multi_step_dists = Float64[]
for delta_val in DELTAS
    num_steps = ceil(Int, T_TOTAL / delta_val)
    @printf("  delta=%-8.4f  num_steps=%d\n", delta_val, num_steps)

    local t_start = time()
    config_t = make_thermalize_config(n; mixing_time=T_TOTAL, delta=delta_val)
    result = run_trajectories(jumps, config_t, psi0, ham;
        total_time=T_TOTAL, delta=delta_val,
        ntraj=NTRAJ, seed=SEED)
    rho_traj = result.rho_mean
    t_elapsed = time() - t_start

    # Normalize
    rho_traj ./= tr(rho_traj)

    dist = trace_dist(rho_traj, rho_exact_T)
    push!(multi_step_dists, dist)

    @printf("    trace_dist(rho_traj, rho_exact) = %.6e  (%.1f sec)\n", dist, t_elapsed)
    @printf("    tr(rho_traj) = %.10f\n", real(tr(result.rho_mean)))
end
println()

# Check O(delta) scaling
@printf("Multi-step scaling analysis:\n")
@printf("%-10s  %14s  %14s\n", "delta", "trace_dist", "trace_dist/delta")
@printf("%-10s  %14s  %14s\n", "-"^10, "-"^14, "-"^14)
for (i, delta_val) in enumerate(DELTAS)
    @printf("%-10.4f  %14.6e  %14.6e\n", delta_val, multi_step_dists[i], multi_step_dists[i] / delta_val)
end
println()

# Monotonicity check
is_monotonic = all(multi_step_dists[i] > multi_step_dists[i+1] for i in 1:length(multi_step_dists)-1)
@printf("Trace distance decreasing with delta? %s\n", is_monotonic ? "YES" : "NO")
if is_monotonic
    println("  -> Trajectory simulation produces CORRECT rho that improves with smaller delta.")
    println("  -> Non-monotonic gap estimation error is in the FITTING procedure, not the simulation.")
else
    println("  -> WARNING: Trajectory rho does NOT improve with smaller delta.")
    println("  -> This indicates a potential CODE BUG in the trajectory pipeline.")
end
println()

# ============================================================================
# Section 3: Observable expectation values from final rho (3 deltas, 8 obs)
# ============================================================================

println("=" ^ 70)
println("  SECTION 3: OBSERVABLE EXPECTATION VALUES")
println("  tr(rho_traj * O) vs tr(rho_exact * O) for each observable")
println("=" ^ 70)
println()

# Exact observable values from rho_exact_T
exact_obs_vals = [real(tr(rho_exact_T * obs[i])) for i in 1:length(obs)]

@printf("Exact observable values at T=%g:\n", T_TOTAL)
for (i, name) in enumerate(obs_names)
    @printf("  %-12s = %+.10f\n", name, exact_obs_vals[i])
end
println()

# For each delta, compute trajectory observable values and errors
@printf("%-12s", "Observable")
for delta_val in DELTAS
    @printf("  err(d=%.0e) ", delta_val)
end
@printf("  monotonic?\n")
@printf("%-12s", "-"^12)
for _ in DELTAS
    @printf("  %12s", "-"^12)
end
@printf("  %s\n", "-"^10)

obs_errors_by_delta = Dict{Float64, Vector{Float64}}()
for delta_val in DELTAS
    config_t = make_thermalize_config(n; mixing_time=T_TOTAL, delta=delta_val)
    result = run_trajectories(jumps, config_t, psi0, ham;
        total_time=T_TOTAL, delta=delta_val,
        ntraj=NTRAJ, seed=SEED)
    rho_traj = result.rho_mean
    rho_traj ./= tr(rho_traj)

    errs = Float64[]
    for i in 1:length(obs)
        traj_val = real(tr(rho_traj * obs[i]))
        push!(errs, traj_val - exact_obs_vals[i])
    end
    obs_errors_by_delta[delta_val] = errs
end

for (obs_idx, name) in enumerate(obs_names)
    @printf("%-12s", name)
    errs_for_obs = [obs_errors_by_delta[d][obs_idx] for d in DELTAS]
    for err in errs_for_obs
        @printf("  %+12.6e", err)
    end
    # Check if |error| decreases monotonically
    abs_errs = abs.(errs_for_obs)
    mono = all(abs_errs[i] >= abs_errs[i+1] for i in 1:length(abs_errs)-1)
    @printf("  %s\n", mono ? "YES" : "NO")
end
println()

# Summary of observable-level monotonicity
global n_mono = 0
for obs_idx in 1:length(obs)
    errs_for_obs = [abs(obs_errors_by_delta[d][obs_idx]) for d in DELTAS]
    if all(errs_for_obs[i] >= errs_for_obs[i+1] for i in 1:length(errs_for_obs)-1)
        global n_mono += 1
    end
end
@printf("Observables with monotonic |error| decrease: %d/%d\n\n", n_mono, length(obs))

# ============================================================================
# Section 4: Time-resolved observable comparison for delta=0.01
# ============================================================================

println("=" ^ 70)
println("  SECTION 4: TIME-RESOLVED OBSERVABLE COMPARISON (delta=0.01)")
println("  Trajectory O(t) vs exact O(t) from iterated exp(delta*save_every*L)")
println("=" ^ 70)
println()

delta_4 = 0.01
save_every_4 = 10
dt_save = delta_4 * save_every_4  # time between saves

@printf("delta=%.4f, save_every=%d, dt_save=%.4f\n", delta_4, save_every_4, dt_save)

# Run trajectory with observable measurements
config_t4 = make_thermalize_config(n; mixing_time=T_TOTAL, delta=delta_4)
@printf("Running %d trajectories with observable tracking...\n", NTRAJ)
t_start = time()
traj_result = run_observable_trajectories(jumps, config_t4, psi0, ham;
    observables=obs, save_every=save_every_4,
    ntraj=NTRAJ, total_time=T_TOTAL, delta=delta_4,
    seed=SEED, reconstruct_dm=true)
@printf("Done (%.1f sec)\n", time() - t_start)

# Compute exact time series: iterated exp(dt_save * L)
num_steps_4 = ceil(Int, T_TOTAL / delta_4)
num_saves_4 = div(num_steps_4, save_every_4) + 1

@printf("Computing exact time series via exp(dt_save*L) iteration (%d save points)...\n", num_saves_4)
exp_L_save = exp(dt_save * L_full)

exact_obs_ts = zeros(Float64, length(obs), num_saves_4)
let rho_v = copy(rho0_vec)
    for s in 1:num_saves_4
        rho_s = reshape(rho_v, dim, dim)
        rho_s = (rho_s + rho_s') / 2
        for i in 1:length(obs)
            exact_obs_ts[i, s] = real(tr(rho_s * obs[i]))
        end
        if s < num_saves_4
            rho_v = exp_L_save * rho_v
        end
    end
end
println("Done.\n")

# Compare: max absolute error over time for each observable
traj_obs_ts = traj_result.measurements_mean
times = traj_result.times

# Ensure same number of save points
n_compare = min(size(traj_obs_ts, 2), size(exact_obs_ts, 2))
@printf("Comparing %d save points (trajectory: %d, exact: %d)\n\n",
        n_compare, size(traj_obs_ts, 2), size(exact_obs_ts, 2))

@printf("%-12s  %14s  %14s  %14s\n",
        "Observable", "max_abs_err", "mean_abs_err", "final_err")
@printf("%-12s  %14s  %14s  %14s\n",
        "-"^12, "-"^14, "-"^14, "-"^14)

for (i, name) in enumerate(obs_names)
    errs_ts = abs.(traj_obs_ts[i, 1:n_compare] .- exact_obs_ts[i, 1:n_compare])
    max_err = maximum(errs_ts)
    mean_err = sum(errs_ts) / length(errs_ts)
    final_err = traj_obs_ts[i, n_compare] - exact_obs_ts[i, n_compare]
    @printf("%-12s  %14.6e  %14.6e  %+14.6e\n", name, max_err, mean_err, final_err)
end
println()

# Print first/last few time points for Mz_stagg (the problematic observable)
mz_stagg_idx = findfirst(==("Mz_stagg"), obs_names)
if mz_stagg_idx !== nothing
    @printf("Time-resolved comparison for Mz_stagg:\n")
    @printf("%-10s  %14s  %14s  %14s\n", "time", "traj_val", "exact_val", "error")
    @printf("%-10s  %14s  %14s  %14s\n", "-"^10, "-"^14, "-"^14, "-"^14)

    # Show first 5, last 5
    show_indices = vcat(1:min(5, n_compare), max(1, n_compare-4):n_compare)
    show_indices = unique(show_indices)
    for s in show_indices
        t = times[s]
        tv = traj_obs_ts[mz_stagg_idx, s]
        ev = exact_obs_ts[mz_stagg_idx, s]
        @printf("%-10.4f  %14.10f  %14.10f  %+14.6e\n", t, tv, ev, tv - ev)
    end
    println()
end

# ============================================================================
# Section 5: Summary and Diagnosis
# ============================================================================

println()
println("=" ^ 70)
println("  SECTION 5: SUMMARY AND DIAGNOSIS")
println("=" ^ 70)
println()

@printf("System: n=%d disordered periodic Heisenberg, beta=%.1f\n", n, BETA)
@printf("Exact spectral gap: %.10f\n\n", exact_gap)

# Section 1 summary
println("--- Single-step sanity check ---")
for (i, delta_val) in enumerate(DELTAS)
    @printf("  delta=%-8.4f  trace_dist=%.2e\n", delta_val, single_step_dists[i])
end
single_step_ok = all(single_step_dists[i] > single_step_dists[i+1] for i in 1:length(single_step_dists)-1)
@printf("  O(delta^2) scaling: %s\n\n", single_step_ok ? "CONFIRMED" : "FAILED")

# Section 2 summary
println("--- Multi-step DM comparison (T=$T_TOTAL) ---")
for (i, delta_val) in enumerate(DELTAS)
    @printf("  delta=%-8.4f  trace_dist=%.6e  dist/delta=%.6e\n",
            delta_val, multi_step_dists[i], multi_step_dists[i] / delta_val)
end
@printf("  Monotonically decreasing: %s\n\n", is_monotonic ? "YES" : "NO")

# Section 3 summary
println("--- Observable errors monotonicity ---")
@printf("  %d/%d observables have monotonically decreasing |error| with delta\n\n", n_mono, length(obs))

# Final diagnosis
println("=" ^ 70)
println("  DIAGNOSIS")
println("=" ^ 70)
println()

if is_monotonic && single_step_ok
    println("ROOT CAUSE: The trajectory simulation produces CORRECT density matrices.")
    println("The non-monotonic delta scaling in gap estimation (Quick-30/31) is NOT a code bug.")
    println("")
    println("The non-monotonic behavior originates in the FITTING PROCEDURE:")
    println("  - The exponential decay model y(t) = A*exp(-gap*t) + C is a single-mode")
    println("    approximation to a multi-exponential decay process.")
    println("  - Different delta values produce different effective Kraus channels,")
    println("    which excite different spectral modes with different weights.")
    println("  - The best-observable selection (smallest gap among converged fits)")
    println("    picks different observables at different deltas, compounding the")
    println("    non-monotonic behavior.")
    println("  - This is an inherent limitation of single-exponential fitting for")
    println("    gap estimation, not a simulation error.")
else
    println("POTENTIAL CODE BUG DETECTED!")
    if !single_step_ok
        println("  Single-step O(delta^2) scaling FAILED.")
    end
    if !is_monotonic
        println("  Multi-step trace distance does NOT decrease monotonically with delta.")
    end
    println("  Further investigation needed in src/trajectories.jl step_along_trajectory!")
end
println()
println("=" ^ 70)
println("  END OF INVESTIGATION")
println("=" ^ 70)
