"""
Trajectory Cross-Validation Tests (Phase 4)
Run via: julia --project test/trajectory_validation/run_trajectory_validation.jl

Tests trajectory-averaged density matrix against Liouvillian-based DM evolution.
NOT included in Pkg.test() -- separate test group.

Requirements validated:
  TVAL-02: EnergyDomain single-step cross-validation
  TVAL-03: TimeDomain single-step cross-validation
  TVAL-04: TrotterDomain (construction=KMS()) single-step cross-validation
  TVAL-06: TrotterDomain coherent convergence to Gibbs state (multi-step)
"""

using Test
using LinearAlgebra
using Random

# Load shared fixtures
include(joinpath(@__DIR__, "..", "test_helpers.jl"))

# ---------------------------------------------------------------------------
# Single-step cross-validation helper
# ---------------------------------------------------------------------------
"""
    single_step_crossval(domain, delta; with_coherent, ntraj, seed) -> Float64

Run single-step trajectory-vs-DM cross-validation for a given domain and delta.

Returns trace distance between trajectory-averaged rho and DM rho.
"""
function single_step_crossval(domain, delta::Float64;
    with_coherent::Bool = (domain isa TrotterDomain),
    ntraj::Int = 50_000,
    seed::Int = 42,
)
    dim = SMALL_DIM  # 8

    # 1. Build Liouvillian and compute exact DM reference
    liouv_config = make_small_liouv_config(domain; construction=with_coherent ? KMS() : GNS())
    trotter_kw = domain isa TrotterDomain ? (; trotter=SMALL_TROTTER) : (;)
    jumps = domain isa TrotterDomain ? SMALL_TROTTER_JUMPS : SMALL_JUMPS
    L = construct_lindbladian(jumps, liouv_config, SMALL_HAM; trotter_kw...)

    psi0 = fill(ComplexF64(1.0), dim) / sqrt(dim)
    rho0 = psi0 * psi0'
    rho_dm = reshape(exp(delta * L) * vec(rho0), dim, dim)
    rho_dm = (rho_dm + rho_dm') / 2

    # 2. Build trajectory framework
    therm_config = make_small_thermalize_config(domain;
        construction=with_coherent ? KMS() : GNS(), delta=delta, mixing_time=Float64(delta))
    ham_or_trott = domain isa TrotterDomain ? SMALL_TROTTER : SMALL_HAM
    precomputed = QuantumFurnace._precompute_data(therm_config, ham_or_trott)
    scratch = QuantumFurnace.ThermalizeScratch(ComplexF64, dim)
    fw = build_trajectoryframework(jumps, ham_or_trott, therm_config,
        precomputed, scratch, delta)

    # 3. Run trajectories and accumulate rho
    rho_traj = zeros(ComplexF64, dim, dim)
    ws = QuantumFurnace.TrajectoryWorkspace(fw)
    rng = Random.Xoshiro(seed)
    for _ in 1:ntraj
        psi = copy(psi0)
        step_along_trajectory!(psi, fw, ws, rng)
        rho_traj .+= psi * psi'
    end
    rho_traj ./= ntraj
    rho_traj = (rho_traj + rho_traj') / 2
    rho_traj ./= tr(rho_traj)  # normalize

    # 4. Compute trace distance
    return trace_distance_h(Hermitian(rho_dm), Hermitian(rho_traj))
end

# ---------------------------------------------------------------------------
# TVAL-02: EnergyDomain single-step cross-validation
# ---------------------------------------------------------------------------
@testset "TVAL-02: EnergyDomain single-step cross-validation" begin
    deltas = [0.2, 0.1, 0.05]
    errors = Float64[]

    for delta in deltas
        dist = single_step_crossval(EnergyDomain(), delta)
        @test dist < 0.01
        push!(errors, dist)
        println("  delta=$delta  trace_dist=$(round(dist; sigdigits=4))")
    end

    # Delta^2 scaling: consecutive ratios should be ~4 (halving delta -> quarter error)
    ratios = [errors[i] / errors[i+1] for i in 1:length(errors)-1]
    for (i, ratio) in enumerate(ratios)
        println("  ratio[$(i)]: $(round(ratio; sigdigits=3))")
        @test 2.0 <= ratio <= 8.0
    end

    # Diagnostic: log-log slope (print only, no assertion)
    log_d = log.(deltas)
    log_e = log.(errors)
    n = length(deltas)
    slope = (n * sum(log_d .* log_e) - sum(log_d) * sum(log_e)) /
            (n * sum(log_d.^2) - sum(log_d)^2)
    println("  log-log slope: $(round(slope; sigdigits=3)) (expect ~2.0)")
end

# ---------------------------------------------------------------------------
# TVAL-03: TimeDomain single-step cross-validation
# ---------------------------------------------------------------------------
@testset "TVAL-03: TimeDomain single-step cross-validation" begin
    deltas = [0.2, 0.1, 0.05]
    errors = Float64[]

    for delta in deltas
        dist = single_step_crossval(TimeDomain(), delta)
        @test dist < 0.01
        push!(errors, dist)
        println("  delta=$delta  trace_dist=$(round(dist; sigdigits=4))")
    end

    # Delta^2 scaling: consecutive ratios should be ~4 (halving delta -> quarter error)
    ratios = [errors[i] / errors[i+1] for i in 1:length(errors)-1]
    for (i, ratio) in enumerate(ratios)
        println("  ratio[$(i)]: $(round(ratio; sigdigits=3))")
        @test 2.0 <= ratio <= 8.0
    end

    # Diagnostic: log-log slope (print only, no assertion)
    log_d = log.(deltas)
    log_e = log.(errors)
    n = length(deltas)
    slope = (n * sum(log_d .* log_e) - sum(log_d) * sum(log_e)) /
            (n * sum(log_d.^2) - sum(log_d)^2)
    println("  log-log slope: $(round(slope; sigdigits=3)) (expect ~2.0)")
end

# ---------------------------------------------------------------------------
# TVAL-04: TrotterDomain (construction=KMS()) single-step cross-validation
# ---------------------------------------------------------------------------
@testset "TVAL-04: TrotterDomain single-step cross-validation" begin
    deltas = [0.2, 0.1, 0.05]
    errors = Float64[]

    for delta in deltas
        dist = single_step_crossval(TrotterDomain(), delta; construction=KMS())
        @test dist < 0.01
        push!(errors, dist)
        println("  delta=$delta  trace_dist=$(round(dist; sigdigits=4))")
    end

    # Delta^2 scaling: consecutive ratios should be ~4 (halving delta -> quarter error)
    ratios = [errors[i] / errors[i+1] for i in 1:length(errors)-1]
    for (i, ratio) in enumerate(ratios)
        println("  ratio[$(i)]: $(round(ratio; sigdigits=3))")
        @test 2.0 <= ratio <= 8.0
    end

    # Diagnostic: log-log slope (print only, no assertion)
    log_d = log.(deltas)
    log_e = log.(errors)
    n = length(deltas)
    slope = (n * sum(log_d .* log_e) - sum(log_d) * sum(log_e)) /
            (n * sum(log_d.^2) - sum(log_d)^2)
    println("  log-log slope: $(round(slope; sigdigits=3)) (expect ~2.0)")
end

# ---------------------------------------------------------------------------
# TVAL-06: TrotterDomain coherent convergence to Gibbs state (multi-step)
# ---------------------------------------------------------------------------
@testset "TVAL-06: TrotterDomain coherent convergence to Gibbs state" begin
    dim = SMALL_DIM  # 8
    delta = 0.01     # Small enough for convergence, large enough for speed

    # Liouvillian for TrotterDomain + coherent
    liouv_config = make_small_liouv_config(TrotterDomain(); construction=KMS())
    L = construct_lindbladian(SMALL_TROTTER_JUMPS, liouv_config, SMALL_HAM; trotter=SMALL_TROTTER)
    exp_L = exp(delta * L)  # Compute once, apply repeatedly

    # Gibbs state in Trotter eigenbasis (following DMTST-02 pattern)
    gibbs_comp = SMALL_HAM.eigvecs * SMALL_GIBBS * SMALL_HAM.eigvecs'
    gibbs_trott = Hermitian(SMALL_TROTTER.eigvecs' * gibbs_comp * SMALL_TROTTER.eigvecs)

    # Liouvillian fixed point (exact steady state of L).
    # The TrotterDomain Liouvillian's fixed point has a negligible domain approximation
    # offset (~1e-6) from the true Gibbs state after the Trotter basis transform fix.
    # We compare against the fixed point for DM convergence (deterministic), and
    # validate trajectory convergence toward the Gibbs state region with a threshold
    # that accounts for statistical noise from finite trajectory count.
    eig = eigen(L)
    ss_idx = argmin(abs.(eig.values))
    ss_vec = eig.vectors[:, ss_idx]
    ss_dm = reshape(ss_vec, dim, dim)
    ss_dm = (ss_dm + ss_dm') / 2
    ss_dm ./= tr(ss_dm)

    # Sanity check: fixed point is close to Gibbs (domain approximation error)
    fp_gibbs_dist = trace_distance_h(Hermitian(ss_dm), gibbs_trott)
    println("  Fixed point -> Gibbs trace distance: $(round(fp_gibbs_dist; sigdigits=4))")
    @test fp_gibbs_dist < 1e-4  # Domain approximation error is negligible after basis fix

    # -- DM evolution via iterated matrix-vector multiply --
    psi0 = fill(ComplexF64(1.0), dim) / sqrt(dim)
    rho_dm_vec = vec(psi0 * psi0')

    max_steps = 5000  # Safety cap
    dm_converged_step = 0
    for step in 1:max_steps
        rho_dm_vec = exp_L * rho_dm_vec
        rho_dm = reshape(copy(rho_dm_vec), dim, dim)
        rho_dm = (rho_dm + rho_dm') / 2
        dist = trace_distance_h(Hermitian(rho_dm), Hermitian(ss_dm))
        if dist < 1e-3
            dm_converged_step = step
            println("  DM converged at step $step (trace_dist=$(round(dist; sigdigits=4)))")
            break
        end
    end
    @test dm_converged_step > 0  # DM must converge within max_steps

    # -- Trajectory evolution for the same total time --
    total_time = dm_converged_step * delta
    therm_config = make_small_thermalize_config(TrotterDomain();
        construction=KMS(), delta=delta, mixing_time=total_time)

    result = run_trajectories(SMALL_TROTTER_JUMPS, therm_config, psi0, SMALL_HAM;
        trotter=SMALL_TROTTER, total_time=total_time, delta=delta, ntraj=10_000, seed=42)
    rho_traj = result.rho_mean

    # -- Assertions --
    # DM converges to the Liouvillian fixed point within 1e-3
    dist_dm = trace_distance_h(
        Hermitian(reshape(copy(rho_dm_vec), dim, dim) |> m -> (m+m')/2), Hermitian(ss_dm))
    @test dist_dm < 1e-3     # DM converged to fixed point

    # Trajectory-averaged rho converges toward the Gibbs state region.
    # With 10,000 trajectories, the 1/sqrt(N) statistical noise floor is ~0.01.
    # Domain approximation is now negligible (~1e-6), so the dominant error source is
    # purely statistical noise. Total expected distance is ~0.01. We assert < 0.015 (safe margin).
    dist_traj_gibbs = trace_distance_h(Hermitian(rho_traj), gibbs_trott)
    @test dist_traj_gibbs < 0.015  # Trajectory near Gibbs (statistical noise dominant)

    # Trajectory must also be closer to the fixed point than the initial state was
    # (confirming thermalization happened, not just noise)
    rho0 = psi0 * psi0'
    dist_init = trace_distance_h(Hermitian(rho0), Hermitian(ss_dm))
    dist_traj = trace_distance_h(Hermitian(rho_traj), Hermitian(ss_dm))
    @test dist_traj < dist_init / 10  # Trajectory moved substantially toward fixed point

    # Report all distances for diagnostics
    println("  DM -> fixed point trace distance:     $(round(dist_dm; sigdigits=4))")
    println("  Traj -> fixed point trace distance:    $(round(dist_traj; sigdigits=4))")
    println("  Traj -> Gibbs trace distance:          $(round(dist_traj_gibbs; sigdigits=4))")
    println("  Initial -> fixed point trace distance: $(round(dist_init; sigdigits=4))")
    println("  Steps used: $dm_converged_step (total_time=$(round(total_time; sigdigits=4)))")
end
