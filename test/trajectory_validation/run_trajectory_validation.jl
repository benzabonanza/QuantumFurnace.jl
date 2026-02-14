"""
Trajectory Cross-Validation Tests (Phase 4)
Run via: julia --project test/trajectory_validation/run_trajectory_validation.jl

Tests trajectory-averaged density matrix against Liouvillian-based DM evolution.
NOT included in Pkg.test() -- separate test group.

Requirements validated:
  TVAL-02: EnergyDomain single-step cross-validation
  TVAL-03: TimeDomain single-step cross-validation
  TVAL-04: TrotterDomain (with_coherent=true) single-step cross-validation
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
    liouv_config = make_small_liouv_config(domain; with_coherent=with_coherent)
    trotter_kw = domain isa TrotterDomain ? (; trotter=SMALL_TROTTER) : (;)
    L = construct_lindbladian(SMALL_JUMPS, liouv_config, SMALL_HAM; trotter_kw...)

    psi0 = fill(ComplexF64(1.0), dim) / sqrt(dim)
    rho0 = psi0 * psi0'
    rho_dm = reshape(exp(delta * L) * vec(rho0), dim, dim)
    rho_dm = (rho_dm + rho_dm') / 2

    # 2. Build trajectory framework
    therm_config = make_small_thermalize_config(domain;
        with_coherent=with_coherent, delta=delta, mixing_time=Float64(delta))
    ham_or_trott = domain isa TrotterDomain ? SMALL_TROTTER : SMALL_HAM
    precomputed = precompute_data(domain, therm_config, ham_or_trott)
    scratch = KrausScratch(ComplexF64, dim)
    fw = build_trajectoryframework(SMALL_JUMPS, ham_or_trott, therm_config,
        precomputed, scratch, delta)

    # 3. Run trajectories and accumulate rho
    rho_traj = zeros(ComplexF64, dim, dim)
    Random.seed!(seed)
    for _ in 1:ntraj
        psi = copy(psi0)
        step_along_trajectory!(psi, fw)
        rho_traj .+= psi * psi'
    end
    rho_traj ./= ntraj
    rho_traj = (rho_traj + rho_traj') / 2
    rho_traj ./= tr(rho_traj)  # normalize

    # 4. Compute trace distance
    return QuantumFurnace.trace_distance_h(Hermitian(rho_dm), Hermitian(rho_traj))
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
# TVAL-04: TrotterDomain (with_coherent=true) single-step cross-validation
# ---------------------------------------------------------------------------
@testset "TVAL-04: TrotterDomain single-step cross-validation" begin
    deltas = [0.2, 0.1, 0.05]
    errors = Float64[]

    for delta in deltas
        dist = single_step_crossval(TrotterDomain(), delta; with_coherent=true)
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
