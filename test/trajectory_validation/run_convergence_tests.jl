"""
Trajectory Convergence Tests (Phase 5, TVAL-05)
Run via: QUANTUMFURNACE_FULL_TESTS=true julia --project test/trajectory_validation/run_convergence_tests.jl

Verifies that trajectory-averaged density matrix error decreases as 1/sqrt(N_traj),
confirming correct Monte Carlo convergence for both EnergyDomain and TrotterDomain.

NOT included in Pkg.test() -- standalone gated test file.
Gated behind QUANTUMFURNACE_FULL_TESTS=true environment variable.

Requirements validated:
  TVAL-05: 1/sqrt(N_traj) convergence rate for EnergyDomain and TrotterDomain

Methodology:
  Compare trajectory-averaged rho at various N_traj against a high-N_traj reference
  (N_ref=500,000 trajectories with independent seed). Use batch averaging (10 batches
  per N_traj point) to reduce noise in the error estimate. The high-N reference
  approximates the true trajectory channel mean, avoiding systematic bias from
  Lie-Trotter splitting that would occur if comparing against Liouvillian DM evolution.
"""

using Test
using LinearAlgebra
using Random
using Statistics

# Load shared fixtures
include(joinpath(@__DIR__, "..", "test_helpers.jl"))

# ---------------------------------------------------------------------------
# Convergence ratio test helper
# ---------------------------------------------------------------------------
"""
    convergence_ratio_test(domain; with_coherent, delta=0.1, n_batches=10) -> NamedTuple

Run a 1/sqrt(N_traj) convergence test for the given domain.

Computes a high-N_traj reference (N_ref=500,000), then for each test N_traj in
[200, 800, 3200, 12800], runs `n_batches` independent batches and averages the
trace distance to the reference. Batch averaging reduces noise so that ratios
of consecutive mean errors reliably reflect the 1/sqrt(N) scaling.

Expected: ratios close to 2.0 (since N increases by 4x, sqrt(4) = 2).
"""
function convergence_ratio_test(domain; with_coherent::Bool, delta::Float64=0.1, n_batches::Int=10)
    dim = SMALL_DIM  # 8
    psi0 = fill(ComplexF64(1.0), dim) / sqrt(dim)

    # 1. Build trajectory framework
    therm_config = make_small_thermalize_config(domain;
        with_coherent=with_coherent, delta=delta, mixing_time=Float64(delta))
    ham_or_trott = domain isa TrotterDomain ? SMALL_TROTTER : SMALL_HAM
    precomputed = QuantumFurnace._precompute_data(therm_config, ham_or_trott)
    scratch = QuantumFurnace.KrausScratch(ComplexF64, dim)
    fw = build_trajectoryframework(SMALL_JUMPS, ham_or_trott, therm_config,
        precomputed, scratch, delta)

    # 2. Compute high-N reference with independent seed
    #    This approximates the true mean of the trajectory channel,
    #    avoiding Lie-Trotter splitting bias from DM comparison.
    N_ref = 500_000
    ws = QuantumFurnace.TrajectoryWorkspace(fw)
    rng = Random.Xoshiro(999)
    rho_ref = zeros(ComplexF64, dim, dim)
    for _ in 1:N_ref
        psi = copy(psi0)
        step_along_trajectory!(psi, fw, ws, rng)
        rho_ref .+= psi * psi'
    end
    rho_ref ./= N_ref
    rho_ref = Matrix{ComplexF64}((rho_ref + rho_ref') / 2)
    rho_ref ./= tr(rho_ref)

    # 3. For each N_traj, run n_batches independent batches and compute mean error
    n_traj_points = [200, 800, 3_200, 12_800]
    mean_errors = Float64[]

    for ntraj in n_traj_points
        batch_errors = Float64[]
        for b in 1:n_batches
            batch_rng = Random.Xoshiro(1000 * b + ntraj)
            rho_traj = zeros(ComplexF64, dim, dim)

            for _ in 1:ntraj
                psi = copy(psi0)
                step_along_trajectory!(psi, fw, ws, batch_rng)
                rho_traj .+= psi * psi'
            end

            rho_traj ./= ntraj
            rho_traj = Matrix{ComplexF64}((rho_traj + rho_traj') / 2)
            rho_traj ./= tr(rho_traj)

            dist = trace_distance_h(Hermitian(rho_ref), Hermitian(rho_traj))
            push!(batch_errors, dist)
        end
        push!(mean_errors, mean(batch_errors))
    end

    # 4. Compute consecutive ratios
    ratios = [mean_errors[i] / mean_errors[i+1] for i in 1:length(mean_errors)-1]

    return (errors=mean_errors, ratios=ratios, n_traj_points=n_traj_points)
end

# ---------------------------------------------------------------------------
# Gated test execution
# ---------------------------------------------------------------------------
if get(ENV, "QUANTUMFURNACE_FULL_TESTS", "false") == "true"

    @testset "TVAL-05: EnergyDomain 1/sqrt(N) convergence" begin
        result = convergence_ratio_test(EnergyDomain(); with_coherent=false)

        # Print diagnostics
        println("  EnergyDomain convergence (delta=0.1, 10 batches):")
        for (n, e) in zip(result.n_traj_points, result.errors)
            println("    N_traj=$n  mean_trace_dist=$(round(e; sigdigits=4))")
        end
        for (i, r) in enumerate(result.ratios)
            println("    ratio[$i]: $(round(r; sigdigits=3)) (expect ~2.0)")
        end

        # Assert errors decrease monotonically (batch-averaged means should be monotone)
        for i in 1:length(result.errors)-1
            @test result.errors[i] > result.errors[i+1]
        end

        # Assert each ratio is in [1.5, 2.5] (1/sqrt(N) scaling)
        for (i, ratio) in enumerate(result.ratios)
            @test 1.5 <= ratio <= 2.5
        end
    end

    @testset "TVAL-05: TrotterDomain 1/sqrt(N) convergence" begin
        result = convergence_ratio_test(TrotterDomain(); with_coherent=true)

        # Print diagnostics
        println("  TrotterDomain convergence (delta=0.1, with_coherent=true, 10 batches):")
        for (n, e) in zip(result.n_traj_points, result.errors)
            println("    N_traj=$n  mean_trace_dist=$(round(e; sigdigits=4))")
        end
        for (i, r) in enumerate(result.ratios)
            println("    ratio[$i]: $(round(r; sigdigits=3)) (expect ~2.0)")
        end

        # Assert errors decrease monotonically (batch-averaged means should be monotone)
        for i in 1:length(result.errors)-1
            @test result.errors[i] > result.errors[i+1]
        end

        # Assert each ratio is in [1.5, 2.5] (1/sqrt(N) scaling)
        for (i, ratio) in enumerate(result.ratios)
            @test 1.5 <= ratio <= 2.5
        end
    end

else
    @info "Skipping TVAL-05 convergence tests (set QUANTUMFURNACE_FULL_TESTS=true to run)"
end
