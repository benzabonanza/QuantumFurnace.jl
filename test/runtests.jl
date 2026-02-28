using Test
using QuantumFurnace
using LinearAlgebra

# Shared fixtures and constants
include("test_helpers.jl")

@testset "QuantumFurnace.jl" begin
    include("test_aqua.jl")
    include("test_compilation.jl")
    include("test_trajectory_fixes.jl")
    include("test_cptp.jl")
    include("test_dm_detailed_balance.jl")
    include("test_dm_scaling.jl")
    include("test_regression.jl")
    include("test_allocation.jl")
    include("test_workspace_independence.jl")
    include("test_threading.jl")
    include("test_gns_trajectory.jl")
    include("test_results.jl")
    include("test_convergence.jl")
    include("test_observable_trajectories.jl")
    include("test_diagnostics.jl")
    include("test_krylov_matvec.jl")
    include("test_krylov_eigsolve.jl")
    include("test_krylov_crossvalidation.jl")

    # Trajectory validation (slow: ~minutes, gated behind env variable)
    if get(ENV, "QUANTUMFURNACE_FULL_TESTS", "false") == "true"
        include("trajectory_validation/run_trajectory_validation.jl")
        include("trajectory_validation/run_convergence_tests.jl")
    else
        @info "Skipping trajectory validation tests (set QUANTUMFURNACE_FULL_TESTS=true to run)"
    end
end
