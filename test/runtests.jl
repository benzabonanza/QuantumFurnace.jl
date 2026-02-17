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
    include("test_fitting.jl")
    include("test_observable_trajectories.jl")
end
