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
end
