using Test
using QuantumFurnace
using LinearAlgebra

# Shared fixtures and constants
include("test_helpers.jl")

@testset "QuantumFurnace.jl" begin
    include("test_compilation.jl")
    include("test_trajectory_fixes.jl")
    # Future phases add test includes here
end
