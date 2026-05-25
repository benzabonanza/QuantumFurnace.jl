#!/usr/bin/env julia
# qf-0fv: smoke test that compute_true_gap kwarg gates Pass-2 cleanly.
#
# Validates:
#   - default (compute_true_gap=false): Pass-1 only, fewer matvecs
#   - opt-in (compute_true_gap=true):   Pass-1 + Pass-2
#   - On the canonical rho_0=|+⟩⟨+|^⊗N + small-n fixture, Pass-1 and Pass-2
#     gaps must match to <1e-8
#   - Trajectory distances are identical between the two (Pass-2 doesn't
#     affect the trajectory)

using QuantumFurnace
using Test, LinearAlgebra
include(joinpath(dirname(@__DIR__), "test", "test_helpers.jl"))

@testset "qf-0fv: Pass-2 gated behind compute_true_gap" begin
    n = 3
    cfg = make_config(Lindbladian(), EnergyDomain(); num_qubits=n)
    ham = N3_HAM
    jumps = N3_JUMPS
    d = 2^n
    psi = ones(ComplexF64, d) ./ sqrt(2.0^n)
    rho_0 = psi * psi'
    t_grid = collect(range(0.0, 80.0; length=21))

    r_p1 = predict_lindbladian_trajectory(cfg, ham, jumps, rho_0, t_grid; krylovdim=40)
    r_p2 = predict_lindbladian_trajectory(cfg, ham, jumps, rho_0, t_grid;
                                          krylovdim=40, compute_true_gap=true)

    @info "L-path matvec cost" pass1_only=r_p1.total_matvecs pass1_plus_pass2=r_p2.total_matvecs
    @info "L-path gap agreement" pass1=r_p1.spectral_gap pass2=r_p2.spectral_gap
    @test r_p2.total_matvecs > r_p1.total_matvecs
    @test isapprox(r_p1.spectral_gap, r_p2.spectral_gap; rtol=1e-8)
    @test isapprox(r_p1.distances, r_p2.distances; rtol=1e-12)

    # Channel path: same fixture but Thermalize config.
    cfg_ch = make_config(Thermalize(), EnergyDomain(); num_qubits=n)
    k_grid = collect(0:10:200)
    rho_init = Matrix{ComplexF64}(rho_0)

    rc_p1 = predict_channel_trajectory(cfg_ch, ham, jumps, rho_init, k_grid; krylovdim=40)
    rc_p2 = predict_channel_trajectory(cfg_ch, ham, jumps, rho_init, k_grid;
                                       krylovdim=40, compute_true_gap=true)
    @info "Channel-path matvec cost" pass1_only=rc_p1.total_matvecs pass1_plus_pass2=rc_p2.total_matvecs
    @info "Channel-path gap (Lindblad units)" pass1=rc_p1.spectral_gap pass2=rc_p2.spectral_gap
    @test rc_p2.total_matvecs > rc_p1.total_matvecs
    @test isapprox(rc_p1.distances, rc_p2.distances; rtol=1e-12)
    # Pass-1 and Pass-2 should match on this clean small-n fixture
    @test isapprox(rc_p1.spectral_gap, rc_p2.spectral_gap; rtol=1e-3)
end
