"""
GNS trajectory path validation tests (Phase 14).

GNS-01: Verifies GNS Lindbladian fixed point is a valid density matrix and is distinct
         from the exact Gibbs state (approximation gap).
GNS-02: Verifies trajectory-averaged density matrix converges to the GNS fixed point
         with trace distance < 0.05.

The GNS (approximate detailed balance) code path uses unshifted transition weights
and omits the coherent B term. Its fixed point approximates but does not equal the
exact Gibbs state -- the approximation gap is documented here as a Phase 18 baseline.
"""

using Random

@testset "GNS-01: Lindbladian fixed point (EnergyDomain)" begin
    config = make_small_liouv_config_gns(EnergyDomain())
    liouv = construct_lindbladian(SMALL_JUMPS, config, SMALL_HAM)

    # Full eigendecomposition (64x64 dense matrix)
    eig = eigen(liouv)

    # Extract fixed point: eigenvalue with smallest |Re(lambda)|
    ss_idx = argmin(abs.(real.(eig.values)))
    ss_vec = eig.vectors[:, ss_idx]
    ss_dm = reshape(ss_vec, SMALL_DIM, SMALL_DIM)
    ss_dm = (ss_dm + ss_dm') / 2   # Hermitianize
    ss_dm ./= tr(ss_dm)            # Normalize

    # Validate fixed point is a valid density matrix
    @test isapprox(tr(ss_dm), 1.0, atol=1e-12)
    @test isapprox(ss_dm, ss_dm', atol=1e-12)   # Hermitian
    @test all(eigvals(Hermitian(ss_dm)) .>= -1e-12)  # PSD

    # GNS fixed point is NOT the exact Gibbs state -- measure the approximation gap
    gap = trace_distance_h(Hermitian(ss_dm), SMALL_GIBBS)
    @info "GNS-01: GNS fixed point to Gibbs trace distance (approximation gap)" gap
    @test gap > 1e-6     # Strictly positive (GNS does not reproduce exact Gibbs)
    @test gap < 0.5      # Sanity bound (should not be wildly far from Gibbs)
end

@testset "GNS-01: CPTP completeness (EnergyDomain)" begin
    config = make_small_thermalize_config_gns(EnergyDomain(); delta=0.01)
    precomputed = QuantumFurnace._precompute_data(config, SMALL_HAM)
    scratch = QuantumFurnace.KrausScratch(ComplexF64, SMALL_DIM)
    fw = build_trajectoryframework(
        SMALL_JUMPS, SMALL_HAM, config, precomputed, scratch, config.delta
    )

    @test fw.n_jumps == length(SMALL_JUMPS)

    identity = Matrix{ComplexF64}(I, SMALL_DIM, SMALL_DIM)
    for (a, per_op) in enumerate(fw.per_operator)
        completeness = per_op.K0' * per_op.K0 + fw.delta_eff * per_op.R + per_op.U_residual' * per_op.U_residual
        @test isapprox(completeness, identity; atol=1e-10)
        @test per_op.U_B === nothing  # No coherent term for GNS
    end
end

@testset "GNS-02: Trajectory convergence to GNS fixed point (EnergyDomain)" begin
    # Compute the GNS reference fixed point
    liouv_config = make_small_liouv_config_gns(EnergyDomain())
    liouv = construct_lindbladian(SMALL_JUMPS, liouv_config, SMALL_HAM)
    eig = eigen(liouv)
    ss_idx = argmin(abs.(real.(eig.values)))
    ss_vec = eig.vectors[:, ss_idx]
    gns_fp = reshape(ss_vec, SMALL_DIM, SMALL_DIM)
    gns_fp = (gns_fp + gns_fp') / 2
    gns_fp ./= tr(gns_fp)

    # Run trajectories with ThermalizeConfigGNS
    config = make_small_thermalize_config_gns(EnergyDomain(); delta=0.01, mixing_time=5.0)
    psi0 = zeros(ComplexF64, SMALL_DIM)
    psi0[1] = 1.0  # computational basis |0>

    result = run_trajectories(SMALL_JUMPS, config, psi0, SMALL_HAM; ntraj=1000, seed=42)

    # Convergence to GNS fixed point (NOT Gibbs)
    dist = trace_distance_h(Hermitian(result.rho_mean), Hermitian(gns_fp))
    @info "GNS-02: Trajectory to GNS fixed point trace distance" dist
    @test dist < 0.05

    # DM validity on the final result
    @test isapprox(tr(result.rho_mean), 1.0, atol=1e-10)
    @test isapprox(result.rho_mean, result.rho_mean', atol=1e-10)  # Hermitian
    @test all(eigvals(Hermitian(result.rho_mean)) .>= -1e-10)      # PSD

    # Log the GNS-to-Gibbs gap for Phase 18 baseline
    gibbs_dist = trace_distance_h(Hermitian(gns_fp), SMALL_GIBBS)
    @info "GNS-02: GNS fixed point to Gibbs distance (Phase 18 baseline)" gibbs_dist
end

@testset "GNS-01: BohrDomain detailed balance" begin
    config = make_small_liouv_config_gns(BohrDomain())
    liouv = construct_lindbladian(SMALL_JUMPS, config, SMALL_HAM)

    # Full eigendecomposition
    eig = eigen(liouv)
    ss_idx = argmin(abs.(real.(eig.values)))
    ss_vec = eig.vectors[:, ss_idx]
    ss_dm_bohr = reshape(ss_vec, SMALL_DIM, SMALL_DIM)
    ss_dm_bohr = (ss_dm_bohr + ss_dm_bohr') / 2
    ss_dm_bohr ./= tr(ss_dm_bohr)

    # Validate fixed point is a valid density matrix
    @test isapprox(tr(ss_dm_bohr), 1.0, atol=1e-12)
    @test isapprox(ss_dm_bohr, ss_dm_bohr', atol=1e-12)   # Hermitian
    @test all(eigvals(Hermitian(ss_dm_bohr)) .>= -1e-12)   # PSD

    # GNS fixed point should be distinct from Gibbs
    gap_bohr = trace_distance_h(Hermitian(ss_dm_bohr), SMALL_GIBBS)
    @info "GNS-01 Bohr: GNS fixed point to Gibbs distance" gap_bohr
    @test gap_bohr > 1e-6  # Strictly positive
end
