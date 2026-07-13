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

@testset "GNS-01: Lindbladian fixed point (TrotterDomain)" begin
    config = make_config(Lindbladian(), TrotterDomain(); num_qubits=3, construction=GNS())
    liouv = construct_lindbladian(N3_TROTTER_JUMPS, config, N3_HAM; trotter=N3_TROTTER)

    # Full eigendecomposition (64x64 dense matrix)
    eig = eigen(liouv)

    # Extract fixed point: eigenvalue with smallest |Re(lambda)|
    ss_idx = argmin(abs.(real.(eig.values)))
    ss_vec = eig.vectors[:, ss_idx]
    ss_dm = reshape(ss_vec, N3_DIM, N3_DIM)
    ss_dm = (ss_dm + ss_dm') / 2   # Hermitianize
    ss_dm ./= tr(ss_dm)            # Normalize

    # Validate fixed point is a valid density matrix (structural checks, atol=1e-12 is machine precision)
    @test isapprox(tr(ss_dm), 1.0, atol=1e-12)
    @test isapprox(ss_dm, ss_dm', atol=1e-12)   # Hermitian
    @test all(eigvals(Hermitian(ss_dm)) .>= -1e-12)  # PSD

    # GNS fixed point is NOT the exact Gibbs state -- measure the approximation gap.
    # The Lindbladian is built in the Trotter eigenbasis, so transform the fixed point
    # to the energy eigenbasis before comparing with N3_GIBBS.
    U_t2e = N3_HAM.eigvecs' * N3_TROTTER.eigvecs  # Trotter-to-energy change of basis
    ss_dm_energy = U_t2e * ss_dm * U_t2e'
    ss_dm_energy = (ss_dm_energy + ss_dm_energy') / 2  # re-Hermitianize after basis change
    gap = trace_distance_h(Hermitian(ss_dm_energy), N3_GIBBS)
    @info "GNS-01: GNS fixed point to Gibbs trace distance (approximation gap)" gap

    # GNS approximate detailed balance: gap is strictly positive because GNS omits the
    # coherent B term and uses unshifted weights. Gap magnitude is system-dependent.
    @test gap > 1e-6     # Strictly positive (GNS does not reproduce exact Gibbs)
    @info "GNS-01: Gap lower bound" gap lower_bound=1e-6

    # Sanity bound: gap should be moderate, not wildly wrong. Empirically ~0.08 for this system.
    @test gap < 0.5      # Sanity bound (should be ~0.08, close to EnergyDomain)
    @info "GNS-01: Gap upper bound" gap upper_bound=0.5
end

@testset "GNS-01: CPTP completeness (TrotterDomain)" begin
    config = make_config(Thermalize(), TrotterDomain(); num_qubits=3, construction=GNS(), delta=0.01)
    ws = QuantumFurnace._build_trajectory_workspace(config, N3_HAM, N3_TROTTER_JUMPS;
        trotter=N3_TROTTER, delta=config.delta)

    @test ws.n_jumps == length(N3_JUMPS)

    identity = Matrix{ComplexF64}(I, N3_DIM, N3_DIM)
    # CPTP completeness: K0'K0 + delta*R + U'U = I (algebraic identity)
    # atol=1e-10 allows for FP accumulation across DIM^2 matrix entries
    max_completeness_err = 0.0
    for a in 1:ws.n_jumps
        completeness = ws.K0s[a]' * ws.K0s[a] + ws.delta * ws.Rs[a] + ws.U_residuals[a]' * ws.U_residuals[a]
        err = norm(completeness - identity)
        max_completeness_err = max(max_completeness_err, err)
        @test isapprox(completeness, identity; atol=1e-10)
        @test ws.U_Bs[a] === nothing  # No coherent term for GNS
    end
    @info "GNS-01: CPTP completeness (TrotterDomain, GNS)" n_jumps=ws.n_jumps max_error=max_completeness_err threshold_atol=1e-10
end

@testset "GNS-02: Trajectory convergence to GNS fixed point (TrotterDomain)" begin
    # Compute the GNS reference fixed point
    liouv_config = make_config(Lindbladian(), TrotterDomain(); num_qubits=3, construction=GNS())
    liouv = construct_lindbladian(N3_TROTTER_JUMPS, liouv_config, N3_HAM; trotter=N3_TROTTER)
    eig = eigen(liouv)
    ss_idx = argmin(abs.(real.(eig.values)))
    ss_vec = eig.vectors[:, ss_idx]
    gns_fp = reshape(ss_vec, N3_DIM, N3_DIM)
    gns_fp = (gns_fp + gns_fp') / 2
    gns_fp ./= tr(gns_fp)

    # Run trajectories with Config{Thermalize, <:Any, GNS}
    # mixing_time=100.0 accounts for the corrected (slower) trajectory evolution rate
    # (the per-step CPTP channel uses bare delta, not delta*n_jumps)
    config = make_config(Thermalize(), TrotterDomain(); num_qubits=3, construction=GNS(), delta=0.01, mixing_time=100.0)
    psi0 = zeros(ComplexF64, N3_DIM)
    psi0[1] = 1.0  # computational basis |0>

    n_traj = 1000
    result = run_trajectories(N3_TROTTER_JUMPS, config, psi0, N3_HAM; trotter=N3_TROTTER, ntraj=n_traj, seed=42)

    # Convergence to GNS fixed point (NOT Gibbs)
    # Statistical noise from trajectory averaging: expected error ~ 1/sqrt(N_traj) ~ 0.03
    # Threshold 0.05 gives ~1.6x margin over expected statistical error
    dist = trace_distance_h(Hermitian(result.rho_mean), Hermitian(gns_fp))
    @info "GNS-02: Trajectory to GNS fixed point trace distance" dist
    @test dist < 0.05
    @info "GNS-02: Trajectory convergence" trace_distance=dist threshold=0.05 n_trajectories=n_traj expected_noise="~1/sqrt($n_traj)≈0.03"

    # DM validity on the final result (structural checks, machine precision tolerance)
    @test isapprox(tr(result.rho_mean), 1.0, atol=1e-10)
    @test isapprox(result.rho_mean, result.rho_mean', atol=1e-10)  # Hermitian
    @test all(eigvals(Hermitian(result.rho_mean)) .>= -1e-10)      # PSD

    # Log the GNS-to-Gibbs gap for Phase 18 baseline.
    # Transform gns_fp from Trotter eigenbasis to energy eigenbasis for comparison.
    U_t2e = N3_HAM.eigvecs' * N3_TROTTER.eigvecs
    gns_fp_energy = U_t2e * gns_fp * U_t2e'
    gns_fp_energy = (gns_fp_energy + gns_fp_energy') / 2
    gibbs_dist = trace_distance_h(Hermitian(gns_fp_energy), N3_GIBBS)
    @info "GNS-02: GNS fixed point to Gibbs distance (Phase 18 baseline)" gibbs_dist
end

@testset "GNS-01: BohrDomain detailed balance" begin
    config = make_config(Lindbladian(), BohrDomain(); num_qubits=3, construction=GNS())
    liouv = construct_lindbladian(N3_JUMPS, config, N3_HAM)

    # Full eigendecomposition
    eig = eigen(liouv)
    ss_idx = argmin(abs.(real.(eig.values)))
    ss_vec = eig.vectors[:, ss_idx]
    ss_dm_bohr = reshape(ss_vec, N3_DIM, N3_DIM)
    ss_dm_bohr = (ss_dm_bohr + ss_dm_bohr') / 2
    ss_dm_bohr ./= tr(ss_dm_bohr)

    # Validate fixed point is a valid density matrix (structural checks, machine precision)
    @test isapprox(tr(ss_dm_bohr), 1.0, atol=1e-12)
    @test isapprox(ss_dm_bohr, ss_dm_bohr', atol=1e-12)   # Hermitian
    @test all(eigvals(Hermitian(ss_dm_bohr)) .>= -1e-12)   # PSD

    # GNS fixed point should be distinct from Gibbs
    # GNS approximate detailed balance: strictly positive gap, system-dependent magnitude
    gap_bohr = trace_distance_h(Hermitian(ss_dm_bohr), N3_GIBBS)
    @info "GNS-01 Bohr: GNS fixed point to Gibbs distance" gap_bohr
    @test gap_bohr > 1e-6  # Strictly positive
    @info "GNS-01 Bohr: Gap lower bound" gap_bohr lower_bound=1e-6
end
