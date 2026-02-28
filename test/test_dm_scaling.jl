using Test
using LinearAlgebra

# DMTST-03: Single DM Euler step error scales as O(delta^2)
# DMTST-04: Multi-step DM Euler error accumulates as O(delta)
#
# Uses the Liouvillian matrix L directly for deterministic DM evolution.
# Euler step: rho(t+delta) = rho(t) + delta * reshape(L * vec(rho(t)), DIM, DIM)
# Reference: exact matrix exponential exp(delta * L).

@testset "DMTST-03: Single-step Euler error O(delta^2)" begin
    # Build Liouvillian for EnergyDomain (simplest non-trivial domain)
    config = make_config(Lindbladian(),EnergyDomain())
    L = construct_lindbladian(TEST_JUMPS, config, TEST_HAM)

    # Initial state: maximally mixed
    rho0 = Matrix{ComplexF64}(I, DIM, DIM) / DIM
    rho0_vec = vec(rho0)

    # Delta sweep
    deltas = [0.1, 0.05, 0.025, 0.0125]
    errors = Float64[]

    for delta in deltas
        # Euler step
        rho_euler = rho0 + delta * reshape(L * rho0_vec, DIM, DIM)
        # Exact step via matrix exponential
        rho_exact = reshape(exp(delta * L) * rho0_vec, DIM, DIM)
        # Error (Frobenius norm of vectorized difference)
        err = norm(vec(rho_euler) - vec(rho_exact))
        push!(errors, err)
    end

    # Compute ratios of consecutive errors
    ratios = [errors[i] / errors[i+1] for i in 1:length(errors)-1]

    # For O(delta^2), when delta halves the error should decrease by factor ~4
    # Bounds [3.0, 5.0] allow for sub-leading O(delta^3) terms at finite delta
    for (i, ratio) in enumerate(ratios)
        @test 3.0 <= ratio <= 5.0
        @info "DMTST-03: Euler single-step ratio" i ratio expected=4.0 lower=3.0 upper=5.0
    end
    @info "DMTST-03: Single-step errors" deltas errors
end

@testset "DMTST-04: Multi-step accumulated error O(delta)" begin
    # Build Liouvillian for EnergyDomain
    config = make_config(Lindbladian(),EnergyDomain())
    L = construct_lindbladian(TEST_JUMPS, config, TEST_HAM)

    # Initial state: maximally mixed
    rho0 = Matrix{ComplexF64}(I, DIM, DIM) / DIM
    rho0_vec = vec(rho0)

    # Fixed total evolution time
    T = 0.5

    # Exact evolution (compute once)
    rho_exact_T = reshape(exp(T * L) * rho0_vec, DIM, DIM)

    # Delta sweep
    deltas = [0.1, 0.05, 0.025, 0.0125]
    errors = Float64[]

    for delta in deltas
        num_steps = Int(round(T / delta))
        rho = copy(rho0)
        for _ in 1:num_steps
            rho .+= delta .* reshape(L * vec(rho), DIM, DIM)
        end
        err = norm(vec(rho) - vec(rho_exact_T))
        push!(errors, err)
    end

    # Compute ratios
    ratios = [errors[i] / errors[i+1] for i in 1:length(errors)-1]

    # For O(delta), when delta halves the error should decrease by factor ~2
    # Bounds [1.5, 2.5] allow for sub-leading terms at finite delta
    for (i, ratio) in enumerate(ratios)
        @test 1.5 <= ratio <= 2.5
        @info "DMTST-04: Multi-step accumulated ratio" i ratio expected=2.0 lower=1.5 upper=2.5
    end
    @info "DMTST-04: Multi-step errors" T deltas errors
end

# DMTST-05: Coherent term B consistency across domains
# B_bohr (exact, Bohr/Energy domain) vs B_time (time quadrature) vs B_trotter (Trotter + time quadrature)
# Expected: B_bohr ~ B_time within TOL_QUADRATURE; B_trotter has additional Trotter error.

@testset "DMTST-05: Coherent term B consistency" begin
    jump = TEST_JUMPS[1]  # X on site 1

    # B_bohr (exact, in Hamiltonian eigenbasis)
    config_bohr = make_config(Lindbladian(),BohrDomain())
    precomputed_bohr = QuantumFurnace._precompute_data(config_bohr, TEST_HAM)
    B_bohr_val = B_bohr(TEST_HAM, JumpOp[jump], config_bohr)
    rmul!(B_bohr_val, precomputed_bohr.gamma_norm_factor)

    # B_time (time quadrature, in Hamiltonian eigenbasis)
    config_time = make_config(Lindbladian(),TimeDomain())
    precomputed_time = QuantumFurnace._precompute_data(config_time, TEST_HAM)
    B_time_val = B_time(JumpOp[jump], TEST_HAM, precomputed_time.b_minus, precomputed_time.b_plus,
                        T0, BETA, SIGMA)
    rmul!(B_time_val, precomputed_time.gamma_norm_factor)

    # B_trotter (Trotter + time quadrature, in Trotter eigenbasis)
    # Use pre-built trotter-basis jump (TEST_TROTTER_JUMPS[1] corresponds to TEST_JUMPS[1])
    config_trott = make_config(Lindbladian(),TrotterDomain())
    precomputed_trott = QuantumFurnace._precompute_data(config_trott, TEST_TROTTER)
    trotter_jump = TEST_TROTTER_JUMPS[1]
    B_trott = B_trotter(JumpOp[trotter_jump], TEST_TROTTER, precomputed_trott.b_minus, precomputed_trott.b_plus,
                         BETA, SIGMA)
    rmul!(B_trott, precomputed_trott.gamma_norm_factor)
    # Transform from Trotter eigenbasis to Hamiltonian eigenbasis
    U_t2e = TEST_TROTTER.eigvecs' * TEST_HAM.eigvecs  # Trotter eigenbasis <- H eigenbasis
    B_trott_in_eigen = U_t2e' * B_trott * U_t2e

    # Compute distances
    dist_bohr_time = norm(B_bohr_val - B_time_val)
    dist_bohr_trott = norm(B_bohr_val - B_trott_in_eigen)
    dist_time_trott = norm(B_time_val - B_trott_in_eigen)

    # B_bohr and B_time agree up to time quadrature tolerance
    # Time quadrature error: O(1/N_time_points) with N=4096, well below 1e-6
    @test dist_bohr_time < TOL_QUADRATURE
    @info "DMTST-05: B_bohr vs B_time" distance=dist_bohr_time threshold=TOL_QUADRATURE

    # Trotter error is at least as large as time quadrature error (with small numerical margin)
    # This is a monotonicity check: Trotter approximation adds error on top of quadrature
    @test dist_bohr_trott >= dist_bohr_time - 1e-10
    @info "DMTST-05: Trotter error >= quadrature error" dist_bohr_trott dist_bohr_time margin=1e-10

    # Trotter error on B term (tightened after basis fix in B_trotter)
    # Measured ~1e-8; threshold 1e-5 gives 1000x margin for system-size variation
    @test dist_bohr_trott < 1e-5
    @info "DMTST-05: B Trotter total error" distance=dist_bohr_trott threshold=1e-5

    @info "DMTST-05: Cross-domain distances" dist_bohr_time dist_bohr_trott dist_time_trott
end

# DMTST-06: OFT consistency across domains
# oft! (analytical, Energy domain) vs time_oft! (time quadrature) vs trotter_oft! (Trotter)
# Expected: energy ~ time within TOL_QUADRATURE; trotter has additional error.

@testset "DMTST-06: OFT consistency" begin
    jump = TEST_JUMPS[1]  # X on site 1
    w = -3 * W0           # Test energy value

    # Energy OFT (analytical, Bohr/Energy domain)
    energy_oft_prefactor = 1 / sqrt(SIGMA * sqrt(2 * pi))
    A_energy = Matrix{ComplexF64}(undef, DIM, DIM)
    oft!(A_energy, jump.in_eigenbasis, TEST_HAM.bohr_freqs, w, 1.0 / (4 * SIGMA^2))
    A_energy .*= energy_oft_prefactor

    # Time OFT (time domain quadrature)
    # Reconstruct oft_time_labels since they are not stored in precomputed_data
    energy_labels = QuantumFurnace._create_energy_labels(NUM_ENERGY_BITS, W0)
    time_labels_full = energy_labels .* (T0 / W0)
    oft_time_labels = QuantumFurnace._truncate_time_labels_for_oft(time_labels_full, SIGMA)

    time_oft_prefactor = T0 * sqrt(SIGMA * sqrt(2 / pi) / (2 * pi))
    caches = QuantumFurnace.OFTCaches{Float64}(DIM)
    A_time = Matrix{ComplexF64}(undef, DIM, DIM)
    QuantumFurnace.time_oft!(A_time, caches, jump, w, TEST_HAM, oft_time_labels, SIGMA)
    A_time .*= time_oft_prefactor

    # Trotter OFT: trotter_oft! needs jump in Trotter eigenbasis
    jump_trott = TEST_TROTTER_JUMPS[1]
    A_trott = Matrix{ComplexF64}(undef, DIM, DIM)
    QuantumFurnace.trotter_oft!(A_trott, caches, jump_trott, w, TEST_TROTTER, oft_time_labels, SIGMA)
    A_trott .*= time_oft_prefactor
    # Transform result from Trotter eigenbasis back to H-eigenbasis
    U_t2e = TEST_TROTTER.eigvecs' * TEST_HAM.eigvecs
    A_trott_in_eigen = U_t2e' * A_trott * U_t2e

    # Compute distances
    dist_energy_time = norm(A_energy - A_time)
    dist_energy_trott = norm(A_energy - A_trott_in_eigen)
    dist_time_trott = norm(A_time - A_trott_in_eigen)

    # Energy and time OFT agree up to quadrature tolerance
    # Time quadrature error: O(1/N_time_points) with N=4096, well below 1e-6
    @test dist_energy_time < TOL_QUADRATURE
    @info "DMTST-06: A_energy vs A_time" distance=dist_energy_time threshold=TOL_QUADRATURE

    # Trotter error is at least as large as time quadrature error (with small numerical margin)
    # Monotonicity: Trotter approximation adds error on top of quadrature
    @test dist_energy_trott >= dist_energy_time - 1e-10
    @info "DMTST-06: Trotter error >= quadrature error" dist_energy_trott dist_energy_time margin=1e-10

    # Trotter OFT error (measured ~1.5e-8, tight threshold)
    # Threshold 1e-5 gives ~1000x margin for system-size variation
    @test dist_energy_trott < 1e-5
    @info "DMTST-06: OFT Trotter total error" distance=dist_energy_trott threshold=1e-5

    @info "DMTST-06: Cross-domain distances" dist_energy_time dist_energy_trott dist_time_trott
end

# DMTST-06b: NUFFT OFT consistency
# Verifies the NUFFT-accelerated OFT matches both the direct summation methods (time_oft!, trotter_oft!)
# and the analytical energy-domain OFT.

@testset "DMTST-06b: NUFFT OFT consistency" begin
    jump = TEST_JUMPS[1]  # X on site 1
    w = -3 * W0           # Test energy value (must be on energy grid)

    # --- Analytical OFT (reference) ---
    energy_oft_prefactor = 1 / sqrt(SIGMA * sqrt(2 * pi))
    A_energy = Matrix{ComplexF64}(undef, DIM, DIM)
    oft!(A_energy, jump.in_eigenbasis, TEST_HAM.bohr_freqs, w, 1.0 / (4 * SIGMA^2))
    A_energy .*= energy_oft_prefactor

    # --- Shared setup ---
    energy_labels = QuantumFurnace._create_energy_labels(NUM_ENERGY_BITS, W0)
    time_labels_full = energy_labels .* (T0 / W0)
    oft_time_labels = QuantumFurnace._truncate_time_labels_for_oft(time_labels_full, SIGMA)
    time_oft_prefactor = T0 * sqrt(SIGMA * sqrt(2 / pi) / (2 * pi))
    caches = QuantumFurnace.OFTCaches{Float64}(DIM)

    # === Time NUFFT OFT ===
    config_time = make_config(Lindbladian(),TimeDomain())
    precomputed_time = QuantumFurnace._precompute_data(config_time, TEST_HAM)

    # Sanity: test energy must be on the NUFFT grid (structural check, no @info)
    @test haskey(precomputed_time.oft_nufft_prefactors.energy_to_index, w)

    nufft_pf_time = QuantumFurnace._prefactor_view(precomputed_time.oft_nufft_prefactors, w)
    A_nufft_time = jump.in_eigenbasis .* nufft_pf_time
    A_nufft_time .*= time_oft_prefactor

    # Direct time_oft! for comparison
    A_time = Matrix{ComplexF64}(undef, DIM, DIM)
    QuantumFurnace.time_oft!(A_time, caches, jump, w, TEST_HAM, oft_time_labels, SIGMA)
    A_time .*= time_oft_prefactor

    # === Trotter NUFFT OFT ===
    config_trott = make_config(Lindbladian(),TrotterDomain())
    precomputed_trott = QuantumFurnace._precompute_data(config_trott, TEST_TROTTER)

    # Sanity: test energy must be on the NUFFT grid (structural check, no @info)
    @test haskey(precomputed_trott.oft_nufft_prefactors.energy_to_index, w)

    jump_trott = TEST_TROTTER_JUMPS[1]
    U_t2e = TEST_TROTTER.eigvecs' * TEST_HAM.eigvecs

    nufft_pf_trott = QuantumFurnace._prefactor_view(precomputed_trott.oft_nufft_prefactors, w)
    A_nufft_trott = jump_trott.in_eigenbasis .* nufft_pf_trott
    A_nufft_trott .*= time_oft_prefactor
    A_nufft_trott_in_eigen = U_t2e' * A_nufft_trott * U_t2e

    # Direct trotter_oft! for comparison
    A_trott = Matrix{ComplexF64}(undef, DIM, DIM)
    QuantumFurnace.trotter_oft!(A_trott, caches, jump_trott, w, TEST_TROTTER, oft_time_labels, SIGMA)
    A_trott .*= time_oft_prefactor
    A_trott_in_eigen = U_t2e' * A_trott * U_t2e

    # Compute distances
    dist_nufft_time_vs_time = norm(A_nufft_time - A_time)
    dist_nufft_trott_vs_trott = norm(A_nufft_trott_in_eigen - A_trott_in_eigen)
    dist_nufft_time_vs_energy = norm(A_nufft_time - A_energy)
    dist_nufft_trott_vs_energy = norm(A_nufft_trott_in_eigen - A_energy)

    # NUFFT time OFT matches direct time_oft! (both compute same sum, NUFFT uses FFT with eps=1e-12)
    # Threshold 1e-10 allows for floating-point reordering differences between FFT and direct sum
    @test dist_nufft_time_vs_time < 1e-10
    @info "DMTST-06b: NUFFT time vs direct time" distance=dist_nufft_time_vs_time threshold=1e-10

    # NUFFT trotter OFT matches direct trotter_oft!
    # Same reasoning: FFT vs direct sum reordering, threshold 1e-10
    @test dist_nufft_trott_vs_trott < 1e-10
    @info "DMTST-06b: NUFFT trotter vs direct trotter" distance=dist_nufft_trott_vs_trott threshold=1e-10

    # NUFFT time OFT matches analytical (same tolerance as DMTST-06 time vs energy)
    # Time quadrature error: O(1/N_time_points), well below TOL_QUADRATURE
    @test dist_nufft_time_vs_energy < TOL_QUADRATURE
    @info "DMTST-06b: NUFFT time vs analytical" distance=dist_nufft_time_vs_energy threshold=TOL_QUADRATURE

    # NUFFT Trotter OFT error (measured ~1.5e-8, tight threshold)
    # Threshold 1e-5 gives ~1000x margin for system-size variation
    @test dist_nufft_trott_vs_energy < 1e-5
    @info "DMTST-06b: NUFFT trotter vs analytical" distance=dist_nufft_trott_vs_energy threshold=1e-5
end
