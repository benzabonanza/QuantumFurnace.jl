"""
DM reference tests: detailed balance (DMTST-01) and domain error hierarchy (DMTST-02).

DMTST-01 verifies that the BohrDomain Lindbladian with coherent term B has the Gibbs state
as its exact fixed point (KMS detailed balance) on a 3-qubit Heisenberg system.

DMTST-02 verifies that the domain approximation hierarchy produces monotonically increasing
errors: dist(Bohr) <= dist(Energy) <= dist(Time) <= dist(Trotter).
"""

@testset "DMTST-01: Bohr detailed balance (3-qubit)" begin
    # Create LiouvConfig for 3-qubit system with BohrDomain and coherent term
    config = LiouvConfig(
        num_qubits = 3,
        with_coherent = true,
        with_linear_combination = true,
        domain = BohrDomain(),
        beta = BETA,
        sigma = SIGMA,
        a = BETA / 30.0,
        b = 0.4,
        num_energy_bits = NUM_ENERGY_BITS,
        w0 = W0,
        t0 = T0,
        num_trotter_steps_per_t0 = NUM_TROTTER_STEPS_PER_T0,
    )

    # Construct the 64x64 Liouvillian
    liouv = construct_lindbladian(SMALL_JUMPS, config, SMALL_HAM)

    # Full eigendecomposition (small dense matrix -- do not use Arpack eigs)
    eig = eigen(liouv)

    # Find the eigenvalue with smallest |Re(lambda)| -- this is the zero eigenvalue (fixed point)
    ss_idx = argmin(abs.(real.(eig.values)))
    ss_vec = eig.vectors[:, ss_idx]

    # Reshape to 8x8 density matrix, Hermitianize, and normalize
    ss_dm = reshape(ss_vec, SMALL_DIM, SMALL_DIM)
    ss_dm = (ss_dm + ss_dm') / 2
    ss_dm ./= tr(ss_dm)

    dist = QuantumFurnace.trace_distance_h(Hermitian(ss_dm), SMALL_GIBBS)
    @info "DMTST-01: Bohr fixed point trace distance to Gibbs" dist
    @test dist < 1e-10
end

@testset "DMTST-02: Domain error hierarchy (4-qubit)" begin
    distances = Dict{Symbol, Float64}()

    for (name, domain) in [(:bohr, BohrDomain()), (:energy, EnergyDomain()),
                            (:time, TimeDomain()), (:trotter, TrotterDomain())]
        config = make_liouv_config(domain)
        trotter_obj = (domain isa TrotterDomain) ? TEST_TROTTER : nothing
        liouv = construct_lindbladian(TEST_JUMPS, config, TEST_HAM; trotter=trotter_obj)

        # Full eigendecomposition (256x256 dense matrix -- fast enough)
        eig = eigen(liouv)
        ss_idx = argmin(abs.(real.(eig.values)))
        ss_vec = eig.vectors[:, ss_idx]
        ss_dm = reshape(ss_vec, DIM, DIM)
        ss_dm = (ss_dm + ss_dm') / 2
        ss_dm ./= tr(ss_dm)

        # TrotterDomain Liouvillian operates in Trotter eigenbasis, so transform
        # Gibbs state: eigenbasis -> computational -> Trotter eigenbasis
        gibbs_ref = if domain isa TrotterDomain
            gibbs_comp = TEST_HAM.eigvecs * TEST_GIBBS * TEST_HAM.eigvecs'
            Hermitian(TEST_TROTTER.eigvecs' * gibbs_comp * TEST_TROTTER.eigvecs)
        else
            TEST_GIBBS
        end
        distances[name] = QuantumFurnace.trace_distance_h(Hermitian(ss_dm), gibbs_ref)
    end

    @info "DMTST-02: Domain distances to Gibbs" bohr=distances[:bohr] energy=distances[:energy] time=distances[:time] trotter=distances[:trotter]

    # Verify monotonic hierarchy with small numerical tolerance
    @test distances[:bohr] <= distances[:energy] + 1e-12
    @test distances[:energy] <= distances[:time] + 1e-12
    @test distances[:time] <= distances[:trotter] + 1e-12
end
