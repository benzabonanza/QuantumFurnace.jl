"""
DM reference tests: detailed balance (DMTST-01) and domain error hierarchy (DMTST-02).

DMTST-01 verifies that the BohrDomain Lindbladian with coherent term B has the Gibbs state
as its exact fixed point (KMS detailed balance) on a 3-qubit Heisenberg system.

DMTST-02 verifies that the domain approximation hierarchy produces monotonically increasing
errors: dist(Bohr) <= dist(Energy) <= dist(Time) <= dist(Trotter).
"""

@testset "DMTST-01: Bohr detailed balance (3-qubit)" begin
    # Create Config{Lindbladian} for 3-qubit system with BohrDomain and KMS (coherent term)
    config = make_config(Lindbladian(), BohrDomain(); num_qubits=3, construction=KMS())

    # Construct the 64x64 Liouvillian
    liouv = construct_lindbladian(N3_JUMPS, config, N3_HAM)

    # Full eigendecomposition (small dense matrix -- do not use Arpack eigs)
    eig = eigen(liouv)

    # Find the eigenvalue with smallest |Re(lambda)| -- this is the zero eigenvalue (fixed point)
    ss_idx = argmin(abs.(real.(eig.values)))
    ss_vec = eig.vectors[:, ss_idx]

    # Reshape to 8x8 density matrix, Hermitianize, and normalize
    ss_dm = reshape(ss_vec, N3_DIM, N3_DIM)
    ss_dm = (ss_dm + ss_dm') / 2
    ss_dm ./= tr(ss_dm)

    dist = trace_distance_h(Hermitian(ss_dm), N3_GIBBS)
    @test dist < 1e-10  # KMS detailed balance: Gibbs is exact fixed point, error is machine precision (N3_DIM * eps ~ 8 * 1e-16 ~ 1e-14)
    @info "DMTST-01: Bohr fixed point trace distance to Gibbs" trace_distance=dist threshold=1e-10
end

@testset "DMTST-02: Domain error hierarchy (4-qubit)" begin
    distances = Dict{Symbol, Float64}()

    for (name, domain) in [(:bohr, BohrDomain()), (:energy, EnergyDomain()),
                            (:time, TimeDomain()), (:trotter, TrotterDomain())]
        config = make_config(Lindbladian(),domain)
        trotter_obj = (domain isa TrotterDomain) ? TEST_TROTTER : nothing
        domain_jumps = (domain isa TrotterDomain) ? TEST_TROTTER_JUMPS : TEST_JUMPS
        liouv = construct_lindbladian(domain_jumps, config, TEST_HAM; trotter=trotter_obj)

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
        distances[name] = trace_distance_h(Hermitian(ss_dm), gibbs_ref)
    end

    @info "DMTST-02: Domain distances to Gibbs" bohr=distances[:bohr] energy=distances[:energy] time=distances[:time] trotter=distances[:trotter]

    # Verify monotonic hierarchy with small numerical tolerance
    # 1e-12 tolerance: eigendecomposition rounding may cause tiny ordering violations at machine precision
    @test distances[:bohr] <= distances[:energy] + 1e-12  # Bohr is exact (KMS), Energy approximates filter function
    @info "DMTST-02: Bohr <= Energy" bohr=distances[:bohr] energy=distances[:energy] tolerance=1e-12
    @test distances[:energy] <= distances[:time] + 1e-12  # Energy uses analytic filter, Time uses quadrature
    @info "DMTST-02: Energy <= Time" energy=distances[:energy] time=distances[:time] tolerance=1e-12
    @test distances[:time] <= distances[:trotter] + 1e-12  # Time uses exact eigenbasis, Trotter uses approximate eigenbasis
    @info "DMTST-02: Time <= Trotter" time=distances[:time] trotter=distances[:trotter] tolerance=1e-12
end
