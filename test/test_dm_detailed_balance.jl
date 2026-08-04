"""
DM reference tests: detailed balance (DMTST-01) and domain error hierarchy (DMTST-02).

DMTST verifies the Gibbs fixed point independently in every construction
domain. No ordering between independent quadrature and Trotter errors is
assumed; those errors can cancel.
"""

@testset "Dense Gibbs fixed point by domain (4-qubit)" begin
    distances = Dict{Symbol, Float64}()
    thresholds = Dict(:bohr => 1e-12, :energy => 1e-12, :time => 1e-12,
        :trotter => 1e-7)

    for (name, domain) in [(:bohr, BohrDomain()), (:energy, EnergyDomain()),
                            (:time, TimeDomain()), (:trotter, TrotterDomain())]
        config = make_config(Lindbladian(),domain)
        trotter_obj = (domain isa TrotterDomain) ? TEST_TROTTER : nothing
        domain_jumps = (domain isa TrotterDomain) ? TEST_TROTTER_JUMPS : TEST_JUMPS
        liouv = construct_lindbladian(domain_jumps, config, TEST_HAM; trotter=trotter_obj)

        # Full eigendecomposition (256x256 dense matrix -- fast enough)
        eig = eigen(liouv)
        ss_idx = argmin(abs.(eig.values))
        ss_vec = eig.vectors[:, ss_idx]
        ss_dm = reshape(ss_vec, DIM, DIM)
        ss_dm = (ss_dm + ss_dm') / 2
        ss_dm ./= tr(ss_dm)
        @test abs(eig.values[ss_idx]) < 1e-10
        @test norm(liouv * vec(ss_dm)) < 1e-10

        # TrotterDomain Liouvillian operates in Trotter eigenbasis, so transform
        # Gibbs state: eigenbasis -> computational -> Trotter eigenbasis
        gibbs_ref = if domain isa TrotterDomain
            gibbs_comp = TEST_HAM.eigvecs * TEST_GIBBS * TEST_HAM.eigvecs'
            Hermitian(TEST_TROTTER.eigvecs' * gibbs_comp * TEST_TROTTER.eigvecs)
        else
            TEST_GIBBS
        end
        distances[name] = trace_distance_h(Hermitian(ss_dm), gibbs_ref)
        @test distances[name] < thresholds[name]
    end

    @info "Dense fixed-point distances" distances thresholds
end
