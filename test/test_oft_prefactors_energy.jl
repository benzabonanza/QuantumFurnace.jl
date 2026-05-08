using Test
using QuantumFurnace
using LinearAlgebra

# Internal handles (struct + builder are not exported; tests reach via fully
# qualified names per project convention for non-API utilities).
const _EnergyPrefactors    = QuantumFurnace.EnergyDomainPrefactors
const _prepare_oft_pre_eng = QuantumFurnace._prepare_oft_prefactors_energy
const _prefactor_view_qf   = QuantumFurnace._prefactor_view

@testset "qf-e60.1: EnergyDomain Gaussian prefactor cache" begin

    @testset "Builder produces correct shape and types (TEST_HAM)" begin
        ham = TEST_HAM
        sigma = SIGMA
        config = make_config(Lindbladian(), EnergyDomain())
        cfg_pd = QuantumFurnace._precompute_data(config, ham)
        energy_labels = cfg_pd.energy_labels

        p = _prepare_oft_pre_eng(ham.bohr_freqs, energy_labels, sigma)

        @test p isa _EnergyPrefactors
        @test eltype(p.data) === Float64
        @test size(p.data) == (size(ham.bohr_freqs, 1), size(ham.bohr_freqs, 2),
                               length(energy_labels))
        @test p.energy_labels === energy_labels
        @test length(p.energy_to_index) == length(energy_labels)
        for (k, w) in pairs(energy_labels)
            @test p.energy_to_index[w] == k
        end
    end

    @testset "Builder matches the closed-form Gaussian formula at every (i, j, k)" begin
        ham = TEST_HAM
        sigma = SIGMA
        inv_4sigma2 = 1.0 / (4 * sigma^2)
        config = make_config(Lindbladian(), EnergyDomain())
        cfg_pd = QuantumFurnace._precompute_data(config, ham)
        energy_labels = cfg_pd.energy_labels

        p = _prepare_oft_pre_eng(ham.bohr_freqs, energy_labels, sigma)
        d = size(ham.bohr_freqs, 1)
        bf = ham.bohr_freqs

        # Direct formula assertion. `oft!` itself uses the same `exp(-Δ² / 4σ²)`
        # instruction (`src/ofts.jl:17`), so this is the gold-standard reference.
        max_diff = 0.0
        for k in 1:length(energy_labels)
            w = energy_labels[k]
            cache_view = _prefactor_view_qf(p, w)
            @inbounds for j in 1:d, i in 1:d
                expected = exp(-(bf[i, j] - w)^2 * inv_4sigma2)
                diff = abs(cache_view[i, j] - expected)
                max_diff = max(max_diff, diff)
            end
        end
        # FP64 rounding budget — tighter than the algorithmic ε≈1e-5 floor by 11 orders.
        @test max_diff ≤ 1e-15
    end

    @testset "Cache · eigenbasis matches oft! on a non-identity input" begin
        ham = TEST_HAM
        sigma = SIGMA
        inv_4sigma2 = 1.0 / (4 * sigma^2)
        config = make_config(Lindbladian(), EnergyDomain())
        cfg_pd = QuantumFurnace._precompute_data(config, ham)
        energy_labels = cfg_pd.energy_labels

        p = _prepare_oft_pre_eng(ham.bohr_freqs, energy_labels, sigma)

        # Use a real jump to exercise the consumer broadcast pattern.
        jump = first(TEST_JUMPS)
        d = size(ham.bohr_freqs, 1)
        A_ref   = zeros(ComplexF64, d, d)
        A_cache = zeros(ComplexF64, d, d)

        for w in energy_labels
            oft!(A_ref, jump.in_eigenbasis, ham.bohr_freqs, w, inv_4sigma2)
            cache_view = _prefactor_view_qf(p, w)
            @. A_cache = jump.in_eigenbasis * cache_view
            @test maximum(abs, A_ref .- A_cache) ≤ 1e-14
        end
    end

    @testset "Determinism: rebuilds are byte-identical" begin
        ham = TEST_HAM
        sigma = SIGMA
        config = make_config(Lindbladian(), EnergyDomain())
        cfg_pd = QuantumFurnace._precompute_data(config, ham)
        energy_labels = cfg_pd.energy_labels

        p1 = _prepare_oft_pre_eng(ham.bohr_freqs, energy_labels, sigma)
        p2 = _prepare_oft_pre_eng(ham.bohr_freqs, energy_labels, sigma)

        @test p1.data == p2.data
        @test p1.energy_labels == p2.energy_labels
        @test p1.energy_to_index == p2.energy_to_index
    end

    @testset "_prefactor_view round-trips every ω" begin
        ham = N3_HAM
        sigma = SIGMA
        config = make_config(Lindbladian(), EnergyDomain(); num_qubits=3)
        cfg_pd = QuantumFurnace._precompute_data(config, ham)
        energy_labels = cfg_pd.energy_labels

        p = _prepare_oft_pre_eng(ham.bohr_freqs, energy_labels, sigma)
        for (k, w) in pairs(energy_labels)
            view_w = _prefactor_view_qf(p, w)
            @test view_w == @view p.data[:, :, k]
        end
    end

    @testset "_precompute_data populates cache for EnergyDomain configs" begin
        for sim in (Lindbladian(), Thermalize())
            config = make_config(sim, EnergyDomain())
            cfg_pd = QuantumFurnace._precompute_data(config, TEST_HAM)
            @test hasproperty(cfg_pd, :oft_prefactors_energy)
            @test cfg_pd.oft_prefactors_energy isa _EnergyPrefactors

            # Must match the formula bit-for-bit; redundant with the standalone
            # builder test, but here we verify the Config-driven pathway.
            inv_4sigma2 = 1.0 / (4 * config.sigma^2)
            d = size(TEST_HAM.bohr_freqs, 1)
            energy_labels = cfg_pd.energy_labels
            cache = cfg_pd.oft_prefactors_energy
            for k in 1:length(energy_labels), j in 1:d, i in 1:d
                expected = exp(-(TEST_HAM.bohr_freqs[i, j] - energy_labels[k])^2 * inv_4sigma2)
                @test abs(cache.data[i, j, k] - expected) ≤ 1e-15
                # Skip the full d^2 N_w on the inner loops in case the test
                # ever runs at larger d; one ω is enough to catch wiring drift.
                k == 1 || break
            end
        end
    end

    @testset "_precompute_data does NOT populate cache for non-EnergyDomain configs" begin
        for domain in (TimeDomain(), TrotterDomain())
            config = make_config(Lindbladian(), domain)
            ham_or_trott = domain isa TrotterDomain ? TEST_TROTTER : TEST_HAM
            cfg_pd = QuantumFurnace._precompute_data(config, ham_or_trott)
            @test !hasproperty(cfg_pd, :oft_prefactors_energy)
        end
        # BohrDomain Lindbladian _precompute_data takes the alpha branch — no
        # energy_labels there at all, so nothing to cache.
        config = make_config(Lindbladian(), BohrDomain())
        cfg_pd = QuantumFurnace._precompute_data(config, TEST_HAM)
        @test !hasproperty(cfg_pd, :oft_prefactors_energy)
    end

    @testset "Workspace constructors propagate oft_prefactors_energy" begin
        # Lindbladian / EnergyDomain: cache populated, NUFFT slot empty.
        config = make_config(Lindbladian(), EnergyDomain())
        ws = QuantumFurnace.Workspace(config, TEST_HAM, TEST_JUMPS)
        @test hasproperty(ws, :oft_prefactors_energy)
        @test ws.oft_prefactors_energy isa _EnergyPrefactors
        @test ws.oft_nufft_prefactors === nothing

        # Thermalize / EnergyDomain: cache populated.
        config_t = make_config(Thermalize(), EnergyDomain())
        ws_t = QuantumFurnace.Workspace(config_t, TEST_HAM, TEST_JUMPS)
        @test ws_t.oft_prefactors_energy isa _EnergyPrefactors

        # Lindbladian / TimeDomain: NUFFT cache populated, energy cache empty.
        config_time = make_config(Lindbladian(), TimeDomain())
        ws_time = QuantumFurnace.Workspace(config_time, TEST_HAM, TEST_JUMPS)
        @test ws_time.oft_prefactors_energy === nothing
        @test ws_time.oft_nufft_prefactors !== nothing

        # Lindbladian / TrotterDomain: NUFFT cache populated, energy cache empty.
        config_trot = make_config(Lindbladian(), TrotterDomain())
        ws_trot = QuantumFurnace.Workspace(config_trot, TEST_HAM, TEST_TROTTER_JUMPS;
                                           trotter=TEST_TROTTER)
        @test ws_trot.oft_prefactors_energy === nothing
        @test ws_trot.oft_nufft_prefactors !== nothing

        # Lindbladian / BohrDomain: both caches empty (alpha-driven path).
        config_bohr = make_config(Lindbladian(), BohrDomain())
        ws_bohr = QuantumFurnace.Workspace(config_bohr, TEST_HAM, TEST_JUMPS)
        @test ws_bohr.oft_prefactors_energy === nothing
    end

    @testset "Lindbladian dense ctor still works post-field-add" begin
        # `construct_lindbladian` builds an internal scratch Workspace with all
        # `nothing` caches; ensure the new positional slot doesn't break it.
        # qf-e60.{2,3,4} migrate the consumer to use the cache.
        config = make_config(Lindbladian(), EnergyDomain())
        L = construct_lindbladian(TEST_JUMPS, config, TEST_HAM)
        @test L isa Matrix
        @test size(L) == (DIM^2, DIM^2)
    end
end
