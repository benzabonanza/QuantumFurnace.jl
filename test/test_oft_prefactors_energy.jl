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
end
