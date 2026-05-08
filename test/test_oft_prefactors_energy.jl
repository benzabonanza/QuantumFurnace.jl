using Test
using QuantumFurnace
using LinearAlgebra
using Random

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

# ---------------------------------------------------------------------------
# qf-e60.5: cross-checks against `oft!` (gold-standard reference)
#
# After qf-e60.{2,3,4} all moved EnergyDomain consumers off `oft!` onto the
# `EnergyDomainPrefactors` cache, the test harness uses `oft!` itself as the
# correctness oracle. The cache and `oft!` evaluate the **same** closed-form
# Gaussian (`exp(-(ω - bohr_freqs)^2 / 4σ²)`) at FP64; absolute equality is
# the right tolerance — any drift is a bug.
# ---------------------------------------------------------------------------
@testset "qf-e60.5: oft! gold-standard cross-checks" begin

    @testset "Per-(jump, ω) bit-equivalence at TEST_HAM (n=4)" begin
        config = make_config(Lindbladian(), EnergyDomain())
        ws = QuantumFurnace.Workspace(config, TEST_HAM, TEST_JUMPS)
        inv_4sigma2 = 1.0 / (4 * config.sigma^2)
        bf = TEST_HAM.bohr_freqs
        d = size(bf, 1)

        A_ref   = zeros(ComplexF64, d, d)
        A_cache = zeros(ComplexF64, d, d)

        max_abs = 0.0
        for jump in TEST_JUMPS, w in ws.energy_labels
            oft!(A_ref, jump.in_eigenbasis, bf, w, inv_4sigma2)
            cache_view = _prefactor_view_qf(ws.oft_prefactors_energy, w)
            @. A_cache = jump.in_eigenbasis * cache_view
            err = maximum(abs, A_ref .- A_cache)
            max_abs = max(max_abs, err)
        end
        # FP64 ULP scale — both paths execute the same `exp(-Δ² · inv_4sigma2)`
        # broadcast.
        @test max_abs ≤ 1e-15
    end

    @testset "Per-(jump, ω) bit-equivalence at N3_HAM (n=3)" begin
        config = make_config(Lindbladian(), EnergyDomain(); num_qubits=3)
        ws = QuantumFurnace.Workspace(config, N3_HAM, N3_JUMPS)
        inv_4sigma2 = 1.0 / (4 * config.sigma^2)
        bf = N3_HAM.bohr_freqs
        d = size(bf, 1)

        A_ref   = zeros(ComplexF64, d, d)
        A_cache = zeros(ComplexF64, d, d)

        max_abs = 0.0
        for jump in N3_JUMPS, w in ws.energy_labels
            oft!(A_ref, jump.in_eigenbasis, bf, w, inv_4sigma2)
            cache_view = _prefactor_view_qf(ws.oft_prefactors_energy, w)
            @. A_cache = jump.in_eigenbasis * cache_view
            err = maximum(abs, A_ref .- A_cache)
            max_abs = max(max_abs, err)
        end
        @test max_abs ≤ 1e-15
    end

    @testset "Lindbladian construct vs apply_lindbladian! parity (n=4)" begin
        # Confirm: dense construct_lindbladian (which uses _jump_contribution!
        # → cache) and matrix-free apply_lindbladian! (which uses cache via
        # _prefactor_view) agree to FP64 GEMM precision.
        config = make_config(Lindbladian(), EnergyDomain())
        L = construct_lindbladian(TEST_JUMPS, config, TEST_HAM)

        ws = QuantumFurnace.Workspace(config, TEST_HAM, TEST_JUMPS)
        d = DIM
        # Probe a few random Hermitian rho's; matvec must equal L * vec(rho).
        max_rel = 0.0
        for seed in 1:3
            rng = Random.MersenneTwister(seed)
            M = randn(rng, ComplexF64, d, d)
            rho = (M + M') / 2
            rho_vec = vec(rho)
            out_dense = reshape(L * rho_vec, d, d)
            apply_lindbladian!(ws, rho, config, TEST_HAM)
            out_mv = ws.scratch.rho_out
            rel = norm(out_mv - out_dense) / max(norm(out_dense), 1e-30)
            max_rel = max(max_rel, rel)
        end
        @test max_rel ≤ 1e-12
    end

    @testset "Krylov spectral gap parity (n=3)" begin
        # `krylov_spectral_gap` matvecs through `apply_lindbladian!` which uses
        # the cache. Cross-check against the dense superop's eigendecomposition.
        config = make_config(Lindbladian(), EnergyDomain(); num_qubits=3)
        L = construct_lindbladian(N3_JUMPS, config, N3_HAM)
        evs = eigvals(L)
        # smallest non-zero |Re(λ)|: gap candidate
        re_evs = sort(abs.(real.(evs)))
        gap_dense = re_evs[2]  # idx 1 is the zero eigenvalue (fixed point)

        result = krylov_spectral_gap(config, N3_HAM, N3_JUMPS; krylovdim=20)
        @test abs(result.spectral_gap - gap_dense) ≤ 1e-9
        # fixed point should be (close to) trace-1 Hermitian PSD
        @test abs(tr(result.fixed_point) - 1.0) ≤ 1e-10
        @test norm(result.fixed_point - result.fixed_point') ≤ 1e-10
    end

    @testset "Channel sandwich (Thermalize/EnergyDomain) parity at n=4" begin
        # `_accumulate_jump_sandwich!` (qf-e60.4 site) vs hand-rolled reference
        # built directly from oft! on the same input. Tighter than 1e-12 because
        # both paths go through identical mul!() chains; only the cache vs
        # oft! source differs.
        config = make_config(Thermalize(), EnergyDomain())
        ws = QuantumFurnace.Workspace(config, TEST_HAM, TEST_JUMPS)
        d = DIM
        rho = let
            rng = Random.MersenneTwister(7)
            M = randn(rng, ComplexF64, d, d)
            ρ = M * M'
            ρ ./= tr(ρ)
            ρ
        end

        out_cache = zeros(ComplexF64, d, d)
        QuantumFurnace._accumulate_jump_sandwich!(
            out_cache, ws, rho, 1.0, config, TEST_HAM)

        # Hand-rolled reference using oft! directly.
        out_ref = zeros(ComplexF64, d, d)
        inv_4sigma2 = 1.0 / (4 * config.sigma^2)
        prefactor = ws.oft_domain_prefactor * ws.gamma_norm_factor
        L = zeros(ComplexF64, d, d)
        tmp = zeros(ComplexF64, d, d)
        for (k, jump) in enumerate(TEST_JUMPS)
            is_herm = ws.jump_hermitian[k]
            for w_raw in ws.energy_labels
                if is_herm
                    w_raw > 1e-12 && continue
                    w = abs(w_raw)
                    oft!(L, jump.in_eigenbasis, TEST_HAM.bohr_freqs, w, inv_4sigma2)
                    rate2 = prefactor * ws.transition(w)
                    mul!(tmp, rho, L')
                    mul!(out_ref, L, tmp, rate2, 1.0)
                    if w > 1e-12
                        rate2_neg = prefactor * ws.transition(-w)
                        mul!(tmp, rho, L)
                        mul!(out_ref, L', tmp, rate2_neg, 1.0)
                    end
                else
                    w = w_raw
                    oft!(L, jump.in_eigenbasis, TEST_HAM.bohr_freqs, w, inv_4sigma2)
                    rate2 = prefactor * ws.transition(w)
                    mul!(tmp, rho, L')
                    mul!(out_ref, L, tmp, rate2, 1.0)
                end
            end
        end

        rel = norm(out_cache - out_ref) / norm(out_ref)
        @test rel ≤ 1e-13
    end

    @testset "Memory cost: cache size matches N_w · d² · 8 bytes (n=3)" begin
        config = make_config(Lindbladian(), EnergyDomain(); num_qubits=3)
        ws = QuantumFurnace.Workspace(config, N3_HAM, N3_JUMPS)
        cache = ws.oft_prefactors_energy
        # Float64 → 8 bytes per element.
        expected_bytes = sizeof(cache.data)
        d = N3_DIM  # 8
        N_w = length(ws.energy_labels)
        @test expected_bytes == d * d * N_w * 8
        # Float64-specific: complex storage would be 16 B/elt.
        @test eltype(cache.data) === Float64
    end

    @testset "DLL has no EnergyDomain dispatches (audit)" begin
        # Acceptance #6: DLL is BohrDomain + TimeDomain only. There is no
        # `_precompute_data(::Config{Lindbladian, EnergyDomain, DLL}, …)`
        # dispatch — Workspace construction must fail. The cache is therefore
        # not relevant for DLL paths.
        config = make_config(Lindbladian(), EnergyDomain(); construction = DLL())
        @test_throws Exception QuantumFurnace.Workspace(config, N3_HAM, N3_JUMPS)
    end

    @testset "Function duplication: only oft! uses the closed-form formula" begin
        # Acceptance #8: `exp(-Δ² / 4σ²)` may appear only in `oft!`
        # (`src/ofts.jl`) and the cache builder (`src/energy_domain.jl`) —
        # never inline in a hot loop again.
        src_dir = joinpath(dirname(@__DIR__), "src")
        offenders = String[]
        for (root, _, files) in walkdir(src_dir)
            for f in files
                endswith(f, ".jl") || continue
                # The two legitimate locations: skip them.
                rel = relpath(joinpath(root, f), src_dir)
                rel in ("ofts.jl", "energy_domain.jl") && continue
                # Skip src/staging — not in the active build.
                startswith(rel, "staging/") && continue
                content = read(joinpath(root, f), String)
                # Look for the live-exp pattern in non-comment lines.
                for (i, line) in enumerate(split(content, "\n"))
                    stripped = strip(line)
                    startswith(stripped, "#") && continue
                    if occursin(r"exp\(-.*4 ?\* ?(?:config\.|self\.)?sigma", line) ||
                       occursin(r"inv_4sigma2.*exp\(", line)
                        push!(offenders, "$(rel):$(i): $(stripped)")
                    end
                end
            end
        end
        @test isempty(offenders)
    end
end
