@testset "sweep_channel_mixing harness (qf-e4z.2)" begin
    using LinearAlgebra
    using QuantumFurnace: predict_channel_trajectory, _load_hamiltonian_bson,
        _load_channel_param_table, _lookup_channel_params, _build_channel_config,
        _jumps_in_basis, _channel_sweep_sidecar_path

    # The smoke cell for P0b: n=3, β=10, ε=1e-3, smooth-Metro KMS, TimeDomain.
    # TimeDomain avoids the qf-d0w shared-δt₀ commensurability issue (downstream
    # sweeps may pick TrotterDomain once the GQSP refactor in qf-e4z.19 lands).
    project_root = dirname(@__DIR__)
    param_table  = joinpath(project_root, "scripts", "output", "channel_param_table.bson")
    ham_path     = joinpath(project_root, "hamiltonians",
                            "heis_xxx_zzdisordered_periodic_n3.bson")

    if !isfile(param_table) || !isfile(ham_path)
        @warn "Skipping sweep_channel_mixing tests: required artefacts missing" param_table ham_path
    else
        @testset "smoke: matches direct predict_channel_trajectory" begin
            results = sweep_channel_mixing(
                [3], [10.0];
                target_epsilons = [1e-3],
                filter_kinds = [:smooth_metro],
                domain = TimeDomain(),
                construction = KMS(),
                seeds = [42],
                krylovdim = 30,
                k_grid_max_log = 4,
                k_grid_length = 40,
                output_dir = nothing,
            )
            @test length(results) == 1
            r = results[1]
            @test r.n == 3
            @test r.beta == 10.0
            @test r.eps == 1e-3
            @test r.filter === :smooth_metro
            @test r.construction == "KMS"
            @test r.domain == "Time"
            @test r.with_gqsp === true
            @test r.gqsp_degree == 1
            # qf-e4y.6: enum is now {:extrapolated, :floor, :nan}. At
            # n=3 β=10 ε=1e-3 smooth-Metro the channel asymptotic floor
            # (~3.94e-3) exceeds ε ⇒ :floor source. tau_mix is the
            # conservative log(d/ε)/λ_gap fallback and stays finite.
            @test r.tau_mix_source in (:extrapolated, :floor)
            @test isfinite(r.tau_mix)
            @test isfinite(r.lambda_gap_channel)
            @test r.lambda_gap_channel > 0.0
            @test isfinite(r.floor_distance) && r.floor_distance >= 0.0
            @test isfinite(r.oft_time_per_step) && r.oft_time_per_step > 0
            @test isfinite(r.b_time_per_step)
            @test r.per_step_time ≈ 2.0 * r.oft_time_per_step + r.b_time_per_step rtol=1e-12
            @test r.n_steps_total > 0
            @test r.total_ham_sim_time ≈ r.n_steps_total * r.per_step_time rtol=1e-12

            # Independent run with the same param-table row reproduces the gap and τ_mix exactly.
            rows = _load_channel_param_table(param_table)
            row  = _lookup_channel_params(rows, 3, 10.0, 1e-3, :smooth_metro)
            ham  = _load_hamiltonian_bson(ham_path, 10.0)
            cfg  = _build_channel_config(row, 3, 10.0, TimeDomain(), KMS())
            jumps = _jumps_in_basis(ham, 3, ham.eigvecs)
            d = size(ham.data, 1)
            rho_0 = Matrix{ComplexF64}(I(d) ./ d)
            k_grid = unique(round.(Int, exp10.(range(0, 4, length=40))))
            res_direct = predict_channel_trajectory(cfg, ham, jumps, rho_0, k_grid;
                                                     krylovdim=30)

            # Mirror the post-processing that sweep_channel_mixing performs
            # (src/lindblad_action.jl::sweep_channel_mixing).  In the
            # :extrapolated branch `r.tau_mix` is the eigenmode bisection
            # crossing, NOT `log(d/ε)/gap` — the latter is only used as the
            # :floor fallback.
            ε = 1e-3
            lambda_eff = log.(res_direct.eigenvalues) ./ res_direct.delta_used
            t_upper_ch = res_direct.spectral_gap > 0 ?
                max(res_direct.t[end], 5.0 * log(d / ε) / res_direct.spectral_gap) :
                res_direct.t[end]
            res_eig = eigenmode_mixing_time(
                lambda_eff, res_direct.c, res_direct.R_modes,
                res_direct.rho_inf, res_direct.sigma_beta, ε;
                t_upper = t_upper_ch,
            )
            tau_mix_direct = if res_eig.source === :extrapolated &&
                                 isfinite(res_eig.mixing_time) && res_eig.mixing_time > 0
                res_eig.mixing_time
            elseif res_direct.spectral_gap > 0
                log(d / ε) / res_direct.spectral_gap
            else
                NaN
            end

            @test r.tau_mix_source === res_eig.source
            @test isapprox(r.lambda_gap_channel, res_direct.spectral_gap, rtol=1e-12)
            @test isapprox(r.tau_mix, tau_mix_direct, rtol=1e-10)
            @test isapprox(r.achieved_dist_at_kmax, res_direct.distances[end], rtol=1e-12)
        end

        @testset "BSON sidecars + skip_existing" begin
            mktempdir() do tmp
                # First run: writes the sidecar.
                r1 = sweep_channel_mixing(
                    [3], [10.0];
                    target_epsilons = [1e-3],
                    filter_kinds = [:smooth_metro],
                    domain = TimeDomain(),
                    construction = KMS(),
                    seeds = [42],
                    krylovdim = 30, k_grid_max_log = 4, k_grid_length = 40,
                    output_dir = tmp,
                )
                @test length(r1) == 1
                sidecar = _channel_sweep_sidecar_path(tmp, 3, 10.0, 42, 1e-3,
                                                       :smooth_metro, "KMS", "Time")
                @test isfile(sidecar)

                # Second run with skip_existing=true: should hit the cache.
                t0 = time()
                r2 = sweep_channel_mixing(
                    [3], [10.0];
                    target_epsilons = [1e-3],
                    filter_kinds = [:smooth_metro],
                    domain = TimeDomain(),
                    construction = KMS(),
                    seeds = [42],
                    krylovdim = 30, k_grid_max_log = 4, k_grid_length = 40,
                    output_dir = tmp, skip_existing = true,
                )
                cached_wall = time() - t0

                @test length(r2) == 1
                @test r2[1].lambda_gap_channel ≈ r1[1].lambda_gap_channel rtol=1e-12
                @test r2[1].tau_mix ≈ r1[1].tau_mix rtol=1e-12
                # Cache hit should be much faster than the original run (~10s).
                @test cached_wall < 5.0
            end
        end

        @testset "(z) floor_distance matches direct ‖ρ_inf - σ_β‖_1 / 2" begin
            # qf-e4y.6: the eigenmode helper's `floor_distance` field must
            # equal a direct svdvals computation on (ρ_inf, σ_β) from the
            # channel predictor. The sweep doesn't expose those matrices on
            # the sidecar (would balloon BSON size), so the test re-runs
            # `predict_channel_trajectory` on the same fixture.
            results = sweep_channel_mixing(
                [3], [10.0];
                target_epsilons = [1e-3],
                filter_kinds = [:smooth_metro],
                domain = TimeDomain(),
                construction = KMS(),
                seeds = [42],
                krylovdim = 30, k_grid_max_log = 4, k_grid_length = 40,
                output_dir = nothing,
            )
            r = results[1]
            rows = _load_channel_param_table(param_table)
            row  = _lookup_channel_params(rows, 3, 10.0, 1e-3, :smooth_metro)
            ham  = _load_hamiltonian_bson(ham_path, 10.0)
            cfg  = _build_channel_config(row, 3, 10.0, TimeDomain(), KMS())
            jumps = _jumps_in_basis(ham, 3, ham.eigvecs)
            d = size(ham.data, 1)
            rho_0 = Matrix{ComplexF64}(I(d) ./ d)
            k_grid = unique(round.(Int, exp10.(range(0, 4, length=40))))
            res_direct = predict_channel_trajectory(cfg, ham, jumps, rho_0, k_grid;
                                                     krylovdim=30)
            floor_direct = sum(svdvals(res_direct.rho_inf .- res_direct.sigma_beta)) / 2
            @test isapprox(r.floor_distance, floor_direct; atol=1e-12, rtol=1e-12)
            @info "(z) floor_distance parity" sweep=r.floor_distance direct=floor_direct
        end

        @testset "(zz) :floor source when target_eps below channel-shift" begin
            # When ε is below ‖ρ_inf - σ_β‖_1 / 2, no t solves d(t) = ε —
            # the harness must signal :floor and populate tau_mix with the
            # conservative log(d/ε)/λ_gap bound (finite, for plotting).
            results = sweep_channel_mixing(
                [3], [10.0];
                target_epsilons = [1e-3],   # below floor at this fixture
                filter_kinds = [:smooth_metro],
                domain = TimeDomain(),
                construction = KMS(),
                seeds = [42],
                krylovdim = 30, k_grid_max_log = 4, k_grid_length = 40,
                output_dir = nothing,
            )
            r = results[1]
            # Sanity: the floor must exceed ε for this assertion to fire.
            if r.floor_distance > r.eps
                @test r.tau_mix_source === :floor
                @test isfinite(r.tau_mix) && r.tau_mix > 0
                @test isapprox(r.tau_mix,
                                log(2^r.n / r.eps) / r.lambda_gap_channel;
                                rtol=1e-12)
                @info "(zz) floor branch" floor=r.floor_distance ε=r.eps τ=r.tau_mix
            else
                @info "(zz) skipped — fixture has floor below ε; assertion moot" floor=r.floor_distance ε=r.eps
            end
        end

        @testset "multi-cell expansion" begin
            # Two filters × one (n, β) × one ε. Both should produce a finite λ.
            results = sweep_channel_mixing(
                [3], [10.0];
                target_epsilons = [1e-3],
                filter_kinds = [:smooth_metro, :gaussian],
                domain = TimeDomain(),
                construction = KMS(),
                seeds = [42],
                krylovdim = 30, k_grid_max_log = 4, k_grid_length = 40,
                output_dir = nothing,
            )
            @test length(results) == 2
            for r in results
                @test r.n == 3
                @test r.beta == 10.0
                @test r.eps == 1e-3
                @test isfinite(r.lambda_gap_channel)
                @test r.lambda_gap_channel > 0.0
                @test isfinite(r.tau_mix)
            end
            @test results[1].filter === :smooth_metro
            @test results[2].filter === :gaussian
        end
    end
end
