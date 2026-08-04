using Test
using LinearAlgebra

# The retained density-matrix and Krylov channel paths share these construction
# helpers. Check every supported domain and both rate-scaling conventions
# directly, without depending on the archived stochastic workspace.
@testset "Retained per-jump CPTP construction" begin
    for (domain, label, trotter) in [
        (BohrDomain(), "Bohr", nothing),
        (EnergyDomain(), "Energy", nothing),
        (TimeDomain(), "Time", nothing),
        (TrotterDomain(), "Trotter", TEST_TROTTER),
    ]
        config = make_config(Thermalize(), domain; delta=TEST_DELTA)
        jumps = domain isa TrotterDomain ? TEST_TROTTER_JUMPS : TEST_JUMPS
        ham_or_trott = trotter === nothing ? TEST_HAM : trotter
        precomputed_data = QuantumFurnace._precompute_data(config, ham_or_trott)
        n_jumps = length(jumps)
        p_jump = 1.0 / n_jumps
        identity = Matrix{ComplexF64}(I, DIM, DIM)

        for (selection, rescale) in ((:sweep, false), (:random, true))
            (; K0s, U_residuals) = QuantumFurnace._precompute_per_jump_channels(
                jumps, ham_or_trott, config, precomputed_data;
                rescale_by_inv_prob=rescale,
            )
            U_coherents = QuantumFurnace._precompute_coherent_unitary(
                jumps, TEST_HAM, config, precomputed_data;
                trotter=trotter,
                delta_scale=rescale ? 1.0 / p_jump : 1.0,
            )

            @test length(K0s) == n_jumps
            @test length(U_residuals) == n_jumps
            @test length(U_coherents) == n_jumps

            builder_scratch = QuantumFurnace.ThermalizeScratch(ComplexF64, DIM)
            max_completeness_err = 0.0
            for a in eachindex(jumps)
                R_bare = copy(@inferred QuantumFurnace._precompute_R(
                    [jumps[a]], ham_or_trott, config, precomputed_data,
                    builder_scratch,
                ))
                @test R_bare isa Matrix{ComplexF64}
                @test isapprox(R_bare, R_bare'; atol=1e-12)

                R = rescale ? R_bare ./ p_jump : R_bare
                built = @inferred QuantumFurnace._build_cptp_channel(R, TEST_DELTA)
                @test isapprox(K0s[a], built.K0; atol=1e-15)
                # The eigendecomposition may choose different signs/phases for
                # degenerate eigenvectors. The channel depends only on U†U.
                @test isapprox(
                    U_residuals[a]' * U_residuals[a],
                    built.U_residual' * built.U_residual;
                    atol=1e-14,
                )
                @test U_coherents[a] !== nothing

                completeness = K0s[a]' * K0s[a] + TEST_DELTA * R +
                    U_residuals[a]' * U_residuals[a]
                err = norm(completeness - identity)
                max_completeness_err = max(max_completeness_err, err)
                @test isapprox(completeness, identity; atol=1e-10)
            end

            # The retained density-matrix substep must stay concretely inferred.
            step_return = Core.Compiler.return_type(
                QuantumFurnace._apply_one_dm_substep!,
                Tuple{
                    Matrix{ComplexF64},
                    QuantumFurnace.ThermalizeScratch{ComplexF64},
                    typeof(jumps[1]),
                    Matrix{ComplexF64},
                    Matrix{ComplexF64},
                    Matrix{ComplexF64},
                    typeof(ham_or_trott),
                    typeof(config),
                    typeof(precomputed_data),
                    Float64,
                },
            )
            @test step_return === Nothing
            @info "Retained CPTP construction ($label, $selection)" n_jumps max_completeness_err threshold_atol=1e-10
        end
    end
end
