"""
TINF-02: Regression tests with frozen BSON reference data.

Compares fresh DM and trajectory computations against frozen reference density matrices
stored in test/reference/*.bson. Any numerical drift from code changes will cause failures.

Always runs as part of Pkg.test() (fast: load BSON + recompute + compare).
"""

using BSON, Random

# Path resolution for Pkg.test() compatibility
source_root = dirname(@__DIR__)
ref_dir = joinpath(source_root, "test", "reference")

@testset "TINF-02: Regression tests" begin

    # Shared initial state for all regression tests
    psi0 = fill(ComplexF64(1.0), SMALL_DIM) / sqrt(SMALL_DIM)
    rho0 = psi0 * psi0'

    # ------------------------------------------------------------------
    # DM regression: EnergyDomain
    # ------------------------------------------------------------------
    @testset "DM regression: EnergyDomain" begin
        ref_data = BSON.load(joinpath(ref_dir, "energy_dm_reference.bson"))
        rho_ref = ref_data[:rho]
        delta = ref_data[:delta]

        liouv_config = make_small_liouv_config(EnergyDomain())
        L = construct_lindbladian(SMALL_JUMPS, liouv_config, SMALL_HAM)
        rho_fresh = reshape(exp(delta * L) * vec(rho0), SMALL_DIM, SMALL_DIM)
        rho_fresh = (rho_fresh + rho_fresh') / 2

        @test isapprox(rho_fresh, rho_ref; atol=1e-10)
    end

    # ------------------------------------------------------------------
    # Trajectory regression: EnergyDomain
    # ------------------------------------------------------------------
    @testset "Trajectory regression: EnergyDomain" begin
        ref_data = BSON.load(joinpath(ref_dir, "energy_traj_reference.bson"))
        rho_ref = ref_data[:rho]
        delta = ref_data[:delta]
        seed = ref_data[:seed]
        ntraj = ref_data[:ntraj]

        therm_config = make_small_thermalize_config(EnergyDomain();
            delta=delta, mixing_time=Float64(delta))
        precomputed = precompute_data(EnergyDomain(), therm_config, SMALL_HAM)
        scratch = KrausScratch(ComplexF64, SMALL_DIM)
        fw = build_trajectoryframework(SMALL_JUMPS, SMALL_HAM, therm_config,
            precomputed, scratch, delta)

        Random.seed!(seed)
        rho_traj = zeros(ComplexF64, SMALL_DIM, SMALL_DIM)
        for _ in 1:ntraj
            psi = copy(psi0)
            step_along_trajectory!(psi, fw)
            rho_traj .+= psi * psi'
        end
        rho_traj ./= ntraj
        rho_traj = (rho_traj + rho_traj') / 2

        @test isapprox(rho_traj, rho_ref; atol=1e-10)
    end

    # ------------------------------------------------------------------
    # DM regression: TrotterDomain (coherent)
    # ------------------------------------------------------------------
    @testset "DM regression: TrotterDomain (coherent)" begin
        ref_data = BSON.load(joinpath(ref_dir, "trotter_coherent_dm_reference.bson"))
        rho_ref = ref_data[:rho]
        delta = ref_data[:delta]

        liouv_config = make_small_liouv_config(TrotterDomain(); with_coherent=true)
        L = construct_lindbladian(SMALL_JUMPS, liouv_config, SMALL_HAM; trotter=SMALL_TROTTER)
        rho_fresh = reshape(exp(delta * L) * vec(rho0), SMALL_DIM, SMALL_DIM)
        rho_fresh = (rho_fresh + rho_fresh') / 2

        @test isapprox(rho_fresh, rho_ref; atol=1e-10)
    end

    # ------------------------------------------------------------------
    # Trajectory regression: TrotterDomain (coherent)
    # ------------------------------------------------------------------
    @testset "Trajectory regression: TrotterDomain (coherent)" begin
        ref_data = BSON.load(joinpath(ref_dir, "trotter_coherent_traj_reference.bson"))
        rho_ref = ref_data[:rho]
        delta = ref_data[:delta]
        seed = ref_data[:seed]
        ntraj = ref_data[:ntraj]

        therm_config = make_small_thermalize_config(TrotterDomain();
            with_coherent=true, delta=delta, mixing_time=Float64(delta))
        precomputed = precompute_data(TrotterDomain(), therm_config, SMALL_TROTTER)
        scratch = KrausScratch(ComplexF64, SMALL_DIM)
        fw = build_trajectoryframework(SMALL_JUMPS, SMALL_TROTTER, therm_config,
            precomputed, scratch, delta)

        Random.seed!(seed)
        rho_traj = zeros(ComplexF64, SMALL_DIM, SMALL_DIM)
        for _ in 1:ntraj
            psi = copy(psi0)
            step_along_trajectory!(psi, fw)
            rho_traj .+= psi * psi'
        end
        rho_traj ./= ntraj
        rho_traj = (rho_traj + rho_traj') / 2

        @test isapprox(rho_traj, rho_ref; atol=1e-10)
    end

end
