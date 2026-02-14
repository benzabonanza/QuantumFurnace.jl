"""
TINF-02: Regression tests with frozen BSON reference data and DM-based trajectory comparison.

DM regression tests compare fresh DM computations against frozen reference density matrices
stored in test/reference/*.bson. Any numerical drift from code changes will cause failures.

Trajectory regression tests compare trajectory-averaged density matrices against the DM
evolution result computed fresh at test time (via exp(delta*L)). This is platform-portable
because DM evolution is deterministic, unlike frozen trajectory BSON data which depends on
BLAS internals (OpenBLAS vs Accelerate) and RNG stream behavior across platforms.

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
    # Trajectory averages are compared against DM evolution (not frozen BSON).
    # This is platform-portable because DM evolution via exp(delta*L) is
    # deterministic across all platforms and Julia versions.
    #
    # Uses delta=0.01 (not delta=0.1) to keep the Lie-Trotter splitting bias
    # small: trajectories apply Kraus operators one-by-one (product formula)
    # while DM uses exp(delta*L) (full Liouvillian), giving O(delta) systematic
    # error. At delta=0.01 the splitting bias is ~0.01 and statistical noise
    # from 1000 trajectories is also ~O(1/sqrt(1000)) ~ 0.03, so the combined
    # error is well within atol=0.05.
    #
    # Real regressions (broken Kraus operators, wrong jump sampling) produce
    # O(1) errors and will be caught easily by the 0.05 tolerance.
    @testset "Trajectory regression: EnergyDomain" begin
        delta = 0.01
        seed = 12345
        ntraj = 1000

        # Compute DM reference fresh (deterministic, platform-portable)
        liouv_config = make_small_liouv_config(EnergyDomain())
        L = construct_lindbladian(SMALL_JUMPS, liouv_config, SMALL_HAM)
        rho_dm = reshape(exp(delta * L) * vec(rho0), SMALL_DIM, SMALL_DIM)
        rho_dm = (rho_dm + rho_dm') / 2

        # Compute trajectory average
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

        @test isapprox(rho_traj, rho_dm; atol=0.05)
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
    # Trajectory averages are compared against DM evolution (not frozen BSON).
    # This is platform-portable because DM evolution via exp(delta*L) is
    # deterministic across all platforms and Julia versions.
    #
    # Uses delta=0.01 (not delta=0.1) to keep the Lie-Trotter splitting bias
    # small: trajectories apply Kraus operators one-by-one (product formula)
    # while DM uses exp(delta*L) (full Liouvillian), giving O(delta) systematic
    # error. At delta=0.01 the splitting bias is ~0.009 and statistical noise
    # from 1000 trajectories is also ~O(1/sqrt(1000)) ~ 0.03, so the combined
    # error is well within atol=0.05.
    #
    # Real regressions (broken Kraus operators, wrong jump sampling) produce
    # O(1) errors and will be caught easily by the 0.05 tolerance.
    @testset "Trajectory regression: TrotterDomain (coherent)" begin
        delta = 0.01
        seed = 12345
        ntraj = 1000

        # Compute DM reference fresh (deterministic, platform-portable)
        liouv_config = make_small_liouv_config(TrotterDomain(); with_coherent=true)
        L = construct_lindbladian(SMALL_JUMPS, liouv_config, SMALL_HAM; trotter=SMALL_TROTTER)
        rho_dm = reshape(exp(delta * L) * vec(rho0), SMALL_DIM, SMALL_DIM)
        rho_dm = (rho_dm + rho_dm') / 2

        # Compute trajectory average
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

        @test isapprox(rho_traj, rho_dm; atol=0.05)
    end

end
