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
    psi0 = fill(ComplexF64(1.0), N3_DIM) / sqrt(N3_DIM)
    rho0 = psi0 * psi0'

    # ------------------------------------------------------------------
    # DM regression: EnergyDomain
    # ------------------------------------------------------------------
    @testset "DM regression: EnergyDomain" begin
        ref_data = BSON.load(joinpath(ref_dir, "energy_dm_reference.bson"))
        rho_ref = ref_data[:rho]
        delta = ref_data[:delta]

        liouv_config = make_config(Lindbladian(), EnergyDomain(); num_qubits=3, construction=GNS())
        L = construct_lindbladian(N3_JUMPS, liouv_config, N3_HAM)
        rho_fresh = reshape(exp(delta * L) * vec(rho0), N3_DIM, N3_DIM)
        rho_fresh = (rho_fresh + rho_fresh') / 2

        @test isapprox(rho_fresh, rho_ref; atol=1e-10)
    end

    # ------------------------------------------------------------------
    # Trajectory regression: EnergyDomain
    # ------------------------------------------------------------------
    @testset "Trajectory regression: EnergyDomain" begin
        delta = 0.01
        seed = 12345
        ntraj = 1000

        # Compute DM reference fresh (deterministic, platform-portable)
        liouv_config = make_config(Lindbladian(), EnergyDomain(); num_qubits=3, construction=GNS())
        L = construct_lindbladian(N3_JUMPS, liouv_config, N3_HAM)
        rho_dm = reshape(exp(delta * L) * vec(rho0), N3_DIM, N3_DIM)
        rho_dm = (rho_dm + rho_dm') / 2

        # Compute trajectory average
        therm_config = make_config(Thermalize(), EnergyDomain();
            num_qubits=3, construction=GNS(), delta=delta, mixing_time=Float64(delta))
        ws = QuantumFurnace._build_trajectory_workspace(therm_config, N3_HAM, N3_JUMPS; delta=delta)

        rng = Random.Xoshiro(seed)
        rho_traj = zeros(ComplexF64, N3_DIM, N3_DIM)
        for _ in 1:ntraj
            psi = copy(psi0)
            step_along_trajectory!(psi, ws, rng)
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

        liouv_config = make_config(Lindbladian(), TrotterDomain(); num_qubits=3, construction=KMS())
        L = construct_lindbladian(N3_TROTTER_JUMPS, liouv_config, N3_HAM; trotter=N3_TROTTER)
        rho_fresh = reshape(exp(delta * L) * vec(rho0), N3_DIM, N3_DIM)
        rho_fresh = (rho_fresh + rho_fresh') / 2

        @test isapprox(rho_fresh, rho_ref; atol=1e-10)
    end

    # ------------------------------------------------------------------
    # Trajectory regression: TrotterDomain (coherent)
    # ------------------------------------------------------------------
    @testset "Trajectory regression: TrotterDomain (coherent)" begin
        delta = 0.01
        seed = 12345
        ntraj = 1000

        # Compute DM reference fresh (deterministic, platform-portable)
        liouv_config = make_config(Lindbladian(), TrotterDomain(); num_qubits=3, construction=KMS())
        L = construct_lindbladian(N3_TROTTER_JUMPS, liouv_config, N3_HAM; trotter=N3_TROTTER)
        rho_dm = reshape(exp(delta * L) * vec(rho0), N3_DIM, N3_DIM)
        rho_dm = (rho_dm + rho_dm') / 2

        # Compute trajectory average
        therm_config = make_config(Thermalize(), TrotterDomain();
            num_qubits=3, construction=KMS(), delta=delta, mixing_time=Float64(delta))
        ws = QuantumFurnace._build_trajectory_workspace(therm_config, N3_HAM, N3_TROTTER_JUMPS;
            trotter=N3_TROTTER, delta=delta)

        rng = Random.Xoshiro(seed)
        rho_traj = zeros(ComplexF64, N3_DIM, N3_DIM)
        for _ in 1:ntraj
            psi = copy(psi0)
            step_along_trajectory!(psi, ws, rng)
            rho_traj .+= psi * psi'
        end
        rho_traj ./= ntraj
        rho_traj = (rho_traj + rho_traj') / 2

        @test isapprox(rho_traj, rho_dm; atol=0.05)
    end

end
