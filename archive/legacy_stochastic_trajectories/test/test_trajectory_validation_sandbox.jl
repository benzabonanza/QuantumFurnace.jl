# test/test_trajectory_validation_sandbox.jl
#
# Sandbox shadow of test/trajectory_validation/*.jl (qf-x56.6). The heavy
# trajectory_validation scripts run a 3-δ slope sweep (50k trajectories at
# each δ) for EnergyDomain / TimeDomain / TrotterDomain plus a 1/√N MC
# convergence test (500k reference + 10 batches across four ntraj points).
# This shadow keeps only the single-step DM ↔ trajectory equivalence check
# at the smallest sandbox-affordable scale:
#
#   - EnergyDomain (the simplest of the three single-step paths)
#   - n = 3 (the smallest available disordered-Heisenberg fixture — there
#     is no n=2 BSON; the issue's "n=2" was an aspirational floor)
#   - ntraj = 5000 (10× smaller than the heavy 50k)
#   - δ = 0.1 (single point — the slope sweep stays NO_SANDBOX)
#   - construction = KMS() (the production target). The heavy TVAL-02 uses
#     GNS() for non-Trotter domains via `with_coherent = false`; KMS adds
#     the −i[B, ·] Lamb-shift kick to the trajectory, so the sandbox covers
#     the coherent-on path required by the project rule "coherent term ON
#     by default" (.claude/rules/julia-code.md).
#
# At ntraj = 5000 the MC error scales as 1/√N ≈ 1.4e-2; combined with the
# O(δ²·‖L‖²) splitting-error bias (~ 1e-2 at δ = 0.1 on the n = 3 fixture)
# the trace distance against the dense DM is ~ 1.5–2e-2 typical. Threshold
# 4e-2 ≈ 2× the typical value gives a regression-guard margin without
# admitting a qualitatively broken pipeline. Empirical at seed = 42:
# dist ≈ 0.030.

using LinearAlgebra: I, Hermitian, tr
using Random
using Test
using QuantumFurnace


@testset "Trajectory ↔ DM equivalence [sandbox shadow] (qf-x56.6)" begin

    # -----------------------------------------------------------------------
    # Single-step DM vs trajectory cross-validation on EnergyDomain.
    # PHYSICS CHECK: at δ = 0.1, the trajectory channel and the dense
    # ρ_dm = exp(δ·L)·ρ_0 differ by O(δ²·‖L‖²) ~ 1e-2 (Lie-Trotter
    # splitting bias) plus MC noise O(1/√ntraj) ~ 1.4e-2 at ntraj = 5000;
    # the trace distance is ~ 1.5–2e-2 typical. Threshold 4e-2 ≈ 2× typical
    # to flag genuine pipeline breaks while staying statistically tolerant
    # at this modest ntraj.
    # -----------------------------------------------------------------------
    @testset "(t1) Single-step DM ↔ trajectory: EnergyDomain @ n=3, δ=0.1" begin
        dim = N3_DIM
        delta = 0.1
        ntraj = 5_000
        seed = 42

        # Liouvillian + dense DM evolution.
        liouv_config = make_config(Lindbladian(), EnergyDomain();
                                    num_qubits = 3, construction = KMS())
        L = construct_lindbladian(N3_JUMPS, liouv_config, N3_HAM)

        psi0 = fill(ComplexF64(1.0), dim) / sqrt(dim)
        rho0 = psi0 * psi0'
        rho_dm = reshape(exp(delta * L) * vec(rho0), dim, dim)
        rho_dm = Hermitian((rho_dm + rho_dm') / 2)

        # Trajectory workspace + ntraj-averaged ρ.
        therm_config = make_config(Thermalize(), EnergyDomain();
                                    num_qubits = 3, construction = KMS(),
                                    delta = delta, mixing_time = delta)
        ws = QuantumFurnace._build_trajectory_workspace(
            therm_config, N3_HAM, N3_JUMPS; trotter = nothing, delta = delta)

        rho_traj = zeros(ComplexF64, dim, dim)
        rng = Random.Xoshiro(seed)
        for _ in 1:ntraj
            psi = copy(psi0)
            step_along_trajectory!(psi, ws, rng)
            rho_traj .+= psi * psi'
        end
        rho_traj ./= ntraj
        rho_traj_h = Hermitian((rho_traj + rho_traj') / 2)
        rho_traj_h = Hermitian(Matrix(rho_traj_h) ./ tr(rho_traj_h))

        dist = trace_distance_h(rho_dm, rho_traj_h)
        @test dist < 4e-2
        @info "(t1) EnergyDomain single-step DM ↔ trajectory" delta ntraj dist threshold=4e-2
    end
end
