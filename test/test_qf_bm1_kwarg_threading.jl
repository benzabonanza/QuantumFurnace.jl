"""
Regression test for `allow_unpaired_nonhermitian` kwarg threading
through `run_trajectory` (qf-bm1.Q1).

The original commit 3b284a2 declared the kwarg in `run_trajectory`'s
signature but did NOT forward it to the four delegated paths
(`_run_trajectory_adaptive`, `_run_trajectory_convergence`,
`_run_trajectory_with_obs`, `_build_framework_and_seed`). The fix in
the qf-bm1 follow-up commit threads the kwarg through every helper
and `run_observable_trajectories` so the user-facing opt-out works
end-to-end. This file is the regression test guarding that fix.
"""

@testset "qf-bm1 run_trajectory kwarg threading regression" begin
    # Build a small fixture
    ham_path = joinpath(dirname(@__DIR__), "hamiltonians",
                        "heis_xxx_zzdisordered_periodic_n3.bson")
    ham = QuantumFurnace._load_hamiltonian_bson(ham_path, 5.0)
    n = 3
    dim = 2^n

    # UNPAIRED non-Hermitian jump
    sigma_plus = ComplexF64[0 1; 0 0]
    norm_fac = 1.0 / sqrt(2)
    op = Matrix(QuantumFurnace.pad_term([sigma_plus], n, 1)) .* norm_fac
    op_eb = ham.eigvecs' * op * ham.eigvecs
    unpaired = JumpOp[JumpOp(op, op_eb, false, false)]

    psi0 = ComplexF64[1.0]
    for _ in 2:dim
        push!(psi0, 0)
    end

    cfg = Config(;
        sim = Thermalize(),
        domain = TimeDomain(),
        construction = KMS(),
        num_qubits = n,
        with_linear_combination = true,
        beta = 5.0,
        sigma = 1.0/5.0,
        a = 5.0/30.0,
        s = 0.4,
        num_energy_bits = 8,
        w0 = 0.1,
        t0 = 2*pi / (2^8 * 0.1),
        num_trotter_steps_per_t0 = 4,
        delta = 0.1,
        mixing_time = 0.5,
    )

    @testset "default path (no observables, no convergence)" begin
        # CURRENT BUG: kwarg is silently ignored, so this throws.
        # When fixed, the call should not throw.
        @test try
            run_trajectory(unpaired, cfg, ham;
                psi0=psi0, ntraj=1, total_time=0.2,
                allow_unpaired_nonhermitian=true)
            true
        catch e
            isa(e, ArgumentError) && occursin("unpaired", string(e)) && false
        end
    end

    @testset "with observables path" begin
        # Pauli Z on site 1, in eigenbasis
        Z_op = Matrix(QuantumFurnace.pad_term([Z], n, 1))
        observables = Matrix{ComplexF64}[ham.eigvecs' * Z_op * ham.eigvecs]

        @test try
            run_trajectory(unpaired, cfg, ham;
                psi0=psi0, ntraj=1, total_time=0.2,
                observables=observables,
                allow_unpaired_nonhermitian=true)
            true
        catch e
            isa(e, ArgumentError) && occursin("unpaired", string(e)) && false
        end
    end

    # The default-true behaviour (validation triggers) is correct — make sure
    # we don't break that.
    @testset "default rejects unpaired NH (sanity check)" begin
        @test_throws ArgumentError run_trajectory(unpaired, cfg, ham;
            psi0=psi0, ntraj=1, total_time=0.2)
    end

    @testset "paired NH path works fine (sanity check)" begin
        # σ⁻ partner
        sigma_minus = ComplexF64[0 0; 1 0]
        op_minus = Matrix(QuantumFurnace.pad_term([sigma_minus], n, 1)) .* norm_fac
        paired = JumpOp[
            JumpOp(op, op_eb, false, false),
            JumpOp(op_minus, ham.eigvecs' * op_minus * ham.eigvecs, false, false),
        ]
        # Should not throw
        result = run_trajectory(paired, cfg, ham;
            psi0=psi0, ntraj=1, total_time=0.2)
        @test result.n_trajectories == 1
    end
end
