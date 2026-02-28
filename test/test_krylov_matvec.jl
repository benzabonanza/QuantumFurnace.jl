using Test
using LinearAlgebra
using Random
using QuantumFurnace

# test_helpers.jl is already included by runtests.jl

# Helper to measure allocations in a hard-scoped function context.
# @testset soft-scope can cause variable boxing (176-byte spurious allocations)
# when the function under test destructures NamedTuples with complex type params.
function _measure_matvec_allocs(ws, rho, config, ham)
    # Warmup forward
    apply_lindbladian!(ws, rho, config, ham)
    apply_lindbladian!(ws, rho, config, ham)
    allocs_fwd = @allocated apply_lindbladian!(ws, rho, config, ham)

    # Warmup adjoint
    apply_adjoint_lindbladian!(ws, rho, config, ham)
    apply_adjoint_lindbladian!(ws, rho, config, ham)
    allocs_adj = @allocated apply_adjoint_lindbladian!(ws, rho, config, ham)

    return allocs_fwd, allocs_adj
end

# ============================================================================
# Round-trip correctness and allocation tests for Krylov matvec
# Phase 27: Core Matvec Infrastructure
# ============================================================================

@testset "Krylov Matvec" begin

    # ========================================================================
    # Testset 1: Workspace construction
    # ========================================================================
    @testset "Workspace construction" begin
        config_kms = make_config(Lindbladian(),EnergyDomain(); construction=KMS())
        ws = Workspace(config_kms, TEST_HAM, TEST_JUMPS)
        @test ws.B_total !== nothing
        @test size(ws.scratch.sandwich_tmp) == (DIM, DIM)
        @test size(ws.scratch.sandwich_out) == (DIM, DIM)
        @test size(ws.scratch.rho_out) == (DIM, DIM)
        @test size(ws.scratch.jump_oft) == (DIM, DIM)

        config_gns = make_config(Lindbladian(), EnergyDomain(); construction=GNS())
        ws_gns = Workspace(config_gns, TEST_HAM, TEST_JUMPS)
        @test ws_gns.B_total === nothing
    end

    # ========================================================================
    # Testset 2: Round-trip matvec vs dense (EnergyDomain KMS, no coherent)
    # ========================================================================
    @testset "Round-trip: matvec vs dense (EnergyDomain KMS, no coherent)" begin
        config = make_config(Lindbladian(),EnergyDomain(); construction=GNS())
        L_dense = construct_lindbladian(TEST_JUMPS, config, TEST_HAM)
        ws = Workspace(config, TEST_HAM, TEST_JUMPS)

        # Threshold rationale: matvec is algebraically exact (same computation, different code path).
        # Error is pure FP accumulation: O(n_jumps * DIM^2 * eps) ~ 12 * 256 * 1e-16 ~ 3e-13.
        # Threshold 1e-12 gives ~3x safety margin.
        max_err = 0.0
        for _ in 1:10
            rho = Matrix(random_density_matrix(NUM_QUBITS))
            v_dense = L_dense * vec(rho)
            L_rho = apply_lindbladian!(ws, rho, config, TEST_HAM)
            err = norm(v_dense - vec(L_rho))
            @test err < 1e-12
            max_err = max(max_err, err)
        end
        @info "Matvec round-trip (EnergyDomain GNS, no coherent)" max_error=max_err n_samples=10 threshold=1e-12
    end

    # ========================================================================
    # Testset 3: Round-trip matvec vs dense (EnergyDomain KMS, with coherent)
    # ========================================================================
    @testset "Round-trip: matvec vs dense (EnergyDomain KMS, with coherent)" begin
        config = make_config(Lindbladian(),EnergyDomain(); construction=KMS())
        L_dense = construct_lindbladian(TEST_JUMPS, config, TEST_HAM)
        ws = Workspace(config, TEST_HAM, TEST_JUMPS)

        # Threshold rationale: same as testset 2 -- algebraically exact, FP accumulation only.
        max_err = 0.0
        for _ in 1:10
            rho = Matrix(random_density_matrix(NUM_QUBITS))
            v_dense = L_dense * vec(rho)
            L_rho = apply_lindbladian!(ws, rho, config, TEST_HAM)
            err = norm(v_dense - vec(L_rho))
            @test err < 1e-12
            max_err = max(max_err, err)
        end
        @info "Matvec round-trip (EnergyDomain KMS, with coherent)" max_error=max_err n_samples=10 threshold=1e-12
    end

    # ========================================================================
    # Testset 4: Round-trip matvec vs dense (EnergyDomain GNS)
    # ========================================================================
    @testset "Round-trip: matvec vs dense (EnergyDomain GNS)" begin
        config = make_config(Lindbladian(), EnergyDomain(); construction=GNS())
        L_dense = construct_lindbladian(TEST_JUMPS, config, TEST_HAM)
        ws = Workspace(config, TEST_HAM, TEST_JUMPS)

        # Threshold rationale: same as testset 2 -- algebraically exact, FP accumulation only.
        max_err = 0.0
        for _ in 1:10
            rho = Matrix(random_density_matrix(NUM_QUBITS))
            v_dense = L_dense * vec(rho)
            L_rho = apply_lindbladian!(ws, rho, config, TEST_HAM)
            err = norm(v_dense - vec(L_rho))
            @test err < 1e-12
            max_err = max(max_err, err)
        end
        @info "Matvec round-trip (EnergyDomain GNS)" max_error=max_err n_samples=10 threshold=1e-12
    end

    # ========================================================================
    # Testset 5: Round-trip adjoint matvec vs dense adjoint (EnergyDomain KMS)
    # ========================================================================
    @testset "Round-trip: adjoint matvec vs dense adjoint (EnergyDomain KMS)" begin
        config = make_config(Lindbladian(),EnergyDomain(); construction=KMS())
        L_dense = construct_lindbladian(TEST_JUMPS, config, TEST_HAM)
        ws = Workspace(config, TEST_HAM, TEST_JUMPS)

        # Threshold rationale: adjoint matvec has same FP accumulation as forward matvec.
        max_err = 0.0
        for _ in 1:10
            rho = Matrix(random_density_matrix(NUM_QUBITS))
            v_adj_dense = L_dense' * vec(rho)
            L_adj_rho = apply_adjoint_lindbladian!(ws, rho, config, TEST_HAM)
            err = norm(v_adj_dense - vec(L_adj_rho))
            @test err < 1e-12
            max_err = max(max_err, err)
        end
        @info "Adjoint round-trip (EnergyDomain KMS)" max_error=max_err n_samples=10 threshold=1e-12
    end

    # ========================================================================
    # Testset 6: Adjoint duality check: tr(X' * L(Y)) == tr(L*(X)' * Y)
    # ========================================================================
    @testset "Adjoint duality check: tr(X' * L(Y)) == tr(L*(X)' * Y)" begin
        config = make_config(Lindbladian(),EnergyDomain(); construction=KMS())
        ws = Workspace(config, TEST_HAM, TEST_JUMPS)

        # Threshold rationale: duality identity tr(X'*L(Y)) == tr(L*(X)'*Y) is exact.
        # Error is FP accumulation from two matvecs + two traces: O(DIM^2 * n_jumps * eps) * 2.
        # 1e-11 gives ~30x margin over expected ~3e-13 per-sample error.
        max_err = 0.0
        for _ in 1:5
            X = Matrix(random_density_matrix(NUM_QUBITS))
            Y = Matrix(random_density_matrix(NUM_QUBITS))

            L_Y = copy(apply_lindbladian!(ws, Y, config, TEST_HAM))
            lhs = tr(X' * L_Y)

            Lstar_X = copy(apply_adjoint_lindbladian!(ws, X, config, TEST_HAM))
            rhs = tr(Lstar_X' * Y)

            err = abs(lhs - rhs)
            @test err < 1e-11
            max_err = max(max_err, err)
        end
        @info "Adjoint duality (EnergyDomain KMS)" max_error=max_err n_samples=5 threshold=1e-11
    end

    # ========================================================================
    # Testset 7: Zero allocations in matvec hot path
    # ========================================================================
    # Transition/alpha values are now computed via dispatched 2-arg methods
    # (pick_transition(config, w) / _pick_alpha(config, nu1, nu2)) instead of
    # stored closures, eliminating the Union{Nothing, Function} boxing overhead.
    MATVEC_ALLOC_BUDGET = 0  # bytes (EnergyDomain / BohrDomain)
    MATVEC_ALLOC_BUDGET_NUFFT = 0  # bytes (TimeDomain / TrotterDomain)

    @testset "Near-zero allocations in matvec hot path" begin
        config = make_config(Lindbladian(),EnergyDomain(); construction=KMS())
        ws = Workspace(config, TEST_HAM, TEST_JUMPS)
        rho = Matrix(random_density_matrix(NUM_QUBITS))

        # Threshold rationale: hot-path matvec must be zero-allocation for Krylov iteration performance.
        allocs, allocs_adj = _measure_matvec_allocs(ws, rho, config, TEST_HAM)
        @test allocs <= MATVEC_ALLOC_BUDGET
        @info "apply_lindbladian! allocations (EnergyDomain)" allocs_bytes=allocs threshold=MATVEC_ALLOC_BUDGET
        @test allocs_adj <= MATVEC_ALLOC_BUDGET
        @info "apply_adjoint_lindbladian! allocations (EnergyDomain)" allocs_bytes=allocs_adj threshold=MATVEC_ALLOC_BUDGET
    end

    # ========================================================================
    # Phase 28: TimeDomain round-trip and allocation tests
    # ========================================================================

    # Testset 8: Round-trip matvec vs dense (TimeDomain KMS, with coherent)
    @testset "Round-trip: matvec vs dense (TimeDomain KMS, with coherent)" begin
        config = make_config(Lindbladian(),TimeDomain(); construction=KMS())
        L_dense = construct_lindbladian(TEST_JUMPS, config, TEST_HAM)
        ws = Workspace(config, TEST_HAM, TEST_JUMPS)

        # Threshold rationale: TimeDomain uses NUFFT but dense reference uses same OFT path.
        # Error is FP accumulation only. Same 1e-12 threshold as EnergyDomain.
        max_err = 0.0
        for _ in 1:10
            rho = Matrix(random_density_matrix(NUM_QUBITS))
            v_dense = L_dense * vec(rho)
            L_rho = apply_lindbladian!(ws, rho, config, TEST_HAM)
            err = norm(v_dense - vec(L_rho))
            @test err < 1e-12
            max_err = max(max_err, err)
        end
        @info "Matvec round-trip (TimeDomain KMS, with coherent)" max_error=max_err n_samples=10 threshold=1e-12
    end

    # Testset 9: Round-trip matvec vs dense (TimeDomain GNS)
    @testset "Round-trip: matvec vs dense (TimeDomain GNS)" begin
        config = make_config(Lindbladian(), TimeDomain(); construction=GNS())
        L_dense = construct_lindbladian(TEST_JUMPS, config, TEST_HAM)
        ws = Workspace(config, TEST_HAM, TEST_JUMPS)

        # Threshold rationale: same as testset 8 -- NUFFT path, FP accumulation only.
        max_err = 0.0
        for _ in 1:10
            rho = Matrix(random_density_matrix(NUM_QUBITS))
            v_dense = L_dense * vec(rho)
            L_rho = apply_lindbladian!(ws, rho, config, TEST_HAM)
            err = norm(v_dense - vec(L_rho))
            @test err < 1e-12
            max_err = max(max_err, err)
        end
        @info "Matvec round-trip (TimeDomain GNS)" max_error=max_err n_samples=10 threshold=1e-12
    end

    # Testset 10: Round-trip adjoint matvec vs dense adjoint (TimeDomain KMS)
    @testset "Round-trip: adjoint matvec vs dense adjoint (TimeDomain KMS)" begin
        config = make_config(Lindbladian(),TimeDomain(); construction=KMS())
        L_dense = construct_lindbladian(TEST_JUMPS, config, TEST_HAM)
        ws = Workspace(config, TEST_HAM, TEST_JUMPS)

        # Threshold rationale: adjoint NUFFT path, same FP accumulation as forward.
        max_err = 0.0
        for _ in 1:10
            rho = Matrix(random_density_matrix(NUM_QUBITS))
            v_adj_dense = L_dense' * vec(rho)
            L_adj_rho = apply_adjoint_lindbladian!(ws, rho, config, TEST_HAM)
            err = norm(v_adj_dense - vec(L_adj_rho))
            @test err < 1e-12
            max_err = max(max_err, err)
        end
        @info "Adjoint round-trip (TimeDomain KMS)" max_error=max_err n_samples=10 threshold=1e-12
    end

    # Testset 11: Near-zero allocations in matvec hot path (TimeDomain)
    @testset "Near-zero allocations in matvec hot path (TimeDomain)" begin
        config = make_config(Lindbladian(),TimeDomain(); construction=KMS())
        ws = Workspace(config, TEST_HAM, TEST_JUMPS)
        rho = Matrix(random_density_matrix(NUM_QUBITS))
        # Threshold rationale: NUFFT hot-path must also be zero-allocation for Krylov performance.
        allocs, allocs_adj = _measure_matvec_allocs(ws, rho, config, TEST_HAM)
        @test allocs <= MATVEC_ALLOC_BUDGET_NUFFT
        @info "apply_lindbladian! allocations (TimeDomain)" allocs_bytes=allocs threshold=MATVEC_ALLOC_BUDGET_NUFFT
        @test allocs_adj <= MATVEC_ALLOC_BUDGET_NUFFT
        @info "apply_adjoint_lindbladian! allocations (TimeDomain)" allocs_bytes=allocs_adj threshold=MATVEC_ALLOC_BUDGET_NUFFT
    end

    # ========================================================================
    # Phase 28: TrotterDomain round-trip and allocation tests
    # ========================================================================

    # Testset 12: Round-trip matvec vs dense (TrotterDomain KMS, with coherent)
    @testset "Round-trip: matvec vs dense (TrotterDomain KMS, with coherent)" begin
        config = make_config(Lindbladian(),TrotterDomain(); construction=KMS())
        L_dense = construct_lindbladian(TEST_TROTTER_JUMPS, config, TEST_HAM; trotter=TEST_TROTTER)
        ws = Workspace(config, TEST_HAM, TEST_TROTTER_JUMPS; trotter=TEST_TROTTER)

        # Threshold rationale: TrotterDomain uses Trotter eigenbasis but same OFT arithmetic.
        # FP accumulation only, same 1e-12 threshold.
        max_err = 0.0
        for _ in 1:10
            rho = Matrix(random_density_matrix(NUM_QUBITS))
            v_dense = L_dense * vec(rho)
            L_rho = apply_lindbladian!(ws, rho, config, TEST_HAM)
            err = norm(v_dense - vec(L_rho))
            @test err < 1e-12
            max_err = max(max_err, err)
        end
        @info "Matvec round-trip (TrotterDomain KMS, with coherent)" max_error=max_err n_samples=10 threshold=1e-12
    end

    # Testset 13: Round-trip matvec vs dense (TrotterDomain GNS)
    @testset "Round-trip: matvec vs dense (TrotterDomain GNS)" begin
        config = make_config(Lindbladian(), TrotterDomain(); construction=GNS())
        L_dense = construct_lindbladian(TEST_TROTTER_JUMPS, config, TEST_HAM; trotter=TEST_TROTTER)
        ws = Workspace(config, TEST_HAM, TEST_TROTTER_JUMPS; trotter=TEST_TROTTER)

        # Threshold rationale: same as testset 12 -- Trotter eigenbasis, FP accumulation only.
        max_err = 0.0
        for _ in 1:10
            rho = Matrix(random_density_matrix(NUM_QUBITS))
            v_dense = L_dense * vec(rho)
            L_rho = apply_lindbladian!(ws, rho, config, TEST_HAM)
            err = norm(v_dense - vec(L_rho))
            @test err < 1e-12
            max_err = max(max_err, err)
        end
        @info "Matvec round-trip (TrotterDomain GNS)" max_error=max_err n_samples=10 threshold=1e-12
    end

    # Testset 14: Round-trip adjoint matvec vs dense adjoint (TrotterDomain KMS)
    @testset "Round-trip: adjoint matvec vs dense adjoint (TrotterDomain KMS)" begin
        config = make_config(Lindbladian(),TrotterDomain(); construction=KMS())
        L_dense = construct_lindbladian(TEST_TROTTER_JUMPS, config, TEST_HAM; trotter=TEST_TROTTER)
        ws = Workspace(config, TEST_HAM, TEST_TROTTER_JUMPS; trotter=TEST_TROTTER)

        # Threshold rationale: adjoint Trotter path, same FP accumulation as forward.
        max_err = 0.0
        for _ in 1:10
            rho = Matrix(random_density_matrix(NUM_QUBITS))
            v_adj_dense = L_dense' * vec(rho)
            L_adj_rho = apply_adjoint_lindbladian!(ws, rho, config, TEST_HAM)
            err = norm(v_adj_dense - vec(L_adj_rho))
            @test err < 1e-12
            max_err = max(max_err, err)
        end
        @info "Adjoint round-trip (TrotterDomain KMS)" max_error=max_err n_samples=10 threshold=1e-12
    end

    # Testset 15: Near-zero allocations in matvec hot path (TrotterDomain)
    @testset "Near-zero allocations in matvec hot path (TrotterDomain)" begin
        config = make_config(Lindbladian(),TrotterDomain(); construction=KMS())
        ws = Workspace(config, TEST_HAM, TEST_TROTTER_JUMPS; trotter=TEST_TROTTER)
        rho = Matrix(random_density_matrix(NUM_QUBITS))
        # Threshold rationale: Trotter NUFFT hot-path must also be zero-allocation.
        allocs, allocs_adj = _measure_matvec_allocs(ws, rho, config, TEST_HAM)
        @test allocs <= MATVEC_ALLOC_BUDGET_NUFFT
        @info "apply_lindbladian! allocations (TrotterDomain)" allocs_bytes=allocs threshold=MATVEC_ALLOC_BUDGET_NUFFT
        @test allocs_adj <= MATVEC_ALLOC_BUDGET_NUFFT
        @info "apply_adjoint_lindbladian! allocations (TrotterDomain)" allocs_bytes=allocs_adj threshold=MATVEC_ALLOC_BUDGET_NUFFT
    end

    # ========================================================================
    # Phase 28: BohrDomain round-trip and duality tests
    # ========================================================================

    # Testset 16: Round-trip matvec vs dense (BohrDomain KMS, with coherent)
    @testset "Round-trip: matvec vs dense (BohrDomain KMS, with coherent)" begin
        config = make_config(Lindbladian(),BohrDomain(); construction=KMS())
        L_dense = construct_lindbladian(TEST_JUMPS, config, TEST_HAM)
        ws = Workspace(config, TEST_HAM, TEST_JUMPS)

        # Threshold rationale: BohrDomain has different loop structure but same algebraic exactness.
        # FP accumulation only, same 1e-12 threshold.
        max_err = 0.0
        for _ in 1:10
            rho = Matrix(random_density_matrix(NUM_QUBITS))
            v_dense = L_dense * vec(rho)
            L_rho = apply_lindbladian!(ws, rho, config, TEST_HAM)
            err = norm(v_dense - vec(L_rho))
            @test err < 1e-12
            max_err = max(max_err, err)
        end
        @info "Matvec round-trip (BohrDomain KMS, with coherent)" max_error=max_err n_samples=10 threshold=1e-12
    end

    # Testset 17: Round-trip matvec vs dense (BohrDomain GNS)
    @testset "Round-trip: matvec vs dense (BohrDomain GNS)" begin
        config = make_config(Lindbladian(), BohrDomain(); construction=GNS())
        L_dense = construct_lindbladian(TEST_JUMPS, config, TEST_HAM)
        ws = Workspace(config, TEST_HAM, TEST_JUMPS)

        # Threshold rationale: same as testset 16 -- BohrDomain, FP accumulation only.
        max_err = 0.0
        for _ in 1:10
            rho = Matrix(random_density_matrix(NUM_QUBITS))
            v_dense = L_dense * vec(rho)
            L_rho = apply_lindbladian!(ws, rho, config, TEST_HAM)
            err = norm(v_dense - vec(L_rho))
            @test err < 1e-12
            max_err = max(max_err, err)
        end
        @info "Matvec round-trip (BohrDomain GNS)" max_error=max_err n_samples=10 threshold=1e-12
    end

    # Testset 18: Round-trip adjoint matvec vs dense adjoint (BohrDomain KMS)
    @testset "Round-trip: adjoint matvec vs dense adjoint (BohrDomain KMS)" begin
        config = make_config(Lindbladian(),BohrDomain(); construction=KMS())
        L_dense = construct_lindbladian(TEST_JUMPS, config, TEST_HAM)
        ws = Workspace(config, TEST_HAM, TEST_JUMPS)

        # Threshold rationale: adjoint BohrDomain path, same FP accumulation as forward.
        max_err = 0.0
        for _ in 1:10
            rho = Matrix(random_density_matrix(NUM_QUBITS))
            v_adj_dense = L_dense' * vec(rho)
            L_adj_rho = apply_adjoint_lindbladian!(ws, rho, config, TEST_HAM)
            err = norm(v_adj_dense - vec(L_adj_rho))
            @test err < 1e-12
            max_err = max(max_err, err)
        end
        @info "Adjoint round-trip (BohrDomain KMS)" max_error=max_err n_samples=10 threshold=1e-12
    end

    # Testset 19: Adjoint duality check (BohrDomain): tr(X' * L(Y)) == tr(L*(X)' * Y)
    @testset "Adjoint duality check (BohrDomain): tr(X' * L(Y)) == tr(L*(X)' * Y)" begin
        config = make_config(Lindbladian(),BohrDomain(); construction=KMS())
        ws = Workspace(config, TEST_HAM, TEST_JUMPS)

        # Threshold rationale: same duality identity as testset 6, BohrDomain loop structure.
        # 1e-11 gives ~30x margin over expected per-sample FP accumulation.
        max_err = 0.0
        for _ in 1:5
            X = Matrix(random_density_matrix(NUM_QUBITS))
            Y = Matrix(random_density_matrix(NUM_QUBITS))

            L_Y = copy(apply_lindbladian!(ws, Y, config, TEST_HAM))
            lhs = tr(X' * L_Y)

            Lstar_X = copy(apply_adjoint_lindbladian!(ws, X, config, TEST_HAM))
            rhs = tr(Lstar_X' * Y)

            err = abs(lhs - rhs)
            @test err < 1e-11
            max_err = max(max_err, err)
        end
        @info "Adjoint duality (BohrDomain KMS)" max_error=max_err n_samples=5 threshold=1e-11
    end

    # ========================================================================
    # Quick-35: Complex non-Hermitian jump operator round-trip tests
    # Validates that Krylov matvec matches dense kron convention for
    # general complex operators (where conj(J) rho J^T != J rho J').
    # ========================================================================

    # Create a single complex non-Hermitian jump operator for the 4-qubit test system
    let rng = MersenneTwister(42)
        raw_jump = randn(rng, ComplexF64, DIM, DIM) ./ sqrt(DIM)
        jump_in_eigen = TEST_HAM.eigvecs' * raw_jump * TEST_HAM.eigvecs
        # Not orthogonal, not Hermitian
        complex_jump = JumpOp(raw_jump, jump_in_eigen, false, false)
        complex_jumps = JumpOp[complex_jump]

        # Testset 20: Round-trip with complex jump (EnergyDomain forward)
        @testset "Round-trip: complex jump forward (EnergyDomain)" begin
            config = make_config(Lindbladian(),EnergyDomain(); construction=KMS())
            L_dense = construct_lindbladian(complex_jumps, config, TEST_HAM)
            ws = Workspace(config, TEST_HAM, complex_jumps)
            # Threshold rationale: complex non-Hermitian jump validates conj(J) convention.
            # Same algebraic exactness as real jumps, FP accumulation only.
            max_err = 0.0
            for _ in 1:10
                rho = Matrix(random_density_matrix(NUM_QUBITS))
                v_dense = L_dense * vec(rho)
                L_rho = apply_lindbladian!(ws, rho, config, TEST_HAM)
                err = norm(v_dense - vec(L_rho))
                @test err < 1e-12
                max_err = max(max_err, err)
            end
            @info "Matvec round-trip (complex jump, EnergyDomain fwd)" max_error=max_err n_samples=10 threshold=1e-12
        end

        # Testset 21: Round-trip with complex jump (EnergyDomain adjoint)
        @testset "Round-trip: complex jump adjoint (EnergyDomain)" begin
            config = make_config(Lindbladian(),EnergyDomain(); construction=KMS())
            L_dense = construct_lindbladian(complex_jumps, config, TEST_HAM)
            ws = Workspace(config, TEST_HAM, complex_jumps)
            # Threshold rationale: adjoint path for complex jumps, same FP accumulation.
            max_err = 0.0
            for _ in 1:10
                rho = Matrix(random_density_matrix(NUM_QUBITS))
                v_adj_dense = L_dense' * vec(rho)
                L_adj_rho = apply_adjoint_lindbladian!(ws, rho, config, TEST_HAM)
                err = norm(v_adj_dense - vec(L_adj_rho))
                @test err < 1e-12
                max_err = max(max_err, err)
            end
            @info "Adjoint round-trip (complex jump, EnergyDomain)" max_error=max_err n_samples=10 threshold=1e-12
        end

        # Testset 22: Adjoint duality with complex jump (EnergyDomain)
        @testset "Adjoint duality: complex jump (EnergyDomain)" begin
            config = make_config(Lindbladian(),EnergyDomain(); construction=KMS())
            ws = Workspace(config, TEST_HAM, complex_jumps)
            # Threshold rationale: duality identity with complex non-Hermitian jumps.
            # Same 1e-11 threshold as real-jump duality tests (testsets 6, 19).
            max_err = 0.0
            for _ in 1:5
                X = Matrix(random_density_matrix(NUM_QUBITS))
                Y = Matrix(random_density_matrix(NUM_QUBITS))
                L_Y = copy(apply_lindbladian!(ws, Y, config, TEST_HAM))
                lhs = tr(X' * L_Y)
                Lstar_X = copy(apply_adjoint_lindbladian!(ws, X, config, TEST_HAM))
                rhs = tr(Lstar_X' * Y)
                err = abs(lhs - rhs)
                @test err < 1e-11
                max_err = max(max_err, err)
            end
            @info "Adjoint duality (complex jump, EnergyDomain)" max_error=max_err n_samples=5 threshold=1e-11
        end

        # Testset 23: Round-trip with complex jump (TimeDomain forward + adjoint)
        @testset "Round-trip: complex jump (TimeDomain)" begin
            config_td = make_config(Lindbladian(),TimeDomain(); construction=KMS())
            L_dense_td = construct_lindbladian(complex_jumps, config_td, TEST_HAM)
            ws_td = Workspace(config_td, TEST_HAM, complex_jumps)
            # Threshold rationale: complex jump + TimeDomain NUFFT path, FP accumulation only.
            max_err_fwd = 0.0
            max_err_adj = 0.0
            for _ in 1:10
                rho = Matrix(random_density_matrix(NUM_QUBITS))
                err_fwd = norm(L_dense_td * vec(rho) - vec(apply_lindbladian!(ws_td, rho, config_td, TEST_HAM)))
                @test err_fwd < 1e-12
                max_err_fwd = max(max_err_fwd, err_fwd)
                err_adj = norm(L_dense_td' * vec(rho) - vec(apply_adjoint_lindbladian!(ws_td, rho, config_td, TEST_HAM)))
                @test err_adj < 1e-12
                max_err_adj = max(max_err_adj, err_adj)
            end
            @info "Matvec round-trip (complex jump, TimeDomain fwd)" max_error=max_err_fwd n_samples=10 threshold=1e-12
            @info "Adjoint round-trip (complex jump, TimeDomain adj)" max_error=max_err_adj n_samples=10 threshold=1e-12
        end
    end

end  # @testset "Krylov Matvec"
