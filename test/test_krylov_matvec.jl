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
    # Testset 1: KrylovWorkspace construction
    # ========================================================================
    @testset "KrylovWorkspace construction" begin
        config_kms = make_liouv_config(EnergyDomain(); construction=KMS())
        ws = KrylovWorkspace(config_kms, TEST_HAM, TEST_JUMPS)
        @test ws.B_total !== nothing
        @test size(ws.scratch.sandwich_tmp) == (DIM, DIM)
        @test size(ws.scratch.sandwich_out) == (DIM, DIM)
        @test size(ws.scratch.rho_out) == (DIM, DIM)
        @test size(ws.scratch.jump_oft) == (DIM, DIM)

        config_gns = make_liouv_config_gns(EnergyDomain())
        ws_gns = KrylovWorkspace(config_gns, TEST_HAM, TEST_JUMPS)
        @test ws_gns.B_total === nothing
    end

    # ========================================================================
    # Testset 2: Round-trip matvec vs dense (EnergyDomain KMS, no coherent)
    # ========================================================================
    @testset "Round-trip: matvec vs dense (EnergyDomain KMS, no coherent)" begin
        config = make_liouv_config(EnergyDomain(); construction=GNS())
        L_dense = construct_lindbladian(TEST_JUMPS, config, TEST_HAM)
        ws = KrylovWorkspace(config, TEST_HAM, TEST_JUMPS)

        for _ in 1:10
            rho = Matrix(random_density_matrix(NUM_QUBITS))
            v_dense = L_dense * vec(rho)
            L_rho = apply_lindbladian!(ws, rho, config, TEST_HAM)
            @test norm(v_dense - vec(L_rho)) < 1e-12
        end
    end

    # ========================================================================
    # Testset 3: Round-trip matvec vs dense (EnergyDomain KMS, with coherent)
    # ========================================================================
    @testset "Round-trip: matvec vs dense (EnergyDomain KMS, with coherent)" begin
        config = make_liouv_config(EnergyDomain(); construction=KMS())
        L_dense = construct_lindbladian(TEST_JUMPS, config, TEST_HAM)
        ws = KrylovWorkspace(config, TEST_HAM, TEST_JUMPS)

        for _ in 1:10
            rho = Matrix(random_density_matrix(NUM_QUBITS))
            v_dense = L_dense * vec(rho)
            L_rho = apply_lindbladian!(ws, rho, config, TEST_HAM)
            @test norm(v_dense - vec(L_rho)) < 1e-12
        end
    end

    # ========================================================================
    # Testset 4: Round-trip matvec vs dense (EnergyDomain GNS)
    # ========================================================================
    @testset "Round-trip: matvec vs dense (EnergyDomain GNS)" begin
        config = make_liouv_config_gns(EnergyDomain())
        L_dense = construct_lindbladian(TEST_JUMPS, config, TEST_HAM)
        ws = KrylovWorkspace(config, TEST_HAM, TEST_JUMPS)

        for _ in 1:10
            rho = Matrix(random_density_matrix(NUM_QUBITS))
            v_dense = L_dense * vec(rho)
            L_rho = apply_lindbladian!(ws, rho, config, TEST_HAM)
            @test norm(v_dense - vec(L_rho)) < 1e-12
        end
    end

    # ========================================================================
    # Testset 5: Round-trip adjoint matvec vs dense adjoint (EnergyDomain KMS)
    # ========================================================================
    @testset "Round-trip: adjoint matvec vs dense adjoint (EnergyDomain KMS)" begin
        config = make_liouv_config(EnergyDomain(); construction=KMS())
        L_dense = construct_lindbladian(TEST_JUMPS, config, TEST_HAM)
        ws = KrylovWorkspace(config, TEST_HAM, TEST_JUMPS)

        for _ in 1:10
            rho = Matrix(random_density_matrix(NUM_QUBITS))
            v_adj_dense = L_dense' * vec(rho)
            L_adj_rho = apply_adjoint_lindbladian!(ws, rho, config, TEST_HAM)
            @test norm(v_adj_dense - vec(L_adj_rho)) < 1e-12
        end
    end

    # ========================================================================
    # Testset 6: Adjoint duality check: tr(X' * L(Y)) == tr(L*(X)' * Y)
    # ========================================================================
    @testset "Adjoint duality check: tr(X' * L(Y)) == tr(L*(X)' * Y)" begin
        config = make_liouv_config(EnergyDomain(); construction=KMS())
        ws = KrylovWorkspace(config, TEST_HAM, TEST_JUMPS)

        for _ in 1:5
            X = Matrix(random_density_matrix(NUM_QUBITS))
            Y = Matrix(random_density_matrix(NUM_QUBITS))

            L_Y = copy(apply_lindbladian!(ws, Y, config, TEST_HAM))
            lhs = tr(X' * L_Y)

            Lstar_X = copy(apply_adjoint_lindbladian!(ws, X, config, TEST_HAM))
            rhs = tr(Lstar_X' * Y)

            @test abs(lhs - rhs) < 1e-11
        end
    end

    # ========================================================================
    # Testset 7: Zero allocations in matvec hot path
    # ========================================================================
    # Budget for Union{Nothing, Function} boxing in the function barrier.
    # The unified Workspace stores transition as Union{Nothing, Function};
    # passing it through the specialization barrier incurs a small fixed cost
    # (~300 bytes) that cannot be eliminated without adding a type parameter.
    # This threshold catches regressions (old pattern was 200K+) while allowing
    # the inherent Union boxing overhead.
    MATVEC_ALLOC_BUDGET = 512  # bytes

    @testset "Near-zero allocations in matvec hot path" begin
        config = make_liouv_config(EnergyDomain(); construction=KMS())
        ws = KrylovWorkspace(config, TEST_HAM, TEST_JUMPS)
        rho = Matrix(random_density_matrix(NUM_QUBITS))

        allocs, allocs_adj = _measure_matvec_allocs(ws, rho, config, TEST_HAM)
        @test allocs <= MATVEC_ALLOC_BUDGET
        @test allocs_adj <= MATVEC_ALLOC_BUDGET
    end

    # ========================================================================
    # Phase 28: TimeDomain round-trip and allocation tests
    # ========================================================================

    # Testset 8: Round-trip matvec vs dense (TimeDomain KMS, with coherent)
    @testset "Round-trip: matvec vs dense (TimeDomain KMS, with coherent)" begin
        config = make_liouv_config(TimeDomain(); construction=KMS())
        L_dense = construct_lindbladian(TEST_JUMPS, config, TEST_HAM)
        ws = KrylovWorkspace(config, TEST_HAM, TEST_JUMPS)

        for _ in 1:10
            rho = Matrix(random_density_matrix(NUM_QUBITS))
            v_dense = L_dense * vec(rho)
            L_rho = apply_lindbladian!(ws, rho, config, TEST_HAM)
            @test norm(v_dense - vec(L_rho)) < 1e-12
        end
    end

    # Testset 9: Round-trip matvec vs dense (TimeDomain GNS)
    @testset "Round-trip: matvec vs dense (TimeDomain GNS)" begin
        config = make_liouv_config_gns(TimeDomain())
        L_dense = construct_lindbladian(TEST_JUMPS, config, TEST_HAM)
        ws = KrylovWorkspace(config, TEST_HAM, TEST_JUMPS)

        for _ in 1:10
            rho = Matrix(random_density_matrix(NUM_QUBITS))
            v_dense = L_dense * vec(rho)
            L_rho = apply_lindbladian!(ws, rho, config, TEST_HAM)
            @test norm(v_dense - vec(L_rho)) < 1e-12
        end
    end

    # Testset 10: Round-trip adjoint matvec vs dense adjoint (TimeDomain KMS)
    @testset "Round-trip: adjoint matvec vs dense adjoint (TimeDomain KMS)" begin
        config = make_liouv_config(TimeDomain(); construction=KMS())
        L_dense = construct_lindbladian(TEST_JUMPS, config, TEST_HAM)
        ws = KrylovWorkspace(config, TEST_HAM, TEST_JUMPS)

        for _ in 1:10
            rho = Matrix(random_density_matrix(NUM_QUBITS))
            v_adj_dense = L_dense' * vec(rho)
            L_adj_rho = apply_adjoint_lindbladian!(ws, rho, config, TEST_HAM)
            @test norm(v_adj_dense - vec(L_adj_rho)) < 1e-12
        end
    end

    # Testset 11: Near-zero allocations in matvec hot path (TimeDomain)
    @testset "Near-zero allocations in matvec hot path (TimeDomain)" begin
        config = make_liouv_config(TimeDomain(); construction=KMS())
        ws = KrylovWorkspace(config, TEST_HAM, TEST_JUMPS)
        rho = Matrix(random_density_matrix(NUM_QUBITS))
        allocs, allocs_adj = _measure_matvec_allocs(ws, rho, config, TEST_HAM)
        @test allocs <= MATVEC_ALLOC_BUDGET
        @test allocs_adj <= MATVEC_ALLOC_BUDGET
    end

    # ========================================================================
    # Phase 28: TrotterDomain round-trip and allocation tests
    # ========================================================================

    # Testset 12: Round-trip matvec vs dense (TrotterDomain KMS, with coherent)
    @testset "Round-trip: matvec vs dense (TrotterDomain KMS, with coherent)" begin
        config = make_liouv_config(TrotterDomain(); construction=KMS())
        L_dense = construct_lindbladian(TEST_TROTTER_JUMPS, config, TEST_HAM; trotter=TEST_TROTTER)
        ws = KrylovWorkspace(config, TEST_HAM, TEST_TROTTER_JUMPS; trotter=TEST_TROTTER)

        for _ in 1:10
            rho = Matrix(random_density_matrix(NUM_QUBITS))
            v_dense = L_dense * vec(rho)
            L_rho = apply_lindbladian!(ws, rho, config, TEST_HAM)
            @test norm(v_dense - vec(L_rho)) < 1e-12
        end
    end

    # Testset 13: Round-trip matvec vs dense (TrotterDomain GNS)
    @testset "Round-trip: matvec vs dense (TrotterDomain GNS)" begin
        config = make_liouv_config_gns(TrotterDomain())
        L_dense = construct_lindbladian(TEST_TROTTER_JUMPS, config, TEST_HAM; trotter=TEST_TROTTER)
        ws = KrylovWorkspace(config, TEST_HAM, TEST_TROTTER_JUMPS; trotter=TEST_TROTTER)

        for _ in 1:10
            rho = Matrix(random_density_matrix(NUM_QUBITS))
            v_dense = L_dense * vec(rho)
            L_rho = apply_lindbladian!(ws, rho, config, TEST_HAM)
            @test norm(v_dense - vec(L_rho)) < 1e-12
        end
    end

    # Testset 14: Round-trip adjoint matvec vs dense adjoint (TrotterDomain KMS)
    @testset "Round-trip: adjoint matvec vs dense adjoint (TrotterDomain KMS)" begin
        config = make_liouv_config(TrotterDomain(); construction=KMS())
        L_dense = construct_lindbladian(TEST_TROTTER_JUMPS, config, TEST_HAM; trotter=TEST_TROTTER)
        ws = KrylovWorkspace(config, TEST_HAM, TEST_TROTTER_JUMPS; trotter=TEST_TROTTER)

        for _ in 1:10
            rho = Matrix(random_density_matrix(NUM_QUBITS))
            v_adj_dense = L_dense' * vec(rho)
            L_adj_rho = apply_adjoint_lindbladian!(ws, rho, config, TEST_HAM)
            @test norm(v_adj_dense - vec(L_adj_rho)) < 1e-12
        end
    end

    # Testset 15: Near-zero allocations in matvec hot path (TrotterDomain)
    @testset "Near-zero allocations in matvec hot path (TrotterDomain)" begin
        config = make_liouv_config(TrotterDomain(); construction=KMS())
        ws = KrylovWorkspace(config, TEST_HAM, TEST_TROTTER_JUMPS; trotter=TEST_TROTTER)
        rho = Matrix(random_density_matrix(NUM_QUBITS))
        allocs, allocs_adj = _measure_matvec_allocs(ws, rho, config, TEST_HAM)
        @test allocs <= MATVEC_ALLOC_BUDGET
        @test allocs_adj <= MATVEC_ALLOC_BUDGET
    end

    # ========================================================================
    # Phase 28: BohrDomain round-trip and duality tests
    # ========================================================================

    # Testset 16: Round-trip matvec vs dense (BohrDomain KMS, with coherent)
    @testset "Round-trip: matvec vs dense (BohrDomain KMS, with coherent)" begin
        config = make_liouv_config(BohrDomain(); construction=KMS())
        L_dense = construct_lindbladian(TEST_JUMPS, config, TEST_HAM)
        ws = KrylovWorkspace(config, TEST_HAM, TEST_JUMPS)

        for _ in 1:10
            rho = Matrix(random_density_matrix(NUM_QUBITS))
            v_dense = L_dense * vec(rho)
            L_rho = apply_lindbladian!(ws, rho, config, TEST_HAM)
            @test norm(v_dense - vec(L_rho)) < 1e-12
        end
    end

    # Testset 17: Round-trip matvec vs dense (BohrDomain GNS)
    @testset "Round-trip: matvec vs dense (BohrDomain GNS)" begin
        config = make_liouv_config_gns(BohrDomain())
        L_dense = construct_lindbladian(TEST_JUMPS, config, TEST_HAM)
        ws = KrylovWorkspace(config, TEST_HAM, TEST_JUMPS)

        for _ in 1:10
            rho = Matrix(random_density_matrix(NUM_QUBITS))
            v_dense = L_dense * vec(rho)
            L_rho = apply_lindbladian!(ws, rho, config, TEST_HAM)
            @test norm(v_dense - vec(L_rho)) < 1e-12
        end
    end

    # Testset 18: Round-trip adjoint matvec vs dense adjoint (BohrDomain KMS)
    @testset "Round-trip: adjoint matvec vs dense adjoint (BohrDomain KMS)" begin
        config = make_liouv_config(BohrDomain(); construction=KMS())
        L_dense = construct_lindbladian(TEST_JUMPS, config, TEST_HAM)
        ws = KrylovWorkspace(config, TEST_HAM, TEST_JUMPS)

        for _ in 1:10
            rho = Matrix(random_density_matrix(NUM_QUBITS))
            v_adj_dense = L_dense' * vec(rho)
            L_adj_rho = apply_adjoint_lindbladian!(ws, rho, config, TEST_HAM)
            @test norm(v_adj_dense - vec(L_adj_rho)) < 1e-12
        end
    end

    # Testset 19: Adjoint duality check (BohrDomain): tr(X' * L(Y)) == tr(L*(X)' * Y)
    @testset "Adjoint duality check (BohrDomain): tr(X' * L(Y)) == tr(L*(X)' * Y)" begin
        config = make_liouv_config(BohrDomain(); construction=KMS())
        ws = KrylovWorkspace(config, TEST_HAM, TEST_JUMPS)

        for _ in 1:5
            X = Matrix(random_density_matrix(NUM_QUBITS))
            Y = Matrix(random_density_matrix(NUM_QUBITS))

            L_Y = copy(apply_lindbladian!(ws, Y, config, TEST_HAM))
            lhs = tr(X' * L_Y)

            Lstar_X = copy(apply_adjoint_lindbladian!(ws, X, config, TEST_HAM))
            rhs = tr(Lstar_X' * Y)

            @test abs(lhs - rhs) < 1e-11
        end
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
            config = make_liouv_config(EnergyDomain(); construction=KMS())
            L_dense = construct_lindbladian(complex_jumps, config, TEST_HAM)
            ws = KrylovWorkspace(config, TEST_HAM, complex_jumps)
            for _ in 1:10
                rho = Matrix(random_density_matrix(NUM_QUBITS))
                v_dense = L_dense * vec(rho)
                L_rho = apply_lindbladian!(ws, rho, config, TEST_HAM)
                @test norm(v_dense - vec(L_rho)) < 1e-12
            end
        end

        # Testset 21: Round-trip with complex jump (EnergyDomain adjoint)
        @testset "Round-trip: complex jump adjoint (EnergyDomain)" begin
            config = make_liouv_config(EnergyDomain(); construction=KMS())
            L_dense = construct_lindbladian(complex_jumps, config, TEST_HAM)
            ws = KrylovWorkspace(config, TEST_HAM, complex_jumps)
            for _ in 1:10
                rho = Matrix(random_density_matrix(NUM_QUBITS))
                v_adj_dense = L_dense' * vec(rho)
                L_adj_rho = apply_adjoint_lindbladian!(ws, rho, config, TEST_HAM)
                @test norm(v_adj_dense - vec(L_adj_rho)) < 1e-12
            end
        end

        # Testset 22: Adjoint duality with complex jump (EnergyDomain)
        @testset "Adjoint duality: complex jump (EnergyDomain)" begin
            config = make_liouv_config(EnergyDomain(); construction=KMS())
            ws = KrylovWorkspace(config, TEST_HAM, complex_jumps)
            for _ in 1:5
                X = Matrix(random_density_matrix(NUM_QUBITS))
                Y = Matrix(random_density_matrix(NUM_QUBITS))
                L_Y = copy(apply_lindbladian!(ws, Y, config, TEST_HAM))
                lhs = tr(X' * L_Y)
                Lstar_X = copy(apply_adjoint_lindbladian!(ws, X, config, TEST_HAM))
                rhs = tr(Lstar_X' * Y)
                @test abs(lhs - rhs) < 1e-11
            end
        end

        # Testset 23: Round-trip with complex jump (TimeDomain forward + adjoint)
        @testset "Round-trip: complex jump (TimeDomain)" begin
            config_td = make_liouv_config(TimeDomain(); construction=KMS())
            L_dense_td = construct_lindbladian(complex_jumps, config_td, TEST_HAM)
            ws_td = KrylovWorkspace(config_td, TEST_HAM, complex_jumps)
            for _ in 1:10
                rho = Matrix(random_density_matrix(NUM_QUBITS))
                @test norm(L_dense_td * vec(rho) - vec(apply_lindbladian!(ws_td, rho, config_td, TEST_HAM))) < 1e-12
                @test norm(L_dense_td' * vec(rho) - vec(apply_adjoint_lindbladian!(ws_td, rho, config_td, TEST_HAM))) < 1e-12
            end
        end
    end

end  # @testset "Krylov Matvec"
