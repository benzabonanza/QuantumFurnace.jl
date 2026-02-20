using Test
using LinearAlgebra
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
        config_kms = make_liouv_config(EnergyDomain(); with_coherent=true)
        ws = KrylovWorkspace(config_kms, TEST_HAM, TEST_JUMPS)
        @test ws.B_total !== nothing
        @test size(ws.tmp1) == (DIM, DIM)
        @test size(ws.tmp2) == (DIM, DIM)
        @test size(ws.LdagL) == (DIM, DIM)
        @test size(ws.rho_out) == (DIM, DIM)
        @test size(ws.jump_oft) == (DIM, DIM)

        config_gns = make_liouv_config_gns(EnergyDomain())
        ws_gns = KrylovWorkspace(config_gns, TEST_HAM, TEST_JUMPS)
        @test ws_gns.B_total === nothing
    end

    # ========================================================================
    # Testset 2: Round-trip matvec vs dense (EnergyDomain KMS, no coherent)
    # ========================================================================
    @testset "Round-trip: matvec vs dense (EnergyDomain KMS, no coherent)" begin
        config = make_liouv_config(EnergyDomain(); with_coherent=false)
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
        config = make_liouv_config(EnergyDomain(); with_coherent=true)
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
        config = make_liouv_config(EnergyDomain(); with_coherent=true)
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
        config = make_liouv_config(EnergyDomain(); with_coherent=true)
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
    @testset "Zero allocations in matvec hot path" begin
        config = make_liouv_config(EnergyDomain(); with_coherent=true)
        ws = KrylovWorkspace(config, TEST_HAM, TEST_JUMPS)
        rho = Matrix(random_density_matrix(NUM_QUBITS))

        # Warmup
        apply_lindbladian!(ws, rho, config, TEST_HAM)
        # Measure
        allocs = @allocated apply_lindbladian!(ws, rho, config, TEST_HAM)
        @test allocs == 0

        # Warmup adjoint
        apply_adjoint_lindbladian!(ws, rho, config, TEST_HAM)
        # Measure adjoint
        allocs_adj = @allocated apply_adjoint_lindbladian!(ws, rho, config, TEST_HAM)
        @test allocs_adj == 0
    end

    # ========================================================================
    # Phase 28: TimeDomain round-trip and allocation tests
    # ========================================================================

    # Testset 8: Round-trip matvec vs dense (TimeDomain KMS, with coherent)
    @testset "Round-trip: matvec vs dense (TimeDomain KMS, with coherent)" begin
        config = make_liouv_config(TimeDomain(); with_coherent=true)
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
        config = make_liouv_config(TimeDomain(); with_coherent=true)
        L_dense = construct_lindbladian(TEST_JUMPS, config, TEST_HAM)
        ws = KrylovWorkspace(config, TEST_HAM, TEST_JUMPS)

        for _ in 1:10
            rho = Matrix(random_density_matrix(NUM_QUBITS))
            v_adj_dense = L_dense' * vec(rho)
            L_adj_rho = apply_adjoint_lindbladian!(ws, rho, config, TEST_HAM)
            @test norm(v_adj_dense - vec(L_adj_rho)) < 1e-12
        end
    end

    # Testset 11: Zero allocations in matvec hot path (TimeDomain)
    @testset "Zero allocations in matvec hot path (TimeDomain)" begin
        config = make_liouv_config(TimeDomain(); with_coherent=true)
        ws = KrylovWorkspace(config, TEST_HAM, TEST_JUMPS)
        rho = Matrix(random_density_matrix(NUM_QUBITS))
        allocs, allocs_adj = _measure_matvec_allocs(ws, rho, config, TEST_HAM)
        @test allocs == 0
        @test allocs_adj == 0
    end

    # ========================================================================
    # Phase 28: TrotterDomain round-trip and allocation tests
    # ========================================================================

    # Testset 12: Round-trip matvec vs dense (TrotterDomain KMS, with coherent)
    @testset "Round-trip: matvec vs dense (TrotterDomain KMS, with coherent)" begin
        config = make_liouv_config(TrotterDomain(); with_coherent=true)
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
        config = make_liouv_config(TrotterDomain(); with_coherent=true)
        L_dense = construct_lindbladian(TEST_TROTTER_JUMPS, config, TEST_HAM; trotter=TEST_TROTTER)
        ws = KrylovWorkspace(config, TEST_HAM, TEST_TROTTER_JUMPS; trotter=TEST_TROTTER)

        for _ in 1:10
            rho = Matrix(random_density_matrix(NUM_QUBITS))
            v_adj_dense = L_dense' * vec(rho)
            L_adj_rho = apply_adjoint_lindbladian!(ws, rho, config, TEST_HAM)
            @test norm(v_adj_dense - vec(L_adj_rho)) < 1e-12
        end
    end

    # Testset 15: Zero allocations in matvec hot path (TrotterDomain)
    @testset "Zero allocations in matvec hot path (TrotterDomain)" begin
        config = make_liouv_config(TrotterDomain(); with_coherent=true)
        ws = KrylovWorkspace(config, TEST_HAM, TEST_TROTTER_JUMPS; trotter=TEST_TROTTER)
        rho = Matrix(random_density_matrix(NUM_QUBITS))
        allocs, allocs_adj = _measure_matvec_allocs(ws, rho, config, TEST_HAM)
        @test allocs == 0
        @test allocs_adj == 0
    end

end  # @testset "Krylov Matvec"
