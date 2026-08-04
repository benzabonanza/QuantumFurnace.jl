using Test
using LinearAlgebra

@testset "Workspace Independence" begin
    rho_gibbs = Matrix{ComplexF64}(N3_GIBBS)
    rho_mixed = Matrix{ComplexF64}(I, N3_DIM, N3_DIM) / N3_DIM

    @testset "Retained channel workspaces own independent scratch" begin
        config = make_config(Thermalize(), EnergyDomain();
            num_qubits=3, construction=KMS(), delta=0.01, mixing_time=0.01)
        ws1 = Workspace(config, N3_HAM, N3_JUMPS)
        ws2 = Workspace(config, N3_HAM, N3_JUMPS)

        @test ws1.scratch !== ws2.scratch
        @test ws1.scratch.rho_next !== ws2.scratch.rho_next
        @test all(iszero, ws2.scratch.rho_next)
        @test length(ws1.scratch.task_scratches) == length(ws2.scratch.task_scratches)
        @test all(ws1.scratch.task_scratches[i] !== ws2.scratch.task_scratches[i]
            for i in eachindex(ws1.scratch.task_scratches))

        out1 = copy(apply_delta_channel!(ws1, copy(rho_gibbs), config, N3_HAM))
        @test all(iszero, ws2.scratch.rho_next)
        out2 = copy(apply_delta_channel!(ws2, copy(rho_mixed), config, N3_HAM))
        @test !isapprox(out1, out2; atol=1e-12)
        @test isapprox(tr(out1), 1.0; atol=1e-12)
        @test isapprox(tr(out2), 1.0; atol=1e-12)

        replay_ws = Workspace(config, N3_HAM, N3_JUMPS)
        replay = copy(apply_delta_channel!(
            replay_ws, copy(rho_gibbs), config, N3_HAM,
        ))
        @test isapprox(replay, out1; atol=1e-14, rtol=0)
    end

    @testset "Retained Lindbladian workspaces own independent scratch" begin
        config = make_config(Lindbladian(), EnergyDomain();
            num_qubits=3, construction=KMS())
        ws1 = Workspace(config, N3_HAM, N3_JUMPS)
        ws2 = Workspace(config, N3_HAM, N3_JUMPS)

        @test ws1.scratch !== ws2.scratch
        @test ws1.scratch.rho_out !== ws2.scratch.rho_out
        @test all(iszero, ws2.scratch.rho_out)

        out1 = copy(apply_lindbladian!(ws1, rho_gibbs, config, N3_HAM))
        @test all(iszero, ws2.scratch.rho_out)
        out2 = copy(apply_lindbladian!(ws2, rho_mixed, config, N3_HAM))
        @test !isapprox(out1, out2; atol=1e-12)

        replay_ws = Workspace(config, N3_HAM, N3_JUMPS)
        replay = copy(apply_lindbladian!(replay_ws, rho_gibbs, config, N3_HAM))
        @test isapprox(replay, out1; atol=1e-14, rtol=0)
    end
end
