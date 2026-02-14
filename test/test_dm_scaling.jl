using Test
using LinearAlgebra

# DMTST-03: Single DM Euler step error scales as O(delta^2)
# DMTST-04: Multi-step DM Euler error accumulates as O(delta)
#
# Uses the Liouvillian matrix L directly for deterministic DM evolution.
# Euler step: rho(t+delta) = rho(t) + delta * reshape(L * vec(rho(t)), DIM, DIM)
# Reference: exact matrix exponential exp(delta * L).

@testset "DMTST-03: Single-step Euler error O(delta^2)" begin
    # Build Liouvillian for EnergyDomain (simplest non-trivial domain)
    config = make_liouv_config(EnergyDomain())
    L = construct_lindbladian(TEST_JUMPS, config, TEST_HAM)

    # Initial state: maximally mixed
    rho0 = Matrix{ComplexF64}(I, DIM, DIM) / DIM
    rho0_vec = vec(rho0)

    # Delta sweep
    deltas = [0.1, 0.05, 0.025, 0.0125]
    errors = Float64[]

    for delta in deltas
        # Euler step
        rho_euler = rho0 + delta * reshape(L * rho0_vec, DIM, DIM)
        # Exact step via matrix exponential
        rho_exact = reshape(exp(delta * L) * rho0_vec, DIM, DIM)
        # Error (Frobenius norm of vectorized difference)
        err = norm(vec(rho_euler) - vec(rho_exact))
        push!(errors, err)
    end

    # Compute ratios of consecutive errors
    ratios = [errors[i] / errors[i+1] for i in 1:length(errors)-1]

    # Diagnostics
    println("DMTST-03 Single-step errors: ", errors)
    println("DMTST-03 Ratios (expect ~4 for O(delta^2)): ", ratios)

    # For O(delta^2), when delta halves the error should decrease by factor ~4
    for (i, ratio) in enumerate(ratios)
        @test 3.0 <= ratio <= 5.0
    end
end

@testset "DMTST-04: Multi-step accumulated error O(delta)" begin
    # Build Liouvillian for EnergyDomain
    config = make_liouv_config(EnergyDomain())
    L = construct_lindbladian(TEST_JUMPS, config, TEST_HAM)

    # Initial state: maximally mixed
    rho0 = Matrix{ComplexF64}(I, DIM, DIM) / DIM
    rho0_vec = vec(rho0)

    # Fixed total evolution time
    T = 0.5

    # Exact evolution (compute once)
    rho_exact_T = reshape(exp(T * L) * rho0_vec, DIM, DIM)

    # Delta sweep
    deltas = [0.1, 0.05, 0.025, 0.0125]
    errors = Float64[]

    for delta in deltas
        num_steps = Int(round(T / delta))
        rho = copy(rho0)
        for _ in 1:num_steps
            rho .+= delta .* reshape(L * vec(rho), DIM, DIM)
        end
        err = norm(vec(rho) - vec(rho_exact_T))
        push!(errors, err)
    end

    # Compute ratios
    ratios = [errors[i] / errors[i+1] for i in 1:length(errors)-1]

    # Diagnostics
    println("DMTST-04 Multi-step errors: ", errors)
    println("DMTST-04 Ratios (expect ~2 for O(delta)): ", ratios)

    # For O(delta), when delta halves the error should decrease by factor ~2
    for (i, ratio) in enumerate(ratios)
        @test 1.5 <= ratio <= 2.5
    end
end
