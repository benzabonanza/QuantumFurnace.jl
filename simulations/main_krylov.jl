# ============================================================================
# Krylov spectrum simulation
# ============================================================================
#
# Demonstrates the run_krylov_spectrum API entry point, which computes the
# leading eigenvalues (and spectral gap) of the Lindbladian via matrix-free
# Krylov iteration (KrylovKit.jl). This avoids constructing the full
# dim^2 x dim^2 Lindbladian matrix.
#
# The script sets up a Heisenberg Hamiltonian with single-site Pauli jump
# operators (same pattern as main_liouv.jl / main_thermalize.jl), creates
# a Config{Lindbladian} with EnergyDomain, and calls run_krylov_spectrum.
#
# Usage:
#   julia --project=@. simulations/main_krylov.jl
#
# Phase 37, Plan 02

using Revise
includet("../src/QuantumFurnace.jl")
using .QuantumFurnace

using LinearAlgebra, Printf

function main()
    #* Config
    num_qubits = 4
    dim = 2^num_qubits
    beta = 10.0
    sigma = 1.0 / beta
    w_gamma = 1 / beta
    sigma_gamma = sqrt(2 * w_gamma / beta - sigma^2)

    # Smooth Metro
    a = 1 / 10
    b = 0.4
    eta = 0.0

    construction = KMS()
    with_linear_combination = true
    domain = EnergyDomain()
    num_energy_bits = 12
    w0 = 0.05
    t0 = 2pi / (2^num_energy_bits * w0)
    num_trotter_steps_per_t0 = 10

    # Krylov uses a Lindbladian-type config (the Lindbladian is never built
    # explicitly -- only its action on density matrices is computed via matvec).
    config = Config(
        sim = Lindbladian(),
        domain = domain,
        construction = construction,
        num_qubits = num_qubits,
        with_linear_combination = with_linear_combination,
        beta = beta,
        sigma = sigma,
        gaussian_parameters = (w_gamma, sigma_gamma),
        a = a,
        b = b,
        num_energy_bits = num_energy_bits,
        w0 = w0,
        t0 = t0,
        eta = eta,
        num_trotter_steps_per_t0 = num_trotter_steps_per_t0,
    )

    #* Hamiltonian
    hamiltonian = load_hamiltonian("heis", num_qubits; beta=beta)

    #* Trotter (only needed for TrotterDomain; nothing for EnergyDomain)
    if domain isa TrotterDomain
        trotter = TrottTrott(hamiltonian, t0, num_trotter_steps_per_t0)
        @printf("Trotter is created.\n")
    else
        trotter = nothing
    end

    #* Jumps
    jump_paulis = [[X], [Y], [Z]]
    num_of_jumps = length(jump_paulis) * num_qubits
    jump_normalization = sqrt(num_of_jumps)
    jumps::Vector{JumpOp} = []
    for pauli in jump_paulis
        for site in 1:num_qubits
            jump_op = pad_term(pauli, num_qubits, site) / jump_normalization
            basis_unitary = (domain isa TrotterDomain) ? trotter.eigvecs : hamiltonian.eigvecs
            jump_op_in_eigenbasis = basis_unitary' * jump_op * basis_unitary
            orthogonal = (jump_op == transpose(jump_op))
            hermitian = (jump_op == jump_op')
            push!(jumps, JumpOp(jump_op, jump_op_in_eigenbasis, orthogonal, hermitian))
        end
    end
    @printf("Jumps created: %d operators\n", length(jumps))

    #* Krylov spectrum
    @printf("\nRunning Krylov spectrum (nev=6, tol=1e-10, krylovdim=30)...\n")
    result = @time run_krylov_spectrum(
        jumps, config, hamiltonian, trotter;
        krylovdim = 30,
        howmany = 6,
        tol = 1e-10,
    )

    #* Results
    @printf("\nKrylov spectrum results:\n")
    @printf("  Spectral gap:    %.8e\n", result.spectral_gap)
    @printf("  Converged:       %s\n", result.converged)
    @printf("  Matvec count:    %d\n", result.matvec_count)
    @printf("  Num restarts:    %d\n", result.num_restarts)
    @printf("  Wall time:       %.2f s\n", result.metadata[:wall_time_seconds])

    @printf("\n  Leading eigenvalues:\n")
    for (i, ev) in enumerate(result.eigenvalues)
        @printf("    lambda_%d = %.8e + %.8e i\n", i, real(ev), imag(ev))
    end

    # Save (uncomment to persist)
    # save_result(result, "results/krylov_spectrum_result.bson")
end

main()
