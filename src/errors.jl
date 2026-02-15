function compute_errors(hamiltonian::HamHam, config::LiouvConfig; trotter::Union{TrottTrott, Nothing} = nothing)

    energy_labels = create_energy_labels(config.num_energy_bits, config.w0)
    truncated_energy_labels = truncate_energy_labels(energy_labels, config)

    energy_error = compute_energy_quadrature_error(config, hamiltonian, truncated_energy_labels; trotter = trotter)
    @printf("Worst quadrature error for the energy integral: %s\n", energy_error)
end

function compute_quadrature_error(integrand::Function, labels::Vector{Float64}, args...)
    integral = quadgk(t->integrand(t, args...), minimum(labels), maximum(labels); atol=1e-10, rtol=1e-10)[1]
    sum = riemann_sum(t->integrand(t, args...), labels)
    return norm(integral - sum)
end

function compute_energy_quadrature_error(config::LiouvConfig, hamiltonian::HamHam, energy_labels::Vector{Float64};
    trotter::Union{TrottTrott, Nothing} = nothing)
    transition = pick_transition(config)
    gaussian_filter(w) =  exp(-config.beta^2 * (w + 1 / (2 * config.beta))^2 / 4) * beta / sqrt(2 * pi) # Worst point at -1/2β
    jump = create_jumpop(["X"], config.num_qubits, Int(round(config.num_qubits / 2)), hamiltonian; trotter = trotter)
    jump_oft(w) = oft(jump, w, hamiltonian, config.beta)
    dm = Matrix{ComplexF64}(I(2^num_qubits) / 2^num_qubits)
    # integrand(w) = transition(w) * gaussian_filter(w)^2
    integrand(w) = transition(w) * tr(jump_oft(w) * dm * jump_oft(w)')
    # energies = [-2:0.01:2.0;]
    # display(plot(energies, integrand.(energies)))

    return compute_quadrature_error(integrand, energy_labels)
end

function compute_time_oft_quadrature_error()
end

function compute_trotter_oft_quadrature_error()
end

function compute_time_B_quadrature_error()
end

function compute_trotter_B_quadrature_error()
end
