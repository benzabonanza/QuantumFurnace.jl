#* TOOLS --------------------------------------------------------------------------------------------------------------------
function pick_transition(config::Union{LiouvConfig, ThermalizeConfig})

    if !(config.with_linear_combination)  # Gaussian case
        @printf("Gaussian\n")
        return w -> begin
            return exp(-(w + config.gaussian_parameters[1])^2 /(2 * config.gaussian_parameters[2]^2))
        end
    end

    # sqrtA = sqrt((4 * a + 1) / 8)
    sqrtA = sqrt(config.beta / 4) * sqrt(4 * config.a + 1)
    if (config.b == 0 && config.a != 0)  # No time singularity but kinky Metro in energy
        return w -> begin
            # sqrtB = beta * abs(w + 1 / (2 * beta)) / sqrt(2)
            sqrtB = sqrt(config.beta / 4) * abs(w + config.beta * config.sigma^2 / 2)
            return exp((- 2 * sqrtA * sqrtB - config.beta * w / 2 - config.beta^2 * config.sigma^2 / 4))
        end
    elseif (config.b != 0 && config.a != 0)  # No time singularity and no kinky Metro (Glauberish)
        @printf("Smooth Metro\n")
        return w -> begin
            sqrtB = sqrt(config.beta / 4) * abs(w + config.beta * config.sigma^2 / 2)
            u_min = sqrt(config.beta * config.sigma^2 * config.b / 2)  # integral lower limit

            transition_b0 = exp((- 2 * sqrtA * sqrtB - config.beta * w / 2 - config.beta^2 * config.sigma^2 / 4))

            return (transition_b0 * (erfc(sqrtA * u_min - sqrtB / u_min) 
                + exp(4 * sqrtA * sqrtB) * erfc(sqrtA * u_min + sqrtB / u_min)) / 2)

            # return (transition_b0 * (erfc(sqrtA * sqrt(b) - sqrtB / sqrt(b)) 
            #     + exp(4 * sqrtA * sqrtB) * erfc(sqrtA * sqrt(b) + sqrtB / sqrt(b))) / 2)
        end
    elseif config.a == 0  # Time singularity and kinky Metro, i.e. simple shifted Metro
        @printf("Kinky Metro\n")
        return w -> begin
            return exp(-config.beta * max(w + config.beta * config.sigma^2 / 2, 0.0))
        end
    end
end

function create_energy_labels(num_energy_bits::Int64, w0::Float64)
    N = 2^(num_energy_bits)
    # N_labels = [0:1:Int(N/2)-1; -Int(N/2):1:-1]  twos complement order
    N_labels = [-Int(N/2):1:Int(N/2)-1;]
    energy_labels = w0 * N_labels
    # @assert maximum(energy_labels) >= 2.0  # For good results
    return energy_labels
end

function truncate_energy_labels(
    energy_labels::Vector{Float64}, 
    config::Union{LiouvConfig, ThermalizeConfig}; 
    cutoff::Float64=1e-12
    )

    transition = pick_transition(config) # we normalize with max(gamma) later
    gaussfilter(w, nu) = exp(- (w - nu)^2 / (4 * config.sigma^2)) * sqrt(1 / (config.sigma * sqrt(2 * pi)))
    integrand_lb(w, nu1, nu2) = transition(w) * gaussfilter(w, nu1) * gaussfilter(w, nu2)
    integrand_ub(w, nu1, nu2) = transition(w) * gaussfilter(w, nu1) * gaussfilter(w, nu2)

    candidate_nus = filter(w -> -0.45 <= w <= (-config.beta * config.sigma^2 / 2), [-0.45:0.05:0.0;])

    start_index = length(energy_labels) + 1
    for (nu1_candidate, nu2_candidate) in Iterators.product(candidate_nus, candidate_nus)
        found_index = findfirst(w -> abs(integrand_lb(w, nu1_candidate, nu2_candidate)) >= cutoff, energy_labels)
        if found_index !== nothing
            start_index = min(start_index, found_index)
        end
    end
    
    if start_index  > length(energy_labels)
        @warn "Lower bound cutoff not found for energies, using default range."
        return energy_labels[abs.(energy_labels) .<= 2.0]
    end

    candidate_nus = Iterators.reverse(filter(w -> (-config.beta * config.sigma^2 / 2) <= w <= 0.45, [-0.1:0.05:0.45;]))
    end_index = 0
    for (nu1_candidate, nu2_candidate) in Iterators.product(candidate_nus, candidate_nus)
        found_index = findlast(w -> abs(integrand_ub(w, nu1_candidate, nu2_candidate)) >= cutoff, energy_labels)
        if found_index !== nothing
            end_index = max(end_index, found_index)
        end
    end

    if end_index === 0
        @warn "Upper bound cutoff not found for energies, using default range."
        return energy_labels[abs.(energy_labels) .<= 2.0]
    end

    if start_index == 1 || end_index == length(energy_labels)
        @warn "No truncation was done, might want more estimating energy range."
    end

    # Symmetrize energy labels around 0.
    sym_limit = max(abs(energy_labels[start_index]), abs(energy_labels[end_index]))
    return energy_labels[abs.(energy_labels) .<= sym_limit]

    return energy_labels[start_index:end_index]
end
#* --------------------------------------------------------------------------------------------------------------------------
#* --------------------------------------------------------------------------------------------------------------------------


#* GAUSS --------------------------------------------------------------------------------------------------------------------
# function transition_gauss_vectorized(jumps::Vector{JumpOp}, hamiltonian::HamHam, energy_labels::Vector{Float64}, 
#     beta::Float64)

#     dim = size(hamiltonian.data, 1)
#     w0 = energy_labels[2] - energy_labels[1]
#     transition_gauss(w) = exp(-beta^2 * (w + 1/beta)^2 /2)

#     T = zeros(ComplexF64, dim^2, dim^2)
#     for jump in jumps
#         for w in energy_labels
#             jump_oft = oft(jump, w, hamiltonian, beta)
#             T .+= transition_gauss(w) * kron(jump_oft, conj(jump_oft))
#         end
#     end
#     return w0 * beta * T / sqrt(2 * pi)  # with OFT normalizations
# end

# function transition_bohr_gauss_from_energy_vectorized(jumps::Vector{JumpOp}, hamiltonian::HamHam, 
#     energy_labels::Vector{Float64}, beta::Float64)

#     dim = size(hamiltonian.data, 1)
    
#     T = zeros(ComplexF64, dim^2, dim^2)
#     for jump in jumps
#         for nu_1 in keys(hamiltonian.bohr_dict)
#             for nu_2 in keys(hamiltonian.bohr_dict)
#                 gaussians(w) = exp(-beta^2 * (w + 1/beta)^2 /2) * exp(-beta^2 * (w - nu_1)^2 / 4) * exp(-beta^2 * (w - nu_2)^2 / 4)
#                 for w in energy_labels
#                     A_nu_1::SparseMatrixCSC{ComplexF64} = spzeros(dim, dim)
#                     A_nu_2::SparseMatrixCSC{ComplexF64} = spzeros(dim, dim)
#                     A_nu_1[hamiltonian.bohr_dict[nu_1]] .= jump.in_eigenbasis[hamiltonian.bohr_dict[nu_1]]
#                     A_nu_2[hamiltonian.bohr_dict[nu_2]] .= jump.in_eigenbasis[bohhamiltonian.bohr_dictr_dict[nu_2]]

#                     T .+= gaussians(w) * kron(A_nu_1, conj(A_nu_2))
#                 end
#             end
#         end
#     end
#     return beta * w0 * T / sqrt(2 * pi)
# end

# function create_alpha_from_gaussians_as_for_real(energy_labels::Vector{Float64}, hamiltonian::HamHam, beta::Float64)

#     w0 = energy_labels[2] - energy_labels[1]

#     gaussian_on_A_nu1(w) = exp.(-beta^2 * (w .- hamiltonian.bohr_freqs).^2 / 4)                  # A_ij
#     gaussian_on_A_nu2_dagger(w) = adjoint(exp.(-beta^2 * (w .- hamiltonian.bohr_freqs).^2 / 4))  # (A_ik)^dagger

#     alpha = zeros(Float64, size(hamiltonian.data))
#     for w in energy_labels
#         # alpha_kj \sim sum_i  alpha_ij_ik
#         alpha .+=  (exp(-beta^2 * (w + 1/beta)^2 /2)) * gaussian_on_A_nu2_dagger(w) * gaussian_on_A_nu1(w) 
#     end
#     return beta * w0 * alpha / sqrt(8 * pi)
# end

# function create_alpha_from_gaussians_as_for_real_but_integrated(energy_labels::Vector{Float64}, 
#     hamiltonian::HamHam, beta::Float64)

#     gaussian_on_A_nu1(w) = exp.(-beta^2 * (w .- hamiltonian.bohr_freqs).^2 / 4)                  # A_ij
#     gaussian_on_A_nu2_dagger(w) = adjoint(exp.(-beta^2 * (w .- hamiltonian.bohr_freqs).^2 / 4))  # (A_ik)^dagger

#     alpha_integrand(w) = (exp(-beta^2 * (w + 1/beta)^2 /2)) * gaussian_on_A_nu2_dagger(w) * gaussian_on_A_nu1(w)
#     alpha, _ = quadgk(alpha_integrand, minimum(energy_labels), maximum(energy_labels); atol=1e-12, rtol=1e-12)
#     return beta * alpha / sqrt(8 * pi)
# end

# function create_alpha_nu1_from_gaussians(nu_2::Float64, hamiltonian::HamHam, energy_labels::Vector{Float64}, beta::Float64)

#     w0 = energy_labels[2] - energy_labels[1]
#     alpha_nu1_matrix = zeros(Float64, size(hamiltonian.data))
#     for w in energy_labels
#         alpha_nu1_matrix .+= beta * exp.(-beta^2 * (w .- hamiltonian.bohr_freqs).^2 / 4) * (exp(-beta^2 * (w + 1/beta)^2 /2)
#                                     * exp(-beta^2 * (w - nu_2)^2 / 4) / sqrt(8 * pi))
#     end
#     return w0 * alpha_nu1_matrix
# end

# function create_alpha_from_gaussians(nu_1::Float64, nu_2::Float64, energy_labels::Vector{Float64}, beta::Float64;
#     energy_cutoff_epsilon::Float64 = 1e-3)  # Truncating with this epsilon is good, any larger and it actually has an effect.

#     w0 = energy_labels[2] - energy_labels[1]
#     alpha = 0.0
#     for w in energy_labels
#         alpha +=  (exp(-beta^2 * (w + 1/beta)^2 / 2)
#                            * exp(-beta^2 * (w - nu_1)^2 / 4) * exp(-beta^2 * (w - nu_2)^2 / 4))
#     end
#     return beta * w0 * alpha / sqrt(8 * pi)
# end

# function create_alpha_from_gaussians_integrated(nu_1::Float64, nu_2::Float64, num_energy_bits::Int64, 
#     w0::Float64, beta::Float64)

#     N = 2^(num_energy_bits)
#     N_labels = [0:1:Int(N/2)-1; -Int(N/2):1:-1]
#     energy_labels = w0 * N_labels
#     energy_domain = (minimum(energy_labels), maximum(energy_labels))
#     alpha_integrand(w) = beta * (exp(-beta^2 * (w + 1/beta)^2 / 2)
#                             * exp(-beta^2 * (w - nu_1)^2 / 4) * exp(-beta^2 * (w - nu_2)^2 / 4) / sqrt(8 * pi))

#     alpha_nu1_nu2, _ = quadgk(alpha_integrand, energy_domain[1], energy_domain[2]; atol=1e-12, rtol=1e-12)

#     return alpha_nu1_nu2
# end


# function integrate_gamma_M(nu_1::Float64, nu_2::Float64, energy_labels::Vector{Float64}, beta::Float64)

#     transition_metro(w) = exp(-beta * max(w + 1/(2 * beta), 0.0))
#     integrand(w) = transition_metro(w) * exp(-beta^2 * (w - nu_1)^2 / 4) * exp(-beta^2 * (w - nu_2)^2 / 4)
#     # integrand(w) = exp(-beta^2 * (w - nu_1)^2 / 4) * exp(-beta^2 * (w - nu_2)^2 / 4)
#     w0 = energy_labels[2] - energy_labels[1]

#     resulting_alpha_M = 0.0
#     for w in energy_labels
#         integrand_w = integrand(w)
#         resulting_alpha_M += integrand_w
#     end

#     return w0 * beta * resulting_alpha_M / sqrt(2*pi)
# end

# function integrate_gamma_gauss(nu_1::Float64, nu_2::Float64, energy_labels::Vector{Float64}, beta::Float64)

#     transition_gauss(w) = exp(-beta^2 * (w + 1/beta)^2 /2)
#     integrand(w) = transition_gauss(w) * exp(-beta^2 * (w - nu_1)^2 / 4) * exp(-beta^2 * (w - nu_2)^2 / 4)
#     w0 = energy_labels[2] - energy_labels[1]

#     resulting_alpha_gauss = 0.0
#     for w in energy_labels
#         integrand_w = integrand(w)
#         resulting_alpha_gauss += integrand_w
#     end

#     return w0 * beta * resulting_alpha_gauss / sqrt(2*pi)
# end