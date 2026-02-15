#! Changed it slightly for speed without debugging
function coherent_bohr(hamiltonian::HamHam, jump::JumpOp, config::Union{LiouvConfig, ThermalizeConfig})

    dim = size(hamiltonian.data, 1)
    unique_freqs = keys(hamiltonian.bohr_dict)

    f = pick_f(config)  # Picks rates for B in Bohr domain

    B = zeros(ComplexF64, dim, dim)
    f_A_nu_1 = zeros(ComplexF64, dim, dim)
    for nu_2 in unique_freqs
        A_nu_2::SparseMatrixCSC{ComplexF64} = spzeros(dim, dim)
        indices = hamiltonian.bohr_dict[nu_2]
        A_nu_2[indices] .= jump.in_eigenbasis[indices]
        @. f_A_nu_1 = f(hamiltonian.bohr_freqs, nu_2) * jump.in_eigenbasis

        B .+= A_nu_2' * f_A_nu_1
    end
    return B
end

function coherent_bohr(hamiltonian::HamHam, jumps::Vector{JumpOp}, config::Union{LiouvConfig, ThermalizeConfig})

    dim = size(hamiltonian.data, 1)
    unique_freqs = keys(hamiltonian.bohr_dict)

    f = pick_f(config)  # Picks rates for B in Bohr domain

    B = zeros(ComplexF64, dim, dim)
    for jump in jumps
        f_A_nu_1 = zeros(ComplexF64, dim, dim)
        for nu_2 in unique_freqs
            A_nu_2::SparseMatrixCSC{ComplexF64} = spzeros(dim, dim)
            indices = hamiltonian.bohr_dict[nu_2]
            A_nu_2[indices] .= jump.in_eigenbasis[indices]
            @. f_A_nu_1 = f(hamiltonian.bohr_freqs, nu_2) * jump.in_eigenbasis

            B .+= A_nu_2' * f_A_nu_1
        end
    end
    return B
end

function pick_f(config::Union{LiouvConfig, ThermalizeConfig})

    beta = config.beta
    sigma = config.sigma
    if config.with_linear_combination
        a = config.a
        b = config.b
        return (nu_1, nu_2) -> create_f(nu_1, nu_2, beta, sigma, a, b)
    else
        gaussian_parameters = config.gaussian_parameters
        return (nu_1, nu_2) -> create_f_gauss(nu_1, nu_2, beta, sigma, gaussian_parameters)
    end 
end

function create_f(nu_1::Float64, nu_2::Float64, beta::Float64, sigma::Float64, a::Float64, b::Float64)
    alpha = create_alpha(nu_1, nu_2, beta, sigma, a, b)
    return tanh(-beta * (nu_1 - nu_2) / 4) * alpha / (2im)
end

function create_f_gauss(nu_1::Float64, nu_2::Float64, beta::Float64, sigma::Float64, 
    gaussian_parameters::Union{Tuple{Float64, Float64}, Tuple{Nothing, Nothing}})
    """Tanh * alpha."""
    alpha_nu1_nu2 = create_alpha_gauss(nu_1, nu_2, sigma, gaussian_parameters)
    return tanh(-beta * (nu_1 - nu_2) / 4) * alpha_nu1_nu2 / (2im)
end

pick_alpha(config::LiouvConfig) = _pick_alpha_kms(config)
pick_alpha(config::ThermalizeConfig) = _pick_alpha_kms(config)
pick_alpha(config::LiouvConfigGNS) = _pick_alpha_gns(config)
pick_alpha(config::ThermalizeConfigGNS) = _pick_alpha_gns(config)

function _pick_alpha_kms(config::Union{LiouvConfig, ThermalizeConfig})

    sigma = config.sigma
    if config.with_linear_combination
        beta = config.beta
        a = config.a
        b = config.b
        return (nu_1, nu_2) -> create_alpha(nu_1, nu_2, beta, sigma, a, b)
    else
        gaussian_parameters = config.gaussian_parameters
        return (nu_1, nu_2) -> create_alpha_gauss(nu_1, nu_2, sigma, gaussian_parameters)
    end 
end

function create_alpha(nu_1::Float64, nu_2::Float64, beta::Float64, sigma::Float64, a::Float64, b::Float64)

    sqrtA = sqrt(beta * (4 * a + 1) / 4)
    sqrtB = sqrt(beta / 16) * abs(nu_1 + nu_2)
    C = beta * (nu_1 + nu_2) / 4
    prefactor = exp(a * beta^2 * sigma^2 / 2) / 2
    u_min = sqrt(beta * sigma^2 * (1 + b) / 2)
    z_plus = sqrtA * u_min + sqrtB / u_min
    z_minus = sqrtA * u_min - sqrtB / u_min

    alpha_nu_1 = (prefactor * exp(-C) * exp(-(nu_1 - nu_2)^2 / (8 * sigma^2)) * exp(- 2 * sqrtA * sqrtB) *
                    (erfc(z_minus) + exp(4 * sqrtA * sqrtB) * erfc(z_plus)))

    return alpha_nu_1
end

function _pick_alpha_gns(config::Union{LiouvConfigGNS, ThermalizeConfigGNS})

    sigma = config.sigma
    if config.with_linear_combination
        beta = config.beta
        a = config.a
        b = config.b
        return (nu_1, nu_2) -> create_alpha_gns(nu_1, nu_2, beta, sigma, a, b)
    else
        gaussian_parameters = config.gaussian_parameters  # but now β = 2ω_γ / σ_γ^2
        return (nu_1, nu_2) -> create_alpha_gauss(nu_1, nu_2, sigma, gaussian_parameters)
    end
end

"""
Coming from unshifted γ(ω) in the energy domain, leads to a partially shifted Kossakowski matrix α.
Difference vs the KMS DB case: |ν1 + ν2| → |ν1 + ν2 + β σ^2 / 2| (and thus not skew symmetric due to the Gaussian filters)
Also: No f and no b-funcs will be constructed because in this setup there is no fine-tuned B in the Lindbladian.
"""
function create_alpha_gns(nu_1::Float64, nu_2::Float64, beta::Float64, sigma::Float64, a::Float64, b::Float64)
    sqrtA = sqrt(beta * (4 * a + 1) / 4)
    sqrtB = sqrt(beta / 16) * abs(nu_1 + nu_2 + beta * sigma^2 / 2)
    C = beta * (nu_1 + nu_2) / 4
    prefactor = exp(a * beta^2 * sigma^2 / 2) / 2
    u_min = sqrt(beta * sigma^2 * (1 + b) / 2)
    z_plus = sqrtA * u_min + sqrtB / u_min
    z_minus = sqrtA * u_min - sqrtB / u_min

    alpha_nu_1 = (prefactor * exp(-C) * exp(-(nu_1 - nu_2)^2 / (8 * sigma^2)) * exp(- 2 * sqrtA * sqrtB) *
                    (erfc(z_minus) + exp(4 * sqrtA * sqrtB) * erfc(z_plus)))

    return alpha_nu_1
end

function create_alpha_gauss(
    nu_1::Float64, 
    nu_2::Float64, 
    sigma::Float64, 
    gaussian_parameters::Union{Tuple{Float64, Float64}, Tuple{Nothing, Nothing}})
    
    (w_gamma, sigma_gamma) = gaussian_parameters
    combined_sigma = sigma^2 + sigma_gamma^2
    prefactor = sigma_gamma / sqrt(combined_sigma)
    alpha_fn(nu_1) = prefactor * (exp(-(nu_1 + nu_2 + 2 * w_gamma)^2 / (8 * combined_sigma)) 
                                    * exp(-(nu_1 - nu_2)^2 / (8 * sigma^2)))
    return alpha_fn(nu_1)
end

function coherent_bohr_gauss(hamiltonian::HamHam, jump::JumpOp, beta::Float64)

    dim = size(hamiltonian.data, 1)
    unique_freqs = keys(hamiltonian.bohr_dict)

    B = zeros(ComplexF64, dim, dim)
    for nu_2 in unique_freqs
        A_nu_2::SparseMatrixCSC{ComplexF64} = spzeros(dim, dim)
        A_nu_2[hamiltonian.bohr_dict[nu_2]] .= jump.in_eigenbasis[hamiltonian.bohr_dict[nu_2]]
        f_nu1_matrix = create_f_gauss.(hamiltonian.bohr_freqs, nu_2, beta)

        B .+= A_nu_2' * (f_nu1_matrix .* jump.in_eigenbasis)
    end
    return B
end

function transition_bohr_gauss_vectorized(jumps::Vector{JumpOp}, hamiltonian::HamHam, beta::Float64; do_adjoint::Bool=false)

    dim = size(hamiltonian.data, 1)

    T = zeros(ComplexF64, dim^2, dim^2)
    for jump in jumps
        for nu_2 in keys(hamiltonian.bohr_dict)
            alpha_nu1_matrix = create_alpha_nu1_matrix(hamiltonian.bohr_freqs, nu_2, beta)

            A_nu_2::SparseMatrixCSC{ComplexF64} = spzeros(dim, dim)
            A_nu_2[hamiltonian.bohr_dict[nu_2]] .= jump.in_eigenbasis[hamiltonian.bohr_dict[nu_2]]
            A_nu_2_dagger::SparseMatrixCSC{ComplexF64} = A_nu_2'

            if !(do_adjoint)
                T .+= kron(alpha_nu1_matrix .* jump.in_eigenbasis, transpose(A_nu_2_dagger))
            else
                T .+= kron(adjoint(alpha_nu1_matrix .* jump.in_eigenbasis), transpose(A_nu_2))
            end
        end
    end
    return T
end

function transition_bohr_gauss_gibbsed_vectorized(jumps::Vector{JumpOp}, hamiltonian::HamHam, beta::Float64)

    dim = size(hamiltonian.data, 1)
    gibbs = gibbs_state_in_eigen(hamiltonian, beta)

    T = zeros(ComplexF64, dim^2, dim^2)
    for jump in jumps
        for nu_2 in keys(hamiltonian.bohr_dict)
            alpha_nu1_matrix = create_alpha_nu1_matrix(hamiltonian.bohr_freqs, nu_2, beta)

            A_nu_2::SparseMatrixCSC{ComplexF64} = spzeros(dim, dim)
            A_nu_2[hamiltonian.bohr_dict[nu_2]] .= jump.in_eigenbasis[hamiltonian.bohr_dict[nu_2]]
            A_nu_2_dagger::SparseMatrixCSC{ComplexF64} = A_nu_2'

            A_nu1s_gibbsed = gibbs^(-1/2) * (alpha_nu1_matrix .* jump.in_eigenbasis) * gibbs^(1/2)
            A_nu_2_dagger_gibbsed = gibbs^(1/2) * A_nu_2_dagger * gibbs^(-1/2)

            T .+= kron(A_nu1s_gibbsed, transpose(A_nu_2_dagger_gibbsed))
        end
    end
    return T
end

function thermalize_bohr_gauss_vectorized(jumps::Vector{JumpOp}, hamiltonian::HamHam, initial_dm::Matrix{ComplexF64}, 
    delta::Float64, mixing_time::Float64, beta::Float64)
    
    dim = size(hamiltonian.data, 1)
    num_liouv_steps = Int(round(mixing_time / delta, digits=0))
    gibbs = gibbs_state_in_eigen(hamiltonian, beta)
    gibbs_vec = vec(gibbs)
    
    time_steps = [0.0:delta:(mixing_time);]
    evolved_dm_vec = vec(copy(initial_dm))
    distances_to_gibbs = [norm(evolved_dm_vec - gibbs_vec)]

    # This implementation applies all jumps at once for one Liouvillian step.
    @showprogress dt=1 desc="Thermalize (Bohr Gaussian)..." for step in 1:num_liouv_steps

        liouv_matrix_for_step = zeros(ComplexF64, dim^2, dim^2)
        for jump in jumps
            # Coherent part
            if with_coherent
                coherent_term = coherent_bohr_gauss(hamiltonian, jump, beta)
                liouv_matrix_for_step .+= vectorize_liouvillian_coherent(coherent_term)
            end

            # Dissipative part
            for nu_2 in keys(hamiltonian.bohr_dict)
                A_nu_2::SparseMatrixCSC{ComplexF64} = spzeros(dim, dim)
                A_nu_2[hamiltonian.bohr_dict[nu_2]] .= jump.in_eigenbasis[hamiltonian.bohr_dict[nu_2]]
                A_nu_2_dagger::SparseMatrixCSC{ComplexF64} = A_nu_2'

                alpha_nu1_matrix = create_alpha_nu1_matrix(hamiltonian.bohr_freqs, nu_2, beta)

                liouv_matrix_for_step .+= vectorize_liouvillian_diss(alpha_nu1_matrix .* jump.in_eigenbasis, A_nu_2_dagger)
            end
        end

        # evolved_dm_vec = exp(delta * liouv_matrix_for_step) * evolved_dm_vec # Perfect Liouvillian evolution
        evolved_dm_vec = evolved_dm_vec + delta * liouv_matrix_for_step * evolved_dm_vec # Trotterized Liouvillian evolution
        dist = norm(evolved_dm_vec - gibbs_vec)
        push!(distances_to_gibbs, dist)
    end
    return HotAlgorithmResults(reshape(evolved_dm_vec, size(hamiltonian.data)), distances_to_gibbs, time_steps)
end

function B_nu_gauss(nu::Float64, nu_2::Float64, hamiltonian::HamHam, jump::JumpOp, beta::Float64)

    dim = size(hamiltonian.data, 1)

    B_nu = zeros(ComplexF64, dim, dim)
    nu_1 = nu + nu_2
    
    A_nu_2::SparseMatrixCSC{ComplexF64} = spzeros(dim, dim)
    A_nu_1::SparseMatrixCSC{ComplexF64} = spzeros(dim, dim)
    A_nu_2[hamiltonian.bohr_dict[nu_2]] .= jump.in_eigenbasis[hamiltonian.bohr_dict[nu_2]]
    A_nu_1[hamiltonian.bohr_dict[nu_1]] .= jump.in_eigenbasis[hamiltonian.bohr_dict[nu_1]]
    f_nu1_nu2 = create_f_gauss(nu_1, nu_2, beta)

    B_nu .= f_nu1_nu2 * A_nu_2' * (A_nu_1)

    return B_nu
end

function R_nu_gauss(nu::Float64, nu_2::Float64, hamiltonian::HamHam, jump::JumpOp, beta::Float64)

    dim = size(hamiltonian.data, 1)

    R_nu = zeros(ComplexF64, dim, dim)
    nu_1 = nu + nu_2
    
    A_nu_2::SparseMatrixCSC{ComplexF64} = spzeros(dim, dim)
    A_nu_1::SparseMatrixCSC{ComplexF64} = spzeros(dim, dim)
    A_nu_2[hamiltonian.bohr_dict[nu_2]] .= jump.in_eigenbasis[hamiltonian.bohr_dict[nu_2]]
    A_nu_1[hamiltonian.bohr_dict[nu_1]] .= jump.in_eigenbasis[hamiltonian.bohr_dict[nu_1]]
    f_nu1_nu2 = create_alpha_gauss(nu_1, nu_2, beta)

    R_nu .= f_nu1_nu2 * A_nu_2' * (A_nu_1)

    return R_nu
end

#* TOOLS --------------------------------------------------------------------------------------------------------------------
function check_alpha_skew_symmetry(alpha::Function, nu_1::Float64, nu_2::Float64, beta::Float64)
    @assert norm(alpha(nu_1, nu_2) - alpha(-nu_2, -nu_1) * exp(-beta * (nu_1 + nu_2) / 2)) < 1e-14
end

function create_bohr_dict(bohr_freqs::Matrix{Float64})
    """Creates a dictionary, where the keys are the Bohr frequencies, and the values are a list of their sparse indices 
    in the Bohr matrix. (With special care on the diagonal elements, that are identically 0.)"""

    bohr_dict = DefaultDict{Float64, Vector{CartesianIndex{2}}}(() -> CartesianIndex{2}[])
    dim = size(bohr_freqs, 1)
    bohr_dict[0.0] = CartesianIndex{2}.(1:dim, 1:dim) # nu = 0.0 is the diagonal and might be other offdiags
    for j in 1:dim
        for i in 1:(j - 1)
            push!(bohr_dict[bohr_freqs[i, j]], CartesianIndex{2}(i, j))
            push!(bohr_dict[-bohr_freqs[i, j]], CartesianIndex{2}(j, i))
        end
    end
    return bohr_dict
end

function find_all_nu1s_to_nu2(nu_2::Float64, nu::Float64, unique_freqs::Set{Float64})
    good_nu1s::Set{Float64} = Set()
    for nu_1 in unique_freqs
        # if round(nu_1 - nu_2, digits=15) == nu
        if (nu_1 - nu_2 == nu)
            push!(good_nu1s, nu_1)
        else
            continue
        end
    end
    return good_nu1s
end
#* --------------------------------------------------------------------------------------------------------------------------

