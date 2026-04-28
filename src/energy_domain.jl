#* TOOLS --------------------------------------------------------------------------------------------------------------------

pick_transition(config::Config{<:Any, <:Any, KMS}) = _pick_transition_kms(config)
pick_transition(config::Config{<:Any, <:Any, GNS}) = _pick_transition_gns(config)

# 2-arg forms: compute transition value directly via dispatch (zero allocation on hot path)
function pick_transition(config::Config{<:Any, <:Any, KMS}, w::Real)
    if !(config.with_linear_combination)
        return exp(-(w + config.gaussian_parameters[1])^2 / (2 * config.gaussian_parameters[2]^2))
    end
    sqrtA = sqrt(config.beta / 4) * sqrt(4 * config.a + 1)
    sqrtB = sqrt(config.beta / 4) * abs(w + config.beta * config.sigma^2 / 2)
    if config.s == 0 && config.a != 0
        return exp(-2 * sqrtA * sqrtB - config.beta * w / 2 - config.beta^2 * config.sigma^2 / 4)
    elseif config.s != 0 && config.a != 0
        u_min = sqrt(config.beta * config.sigma^2 * config.s / 2)
        transition_b0 = exp(-2 * sqrtA * sqrtB - config.beta * w / 2 - config.beta^2 * config.sigma^2 / 4)
        return transition_b0 * (erfc(sqrtA * u_min - sqrtB / u_min) + exp(4 * sqrtA * sqrtB) * erfc(sqrtA * u_min + sqrtB / u_min)) / 2
    elseif config.a == 0
        return exp(-config.beta * max(w + config.beta * config.sigma^2 / 2, 0.0))
    end
end

function pick_transition(config::Config{<:Any, <:Any, GNS}, w::Real)
    if !(config.with_linear_combination)
        w_gamma = config.gaussian_parameters[1]
        sigma_gamma = config.gaussian_parameters[2]
        return exp(-(w + w_gamma)^2 / (2 * sigma_gamma^2))
    end
    sqrtA = sqrt(config.beta / 4) * sqrt(4 * config.a + 1)
    if config.s == 0 && config.a != 0
        sqrtB = sqrt(config.beta / 4) * abs(w)
        return exp(-2 * sqrtA * sqrtB - config.beta * w / 2)
    elseif config.s != 0 && config.a != 0
        sqrtB = sqrt(config.beta / 4) * abs(w)
        u_min = sqrt(config.beta * config.sigma^2 * config.s / 2)
        transition_b0 = exp(-2 * sqrtA * sqrtB - config.beta * w / 2)
        return transition_b0 * (erfc(sqrtA * u_min - sqrtB / u_min) + exp(4 * sqrtA * sqrtB) * erfc(sqrtA * u_min + sqrtB / u_min)) / 2
    elseif config.a == 0
        return exp(-config.beta * max(w, 0.0))
    end
end


function _pick_transition_kms(config::Config{<:Any, <:Any, KMS})

    if !(config.with_linear_combination)  # Gaussian case
        # @printf("Gaussian\n")
        return w -> begin
            return exp(-(w + config.gaussian_parameters[1])^2 /(2 * config.gaussian_parameters[2]^2))
        end
    end

    # sqrtA = sqrt((4 * a + 1) / 8)
    sqrtA = sqrt(config.beta / 4) * sqrt(4 * config.a + 1)
    if (config.s == 0 && config.a != 0)  # a-regularized, no smoothing (s=0)
        return w -> begin
            # sqrtB = beta * abs(w + 1 / (2 * beta)) / sqrt(2)
            sqrtB = sqrt(config.beta / 4) * abs(w + config.beta * config.sigma^2 / 2)
            return exp((- 2 * sqrtA * sqrtB - config.beta * w / 2 - config.beta^2 * config.sigma^2 / 4))
        end
    elseif (config.s != 0 && config.a != 0)  # a-regularized + smoothed (s>0)
        # @printf("Smooth Metro\n")
        return w -> begin
            sqrtB = sqrt(config.beta / 4) * abs(w + config.beta * config.sigma^2 / 2)
            u_min = sqrt(config.beta * config.sigma^2 * config.s / 2)  # integral lower limit

            transition_b0 = exp((- 2 * sqrtA * sqrtB - config.beta * w / 2 - config.beta^2 * config.sigma^2 / 4))

            return (transition_b0 * (erfc(sqrtA * u_min - sqrtB / u_min)
                + exp(4 * sqrtA * sqrtB) * erfc(sqrtA * u_min + sqrtB / u_min)) / 2)

            # return (transition_b0 * (erfc(sqrtA * sqrt(s) - sqrtB / sqrt(s))
            #     + exp(4 * sqrtA * sqrtB) * erfc(sqrtA * sqrt(s) + sqrtB / sqrt(s))) / 2)
        end
    elseif config.a == 0  # Time singularity and kinky Metro, i.e. simple shifted Metro
        # @printf("Kinky Metro\n")
        return w -> begin
            return exp(-config.beta * max(w + config.beta * config.sigma^2 / 2, 0.0))
        end
    end
end

"""
    _pick_transition_gns(config) -> (ω::Real -> Real)

    Return the (approx.) GNS-detailed-balance transition weight (\tilde{γ}(ω)).

    KMS-conditioned rates (CKBG23, Eq. (1.12)):

        \tilde{γ}(ω) = \tilde{γ}(-ω) e^{- β ω).

    In contrast, in the exact KMS-DB (CKG) construction the weight used in the ω-integral is a
    βσ_E²/2-shifted version of such a (\tilde{ω}) (Ramkumar–Soleimanifar Lemma 7.1):

        γ(ω) = \tilde{γ}(ω + βσ_E^2/2).

    This helper returns the unshifted (\tilde{γ}) (i.e., the version that itself satisfies KMS condition).
"""
function _pick_transition_gns(config::Config{<:Any, <:Any, GNS})

    # Gaussian case
    # KMS condition satisfied at inverse temperature β requires β = 2\omega_\gamma/\sigma_\gamma^2.
    if !(config.with_linear_combination)
        # @printf("Gaussian approx GNS gamma\n")
        return w -> begin
            # gaussian_parameters = [ωγ, σγ]
            w_gamma = config.gaussian_parameters[1]
            sigma_gamma = config.gaussian_parameters[2]
            return exp(-(w + w_gamma)^2 / (2 * sigma_gamma^2))
        end
    end

    sqrtA = sqrt(config.beta / 4) * sqrt(4 * config.a + 1)
    if (config.s == 0 && config.a != 0)
        # Smooth Metropolis (no time singularity) — UN-SHIFTED.
        return w -> begin
            sqrtB = sqrt(config.beta / 4) * abs(w)
            return exp((-2 * sqrtA * sqrtB - config.beta * w / 2))
        end
    elseif (config.s != 0 && config.a != 0)
        # Smooth Metropolis (Glauber-like smoothing) — UN-SHIFTED.
        # @printf("Smooth Metro approx GNS gamma \n")
        return w -> begin
            sqrtB = sqrt(config.beta / 4) * abs(w)
            u_min = sqrt(config.beta * config.sigma^2 * config.s / 2)  # integral lower limit

            transition_b0 = exp((-2 * sqrtA * sqrtB - config.beta * w / 2))

            return (transition_b0 * (erfc(sqrtA * u_min - sqrtB / u_min)
                + exp(4 * sqrtA * sqrtB) * erfc(sqrtA * u_min + sqrtB / u_min)) / 2)
        end
    elseif config.a == 0
        # Kinky Metropolis — UN-SHIFTED.
        # @printf("Kinky Metro approx GNS gamma\n")
        return w -> begin
            return exp(-config.beta * max(w, 0.0))
        end
    end
end


function _create_energy_labels(num_energy_bits::Integer, w0::Real)
    N = 2^(num_energy_bits)
    # N_labels = [0:1:Int(N/2)-1; -Int(N/2):1:-1]  twos complement order
    N_labels = [-Int(N/2):1:Int(N/2)-1;]
    energy_labels = w0 * N_labels
    # @assert maximum(energy_labels) >= 2.0  # For good results
    return energy_labels
end

function _truncate_energy_labels(
    energy_labels::AbstractVector{<:Real},
    config::Config;
    cutoff::Real=1e-12
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