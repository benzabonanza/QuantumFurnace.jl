function precompute_labels(::Union{BohrDomain, EnergyDomain}, config::AbstractConfig)
    energy_labels = create_energy_labels(config.num_energy_bits, config.w0)
    truncated_energy_labels = truncate_energy_labels(energy_labels, config)
    return (truncated_energy_labels,)  # Energy labels
end

function precompute_labels(::Union{TimeDomain, TrotterDomain}, config::AbstractConfig)
    energy_labels = create_energy_labels(config.num_energy_bits, config.w0)
    truncated_energy_labels = truncate_energy_labels(energy_labels, config)
    time_labels = energy_labels .* (config.t0 / config.w0)
    return (truncated_energy_labels, time_labels) # Energy and time labels
end  

function precompute_data(::BohrDomain,
    config::AbstractConfig,
    ham_or_trott::Union{HamHam, TrottTrott}
)

    alpha = pick_alpha(config)
    # Was the only way to bring in the normalizing factor 1 / ||γ||_∞
    energy_labels, = precompute_labels(config.domain, config)
    transition = pick_transition(config)
    gamma_norm_factor =  1.0 / maximum(transition.(energy_labels))
    return (
        alpha = alpha,
        gamma_norm_factor = gamma_norm_factor
    )
end

function precompute_data(::EnergyDomain, 
    config::AbstractConfig, 
    ham_or_trott::Union{HamHam, TrottTrott}
)
    energy_labels, = precompute_labels(config.domain, config)
    transition = pick_transition(config)
    gamma_norm_factor =  1.0 / maximum(transition.(energy_labels))
    
    return (
        transition = transition,
        gamma_norm_factor = gamma_norm_factor,
        energy_labels = energy_labels
    )
end

function precompute_data(::Union{TimeDomain, TrotterDomain},
    config::AbstractConfig,
    ham_or_trott::Union{HamHam, TrottTrott}
)
    energy_labels, time_labels = precompute_labels(config.domain, config)
    oft_time_labels = truncate_time_labels_for_oft(time_labels, config.sigma)

    transition = pick_transition(config)
    gamma_norm_factor =  1.0 / maximum(transition.(energy_labels))

    # Coherent term B
    b_minus, b_plus = if config.with_coherent
        _b_minus = compute_truncated_func(compute_b_minus, time_labels, config.beta, config.sigma)
        chosen_b_plus, b_plus_args = select_b_plus_calculator(config)
        _b_plus = compute_truncated_func(chosen_b_plus, time_labels, b_plus_args...)
        (_b_minus, _b_plus)
    else
        (nothing, nothing)
    end

    # OFT NUFFT prefactors
    if config.domain isa TimeDomain
        oft_nufft_prefactors = prepare_oft_nufft_prefactors(
            ham_or_trott.bohr_freqs,
            oft_time_labels,
            energy_labels,
            config.sigma;
            eps=1e-12,
            nthreads=1,
            use_shared_array=(nprocs() > 1),
        )
    elseif config.domain isa TrotterDomain
        oft_nufft_prefactors = prepare_oft_nufft_prefactors(
            ham_or_trott.bohr_freqs,
            oft_time_labels,
            energy_labels,
            config.sigma;
            eps=1e-12,
            nthreads=1,
            use_shared_array=(nprocs() > 1),
        )
    end
    
    return (
        transition = transition,
        gamma_norm_factor = gamma_norm_factor,
        energy_labels = energy_labels,
        oft_nufft_prefactors = oft_nufft_prefactors,
        b_minus = b_minus,
        b_plus = b_plus,
    )
end

function select_b_plus_calculator(config::Union{LiouvConfig, ThermalizeConfig})
    if !config.with_linear_combination
        # Gaussian
        return (compute_b_plus, (config.beta, config.gaussian_parameters[1], config.gaussian_parameters[2]))
    else
        if config.a != 0.0
            # Improved Metro / Glauber
            return (compute_b_plus_smooth, (config.beta, config.sigma, config.a, config.b))
        else
            # Metro
            return (compute_b_plus_metro, (config.beta, config.sigma, config.eta,))
        end
    end
end