function _precompute_labels(config::Config{<:Any, D}) where {D<:Union{BohrDomain, EnergyDomain}}
    energy_labels = _create_energy_labels(config.num_energy_bits, config.w0)
    truncated_energy_labels = _truncate_energy_labels(energy_labels, config)
    return (truncated_energy_labels,)  # Energy labels
end

function _precompute_labels(config::Config{<:Any, D}) where {D<:Union{TimeDomain, TrotterDomain}}
    energy_labels = _create_energy_labels(config.num_energy_bits, config.w0)
    truncated_energy_labels = _truncate_energy_labels(energy_labels, config)
    time_labels = energy_labels .* (config.t0 / config.w0)
    return (truncated_energy_labels, time_labels) # Energy and time labels
end  

function _precompute_data(
    config::Config{Lindbladian, BohrDomain},
    ham_or_trott::Union{HamHam, TrottTrott}
)

    alpha = _pick_alpha(config)
    # Was the only way to bring in the normalizing factor 1 / ||γ||_∞
    energy_labels, = _precompute_labels(config)
    transition = pick_transition(config)
    gamma_norm_factor =  1.0 / maximum(transition.(energy_labels))
    return (
        alpha = alpha,
        gamma_norm_factor = gamma_norm_factor
    )
end

function _precompute_data(
    config::Config{Thermalize, BohrDomain},
    hamiltonian::HamHam
)
    alpha = _pick_alpha(config)

    energy_labels, = _precompute_labels(config)
    transition = pick_transition(config)
    gamma_norm_factor = 1.0 / maximum(transition.(energy_labels))
    
    # Cache the Bohr buckets as plain Int index pairs to avoid CartesianIndex overhead
    # and avoid rebuilding any per-frequency index lists inside jump_contribution!.
    bohr_keys = collect(keys(hamiltonian.bohr_dict))
    bohr_is = Vector{Vector{Int}}(undef, length(bohr_keys))
    bohr_js = Vector{Vector{Int}}(undef, length(bohr_keys))
    @inbounds for (k, nu) in pairs(bohr_keys)
        idxs = hamiltonian.bohr_dict[nu]
        is = Vector{Int}(undef, length(idxs))
        js = Vector{Int}(undef, length(idxs))
        @inbounds for t in eachindex(idxs)
            is[t] = idxs[t][1]
            js[t] = idxs[t][2]
        end
        bohr_is[k] = is
        bohr_js[k] = js
    end

    return (
        alpha = alpha,
        gamma_norm_factor = gamma_norm_factor,
        bohr_keys = bohr_keys,
        bohr_is = bohr_is,
        bohr_js = bohr_js,
    )

end

function _precompute_data(
    config::Config{<:Any, EnergyDomain},
    ham_or_trott::Union{HamHam, TrottTrott}
)
    energy_labels, = _precompute_labels(config)
    transition = pick_transition(config)
    gamma_norm_factor =  1.0 / maximum(transition.(energy_labels))
    
    return (
        transition = transition,
        gamma_norm_factor = gamma_norm_factor,
        energy_labels = energy_labels
    )
end

function _precompute_data(
    config::Config{<:Any, D},
    ham_or_trott::Union{HamHam, TrottTrott}
) where {D<:Union{TimeDomain, TrotterDomain}}
    energy_labels, time_labels = _precompute_labels(config)
    oft_time_labels = _truncate_time_labels_for_oft(time_labels, config.sigma)

    transition = pick_transition(config)
    gamma_norm_factor =  1.0 / maximum(transition.(energy_labels))

    # Coherent term B
    b_minus, b_plus = if with_coherent(config.construction)
        _b_minus = _compute_truncated_func(_compute_b_minus, time_labels, config.beta, config.sigma)
        chosen_b_plus, b_plus_args = _select_b_plus_calculator(config)
        _b_plus = _compute_truncated_func(chosen_b_plus, time_labels, b_plus_args...)
        (_b_minus, _b_plus)
    else
        (nothing, nothing)
    end

    # OFT NUFFT prefactors (same call for both TimeDomain and TrotterDomain)
    oft_nufft_prefactors = _prepare_oft_nufft_prefactors(
        ham_or_trott.bohr_freqs,
        oft_time_labels,
        energy_labels,
        config.sigma;
        eps=1e-12,
        nthreads=1,
        use_shared_array=(nprocs() > 1),
    )
    
    return (
        transition = transition,
        gamma_norm_factor = gamma_norm_factor,
        energy_labels = energy_labels,
        oft_nufft_prefactors = oft_nufft_prefactors,
        b_minus = b_minus,
        b_plus = b_plus,
    )
end

function _select_b_plus_calculator(config::Config{<:Any, <:Any, KMS})
    if !config.with_linear_combination
        # Gaussian
        return (_compute_b_plus, (config.beta, config.gaussian_parameters[1], config.gaussian_parameters[2]))
    else
        if config.a != 0.0
            # Improved Metro / Glauber
            return (_compute_b_plus_smooth, (config.beta, config.sigma, config.a, config.b))
        else
            # Metro
            return (_compute_b_plus_metro, (config.beta, config.sigma, config.eta,))
        end
    end
end