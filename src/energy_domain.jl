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
    if config.s == 0 && config.a == 0
        return exp(-config.beta * max(w + config.beta * config.sigma^2 / 2, 0.0))
    else
        # Smooth Metropolis (thesis eq:smooth-metro). Handles any a, including a == 0
        # (the thesis-main case). At a == 0 reduces to γ_M^{(0)} × (1/2)[erfc(z_-) + e^{β|ω̃|}erfc(z_+)].
        # The (s=0, a>0) "a-regularized, no smoothing" case is rejected by
        # validate_config! and therefore unreachable here.
        u_min = sqrt(config.beta * config.sigma^2 * config.s / 2)
        transition_b0 = exp(-2 * sqrtA * sqrtB - config.beta * w / 2 - config.beta^2 * config.sigma^2 / 4)
        return transition_b0 * (erfc(sqrtA * u_min - sqrtB / u_min) + exp(4 * sqrtA * sqrtB) * erfc(sqrtA * u_min + sqrtB / u_min)) / 2
    end
end

function pick_transition(config::Config{<:Any, <:Any, GNS}, w::Real)
    if !(config.with_linear_combination)
        w_gamma = config.gaussian_parameters[1]
        sigma_gamma = config.gaussian_parameters[2]
        return exp(-(w + w_gamma)^2 / (2 * sigma_gamma^2))
    end
    sqrtA = sqrt(config.beta / 4) * sqrt(4 * config.a + 1)
    sqrtB = sqrt(config.beta / 4) * abs(w)
    if config.s == 0 && config.a == 0
        return exp(-config.beta * max(w, 0.0))
    else
        # Smooth Metropolis (un-shifted GNS form). Handles any a, including a == 0.
        # The (s=0, a>0) "a-regularized, no smoothing" case is rejected by
        # validate_config! and therefore unreachable here.
        u_min = sqrt(config.beta * config.sigma^2 * config.s / 2)
        transition_b0 = exp(-2 * sqrtA * sqrtB - config.beta * w / 2)
        return transition_b0 * (erfc(sqrtA * u_min - sqrtB / u_min) + exp(4 * sqrtA * sqrtB) * erfc(sqrtA * u_min + sqrtB / u_min)) / 2
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
    # The (s=0, a>0) "a-regularized, no smoothing" case is rejected by
    # validate_config! and therefore unreachable here.
    if (config.s == 0 && config.a == 0)  # Kinky Metropolis (γ_M^{(0)})
        return w -> exp(-config.beta * max(w + config.beta * config.sigma^2 / 2, 0.0))
    else  # config.s != 0 — smooth Metropolis (thesis eq:smooth-metro), any a (incl. a == 0)
        return w -> begin
            sqrtB = sqrt(config.beta / 4) * abs(w + config.beta * config.sigma^2 / 2)
            u_min = sqrt(config.beta * config.sigma^2 * config.s / 2)
            transition_b0 = exp((- 2 * sqrtA * sqrtB - config.beta * w / 2 - config.beta^2 * config.sigma^2 / 4))
            return (transition_b0 * (erfc(sqrtA * u_min - sqrtB / u_min)
                + exp(4 * sqrtA * sqrtB) * erfc(sqrtA * u_min + sqrtB / u_min)) / 2)
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
    # The (s=0, a>0) "a-regularized, no smoothing" case is rejected by
    # validate_config! and therefore unreachable here.
    if (config.s == 0 && config.a == 0)
        # Kinky Metropolis — UN-SHIFTED.
        return w -> exp(-config.beta * max(w, 0.0))
    else  # config.s != 0 — smooth Metropolis (un-shifted GNS form), any a (incl. a == 0)
        return w -> begin
            sqrtB = sqrt(config.beta / 4) * abs(w)
            u_min = sqrt(config.beta * config.sigma^2 * config.s / 2)
            transition_b0 = exp((-2 * sqrtA * sqrtB - config.beta * w / 2))
            return (transition_b0 * (erfc(sqrtA * u_min - sqrtB / u_min)
                + exp(4 * sqrtA * sqrtB) * erfc(sqrtA * u_min + sqrtB / u_min)) / 2)
        end
    end
end


"""
    pick_gamma_sup(config::Config) -> Real

Closed-form continuum supremum `‖γ‖_∞` for the rate function `γ(ω)` selected by
`config`. Returns the analytical sup independent of any energy grid — this is
the grid-independent replacement for the prior `1.0 / maximum(transition.(energy_labels))`
used to populate `gamma_norm_factor`.

For every standard rate family currently supported (KMS / GNS Gaussian, kinky
Metropolis at `s = a = 0`, smooth Metropolis at any `s > 0, a ≥ 0`), the
continuum supremum is exactly `1.0`:
- Gaussian: `γ(ω) = exp(-(ω+ω_γ)²/(2σ_γ²))`, sup = 1 at `ω = -ω_γ`.
- Kinky Metropolis (KMS): `γ(ω) = exp(-β·max(ω + βσ²/2, 0))`, sup = 1 on `ω ≤ -βσ²/2`.
- Kinky Metropolis (GNS, un-shifted): `γ(ω) = exp(-β·max(ω, 0))`, sup = 1 on `ω ≤ 0`.
- Smooth Metropolis (any `a ≥ 0`, `s > 0`): sup = 1, attained as `ω → -∞`.

The `(s = 0, a > 0)` "a-regularized, no smoothing" case is rejected by
`validate_config!` (it is outside the thesis's `(a, s)` taxonomy: smoothing
requires `s > 0`; kinky requires `s = a = 0`).

The reciprocal of this value is what populates `gamma_norm_factor` in
`_precompute_data` (so `gamma_norm_factor = 1 / pick_gamma_sup(config) = 1.0`
for all current families).
"""
pick_gamma_sup(config::Config{<:Any, <:Any, KMS}) = 1.0
pick_gamma_sup(config::Config{<:Any, <:Any, GNS}) = 1.0


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
end
#* --------------------------------------------------------------------------------------------------------------------------