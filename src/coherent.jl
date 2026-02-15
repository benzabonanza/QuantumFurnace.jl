#* COHERENT TERMS -----------------------------------------------------------------------------------------------------------
"""
    _precompute_coherent_total_B(
        jumps,
        hamiltonian,
        config,
        precomputed_data;
        trotter=nothing,
    ) -> Union{Nothing, Matrix{<:Complex}}

    Returns the total coherent operator B = sum_k B_k, already scaled by gamma_norm_factor.
    Returns nothing if config.with_coherent == false.
"""
function _precompute_coherent_total_B(
    jumps::AbstractVector{<:JumpOp},
    ham_or_trott::Union{HamHam, TrottTrott},
    config::AbstractConfig,
    precomputed_data;
    )

    config.with_coherent || return nothing

    if config.domain isa TimeDomain
        (; b_minus, b_plus, gamma_norm_factor) = precomputed_data
        B = B_time(jumps, ham_or_trott, b_minus, b_plus, config.t0, config.beta, config.sigma)

    elseif config.domain isa TrotterDomain
        (; b_minus, b_plus, gamma_norm_factor) = precomputed_data
        @assert ham_or_trott !== nothing
        trotter_jumps = transform_jumps_to_basis(jumps, ham_or_trott.eigvecs)
        B = B_trotter(trotter_jumps, ham_or_trott, b_minus, b_plus, config.beta, config.sigma)

    else
        # BohrDomain / EnergyDomain
        (; gamma_norm_factor) = precomputed_data
        B = B_bohr(ham_or_trott, jumps, config)
    end

    rmul!(B, gamma_norm_factor)
    return B
end

"""
    _precompute_coherent_unitary_terms(
        jumps::AbstractVector{<:JumpOp},
        hamiltonian::HamHam,
        config::AbstractThermalizeConfig,
        precomputed_data;
        trotter::Union{Nothing, TrottTrott}=nothing,
    ) -> Union{Nothing, Vector{Matrix{<:Complex}}}

    Precompute per-jump coherent unitaries for Kraus thermalization:
        U_k = exp(-1im * config.delta * B_k)

    Each B_k is constructed exactly as in the coherent-term definitions (domain-dependent),
    scaled by `gamma_norm_factor` (same convention as Liouvillian construction).
    Returns `nothing` if `config.with_coherent == false`.
"""
function _precompute_coherent_unitary_terms(
    jumps::AbstractVector{<:JumpOp},
    hamiltonian::HamHam,
    config::AbstractThermalizeConfig,
    precomputed_data;
    trotter::Union{Nothing, TrottTrott}=nothing,
    delta_scale::Real = 1.0  # for randomized channels
    )

    config.with_coherent || return nothing

    delta = delta_scale * config.delta
    CT = Complex{eltype(hamiltonian.eigvals)}
    U_terms = Vector{Matrix{CT}}(undef, length(jumps))

    if config.domain isa TimeDomain
        (; b_minus, b_plus, gamma_norm_factor) = precomputed_data
        @inbounds for (k, jump) in pairs(jumps)
            B = B_time(jump, hamiltonian, b_minus, b_plus, config.t0, config.beta, config.sigma)
            rmul!(B, gamma_norm_factor)
            U_terms[k] = exp(-1im * delta * Hermitian(B))
        end

    elseif config.domain isa TrotterDomain
        (; b_minus, b_plus, gamma_norm_factor) = precomputed_data
        @assert trotter !== nothing
        trotter_jumps = transform_jumps_to_basis(jumps, trotter.eigvecs)
        @inbounds for (k, jump) in pairs(trotter_jumps)
            B = B_trotter(jump, trotter, b_minus, b_plus, config.beta, config.sigma)
            rmul!(B, gamma_norm_factor)
            U_terms[k] = exp(-1im * delta * Hermitian(B))
        end

    else
        # BohrDomain / EnergyDomain
        (; gamma_norm_factor) = precomputed_data
        @inbounds for (k, jump) in pairs(jumps)
            B = B_bohr(hamiltonian, jump, config)
            rmul!(B, gamma_norm_factor)
            U_terms[k] = exp(-1im * delta * Hermitian(B))
        end
    end
    return U_terms
end


"""
    precompute_coherent_terms(
        jumps::AbstractVector{<:JumpOp},
        hamiltonian::HamHam,
        config::AbstractConfig,
        precomputed_data;
        trotter::Union{Nothing, TrottTrott}=nothing,
    ) -> Union{Nothing, Vector{Matrix{<:Complex}}}

    Precompute and cache the coherent B term for each `JumpOp`, already scaled by `gamma_norm_factor`.
    Returns `nothing` if `config.with_coherent == false`.
"""
function _precompute_coherent_terms(
    jumps::AbstractVector{<:JumpOp},
    hamiltonian::HamHam,
    config::AbstractConfig,
    precomputed_data;
    trotter::Union{Nothing, TrottTrott}=nothing,
    )

    config.with_coherent || return nothing

    CT = Complex{eltype(hamiltonian.eigvals)}
    coherent_terms = Vector{Matrix{CT}}(undef, length(jumps))

    if config.domain isa TimeDomain
        (; b_minus, b_plus, gamma_norm_factor) = precomputed_data
        @inbounds for (k, jump) in pairs(jumps)
            B = B_time(jump, hamiltonian, b_minus, b_plus, config.t0, config.beta, config.sigma)
            rmul!(B, gamma_norm_factor)
            coherent_terms[k] = B
        end

    elseif config.domain isa TrotterDomain
        (; b_minus, b_plus, gamma_norm_factor) = precomputed_data
        @assert trotter !== nothing
        trotter_jumps = transform_jumps_to_basis(jumps, trotter.eigvecs)
        @inbounds for (k, jump) in pairs(trotter_jumps)
            B = B_trotter(jump, trotter, b_minus, b_plus, config.beta, config.sigma)
            rmul!(B, gamma_norm_factor)
            coherent_terms[k] = B
        end

    else
        # BohrDomain / EnergyDomain: coherent term is the BohrDomain B operator.
        (; gamma_norm_factor) = precomputed_data
        @inbounds for (k, jump) in pairs(jumps)
            B = B_bohr(hamiltonian, jump, config)
            rmul!(B, gamma_norm_factor)
            coherent_terms[k] = B
        end
    end

    return coherent_terms
end

# (3.1) and Proposition III.1
# Has to be on a symmetric time domain, otherwise it can't be Hermitian.
function B_time(jump::JumpOp, hamiltonian::HamHam, b_minus, b_plus, t0, beta, sigma)

    d = size(hamiltonian.data, 1)
    CT = Complex{eltype(hamiltonian.eigvals)}
    eigvals = hamiltonian.eigvals
    jump_eig = jump.in_eigenbasis

    # Pre-allocated diagonal vector buffers (replace Diagonal wrappers)
    diag_u = Vector{CT}(undef, d)
    diag_u2 = Vector{CT}(undef, d)

    # Pre-allocated matrix scratch buffers
    b_plus_summand = zeros(CT, d, d)
    tmp = Matrix{CT}(undef, d, d)
    M = Matrix{CT}(undef, d, d)

    # Inner summand b_plus
    # Product: diag(u) * A' * diag(u2) * A * diag(u)
    for s in keys(b_plus)
        t_s = s * beta
        @. diag_u = exp(1im * eigvals * t_s)
        @. diag_u2 = exp(-2im * eigvals * t_s)

        # tmp = diag(u2) * A  (row-scale A by u2)
        @. tmp = diag_u2 * jump_eig
        # M = A' * tmp = A' * diag(u2) * A
        mul!(M, jump_eig', tmp)
        # b_plus_summand += b_s * diag(u) * M * diag(u) = b_s * (u .* M .* u')
        diag_u_row = transpose(diag_u)  # (1,d) row view for column-scaling
        b_plus_summand .+= b_plus[s] .* diag_u .* M .* diag_u_row
    end

    # Outer summand b_minus
    # B += b_t * diag(u)' * b_plus_summand * diag(u) = b_t * conj(u) .* bps .* u'
    B = zeros(CT, d, d)
    for t in keys(b_minus)
        @. diag_u = exp(1im * eigvals * (t / sigma))
        diag_u_row = transpose(diag_u)
        B .+= b_minus[t] .* conj.(diag_u) .* b_plus_summand .* diag_u_row
    end

    return B .* t0^2
end

function B_time(jumps::Vector{JumpOp}, hamiltonian::HamHam, b_minus, b_plus, t0, beta, sigma)

    d = size(hamiltonian.data, 1)
    CT = Complex{eltype(hamiltonian.eigvals)}
    eigvals = hamiltonian.eigvals

    # Pre-allocated diagonal vector buffers (replace Diagonal wrappers)
    diag_u = Vector{CT}(undef, d)
    diag_u2 = Vector{CT}(undef, d)

    # Pre-allocated matrix scratch buffers
    b_plus_summand = zeros(CT, d, d)
    tmp = Matrix{CT}(undef, d, d)
    M = Matrix{CT}(undef, d, d)

    # Inner summand b_plus (all jumps A^a)
    # Product per jump: diag(u) * A' * diag(u2) * A * diag(u)
    for s in keys(b_plus)
        t_s = s * beta
        @. diag_u = exp(1im * eigvals * t_s)
        @. diag_u2 = exp(-2im * eigvals * t_s)
        diag_u_row = transpose(diag_u)

        for jump_a in jumps
            jump_eig = jump_a.in_eigenbasis
            # tmp = diag(u2) * A  (row-scale A by u2)
            @. tmp = diag_u2 * jump_eig
            # M = A' * tmp = A' * diag(u2) * A
            mul!(M, jump_eig', tmp)
            # b_plus_summand += b_s * diag(u) * M * diag(u)
            b_plus_summand .+= b_plus[s] .* diag_u .* M .* diag_u_row
        end
    end

    # Outer summand b_minus
    # B += b_t * diag(u)' * b_plus_summand * diag(u)
    B = zeros(CT, d, d)
    for t in keys(b_minus)
        @. diag_u = exp(1im * eigvals * (t / sigma))
        diag_u_row = transpose(diag_u)
        B .+= b_minus[t] .* conj.(diag_u) .* b_plus_summand .* diag_u_row
    end

    return B .* t0^2
end

function B_trotter(jump::JumpOp, trotter::TrottTrott, b_minus, b_plus, beta, sigma)

    d = size(trotter.eigvecs, 1)
    CT = Complex{eltype(trotter.bohr_freqs)}

    # Use jump.in_eigenbasis directly -- callers ensure it contains
    # the Trotter-basis representation via transform_jumps_to_basis
    jump_eig = jump.in_eigenbasis

    # Pre-allocated diagonal vector buffers (replace Diagonal wrappers)
    diag_u = Vector{CT}(undef, d)
    diag_u2 = Vector{CT}(undef, d)

    # Pre-allocated matrix scratch buffers
    b_plus_summand = zeros(CT, d, d)
    tmp = Matrix{CT}(undef, d, d)
    M = Matrix{CT}(undef, d, d)

    # Inner summand b_plus
    # Product: diag(u) * A' * diag(u2) * A * diag(u)
    for (s, b_s) in b_plus
        num_t0_steps = Int(round(s * beta / trotter.t0))

        @. diag_u = trotter.eigvals_t0 ^ num_t0_steps
        @. diag_u2 = trotter.eigvals_t0 ^ (-2 * num_t0_steps)

        # tmp = diag(u2) * A  (row-scale A by u2)
        @. tmp = diag_u2 * jump_eig
        # M = A' * tmp = A' * diag(u2) * A
        mul!(M, jump_eig', tmp)
        # b_plus_summand += b_s * diag(u) * M * diag(u)
        diag_u_row = transpose(diag_u)
        b_plus_summand .+= b_s .* diag_u .* M .* diag_u_row
    end

    B = zeros(CT, d, d)
    for (t, b_t) in b_minus
        num_t0_steps = Int(round(t / (sigma * trotter.t0)))
        @. diag_u = trotter.eigvals_t0 ^ num_t0_steps

        # B += b_t * diag(u)' * b_plus_summand * diag(u)
        diag_u_row = transpose(diag_u)
        B .+= b_t .* conj.(diag_u) .* b_plus_summand .* diag_u_row
    end

    return B .* trotter.t0^2  # B in Trotter basis
end

function B_trotter(jumps::Vector{JumpOp}, trotter::TrottTrott, b_minus, b_plus, beta, sigma)

    d = size(trotter.eigvecs, 1)
    CT = Complex{eltype(trotter.bohr_freqs)}

    # Pre-allocated diagonal vector buffers (replace Diagonal wrappers)
    diag_u = Vector{CT}(undef, d)
    diag_u2 = Vector{CT}(undef, d)

    # Pre-allocated matrix scratch buffers
    b_plus_summand = zeros(CT, d, d)
    tmp = Matrix{CT}(undef, d, d)
    M = Matrix{CT}(undef, d, d)

    # Inner summand b_plus
    # Product per jump: diag(u) * A' * diag(u2) * A * diag(u)
    for (s, b_s) in b_plus
        num_t0_steps = Int(round(s * beta / trotter.t0))

        @. diag_u = trotter.eigvals_t0 ^ num_t0_steps
        @. diag_u2 = trotter.eigvals_t0 ^ (-2 * num_t0_steps)
        diag_u_row = transpose(diag_u)

        for jump_a in jumps
            # Use jump_a.in_eigenbasis directly -- callers ensure it contains
            # the Trotter-basis representation via transform_jumps_to_basis
            jump_a_eig = jump_a.in_eigenbasis
            # tmp = diag(u2) * A  (row-scale A by u2)
            @. tmp = diag_u2 * jump_a_eig
            # M = A' * tmp = A' * diag(u2) * A
            mul!(M, jump_a_eig', tmp)
            # b_plus_summand += b_s * diag(u) * M * diag(u)
            b_plus_summand .+= b_s .* diag_u .* M .* diag_u_row
        end
    end

    B = zeros(CT, d, d)
    for (t, b_t) in b_minus
        num_t0_steps = Int(round(t / (sigma * trotter.t0)))
        @. diag_u = trotter.eigvals_t0 ^ num_t0_steps

        # B += b_t * diag(u)' * b_plus_summand * diag(u)
        diag_u_row = transpose(diag_u)
        B .+= b_t .* conj.(diag_u) .* b_plus_summand .* diag_u_row
    end

    return B .* trotter.t0^2  # B in Trotter basis
end

#* B1 AND B2 ---------------------------------------------------------------------------------------------------------------
#TODO: Reintroduce sigmas here
function _compute_b_minus(t::Real, beta::Real, sigma::Real)  # 2pi sqrt(pi) * f_minus(t / sigma_E)
    f1(t) = 1 / cosh(2 * pi * t / (beta * sigma))
    f2(t) = sin(-t * beta * sigma) * exp(-2 * t^2)
    return 2 * sqrt(pi) * exp(beta^2 * sigma^2 / 8) * _convolute(f1, f2, t) / (beta * sigma)
end

function _compute_b_plus(t::Real, beta::Real, w_gamma::Real, sigma_gamma::Real)  # f_plus(t * beta) / (2pi sqrt(pi))
    return beta * sigma_gamma * exp(- 2 * beta * w_gamma * (2 * t^2 + im * t)) / sqrt(pi^3)
end

function _compute_b_plus_metro(t::Real, beta::Real, sigma::Real, eta::Real)
    if abs(t) < 1e-12  # Handle t=0
        return complex(1 / (2 * sqrt(2) * pi^2))
    elseif abs(t) <= eta
        numerator = exp(- sigma^2 * beta^2 * (2 * t^2 + 1im * t)) + 1im * (2 * t + 1im)
    else
        numerator = exp(- sigma^2 * beta^2 * (2 * t^2 + 1im * t))
    end
    denominator = t * (2 * t + 1im)
    return (1 / (2 * sqrt(2) * pi^2)) * numerator / denominator
end

function _compute_b_plus_smooth(t::Real, beta::Real, sigma::Real, a::Real, b::Real)
    b_vals = exp(- a * b / 2) * exp(- sigma^2 * beta^2 * t * (2 * t + 1im) * (1 + b)) / (4 * t^2 + a + 2im * t)
    return sqrt(4 * a + 1) * b_vals / (sqrt(2) * pi^2)
end

function _compute_truncated_func(target_func::Function, time_labels::AbstractVector{<:Real}, fixed_args...; atol::Real = 1e-12)
    f_vals = Vector{ComplexF64}(target_func.(time_labels, fixed_args...))
    indices_to_keep = _get_truncated_indices(f_vals; atol=atol)
    return Dict(zip(time_labels[indices_to_keep], f_vals[indices_to_keep]))
end

#* TOOLS --------------------------------------------------------------------------------------------------------------------
function _get_truncated_indices(fvals::AbstractVector{<:Real}; atol::Real = 1e-12)
   """Find elements in `fvals` that are larger than `atol`"""
    return findall(abs.(fvals) .>= atol)
end

function _get_truncated_indices(fvals::AbstractVector{<:Complex}; atol::Real = 1e-12)
    """Find elements in `fvals` that are larger than `atol`"""
     return findall(abs.(fvals) .>= atol)
 end

function _convolute(f::Function, g::Function, t::Real; atol=1e-12, rtol=1e-12)
    integrand(s) = f(s) * g(t - s)
    result, _ = quadgk(integrand, -Inf, Inf; atol=atol, rtol=rtol)
    return result
end
