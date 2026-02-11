# struct KrausFramework{T}
#     kraus_H_eff::Matrix{T}                   # non-Hermitian no jump evolution
#     kraus_jumps::Vector{Matrix{T}}      # jump Kraus
#     R::Matrix{T}                    # sum(M^\dagger M) (without delta)
#     psi_temp::Vector{T}             # buffer for less allocations
#     delta::Float64                  # delta steps for trajectory
# end

# function krausframework(
#     H::AbstractMatrix{T}, 
#     kraus_jumps::Vector{Tuple{Float64, <:AbstractMatrix{T}}}, 
#     R::AbstractMatrix{T},
#     delta::Float64) where T

#     dim = size(H, 1)
    
#     kraus_jumps = [Matrix{T}(sqrt(delta) * rate * op) for (rate, op) in kraus_jumps]

#     # Construct H_eff and the 0th Kraus operator
#     H_eff = 1im * H + 0.5 * R
#     kraus_H_eff = Matrix{T}(I, dim, dim) - delta * H_eff

#     # Buffer for statevector
#     psi_temp = zeros(T, dim)

#     return KrausFramework(kraus_H_eff, kraus_jumps, R, psi_temp, delta)
# end

"""
Run-time cache for trajectory simulations
"""
struct TrajectoryWorkspace{T}
    jump_oft::Matrix{T}   # buffer for A_ω (or its dagger-partner is handled by swapping mul order)
    psi_tmp::Vector{T}    # generic tmp for matvec results (K0ψ, Aωψ, Uresψ, etc.)
    Rpsi::Vector{T}       # cache Rψ so K0ψ can be formed as ψ - α(Rψ) without a second matvec
end

struct TrajectoryFramework{T,C,H,PD,D<:AbstractDomain}
    domain::D
    jumps::Vector{JumpOp}
    ham_or_trott::H
    config::C
    precomputed_data::PD

    # Coherent term (optional)
    B::Union{Nothing, Matrix{T}}
    U_B::Union{Nothing, Matrix{T}}

    # Dissipative Kraus skeleton
    R::Matrix{T}
    K0::Matrix{T}
    U_residual::Matrix{T}   # Cholesky factor U with S ≈ U' U (used as residual Kraus)

    # Step parameters
    delta::Float64
    alpha::Float64           # α = 1 - sqrt(1-δ)

    # Runtime workspace (mutated during stepping)
    ws::TrajectoryWorkspace{T}
end

function build_trajectoryframework(
    jumps::AbstractVector{<:JumpOp},
    ham_or_trott::Union{HamHam, TrottTrott},
    config::AbstractThermalizeConfig,
    precomputed_data,
    scratch::KrausScratch{ComplexF64},
    delta::Float64) where T

    dim = size(H, 1)
    
    B_total = precompute_coherent_total_B(jumps, hamiltonian, config, precomputed_data; trotter=trotter)
    B_total .= 0.5 .* (B_total .+ B_total')
    U_B = exp(-1im * delta * Hermitian(B_total))

    R = precompute_R(config.domain, jumps, ham_or_trott, config, precomputed_data, scratch)
    R = copy(scratch.R)

    alpha = 1 - sqrt(1 - delta)
    K0 = copy(R)
    K0 .*= (-alpha)
    @inbounds for i in 1:dim
        K0[i,i] += 1
    end

    mul!(scratch.tmp1, R, R)  # tmp1 := R²
    s1 = 2 * alpha - delta
    s2 = alpha * alpha
    @. scratch.tmp2 = s1 * R - s2 * scratch.tmp1   # tmp2 := S (up to roundoff)

    # Numerical symmetrization and tiny diagonal shift (matches thermalization code logic)
    scratch.tmp2 .= 0.5 .* (scratch.tmp2 .+ scratch.tmp2')
    eps_shift = 10 * eps(Float64)
    @inbounds for i in 1:dim
        scratch.tmp2[i,i] += eps_shift
    end

    S = copy(scratch.tmp2)  # store S explicitly (post-shift)

    # Factor for applying the residual Kraus to state vectors (K_res ≈ U where S ≈ U' U)
    cholesky_S = cholesky!(Hermitian(scratch.tmp2), check=false)
    U_residual = Matrix{ComplexF64}(cholesky_S.U)  # detach from scratch.tmp2

    ws = TrajectoryWorkspace(ComplexF64, dim)

    return TrajectoryFramework(
        config.domain,
        collect(jumps),
        ham_or_trott,
        config,
        precomputed_data,
        B_total,
        U_B,
        R,
        K0,
        U_residual,
        delta,
        alpha,
        ws,
    )
    return
end

"""
    precompute_R(domain, jumps, ham_or_trott, config, precomputed_data, scratch) -> Matrix{ComplexF64}

    Compute
        R = ∑_{k>0} L_k† L_k
    in the same basis as `jump.in_eigenbasis` (Hamiltonian eigenbasis for `HamHam`, Trotter basis for `TrottTrott`).

    Conventions are matched to `jump_contribution!(domain, ::AbstractThermalizeConfig, ...)`:
    - the weights are `rate2(ω) = base_prefactor * transition(ω)` (no extra `δ` factor),
    - for Hermitian jumps we iterate half-grid and add the mirrored negative-frequency partner explicitly.

    This returns `scratch.R` (Hermitianized).
"""
function precompute_R(
    ::EnergyDomain,
    jumps::AbstractVector{<:JumpOp},
    hamiltonian::HamHam,
    config::AbstractThermalizeConfig,
    precomputed_data,
    scratch::KrausScratch{ComplexF64},
)
    dim = size(hamiltonian.data, 1)
    (; transition, gamma_norm_factor, energy_labels) = precomputed_data


    base_prefactor = config.w0 / (config.sigma * sqrt(2 * pi)) * gamma_norm_factor

    fill!(scratch.R, 0)

    @inbounds for jump in jumps
        if jump.hermitian
            # Half-grid (w_raw <= 0) and mirror partner at -w using Aω† as the Lindblad operator.
            for w_raw in energy_labels
                w_raw > 1e-12 && continue
                w = abs(w_raw)

                # Aω := A ∘ exp(-(w-ν)^2/(4σ^2))   (elementwise in eigenbasis)
                oft!(scratch.jump_oft, jump, w, hamiltonian, config.sigma)

                # Positive-frequency contribution: rate2(w) * (Aω† Aω)
                rate2_pos = base_prefactor * transition(w)
                mul!(scratch.LdagL, scratch.jump_oft', scratch.jump_oft)
                @. scratch.R += rate2_pos * scratch.LdagL

                if w > 1e-12
                    # Negative-frequency partner uses Lindblad op = Aω†:
                    # contribution is rate2(-w) * (Aω Aω†)
                    rate2_neg = base_prefactor * transition(-w)
                    mul!(scratch.LdagL, scratch.jump_oft, scratch.jump_oft')
                    @. scratch.R += rate2_neg * scratch.LdagL
                end
            end
        else
            # Non-Hermitian jump: full grid, no mirroring shortcut.
            for w in energy_labels
                oft!(scratch.jump_oft, jump, w, hamiltonian, config.sigma)

                rate2 = base_prefactor * transition(w)
                mul!(scratch.LdagL, scratch.jump_oft', scratch.jump_oft)
                @. scratch.R += rate2 * scratch.LdagL
            end
        end
    end

    # Numerical Hermitianization (R should be Hermitian PSD by construction).
    scratch.R .= 0.5 .* (scratch.R .+ scratch.R')
    return scratch.R
end


function precompute_R(
    ::Union{TimeDomain, TrotterDomain},
    jumps::AbstractVector{<:JumpOp},
    ham_or_trott::Union{HamHam, TrottTrott},
    config::AbstractThermalizeConfig,
    precomputed_data,
    scratch::KrausScratch{ComplexF64},
)
    dim = size(jumps[1].in_eigenbasis, 1)
    (; transition, gamma_norm_factor, energy_labels, oft_nufft_prefactors) = precomputed_data

    # Same weight as in jump_contribution!(::Union{TimeDomain,TrotterDomain}, ...), but without δ.
    base_prefactor = config.w0 * config.t0^2 * (config.sigma * sqrt(2 / pi)) / (2 * pi) * gamma_norm_factor

    fill!(scratch.R, 0)

    @inbounds for jump in jumps
        if jump.hermitian
            # Half-grid and explicit mirror, allocation-free (no filter/abs vector creation).
            for w_raw in energy_labels
                w_raw > 1e-12 && continue
                w = abs(w_raw)

                nufft_prefactor_matrix = prefactor_view(oft_nufft_prefactors, w)

                # Aω := A ∘ prefactor_matrix(ω)  (elementwise)
                @. scratch.jump_oft = jump.in_eigenbasis * nufft_prefactor_matrix

                rate2_pos = base_prefactor * transition(w)
                mul!(scratch.LdagL, scratch.jump_oft', scratch.jump_oft)
                @. scratch.R += rate2_pos * scratch.LdagL

                if w > 1e-12
                    rate2_neg = base_prefactor * transition(-w)
                    mul!(scratch.LdagL, scratch.jump_oft, scratch.jump_oft')
                    @. scratch.R += rate2_neg * scratch.LdagL
                end
            end
        else
            for w in energy_labels
                nufft_prefactor_matrix = prefactor_view(oft_nufft_prefactors, w)
                @. scratch.jump_oft = jump.in_eigenbasis * nufft_prefactor_matrix

                rate2 = base_prefactor * transition(w)
                mul!(scratch.LdagL, scratch.jump_oft', scratch.jump_oft)
                @. scratch.R += rate2 * scratch.LdagL
            end
        end
    end

    scratch.R .= 0.5 .* (scratch.R .+ scratch.R')
    return scratch.R
end


function precompute_kraus_jump!(::EnergyDomain,
    fw::TrajectoryScratch,
    jump::JumpOp,
    hamiltonian::HamHam,
    config::AbstractThermalizeConfig,
    precomputed_data,
    )
end


function evolve_along_trajectory(psi0::Vector{ComplexF64}, fw::KrausFramework, total_time::Float64)::Vector{ComplexF64}

    num_steps = round(Int, total_time / fw.delta)  # Number of trajectory steps

    psi = copy(psi0)
    for _ in 1:num_steps
        step_along_the_trajectory!(psi, fw)
    end

    return psi
end

function evolve_and_measure_along_trajectory(
    psi0::Vector{ComplexF64}, 
    fw::KrausFramework, 
    T::Float64, 
    observables::Vector{Matrix{ComplexF64}};
    save_every::Int = 1
    )

    num_steps = round(Int, T / fw.delta)  # Number of trajectory steps
    num_saves = div(num_steps, save_every) + 1
    num_obs = length(observables)
    data = zeros(Float64, num_obs, num_saves)
    times = zeros(Float64, num_saves)

    psi = copy(psi0)
    # Measure initial state
    save_index = 1
    times[save_index] = 0.0
    measure!(view(data, :, save_index), psi, observables)

    for step in 1:num_steps
        step_along_the_trajectory!(psi, fw)

        if step % save_every == 0
            save_index += 1
            times[save_index] = step * fw.delta
            measure!(view(data, :, save_index), psi, observables)
        end
    end

    return (psi = psi, times = times, measurements = data)
end

function step_along_the_trajectory!(psi::Vector{ComplexF64}, fw::KrausFramework)

    # No jump
    mul!(fw.psi_temp, fw.M0, psi)
    prob_no_jump = norm(fw.psi_temp)^2

    # Total jump probability
    p_jump_total = fw.delta * real(dot(psi, fw.R, psi))

    total_weight = prob_no_jump + p_jump_total

    r = rand() * total_weight

    if r < prob_no_jump
        # No jump
        copyto!(psi, fw.psi_temp)

        # Force normalize
        rmul!(psi, 1.0 / sqrt(prob_no_jump))
    else  # Jump.
        # Iterate through jumps and their probabilites till we find the winner
        target_cummulative = r - prob_no_jump
        current_cummulative = 0.0
        
        for k in 1:length(fw.M_jumps)
            mul!(fw.psi_temp, fw.M_jumps[k], psi)
            prob_jump_k = norm(fw.psi_temp)^2
            # println("Prob of jump k")
            # println(prob_jump_k)

            current_cummulative += prob_jump_k

            if current_cummulative >= target_cummulative
                copyto!(psi, fw.psi_temp)

                # Force normalize
                rmul!(psi, 1.0 / sqrt(prob_jump_k))
                return
            end
        end

        # If somehow we haven't picked any jumps, then use the last one.
        copyto!(psi, fw.psi_temp)
        normalize!(psi)
    end
end

function measure!(measured_values, state, observable_list)
    for i in eachindex(observable_list)
        measured_values[i] = real(dot(state, observable_list[i], state))
    end
end

function construct_gksl_lindbladian(
    H::AbstractMatrix{ComplexF64}, 
    kraus_jumps::Vector{Matrix{ComplexF64}}
    )
    dim = size(H, 1)

    lindblad = zeros(ComplexF64, dim^2, dim^2)

    # Add coherent part
    vectorize_liouvillian_coherent!(lindblad, H)

    # Add dissipative part
    for kraus_jump in kraus_jumps
        vectorize_liouv_diss_and_add!(lindblad, kraus_jump, 1.0)
    end

    return lindblad
end