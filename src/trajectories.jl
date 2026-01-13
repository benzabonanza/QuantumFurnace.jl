#TODO: Test this code!
struct KrausFramework{T}
    M0::Matrix{T}                   # non-Hermitian no jump evolution
    M_jumps::Vector{Matrix{T}}      # jump Kraus
    R::Matrix{T}                    # sum(M^\dagger M) (without delta)
    psi_temp::Vector{T}             # buffer for less allocations
    delta::Float64                  # delta steps for trajectory
end

function build_krausframework(H::AbstractMatrix{T}, 
                              kraus_jumps::Vector{<:AbstractMatrix{T}}, 
                              R::AbstractMatrix{T},
                              delta::Float64) where T
    dim = size(H, 1)
    
    M_jumps = [Matrix{T}(sqrt(delta) * M) for M in kraus_jumps]

    # Construct H_eff and the 0th Kraus operator
    H_eff = 1im * H + 0.5 * R
    M0 = Matrix{T}(I, dim, dim) - delta * H_eff

    # Buffer for statevector
    psi_temp = zeros(T, dim)

    return KrausFramework(M0, M_jumps, R, psi_temp, delta)
end


function evolve_along_trajectory(psi0::Vector{ComplexF64}, fw::KrausFramework, T::Float64)

    num_steps = round(Int, T / fw.delta)  # Number of trajectory steps

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
    p_jump_total = delta * real(dot(psi, fw.R, psi))

    total_weight = prob_no_jump + p_jump_total

    r = rand() * total_weight

    if r < prob_no_jump
        # No jump
        copyto!(psi, fw.psi_temp)

        # Force normalize
        rmul!(psi, 1.0 / sqrt(prob_no_jump))
    else
        # Jump, but which jump?
        # Iterate through jumps and their probabilites till we find the winner
        target_cummulative = r - prob_no_jump
        current_cummulative = 0.0
        
        for k in 1:length(fw.M_jumps)
            mul!(fw.psi_temp, fw.M_jumps[k], psi)
            prob_jump_k = norm(fw.psi_temp)^2

            current_cummulative += prob_jump_k

            if current_cummulative >= target_cummulative
                copyto!(psi,  fw.psi_temp)

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

    for kraus_jump in kraus_jumps
        vectorize_liouv_diss_and_add!(lindblad, kraus_jump, 1.0)
    end

    return lindblad
end