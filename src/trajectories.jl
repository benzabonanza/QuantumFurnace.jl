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