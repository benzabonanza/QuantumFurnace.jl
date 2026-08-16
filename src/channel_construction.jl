# Retained construction of the per-jump rate operator R.
#
# This code is shared by the full-density-matrix channel and the matrix-free
# channel workspace.

"""
    _precompute_R(jumps, ham_or_trott, config, precomputed_data, scratch)

Accumulate the channel rate operator in the active spectral basis.

Hermitian jumps use the positive half-grid plus the explicit negative-frequency
partner. Weights exclude the channel step size.

# Returns
Hermitianized `scratch.R`.
"""
function _precompute_R(
    jumps::AbstractVector{<:JumpOp},
    hamiltonian::HamHam,
    config::Config{Thermalize, EnergyDomain},
    precomputed_data,
    scratch::ThermalizeScratch{<:Complex},
    )
    dim = size(hamiltonian.data, 1)
    (; transition, gamma_norm_factor, energy_labels) = precomputed_data

    n_labels = length(energy_labels)
    if Threads.nthreads() > 1 && n_labels >= OMEGA_THREAD_THRESHOLD
        return _precompute_R_threaded_energy!(jumps, hamiltonian, config, precomputed_data, scratch)
    end

    base_prefactor = precomputed_data.oft_domain_prefactor * gamma_norm_factor
    inv_4sigma2 = 1.0 / (4 * config.sigma^2)

    fill!(scratch.R, 0)

    @inbounds for jump in jumps
        if jump.hermitian
            # Half-grid (w_raw <= 0) and mirror partner at -w using Aw^dagger as the Lindblad operator.
            for w_raw in energy_labels
                w_raw > 1e-12 && continue
                w = abs(w_raw)

                # Aw := A .* exp(-(w-nu)^2/(4sigma^2))   (elementwise in eigenbasis)
                oft!(scratch.jump_oft, jump.in_eigenbasis, hamiltonian.bohr_freqs, w, inv_4sigma2)

                # Positive-frequency contribution: rate2(w) * (Aw^dagger Aw)
                rate2_pos = base_prefactor * transition(w)
                mul!(scratch.LdagL, scratch.jump_oft', scratch.jump_oft)
                @. scratch.R += rate2_pos * scratch.LdagL

                if w > 1e-12
                    # Negative-frequency partner uses Lindblad op = Aw^dagger:
                    # contribution is rate2(-w) * (Aw Aw^dagger)
                    rate2_neg = base_prefactor * transition(-w)
                    mul!(scratch.LdagL, scratch.jump_oft, scratch.jump_oft')
                    @. scratch.R += rate2_neg * scratch.LdagL
                end
            end
        else
            # Non-Hermitian jump: full grid, no mirroring shortcut.
            for w in energy_labels
                oft!(scratch.jump_oft, jump.in_eigenbasis, hamiltonian.bohr_freqs, w, inv_4sigma2)

                rate2 = base_prefactor * transition(w)
                mul!(scratch.LdagL, scratch.jump_oft', scratch.jump_oft)
                @. scratch.R += rate2 * scratch.LdagL
            end
        end
    end

    # Numerical Hermitianization (R should be Hermitian PSD by construction).
    hermitianize!(scratch.R)
    return scratch.R
end


function _precompute_R(
    jumps::AbstractVector{<:JumpOp},
    ham_or_trott::Union{HamHam, AbstractTrotter},
    config::Config{Thermalize, D},
    precomputed_data,
    scratch::ThermalizeScratch{<:Complex},
    ) where {D<:Union{TimeDomain, TrotterDomain}}
    dim = size(jumps[1].in_eigenbasis, 1)
    (; transition, gamma_norm_factor, energy_labels, oft_nufft_prefactors) = precomputed_data

    n_labels = length(energy_labels)
    if Threads.nthreads() > 1 && n_labels >= OMEGA_THREAD_THRESHOLD
        return _precompute_R_threaded_timetrot!(jumps, ham_or_trott, config, precomputed_data, scratch)
    end

    # Same weight as in jump_contribution!(::Union{TimeDomain,TrotterDomain}, ...), but without delta.
    base_prefactor = precomputed_data.oft_domain_prefactor * gamma_norm_factor

    fill!(scratch.R, 0)

    @inbounds for jump in jumps
        if jump.hermitian
            # Half-grid and explicit mirror, allocation-free (no filter/abs vector creation).
            for w_raw in energy_labels
                w_raw > 1e-12 && continue
                w = abs(w_raw)

                nufft_prefactor_matrix = _prefactor_view(oft_nufft_prefactors, w)

                # Aw := A .* prefactor_matrix(w)  (elementwise)
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
                nufft_prefactor_matrix = _prefactor_view(oft_nufft_prefactors, w)
                @. scratch.jump_oft = jump.in_eigenbasis * nufft_prefactor_matrix

                rate2 = base_prefactor * transition(w)
                mul!(scratch.LdagL, scratch.jump_oft', scratch.jump_oft)
                @. scratch.R += rate2 * scratch.LdagL
            end
        end
    end

    hermitianize!(scratch.R)
    return scratch.R
end

"""
    _precompute_R(jumps, hamiltonian, config::Config{Thermalize, BohrDomain}, precomputed_data, scratch)

Accumulate the Bohr-domain rate operator.

Math: \$R = sum_(nu_2) A_(nu_2)^dagger sum_(nu_1) alpha_(nu_1,nu_2) A_(nu_1)\$.
Per-jump probability rescaling is left to the caller.

# Returns
Hermitianized `scratch.R`.
"""
function _precompute_R(
    jumps::AbstractVector{<:JumpOp},
    hamiltonian::HamHam,
    config::Config{Thermalize, BohrDomain},
    precomputed_data,
    scratch::ThermalizeScratch{<:Complex},
)
    dim = size(hamiltonian.data, 1)
    (; alpha, gamma_norm_factor) = precomputed_data

    bohr_keys = hasproperty(precomputed_data, :bohr_keys) ? precomputed_data.bohr_keys : collect(keys(hamiltonian.bohr_dict))
    bohr_is   = hasproperty(precomputed_data, :bohr_is)   ? precomputed_data.bohr_is   : nothing
    bohr_js   = hasproperty(precomputed_data, :bohr_js)   ? precomputed_data.bohr_js   : nothing

    n_keys = length(bohr_keys)
    if Threads.nthreads() > 1 && n_keys >= OMEGA_THREAD_THRESHOLD
        return _precompute_R_threaded_bohr!(jumps, hamiltonian, config, precomputed_data, scratch, bohr_keys, bohr_is, bohr_js)
    end

    fill!(scratch.R, 0)

    # Keep the hot-loop matrix view concretely typed.
    CT = eltype(scratch.jump_oft)

    @inbounds for jump in jumps
        in_eb = jump.in_eigenbasis::Matrix{CT}
        for (k, nu_2) in pairs(bohr_keys)
            # B_{nu_2} = sum_{nu_1} alpha(nu_1, nu_2) * A^a
            @. scratch.jump_oft = alpha(hamiltonian.bohr_freqs, nu_2) * in_eb

            # R += gamma_norm_factor * A_{nu_2}^dagger * B_{nu_2}
            if bohr_is !== nothing
                is = bohr_is[k]
                js = bohr_js[k]
                @inbounds for t in eachindex(is)
                    i = is[t]
                    j = js[t]
                    v = conj(in_eb[i, j]) * gamma_norm_factor
                    @inbounds for q in 1:dim
                        scratch.R[j, q] += v * scratch.jump_oft[i, q]
                    end
                end
            else
                indices = hamiltonian.bohr_dict[nu_2]
                @inbounds for idx in indices
                    i = idx[1]
                    j = idx[2]
                    v = conj(in_eb[i, j]) * gamma_norm_factor
                    @inbounds for q in 1:dim
                        scratch.R[j, q] += v * scratch.jump_oft[i, q]
                    end
                end
            end
        end
    end

    hermitianize!(scratch.R)
    return scratch.R
end

#* Threaded _precompute_R variants (THREAD-02) -------------------------------------------------------

# --- EnergyDomain threaded _precompute_R ---

function _precompute_R_threaded_energy!(
    jumps::AbstractVector{<:JumpOp},
    hamiltonian::HamHam,
    config::Config{Thermalize, EnergyDomain},
    precomputed_data,
    scratch::ThermalizeScratch{CT},
) where {CT<:Complex}
    dim = size(hamiltonian.data, 1)
    (; transition, gamma_norm_factor, energy_labels) = precomputed_data
    base_prefactor = precomputed_data.oft_domain_prefactor * gamma_norm_factor
    inv_4sigma2 = 1.0 / (4 * config.sigma^2)

    fill!(scratch.R, 0)

    @inbounds for jump in jumps
        # For Hermitian jumps, pre-filter to half-grid indices
        if jump.hermitian
            half_indices = [i for i in eachindex(energy_labels) if energy_labels[i] <= 1e-12]
        else
            half_indices = collect(eachindex(energy_labels))
        end

        n_work = length(half_indices)
        nt = min(Threads.nthreads(), n_work)
        chunks = _partition_range(1:n_work, nt)

        task_scratches = [ThermalizeScratch(CT, dim) for _ in 1:length(chunks)]

        @sync for (idx, chunk) in enumerate(chunks)
            Threads.@spawn _precompute_R_chunk_energy!(
                task_scratches[idx], jump, hamiltonian, precomputed_data,
                half_indices[chunk];
                base_prefactor=base_prefactor, inv_4sigma2=inv_4sigma2)
        end

        # Reduce: sum per-task R into scratch.R (additive across jumps)
        for ts in task_scratches
            scratch.R .+= ts.R
        end
    end

    hermitianize!(scratch.R)
    return scratch.R
end

function _precompute_R_chunk_energy!(
    scratch::ThermalizeScratch{CT},
    jump::JumpOp,
    hamiltonian::HamHam,
    precomputed_data,
    label_indices::AbstractVector{Int};
    base_prefactor::Real,
    inv_4sigma2::Real,
) where {CT<:Complex}
    (; transition, energy_labels) = precomputed_data
    fill!(scratch.R, 0)

    @inbounds for li in label_indices
        w_raw = energy_labels[li]
        w = jump.hermitian ? abs(w_raw) : w_raw

        oft!(scratch.jump_oft, jump.in_eigenbasis, hamiltonian.bohr_freqs, w, inv_4sigma2)

        rate2_pos = base_prefactor * transition(w)
        mul!(scratch.LdagL, scratch.jump_oft', scratch.jump_oft)
        @. scratch.R += rate2_pos * scratch.LdagL

        if jump.hermitian && w > 1e-12
            rate2_neg = base_prefactor * transition(-w)
            mul!(scratch.LdagL, scratch.jump_oft, scratch.jump_oft')
            @. scratch.R += rate2_neg * scratch.LdagL
        end
    end

    return nothing
end

# --- TimeDomain/TrotterDomain threaded _precompute_R ---

function _precompute_R_threaded_timetrot!(
    jumps::AbstractVector{<:JumpOp},
    ham_or_trott::Union{HamHam, AbstractTrotter},
    config::Config{Thermalize, D},
    precomputed_data,
    scratch::ThermalizeScratch{CT},
) where {CT<:Complex, D<:Union{TimeDomain, TrotterDomain}}
    dim = size(jumps[1].in_eigenbasis, 1)
    (; transition, gamma_norm_factor, energy_labels, oft_nufft_prefactors) = precomputed_data
    base_prefactor = precomputed_data.oft_domain_prefactor * gamma_norm_factor

    fill!(scratch.R, 0)

    @inbounds for jump in jumps
        if jump.hermitian
            half_indices = [i for i in eachindex(energy_labels) if energy_labels[i] <= 1e-12]
        else
            half_indices = collect(eachindex(energy_labels))
        end

        n_work = length(half_indices)
        nt = min(Threads.nthreads(), n_work)
        chunks = _partition_range(1:n_work, nt)

        task_scratches = [ThermalizeScratch(CT, dim) for _ in 1:length(chunks)]

        @sync for (idx, chunk) in enumerate(chunks)
            Threads.@spawn _precompute_R_chunk_timetrot!(
                task_scratches[idx], jump, precomputed_data,
                half_indices[chunk];
                base_prefactor=base_prefactor)
        end

        for ts in task_scratches
            scratch.R .+= ts.R
        end
    end

    hermitianize!(scratch.R)
    return scratch.R
end

function _precompute_R_chunk_timetrot!(
    scratch::ThermalizeScratch{CT},
    jump::JumpOp,
    precomputed_data,
    label_indices::AbstractVector{Int};
    base_prefactor::Real,
) where {CT<:Complex}
    (; transition, energy_labels, oft_nufft_prefactors) = precomputed_data
    fill!(scratch.R, 0)

    @inbounds for li in label_indices
        w_raw = energy_labels[li]
        w = jump.hermitian ? abs(w_raw) : w_raw

        nufft_prefactor_matrix = _prefactor_view(oft_nufft_prefactors, w)
        @. scratch.jump_oft = jump.in_eigenbasis * nufft_prefactor_matrix

        rate2_pos = base_prefactor * transition(w)
        mul!(scratch.LdagL, scratch.jump_oft', scratch.jump_oft)
        @. scratch.R += rate2_pos * scratch.LdagL

        if jump.hermitian && w > 1e-12
            rate2_neg = base_prefactor * transition(-w)
            mul!(scratch.LdagL, scratch.jump_oft, scratch.jump_oft')
            @. scratch.R += rate2_neg * scratch.LdagL
        end
    end

    return nothing
end

# --- BohrDomain threaded _precompute_R ---

function _precompute_R_threaded_bohr!(
    jumps::AbstractVector{<:JumpOp},
    hamiltonian::HamHam,
    config::Config{Thermalize, BohrDomain},
    precomputed_data,
    scratch::ThermalizeScratch{CT},
    bohr_keys::AbstractVector,
    bohr_is::Union{Nothing, Vector{Vector{Int}}},
    bohr_js::Union{Nothing, Vector{Vector{Int}}},
) where {CT<:Complex}
    dim = size(hamiltonian.data, 1)
    (; alpha, gamma_norm_factor) = precomputed_data

    fill!(scratch.R, 0)

    @inbounds for jump in jumps
        n_keys = length(bohr_keys)
        nt = min(Threads.nthreads(), n_keys)
        chunks = _partition_range(1:n_keys, nt)

        task_scratches = [ThermalizeScratch(CT, dim) for _ in 1:length(chunks)]

        @sync for (idx, chunk) in enumerate(chunks)
            Threads.@spawn _precompute_R_chunk_bohr!(
                task_scratches[idx], jump, hamiltonian, precomputed_data,
                bohr_keys, bohr_is, bohr_js, chunk;
                gamma_norm_factor=gamma_norm_factor)
        end

        for ts in task_scratches
            scratch.R .+= ts.R
        end
    end

    hermitianize!(scratch.R)
    return scratch.R
end

function _precompute_R_chunk_bohr!(
    scratch::ThermalizeScratch{CT},
    jump::JumpOp,
    hamiltonian::HamHam,
    precomputed_data,
    bohr_keys::AbstractVector,
    bohr_is::Union{Nothing, Vector{Vector{Int}}},
    bohr_js::Union{Nothing, Vector{Vector{Int}}},
    key_indices::UnitRange{Int};
    gamma_norm_factor::Real,
) where {CT<:Complex}
    dim = size(hamiltonian.data, 1)
    (; alpha) = precomputed_data
    fill!(scratch.R, 0)

    # Keep the hot-loop matrix view concretely typed.
    in_eb = jump.in_eigenbasis::Matrix{CT}

    @inbounds for k in key_indices
        nu_2 = bohr_keys[k]
        @. scratch.jump_oft = alpha(hamiltonian.bohr_freqs, nu_2) * in_eb

        if bohr_is !== nothing
            is = bohr_is[k]
            js = bohr_js[k]
            @inbounds for t in eachindex(is)
                i = is[t]
                j = js[t]
                v = conj(in_eb[i, j]) * gamma_norm_factor
                @inbounds for q in 1:dim
                    scratch.R[j, q] += v * scratch.jump_oft[i, q]
                end
            end
        else
            indices = hamiltonian.bohr_dict[nu_2]
            @inbounds for idx in indices
                i = idx[1]
                j = idx[2]
                v = conj(in_eb[i, j]) * gamma_norm_factor
                @inbounds for q in 1:dim
                    scratch.R[j, q] += v * scratch.jump_oft[i, q]
                end
            end
        end
    end

    return nothing
end
