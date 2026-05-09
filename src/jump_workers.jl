#* Liouvillian (vectorized) jump contributions ------------------------------------------------------------------------------
"""
    Accumulate the Liouvillian contribution of a single jump operator in-place.

    This avoids allocating a full `dim^2 x dim^2` matrix per jump. Call with a
    preallocated `L_target` (dense) and a `Workspace{Lindbladian}`.

    If `with_coherent(config.construction)==true`, pass `coherent_term` already scaled by
    `gamma_norm_factor` to avoid modifying cached matrices.
"""
function _jump_contribution!(
    L_target::AbstractMatrix{<:Complex},
    jump::JumpOp,
    hamiltonian::HamHam,
    config::Config{Lindbladian, BohrDomain},
    precomputed_data,
    ws::Workspace{Lindbladian};
    coherent_term::Union{Nothing, AbstractMatrix{<:Complex}} = nothing,
    )
    dim = size(hamiltonian.data, 1)
    unique_freqs = keys(hamiltonian.bohr_dict)
    (; alpha, gamma_norm_factor) = precomputed_data

    B = coherent_term
    if B !== nothing
        _vectorize_liouvillian_coherent!(L_target, B, ws)
    end

    alpha_A_nu1 = ws.scratch.jump_tmp
    for nu_2 in unique_freqs
        @. alpha_A_nu1 = alpha(hamiltonian.bohr_freqs, nu_2) * jump.in_eigenbasis

        indices = hamiltonian.bohr_dict[nu_2]
        A_nu_2_vals = view(jump.in_eigenbasis, indices)

        # swapped for dagger
        rows_dag = getindex.(indices, 2)
        cols_dag = getindex.(indices, 1)

        A_nu_2_dag = sparse(rows_dag, cols_dag, conj.(A_nu_2_vals), dim, dim)

        _vectorize_liouv_diss_and_add!(L_target, alpha_A_nu1, A_nu_2_dag, gamma_norm_factor, ws)
    end
    return L_target
end

function _jump_contribution!(
    L_target::AbstractMatrix{<:Complex},
    jump::JumpOp,
    hamiltonian::HamHam,
    config::Config{Lindbladian, EnergyDomain},
    precomputed_data,
    ws::Workspace{Lindbladian};
    coherent_term::Union{Nothing, AbstractMatrix{<:Complex}} = nothing,
    )

    (; transition, gamma_norm_factor, energy_labels) = precomputed_data

    B = coherent_term
    if B !== nothing
        _vectorize_liouvillian_coherent!(L_target, B, ws)
    end

    jump_oft = ws.scratch.jump_tmp
    prefactor = precomputed_data.oft_domain_prefactor * gamma_norm_factor
    inv_4sigma2 = 1.0 / (4 * config.sigma^2)

    if jump.hermitian
        for w_raw in energy_labels
            # iterate only half-grid (w<=0) and mirror manually
            w_raw > 1e-12 && continue
            w = abs(w_raw)
            oft!(jump_oft, jump.in_eigenbasis, hamiltonian.bohr_freqs, w, inv_4sigma2)
            scalar_w = prefactor * transition(w)
            _vectorize_liouv_diss_and_add!(L_target, jump_oft, scalar_w, ws)
            if w > 1e-12
                scalar_negative_w = prefactor * transition(-w)
                _vectorize_liouv_diss_and_add!(L_target, jump_oft', scalar_negative_w, ws)
            end
        end
    else
        for w in energy_labels
            oft!(jump_oft, jump.in_eigenbasis, hamiltonian.bohr_freqs, w, inv_4sigma2)
            scalar_w = prefactor * transition(w)
            _vectorize_liouv_diss_and_add!(L_target, jump_oft, scalar_w, ws)
        end
    end

    return L_target
end

"""
    _jump_contribution!(
        L_target, jump, hamiltonian,
        config::Config{Lindbladian, BohrDomain, DLL}, precomputed_data, ws;
        coherent_term=nothing,
    )

DLL Bohr-domain Liouvillian contribution (Ding–Li–Lin 2024, Eqs. 3.4 / 3.8).
Builds the single Lindblad operator `L_a = Σ_ν freq_kernel(filter, ν) A^a_ν`
once via `dll_lindblad_op_bohr` and accumulates the dissipator
`L_a ρ L_a† − (1/2){L_a† L_a, ρ}` into `L_target` with no outer ω-loop and
no γ(ω) prefactor.

For multi-channel DLL filters (qf-7go.1), `dll_lindblad_op_bohr` returns a
length-`k` vector of per-channel operators; this loop sums the per-channel
dissipators (no cross terms in the multi-channel α).

The optional `coherent_term` plumbing matches the CKG signature; the actual
DLL coherent operator (G in Eq. 3.8) is built in DLL-3, so callers currently
pass `coherent_term=nothing`.
"""
function _jump_contribution!(
    L_target::AbstractMatrix{<:Complex},
    jump::JumpOp,
    hamiltonian::HamHam,
    config::Config{Lindbladian, BohrDomain, DLL},
    precomputed_data,
    ws::Workspace{Lindbladian};
    coherent_term::Union{Nothing, AbstractMatrix{<:Complex}} = nothing,
    )
    (; filter) = precomputed_data

    if coherent_term !== nothing
        _vectorize_liouvillian_coherent!(L_target, coherent_term, ws)
    end

    _accumulate_dll_bohr_dissipator!(L_target, jump, hamiltonian, filter, ws)
    return L_target
end

# Single-channel DLL filter: one Lindblad operator per coupling.
@inline function _accumulate_dll_bohr_dissipator!(
    L_target::AbstractMatrix{<:Complex},
    jump::JumpOp,
    hamiltonian::HamHam,
    filter::AbstractFilter,
    ws::Workspace{Lindbladian},
)
    L_a = dll_lindblad_op_bohr(jump, hamiltonian, filter)
    _vectorize_liouv_diss_and_add!(L_target, L_a, 1.0, ws)
    return L_target
end

# Multi-channel DLL filter (qf-7go.1): see `src/dll_multichannel.jl` for the
# `_accumulate_dll_bohr_dissipator!(::DLLMultiChannelFilter, ...)` overload.

"""
    _jump_contribution!(
        L_target, jump, hamiltonian,
        config::Config{Lindbladian, TimeDomain, DLL}, precomputed_data, ws;
        coherent_term=nothing,
    )

DLL Time-domain Liouvillian contribution (Ding–Li–Lin 2024, Eqs. 3.4 / 3.8 /
3.13). The single Lindblad operator on the simulator's truncated time grid is

    L_a[i, j] = A_eb[i, j] · τ · Σ_m time_kernel(filter, t_m) · cis((λ_i − λ_j) · t_m)

The DFT factor `Σ_m … cis((λ_i − λ_j) · t_m)` is precomputed once per
Liouvillian build via FINUFFT in `_precompute_data`
(`oft_nufft_at_zero_list`), collapsing the per-jump cost to a single
elementwise multiply. Same Riemann sum, same Eq. 3.15 quadrature error
structure as the explicit `dll_lindblad_op_time`; only the DFT is now
`O(Nt log Nt + n² log(1/ε))`.

For multi-channel DLL filters (qf-7go.1) the list has one prefactor matrix
per channel; the dissipator sums `Σ_ℓ L^(ℓ) ρ (L^(ℓ))† − …` with no cross
terms.

The optional `coherent_term` plumbing matches the CKG signature; DLL coherent
`G` is wired through `_precompute_coherent_B`.
"""
function _jump_contribution!(
    L_target::AbstractMatrix{<:Complex},
    jump::JumpOp,
    hamiltonian::HamHam,
    config::Config{Lindbladian, TimeDomain, DLL},
    precomputed_data,
    ws::Workspace{Lindbladian};
    coherent_term::Union{Nothing, AbstractMatrix{<:Complex}} = nothing,
    )
    (; t0, oft_nufft_at_zero_list) = precomputed_data

    if coherent_term !== nothing
        _vectorize_liouvillian_coherent!(L_target, coherent_term, ws)
    end

    L_a = ws.scratch.jump_tmp
    @inbounds for nufft_at_zero in oft_nufft_at_zero_list
        @. L_a = jump.in_eigenbasis * nufft_at_zero * t0
        _vectorize_liouv_diss_and_add!(L_target, L_a, 1.0, ws)
    end
    return L_target
end

function _jump_contribution!(
    L_target::AbstractMatrix{<:Complex},
    jump::JumpOp,
    ham_or_trott::Union{HamHam, TrottTrott},
    config::Config{Lindbladian, D},
    precomputed_data,
    ws::Workspace{Lindbladian};
    coherent_term::Union{Nothing, AbstractMatrix{<:Complex}} = nothing,
    ) where {D<:Union{TimeDomain, TrotterDomain}}

    (; transition, gamma_norm_factor, energy_labels, oft_nufft_prefactors, b_minus, b_plus) = precomputed_data

    B = coherent_term
    if B !== nothing
        _vectorize_liouvillian_coherent!(L_target, B, ws)
    end

    jump_oft = ws.scratch.jump_tmp
    prefactor = precomputed_data.oft_domain_prefactor * gamma_norm_factor

    if jump.hermitian
        for w_raw in energy_labels
            w_raw > 1e-12 && continue
            w = abs(w_raw)
            nufft_prefactor_matrix = _prefactor_view(oft_nufft_prefactors, w)
            @. jump_oft = jump.in_eigenbasis * nufft_prefactor_matrix

            scalar_w = prefactor * transition(w)
            _vectorize_liouv_diss_and_add!(L_target, jump_oft, scalar_w, ws)
            if w > 1e-12
                scalar_negative_w = prefactor * transition(-w)
                _vectorize_liouv_diss_and_add!(L_target, jump_oft', scalar_negative_w, ws)
            end
        end
    else
        for w in energy_labels
            nufft_prefactor_matrix = _prefactor_view(oft_nufft_prefactors, w)
            @. jump_oft = jump.in_eigenbasis * nufft_prefactor_matrix
            scalar_w = prefactor * transition(w)
            _vectorize_liouv_diss_and_add!(L_target, jump_oft, scalar_w, ws)
        end
    end

    return L_target
end

#* Algorithmic jump contributions -------------------------------------------------------------------------------------------

"""
    _apply_coherent_unitary!(evolving_dm, U_B, scratch) -> nothing

Apply coherent unitary evolution: rho -> U_B * rho * U_B'.
No-op if U_B is nothing.
"""
@inline function _apply_coherent_unitary!(
    evolving_dm::Matrix{<:Complex},
    U_B::Union{Nothing,Matrix{<:Complex}},
    scratch::ThermalizeScratch{<:Complex},
)
    U_B === nothing && return nothing
    mul!(scratch.sandwich_tmp, U_B, evolving_dm)
    mul!(scratch.rho_next, scratch.sandwich_tmp, U_B')
    copyto!(evolving_dm, scratch.rho_next)
    return nothing
end

"""
    _finalize_kraus_step!(evolving_dm, delta, scratch) -> evolving_dm

Apply the CPTP weak-measurement channel after R and rho_jump have been accumulated.

Implements Chen Eq. 3.2:
  K0 = I - alpha * R,  alpha = 1 - sqrt(1 - delta)
  S  = (2*alpha - delta)*R - alpha^2 * R^2  (residual, O(delta^2))
  U_residual = sqrt_psd(S)  (eigendecomposition with clamped eigenvalues)
  rho_next = K0 * rho * K0' + rho_jump + U_res * rho * U_res'

Expects scratch.R (Hermitianized) and scratch.rho_jump to be pre-filled by the
domain-specific dissipative accumulation loop.
"""
function _finalize_kraus_step!(
    evolving_dm::Matrix{<:Complex},
    delta::Real,
    scratch::ThermalizeScratch{<:Complex},
)
    # Build CPTP channel from accumulated R (Chen Eq. 3.2)
    (; K0, U_residual) = _build_cptp_channel(scratch.R, delta)

    # rho_next = K0 * rho * K0' + rho_jump + U_res * rho * U_res'
    mul!(scratch.sandwich_tmp, K0, evolving_dm)
    mul!(scratch.rho_next, scratch.sandwich_tmp, K0')
    scratch.rho_next .+= scratch.rho_jump

    mul!(scratch.sandwich_tmp, U_residual, evolving_dm)
    mul!(scratch.rho_next, scratch.sandwich_tmp, U_residual', 1.0, 1.0)

    # Keep it a density matrix numerically
    hermitianize!(scratch.rho_next)
    copyto!(evolving_dm, scratch.rho_next)
    return evolving_dm
end

function _jump_contribution!(
    evolving_dm::Matrix{<:Complex},
    jump::JumpOp,
    hamiltonian::HamHam,
    config::Config{Thermalize, BohrDomain},
    precomputed_data,
    scratch::ThermalizeScratch{<:Complex};
    coherent_unitary_cache::Union{Nothing,Matrix{<:Complex}} = nothing,
    jump_prob::Real = 1.0,
    rescale_by_inv_prob::Bool = false
    )

    dim = size(evolving_dm, 1)
    (; alpha, gamma_norm_factor) = precomputed_data

    bohr_keys = hasproperty(precomputed_data, :bohr_keys) ? precomputed_data.bohr_keys : collect(keys(hamiltonian.bohr_dict))
    bohr_is   = hasproperty(precomputed_data, :bohr_is)   ? precomputed_data.bohr_is   : nothing
    bohr_js   = hasproperty(precomputed_data, :bohr_js)   ? precomputed_data.bohr_js   : nothing

    jump_weight_scaling = rescale_by_inv_prob ? (gamma_norm_factor / jump_prob) : gamma_norm_factor
    scaled_delta = config.delta * jump_weight_scaling

    _apply_coherent_unitary!(evolving_dm, coherent_unitary_cache, scratch)

    fill!(scratch.R, 0)
    fill!(scratch.rho_jump, 0)

    # For each fixed "right" Bohr label v2 build the composite
    #   B_{v2} = \sum_{v1} alpha_{v1,v2} A_{v1}
    # and accumulate
    #   rho_jump += delta * B_{v2} rho A_{v2}dag
    #   R        +=     A_{v2}dag B_{v2}
    @inbounds for (k, nu_2) in pairs(bohr_keys)
        # B_{v2}
        @. scratch.jump_oft = alpha(hamiltonian.bohr_freqs, nu_2) * jump.in_eigenbasis

        # sandwich_tmp := rho A_{v2}dag without explicitly building A_{v2}dag.
        # If (i,j) is in the v2 bucket, then (A_{v2}dag)_{j,i} = conj(A_{i,j}).
        fill!(scratch.sandwich_tmp, 0)
        if bohr_is !== nothing
            is = bohr_is[k]
            js = bohr_js[k]
            @inbounds for t in eachindex(is)
                i = is[t]
                j = js[t]
                v = conj(jump.in_eigenbasis[i, j])
                @inbounds for p in 1:dim
                    scratch.sandwich_tmp[p, i] += evolving_dm[p, j] * v
                end
            end
        else
            indices = hamiltonian.bohr_dict[nu_2]
            @inbounds for idx in indices
                i = idx[1]
                j = idx[2]
                v = conj(jump.in_eigenbasis[i, j])
                @inbounds for p in 1:dim
                    scratch.sandwich_tmp[p, i] += evolving_dm[p, j] * v
                end
            end
        end

        # rho_jump += delta * B_{v2} * (rho A_{v2}dag)
        mul!(scratch.rho_jump, scratch.jump_oft, scratch.sandwich_tmp, scaled_delta, 1.0)

        # R += A_{v2}dag * B_{v2}  (no delta factor)
        if bohr_is !== nothing
            is = bohr_is[k]
            js = bohr_js[k]
            @inbounds for t in eachindex(is)
                i = is[t]
                j = js[t]
                v = conj(jump.in_eigenbasis[i, j]) * jump_weight_scaling
                @inbounds for q in 1:dim
                    scratch.R[j, q] += v * scratch.jump_oft[i, q]
                end
            end
        else
            indices = hamiltonian.bohr_dict[nu_2]
            @inbounds for idx in indices
                i = idx[1]
                j = idx[2]
                v = conj(jump.in_eigenbasis[i, j]) * jump_weight_scaling
                @inbounds for q in 1:dim
                    scratch.R[j, q] += v * scratch.jump_oft[i, q]
                end
            end
        end
    end

    # Hermitianize R (numerical)
    hermitianize!(scratch.R)

    # Apply R, K0, U_residual
    _finalize_kraus_step!(evolving_dm, config.delta, scratch)
    return evolving_dm
end

function _jump_contribution!(
    evolving_dm::Matrix{<:Complex},
    jump::JumpOp,
    hamiltonian::HamHam,
    config::Config{Thermalize, EnergyDomain},
    precomputed_data,
    scratch::ThermalizeScratch{<:Complex};
    coherent_unitary_cache::Union{Nothing,Matrix{<:Complex}} = nothing,
    jump_prob::Real = 1.0,
    rescale_by_inv_prob::Bool = false
    )

    dim = size(evolving_dm, 1)
    (; transition, gamma_norm_factor, energy_labels) = precomputed_data

    jump_weight_scaling = rescale_by_inv_prob ? (gamma_norm_factor / jump_prob) : gamma_norm_factor

    _apply_coherent_unitary!(evolving_dm, coherent_unitary_cache, scratch)

    # --- Dissipative part ---
    base_prefactor = precomputed_data.oft_domain_prefactor * jump_weight_scaling
    inv_4sigma2 = 1.0 / (4 * config.sigma^2)

    fill!(scratch.R, 0)
    fill!(scratch.rho_jump, 0)

    if jump.hermitian
        @inbounds for w_raw in energy_labels
            w_raw > 1e-12 && continue
            w = abs(w_raw)

            # Aw
            oft!(scratch.jump_oft, jump.in_eigenbasis, hamiltonian.bohr_freqs, w, inv_4sigma2)

            rate2_pos = base_prefactor * transition(w)

            # R += rate^2 * (Aw_dag Aw)
            mul!(scratch.LdagL, scratch.jump_oft', scratch.jump_oft)
            @. scratch.R += rate2_pos * scratch.LdagL

            # rho_jump += delta * rate^2 * (Aw rho Aw_dag)
            mul!(scratch.sandwich_tmp, evolving_dm, scratch.jump_oft')  # rho Aw_dag
            mul!(scratch.rho_jump, scratch.jump_oft, scratch.sandwich_tmp, config.delta * rate2_pos, 1.0)

            if w > 1e-12
                rate2_neg = base_prefactor * transition(-w)

                # Negative-frequency partner uses (Aw)_dag as Lindblad operator.
                mul!(scratch.LdagL, scratch.jump_oft, scratch.jump_oft')
                @. scratch.R += rate2_neg * scratch.LdagL

                mul!(scratch.sandwich_tmp, evolving_dm, scratch.jump_oft)  # rho Aw
                mul!(scratch.rho_jump, scratch.jump_oft', scratch.sandwich_tmp, config.delta * rate2_neg, 1.0)
            end
        end
    else
        @inbounds for w in energy_labels
            oft!(scratch.jump_oft, jump.in_eigenbasis, hamiltonian.bohr_freqs, w, inv_4sigma2)

            rate2 = base_prefactor * transition(w)

            mul!(scratch.LdagL, scratch.jump_oft', scratch.jump_oft)
            @. scratch.R += rate2 * scratch.LdagL

            mul!(scratch.sandwich_tmp, evolving_dm, scratch.jump_oft')
            mul!(scratch.rho_jump, scratch.jump_oft, scratch.sandwich_tmp, config.delta * rate2, 1.0)
        end
    end

    # Hermitianize R (numerical)
    hermitianize!(scratch.R)

    # Apply R, K0, U_residual
    _finalize_kraus_step!(evolving_dm, config.delta, scratch)
    return evolving_dm
end

function _jump_contribution!(
    evolving_dm::Matrix{<:Complex},
    jump::JumpOp,
    ham_or_trott,              # HamHam or TrottTrott depending on domain
    config::Config{Thermalize, D},
    precomputed_data,
    scratch::ThermalizeScratch{<:Complex};
    coherent_unitary_cache::Union{Nothing,Matrix{<:Complex}} = nothing,
    jump_prob::Real = 1.0,
    rescale_by_inv_prob::Bool = false
    ) where {D<:Union{TimeDomain, TrotterDomain}}

    dim = size(evolving_dm, 1)
    (; transition, gamma_norm_factor, energy_labels, oft_nufft_prefactors, b_minus, b_plus) = precomputed_data

    jump_weight_scaling = rescale_by_inv_prob ? (gamma_norm_factor / jump_prob) : gamma_norm_factor

    _apply_coherent_unitary!(evolving_dm, coherent_unitary_cache, scratch)

    base_prefactor = precomputed_data.oft_domain_prefactor * jump_weight_scaling

    fill!(scratch.R, 0)
    fill!(scratch.rho_jump, 0)

    if jump.hermitian
        @inbounds for w_raw in energy_labels
            w_raw > 1e-12 && continue
            w = abs(w_raw)

            nufft_prefactor_matrix = _prefactor_view(oft_nufft_prefactors, w)
            @. scratch.jump_oft = jump.in_eigenbasis * nufft_prefactor_matrix

            rate2_pos = base_prefactor * transition(w)

            mul!(scratch.LdagL, scratch.jump_oft', scratch.jump_oft)
            @. scratch.R += rate2_pos * scratch.LdagL

            mul!(scratch.sandwich_tmp, evolving_dm, scratch.jump_oft')
            mul!(scratch.rho_jump, scratch.jump_oft, scratch.sandwich_tmp, config.delta*rate2_pos, 1.0)

            if w > 1e-12
                rate2_neg = base_prefactor * transition(-w)

                mul!(scratch.LdagL, scratch.jump_oft, scratch.jump_oft')
                @. scratch.R += rate2_neg * scratch.LdagL

                mul!(scratch.sandwich_tmp, evolving_dm, scratch.jump_oft)
                mul!(scratch.rho_jump, scratch.jump_oft', scratch.sandwich_tmp, config.delta*rate2_neg, 1.0)
            end
        end
    else
        for w in energy_labels
            nufft_prefactor_matrix = _prefactor_view(oft_nufft_prefactors, w)
            @. scratch.jump_oft = jump.in_eigenbasis * nufft_prefactor_matrix

            rate2_pos = base_prefactor * transition(w)

            mul!(scratch.LdagL, scratch.jump_oft', scratch.jump_oft)
            @. scratch.R += rate2_pos * scratch.LdagL

            mul!(scratch.sandwich_tmp, evolving_dm, scratch.jump_oft')
            mul!(scratch.rho_jump, scratch.jump_oft, scratch.sandwich_tmp, config.delta*rate2_pos, 1.0)
        end
    end

    # Hermitianize R
    hermitianize!(scratch.R)

    # Apply R, K0, U_residual
    _finalize_kraus_step!(evolving_dm, config.delta, scratch)
    return evolving_dm
end

#* Omega-loop threading infrastructure ---------------------------------------------------------------

# Minimum number of frequency labels to enable omega-loop parallelism.
# Below this threshold, serial execution is faster due to task spawn overhead.
#
# 2026-05-06 (qf-in3.4) — set by empirical sweep
# (`scripts/scratch_omega_threading_threshold.jl`) at n ∈ {3, 4, 5},
# nthreads=4, EnergyDomain & TimeDomain:
#  - n=3 (dim=8): serial wins at N=5 by ~20%; thread wins from N=10.
#  - n=4 (dim=16): thread wins from N=5 by 2×, climbing to 3.3× at N≥75.
#  - n=5 (dim=32): thread wins from N=5 by 3.3×, climbing to 3.7× at N≥75.
# The new value of 10 captures the n=3 crossover where overhead is largest
# in absolute terms; n ≥ 4 is then strictly faster threaded. The threading
# dispatch costs ~1 kB / matvec for `Threads.@spawn` Task objects, which is
# accepted in exchange for the 2–3× wall-time speedup for n_labels ∈ [10, 49]
# that the previous (folkloric) value of 50 was forfeiting.
const OMEGA_THREAD_THRESHOLD = 10

"""
    _partition_range(range, n_chunks) -> Vector{UnitRange{Int}}

Partition a range into approximately equal chunks for parallel execution.
Same algorithm as _partition_trajectories in trajectories.jl.
"""
function _partition_range(range::UnitRange{Int}, n_chunks::Int)
    len = length(range)
    n_chunks = min(n_chunks, len)
    base = div(len, n_chunks)
    remainder = rem(len, n_chunks)
    chunks = Vector{UnitRange{Int}}(undef, n_chunks)
    start = first(range)
    for i in 1:n_chunks
        chunk_size = base + (i <= remainder ? 1 : 0)
        chunks[i] = start:(start + chunk_size - 1)
        start += chunk_size
    end
    return chunks
end

#* Per-jump rho_jump-only accumulation (precomputed channel path) ----------------------------------------

"""
    _accumulate_rho_jump!(scratch, evolving_dm, jump, hamiltonian, config::Config{Thermalize, EnergyDomain},
                          precomputed_data; jump_weight_scaling)

Accumulate rho_jump = delta * sum_w rate^2(w) * L_{a,w} * rho * L_{a,w}^dagger for EnergyDomain.

Extracts ONLY the rho_jump accumulation from the EnergyDomain `_jump_contribution!`,
with no R or LdagL computation. Used with precomputed K0/U_residual channels.
"""
function _accumulate_rho_jump!(
    scratch::ThermalizeScratch{<:Complex},
    evolving_dm::Matrix{<:Complex},
    jump::JumpOp,
    hamiltonian::HamHam,
    config::Config{Thermalize, EnergyDomain},
    precomputed_data;
    jump_weight_scaling::Real,
)
    (; transition, energy_labels) = precomputed_data

    n_labels = length(energy_labels)
    if Threads.nthreads() > 1 && n_labels >= OMEGA_THREAD_THRESHOLD
        return _accumulate_rho_jump_threaded_energy!(scratch, evolving_dm, jump, hamiltonian, config, precomputed_data; jump_weight_scaling=jump_weight_scaling)
    end

    base_prefactor = precomputed_data.oft_domain_prefactor * jump_weight_scaling
    inv_4sigma2 = 1.0 / (4 * config.sigma^2)

    fill!(scratch.rho_jump, 0)

    if jump.hermitian
        @inbounds for w_raw in energy_labels
            w_raw > 1e-12 && continue
            w = abs(w_raw)

            oft!(scratch.jump_oft, jump.in_eigenbasis, hamiltonian.bohr_freqs, w, inv_4sigma2)

            rate2_pos = base_prefactor * transition(w)

            # rho_jump += delta * rate^2 * (Aw rho Aw_dag)
            mul!(scratch.sandwich_tmp, evolving_dm, scratch.jump_oft')
            mul!(scratch.rho_jump, scratch.jump_oft, scratch.sandwich_tmp, config.delta * rate2_pos, 1.0)

            if w > 1e-12
                rate2_neg = base_prefactor * transition(-w)

                # Negative-frequency partner uses (Aw)_dag as Lindblad operator.
                mul!(scratch.sandwich_tmp, evolving_dm, scratch.jump_oft)
                mul!(scratch.rho_jump, scratch.jump_oft', scratch.sandwich_tmp, config.delta * rate2_neg, 1.0)
            end
        end
    else
        @inbounds for w in energy_labels
            oft!(scratch.jump_oft, jump.in_eigenbasis, hamiltonian.bohr_freqs, w, inv_4sigma2)

            rate2 = base_prefactor * transition(w)

            mul!(scratch.sandwich_tmp, evolving_dm, scratch.jump_oft')
            mul!(scratch.rho_jump, scratch.jump_oft, scratch.sandwich_tmp, config.delta * rate2, 1.0)
        end
    end

    return nothing
end

"""
    _accumulate_rho_jump!(scratch, evolving_dm, jump, ham_or_trott, config::Config{Thermalize, D},
                          precomputed_data; jump_weight_scaling) where D<:Union{TimeDomain, TrotterDomain}

Accumulate rho_jump for TimeDomain/TrotterDomain. Uses NUFFT prefactors for jump_oft
computation. No R or LdagL accumulation.
"""
function _accumulate_rho_jump!(
    scratch::ThermalizeScratch{<:Complex},
    evolving_dm::Matrix{<:Complex},
    jump::JumpOp,
    ham_or_trott,
    config::Config{Thermalize, D},
    precomputed_data;
    jump_weight_scaling::Real,
) where {D<:Union{TimeDomain, TrotterDomain}}
    (; transition, energy_labels, oft_nufft_prefactors) = precomputed_data

    n_labels = length(energy_labels)
    if Threads.nthreads() > 1 && n_labels >= OMEGA_THREAD_THRESHOLD
        return _accumulate_rho_jump_threaded_timetrot!(scratch, evolving_dm, jump, ham_or_trott, config, precomputed_data; jump_weight_scaling=jump_weight_scaling)
    end

    base_prefactor = precomputed_data.oft_domain_prefactor * jump_weight_scaling

    fill!(scratch.rho_jump, 0)

    if jump.hermitian
        @inbounds for w_raw in energy_labels
            w_raw > 1e-12 && continue
            w = abs(w_raw)

            nufft_prefactor_matrix = _prefactor_view(oft_nufft_prefactors, w)
            @. scratch.jump_oft = jump.in_eigenbasis * nufft_prefactor_matrix

            rate2_pos = base_prefactor * transition(w)

            mul!(scratch.sandwich_tmp, evolving_dm, scratch.jump_oft')
            mul!(scratch.rho_jump, scratch.jump_oft, scratch.sandwich_tmp, config.delta * rate2_pos, 1.0)

            if w > 1e-12
                rate2_neg = base_prefactor * transition(-w)

                mul!(scratch.sandwich_tmp, evolving_dm, scratch.jump_oft)
                mul!(scratch.rho_jump, scratch.jump_oft', scratch.sandwich_tmp, config.delta * rate2_neg, 1.0)
            end
        end
    else
        for w in energy_labels
            nufft_prefactor_matrix = _prefactor_view(oft_nufft_prefactors, w)
            @. scratch.jump_oft = jump.in_eigenbasis * nufft_prefactor_matrix

            rate2_pos = base_prefactor * transition(w)

            mul!(scratch.sandwich_tmp, evolving_dm, scratch.jump_oft')
            mul!(scratch.rho_jump, scratch.jump_oft, scratch.sandwich_tmp, config.delta * rate2_pos, 1.0)
        end
    end

    return nothing
end

"""
    _accumulate_rho_jump!(scratch, evolving_dm, jump, hamiltonian, config::Config{Thermalize, BohrDomain},
                          precomputed_data; jump_weight_scaling)

Accumulate rho_jump for BohrDomain. Iterates over Bohr frequency buckets, computing
rho_jump += scaled_delta * B_{nu_2} * (rho * A_{nu_2}^dagger) for each bucket.
No R accumulation.
"""
function _accumulate_rho_jump!(
    scratch::ThermalizeScratch{<:Complex},
    evolving_dm::Matrix{<:Complex},
    jump::JumpOp,
    hamiltonian::HamHam,
    config::Config{Thermalize, BohrDomain},
    precomputed_data;
    jump_weight_scaling::Real,
)
    dim = size(evolving_dm, 1)
    (; alpha) = precomputed_data

    bohr_keys = hasproperty(precomputed_data, :bohr_keys) ? precomputed_data.bohr_keys : collect(keys(hamiltonian.bohr_dict))
    bohr_is   = hasproperty(precomputed_data, :bohr_is)   ? precomputed_data.bohr_is   : nothing
    bohr_js   = hasproperty(precomputed_data, :bohr_js)   ? precomputed_data.bohr_js   : nothing

    n_keys = length(bohr_keys)
    if Threads.nthreads() > 1 && n_keys >= OMEGA_THREAD_THRESHOLD
        return _accumulate_rho_jump_threaded_bohr!(scratch, evolving_dm, jump, hamiltonian, config, precomputed_data, bohr_keys, bohr_is, bohr_js; jump_weight_scaling=jump_weight_scaling)
    end

    scaled_delta = config.delta * jump_weight_scaling

    fill!(scratch.rho_jump, 0)

    @inbounds for (k, nu_2) in pairs(bohr_keys)
        # B_{v2} = sum_{v1} alpha(v1, v2) * A
        @. scratch.jump_oft = alpha(hamiltonian.bohr_freqs, nu_2) * jump.in_eigenbasis

        # sandwich_tmp := rho A_{v2}dag
        fill!(scratch.sandwich_tmp, 0)
        if bohr_is !== nothing
            is = bohr_is[k]
            js = bohr_js[k]
            @inbounds for t in eachindex(is)
                i = is[t]
                j = js[t]
                v = conj(jump.in_eigenbasis[i, j])
                @inbounds for p in 1:dim
                    scratch.sandwich_tmp[p, i] += evolving_dm[p, j] * v
                end
            end
        else
            indices = hamiltonian.bohr_dict[nu_2]
            @inbounds for idx in indices
                i = idx[1]
                j = idx[2]
                v = conj(jump.in_eigenbasis[i, j])
                @inbounds for p in 1:dim
                    scratch.sandwich_tmp[p, i] += evolving_dm[p, j] * v
                end
            end
        end

        # rho_jump += scaled_delta * B_{v2} * (rho A_{v2}dag)
        mul!(scratch.rho_jump, scratch.jump_oft, scratch.sandwich_tmp, scaled_delta, 1.0)
    end

    return nothing
end

#* Threaded omega-loop variants (THREAD-01, THREAD-04) ---------------------------------------------------

# --- EnergyDomain threaded variant ---

function _accumulate_rho_jump_threaded_energy!(
    scratch::ThermalizeScratch{CT},
    evolving_dm::Matrix{CT},
    jump::JumpOp,
    hamiltonian::HamHam,
    config::Config{Thermalize, EnergyDomain},
    precomputed_data;
    jump_weight_scaling::Real,
    task_scratches::Union{Nothing, Vector{ThermalizeScratch{CT}}}=nothing,
) where {CT<:Complex}
    (; transition, energy_labels) = precomputed_data
    base_prefactor = precomputed_data.oft_domain_prefactor * jump_weight_scaling
    inv_4sigma2 = 1.0 / (4 * config.sigma^2)

    # For Hermitian jumps, pre-filter to half-grid indices for balanced partitioning
    if jump.hermitian
        half_indices = [i for i in eachindex(energy_labels) if energy_labels[i] <= 1e-12]
    else
        half_indices = collect(eachindex(energy_labels))
    end

    n_work = length(half_indices)
    nt = min(Threads.nthreads(), n_work)
    chunks = _partition_range(1:n_work, nt)
    dim = size(evolving_dm, 1)

    # Per-task scratch: each needs rho_jump, jump_oft, sandwich_tmp.
    # If a pool is supplied (qf-po5), reuse it; otherwise allocate fresh
    # (backward-compat path used by run_thermalize / run_trajectory today).
    local_pool = if task_scratches === nothing
        [ThermalizeScratch(CT, dim) for _ in 1:length(chunks)]
    else
        @assert length(task_scratches) >= length(chunks)
        task_scratches
    end

    old_blas = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    try
        @sync for (idx, chunk) in enumerate(chunks)
            Threads.@spawn _accumulate_rho_jump_chunk_energy!(
                local_pool[idx], evolving_dm, jump, hamiltonian,
                config, precomputed_data, half_indices[chunk];
                base_prefactor=base_prefactor, inv_4sigma2=inv_4sigma2)
        end
    finally
        BLAS.set_num_threads(old_blas)
    end

    # Reduce: sum per-task rho_jump into scratch.rho_jump
    fill!(scratch.rho_jump, 0)
    for idx in 1:length(chunks)
        scratch.rho_jump .+= local_pool[idx].rho_jump
    end

    return nothing
end

function _accumulate_rho_jump_chunk_energy!(
    scratch::ThermalizeScratch{CT},
    evolving_dm::Matrix{CT},
    jump::JumpOp,
    hamiltonian::HamHam,
    config::Config{Thermalize, EnergyDomain},
    precomputed_data,
    label_indices::AbstractVector{Int};
    base_prefactor::Real,
    inv_4sigma2::Real,
) where {CT<:Complex}
    (; transition, energy_labels) = precomputed_data
    fill!(scratch.rho_jump, 0)

    @inbounds for li in label_indices
        w_raw = energy_labels[li]
        w = abs(w_raw)

        oft!(scratch.jump_oft, jump.in_eigenbasis, hamiltonian.bohr_freqs, w, inv_4sigma2)

        rate2_pos = base_prefactor * transition(w)
        mul!(scratch.sandwich_tmp, evolving_dm, scratch.jump_oft')
        mul!(scratch.rho_jump, scratch.jump_oft, scratch.sandwich_tmp, config.delta * rate2_pos, 1.0)

        if jump.hermitian && w > 1e-12
            rate2_neg = base_prefactor * transition(-w)
            mul!(scratch.sandwich_tmp, evolving_dm, scratch.jump_oft)
            mul!(scratch.rho_jump, scratch.jump_oft', scratch.sandwich_tmp, config.delta * rate2_neg, 1.0)
        end
    end

    return nothing
end

# --- TimeDomain/TrotterDomain threaded variant ---

function _accumulate_rho_jump_threaded_timetrot!(
    scratch::ThermalizeScratch{CT},
    evolving_dm::Matrix{CT},
    jump::JumpOp,
    ham_or_trott,
    config::Config{Thermalize, D},
    precomputed_data;
    jump_weight_scaling::Real,
    task_scratches::Union{Nothing, Vector{ThermalizeScratch{CT}}}=nothing,
) where {CT<:Complex, D<:Union{TimeDomain, TrotterDomain}}
    (; transition, energy_labels, oft_nufft_prefactors) = precomputed_data
    base_prefactor = precomputed_data.oft_domain_prefactor * jump_weight_scaling

    # For Hermitian jumps, pre-filter to half-grid indices for balanced partitioning
    if jump.hermitian
        half_indices = [i for i in eachindex(energy_labels) if energy_labels[i] <= 1e-12]
    else
        half_indices = collect(eachindex(energy_labels))
    end

    n_work = length(half_indices)
    nt = min(Threads.nthreads(), n_work)
    chunks = _partition_range(1:n_work, nt)
    dim = size(evolving_dm, 1)

    local_pool = if task_scratches === nothing
        [ThermalizeScratch(CT, dim) for _ in 1:length(chunks)]
    else
        @assert length(task_scratches) >= length(chunks)
        task_scratches
    end

    old_blas = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    try
        @sync for (idx, chunk) in enumerate(chunks)
            Threads.@spawn _accumulate_rho_jump_chunk_timetrot!(
                local_pool[idx], evolving_dm, jump,
                config, precomputed_data, half_indices[chunk];
                base_prefactor=base_prefactor)
        end
    finally
        BLAS.set_num_threads(old_blas)
    end

    # Reduce: sum per-task rho_jump into scratch.rho_jump
    fill!(scratch.rho_jump, 0)
    for idx in 1:length(chunks)
        scratch.rho_jump .+= local_pool[idx].rho_jump
    end

    return nothing
end

function _accumulate_rho_jump_chunk_timetrot!(
    scratch::ThermalizeScratch{CT},
    evolving_dm::Matrix{CT},
    jump::JumpOp,
    config::Config{Thermalize, D},
    precomputed_data,
    label_indices::AbstractVector{Int};
    base_prefactor::Real,
) where {CT<:Complex, D<:Union{TimeDomain, TrotterDomain}}
    (; transition, energy_labels, oft_nufft_prefactors) = precomputed_data
    fill!(scratch.rho_jump, 0)

    @inbounds for li in label_indices
        w_raw = energy_labels[li]
        w = abs(w_raw)

        nufft_prefactor_matrix = _prefactor_view(oft_nufft_prefactors, w)
        @. scratch.jump_oft = jump.in_eigenbasis * nufft_prefactor_matrix

        rate2_pos = base_prefactor * transition(w)
        mul!(scratch.sandwich_tmp, evolving_dm, scratch.jump_oft')
        mul!(scratch.rho_jump, scratch.jump_oft, scratch.sandwich_tmp, config.delta * rate2_pos, 1.0)

        if jump.hermitian && w > 1e-12
            rate2_neg = base_prefactor * transition(-w)
            mul!(scratch.sandwich_tmp, evolving_dm, scratch.jump_oft)
            mul!(scratch.rho_jump, scratch.jump_oft', scratch.sandwich_tmp, config.delta * rate2_neg, 1.0)
        end
    end

    return nothing
end

# --- BohrDomain threaded variant ---

function _accumulate_rho_jump_threaded_bohr!(
    scratch::ThermalizeScratch{CT},
    evolving_dm::Matrix{CT},
    jump::JumpOp,
    hamiltonian::HamHam,
    config::Config{Thermalize, BohrDomain},
    precomputed_data,
    bohr_keys::AbstractVector,
    bohr_is::Union{Nothing, Vector{Vector{Int}}},
    bohr_js::Union{Nothing, Vector{Vector{Int}}};
    jump_weight_scaling::Real,
    task_scratches::Union{Nothing, Vector{ThermalizeScratch{CT}}}=nothing,
) where {CT<:Complex}
    dim = size(evolving_dm, 1)
    (; alpha) = precomputed_data
    scaled_delta = config.delta * jump_weight_scaling

    n_keys = length(bohr_keys)
    nt = min(Threads.nthreads(), n_keys)
    chunks = _partition_range(1:n_keys, nt)

    local_pool = if task_scratches === nothing
        [ThermalizeScratch(CT, dim) for _ in 1:length(chunks)]
    else
        @assert length(task_scratches) >= length(chunks)
        task_scratches
    end

    old_blas = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    try
        @sync for (idx, chunk) in enumerate(chunks)
            Threads.@spawn _accumulate_rho_jump_chunk_bohr!(
                local_pool[idx], evolving_dm, jump, hamiltonian,
                precomputed_data, bohr_keys, bohr_is, bohr_js, chunk;
                scaled_delta=scaled_delta)
        end
    finally
        BLAS.set_num_threads(old_blas)
    end

    # Reduce: sum per-task rho_jump into scratch.rho_jump
    fill!(scratch.rho_jump, 0)
    for idx in 1:length(chunks)
        scratch.rho_jump .+= local_pool[idx].rho_jump
    end

    return nothing
end

function _accumulate_rho_jump_chunk_bohr!(
    scratch::ThermalizeScratch{CT},
    evolving_dm::Matrix{CT},
    jump::JumpOp,
    hamiltonian::HamHam,
    precomputed_data,
    bohr_keys::AbstractVector,
    bohr_is::Union{Nothing, Vector{Vector{Int}}},
    bohr_js::Union{Nothing, Vector{Vector{Int}}},
    key_indices::UnitRange{Int};
    scaled_delta::Real,
) where {CT<:Complex}
    dim = size(evolving_dm, 1)
    (; alpha) = precomputed_data
    fill!(scratch.rho_jump, 0)

    @inbounds for k in key_indices
        nu_2 = bohr_keys[k]

        # B_{v2} = sum_{v1} alpha(v1, v2) * A
        @. scratch.jump_oft = alpha(hamiltonian.bohr_freqs, nu_2) * jump.in_eigenbasis

        # sandwich_tmp := rho A_{v2}dag
        fill!(scratch.sandwich_tmp, 0)
        if bohr_is !== nothing
            is = bohr_is[k]
            js = bohr_js[k]
            @inbounds for t in eachindex(is)
                i = is[t]
                j = js[t]
                v = conj(jump.in_eigenbasis[i, j])
                @inbounds for p in 1:dim
                    scratch.sandwich_tmp[p, i] += evolving_dm[p, j] * v
                end
            end
        else
            indices = hamiltonian.bohr_dict[nu_2]
            @inbounds for idx in indices
                i = idx[1]
                j = idx[2]
                v = conj(jump.in_eigenbasis[i, j])
                @inbounds for p in 1:dim
                    scratch.sandwich_tmp[p, i] += evolving_dm[p, j] * v
                end
            end
        end

        # rho_jump += scaled_delta * B_{v2} * (rho A_{v2}dag)
        mul!(scratch.rho_jump, scratch.jump_oft, scratch.sandwich_tmp, scaled_delta, 1.0)
    end

    return nothing
end
