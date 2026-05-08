"""
    Workspace(config::Config{Lindbladian}, hamiltonian, jumps; trotter=nothing)

Construct a `Workspace{KrylovSpectrum}` pre-allocating all scratch matrices for the given
(config, hamiltonian) pair. Mirrors `construct_lindbladian` setup in `furnace.jl`.

Returns `Workspace{KrylovSpectrum,D,C,T}`.
"""
function Workspace(
    config::Config{Lindbladian},
    hamiltonian::HamHam,
    jumps::Vector{JumpOp};
    trotter::Union{TrottTrott, Nothing}=nothing,
)
    # Determine ham_or_trott (mirrors construct_lindbladian in furnace.jl)
    ham_or_trott = if config.domain isa TrotterDomain
        trotter === nothing && error("A Trotter object must be provided for the TrotterDomain")
        trotter
    else
        hamiltonian
    end

    # Precompute domain-specific data (transition, energy_labels, gamma_norm_factor, ...)
    precomputed_data = _precompute_data(config, ham_or_trott)

    # Precompute coherent B_total (returns nothing for GNS / with_coherent=false)
    B_total = _precompute_coherent_B(jumps, ham_or_trott, config, precomputed_data)

    # Determine dimensions and element type
    dim = size(hamiltonian.data, 1)
    T = eltype(hamiltonian.eigvals)
    CT = Complex{T}

    # Extract concrete-typed eigenbasis matrices and hermitian flags
    jump_eigenbases = [Matrix{CT}(j.in_eigenbasis) for j in jumps]
    jump_hermitian  = [j.hermitian for j in jumps]

    # Allocate KrylovScratch (no channel_rho_jump for Lindbladian)
    sc = KrylovScratch(CT, dim; with_channel_rho_jump=false)

    # Precompute G_left/G_right for optimized Lindbladian matvec (Phase 32)
    R_total = zeros(CT, dim, dim)
    _accumulate_R_total!(R_total, jump_eigenbases, jump_hermitian,
                         precomputed_data, config, ham_or_trott)
    hermitianize!(R_total)

    if B_total !== nothing
        B_T = Matrix{CT}(transpose(B_total))
        G_left  = Matrix{CT}(1im .* B_T .- 0.5 .* R_total)
        G_right = Matrix{CT}(-1im .* B_T .- 0.5 .* R_total)
    else
        G_left  = Matrix{CT}(-0.5 .* R_total)
        G_right = Matrix{CT}(-0.5 .* R_total)
    end

    G_left_adj  = G_right
    G_right_adj = G_left

    # Absorb precomputed_data fields into flat workspace fields
    pd_transition = hasproperty(precomputed_data, :transition) ? precomputed_data.transition : nothing
    pd_gnf = hasproperty(precomputed_data, :gamma_norm_factor) ? precomputed_data.gamma_norm_factor : nothing
    pd_el = hasproperty(precomputed_data, :energy_labels) ? precomputed_data.energy_labels : nothing
    pd_odp = hasproperty(precomputed_data, :oft_domain_prefactor) ? precomputed_data.oft_domain_prefactor : nothing
    pd_nufft = hasproperty(precomputed_data, :oft_nufft_prefactors) ? precomputed_data.oft_nufft_prefactors : nothing
    pd_oft_pre_energy = hasproperty(precomputed_data, :oft_prefactors_energy) ? precomputed_data.oft_prefactors_energy : nothing
    pd_alpha = hasproperty(precomputed_data, :alpha) ? precomputed_data.alpha : nothing
    pd_bkeys = hasproperty(precomputed_data, :bohr_keys) ? precomputed_data.bohr_keys : nothing
    pd_bis = hasproperty(precomputed_data, :bohr_is) ? precomputed_data.bohr_is : nothing
    pd_bjs = hasproperty(precomputed_data, :bohr_js) ? precomputed_data.bohr_js : nothing
    pd_bminus = hasproperty(precomputed_data, :b_minus) ? precomputed_data.b_minus : nothing
    pd_bplus = hasproperty(precomputed_data, :b_plus) ? precomputed_data.b_plus : nothing

    # qf-in3.4: pre-build the (jump_idx, label_idx) work list for the threaded
    # ω-loop dispatch in apply_lindbladian! (EnergyDomain / TimeDomain /
    # TrotterDomain). BohrDomain leaves the list empty.
    if pd_el !== nothing
        _populate_lindblad_work_list!(sc.work_list, jump_hermitian, pd_el)
    end

    D = typeof(config.domain)
    C = typeof(config.construction)

    return Workspace{KrylovSpectrum, D, C, T}(
        jump_eigenbases, jump_hermitian, jumps, B_total,
        nothing,  # dll_lindblads (CKG/GNS path; populated by DLL specialised constructor)
        G_left, G_right, G_left_adj, G_right_adj,
        nothing, nothing, nothing, nothing, nothing,  # channel fields
        pd_transition, pd_gnf, pd_el, pd_odp, pd_nufft, pd_oft_pre_energy,
        pd_alpha, pd_bkeys, pd_bis, pd_bjs, pd_bminus, pd_bplus,
        nothing,  # coherent_unitaries
        nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing,  # trajectory fields
        nothing,  # jump_selection (KrylovSpectrum path: no per-step jump selection)
        nothing,  # Id
        sc,
    )
end

# ---------------------------------------------------------------------------
# R_total accumulation helpers (physics convention: R = sum rate^2 * L' * L)
# ---------------------------------------------------------------------------

"""
    _accumulate_R_total!(R, ws, config, hamiltonian) -> nothing

Accumulate R_total = sum over all jumps and frequencies of rate^2 * (L' * L)
in physics convention. Used at workspace construction time (not per-matvec).
"""
function _accumulate_R_total!(
    R::Matrix{T},
    ws_eigenbases::Vector{Matrix{T}},
    ws_hermitian::Vector{Bool},
    precomputed_data,
    config::Config{<:Any, EnergyDomain},
    hamiltonian::HamHam,
) where {T<:Complex}
    (; transition, gamma_norm_factor, energy_labels, oft_prefactors_energy) = precomputed_data
    prefactor = precomputed_data.oft_domain_prefactor * gamma_norm_factor

    n_labels = length(energy_labels)
    n_jumps  = length(ws_eigenbases)
    if Threads.nthreads() > 1 && n_jumps * n_labels >= OMEGA_THREAD_THRESHOLD
        return _accumulate_R_total_threaded_energy!(
            R, ws_eigenbases, ws_hermitian, oft_prefactors_energy, energy_labels,
            transition, prefactor)
    end

    dim = size(R, 1)
    jump_oft = zeros(T, dim, dim)
    LdagL = zeros(T, dim, dim)

    for (k, eigenbasis) in enumerate(ws_eigenbases)
        is_herm = ws_hermitian[k]
        if is_herm
            for w_raw in energy_labels
                w_raw > 1e-12 && continue
                w = abs(w_raw)
                pref_view = _prefactor_view(oft_prefactors_energy, w)
                @. jump_oft = eigenbasis * pref_view
                rate2 = prefactor * transition(w)
                mul!(LdagL, jump_oft', jump_oft)
                @. R += rate2 * LdagL
                if w > 1e-12
                    rate2_neg = prefactor * transition(-w)
                    mul!(LdagL, jump_oft, jump_oft')
                    @. R += rate2_neg * LdagL
                end
            end
        else
            for w in energy_labels
                pref_view = _prefactor_view(oft_prefactors_energy, w)
                @. jump_oft = eigenbasis * pref_view
                rate2 = prefactor * transition(w)
                mul!(LdagL, jump_oft', jump_oft)
                @. R += rate2 * LdagL
            end
        end
    end
    return nothing
end

function _accumulate_R_total!(
    R::Matrix{T},
    ws_eigenbases::Vector{Matrix{T}},
    ws_hermitian::Vector{Bool},
    precomputed_data,
    config::Config{<:Any, D},
    ham_or_trott::Union{HamHam, TrottTrott},
) where {T<:Complex, D<:Union{TimeDomain, TrotterDomain}}
    (; transition, gamma_norm_factor, energy_labels, oft_nufft_prefactors) = precomputed_data
    prefactor = precomputed_data.oft_domain_prefactor * gamma_norm_factor

    n_labels = length(energy_labels)
    n_jumps  = length(ws_eigenbases)
    if Threads.nthreads() > 1 && n_jumps * n_labels >= OMEGA_THREAD_THRESHOLD
        return _accumulate_R_total_threaded_timetrot!(
            R, ws_eigenbases, ws_hermitian, oft_nufft_prefactors,
            energy_labels, transition, prefactor)
    end

    dim = size(R, 1)
    jump_oft = zeros(T, dim, dim)
    LdagL = zeros(T, dim, dim)

    for (k, eigenbasis) in enumerate(ws_eigenbases)
        is_herm = ws_hermitian[k]
        if is_herm
            for w_raw in energy_labels
                w_raw > 1e-12 && continue
                w = abs(w_raw)
                nufft_pf = _prefactor_view(oft_nufft_prefactors, w)
                @. jump_oft = eigenbasis * nufft_pf
                rate2 = prefactor * transition(w)
                mul!(LdagL, jump_oft', jump_oft)
                @. R += rate2 * LdagL
                if w > 1e-12
                    rate2_neg = prefactor * transition(-w)
                    mul!(LdagL, jump_oft, jump_oft')
                    @. R += rate2_neg * LdagL
                end
            end
        else
            for w in energy_labels
                nufft_pf = _prefactor_view(oft_nufft_prefactors, w)
                @. jump_oft = eigenbasis * nufft_pf
                rate2 = prefactor * transition(w)
                mul!(LdagL, jump_oft', jump_oft)
                @. R += rate2 * LdagL
            end
        end
    end
    return nothing
end

function _accumulate_R_total!(
    R::Matrix{T},
    ws_eigenbases::Vector{Matrix{T}},
    ws_hermitian::Vector{Bool},
    precomputed_data,
    config::Config{<:Any, BohrDomain},
    hamiltonian::HamHam,
) where {T<:Complex}
    (; alpha, gamma_norm_factor) = precomputed_data
    bohr_keys = collect(keys(hamiltonian.bohr_dict))
    n_jumps = length(ws_eigenbases)
    n_keys  = length(bohr_keys)
    if Threads.nthreads() > 1 && n_jumps * n_keys >= OMEGA_THREAD_THRESHOLD
        return _accumulate_R_total_threaded_bohr!(
            R, ws_eigenbases, alpha, gamma_norm_factor,
            hamiltonian.bohr_freqs, hamiltonian.bohr_dict, bohr_keys)
    end

    dim = size(R, 1)
    jump_oft = zeros(T, dim, dim)
    A_nu2_dag = zeros(T, dim, dim)

    for (k, eigenbasis) in enumerate(ws_eigenbases)
        for nu_2 in bohr_keys
            @. jump_oft = alpha(hamiltonian.bohr_freqs, nu_2) * eigenbasis

            fill!(A_nu2_dag, 0)
            indices = hamiltonian.bohr_dict[nu_2]
            @inbounds for idx in indices
                i = idx[1]; j = idx[2]
                A_nu2_dag[j, i] = conj(eigenbasis[i, j])
            end

            mul!(R, A_nu2_dag, jump_oft, gamma_norm_factor, 1.0)
        end
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Threaded R_total accumulation (qf-6af) — mirrors the matvec ω-loop pattern
# from `src/krylov_matvec.jl`, but for the construction-time additive
# accumulator. Allocates per-task buffers locally (one-time at Workspace
# construction; trivial overhead vs. the BLAS work it covers). Each task
# accumulates a private `R_partial`; final reduction sums into the caller's
# `R` after `@sync`.
# ---------------------------------------------------------------------------

# Build a flat (k, li) work list with the same Hermitian-fold convention as
# `_populate_lindblad_work_list!` in krylov_matvec.jl. Returns a fresh vector
# (no scratch reuse — this fires only at Workspace construction).
function _build_R_total_work_list(
    jump_hermitian::Vector{Bool},
    energy_labels::AbstractVector{<:Real},
)
    n_jumps  = length(jump_hermitian)
    n_labels = length(energy_labels)
    work = Tuple{Int, Int}[]
    sizehint!(work, n_jumps * n_labels)
    @inbounds for k in 1:n_jumps
        is_herm = jump_hermitian[k]
        for li in 1:n_labels
            if is_herm && energy_labels[li] > 1e-12
                continue
            end
            push!(work, (k, li))
        end
    end
    return work
end

function _accumulate_R_total_threaded_energy!(
    R::Matrix{T},
    ws_eigenbases::Vector{Matrix{T}},
    ws_hermitian::Vector{Bool},
    oft_prefactors_energy::EnergyDomainPrefactors{Float64, Array{Float64, 3}},
    energy_labels::AbstractVector{Float64},
    transition,
    prefactor::Float64,
) where {T<:Complex}
    work = _build_R_total_work_list(ws_hermitian, energy_labels)
    n_work = length(work)
    n_work == 0 && return nothing

    dim = size(R, 1)
    nt = min(Threads.nthreads(), n_work)
    chunks = _partition_range(1:n_work, nt)
    n_chunks = length(chunks)

    R_partials = [zeros(T, dim, dim) for _ in 1:n_chunks]
    jump_ofts  = [Matrix{T}(undef, dim, dim) for _ in 1:n_chunks]
    LdagLs     = [Matrix{T}(undef, dim, dim) for _ in 1:n_chunks]

    old_blas = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    try
        @sync for (idx, chunk) in enumerate(chunks)
            Threads.@spawn _accumulate_R_total_chunk_energy!(
                R_partials[idx], jump_ofts[idx], LdagLs[idx],
                ws_eigenbases, ws_hermitian, oft_prefactors_energy, energy_labels,
                work, chunk, transition, prefactor)
        end
    finally
        BLAS.set_num_threads(old_blas)
    end

    @inbounds for idx in 1:n_chunks
        R .+= R_partials[idx]
    end
    return nothing
end

function _accumulate_R_total_chunk_energy!(
    R_partial::Matrix{T},
    jump_oft::Matrix{T},
    LdagL::Matrix{T},
    ws_eigenbases::Vector{Matrix{T}},
    ws_hermitian::Vector{Bool},
    oft_prefactors_energy::EnergyDomainPrefactors{Float64, Array{Float64, 3}},
    energy_labels::AbstractVector{Float64},
    work::Vector{Tuple{Int, Int}},
    chunk::UnitRange{Int},
    transition,
    prefactor::Float64,
) where {T<:Complex}
    @inbounds for w_idx in chunk
        (k, li) = work[w_idx]
        eigenbasis = ws_eigenbases[k]
        is_herm    = ws_hermitian[k]

        w_raw = energy_labels[li]
        w = is_herm ? abs(w_raw) : w_raw

        pref_view = _prefactor_view(oft_prefactors_energy, w)
        @. jump_oft = eigenbasis * pref_view
        rate2 = prefactor * transition(w)
        mul!(LdagL, jump_oft', jump_oft)
        @. R_partial += rate2 * LdagL

        if is_herm && w > 1e-12
            rate2_neg = prefactor * transition(-w)
            mul!(LdagL, jump_oft, jump_oft')
            @. R_partial += rate2_neg * LdagL
        end
    end
    return nothing
end

function _accumulate_R_total_threaded_timetrot!(
    R::Matrix{T},
    ws_eigenbases::Vector{Matrix{T}},
    ws_hermitian::Vector{Bool},
    oft_nufft_prefactors,
    energy_labels::AbstractVector{Float64},
    transition,
    prefactor::Float64,
) where {T<:Complex}
    work = _build_R_total_work_list(ws_hermitian, energy_labels)
    n_work = length(work)
    n_work == 0 && return nothing

    dim = size(R, 1)
    nt = min(Threads.nthreads(), n_work)
    chunks = _partition_range(1:n_work, nt)
    n_chunks = length(chunks)

    R_partials = [zeros(T, dim, dim) for _ in 1:n_chunks]
    jump_ofts  = [Matrix{T}(undef, dim, dim) for _ in 1:n_chunks]
    LdagLs     = [Matrix{T}(undef, dim, dim) for _ in 1:n_chunks]

    old_blas = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    try
        @sync for (idx, chunk) in enumerate(chunks)
            Threads.@spawn _accumulate_R_total_chunk_timetrot!(
                R_partials[idx], jump_ofts[idx], LdagLs[idx],
                ws_eigenbases, ws_hermitian, oft_nufft_prefactors,
                energy_labels, work, chunk, transition, prefactor)
        end
    finally
        BLAS.set_num_threads(old_blas)
    end

    @inbounds for idx in 1:n_chunks
        R .+= R_partials[idx]
    end
    return nothing
end

function _accumulate_R_total_chunk_timetrot!(
    R_partial::Matrix{T},
    jump_oft::Matrix{T},
    LdagL::Matrix{T},
    ws_eigenbases::Vector{Matrix{T}},
    ws_hermitian::Vector{Bool},
    oft_nufft_prefactors,
    energy_labels::AbstractVector{Float64},
    work::Vector{Tuple{Int, Int}},
    chunk::UnitRange{Int},
    transition,
    prefactor::Float64,
) where {T<:Complex}
    @inbounds for w_idx in chunk
        (k, li) = work[w_idx]
        eigenbasis = ws_eigenbases[k]
        is_herm    = ws_hermitian[k]

        w_raw = energy_labels[li]
        # Hermitian fold: only `w_raw <= 1e-12` queued; index NUFFT slice via
        # the |w_raw| key. Non-Hermitian: index by label position `li`.
        w = is_herm ? abs(w_raw) : w_raw
        nufft_pf = is_herm ?
            _prefactor_view(oft_nufft_prefactors, w) :
            (@view oft_nufft_prefactors.data[:, :, li])
        @. jump_oft = eigenbasis * nufft_pf

        rate2 = prefactor * transition(w)
        mul!(LdagL, jump_oft', jump_oft)
        @. R_partial += rate2 * LdagL

        if is_herm && w > 1e-12
            rate2_neg = prefactor * transition(-w)
            mul!(LdagL, jump_oft, jump_oft')
            @. R_partial += rate2_neg * LdagL
        end
    end
    return nothing
end

function _accumulate_R_total_threaded_bohr!(
    R::Matrix{T},
    ws_eigenbases::Vector{Matrix{T}},
    alpha,
    gamma_norm_factor::Real,
    bohr_freqs::AbstractMatrix{<:Real},
    bohr_dict,
    bohr_keys::Vector,
) where {T<:Complex}
    n_jumps = length(ws_eigenbases)
    n_keys  = length(bohr_keys)
    n_work  = n_jumps * n_keys
    n_work == 0 && return nothing

    dim = size(R, 1)
    nt = min(Threads.nthreads(), n_work)
    chunks = _partition_range(1:n_work, nt)
    n_chunks = length(chunks)

    R_partials  = [zeros(T, dim, dim) for _ in 1:n_chunks]
    jump_ofts   = [Matrix{T}(undef, dim, dim) for _ in 1:n_chunks]
    A_nu2_dags  = [zeros(T, dim, dim) for _ in 1:n_chunks]

    old_blas = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    try
        @sync for (idx, chunk) in enumerate(chunks)
            Threads.@spawn _accumulate_R_total_chunk_bohr!(
                R_partials[idx], jump_ofts[idx], A_nu2_dags[idx],
                ws_eigenbases, alpha, gamma_norm_factor,
                bohr_freqs, bohr_dict, bohr_keys, n_keys, chunk)
        end
    finally
        BLAS.set_num_threads(old_blas)
    end

    @inbounds for idx in 1:n_chunks
        R .+= R_partials[idx]
    end
    return nothing
end

function _accumulate_R_total_chunk_bohr!(
    R_partial::Matrix{T},
    jump_oft::Matrix{T},
    A_nu2_dag::Matrix{T},
    ws_eigenbases::Vector{Matrix{T}},
    alpha,
    gamma_norm_factor::Real,
    bohr_freqs::AbstractMatrix{<:Real},
    bohr_dict,
    bohr_keys::Vector,
    n_keys::Int,
    chunk::UnitRange{Int},
) where {T<:Complex}
    @inbounds for w_idx in chunk
        # Linear (k, key_idx) decoding: outer over jumps, inner over keys.
        k       = ((w_idx - 1) ÷ n_keys) + 1
        key_idx = ((w_idx - 1) % n_keys) + 1
        eigenbasis = ws_eigenbases[k]
        nu_2 = bohr_keys[key_idx]

        @. jump_oft = alpha(bohr_freqs, nu_2) * eigenbasis

        fill!(A_nu2_dag, 0)
        indices = bohr_dict[nu_2]
        @inbounds for idx in indices
            i = idx[1]; j = idx[2]
            A_nu2_dag[j, i] = conj(eigenbasis[i, j])
        end

        mul!(R_partial, A_nu2_dag, jump_oft, gamma_norm_factor, 1.0)
    end
    return nothing
end

"""
    _accumulate_R_total_dll!(R, dll_lindblads_out, jumps, precomputed_data, config, hamiltonian) -> nothing

DLL Bohr-domain accumulation of `R = Σ_a L_a' L_a` plus capture of the per-jump
`L_a = dll_lindblad_op_bohr(jump, hamiltonian, filter)` matrices into
`dll_lindblads_out` (mutated in-place). Distinct from the
`_accumulate_R_total!` family because (i) the dispatch key here is
`Config{<:Any, BohrDomain, DLL}` and (ii) we need the per-jump matrices
themselves for the matrix-free dispatch in `apply_lindbladian!`.
"""
function _accumulate_R_total_dll!(
    R::Matrix{T},
    dll_lindblads_out::Vector{Matrix{T}},
    jumps::Vector{JumpOp},
    precomputed_data,
    config::Config{<:Any, BohrDomain, DLL},
    hamiltonian::HamHam,
) where {T<:Complex}
    (; filter) = precomputed_data
    n_jumps = length(jumps)

    # Threaded path (qf-6af.4): each task builds its (per-channel) DLL Lindblad
    # operator(s) and accumulates a private R_partial. Reduce after `@sync`.
    # `dll_lindblads_out` retains deterministic per-jump order.
    if Threads.nthreads() > 1 && n_jumps >= 2
        nt = min(Threads.nthreads(), n_jumps)
        chunks = _partition_range(1:n_jumps, nt)
        n_chunks = length(chunks)
        dim = size(R, 1)

        # Each task accumulates into its own R_partial and a per-jump-vector
        # of operator lists, sized to the original jump count for ordered merge.
        R_partials = [zeros(T, dim, dim) for _ in 1:n_chunks]
        per_jump_ops = Vector{Vector{Matrix{T}}}(undef, n_jumps)

        old_blas = BLAS.get_num_threads()
        BLAS.set_num_threads(1)
        try
            @sync for (idx, chunk) in enumerate(chunks)
                Threads.@spawn _accumulate_R_total_dll_chunk!(
                    R_partials[idx], per_jump_ops, jumps, hamiltonian,
                    filter, chunk)
            end
        finally
            BLAS.set_num_threads(old_blas)
        end

        @inbounds for idx in 1:n_chunks
            R .+= R_partials[idx]
        end
        @inbounds for k in 1:n_jumps
            for L_a in per_jump_ops[k]
                push!(dll_lindblads_out, L_a)
            end
        end
        return nothing
    end

    for jump in jumps
        L_or_Ls = dll_lindblad_op_bohr(jump, hamiltonian, filter)
        # Single-channel filters return a Matrix; multi-channel filters
        # (qf-7go.1) return a Vector{Matrix} of length k. Flatten into the
        # output `dll_lindblads_out` — the matrix-free hot path
        # `apply_lindbladian!` iterates `for L_a in dll_lindblads` and the
        # dissipator `Σ_a L_a ρ L_a†` is a flat sum over all per-channel
        # operators (no cross terms in the multi-channel α).
        if L_or_Ls isa AbstractMatrix
            L_a = Matrix{T}(L_or_Ls)
            push!(dll_lindblads_out, L_a)
            mul!(R, L_a', L_a, 1.0, 1.0)
        else
            for L_one in L_or_Ls
                L_a = Matrix{T}(L_one)
                push!(dll_lindblads_out, L_a)
                mul!(R, L_a', L_a, 1.0, 1.0)
            end
        end
    end
    return nothing
end

function _accumulate_R_total_dll_chunk!(
    R_partial::Matrix{T},
    per_jump_ops::Vector{Vector{Matrix{T}}},
    jumps::Vector{JumpOp},
    hamiltonian::HamHam,
    filter,
    chunk::UnitRange{Int},
) where {T<:Complex}
    @inbounds for k in chunk
        L_or_Ls = dll_lindblad_op_bohr(jumps[k], hamiltonian, filter)
        ops = Vector{Matrix{T}}()
        if L_or_Ls isa AbstractMatrix
            L_a = Matrix{T}(L_or_Ls)
            push!(ops, L_a)
            mul!(R_partial, L_a', L_a, 1.0, 1.0)
        else
            for L_one in L_or_Ls
                L_a = Matrix{T}(L_one)
                push!(ops, L_a)
                mul!(R_partial, L_a', L_a, 1.0, 1.0)
            end
        end
        per_jump_ops[k] = ops
    end
    return nothing
end

# ---------------------------------------------------------------------------
# DLL BohrDomain Lindbladian workspace constructor
# ---------------------------------------------------------------------------

"""
    Workspace(config::Config{Lindbladian, BohrDomain, DLL}, hamiltonian, jumps; trotter=nothing)

DLL-specialised `Workspace{KrylovSpectrum, BohrDomain, DLL, T}` for Bohr-domain
matrix-free `apply_lindbladian!` (qf-lkb.9). Precomputes:

- `dll_lindblads = [dll_lindblad_op_bohr(jump, hamiltonian, filter) for jump in jumps]`
  (Ding–Li–Lin 2024 Eq. 3.4 first form),
- `R_total = Σ_a L_a' L_a`, hermitised post-accumulation,
- `G = dll_coherent_op_bohr(jumps, hamiltonian, filter, β)` (Eq. 3.5 + 3.7),
  stored in the existing `B_total` slot.

Sign convention for `G_left`/`G_right` matches the CKG path
(`_vectorize_liouvillian_coherent!`): with `B_T = transpose(G)`,

    G_left  = +1im · B_T − 0.5 · R_total
    G_right = −1im · B_T − 0.5 · R_total

so that `apply_lindbladian!`'s gemm pattern `G_left·ρ + ρ·G_right` realises
`+i [B_T, ρ] − 0.5 {R, ρ}`, byte-for-byte agreement with
`Matrix(construct_lindbladian(jumps, config, hamiltonian)) · vec(ρ)`.
"""
function Workspace(
    config::Config{Lindbladian, BohrDomain, DLL},
    hamiltonian::HamHam,
    jumps::Vector{JumpOp};
    trotter::Union{TrottTrott, Nothing}=nothing,
)
    @assert trotter === nothing  "DLL BohrDomain does not use Trotter"

    precomputed_data = _precompute_data(config, hamiltonian)
    (; filter) = precomputed_data

    dim = size(hamiltonian.data, 1)
    T = eltype(hamiltonian.eigvals)
    CT = Complex{T}

    # Per-jump DLL Lindblad operators + accumulated R_total. For
    # multi-channel filters (qf-7go.1) this stores k×|jumps| flattened
    # operators (one per channel per coupling).
    n_channels = filter isa DLLMultiChannelFilter ? length(filter.channels) : 1
    dll_lindblads = Vector{Matrix{CT}}()
    sizehint!(dll_lindblads, length(jumps) * n_channels)
    R_total = zeros(CT, dim, dim)
    _accumulate_R_total_dll!(R_total, dll_lindblads, jumps,
                              precomputed_data, config, hamiltonian)
    hermitianize!(R_total)

    # Coherent G via dll_coherent_op_bohr (Hermitian by Theorem 10).
    G = Matrix{CT}(dll_coherent_op_bohr(jumps, hamiltonian, filter, config.beta))

    # Match CKG sign convention: dense uses kron(B, I)·(-1im) + kron(I, B^T)·(+1im) = +i[B^T, ρ].
    B_T = Matrix{CT}(transpose(G))
    G_left  = Matrix{CT}( 1im .* B_T .- 0.5 .* R_total)
    G_right = Matrix{CT}(-1im .* B_T .- 0.5 .* R_total)
    G_left_adj  = G_right
    G_right_adj = G_left

    jump_eigenbases = [Matrix{CT}(j.in_eigenbasis) for j in jumps]
    jump_hermitian  = [j.hermitian for j in jumps]

    sc = KrylovScratch(CT, dim; with_channel_rho_jump=false)

    return Workspace{KrylovSpectrum, BohrDomain, DLL, T}(
        jump_eigenbases, jump_hermitian, jumps, G,  # B_total slot stores G
        dll_lindblads,
        G_left, G_right, G_left_adj, G_right_adj,
        nothing, nothing, nothing, nothing, nothing,  # channel fields
        nothing, nothing, nothing, nothing, nothing, nothing,  # transition/gnf/energy_labels/odp/nufft/oft_prefactors_energy
        nothing, nothing, nothing, nothing, nothing, nothing,  # bohr_alpha/bohr_keys/bohr_is/bohr_js/b_minus/b_plus
        nothing,  # coherent_unitaries
        nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing,  # trajectory fields
        nothing,  # jump_selection
        nothing,  # Id
        sc,
    )
end

# ---------------------------------------------------------------------------
# Config{Thermalize} workspace constructor
# ---------------------------------------------------------------------------

"""
    Workspace(config::Config{Thermalize}, hamiltonian, jumps; trotter=nothing)

Construct a `Workspace{KrylovSpectrum}` with precomputed CPTP channel matrices for
the faithful Chen channel (Eq. 3.2).

Returns `Workspace{KrylovSpectrum,D,C,T}`.
"""
function Workspace(
    config::Config{Thermalize},
    hamiltonian::HamHam,
    jumps::Vector{JumpOp};
    trotter::Union{TrottTrott, Nothing}=nothing,
)
    # Determine ham_or_trott
    ham_or_trott = if config.domain isa TrotterDomain
        trotter === nothing && error("A Trotter object must be provided for the TrotterDomain")
        trotter
    else
        hamiltonian
    end

    # Precompute domain-specific data
    precomputed_data = _precompute_data(config, ham_or_trott)

    # Precompute coherent B_total (returns nothing for GNS / with_coherent=false)
    B_total = _precompute_coherent_B(jumps, ham_or_trott, config, precomputed_data)

    # Determine dimensions and element type
    dim = size(hamiltonian.data, 1)
    T = eltype(hamiltonian.eigvals)
    CT = Complex{T}

    # Extract concrete-typed eigenbasis matrices and hermitian flags
    jump_eigenbases = [Matrix{CT}(j.in_eigenbasis) for j in jumps]
    jump_hermitian  = [j.hermitian for j in jumps]

    # Allocate KrylovScratch (with channel_rho_jump for Thermalize)
    sc = KrylovScratch(CT, dim; with_channel_rho_jump=true)

    # --- Precompute channel matrices ---
    delta = config.delta

    # 1. Compute R_total (physics convention)
    R_total = zeros(CT, dim, dim)
    _accumulate_R_total!(R_total, jump_eigenbases, jump_hermitian,
                         precomputed_data, config, ham_or_trott)
    hermitianize!(R_total)

    # Precompute G_left/G_right
    if B_total !== nothing
        B_T = Matrix{CT}(transpose(B_total))
        G_left  = Matrix{CT}(1im .* B_T .- 0.5 .* R_total)
        G_right = Matrix{CT}(-1im .* B_T .- 0.5 .* R_total)
    else
        G_left  = Matrix{CT}(-0.5 .* R_total)
        G_right = Matrix{CT}(-0.5 .* R_total)
    end

    G_left_adj  = G_right
    G_right_adj = G_left

    # 2. Compute K0, U_residual (Chen Eq. 3.2)
    channel = _build_cptp_channel(R_total, delta)

    # 3. Compute U_coherent = exp(-i*delta*B_total) if coherent
    U_coherent = if B_total !== nothing
        Matrix{CT}(exp(-1im * delta * Hermitian(B_total)))
    else
        nothing
    end

    # Absorb precomputed_data fields
    pd_transition = hasproperty(precomputed_data, :transition) ? precomputed_data.transition : nothing
    pd_gnf = hasproperty(precomputed_data, :gamma_norm_factor) ? precomputed_data.gamma_norm_factor : nothing
    pd_el = hasproperty(precomputed_data, :energy_labels) ? precomputed_data.energy_labels : nothing
    pd_odp = hasproperty(precomputed_data, :oft_domain_prefactor) ? precomputed_data.oft_domain_prefactor : nothing
    pd_nufft = hasproperty(precomputed_data, :oft_nufft_prefactors) ? precomputed_data.oft_nufft_prefactors : nothing
    pd_oft_pre_energy = hasproperty(precomputed_data, :oft_prefactors_energy) ? precomputed_data.oft_prefactors_energy : nothing
    pd_alpha = hasproperty(precomputed_data, :alpha) ? precomputed_data.alpha : nothing
    pd_bkeys = hasproperty(precomputed_data, :bohr_keys) ? precomputed_data.bohr_keys : nothing
    pd_bis = hasproperty(precomputed_data, :bohr_is) ? precomputed_data.bohr_is : nothing
    pd_bjs = hasproperty(precomputed_data, :bohr_js) ? precomputed_data.bohr_js : nothing
    pd_bminus = hasproperty(precomputed_data, :b_minus) ? precomputed_data.b_minus : nothing
    pd_bplus = hasproperty(precomputed_data, :b_plus) ? precomputed_data.b_plus : nothing

    # qf-in3.4: pre-build the (jump_idx, label_idx) work list for the threaded
    # ω-loop dispatch in apply_delta_channel! (EnergyDomain / TimeDomain /
    # TrotterDomain). BohrDomain leaves the list empty.
    if pd_el !== nothing
        _populate_lindblad_work_list!(sc.work_list, jump_hermitian, pd_el)
    end

    D = typeof(config.domain)
    C = typeof(config.construction)

    return Workspace{KrylovSpectrum, D, C, T}(
        jump_eigenbases, jump_hermitian, jumps, B_total,
        nothing,  # dll_lindblads (Thermalize path)
        G_left, G_right, G_left_adj, G_right_adj,
        channel.K0, channel.U_residual, U_coherent, nothing, Float64(delta),
        pd_transition, pd_gnf, pd_el, pd_odp, pd_nufft, pd_oft_pre_energy,
        pd_alpha, pd_bkeys, pd_bis, pd_bjs, pd_bminus, pd_bplus,
        nothing,  # coherent_unitaries
        nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing,  # trajectory fields
        nothing,  # jump_selection
        nothing,  # Id
        sc,
    )
end
