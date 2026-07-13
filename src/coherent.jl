#* COHERENT TERMS -----------------------------------------------------------------------------------------------------------
"""
    _precompute_coherent_B(
        jumps,
        hamiltonian,
        config,
        precomputed_data;
        trotter=nothing,
    ) -> Union{Nothing, Matrix{<:Complex}}

    Returns the total coherent operator B = sum_k B_k, already scaled by gamma_norm_factor.
    Returns nothing if with_coherent(config.construction) == false.
"""
function _precompute_coherent_B(
    jumps::AbstractVector{<:JumpOp},
    ham_or_trott::Union{HamHam, AbstractTrotter},
    config::Config,
    precomputed_data;
    )

    with_coherent(config.construction) || return nothing

    if config.domain isa TimeDomain
        (; b_minus, b_plus, gamma_norm_factor) = precomputed_data
        # qf-9z0.3: outer integration over `b_-(t)` runs on the `b_minus`
        # register; inner over `b_+(τ)` on the `b_plus` register.
        B = B_time(jumps, ham_or_trott, b_minus, b_plus,
            register_t0_b_minus(config), register_t0_b_plus(config),
            config.beta, config.sigma)

    elseif config.domain isa TrotterDomain
        (; b_minus, b_plus, gamma_norm_factor) = precomputed_data
        @assert ham_or_trott !== nothing
        B = B_trotter(jumps, ham_or_trott, b_minus, b_plus,
            register_t0_b_minus(config), register_t0_b_plus(config),
            config.beta, config.sigma)

    else
        # BohrDomain / EnergyDomain
        (; gamma_norm_factor) = precomputed_data
        B = B_bohr(ham_or_trott, jumps, config)
    end

    rmul!(B, gamma_norm_factor)
    return B
end

"""
    _precompute_coherent_B(jumps, ham_or_trott, config::Config{<:Any, BohrDomain, DLL}, precomputed_data)
    _precompute_coherent_B(jumps, ham_or_trott, config::Config{<:Any, TimeDomain, DLL}, precomputed_data)

DLL coherent operator `G` (Ding–Li–Lin 2024, Eqs. 3.5–3.7) computed via the
`dll_coherent_op_bohr` / `dll_coherent_op_time` helpers in `src/dll.jl`. DLL
has no `gamma_norm_factor` — the filter `f̂(ν)` already encodes the KMS
weight — so the result is returned without rescaling (cf. the CKG branch
above which applies `rmul!(B, gamma_norm_factor)`).
"""
function _precompute_coherent_B(
    jumps::AbstractVector{<:JumpOp},
    hamiltonian::HamHam,
    config::Config{<:AbstractSimulation, BohrDomain, DLL},
    precomputed_data,
    )
    (; filter) = precomputed_data
    return dll_coherent_op_bohr(jumps, hamiltonian, filter, config.beta)
end

function _precompute_coherent_B(
    jumps::AbstractVector{<:JumpOp},
    hamiltonian::HamHam,
    config::Config{<:AbstractSimulation, TimeDomain, DLL},
    precomputed_data,
    )
    (; filter, time_labels, t0) = precomputed_data
    return dll_coherent_op_time(jumps, hamiltonian, time_labels, filter, config.beta, t0)
end

"""
    _precompute_coherent_unitary(
        jumps::AbstractVector{<:JumpOp},
        hamiltonian::HamHam,
        config::Config{Thermalize},
        precomputed_data;
        trotter::Union{Nothing, AbstractTrotter}=nothing,
    ) -> Union{Nothing, Vector{Matrix{<:Complex}}}

    Precompute per-jump coherent unitaries for Kraus thermalization:
        U_k = exp(-1im * config.delta * B_k)

    Each B_k is constructed exactly as in the coherent-term definitions (domain-dependent),
    scaled by `gamma_norm_factor` (same convention as Liouvillian construction).
    Returns `nothing` if `with_coherent(config.construction) == false`.
"""
function _precompute_coherent_unitary(
    jumps::AbstractVector{<:JumpOp},
    hamiltonian::HamHam,
    config::Config{Thermalize},
    precomputed_data;
    trotter::Union{Nothing, AbstractTrotter}=nothing,
    delta_scale::Real = 1.0  # for randomized channels
    )

    with_coherent(config.construction) || return nothing

    delta = delta_scale * config.delta
    CT = Complex{eltype(hamiltonian.eigvals)}
    U_terms = Vector{Matrix{CT}}(undef, length(jumps))

    if config.domain isa TimeDomain
        (; b_minus, b_plus, gamma_norm_factor) = precomputed_data
        # Outer/inner integration registers (qf-9z0.3).
        t0_outer = register_t0_b_minus(config)
        t0_inner = register_t0_b_plus(config)
        @inbounds for (k, jump) in pairs(jumps)
            B = B_time([jump], hamiltonian, b_minus, b_plus, t0_outer, t0_inner, config.beta, config.sigma)
            rmul!(B, gamma_norm_factor)
            U_terms[k] = _coherent_unitary_step(jump, B, precomputed_data,
                t0_outer, t0_inner, delta, config.with_gqsp, config.gqsp_degree)
        end

    elseif config.domain isa TrotterDomain
        (; b_minus, b_plus, gamma_norm_factor) = precomputed_data
        @assert trotter !== nothing
        # Outer/inner integration registers (qf-9z0.3); the Trotter step
        # `trotter.t0` is independent and stays as the per-substep duration.
        t0_outer = register_t0_b_minus(config)
        t0_inner = register_t0_b_plus(config)
        @inbounds for (k, jump) in pairs(jumps)
            B = B_trotter([jump], trotter, b_minus, b_plus, t0_outer, t0_inner, config.beta, config.sigma)
            rmul!(B, gamma_norm_factor)
            U_terms[k] = _coherent_unitary_step(jump, B, precomputed_data,
                t0_outer, t0_inner, delta, config.with_gqsp, config.gqsp_degree)
        end

    else
        # BohrDomain / EnergyDomain — validate_config! prevents with_gqsp here
        (; gamma_norm_factor) = precomputed_data
        @inbounds for (k, jump) in pairs(jumps)
            B = B_bohr(hamiltonian, [jump], config)
            rmul!(B, gamma_norm_factor)
            U_terms[k] = exp(-1im * delta * Hermitian(B))
        end
    end
    return U_terms
end

# (3.1) and Proposition III.1
# Has to be on a symmetric time domain, otherwise it can't be Hermitian.
# qf-9z0.3: outer integration spacing `t0_outer` (b_-(t)) and inner spacing
# `t0_inner` (b_+(τ)) are now independent, mapping each leg to its own QPE
# register on the quantum side. Final scaling is `t0_outer * t0_inner` from
# the nested Riemann sums (thesis Eq. on line 744 of `2_methods.tex`).
function B_time(jumps::AbstractVector{<:JumpOp}, hamiltonian::HamHam,
        b_minus, b_plus, t0_outer::Real, t0_inner::Real, beta, sigma)

    d = size(hamiltonian.data, 1)
    CT = Complex{eltype(hamiltonian.eigvals)}
    eigvals = hamiltonian.eigvals

    # qf-6af.5: at construction time, the dominant cost in `_precompute_coherent_B`
    # is this nested b_+(τ) × jumps × b_-(t) loop. Parallelise the inner τ-loop
    # over (jump, τ) pairs, then the outer t-loop over `t` values, with private
    # partials reduced afterwards. Each task pins BLAS to 1 thread (Julia tasks
    # do the work; BLAS multithreading inside a per-thread mul! actively hurt
    # in measurements). Falls back to the serial inline body when nthreads()=1
    # or work below threshold (~OMEGA_THREAD_THRESHOLD).
    tau_keys = collect(keys(b_plus))
    n_jumps  = length(jumps)
    n_inner_work = length(tau_keys) * n_jumps

    if Threads.nthreads() > 1 && n_inner_work >= OMEGA_THREAD_THRESHOLD
        b_plus_summand = _b_time_inner_threaded(jumps, eigvals, b_plus, tau_keys, beta, d, CT)
    else
        b_plus_summand = zeros(CT, d, d)
        diag_u  = Vector{CT}(undef, d)
        diag_u2 = Vector{CT}(undef, d)
        tmp     = Matrix{CT}(undef, d, d)
        M       = Matrix{CT}(undef, d, d)
        for tau in tau_keys
            t_tau = tau * beta
            @. diag_u = exp(1im * eigvals * t_tau)
            @. diag_u2 = exp(-2im * eigvals * t_tau)
            diag_u_row = transpose(diag_u)
            for jump_a in jumps
                jump_eig = jump_a.in_eigenbasis
                @. tmp = diag_u2 * jump_eig
                mul!(M, jump_eig', tmp)
                b_plus_summand .+= b_plus[tau] .* diag_u .* M .* diag_u_row
            end
        end
    end

    # Outer summand b_minus — uses the outer register grid t0_outer.
    t_keys = collect(keys(b_minus))
    if Threads.nthreads() > 1 && length(t_keys) >= OMEGA_THREAD_THRESHOLD
        B = _b_time_outer_threaded(eigvals, b_plus_summand, b_minus, t_keys, sigma, d, CT)
    else
        B = zeros(CT, d, d)
        diag_u = Vector{CT}(undef, d)
        for t in t_keys
            @. diag_u = exp(1im * eigvals * (t / sigma))
            diag_u_row = transpose(diag_u)
            B .+= b_minus[t] .* conj.(diag_u) .* b_plus_summand .* diag_u_row
        end
    end

    return B .* (t0_outer * t0_inner)
end

# Threaded inner τ × jumps accumulation for `B_time`. Each task builds a
# private partial b_+(τ)-summand (d×d) using its own diag_u/diag_u2/tmp/M
# buffers. Reduction sums into a final result without holding the GIL.
function _b_time_inner_threaded(jumps, eigvals, b_plus, tau_keys, beta, d, ::Type{CT}) where {CT}
    n_inner = length(tau_keys) * length(jumps)
    nt = min(Threads.nthreads(), n_inner)
    chunks = _partition_range(1:n_inner, nt)
    n_chunks = length(chunks)
    n_jumps = length(jumps)

    partials = [zeros(CT, d, d) for _ in 1:n_chunks]

    old_blas = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    try
        @sync for (idx, chunk) in enumerate(chunks)
            Threads.@spawn _b_time_inner_chunk!(
                partials[idx], jumps, eigvals, b_plus, tau_keys,
                beta, d, n_jumps, chunk, CT)
        end
    finally
        BLAS.set_num_threads(old_blas)
    end

    summand = zeros(CT, d, d)
    @inbounds for idx in 1:n_chunks
        summand .+= partials[idx]
    end
    return summand
end

function _b_time_inner_chunk!(partial::Matrix{CT}, jumps, eigvals, b_plus,
    tau_keys, beta, d, n_jumps, chunk, ::Type{CT}) where {CT}

    diag_u  = Vector{CT}(undef, d)
    diag_u2 = Vector{CT}(undef, d)
    tmp     = Matrix{CT}(undef, d, d)
    M       = Matrix{CT}(undef, d, d)

    last_tau_idx = 0
    @inbounds for w_idx in chunk
        # Linear (tau_idx, jump_idx) decoding: outer tau, inner jump.
        tau_idx  = ((w_idx - 1) ÷ n_jumps) + 1
        jump_idx = ((w_idx - 1) % n_jumps) + 1
        if tau_idx != last_tau_idx
            tau   = tau_keys[tau_idx]
            t_tau = tau * beta
            @. diag_u  = exp(1im * eigvals * t_tau)
            @. diag_u2 = exp(-2im * eigvals * t_tau)
            last_tau_idx = tau_idx
        end
        tau    = tau_keys[tau_idx]
        b_tau  = b_plus[tau]
        jump_a = jumps[jump_idx]
        jump_eig = jump_a.in_eigenbasis
        diag_u_row = transpose(diag_u)
        @. tmp = diag_u2 * jump_eig
        mul!(M, jump_eig', tmp)
        partial .+= b_tau .* diag_u .* M .* diag_u_row
    end
    return nothing
end

# Threaded outer t-loop accumulation for `B_time`. Reads `b_plus_summand` as
# a constant (already final at this point) and accumulates B partials.
function _b_time_outer_threaded(eigvals, b_plus_summand, b_minus, t_keys,
    sigma, d, ::Type{CT}) where {CT}
    n_t = length(t_keys)
    nt  = min(Threads.nthreads(), n_t)
    chunks = _partition_range(1:n_t, nt)
    n_chunks = length(chunks)

    partials = [zeros(CT, d, d) for _ in 1:n_chunks]

    old_blas = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    try
        @sync for (idx, chunk) in enumerate(chunks)
            Threads.@spawn _b_time_outer_chunk!(
                partials[idx], eigvals, b_plus_summand, b_minus,
                t_keys, sigma, d, chunk, CT)
        end
    finally
        BLAS.set_num_threads(old_blas)
    end

    B = zeros(CT, d, d)
    @inbounds for idx in 1:n_chunks
        B .+= partials[idx]
    end
    return B
end

function _b_time_outer_chunk!(partial::Matrix{CT}, eigvals, b_plus_summand,
    b_minus, t_keys, sigma, d, chunk, ::Type{CT}) where {CT}

    diag_u = Vector{CT}(undef, d)
    @inbounds for w_idx in chunk
        t = t_keys[w_idx]
        @. diag_u = exp(1im * eigvals * (t / sigma))
        diag_u_row = transpose(diag_u)
        partial .+= b_minus[t] .* conj.(diag_u) .* b_plus_summand .* diag_u_row
    end
    return nothing
end

# Legacy 6-arg form (single t0); forwards to the explicit-outer-inner form
# with `t0_outer = t0_inner = t0`. Retained so callers that did not migrate
# still see byte-identical behaviour (B .* t0² ≡ B .* t0 .* t0).
B_time(jumps::AbstractVector{<:JumpOp}, hamiltonian::HamHam, b_minus, b_plus,
    t0::Real, beta::Real, sigma::Real) =
    B_time(jumps, hamiltonian, b_minus, b_plus, t0, t0, beta, sigma)

# qf-9z0.3: Trotter coherent gains (t0_outer, t0_inner) for the nested Riemann
# sum integration weight. This single-cache `B_trotter(::TrottTrott)` method runs
# both the inner b_+(τ) loop and the outer b_-(t) loop against the SAME Trotter
# cache (`trotter.eigvals_t0` at step `trotter.t0`). The canonical KMS coherent
# path uses `B_trotter(::TrotterTriple)` instead, which gives each leg its own
# cache and eigenbasis. Per-grid step counts:
#   - inner b_+(τ) loop: round(τ·β / t0), since τ·β advances the evolution by
#     k·t0_grid·β per τ-grid increment.
#   - outer b_-(t) loop: round(t / (σ·t0)), correct for general σ (coincides
#     with the β form when σ = 1/β).
function B_trotter(jumps::AbstractVector{<:JumpOp}, trotter::TrottTrott,
        b_minus, b_plus, t0_outer::Real, t0_inner::Real, beta, sigma)

    d = size(trotter.eigvecs, 1)
    CT = Complex{eltype(trotter.bohr_freqs)}

    # Single-cache: both legs run against the same Trotter cache, diagonal in
    # `trotter.eigvecs`.
    eigvals_outer = trotter.eigvals_t0
    eigvals_inner = trotter.eigvals_t0
    t0_step_outer = trotter.t0
    t0_step_inner = trotter.t0

    # qf-6af.5: thread the inner τ × jumps loop and the outer t-loop. See the
    # `B_time` parallelisation for rationale (Julia tasks; BLAS pinned to 1).
    tau_keys = collect(keys(b_plus))
    n_jumps  = length(jumps)
    n_inner_work = length(tau_keys) * n_jumps

    if Threads.nthreads() > 1 && n_inner_work >= OMEGA_THREAD_THRESHOLD
        b_plus_summand = _b_trotter_inner_threaded(
            jumps, eigvals_inner, b_plus, tau_keys, beta, t0_step_inner, d, CT)
    else
        b_plus_summand = zeros(CT, d, d)
        diag_u  = Vector{CT}(undef, d)
        diag_u2 = Vector{CT}(undef, d)
        tmp     = Matrix{CT}(undef, d, d)
        M       = Matrix{CT}(undef, d, d)
        for (tau, b_tau) in b_plus
            num_t0_steps = Int(round(tau * beta / t0_step_inner))
            @. diag_u  = eigvals_inner ^ num_t0_steps
            @. diag_u2 = eigvals_inner ^ (-2 * num_t0_steps)
            diag_u_row = transpose(diag_u)
            for jump_a in jumps
                jump_a_eig = jump_a.in_eigenbasis
                @. tmp = diag_u2 * jump_a_eig
                mul!(M, jump_a_eig', tmp)
                b_plus_summand .+= b_tau .* diag_u .* M .* diag_u_row
            end
        end
    end

    t_keys = collect(keys(b_minus))
    if Threads.nthreads() > 1 && length(t_keys) >= OMEGA_THREAD_THRESHOLD
        B = _b_trotter_outer_threaded(
            eigvals_outer, b_plus_summand, b_minus, t_keys, sigma, t0_step_outer, d, CT)
    else
        B = zeros(CT, d, d)
        diag_u = Vector{CT}(undef, d)
        for (t, b_t) in b_minus
            num_t0_steps = Int(round(t / (sigma * t0_step_outer)))
            @. diag_u = eigvals_outer ^ num_t0_steps
            diag_u_row = transpose(diag_u)
            B .+= b_t .* conj.(diag_u) .* b_plus_summand .* diag_u_row
        end
    end

    return B .* (t0_outer * t0_inner)  # B in Trotter basis
end

function _b_trotter_inner_threaded(jumps, eigvals_inner, b_plus, tau_keys,
    beta, t0_step_inner, d, ::Type{CT}) where {CT}
    n_inner = length(tau_keys) * length(jumps)
    nt = min(Threads.nthreads(), n_inner)
    chunks = _partition_range(1:n_inner, nt)
    n_chunks = length(chunks)
    n_jumps  = length(jumps)

    partials = [zeros(CT, d, d) for _ in 1:n_chunks]

    old_blas = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    try
        @sync for (idx, chunk) in enumerate(chunks)
            Threads.@spawn _b_trotter_inner_chunk!(
                partials[idx], jumps, eigvals_inner, b_plus, tau_keys,
                beta, t0_step_inner, d, n_jumps, chunk, CT)
        end
    finally
        BLAS.set_num_threads(old_blas)
    end

    summand = zeros(CT, d, d)
    @inbounds for idx in 1:n_chunks
        summand .+= partials[idx]
    end
    return summand
end

function _b_trotter_inner_chunk!(partial::Matrix{CT}, jumps, eigvals_inner,
    b_plus, tau_keys, beta, t0_step_inner, d, n_jumps, chunk, ::Type{CT}) where {CT}

    diag_u  = Vector{CT}(undef, d)
    diag_u2 = Vector{CT}(undef, d)
    tmp     = Matrix{CT}(undef, d, d)
    M       = Matrix{CT}(undef, d, d)

    last_tau_idx = 0
    @inbounds for w_idx in chunk
        tau_idx  = ((w_idx - 1) ÷ n_jumps) + 1
        jump_idx = ((w_idx - 1) % n_jumps) + 1
        if tau_idx != last_tau_idx
            tau   = tau_keys[tau_idx]
            num_t0_steps = Int(round(tau * beta / t0_step_inner))
            @. diag_u  = eigvals_inner ^ num_t0_steps
            @. diag_u2 = eigvals_inner ^ (-2 * num_t0_steps)
            last_tau_idx = tau_idx
        end
        tau   = tau_keys[tau_idx]
        b_tau = b_plus[tau]
        jump_a = jumps[jump_idx]
        jump_a_eig = jump_a.in_eigenbasis
        diag_u_row = transpose(diag_u)
        @. tmp = diag_u2 * jump_a_eig
        mul!(M, jump_a_eig', tmp)
        partial .+= b_tau .* diag_u .* M .* diag_u_row
    end
    return nothing
end

function _b_trotter_outer_threaded(eigvals_outer, b_plus_summand, b_minus, t_keys,
    sigma, t0_step_outer, d, ::Type{CT}) where {CT}
    n_t = length(t_keys)
    nt  = min(Threads.nthreads(), n_t)
    chunks = _partition_range(1:n_t, nt)
    n_chunks = length(chunks)

    partials = [zeros(CT, d, d) for _ in 1:n_chunks]

    old_blas = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    try
        @sync for (idx, chunk) in enumerate(chunks)
            Threads.@spawn _b_trotter_outer_chunk!(
                partials[idx], eigvals_outer, b_plus_summand, b_minus,
                t_keys, sigma, t0_step_outer, d, chunk, CT)
        end
    finally
        BLAS.set_num_threads(old_blas)
    end

    B = zeros(CT, d, d)
    @inbounds for idx in 1:n_chunks
        B .+= partials[idx]
    end
    return B
end

function _b_trotter_outer_chunk!(partial::Matrix{CT}, eigvals_outer,
    b_plus_summand, b_minus, t_keys, sigma, t0_step_outer, d, chunk,
    ::Type{CT}) where {CT}

    diag_u = Vector{CT}(undef, d)
    @inbounds for w_idx in chunk
        t = t_keys[w_idx]
        b_t = b_minus[t]
        num_t0_steps = Int(round(t / (sigma * t0_step_outer)))
        @. diag_u = eigvals_outer ^ num_t0_steps
        diag_u_row = transpose(diag_u)
        partial .+= b_t .* conj.(diag_u) .* b_plus_summand .* diag_u_row
    end
    return nothing
end

# Legacy form: forwards to the explicit-outer-inner form using `trotter.t0`
# for both legs (matches the pre-qf-9z0 behaviour `B .* trotter.t0²`).
B_trotter(jumps::AbstractVector{<:JumpOp}, trotter::TrottTrott, b_minus, b_plus,
    beta::Real, sigma::Real) =
    B_trotter(jumps, trotter, b_minus, b_plus, trotter.t0, trotter.t0, beta, sigma)

# ---------------------------------------------------------------------------
# qf-e4z.20.3 — B_trotter on TrotterTriple (independent per-leg caches).
#
# Pipeline (jumps arrive in V_D basis; final B is returned in V_D):
#
#   1. Rotate jumps V_D → V_bp:   J_bp = R_bp_in_D · J_D · R_bp_in_D'.
#   2. Inner τ-loop in V_bp:      runs against `triple.b_plus.eigvals_t0`,
#                                  `triple.b_plus.t0`. Output: b_+_summand in V_bp.
#   3. Rotate summand V_bp → V_bm: summand_bm = R_bm_in_bp · summand_bp · R_bm_in_bp'.
#   4. Outer t-loop in V_bm:      runs against `triple.b_minus.eigvals_t0`,
#                                  `triple.b_minus.t0`. Output: B in V_bm.
#   5. Rotate B V_bm → V_D:       B_D = R_bm_in_D' · B_bm · R_bm_in_D.
#
# Convention: R_{Y←X} := V_Y' · V_X, so M_Y = R_{Y←X} · M_X · R_{Y←X}'.
# The rotations are O(d^3) per call (6 GEMMs) — sub-dominant to the
# threaded τ × jumps inner loop at d ≤ 64 and competitive at d ≤ 512.
# ---------------------------------------------------------------------------
function B_trotter(
    jumps::AbstractVector{<:JumpOp},
    triple::TrotterTriple,
    b_minus, b_plus,
    t0_outer::Real, t0_inner::Real,
    beta, sigma,
)
    d  = size(triple.D.eigvecs, 1)
    CT = Complex{eltype(triple.D.bohr_freqs)}

    # Step 1: rotate jumps V_D → V_bp.  J_bp = R_bp_in_D · J_D · R_bp_in_D'.
    R_bp_in_D = triple.R_bp_in_D
    jumps_bp  = Vector{JumpOp}(undef, length(jumps))
    @inbounds for (k, j) in pairs(jumps)
        j_bp = Matrix{CT}(R_bp_in_D * j.in_eigenbasis * R_bp_in_D')
        jumps_bp[k] = JumpOp(j.data, j_bp, j.orthogonal, j.hermitian)
    end

    # Step 2: inner τ-loop in V_bp. Reuses the existing threaded helper —
    # `eigvals_inner = triple.b_plus.eigvals_t0` is diagonal in V_bp.
    eigvals_inner = triple.b_plus.eigvals_t0
    t0_step_inner = triple.b_plus.t0
    tau_keys = collect(keys(b_plus))
    n_jumps  = length(jumps_bp)
    n_inner_work = length(tau_keys) * n_jumps

    if Threads.nthreads() > 1 && n_inner_work >= OMEGA_THREAD_THRESHOLD
        b_plus_summand_bp = _b_trotter_inner_threaded(
            jumps_bp, eigvals_inner, b_plus, tau_keys, beta, t0_step_inner, d, CT)
    else
        b_plus_summand_bp = zeros(CT, d, d)
        diag_u  = Vector{CT}(undef, d)
        diag_u2 = Vector{CT}(undef, d)
        tmp     = Matrix{CT}(undef, d, d)
        M       = Matrix{CT}(undef, d, d)
        for (tau, b_tau) in b_plus
            num_t0_steps = Int(round(tau * beta / t0_step_inner))
            @. diag_u  = eigvals_inner ^ num_t0_steps
            @. diag_u2 = eigvals_inner ^ (-2 * num_t0_steps)
            diag_u_row = transpose(diag_u)
            for jump_a in jumps_bp
                jump_a_eig = jump_a.in_eigenbasis
                @. tmp = diag_u2 * jump_a_eig
                mul!(M, jump_a_eig', tmp)
                b_plus_summand_bp .+= b_tau .* diag_u .* M .* diag_u_row
            end
        end
    end

    # Step 3: rotate summand V_bp → V_bm.  summand_bm = R_bm_in_bp · summand_bp · R_bm_in_bp'.
    R_bm_in_bp        = triple.R_bm_in_bp
    b_plus_summand_bm = Matrix{CT}(R_bm_in_bp * b_plus_summand_bp * R_bm_in_bp')

    # Step 4: outer t-loop in V_bm. `eigvals_outer = triple.b_minus.eigvals_t0`
    # is diagonal in V_bm.
    eigvals_outer = triple.b_minus.eigvals_t0
    t0_step_outer = triple.b_minus.t0
    t_keys = collect(keys(b_minus))
    if Threads.nthreads() > 1 && length(t_keys) >= OMEGA_THREAD_THRESHOLD
        B_bm = _b_trotter_outer_threaded(
            eigvals_outer, b_plus_summand_bm, b_minus, t_keys, sigma, t0_step_outer, d, CT)
    else
        B_bm = zeros(CT, d, d)
        diag_u = Vector{CT}(undef, d)
        for (t, b_t) in b_minus
            num_t0_steps = Int(round(t / (sigma * t0_step_outer)))
            @. diag_u = eigvals_outer ^ num_t0_steps
            diag_u_row = transpose(diag_u)
            B_bm .+= b_t .* conj.(diag_u) .* b_plus_summand_bm .* diag_u_row
        end
    end

    # Step 5: rotate V_bm → V_D.  B_D = R_bm_in_D' · B_bm · R_bm_in_D.
    R_bm_in_D = triple.R_bm_in_D
    B_D = Matrix{CT}(R_bm_in_D' * B_bm * R_bm_in_D)

    return B_D .* (t0_outer * t0_inner)  # B in V_D
end

# Legacy form for TrotterTriple: forwards using `triple.D.t0` for both grid
# spacings (matches the pre-qf-9z0 single-cache convention).
B_trotter(jumps::AbstractVector{<:JumpOp}, triple::TrotterTriple, b_minus, b_plus,
    beta::Real, sigma::Real) =
    B_trotter(jumps, triple, b_minus, b_plus, triple.D.t0, triple.D.t0, beta, sigma)

#* GQSP POLYNOMIAL APPROXIMATION ------------------------------------------------------------------------------------------
"""
    _coherent_unitary_step(jump, B, precomputed_data, t0_outer, t0_inner, delta_eff,
                           with_gqsp, gqsp_degree) -> Matrix{<:Complex}

Single-jump coherent step `U_a ≈ exp(-i·delta_eff·B_a)`. `B` is assumed already scaled by
`gamma_norm_factor` (same convention as the dissipator construction). It is hermitised
in-place to absorb numerical noise; both branches then operate on the same Hermitian
matrix:
- `with_gqsp = false`: returns `exp(-i·delta_eff·Hermitian(B))` (exact matrix exp).
- `with_gqsp = true`: returns the post-selected GQSP polynomial `f_{gqsp_degree}(B/α)`
  approximating `exp(-i·delta_eff·B_a)` with truncation error
  `O((delta_eff·α)^{gqsp_degree+1})` (Bessel-tail bound; MW 2024 Eq. 62–63).
  `α = _gqsp_block_encoding_alpha(...)` is built from the **outer** and **inner**
  Riemann-sum spacings independently (qf-9z0.3). Reads `b_minus`, `b_plus`,
  `gamma_norm_factor` from `precomputed_data` only in this branch — `validate_config!`
  guarantees those fields exist when `with_gqsp = true` (Time/TrotterDomain only).
"""
function _coherent_unitary_step(
    jump::JumpOp,
    B::AbstractMatrix{<:Complex},
    precomputed_data,
    t0_outer::Real,
    t0_inner::Real,
    delta_eff::Real,
    with_gqsp::Bool,
    gqsp_degree::Int,
)
    hermitianize!(B)
    if with_gqsp
        α_be = _gqsp_block_encoding_alpha(jump,
            precomputed_data.b_minus, precomputed_data.b_plus,
            t0_outer, t0_inner, precomputed_data.gamma_norm_factor)
        return _gqsp_apply_polynomial(B, α_be, delta_eff, gqsp_degree)
    else
        return exp(-1im * delta_eff * Hermitian(B))
    end
end

# Legacy 7-arg form: forwards to the explicit-outer-inner form with
# `t0_outer = t0_inner = t0_sim`.
_coherent_unitary_step(jump, B, precomputed_data, t0_sim::Real, delta_eff::Real,
    with_gqsp::Bool, gqsp_degree::Int) =
    _coherent_unitary_step(jump, B, precomputed_data, t0_sim, t0_sim, delta_eff,
        with_gqsp, gqsp_degree)

"""
    _gqsp_block_encoding_alpha(jump, b_minus, b_plus, t0_outer, t0_inner,
                               gamma_norm_factor) -> Real

Block-encoding norm `α_a` of the simulator-side `B_a` operator (Time/TrotterDomain),
with **independent outer/inner** Riemann-sum spacings (qf-9z0.3).

Following Alg. `alg:coh` in the thesis, with `B_a` already scaled by `gamma_norm_factor`:

    α_a = γ_nf · t0_outer · t0_inner · ‖b_-‖_{ℓ¹} · ‖b_+‖_{ℓ¹} · ‖A_a‖²_op

where the ℓ¹ norms are the truncated-grid sums of `|b_minus[t]|` and `|b_plus[τ]|`.
The `‖A_a‖²` factor (not `‖A_a‖`) reflects the two `A_a` factors in the inner
integrand `A_a^† · e^{-2iHβτ} · A_a` (the inner unitary has operator norm 1).
This guarantees `‖B_a / α_a‖_op ≤ 1` so the GQSP polynomial of `B_a / α_a` is faithful.
"""
function _gqsp_block_encoding_alpha(
    jump::JumpOp,
    b_minus,
    b_plus,
    t0_outer::Real,
    t0_inner::Real,
    gamma_norm_factor::Real,
)
    l1_minus = sum(abs, values(b_minus))
    l1_plus  = sum(abs, values(b_plus))
    A_norm_sq = opnorm(jump.data)^2
    return gamma_norm_factor * t0_outer * t0_inner * l1_minus * l1_plus * A_norm_sq
end

# Legacy 5-arg form: forwards to the explicit-outer-inner form.
_gqsp_block_encoding_alpha(jump, b_minus, b_plus, t0_sim::Real, gamma_norm_factor::Real) =
    _gqsp_block_encoding_alpha(jump, b_minus, b_plus, t0_sim, t0_sim, gamma_norm_factor)

"""
    _gqsp_apply_polynomial(B::AbstractMatrix{<:Complex}, alpha::Real, delta::Real, d::Int)

Compute the post-selected anc=|0⟩ block of the GQSP circuit at degree `d`, i.e. the
Chebyshev expansion produced by qubitization + Jacobi-Anger truncation:

    f_d(B/α) = J_0(δα) I + Σ_{k=1}^{d} 2 (-i)^k J_k(δα) T_k(B/α)

This is the post-selected anc=|0⟩ block of the GQSP Laurent polynomial
`L_d(W) = Σ_{n=-d}^{d} (-i)^n J_n(δα) W^n` evaluated on the qubitization walk `W`
of `B/α` (Motlagh & Wiebe 2024, Theorem 7 + Cor. 8, Eq. 62–66). The `(-i)^n`
phase comes from MW Eq. 62 `e^{it cos θ} = Σ_n i^n J_n(t) e^{inθ}` with `t = -δα`
combined with `J_n(-x) = (-1)^n J_n(x)`; the factor `2` on `k ≥ 1` collapses the
two-sided sum onto Chebyshev `T_k(cos θ) = cos(kθ)`.

To `O((δα)^{d+1})` this approximates `exp(-iδ B)` (Bessel-tail bound,
MW Eq. 63). The d=1 special case `f_1(B/α) = J_0(δα) I − 2i J_1(δα) (B/α)` requires
no matmul. For `d ≥ 2`, the function uses the Clenshaw recurrence on `T_k(B/α)`
with three n×n scratch buffers reused across iterations (allocation is `O(1)` in d).

Returns a `Matrix{eltype(B)}`. `f_d` is unitary up to `O((δα)^{d+1})` but in general not
Hermitian; do not Hermitianize the output.

# Circuit form (cost-model only — has no effect on this evaluator)

This routine is a circuit-form-agnostic Clenshaw evaluator of `f_d(B/α)`; it
returns the post-selected QSP=|0⟩ block of the GQSP unitary, which is the
*same* operator for both Form B and Form C (algebraic identity from MW2024
Eqs. 49→53). The implementation target on hardware is **Form B (MW2024 Eq. 46)**:
`d` open-controlled-`W` slots + `d` closed-controlled-`W†` slots (the `A'`
of MW Eq. 45) = `2d` block-encoding queries — 1.5× cheaper than Form C
(Eq. 52: all controlled-`W` + uncontrolled `W^{-d}` tail = `3d` queries).
The Hamiltonian-simulation cost model in `src/simulation_time.jl` charges
the Form-B count. See `src/python/tests/test_gqsp.py::test_form_b_equivalent_to_form_c`
for the numerical confirmation that Form B ≡ Form C on the post-selected
block (≤ 1e-12, qf-e4z.19).

Reference: scripts/scratch_gqsp_B_n3.jl, scripts/scratch_gqsp_random_h.jl (qf-0x6 POC).
"""
function _gqsp_apply_polynomial(
    B::AbstractMatrix{<:Complex},
    alpha::Real,
    delta::Real,
    d::Int,
)
    @assert d ≥ 1 "gqsp_degree must be ≥ 1 (got $d)"
    n = LinearAlgebra.checksquare(B)
    CT = eltype(B)

    delta_alpha = delta * alpha
    a0 = CT(besselj(0, delta_alpha))

    # d = 1 fast path: f_1 = J_0 I − 2i J_1 (B/α)  -- scalar-axpy, no matmul
    if d == 1
        a1_over_alpha = CT(-2im * besselj(1, delta_alpha) / alpha)
        result = Matrix{CT}(undef, n, n)
        @inbounds for j in 1:n, i in 1:n
            result[i, j] = a1_over_alpha * B[i, j]
        end
        @inbounds for i in 1:n
            result[i, i] += a0
        end
        return result
    end

    # d ≥ 2: Clenshaw recurrence on T_k(x), x = B/α
    # f(x) = a_0 + Σ_{k=1}^d a_k T_k(x), a_k = 2 (-i)^k J_k(δα)
    # b_{d+1} = b_{d+2} = 0;  for k = d, d-1, …, 1: b_k = a_k I + 2 x b_{k+1} - b_{k+2}
    # f(x) = a_0 I + x b_1 - b_2
    inv_alpha = inv(alpha)
    x = Matrix{CT}(undef, n, n)
    @inbounds for j in 1:n, i in 1:n
        x[i, j] = B[i, j] * inv_alpha
    end

    b_kp2 = zeros(CT, n, n)
    b_kp1 = zeros(CT, n, n)
    b_k   = Matrix{CT}(undef, n, n)
    tmp   = Matrix{CT}(undef, n, n)

    @inbounds for k in d:-1:1
        ak = CT(2 * cis(-π/2 * k) * besselj(k, delta_alpha))
        mul!(tmp, x, b_kp1)            # tmp = x · b_{k+1}
        @. b_k = 2 * tmp - b_kp2       # b_k = 2 x b_{k+1} − b_{k+2}
        for i in 1:n
            b_k[i, i] += ak            # add a_k I on diagonal
        end
        # Rotate buffers: next iter's (b_{k+2}, b_{k+1}) ← (b_{k+1}, b_k); old b_{k+2} freed
        b_kp2, b_kp1, b_k = b_kp1, b_k, b_kp2
    end

    # b_kp1 holds b_1, b_kp2 holds b_2
    result = Matrix{CT}(undef, n, n)
    mul!(result, x, b_kp1)
    @. result = result - b_kp2
    @inbounds for i in 1:n
        result[i, i] += a0
    end
    return result
end

#* B1 AND B2 ---------------------------------------------------------------------------------------------------------------
function _compute_b_minus(t::Real, beta::Real, sigma::Real)  # 2pi sqrt(pi) * f_minus(t / sigma_E)
    f1(t) = 1 / cosh(2 * pi * t / (beta * sigma))
    f2(t) = sin(-t * beta * sigma) * exp(-2 * t^2)
    return 2 * sqrt(pi) * exp(beta^2 * sigma^2 / 8) * _convolute(f1, f2, t) / (beta * sigma)
end

function _compute_b_plus(t::Real, beta::Real, w_gamma::Real, sigma_gamma::Real)  # f_plus(t * beta) / (2pi sqrt(pi))
    return beta * sigma_gamma * exp(- 2 * beta * w_gamma * (2 * t^2 + im * t)) / sqrt(pi^3)
end

function _compute_b_plus_metro(t::Real, beta::Real, sigma::Real, eta::Real, s::Real=0.0)
    if abs(t) < 1e-12  # Handle t=0 (L'Hopital limit; reduces to 1/(2√2 π²) at σβ=1, s=0)
        return complex((2 - sigma^2 * beta^2 * (1 + s)) / (2 * sqrt(2) * pi^2))
    elseif abs(t) <= eta
        numerator = exp(- sigma^2 * beta^2 * (2 * t^2 + 1im * t) * (1 + s)) + 1im * (2 * t + 1im)
    else
        numerator = exp(- sigma^2 * beta^2 * (2 * t^2 + 1im * t) * (1 + s))
    end
    denominator = t * (2 * t + 1im)
    return (1 / (2 * sqrt(2) * pi^2)) * numerator / denominator
end

function _compute_b_plus_smooth(t::Real, beta::Real, sigma::Real, a::Real, s::Real)
    b_vals = exp(- a * s / 2) * exp(- sigma^2 * beta^2 * t * (2 * t + 1im) * (1 + s)) / (4 * t^2 + a + 2im * t)
    return sqrt(4 * a + 1) * b_vals / (sqrt(2) * pi^2)
end

function _compute_truncated_func(target_func::Function, time_labels::AbstractVector{<:Real}, fixed_args...; atol::Real = 1e-12)
    f_vals = Vector{ComplexF64}(target_func.(time_labels, fixed_args...))
    indices_to_keep = _get_truncated_indices(f_vals; atol=atol)
    return Dict(zip(time_labels[indices_to_keep], f_vals[indices_to_keep]))
end

#* TOOLS --------------------------------------------------------------------------------------------------------------------
function _get_truncated_indices(fvals::AbstractVector{<:Number}; atol::Real = 1e-12)
    return findall(abs.(fvals) .>= atol)
end

function _convolute(f::Function, g::Function, t::Real; atol=1e-12, rtol=1e-12)
    integrand(s) = f(s) * g(t - s)
    result, _ = quadgk(integrand, -Inf, Inf; atol=atol, rtol=rtol)
    return result
end
