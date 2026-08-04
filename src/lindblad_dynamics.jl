# Matrix-free Krylov integration of Lindbladian and discriminant flows.
# Math: $dot(rho) = L(rho)$ and $psi = sigma^(-1/4) rho sigma^(-1/4)$.

"""
    lindblad_action_integrate(L_apply!, rho_0, sigma_beta, t_grid;
                              krylovdim=30, tol=1e-10, save_states=false)

Integrate a matrix-free Lindblad equation over an ordered time grid.

# Arguments
- `L_apply!`: In-place closure writing `L(input)` into its first argument.
- `rho_0`: Initial density matrix.
- `sigma_beta`: Reference state for trace distance.
- `t_grid`: Ordered sample times.

# Keywords
- `krylovdim`: Arnoldi subspace size.
- `tol`: Per-step Krylov tolerance.
- `save_states`: Retain every propagated state.

# Returns
A named tuple with times, trace distances, final state, matvec count,
convergence status, and optional states.
"""
function lindblad_action_integrate(
    L_apply!::F,
    rho_0::Matrix{T},
    sigma_beta::Matrix{T},
    t_grid::AbstractVector{<:Real};
    krylovdim::Int = 30,
    tol::Real = 1e-10,
    save_states::Bool = false,
)::NamedTuple where {T<:Complex, F}
    d = size(rho_0, 1)
    @assert size(rho_0, 2) == d  "rho_0 must be square"

    # Persistent buffers captured by the closure (no allocations in the hot loop).
    rho_buf = Matrix{T}(undef, d, d)
    out_buf = Matrix{T}(undef, d, d)

    # KrylovKit-side closure on flat vectors. KrylovKit may overwrite the
    # buffer `out_buf` on the next call, so we MUST `copy(vec(out_buf))`
    # before returning; we also `copyto!` the input into our private buffer
    # in case `v` aliases internal state.
    function L_vec_apply(v::AbstractVector)
        copyto!(rho_buf, reshape(v, d, d))
        L_apply!(out_buf, rho_buf)
        return copy(vec(out_buf))
    end

    n_steps   = length(t_grid)
    distances = Vector{Float64}(undef, n_steps)
    states    = save_states ? Vector{Matrix{T}}(undef, n_steps) : Matrix{T}[]

    # Initial state: defensive copy + trace-distance.
    rho     = copy(rho_0)
    distances[1] = sum(svdvals(rho - sigma_beta)) / 2
    save_states && (states[1] = copy(rho))

    # Working flat vector for Krylov.
    v_rho   = copy(vec(rho))

    total_matvecs = 0
    all_converged = true

    @inbounds for i in 1:(n_steps - 1)
        dt = float(t_grid[i + 1] - t_grid[i])

        v_next, info = exponentiate(L_vec_apply, dt, v_rho;
                                    krylovdim = krylovdim,
                                    tol = tol,
                                    ishermitian = false)
        total_matvecs += info.numops
        if info.converged == 0
            all_converged = false
            @warn "L-mode exponentiate did not converge at step" i numops=info.numops
        end
        copyto!(v_rho, v_next)

        # Reshape into rho, then defensively re-Hermitise + re-trace-normalise:
        # Davies preserves Hermiticity and trace exactly in the continuum, but
        # Krylov truncation introduces O(tol)-level violations.
        copyto!(rho, reshape(v_rho, d, d))
        @inbounds for j in 1:d, k in 1:d
            rho[k, j] = (rho[k, j] + conj(rho[j, k])) / 2
        end
        # rho is now Hermitian; re-vec and renormalise trace (only the real diagonal).
        tr_now = real(tr(rho))
        if tr_now != 0
            rho ./= tr_now
        end
        copyto!(v_rho, vec(rho))

        distances[i + 1] = sum(svdvals(rho - sigma_beta)) / 2
        save_states && (states[i + 1] = copy(rho))
    end

    return (
        t              = collect(t_grid),
        distances      = distances,
        rho_final      = copy(rho),
        total_matvecs  = total_matvecs,
        all_converged  = all_converged,
        states         = states,
    )
end


"""
    discriminant_action_integrate(K_apply!, psi_0, psi_eq, t_grid;
                                  krylovdim=30, tol=1e-10,
                                  is_hermitian=true, save_states=false)

Integrate a matrix-free discriminant equation over an ordered time grid.

Set `is_hermitian=true` only when KMS detailed balance makes the discriminant
Hilbert--Schmidt self-adjoint; otherwise Arnoldi is required.

# Returns
A named tuple with times, Frobenius distances, final state, matvec count,
convergence status, and optional states.
"""
function discriminant_action_integrate(
    K_apply!::F,
    psi_0::Matrix{T},
    psi_eq::Matrix{T},
    t_grid::AbstractVector{<:Real};
    krylovdim::Int = 30,
    tol::Real = 1e-10,
    is_hermitian::Bool = true,
    save_states::Bool = false,
)::NamedTuple where {T<:Complex, F}
    d = size(psi_0, 1)
    @assert size(psi_0, 2) == d  "psi_0 must be square"

    in_buf  = Matrix{T}(undef, d, d)
    out_buf = Matrix{T}(undef, d, d)

    # KrylovKit may overwrite `out_buf` on the next call, so we
    # `copy(vec(out_buf))` before returning; we also `copyto!` the input
    # into our private buffer in case `v` aliases internal state.
    function K_vec_apply(v::AbstractVector)
        copyto!(in_buf, reshape(v, d, d))
        K_apply!(out_buf, in_buf)
        return copy(vec(out_buf))
    end

    n_steps   = length(t_grid)
    distances = Vector{Float64}(undef, n_steps)
    states    = save_states ? Vector{Matrix{T}}(undef, n_steps) : Matrix{T}[]

    psi = copy(psi_0)
    distances[1] = norm(psi - psi_eq)             # Frobenius distance (chi metric)
    save_states && (states[1] = copy(psi))

    v_psi = copy(vec(psi))
    total_matvecs = 0
    all_converged = true

    @inbounds for i in 1:(n_steps - 1)
        dt = float(t_grid[i + 1] - t_grid[i])

        v_next, info = exponentiate(K_vec_apply, dt, v_psi;
                                    krylovdim = krylovdim,
                                    tol = tol,
                                    ishermitian = is_hermitian)
        total_matvecs += info.numops
        if info.converged == 0
            all_converged = false
            @warn "K-mode exponentiate did not converge at step" i numops=info.numops
        end
        copyto!(v_psi, v_next)
        copyto!(psi, reshape(v_psi, d, d))

        # KMS-DB preserves Hermiticity and $<psi_eq,psi>_F = tr(rho) = 1$;
        # correct accumulated Krylov round-off in both invariants.
        @inbounds for j in 1:d, k in 1:d
            psi[k, j] = (psi[k, j] + conj(psi[j, k])) / 2
        end
        c_now = real(dot(psi_eq, psi))
        if c_now != 0
            psi ./= c_now
        end
        copyto!(v_psi, vec(psi))

        distances[i + 1] = norm(psi - psi_eq)
        save_states && (states[i + 1] = copy(psi))
    end

    return (
        t              = collect(t_grid),
        distances      = distances,
        psi_final      = copy(psi),
        total_matvecs  = total_matvecs,
        all_converged  = all_converged,
        states         = states,
    )
end


"""
    integrate_to_gibbs(config, hamiltonian, jumps, rho_0, t_grid;
                       mode=:L, krylovdim=30, tol=1e-10, save_states=false)

Integrate a configured Lindbladian or its KMS discriminant toward Gibbs.

# Arguments
- `config`: Bohr- or energy-domain Lindbladian configuration.
- `hamiltonian`: Hamiltonian with cached Gibbs state.
- `jumps`: Coupling operators.
- `rho_0`: Initial density matrix.
- `t_grid`: Ordered sample times.

# Keywords
- `mode`: `:L` for density matrices or `:K` for the discriminant representation.
- `krylovdim`, `tol`, `save_states`: Forwarded to the selected integrator.
- `allow_unpaired_nonhermitian`: Opt out of adjoint-pair validation.

# Returns
The selected integrator's named tuple. Configuration and DLL filter constraints
are validated before constructing the workspace.
"""
function integrate_to_gibbs(
    config::Config{Lindbladian, <:Union{BohrDomain, EnergyDomain}},
    hamiltonian::HamHam,
    jumps::Vector{JumpOp},
    rho_0::Matrix{T},
    t_grid::AbstractVector{<:Real};
    mode::Symbol = :L,
    krylovdim::Int = 30,
    tol::Real = 1e-10,
    save_states::Bool = false,
    allow_unpaired_nonhermitian::Bool = false,
)::NamedTuple where {T<:Complex}
    mode in (:L, :K) || throw(ArgumentError("mode must be :L or :K (got :$mode)"))
    d = size(rho_0, 1)
    @assert size(rho_0, 2) == d  "rho_0 must be square"

    # validate_config! is invoked by run_lindblad/run_thermalize; we call it
    # explicitly here since this entry point bypasses those.
    validate_config!(config)
    validate_jump_pairing(jumps; allow_unpaired_nonhermitian=allow_unpaired_nonhermitian)

    # Build the L_apply!(out, in) closure: matrix-free for both KMS and DLL.
    # `let` scope binds the captured state directly to dodge Box wrapping under
    # Julia 1.11+ closure capture rules.
    ws = Workspace(config, hamiltonian, jumps)
    L_apply! = let ws = ws, config = config, ham = hamiltonian
        (out::AbstractMatrix, x::AbstractMatrix) -> begin
            apply_lindbladian!(ws, x, config, ham)
            copyto!(out, ws.scratch.rho_out)
            return out
        end
    end

    if mode == :L
        sigma_beta = Matrix{T}(hamiltonian.gibbs)
        return lindblad_action_integrate(
            L_apply!, rho_0, sigma_beta, t_grid;
            krylovdim = krylovdim, tol = tol, save_states = save_states,
        )
    else  # mode == :K
        powers = gibbs_fractional_powers(hamiltonian.gibbs)
        sq, sq_inv, sh = powers.sigma_quarter, powers.sigma_inv_quarter, powers.sigma_half

        # psi_0 = sigma^{-1/4} rho_0 sigma^{-1/4} (diagonal multiply, BohrDomain).
        psi_0 = Matrix{T}(undef, d, d)
        @inbounds for j in 1:d, i in 1:d
            psi_0[i, j] = sq_inv[i] * rho_0[i, j] * sq_inv[j]
        end
        # psi_eq = sigma^{1/2} as a full Matrix (the integrator wants a Matrix, not Diagonal).
        psi_eq = Matrix{T}(Diagonal(complex.(sh)))

        bufs = DiscriminantBuffers{T}(d)
        K_apply! = let L = L_apply!, sq = sq, sq_inv = sq_inv, bufs = bufs
            (out::AbstractMatrix, x::AbstractMatrix) -> begin
                apply_discriminant!(out, x, L, sq, sq_inv, bufs)
                return out
            end
        end

        return discriminant_action_integrate(
            K_apply!, psi_0, psi_eq, t_grid;
            krylovdim = krylovdim, tol = tol,
            is_hermitian = true,  # KMS-DB ⇒ K is HS-self-adjoint (Lanczos OK)
            save_states = save_states,
        )
    end
end
