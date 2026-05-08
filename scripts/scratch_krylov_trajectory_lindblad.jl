# scripts/scratch_krylov_trajectory_lindblad.jl
#
# qf-ev5.4 prototype: Krylov-spectral-expansion trajectory predictor for
# Config{Lindbladian}. Validates the spectral-expansion approach against the
# existing ODE-based `lindblad_action_integrate` and a dense reference at n=3.
#
# The Lindbladian L is a non-normal generator on operator space. We compute:
#
#   L  R_i = lambda_i  R_i           (right eigenpairs, Krylov forward)
#   L^* L_i = conj(lambda_i)  L_i    (left eigenpairs, Krylov adjoint)
#
# biorthogonalise (L_i^* R_j = delta_ij), and reconstruct
#
#   rho(t) = rho_inf + sum_{i != 0} c_i  exp(lambda_i  t)  R_i,
#   c_i = <L_i, rho_0 - rho_inf>_HS.
#
# Trace distance to sigma_beta is then evaluated on a t-grid in O(h * d^2)
# per grid point.
#
# Run:  julia --project scripts/scratch_krylov_trajectory_lindblad.jl

using QuantumFurnace
using LinearAlgebra
using KrylovKit
using Printf

# Ensure test_helpers is on path for fixture access; fall back to inlined builder otherwise.
include(joinpath(dirname(@__DIR__), "test", "test_helpers.jl"))

# ---------------------------------------------------------------------------
# Spectral-expansion engine (prototype, scratch)
# ---------------------------------------------------------------------------

"""
    spectral_decompose_lindbladian(forward_apply!, adjoint_apply!, rho_0, dim;
                                   howmany=20, krylovdim=40, tol=1e-10)

Compute the leading `howmany` right + left eigenpairs of L matrix-free,
biorthogonalise, and project `rho_0 - rho_inf` onto the slow-mode subspace.

Returns a NamedTuple with
  lambdas::Vector{ComplexF64}            -- eigenvalues sorted by |Re| ascending
  R_modes::Vector{Matrix{ComplexF64}}    -- right eigenvectors, R_modes[1] = rho_inf
  L_modes::Vector{Matrix{ComplexF64}}    -- biorthogonalised left eigvecs
  c::Vector{ComplexF64}                  -- coefficients (c[1] should be ~0)
  rho_inf::Matrix{ComplexF64}            -- hermitised, trace-normalised steady state
  matvec_count::Int                      -- total matvecs (forward + adjoint)
  converged::Int                         -- min(forward_converged, adjoint_converged)
"""
function spectral_decompose_lindbladian(
    forward_apply!::F1,
    adjoint_apply!::F2,
    rho_0::Matrix{T},
    dim::Int;
    howmany::Int = 20,
    krylovdim::Int = 40,
    tol::Real = 1e-10,
) where {T<:Complex, F1, F2}

    # KrylovKit closures on flat vectors
    rho_buf_f = Matrix{T}(undef, dim, dim)
    out_buf_f = Matrix{T}(undef, dim, dim)
    function fwd_vec(v::AbstractVector)
        copyto!(rho_buf_f, reshape(v, dim, dim))
        forward_apply!(out_buf_f, rho_buf_f)
        return copy(vec(out_buf_f))
    end

    rho_buf_a = Matrix{T}(undef, dim, dim)
    out_buf_a = Matrix{T}(undef, dim, dim)
    function adj_vec(v::AbstractVector)
        copyto!(rho_buf_a, reshape(v, dim, dim))
        adjoint_apply!(out_buf_a, rho_buf_a)
        return copy(vec(out_buf_a))
    end

    # Initial vectors. For the FORWARD problem, x0 = vec(I/d) is fine: L(I/d) is
    # generically nonzero and Arnoldi converges normally. For the ADJOINT, however,
    # L*(I) = 0 exactly (trace preservation), so vec(I/d) is an exact eigenvector
    # and Arnoldi terminates after one matvec. Use a generic Hermitian seed for L*.
    x0_fwd = vec(Matrix{T}(I(dim) / dim))
    rng_seed = ComplexF64.(reshape(collect(range(0.1, 1.0, length=dim^2)), dim, dim))
    rng_seed .= (rng_seed .+ rng_seed') ./ 2
    rng_seed ./= norm(rng_seed)
    x0_adj = vec(rng_seed)

    # Forward eigsolve: largest real part = steady state (~0) + slow modes
    valsR, vecsR_flat, infoR = eigsolve(fwd_vec, x0_fwd, howmany, :LR,
        Arnoldi(; krylovdim=krylovdim, tol=tol, maxiter=100, verbosity=0))
    # Adjoint eigsolve: same targeting (eigenvalues are conj(lambda_i), Re unchanged)
    valsL, vecsL_flat, infoL = eigsolve(adj_vec, x0_adj, howmany, :LR,
        Arnoldi(; krylovdim=krylovdim, tol=tol, maxiter=100, verbosity=0))

    h = min(length(valsR), length(valsL))
    h >= 2 || error("need at least 2 eigenpairs (got $h)")

    # Sort forward by |Re(lambda)| ascending (steady state at index 1)
    permR = sortperm(valsR; by = v -> abs(real(v)))[1:h]
    valsR_sorted = valsR[permR]
    R_modes = [reshape(vecsR_flat[i], dim, dim) for i in permR]

    # Match adjoint eigenvalues to forward by closeness (handles near-degeneracy
    # and accidental swaps from independent Arnoldi runs).
    L_modes = Vector{Matrix{T}}(undef, h)
    used = falses(length(valsL))
    for i in 1:h
        target = conj(valsR_sorted[i])
        best_j = 0
        best_d = Inf
        for j in 1:length(valsL)
            used[j] && continue
            d = abs(valsL[j] - target)
            if d < best_d
                best_d = d
                best_j = j
            end
        end
        best_j == 0 && error("could not match adjoint eigenvalue $i (target=$target)")
        used[best_j] = true
        L_modes[i] = reshape(vecsL_flat[best_j], dim, dim)
    end

    # Steady state: hermitise R_modes[1], trace-normalise
    rho_inf = (R_modes[1] + R_modes[1]') / 2
    tr_inf = real(tr(rho_inf))
    if tr_inf == 0
        error("steady-state mode has zero trace; cannot normalise")
    end
    rho_inf ./= tr_inf
    R_modes[1] = rho_inf

    # Build h x h overlap matrix M[i, j] = tr(L_i' * R_j)
    M = Matrix{T}(undef, h, h)
    for j in 1:h, i in 1:h
        M[i, j] = dot(L_modes[i], R_modes[j])  # = tr(L_i' R_j) for matrices via column-stacking
    end

    # Solve M c = b where b[i] = <L_i, rho_0 - rho_inf>_HS = tr(L_i' (rho_0 - rho_inf))
    delta_rho = rho_0 - rho_inf
    b = Vector{T}(undef, h)
    for i in 1:h
        b[i] = dot(L_modes[i], delta_rho)
    end
    c = M \ b

    matvec_count = infoR.numops + infoL.numops
    converged = min(infoR.converged, infoL.converged)

    return (
        lambdas = Complex{Float64}.(valsR_sorted),
        R_modes = R_modes,
        L_modes = L_modes,
        c = Complex{Float64}.(c),
        rho_inf = rho_inf,
        matvec_count = matvec_count,
        converged = converged,
        M = M,
    )
end

"""
    predict_trajectory_lindblad(decomp, t_grid, sigma_beta;
                                skip_steady=true, eps_lambda=1e-10)

Reconstruct rho(t) = rho_inf + sum_i c_i  exp(lambda_i  t)  R_i and
return trace distance to `sigma_beta` on the time grid.
"""
function predict_trajectory_lindblad(
    decomp::NamedTuple,
    t_grid::AbstractVector{<:Real},
    sigma_beta::Matrix{T};
    skip_steady::Bool = true,
    eps_lambda::Real = 1e-10,
) where {T<:Complex}
    n_t = length(t_grid)
    distances = Vector{Float64}(undef, n_t)
    rho_inf = decomp.rho_inf
    rho_t = Matrix{T}(undef, size(rho_inf)...)
    for (k, t) in enumerate(t_grid)
        copyto!(rho_t, rho_inf)
        for i in eachindex(decomp.lambdas)
            # Skip the steady-state mode (its c should be ~0 for trace-pres dynamics)
            skip_steady && abs(decomp.lambdas[i]) < eps_lambda && continue
            phase = exp(decomp.lambdas[i] * t)
            rho_t .+= (decomp.c[i] * phase) .* decomp.R_modes[i]
        end
        # Defensive Hermitisation
        rho_t .= (rho_t .+ rho_t') ./ 2
        distances[k] = sum(svdvals(rho_t .- sigma_beta)) / 2
    end
    return distances
end

# ---------------------------------------------------------------------------
# Main: validation on n=3 disordered Heisenberg, beta=10
# ---------------------------------------------------------------------------

function main()
    n = 3
    beta = 10.0
    sigma = 1.0 / beta

    println("="^72)
    println("qf-ev5.4 prototype: Krylov spectral trajectory vs ODE vs dense")
    println("n=$n  beta=$beta  CKG smooth-Metropolis (a=0, s=0.25, EnergyDomain)")
    println("="^72)

    fix = make_dll_n3_system(beta)
    ham = fix.ham
    jumps = fix.jumps
    dim = size(ham.data, 1)

    cfg = Config(
        sim = Lindbladian(),
        domain = EnergyDomain(),
        construction = KMS(),
        num_qubits = n,
        with_linear_combination = true,
        beta = beta,
        sigma = sigma,
        a = 0.0, s = 0.25,
        num_energy_bits = 12,
        w0 = 0.05,
        t0 = 2π / (2^12 * 0.05),
        num_trotter_steps_per_t0 = 10,
    )

    rho_0 = Matrix{ComplexF64}(I(dim) / dim)
    sigma_beta = Matrix{ComplexF64}(ham.gibbs)

    # Spectral gap estimate via existing krylov_spectral_gap
    gap_t0 = time()
    gap_res = krylov_spectral_gap(cfg, ham, jumps; krylovdim=30, howmany=4, tol=1e-10)
    gap = gap_res.spectral_gap
    println(@sprintf("  spectral gap (Krylov): %.6e   (in %.3fs, %d matvecs)",
                     gap, time()-gap_t0, gap_res.matvec_count))

    t_max = 5.0 / gap
    t_grid = collect(range(0.0, t_max, length=41))

    # ---- Dense reference: build full L, exact e^{tL}rho_0 ----
    L_dense = construct_lindbladian(jumps, cfg, ham)
    rho_0_vec = vec(rho_0)
    sigma_beta_vec = vec(sigma_beta)
    distances_dense = Vector{Float64}(undef, length(t_grid))
    t0_d = time()
    for (k, t) in enumerate(t_grid)
        rho_t_vec = exp(t * L_dense) * rho_0_vec
        rho_t = reshape(rho_t_vec, dim, dim)
        rho_t .= (rho_t .+ rho_t') ./ 2
        distances_dense[k] = sum(svdvals(rho_t - sigma_beta)) / 2
    end
    t_dense = time() - t0_d
    println(@sprintf("  dense reference           done in %.3fs (d^2=%d)", t_dense, dim^2))

    # ---- ODE-based: lindblad_action_integrate ----
    t0_o = time()
    res_ode = integrate_to_gibbs(cfg, ham, jumps, rho_0, t_grid;
                                  mode=:L, krylovdim=30, tol=1e-10)
    t_ode = time() - t0_o
    println(@sprintf("  lindblad_action_integrate done in %.3fs (matvecs=%d)",
                     t_ode, res_ode.total_matvecs))

    # ---- Krylov spectral: prototype ----
    ws = Workspace(cfg, ham, jumps)
    fwd! = (out, x) -> begin
        apply_lindbladian!(ws, x, cfg, ham)
        copyto!(out, ws.scratch.rho_out)
    end
    adj! = (out, x) -> begin
        apply_adjoint_lindbladian!(ws, x, cfg, ham)
        copyto!(out, ws.scratch.rho_out)
    end

    t0_k = time()
    decomp = spectral_decompose_lindbladian(fwd!, adj!, rho_0, dim;
                                             howmany=20, krylovdim=40, tol=1e-10)
    distances_krylov = predict_trajectory_lindblad(decomp, t_grid, sigma_beta)
    t_krylov = time() - t0_k
    println(@sprintf("  Krylov spectral predict   done in %.3fs (matvecs=%d, conv=%d)",
                     t_krylov, decomp.matvec_count, decomp.converged))

    # ---- Comparisons ----
    println("\nDistances at sample t values:")
    println("  k     t           dense          ODE            Krylov         |Krylov-dense|/dense")
    sample_idx = [1, length(t_grid) ÷ 4, length(t_grid) ÷ 2, 3 * length(t_grid) ÷ 4, length(t_grid)]
    for k in sample_idx
        rel_diff = distances_dense[k] > 0 ?
            abs(distances_krylov[k] - distances_dense[k]) / distances_dense[k] : 0.0
        println(@sprintf("  %3d  %8.3f  %.6e  %.6e  %.6e  %.3e",
                          k, t_grid[k], distances_dense[k], res_ode.distances[k],
                          distances_krylov[k], rel_diff))
    end

    println("\n--- error norms ---")
    err_dense_ode    = maximum(abs.(res_ode.distances .- distances_dense))
    err_dense_krylov = maximum(abs.(distances_krylov .- distances_dense))
    err_ode_krylov   = maximum(abs.(distances_krylov .- res_ode.distances))
    println(@sprintf("  max |ODE - dense|:    %.3e", err_dense_ode))
    println(@sprintf("  max |Krylov - dense|: %.3e", err_dense_krylov))
    println(@sprintf("  max |Krylov - ODE|:   %.3e", err_ode_krylov))

    println("\n--- spectrum check (Krylov vs dense eigen) ---")
    eigs_dense = sort(eigvals(L_dense); by = v -> abs(real(v)))
    println("  dense eigenvalues (5 smallest by |Re|):")
    for i in 1:5
        println(@sprintf("    [%d]  %12.6e + %12.6e i", i, real(eigs_dense[i]), imag(eigs_dense[i])))
    end
    println("  Krylov eigenvalues (5 smallest by |Re|):")
    for i in 1:min(5, length(decomp.lambdas))
        println(@sprintf("    [%d]  %12.6e + %12.6e i", i, real(decomp.lambdas[i]), imag(decomp.lambdas[i])))
    end

    println("\n--- biorthogonality check (M = L'*R) ---")
    M_diag_err  = maximum(abs.(diag(decomp.M) .- 1)) # post-solve diag should be ~1 only after biorth
    println(@sprintf("  max |off-diag M|:    %.3e", maximum(abs.(decomp.M - Diagonal(diag(decomp.M))))))
    println(@sprintf("  range |diag M|:      [%.3e, %.3e]",
                     minimum(abs.(diag(decomp.M))), maximum(abs.(diag(decomp.M)))))

    println("\n--- summary ---")
    speedup_ode_krylov = res_ode.total_matvecs / decomp.matvec_count
    println(@sprintf("  ODE matvecs / Krylov matvecs = %d / %d = %.2fx",
                     res_ode.total_matvecs, decomp.matvec_count, speedup_ode_krylov))
    println(@sprintf("  PASS criterion (max |Krylov - dense| < 1e-6): %s",
                     err_dense_krylov < 1e-6 ? "PASS" : "FAIL"))
    return nothing
end

main()
