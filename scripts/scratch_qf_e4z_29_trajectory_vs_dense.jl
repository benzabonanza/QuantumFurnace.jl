"""
qf-e4z.29 — Direct measurement of `predict_lindbladian_trajectory`'s
trajectory reconstruction accuracy vs a DENSE `exp(L·t)` reference.

The qf-e4z.27 design assumes (by symmetry argument) that single-pass
Arnoldi from `vec(rho_0=I/d)` reconstructs the true trajectory of L
acting on I/d, even on parity-symmetric L where the captured spectrum
is restricted to the parity-EVEN sector and the reported "spectral_gap"
is wrong by 17–24%. Nobody has measured this directly — the qf-e4z.27
and qf-e4z.28 verification scripts only compared spectral-gap values.

This script does the trajectory check at n=5 only (where the dense
1024×1024 Lindbladian fits in memory and `exp(L·t)` is tractable).
Compares the predictor's reconstructed states[k] against the dense
reference at every t-grid point.

Acceptance:
  • Z+ZZ (exact parity): error ≤ 1e-7 at every t (symmetry argument
    predicts byte-equivalence up to Krylov truncation).
  • X+ZZ (parity broken at 0.1): error ≤ disorder_strength × leakage
    ~ 1e-5 ballpark; should grow slowly with t.

Also measures the |+⟩⟨+|^⊗N Arnoldi-seed alternative (qf-e4z.28
Test D extension) at multiple krylovdim values to see if a single-pass
strategy with a different seed converges to 1e-7 at krylovdim ≤ 80.

Coherent term ON throughout (default `include_coherent=true`).

Run:
    julia --project scripts/scratch_qf_e4z_29_trajectory_vs_dense.jl

Cost: ~3 min wall at n=5 (dense exp is the dominant cost).
"""

using LinearAlgebra
using Printf
using BSON
using QuantumFurnace
using QuantumFurnace: X, Y, Z, _parse_hamiltonian_bson, HamHam, beta_alg,
                     _build_jump_set

const FIXTURE_DIR = joinpath(@__DIR__, "output", "multiseed_fixtures")
const OUT_DIR     = joinpath(@__DIR__, "output", "qf_e4z_29")
isdir(OUT_DIR) || mkpath(OUT_DIR)

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

function build_cfg(ham, n, β_phys)
    β_alg = beta_alg(ham, β_phys)
    σ = 1.0 / β_alg
    H_norm = maximum(abs, ham.eigvals)
    omega_range = 2.0 * (H_norm + 8 * σ)
    r_D = 7
    w0_D = omega_range / 2.0^r_D
    t0_D = 2π / (2.0^r_D * w0_D)
    return Config(
        sim = Lindbladian(), domain = EnergyDomain(), construction = KMS(),
        num_qubits = n, with_linear_combination = true,
        beta = β_alg, beta_phys = β_phys, sigma = σ,
        a = 0.0, s = 0.25, gaussian_parameters = (nothing, nothing),
        num_energy_bits_D = r_D, w0_D = w0_D, t0_D = t0_D,
        num_trotter_steps_per_t0 = 10, filter = nothing,
    )
end

function rho_plus_tensor(n)
    plus = ComplexF64[0.5 0.5; 0.5 0.5]
    rho = plus
    for _ in 2:n
        rho = kron(rho, plus)
    end
    return rho
end

"""
    dense_trajectory(L_dense, rho_0, t_grid) -> Vector{Matrix{ComplexF64}}

Reference trajectory at every t in `t_grid`, via dense matrix exponential.
The Lindbladian acts as a d²×d² matrix on `vec(rho)`; exp(L·t)·vec(rho_0)
gives vec(rho(t)), which we reshape back to a d×d operator.

Note: this is O(d^6) per t-grid point via `exp(t·L)` — manageable at d=32
(n=5) where d² = 1024. At d=128 (n=7) the dense L is 16384×16384, too
big for sandbox (would need ~32 GB).
"""
function dense_trajectory(L_dense::Matrix{ComplexF64}, rho_0::Matrix{ComplexF64},
                          t_grid::AbstractVector{<:Real})
    d = size(rho_0, 1)
    states = Vector{Matrix{ComplexF64}}(undef, length(t_grid))
    vec_rho_0 = vec(rho_0)
    for (k, t) in enumerate(t_grid)
        vec_rho_t = exp(L_dense .* t) * vec_rho_0
        states[k] = reshape(vec_rho_t, d, d)
    end
    return states
end

function trace_distance(A::Matrix, B::Matrix)
    diff = A .- B
    diff_h = (diff .+ diff') ./ 2  # Hermitise
    return sum(svdvals(diff_h)) / 2
end

"""
    operator_l1(A, B) -> Float64

Trace-norm distance of two operators (NOT halved). Useful for measuring
predictor-vs-dense error without imposing the trace-distance 'physical
probability' factor of 1/2.
"""
operator_l1(A, B) = sum(svdvals(A .- B))

# -------------------------------------------------------------------------
# Per-fixture / per-rho_0-kind runner
# -------------------------------------------------------------------------

"""
    run_cell(fixture_stem, n, β_phys, seed; rho_0_kind, krylovdim) -> NamedTuple

Build the dense Lindbladian, compute the reference trajectory via
exp(L·t)·vec(rho_0), and compare against `predict_lindbladian_trajectory`'s
reconstructed states[k] at every t-grid point.

Returns a NamedTuple with per-t pointwise errors plus aggregate
diagnostics (gap_traj, gap_ref, total_matvecs).
"""
function run_cell(fixture_stem, n, β_phys, seed;
                  rho_0_kind::Symbol = :identity,
                  krylovdim::Int = 40)
    ham_path = joinpath(FIXTURE_DIR,
        @sprintf("%s_periodic_n%d_seed%d.bson", fixture_stem, n, seed))
    raw = _parse_hamiltonian_bson(ham_path)
    ham = HamHam(raw; beta_phys=β_phys)
    jumps = _build_jump_set(ham, n)
    cfg = build_cfg(ham, n, β_phys)

    d = 2^n
    rho_0 = rho_0_kind === :identity ?
        Matrix{ComplexF64}(I(d) / d) :
        rho_plus_tensor(n)

    # Dense Lindbladian (for the reference trajectory).
    L_dense = construct_lindbladian(jumps, cfg, ham)

    # Reference gap = smallest |Re(λ)| eigenvalue with |Re(λ)| > 1e-12.
    eigvals_dense = eigvals(L_dense)
    re_eigs = real.(eigvals_dense)
    gap_dense = minimum(abs.(re_eigs[abs.(re_eigs) .> 1e-10]))

    # t-grid: range 0..10/gap_dense. 11 points (default Krylov pred range).
    t_max = 10.0 / gap_dense
    t_grid = collect(range(0.0, t_max, length=11))

    # Dense reference trajectory.
    rho_ref = dense_trajectory(L_dense, rho_0, t_grid)

    # Predictor trajectory (current qf-e4z.27 two-pass).
    traj = predict_lindbladian_trajectory(cfg, ham, jumps, rho_0, t_grid;
                                          krylovdim=krylovdim, tol=1e-10,
                                          save_states=true)

    # Pointwise errors (L1 of operator difference, NOT halved).
    err_ops = [operator_l1(traj.states[k], rho_ref[k]) for k in eachindex(t_grid)]
    err_max = maximum(err_ops)
    err_final = err_ops[end]

    # Sanity: trace preservation of reference (should be 1 throughout).
    tr_drift = maximum(abs(real(tr(rho_ref[k])) - 1.0) for k in eachindex(t_grid))

    return (
        t_grid          = t_grid,
        err_ops         = err_ops,
        err_max         = err_max,
        err_final       = err_final,
        gap_traj        = abs(real(traj.eigenvalues[2])),
        gap_2pass       = traj.spectral_gap,
        gap_dense       = gap_dense,
        ref_trace_drift = tr_drift,
    )
end

# -------------------------------------------------------------------------
# Driver
# -------------------------------------------------------------------------

function main()
    println("=" ^ 78)
    println("qf-e4z.29 — Trajectory-vs-dense check for predict_lindbladian_trajectory")
    println("=" ^ 78)

    n = 5
    β_phys = 2.5
    seed = 42

    summary = String[]
    all_data = Dict{Symbol, Any}()

    # --- Test 1: I/d on both fixtures (qf-e4z.27 symmetry-argument check) ---
    println("\n--- Test 1: rho_0 = I/d, qf-e4z.27 two-pass at krylovdim=40 ---")
    println("Checks whether single-pass Arnoldi from vec(I/d), which gives the")
    println("WRONG gap on parity-symmetric L, nonetheless reconstructs the")
    println("trajectory rho(t) = e^{tL}(I/d) correctly.")
    @printf("%-5s %-14s %-14s %-14s %-12s %-12s %-12s\n",
            "fix", "err_max", "err_final", "tr_drift_ref",
            "gap_traj", "gap_2pass", "gap_dense")
    println("-" ^ 95)
    for fix in (:zzdis, :xzdis)
        stem = fix === :zzdis ? "heis_xxx_disordered" : "heis_xxx_XZdisordered"
        r = run_cell(stem, n, β_phys, seed; rho_0_kind=:identity, krylovdim=40)
        all_data[Symbol("test1_", fix)] = r
        @printf("%-5s %-14.4e %-14.4e %-14.4e %-12.4e %-12.4e %-12.4e\n",
            String(fix), r.err_max, r.err_final, r.ref_trace_drift,
            r.gap_traj, r.gap_2pass, r.gap_dense)
        push!(summary, @sprintf("T1 [%s I/d]: err_max=%.2e gap_traj=%.4f gap_dense=%.4f",
            String(fix), r.err_max, r.gap_traj, r.gap_dense))
    end

    # --- Test 2: |+⟩⟨+|^⊗N seed at krylovdim ∈ {40, 60, 80} ---
    # Tests whether single-pass with this non-symmetric seed converges
    # to the dense reference for both trajectory AND gap.
    println("\n--- Test 2: rho_0 = |+⟩⟨+|^⊗N, krylovdim ∈ {40, 60, 80} ---")
    println("With parity-breaking seed, single-pass Arnoldi captures both")
    println("parity sectors. Does it converge to dense reference at larger")
    println("krylovdim, and from what dim?")
    @printf("%-5s %-6s %-14s %-14s %-12s %-12s\n",
            "fix", "kdim", "err_max", "err_final", "gap_traj", "gap_dense")
    println("-" ^ 75)
    for fix in (:zzdis, :xzdis)
        stem = fix === :zzdis ? "heis_xxx_disordered" : "heis_xxx_XZdisordered"
        for kdim in (40, 60, 80)
            r = run_cell(stem, n, β_phys, seed; rho_0_kind=:plus_tensor, krylovdim=kdim)
            all_data[Symbol("test2_", fix, "_kdim", kdim)] = r
            @printf("%-5s %-6d %-14.4e %-14.4e %-12.4e %-12.4e\n",
                String(fix), kdim, r.err_max, r.err_final, r.gap_traj, r.gap_dense)
            push!(summary, @sprintf("T2 [%s |+⟩^N kdim=%d]: err_max=%.2e gap_traj=%.4f gap_dense=%.4f",
                String(fix), kdim, r.err_max, r.gap_traj, r.gap_dense))
        end
    end

    # --- Save full per-t data to BSON for plotting ---
    bson_path = joinpath(OUT_DIR, @sprintf("trajectory_vs_dense_n%d_beta%.2f_seed%d.bson",
                                            n, β_phys, seed))
    BSON.bson(bson_path, all_data)
    println("\nFull per-t error arrays saved to: ", bson_path)

    # --- Summary ---
    println("\n" * "=" ^ 78)
    println("Summary:")
    println("=" ^ 78)
    for line in summary
        println("  ", line)
    end

    # Interpretive verdict.
    println("\n--- Verdict ---")
    t1_zz_err = all_data[:test1_zzdis].err_max
    t1_xz_err = all_data[:test1_xzdis].err_max
    println(@sprintf("Test 1 (I/d): Z+ZZ err_max = %.2e, X+ZZ err_max = %.2e",
                     t1_zz_err, t1_xz_err))
    if t1_zz_err < 1e-7
        println("✓ qf-e4z.27 trajectory claim VERIFIED on Z+ZZ: trajectory is correct")
        println("  to 1e-7 despite the 17 % gap reporting error.")
    elseif t1_zz_err < 1e-4
        println("≈ qf-e4z.27 trajectory claim partially verified on Z+ZZ: trajectory")
        println(@sprintf("  is correct to %.0e (Krylov-truncation regime).", t1_zz_err))
    else
        println("✗ qf-e4z.27 trajectory claim FAILS on Z+ZZ: trajectory error")
        println(@sprintf("  reaches %.2e — the symmetry argument is wrong or", t1_zz_err))
        println("  Krylov truncation is much worse than expected. INVESTIGATE.")
    end

    # Test 2 convergence check.
    println()
    for fix in (:zzdis, :xzdis)
        e40 = all_data[Symbol("test2_", fix, "_kdim40")].err_max
        e80 = all_data[Symbol("test2_", fix, "_kdim80")].err_max
        ratio = e40 / max(e80, eps())
        println(@sprintf("Test 2 |+⟩ on %s: err_max kdim=40 → 80 : %.2e → %.2e  (×%.2f)",
                          String(fix), e40, e80, ratio))
    end
end

main()
