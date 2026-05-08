#=
Benchmark L-mode (Arnoldi on ρ) vs K-mode (Lanczos on ψ) integrator paths
for the KMS BohrDomain Lindbladian, sweeping n ∈ {3, 4} and β ∈ {1, 5, 10}.

For each (n, β) cell we:
  1. Build the disordered Heisenberg system (load BSON fixture, attach
     XYZ-on-each-site jumps normalised by sqrt(num_jumps)).
  2. Construct the dense KMS Lindbladian once and pass it through
     `discriminant_spectrum` to read `H_gap` (the parent-Hamiltonian gap
     of the discriminant; KMS-DB asymptotic decay rate of both ρ and ψ).
  3. Set t_max = T_MAX_FACTOR / H_gap and t_grid as a linspace of length N_GRID.
  4. Run both modes via `integrate_to_gibbs` with identical config/jumps/rho_0.
  5. Fit bi-exponential decay via `estimate_mixing_time` to extract τ_mix
     (extrapolated to target_epsilon) and the slow `fitted_gap`.
  6. Record wall time (time_ns deltas), matvec count, R², convergence flag.

A warm-up call at n=3, β=1, mode=:L compiles the integrator method body
before we read time_ns; per-cell timings exclude compilation.

Why DLL is excluded here: at n=4 the dense L-superoperator is 256² ComplexF64
≈ 1 MB (fine) but DLL BohrDomain still goes through the dense
construct_lindbladian path, while CKG uses a matrix-free Workspace closure.
Comparing modes head-to-head is cleanest on the matrix-free CKG path; DLL
mode comparison is a separate question for qf-lkb later.

PHYSICS CHECK summary:
- T_MAX_FACTOR = 5.0 → exp(-5) ≈ 0.007 of slow-mode amplitude remains; gives
  bi-exp fit enough decay to push offset C below target_eps for extrapolation.
- target_epsilon: TARGET_EPS_L = TARGET_EPS_K = 1e-3. The K integrator returns
  ‖ψ−ψ_eq‖_F = √χ², which is the natural Lyapunov function in the chi-metric
  representation and is metric-equivalent to the trace distance up to O(1)
  factors. Both modes therefore target the same Gibbs-proximity scale, and
  both bi-exp **slow gaps** should converge to the same value ≈ H_gap (the
  spectral gap of the parent Hamiltonian H_part = D for KMS-DB).
- BLAS oversubscription mitigated by serial outer loop (no @threads); BLAS
  internal threading kept at the system default.

Run: julia --project scripts/benchmark_lmode_vs_kmode.jl
=#

using QuantumFurnace
using LinearAlgebra
using LinearAlgebra: BLAS
using Printf
using BSON
using Dates

include(joinpath(@__DIR__, "..", "test", "test_helpers.jl"))

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

const NS              = [3, 4]
const BETAS           = [1.0, 5.0, 10.0]
const KRYLOVDIM       = 30
const TOL             = 1e-10
const N_GRID          = 81
# PHYSICS CHECK: 5/H_gap ≈ time at which slow mode is ~0.7% of initial → bi-exp
# fit offset C drops below tight target_eps (esp. K-mode 1e-6) for extrapolation.
const T_MAX_FACTOR    = 5.0
# PHYSICS CHECK: TD target for L-mode trajectories (matches existing tests).
const TARGET_EPS_L    = 1e-3
# PHYSICS CHECK: K integrator returns ‖ψ−ψ_eq‖_F = √χ², which is the natural
# Lyapunov function in the chi-metric representation and is metric-equivalent
# to the trace distance up to O(1) factors. Targeting `√χ² = 1e-3` matches
# the same Gibbs-proximity scale as TARGET_EPS_L = 1e-3 (corresponds to
# χ² = 1e-6) and lets the bi-exp extrapolator complete in the same horizon
# range as L-mode. The earlier choice 1e-6 — interpreting the integrator
# output as χ² rather than √χ² — was tighter than 1e-3 by an extra factor
# of 1e-3 squared and caused the n=4, β ∈ {5,10} cells to NaN out in the
# pre-fix benchmark run.
const TARGET_EPS_K    = 1e-3
const OUTPUT_PATH     = joinpath(@__DIR__, "output", "benchmark_lmode_vs_kmode.bson")
const HAM_DIR         = joinpath(@__DIR__, "..", "hamiltonians")
const HAM_PATH(n)     = joinpath(HAM_DIR, "heis_disordered_periodic_n$n.bson")

# ---------------------------------------------------------------------------
# System builders
# ---------------------------------------------------------------------------

"""
    build_n3_system(beta) -> (; ham, jumps, gibbs)

Delegate to the test fixture factory used by all DLL/CKG tests at n=3.
"""
build_n3_system(beta::Real) = make_dll_n3_system(beta)

"""
    build_n4_system(beta) -> (; ham, jumps, gibbs) | nothing

Load the n=4 cached Hamiltonian (if present) and attach the same XYZ-on-each-
site jump set used at n=3 (3 paulis × 4 sites = 12 jumps, normalised by √12).
Returns `nothing` and emits an `@warn` if the BSON fixture is missing.
"""
function build_n4_system(beta::Real)
    path = HAM_PATH(4)
    if !isfile(path)
        @warn "Missing fixture; skipping n=4 cells" path
        return nothing
    end
    ham = QuantumFurnace._load_hamiltonian_bson(path, Float64(beta))
    jump_paulis = [[X], [Y], [Z]]
    n = 4
    num_jumps = length(jump_paulis) * n
    jump_norm = sqrt(num_jumps)
    jumps = JumpOp[]
    for pauli in jump_paulis
        for site in 1:n
            op = Matrix(pad_term(pauli, n, site)) ./ jump_norm
            op_eb = ham.eigvecs' * op * ham.eigvecs
            push!(jumps, JumpOp(op, op_eb, op == transpose(op), op == op'))
        end
    end
    return (; ham, jumps, gibbs = ham.gibbs)
end

"""
    build_system(n, beta) -> (; ham, jumps, gibbs) | nothing
"""
function build_system(n::Int, beta::Real)
    n == 3 && return build_n3_system(beta)
    n == 4 && return build_n4_system(beta)
    error("Unsupported n=$n (expected 3 or 4)")
end

# ---------------------------------------------------------------------------
# Per-cell config + run
# ---------------------------------------------------------------------------

"""
    make_kms_bohr_config(n, beta) -> Config

KMS BohrDomain config; the trotter/time grid kwargs are required by the
`Config` validator but ignored by BohrDomain. Mirrors the parameters used in
test_helpers `make_config(Lindbladian(), BohrDomain(); construction=KMS())`.
"""
function make_kms_bohr_config(n::Int, beta::Real)
    return Config(;
        sim                      = Lindbladian(),
        domain                   = BohrDomain(),
        construction             = KMS(),
        num_qubits               = n,
        with_linear_combination  = true,
        beta                     = float(beta),
        sigma                    = 1.0 / float(beta),
        a                        = float(beta) / 30.0,
        s                        = 0.4,
        num_energy_bits          = 12,
        w0                       = 0.05,
        t0                       = 2π / (2^12 * 0.05),
        num_trotter_steps_per_t0 = 10,
    )
end

"""
    compute_H_gap(config, ham, jumps) -> Float64

Build the dense KMS Lindbladian once and read the parent-Hamiltonian gap
from its discriminant spectrum. Note: `discriminant_spectrum` requires a
regular `Matrix`, not a `SparseMatrixCSC`, so we call `Matrix(...)` defensively.
"""
function compute_H_gap(config::Config, ham, jumps::Vector{JumpOp})
    L_dense = Matrix(construct_lindbladian(jumps, config, ham))
    spec = discriminant_spectrum(L_dense, ham.gibbs; n_modes = 4)
    return spec.H_gap
end

"""
    run_cell(n, beta, mode, ham, jumps, H_gap) -> NamedTuple

Build t_grid from H_gap, run `integrate_to_gibbs` with the given mode,
fit bi-exp, return the per-cell metrics record.
"""
function run_cell(n::Int, beta::Real, mode::Symbol,
                  ham, jumps::Vector{JumpOp}, H_gap::Float64,
                  t_max_factor::Float64 = T_MAX_FACTOR,
                  n_grid::Int = N_GRID)
    config = make_kms_bohr_config(n, beta)
    t_max = t_max_factor / H_gap
    t_grid = collect(range(0.0, t_max; length = n_grid))
    d = size(ham.data, 1)
    rho_0 = Matrix{ComplexF64}(I(d) / d)

    target_eps = mode == :L ? TARGET_EPS_L : TARGET_EPS_K

    t0_ns = time_ns()
    res = integrate_to_gibbs(config, ham, jumps, rho_0, t_grid;
                             mode = mode, krylovdim = KRYLOVDIM, tol = TOL)
    wall_s = (time_ns() - t0_ns) / 1e9

    est = estimate_mixing_time(res; model = :biexp,
                               target_epsilon = target_eps,
                               extrapolate = true)

    return (
        n              = n,
        beta           = float(beta),
        mode           = mode,
        H_gap          = H_gap,
        t_max          = t_max,
        n_grid         = n_grid,
        wall_s         = wall_s,
        matvecs        = res.total_matvecs,
        all_converged  = res.all_converged,
        dist_end       = res.distances[end],
        fitted_gap     = est.fitted_gap,
        r_squared      = est.r_squared,
        tau_mix        = est.mixing_time,
        converged_fit  = est.converged,
        target_eps     = target_eps,
        failed         = false,
        error_message  = "",
    )
end

"""
    run_cell_safe(...) -> NamedTuple

Wrap `run_cell` so a single cell failure (e.g. fit divergence, OOM) does
not abort the sweep. On error returns a placeholder NamedTuple with NaN
metrics + the exception message.
"""
function run_cell_safe(n::Int, beta::Real, mode::Symbol,
                       ham, jumps::Vector{JumpOp}, H_gap::Float64;
                       t_max_factor::Float64 = T_MAX_FACTOR,
                       n_grid::Int = N_GRID)
    try
        return run_cell(n, beta, mode, ham, jumps, H_gap, t_max_factor, n_grid)
    catch err
        @warn "Cell failed" n beta mode error=err
        return (
            n              = n,
            beta           = float(beta),
            mode           = mode,
            H_gap          = H_gap,
            t_max          = t_max_factor / H_gap,
            n_grid         = n_grid,
            wall_s         = NaN,
            matvecs        = 0,
            all_converged  = false,
            dist_end       = NaN,
            fitted_gap     = NaN,
            r_squared      = NaN,
            tau_mix        = NaN,
            converged_fit  = false,
            target_eps     = mode == :L ? TARGET_EPS_L : TARGET_EPS_K,
            failed         = true,
            error_message  = sprint(showerror, err),
        )
    end
end

# ---------------------------------------------------------------------------
# Pretty-print helpers
# ---------------------------------------------------------------------------

function print_table(rows::Vector{<:NamedTuple})
    @printf "\n%-3s %-5s %-4s %-9s %-6s %-8s %-9s %-7s %-9s %-12s %-9s\n" "n" "β" "mode" "t_max" "Ngrid" "wall(s)" "matvecs" "R²" "τ_mix" "fitted_gap" "H_gap"
    println(repeat('-', 95))
    for r in rows
        if r.failed
            @printf "%-3d %-5.1f %-4s %-9.2f %-6d %-8s %-9s %-7s %-9s %-12s %-9.3e\n" r.n r.beta string(r.mode) r.t_max r.n_grid "FAIL" "FAIL" "FAIL" "FAIL" "FAIL" r.H_gap
        else
            @printf "%-3d %-5.1f %-4s %-9.2f %-6d %-8.3f %-9d %-7.3f %-9.3f %-12.3e %-9.3e\n" r.n r.beta string(r.mode) r.t_max r.n_grid r.wall_s r.matvecs r.r_squared r.tau_mix r.fitted_gap r.H_gap
        end
    end
    println(repeat('-', 95))
end

# ---------------------------------------------------------------------------
# Summary statistics + recommendation rule
# ---------------------------------------------------------------------------

_safemedian(xs) = isempty(xs) ? NaN : begin
    cleaned = filter(x -> isfinite(x), xs)
    isempty(cleaned) ? NaN : median(cleaned)
end

# Bare-bones median (avoid pulling in Statistics dep just for one call).
function median(xs::AbstractVector{<:Real})
    sorted = sort(collect(xs))
    n = length(sorted)
    n == 0 && return NaN
    isodd(n) ? float(sorted[(n + 1) ÷ 2]) : (sorted[n ÷ 2] + sorted[n ÷ 2 + 1]) / 2
end

function recommend(rows::Vector{<:NamedTuple})
    L_rows = filter(r -> r.mode == :L && !r.failed, rows)
    K_rows = filter(r -> r.mode == :K && !r.failed, rows)
    isempty(L_rows) && return ("tied", "no L-mode cells succeeded")
    isempty(K_rows) && return ("tied", "no K-mode cells succeeded")

    med_wall_L = _safemedian([r.wall_s for r in L_rows])
    med_wall_K = _safemedian([r.wall_s for r in K_rows])
    med_mv_L   = _safemedian([float(r.matvecs) for r in L_rows])
    med_mv_K   = _safemedian([float(r.matvecs) for r in K_rows])
    med_r2_L   = _safemedian([r.r_squared for r in L_rows])
    med_r2_K   = _safemedian([r.r_squared for r in K_rows])

    rec, reason = if med_wall_K < 0.9 * med_wall_L && med_r2_K >= med_r2_L - 0.01
        ("K", @sprintf("med(wall_K)=%.3fs vs med(wall_L)=%.3fs (≥10%% faster), R² parity (%.3f vs %.3f)", med_wall_K, med_wall_L, med_r2_K, med_r2_L))
    elseif med_wall_L < 0.9 * med_wall_K && med_r2_L >= med_r2_K - 0.01
        ("L", @sprintf("med(wall_L)=%.3fs vs med(wall_K)=%.3fs (≥10%% faster), R² parity (%.3f vs %.3f)", med_wall_L, med_wall_K, med_r2_L, med_r2_K))
    else
        ("tied", @sprintf("walls within 10%% (L=%.3fs, K=%.3fs); R²: L=%.3f, K=%.3f", med_wall_L, med_wall_K, med_r2_L, med_r2_K))
    end

    return (rec, reason, (; med_wall_L, med_wall_K, med_mv_L, med_mv_K, med_r2_L, med_r2_K))
end

# ---------------------------------------------------------------------------
# Main sweep
# ---------------------------------------------------------------------------

function main()
    println("# Benchmark L-mode vs K-discriminant integrator paths")
    println("# qf-lkb.7 — KMS BohrDomain only (matrix-free CKG path)")
    @printf "# Julia threads: %d   BLAS threads: %d   Date: %s\n" Threads.nthreads() BLAS.get_num_threads() Dates.now()
    println("#")
    @printf "# Sweep: n ∈ %s, β ∈ %s, mode ∈ {:L, :K}\n" repr(NS) repr(BETAS)
    @printf "# t_max = %.2f / H_gap, n_grid = %d, krylovdim = %d, tol = %g\n" T_MAX_FACTOR N_GRID KRYLOVDIM TOL
    @printf "# target_eps_L = %g, target_eps_K = %g\n" TARGET_EPS_L TARGET_EPS_K
    println()

    t_total_0 = time_ns()

    # Warm-up at n=3, β=1, mode=:L to compile method specialisations
    # (so per-cell `time_ns` excludes integrator compile time).
    println("Warm-up: n=3, β=1, mode=:L (3-point t_grid) ...")
    sys_warm = build_n3_system(1.0)
    cfg_warm = make_kms_bohr_config(3, 1.0)
    rho0_warm = Matrix{ComplexF64}(I(8) / 8)
    _ = integrate_to_gibbs(cfg_warm, sys_warm.ham, sys_warm.jumps, rho0_warm,
                           collect(range(0.0, 1.0; length = 3));
                           mode = :L, krylovdim = KRYLOVDIM, tol = TOL)
    _ = integrate_to_gibbs(cfg_warm, sys_warm.ham, sys_warm.jumps, rho0_warm,
                           collect(range(0.0, 1.0; length = 3));
                           mode = :K, krylovdim = KRYLOVDIM, tol = TOL)
    println("Warm-up done.")

    rows = NamedTuple[]
    for n in NS
        for beta in BETAS
            sys = build_system(n, beta)
            sys === nothing && continue  # n=4 missing fixture: skip both modes

            cfg = make_kms_bohr_config(n, beta)
            print("Computing H_gap @ n=$n, β=$beta ... ")
            flush(stdout)
            t_g0 = time_ns()
            H_gap = compute_H_gap(cfg, sys.ham, sys.jumps)
            t_g1 = (time_ns() - t_g0) / 1e9
            @printf "H_gap = %.3e (%.2fs)\n" H_gap t_g1
            flush(stdout)

            for mode in (:L, :K)
                print("  cell n=$n β=$beta mode=$mode ... ")
                flush(stdout)
                row = run_cell_safe(n, beta, mode, sys.ham, sys.jumps, H_gap)
                if row.failed
                    @printf "FAILED (%s)\n" row.error_message
                else
                    @printf "wall=%.3fs matvecs=%d τ_mix=%.3f gap=%.3e R²=%.3f\n" row.wall_s row.matvecs row.tau_mix row.fitted_gap row.r_squared
                end
                flush(stdout)
                push!(rows, row)
            end
        end
    end

    print_table(rows)

    # Per-cell L vs K fitted_gap parity check.
    # PHYSICS CHECK: For KMS-DB, K = D is HS-self-adjoint with the same spectrum
    # as the parent Hamiltonian H_part (similarity transform). The slowest
    # operator-space mode decays at rate λ_1 = H_gap. Then ‖ψ−ψ_eq‖_F decays
    # at rate H_gap, and ‖ρ−σ‖_TD also decays at rate H_gap (both inherit the
    # spectral gap of D). So gap_L ≈ gap_K ≈ H_gap is the right expectation;
    # the bi-exp fit may not extract the asymptotic slope perfectly when
    # multiple modes contribute non-trivially — print all three and flag
    # ratios well outside ±20% as suspicious.
    println("\nL vs K fitted_gap parity (KMS-DB asymptotics):")
    println("  Expectation: gap_L ≈ gap_K ≈ H_gap (both modes inherit the spectral gap of D = H_part)")
    @printf "  %-3s %-5s %-12s %-12s %-12s %-12s %-10s\n" "n" "β" "gap_L" "gap_K" "H_gap" "gap_K/gap_L" "match?"
    println(repeat('-', 80))
    n_match = 0
    n_pairs = 0
    for n in NS
        for beta in BETAS
            L = filter(r -> r.n == n && r.beta == beta && r.mode == :L && !r.failed, rows)
            K = filter(r -> r.n == n && r.beta == beta && r.mode == :K && !r.failed, rows)
            (isempty(L) || isempty(K)) && continue
            n_pairs += 1
            gL, gK, Hg = L[1].fitted_gap, K[1].fitted_gap, L[1].H_gap
            ratio_L = abs(gL - Hg) / Hg
            ratio_K = abs(gK - Hg) / Hg
            match = ratio_L < 0.20 && ratio_K < 0.20
            match && (n_match += 1)
            @printf "  %-3d %-5.1f %-12.3e %-12.3e %-12.3e %-12.3f %-10s\n" n beta gL gK Hg gK / gL (match ? "yes" : "NO")
        end
    end
    @printf "  parity matches (within 20%%): %d / %d\n" n_match n_pairs

    # Recommendation.
    rec_tuple = recommend(rows)
    rec, reason, summary = rec_tuple
    println("\n=== Recommendation ===")
    @printf "  Use mode: %s\n" rec
    @printf "  Reason  : %s\n" reason
    @printf "  Summary : med wall L=%.3fs / K=%.3fs; med matvecs L=%.0f / K=%.0f; med R² L=%.3f / K=%.3f\n" summary.med_wall_L summary.med_wall_K summary.med_mv_L summary.med_mv_K summary.med_r2_L summary.med_r2_K

    # Save BSON sidecar (run artefact, not committed).
    mkpath(dirname(OUTPUT_PATH))
    BSON.bson(OUTPUT_PATH, Dict(
        :rows     => rows,
        :summary  => summary,
        :recommendation => rec,
        :reason   => reason,
        :metadata => (
            julia_threads = Threads.nthreads(),
            blas_threads  = BLAS.get_num_threads(),
            timestamp_iso = string(Dates.now()),
            ns            = NS,
            betas         = BETAS,
            n_grid        = N_GRID,
            t_max_factor  = T_MAX_FACTOR,
            krylovdim     = KRYLOVDIM,
            tol           = TOL,
            target_eps_L  = TARGET_EPS_L,
            target_eps_K  = TARGET_EPS_K,
        ),
    ))
    @printf "\nSidecar saved: %s\n" OUTPUT_PATH

    t_total = (time_ns() - t_total_0) / 1e9
    @printf "Total wall time: %.1fs\n" t_total

    # Convergence summary
    n_total = length(rows)
    n_ok = count(r -> !r.failed && r.all_converged && isfinite(r.tau_mix), rows)
    @printf "Cells converged: %d / %d\n" n_ok n_total

    return rows
end

main()
