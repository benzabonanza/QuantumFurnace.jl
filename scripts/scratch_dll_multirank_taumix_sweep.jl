#!/usr/bin/env julia
#
# Multi-rank DLL τ_mix sweep (qf-7go.6).
#
# Headline question
# -----------------
# Does adding rank to the DLL Lindbladian — k > 1 symmetrised-translate
# channels per coupling — close the τ_mix gap to CKG smooth-Metropolis,
# and where does τ_mix saturate as k grows?
#
# From `drafts/ckg-vs-dll-comparison-findings.md` (memory:
# `ckg_vs_dll_first_findings.md` + `fair_comparison_dirichlet_qf_mto.md`):
# - Standard rank-1 DLL Metropolis tracks CKG smooth-Metropolis to ±10–25%
#   in τ_mix at n=3, β ≤ 10.
# - DLL Gaussian collapses 8× at n=3, β=20 (single narrow channel; no
#   parallelism).
# - Dirichlet-form analysis (qf-mto) confirmed Metropolis-shape α_diagonal
#   is the τ_mix driver (H2 hypothesis confirmed); rank had no effect on
#   the diagonal at the comparison cells. Predicts: multi-rank DLL Metro
#   gives null result on τ_mix at moderate (n, β); the test is whether the
#   DLL Gaussian collapse at large β is recovered by adding channels.
#
# Sweep grid
# ----------
#   n      ∈ {3}          (only n=3 here; n≥4 deferred to follow-up)
#   β      ∈ {1, 5, 10, 20}
#   k      ∈ {1, 2, 4, 8} (k=1 reproduces single-channel baseline)
#   filter ∈ {DLLMetropolis, DLLGaussian} (both base shapes)
#
# Centers for k > 1: uniform spacing on [0, c_max] with c_max ≤ S/2 = 1.0
# (Metropolis bump flat-top constraint, validated by
# `dll_multichannel_translates`). For Gaussian (no compact support), use
# the same centers — the per-channel cosh envelope only grows by O(1) at
# the largest center, so the construction stays well-behaved.
#
# PHYSICS CHECK: target_ε = 1e-3; t_max_factor = :auto (qf-lkb.10
# adaptive). Smooth-Metropolis a/s knobs are inert for DLL filters.
#
# Output
# ------
# - BSON file at `scripts/output/dll_multirank/taumix_vs_k.bson` with
#   per-cell (β, k, base_shape) τ_mix, fitted gap, etc.
# - Summary table printed to stdout (k → τ_mix at each β).
#
# Acceptance
# ----------
# - All 32 cells (4 β × 4 k × 2 shapes) complete without error.
# - k=1 column matches the existing single-channel DLL τ_mix to ≤1e-3
#   relative on the same fixture (sanity).
# - Plot/table shows τ_mix vs k at each (β, shape).

using Printf
using Random
using LinearAlgebra
using BSON
using QuantumFurnace
using QuantumFurnace: _load_hamiltonian_bson, _build_jump_set, _make_init_state

const N_QUBITS    = 3
const BETA_VALUES = [1.0, 5.0, 10.0, 20.0]
const K_VALUES    = [1, 2, 4, 8]
const BASE_SHAPES = (:metropolis, :gaussian)
const SEED        = 42
const TARGET_EPS  = 1e-3
const T_GRID_LEN  = 81
const KRYLOV_DIM  = 30
const TOL         = 1e-10
const T_MAX_FACTOR = max(5.0, 1.5 * log10(1.0 / TARGET_EPS))  # qf-lkb.10

const OUT_DIR  = joinpath(@__DIR__, "output", "dll_multirank")
const BSON_OUT = joinpath(OUT_DIR, "taumix_vs_k.bson")
mkpath(OUT_DIR)

# Centers for k channels, uniform on [0, c_max]. k=1 ⇒ [0.0] (single
# channel, equivalent to the standard rank-1 DLL filter).
function _centers_for_k(k::Int, c_max::Float64)
    k == 1 && return [0.0]
    return collect(range(0.0, c_max; length = k))
end

function _build_filter(base_shape::Symbol, beta::Float64, k::Int; S::Float64 = 2.0)
    if base_shape === :metropolis
        base = DLLMetropolisFilter(beta; S = S)
        c_max = S / 2 - 1e-9   # stay strictly inside the Metropolis flat top
        return dll_multichannel_translates(base; centers = _centers_for_k(k, c_max))
    elseif base_shape === :gaussian
        base = DLLGaussianFilter(beta)
        # Gaussian has no compact support; pick c_max comparable to the
        # Metropolis flat top so both shapes share the same center grid.
        c_max = 1.0
        return dll_multichannel_translates(base; centers = _centers_for_k(k, c_max))
    else
        throw(ArgumentError("base_shape must be :metropolis or :gaussian (got $base_shape)"))
    end
end

# Run one sweep cell. Returns NamedTuple with τ_mix + diagnostics.
function _run_cell(ham, jumps, beta::Float64, k::Int, base_shape::Symbol)
    filter = _build_filter(base_shape, beta, k)
    config = Config(;
        sim                       = Lindbladian(),
        domain                    = BohrDomain(),
        construction              = DLL(),
        num_qubits                = N_QUBITS,
        with_linear_combination   = true,
        beta                      = beta,
        sigma                     = 1.0 / beta,
        a                         = 0.0,
        s                         = 0.25,
        num_energy_bits           = 12,
        w0                        = 0.05,
        t0                        = 2π / (2^12 * 0.05),
        num_trotter_steps_per_t0  = 10,
        filter                    = filter,
    )

    # Spectral-gap pre-pass via matrix-free Krylov.
    gap_est = try
        kr = krylov_spectral_gap(config, ham, jumps;
                                  krylovdim = 30, howmany = 2, tol = 1e-8)
        kr.spectral_gap > 0 ? kr.spectral_gap : 1.0 / beta
    catch err
        @warn "krylov_spectral_gap failed; falling back to 1/β" beta k base_shape err
        1.0 / beta
    end
    t_max = T_MAX_FACTOR / max(gap_est, 1e-12)
    t_grid = collect(range(0.0, t_max; length = T_GRID_LEN))

    rho_0 = _make_init_state(:maximally_mixed,
                             size(ham.data, 1),
                             Matrix{ComplexF64}(ham.gibbs),
                             SEED)

    t0 = time()
    res = integrate_to_gibbs(config, ham, jumps, rho_0, t_grid;
                              mode = :L, krylovdim = KRYLOV_DIM, tol = TOL)
    est = estimate_mixing_time(res; model = :biexp,
                                target_epsilon = TARGET_EPS,
                                extrapolate = true)
    wall = time() - t0

    mt_source = if isfinite(est.mixing_time)
        :extrapolated
    elseif est.mixing_time_actual !== nothing && isfinite(est.mixing_time_actual)
        :observed
    else
        :nan
    end
    mt = if mt_source === :extrapolated
        est.mixing_time
    elseif mt_source === :observed
        est.mixing_time_actual::Float64
    else
        NaN
    end

    return (
        beta              = beta,
        k                 = k,
        base_shape        = base_shape,
        gap_est           = gap_est,
        t_max             = t_max,
        n_grid            = length(t_grid),
        fitted_gap        = est.fitted_gap,
        mixing_time       = mt,
        mixing_time_source = mt_source,
        r_squared         = est.r_squared,
        converged_fit     = est.converged,
        wall_time         = wall,
    )
end

# ── Driver ────────────────────────────────────────────────────────────────────
println("="^72)
println("DLL multi-rank τ_mix sweep (qf-7go.6)")
println("="^72)
@printf("n            : %d\n", N_QUBITS)
@printf("beta_values  : %s\n", string(BETA_VALUES))
@printf("k_values     : %s\n", string(K_VALUES))
@printf("base_shapes  : %s\n", string(BASE_SHAPES))
@printf("target ε     : %.1e\n", TARGET_EPS)
@printf("t_max_factor : %.2f (qf-lkb.10 derived)\n", T_MAX_FACTOR)
@printf("Julia threads: %d, BLAS threads: %d\n",
        Threads.nthreads(), BLAS.get_num_threads())
println("="^72)
flush(stdout)

ham_path = joinpath(dirname(@__DIR__), "hamiltonians",
                     "heis_disordered_periodic_n$(N_QUBITS).bson")
isfile(ham_path) || error("Hamiltonian fixture missing: $ham_path")
ham = _load_hamiltonian_bson(ham_path, BETA_VALUES[1])  # placeholder; rebuilt per β
jumps_template = _build_jump_set(ham, N_QUBITS)

results = NamedTuple[]
t_total = time()
for shape in BASE_SHAPES
    for beta in BETA_VALUES
        # Reload hamiltonian at the right β (the .gibbs cache depends on β).
        ham_β = _load_hamiltonian_bson(ham_path, beta)
        jumps_β = _build_jump_set(ham_β, N_QUBITS)
        for k in K_VALUES
            t_cell0 = time()
            res = _run_cell(ham_β, jumps_β, beta, k, shape)
            push!(results, res)
            @printf("[%s] β=%-5.1f k=%-2d  τ_mix=%-9.3g  src=%-12s  fit_gap=%-9.3g  wall=%.1fs\n",
                    string(shape), beta, k, res.mixing_time,
                    string(res.mixing_time_source), res.fitted_gap, res.wall_time)
            flush(stdout)
        end
    end
end

@printf("\nTotal wall time: %.1fs\n", time() - t_total)

# Persist BSON for the plotter (qf-7go.7).
bson_payload = Dict{Symbol, Any}(
    :results       => results,
    :n             => N_QUBITS,
    :beta_values   => BETA_VALUES,
    :k_values      => K_VALUES,
    :base_shapes   => collect(BASE_SHAPES),
    :target_epsilon => TARGET_EPS,
    :t_max_factor  => T_MAX_FACTOR,
    :seed          => SEED,
    :git_sha       => try
        read(`git -C $(dirname(@__DIR__)) rev-parse HEAD`, String) |> strip
    catch
        "unknown"
    end,
)
BSON.bson(BSON_OUT, bson_payload)
@printf("\nWrote BSON: %s\n", BSON_OUT)

# ── Summary table (k vs τ_mix at each (β, shape)) ────────────────────────────
println("\n" * "="^72)
println("SUMMARY — τ_mix vs k")
println("="^72)
for shape in BASE_SHAPES
    @printf("\n[base = %s]\n", shape)
    @printf("  β        ")
    for k in K_VALUES
        @printf("k=%-9d ", k)
    end
    println()
    @printf("  ---------")
    for _ in K_VALUES
        @printf("-----------")
    end
    println()
    for beta in BETA_VALUES
        @printf("  %-9.1f", beta)
        for k in K_VALUES
            row = findfirst(r -> r.beta == beta && r.k == k && r.base_shape == shape,
                            results)
            mt = row === nothing ? NaN : results[row].mixing_time
            @printf("%-11.3g", mt)
        end
        println()
    end
end
println()
