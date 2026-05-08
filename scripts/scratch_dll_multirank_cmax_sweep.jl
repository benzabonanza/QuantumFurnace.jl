#!/usr/bin/env julia
#
# qf-9ld follow-up — multi-rank DLL τ_mix and ‖L‖_be vs c_max.
#
# What this answers
# -----------------
# qf-7go.6 chose c_max ≈ 1.0 (Metropolis bump flat-top constraint), but the
# rescaled fixture has ν_max(H) = 0.45 — so off-center channels were
# placing Kossakowski mass at frequencies outside the Hamiltonian's actual
# Bohr support. The cosh(β·ν_ℓ/4)² blow-up of the LCU 1-norm penalises
# this catastrophically (see drafts/ckg-vs-dll-comparison-findings.md §9.7).
#
# Question: does the qf-7go.6 τ_mix speedup persist when we constrain
# c_max ≤ ν_max(H)? And does the algorithmic figure of merit
# ‖L‖_be · τ_mix at any (k, c_max) cell beat the k=1 baseline?
#
# Sweep
# -----
#   n          = 3  (rescaled fixture; ν_max(H) = 0.45)
#   β          ∈ {5, 10, 20}
#   c_max      ∈ {0.10, 0.20, 0.45, 1.00}  (last is the original qf-7go.6)
#   k          ∈ {1, 2, 4, 8}
#   base_shape = :metropolis (Gaussian inherits the same logic; skip for runtime)
#
# Per cell we compute:
#   • τ_mix (matrix-free integrator + bi-exp extrapolation, target ε=1e-3)
#   • Z_f, Z_g (LCU 1-norms of time-domain weights)
#   • α_G, α_L (block-encoding factors per §9.3)
#   • ‖L‖_be   (LW23 aggregator, Eq. dll-Lbe)
#   • Product ‖L‖_be · τ_mix (the algorithmic figure of merit)
#
# k=1 cell is used as the baseline; ratios reported relative to it.

using Printf
using Random
using LinearAlgebra
using BSON
using QuantumFurnace
using QuantumFurnace: _load_hamiltonian_bson, _build_jump_set, _make_init_state

const N_QUBITS    = 3
const BETA_VALUES = [5.0, 10.0, 20.0]
const C_MAX_VALUES = [0.10, 0.20, 0.45, 1.0 - 1e-9]
const K_VALUES    = [1, 2, 4, 8]
const SEED        = 42
const TARGET_EPS  = 1e-3
const T_GRID_LEN  = 81
const KRYLOV_DIM  = 30
const TOL         = 1e-10
const T_MAX_FACTOR = max(5.0, 1.5 * log10(1.0 / TARGET_EPS))   # qf-lkb.10

# LCU-norm time grid (matches scripts/scratch_dll_norm_audit.jl).
const TAU_LCU = 0.2
const TMAX_FAC = 6.0
const NU_RADIUS_GRID = 192   # Nν for the (ν,ν') quadrature in Z_g

const OUT_DIR  = joinpath(@__DIR__, "output", "dll_multirank")
const BSON_OUT = joinpath(OUT_DIR, "taumix_vs_cmax.bson")
mkpath(OUT_DIR)

# ── Filter construction (matches qf-7go.6 except c_max is a parameter) ────────
function _centers_for_k(k::Int, c_max::Float64)
    k == 1 && return [0.0]
    return collect(range(0.0, c_max; length = k))
end

function _build_filter(beta::Float64, k::Int, c_max::Float64; S::Float64 = 2.0)
    base = DLLMetropolisFilter(beta; S = S)
    centers = _centers_for_k(k, c_max)
    return k == 1 ? base : dll_multichannel_translates(base; centers = centers)
end

# ── LCU-norm computation (per-channel sum for multi-channel) ──────────────────
function compute_Z_f(filter::AbstractFilter, β::Real)
    Tmax = TMAX_FAC * β
    Nt = 2 * round(Int, Tmax / TAU_LCU) + 1
    ts = range(-Tmax, Tmax; length = Nt)
    τ = step(ts)
    return sum(abs(time_kernel(filter, t)) * τ for t in ts)
end

# Single-channel Z_g via 2D-FT of ĝ(ν, ν') on a uniform (ν, ν') grid.
function _Z_g_single(filter, β; Nν = NU_RADIUS_GRID)
    ν_radius = max(2.0, 4.0 / β)   # Metropolis S=2 ⇒ bump on [-2, 2]
    νs = range(-ν_radius, ν_radius; length = Nν)
    Δν = step(νs)
    fk = ComplexF64[freq_kernel(filter, ν) for ν in νs]
    ĝ = Matrix{ComplexF64}(undef, Nν, Nν)
    pref = ComplexF64(1) / 2im
    for q in 1:Nν, p in 1:Nν
        th = tanh(β * (νs[q] - νs[p]) / 4)
        ĝ[p, q] = pref * th * fk[p] * conj(fk[q])
    end
    Tmax = TMAX_FAC * β
    Nt = 2 * round(Int, Tmax / TAU_LCU) + 1
    ts = range(-Tmax, Tmax; length = Nt)
    τ = step(ts)
    Nsrc = Nν * Nν
    src_x = Float64[-νs[p] for q in 1:Nν, p in 1:Nν][:]
    src_y = Float64[ νs[q] for q in 1:Nν, p in 1:Nν][:]
    src_c = ComplexF64[ĝ[p, q] for q in 1:Nν, p in 1:Nν][:]
    Ntgt = Nt * Nt
    tgt_x = Float64[ts[m] for nn in 1:Nt, m in 1:Nt][:]
    tgt_y = Float64[ts[nn] for nn in 1:Nt, m in 1:Nt][:]
    plan = QuantumFurnace.FINUFFT.finufft_makeplan(3, 2, +1, 1, 1e-10; dtype=Float64, nthreads=1)
    QuantumFurnace.FINUFFT.finufft_setpts!(plan, src_x, src_y, Float64[], tgt_x, tgt_y, Float64[])
    out = Vector{ComplexF64}(undef, Ntgt)
    QuantumFurnace.FINUFFT.finufft_exec!(plan, src_c, out)
    QuantumFurnace.FINUFFT.finufft_destroy!(plan)
    norm_factor = (Δν / (2π))^2
    return sum(abs(v * norm_factor) for v in out) * τ^2
end

compute_Z_g(filter::AbstractFilter, β) = _Z_g_single(filter, β)
compute_Z_g(filter::DLLMultiChannelFilter, β) =
    sum(_Z_g_single(c, β) for c in filter.channels)

# ── Full-Lindbladian norms via dense superop materialisation ──────────────────
# Materialises L_S as a dim²×dim² matrix from `apply_lindbladian!`, then
# extracts ‖L‖_op (HS-induced 2-2 superop norm) and the scale-invariant
# intrinsic mixing ratio ρ = λ / Λ_max (KMS Dirichlet form).
function _build_L_super(config, ham, jumps)
    dim = size(ham.data, 1)
    ws  = Workspace(config, ham, jumps)
    L_apply! = function(out, X)
        apply_lindbladian!(ws, X, config, ham)
        copyto!(out, ws.scratch.rho_out)
        return out
    end
    return build_dense_superoperator(L_apply!, dim)
end

function dll_norms(filter, ham, jumps, beta)
    cfg = Config(;
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
    L_super = _build_L_super(cfg, ham, jumps)
    op = hs_operator_norm(L_super)
    ρ  = intrinsic_mixing_ratio(L_super, Matrix{ComplexF64}(ham.gibbs))
    return (op = op, rho = ρ)
end

function ckg_norms(ham, jumps, beta)
    cfg = Config(;
        sim                       = Lindbladian(),
        domain                    = BohrDomain(),
        construction              = KMS(),
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
        filter                    = nothing,
    )
    L_super = _build_L_super(cfg, ham, jumps)
    op = hs_operator_norm(L_super)
    ρ  = intrinsic_mixing_ratio(L_super, Matrix{ComplexF64}(ham.gibbs))
    return (op = op, rho = ρ)
end

# ── Per-cell driver ──────────────────────────────────────────────────────────
function _run_cell(ham, jumps, beta::Float64, k::Int, c_max::Float64)
    filter = _build_filter(beta, k, c_max)

    # τ_mix from the integrator + bi-exp pipeline.
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

    gap_est = try
        kr = krylov_spectral_gap(config, ham, jumps;
                                  krylovdim = 30, howmany = 2, tol = 1e-8)
        kr.spectral_gap > 0 ? kr.spectral_gap : 1.0 / beta
    catch err
        @warn "krylov_spectral_gap failed; falling back to 1/β" beta k c_max err
        1.0 / beta
    end
    t_max = T_MAX_FACTOR / max(gap_est, 1e-12)
    t_grid = collect(range(0.0, t_max; length = T_GRID_LEN))

    rho_0 = _make_init_state(:maximally_mixed,
                             size(ham.data, 1),
                             Matrix{ComplexF64}(ham.gibbs),
                             SEED)

    res = integrate_to_gibbs(config, ham, jumps, rho_0, t_grid;
                              mode = :L, krylovdim = KRYLOV_DIM, tol = TOL)
    est = estimate_mixing_time(res; model = :biexp,
                                target_epsilon = TARGET_EPS,
                                extrapolate = true)

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

    # Block-encoding factors. The CORRECT LW23 multi-channel aggregation:
    #   α_G = (Σ_ℓ Z_g^(ℓ)) · Z_A² · √|𝒜|     (sum: linear in k for identical channels)
    #   A_L = max_ℓ (Z_f^(ℓ) · Z_A)            (max, NOT sum: per-Lindblad BE factor)
    #   ‖L‖_be = α_G + ½ · A_L² · (|𝒜|·k)     (J = |𝒜|·k total Lindblad operators)
    # For uniform rescale (k identical channels at c=0) this gives ‖L‖_be ∝ k,
    # which combined with τ_mix ∝ 1/k makes the product scale-invariant. ✓
    #
    # The earlier "naive" formula (α_L = Σ_ℓ Z_f^(ℓ) · Z_A, no extra k factor)
    # double-counts: it summed per-channel Z_f's into a single block-encoding
    # factor as if the k channels were one giant LCU prep, instead of k separate
    # diagonal blocks each with its own BE factor.
    Z_g_per_channel, Z_f_per_channel = if k == 1
        ([compute_Z_g(filter, beta)], [compute_Z_f(filter, beta)])
    else
        ch = filter.channels
        ([compute_Z_g(c, beta) for c in ch], [compute_Z_f(c, beta) for c in ch])
    end
    Z_A = maximum(opnorm(j.in_eigenbasis) for j in jumps)
    A_card = length(jumps)
    Z_g_total = sum(Z_g_per_channel)
    Z_f_max   = maximum(Z_f_per_channel)
    α_G = Z_g_total * Z_A^2 * sqrt(A_card)
    α_L_per = Z_f_max * Z_A
    L_be = α_G + 0.5 * α_L_per^2 * A_card * k     # J = |𝒜|·k

    # Full-Lindbladian norms (HS-induced 2-2, and ρ_intrinsic = λ/Λ_max).
    norms = dll_norms(filter, ham, jumps, beta)

    return (
        beta = beta, k = k, c_max = c_max,
        gap_est = gap_est, t_max = t_max,
        fitted_gap = est.fitted_gap,
        mixing_time = mt, mixing_time_source = mt_source,
        Z_g_total = Z_g_total, Z_f_max = Z_f_max,
        α_L_per = α_L_per, α_G = α_G, L_be = L_be,
        product = L_be * mt,
        L_op = norms.op, ρ_intr = norms.rho,
    )
end

# ── Driver ───────────────────────────────────────────────────────────────────
println("="^110)
println("Multi-rank DLL τ_mix + ‖L‖_be sweep, qf-9ld follow-up to qf-7go.6")
println("="^110)
@printf("n=%d, β ∈ %s, k ∈ %s, c_max ∈ %s\n",
        N_QUBITS, BETA_VALUES, K_VALUES, C_MAX_VALUES)
println("="^110)
flush(stdout)

ham_path = joinpath(dirname(@__DIR__), "hamiltonians",
                     "heis_disordered_periodic_n$(N_QUBITS).bson")
isfile(ham_path) || error("Hamiltonian fixture missing: $ham_path")

results = NamedTuple[]
ckg_baseline = Dict{Float64, NamedTuple}()
t0_run = time()
for beta in BETA_VALUES
    ham_β = _load_hamiltonian_bson(ham_path, beta)
    jumps_β = _build_jump_set(ham_β, N_QUBITS)

    # CKG baseline at this β
    ckg = ckg_norms(ham_β, jumps_β, beta)
    ckg_baseline[beta] = ckg
    @printf("\n--- β = %.1f --- CKG baseline: ‖L‖_op = %.4e   ρ_intrinsic = %.4f\n",
            beta, ckg.op, ckg.rho)

    for c_max in C_MAX_VALUES
        @printf("  c_max = %.3f\n", c_max)
        for k in K_VALUES
            tcell = time()
            res = _run_cell(ham_β, jumps_β, beta, k, c_max)
            push!(results, res)
            @printf("    k=%-2d  τ_mix=%-9.3g  ‖L‖_op=%-10.3e  ρ_intr=%-7.4f  ‖L‖_be=%-10.3e  product=%-10.3e   wall=%.1fs\n",
                    k, res.mixing_time, res.L_op, res.ρ_intr,
                    res.L_be, res.product, time() - tcell)
            flush(stdout)
        end
    end
end
@printf("\nTotal wall time: %.1fs\n", time() - t0_run)

# Save BSON.
bson_payload = Dict{Symbol, Any}(
    :results        => results,
    :ckg_baseline   => ckg_baseline,
    :n              => N_QUBITS,
    :beta_values    => BETA_VALUES,
    :k_values       => K_VALUES,
    :c_max_values   => C_MAX_VALUES,
    :target_epsilon => TARGET_EPS,
    :tau_lcu        => TAU_LCU,
    :tmax_fac       => TMAX_FAC,
    :seed           => SEED,
)
BSON.bson(BSON_OUT, bson_payload)
@printf("\nWrote BSON: %s\n", BSON_OUT)

# ── Summary tables ───────────────────────────────────────────────────────────
println("\n" * "="^110)
println("SUMMARY — τ_mix vs (k, c_max) at each β")
println("="^110)
for beta in BETA_VALUES
    @printf("\n[β = %.1f]\n", beta)
    @printf("  c_max      ")
    for k in K_VALUES
        @printf("k=%-9d ", k)
    end
    println()
    @printf("  -----------")
    for _ in K_VALUES
        @printf("-----------")
    end
    println()
    for c_max in C_MAX_VALUES
        @printf("  %-11.3f", c_max)
        for k in K_VALUES
            row = findfirst(r -> r.beta == beta && r.k == k && r.c_max == c_max,
                            results)
            mt = row === nothing ? NaN : results[row].mixing_time
            @printf("%-11.3g", mt)
        end
        println()
    end
end

println("\n" * "="^110)
println("‖L‖_be vs (k, c_max) at each β")
println("="^110)
for beta in BETA_VALUES
    @printf("\n[β = %.1f]\n", beta)
    @printf("  c_max      ")
    for k in K_VALUES
        @printf("k=%-9d ", k)
    end
    println()
    @printf("  -----------")
    for _ in K_VALUES
        @printf("-----------")
    end
    println()
    for c_max in C_MAX_VALUES
        @printf("  %-11.3f", c_max)
        for k in K_VALUES
            row = findfirst(r -> r.beta == beta && r.k == k && r.c_max == c_max,
                            results)
            v = row === nothing ? NaN : results[row].L_be
            @printf("%-11.3e", v)
        end
        println()
    end
end

println("\n" * "="^110)
println("PRODUCT ‖L‖_be · τ_mix vs (k, c_max) at each β  ← the algorithmic figure of merit")
println("="^110)
for beta in BETA_VALUES
    @printf("\n[β = %.1f]\n", beta)
    @printf("  c_max      ")
    for k in K_VALUES
        @printf("k=%-9d ", k)
    end
    println()
    @printf("  -----------")
    for _ in K_VALUES
        @printf("-----------")
    end
    println()
    for c_max in C_MAX_VALUES
        @printf("  %-11.3f", c_max)
        for k in K_VALUES
            row = findfirst(r -> r.beta == beta && r.k == k && r.c_max == c_max,
                            results)
            v = row === nothing ? NaN : results[row].product
            @printf("%-11.3e", v)
        end
        println()
    end
end

# Ratio table: product vs k=1 baseline at the same (β, c_max).
println("\n" * "="^110)
println("RATIO: product / product(k=1, same c_max)  ← < 1.0 means multi-rank wins")
println("="^110)
for beta in BETA_VALUES
    @printf("\n[β = %.1f]\n", beta)
    @printf("  c_max      ")
    for k in K_VALUES
        @printf("k=%-9d ", k)
    end
    println()
    @printf("  -----------")
    for _ in K_VALUES
        @printf("-----------")
    end
    println()
    for c_max in C_MAX_VALUES
        @printf("  %-11.3f", c_max)
        baseline_row = findfirst(r -> r.beta == beta && r.k == 1 && r.c_max == c_max, results)
        baseline = baseline_row === nothing ? NaN : results[baseline_row].product
        for k in K_VALUES
            row = findfirst(r -> r.beta == beta && r.k == k && r.c_max == c_max,
                            results)
            v = row === nothing ? NaN : results[row].product / baseline
            @printf("%-11.3f", v)
        end
        println()
    end
end

println()

# ── New tables: ‖L‖_op vs CKG, and ρ_intrinsic vs CKG ─────────────────────────
println("\n" * "="^110)
println("‖L‖_op (HS-induced 2-2 superop norm) vs (k, c_max) at each β")
println("Reported as DLL/CKG ratio — > 1 means DLL Lindbladian is bigger than CKG's")
println("="^110)
for beta in BETA_VALUES
    @printf("\n[β = %.1f]   CKG baseline ‖L‖_op = %.4e\n", beta, ckg_baseline[beta].op)
    @printf("  c_max      ")
    for k in K_VALUES
        @printf("k=%-9d ", k)
    end
    println()
    @printf("  -----------")
    for _ in K_VALUES
        @printf("-----------")
    end
    println()
    for c_max in C_MAX_VALUES
        @printf("  %-11.3f", c_max)
        for k in K_VALUES
            row = findfirst(r -> r.beta == beta && r.k == k && r.c_max == c_max, results)
            v = row === nothing ? NaN : results[row].L_op / ckg_baseline[beta].op
            @printf("%-11.3f", v)
        end
        println()
    end
end

println("\n" * "="^110)
println("ρ_intrinsic (= λ/Λ_max, scale-invariant) — closer to 1.0 means structurally faster")
println("Compared to CKG baseline (printed in header).")
println("="^110)
for beta in BETA_VALUES
    @printf("\n[β = %.1f]   CKG baseline ρ_intrinsic = %.4f\n", beta, ckg_baseline[beta].rho)
    @printf("  c_max      ")
    for k in K_VALUES
        @printf("k=%-9d ", k)
    end
    println()
    @printf("  -----------")
    for _ in K_VALUES
        @printf("-----------")
    end
    println()
    for c_max in C_MAX_VALUES
        @printf("  %-11.3f", c_max)
        for k in K_VALUES
            row = findfirst(r -> r.beta == beta && r.k == k && r.c_max == c_max, results)
            v = row === nothing ? NaN : results[row].ρ_intr
            @printf("%-11.4f", v)
        end
        println()
    end
end

println()
