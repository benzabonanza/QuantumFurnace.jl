#!/usr/bin/env julia
#
# Smooth-Metropolis qubit-cost analysis (qf-3il, corrected).
#
# Headline question
# -----------------
# What is the smallest energy-register spacing `w0` and time-register spacing
# `t0` such that the simulator's discrete Lindbladian matches the BohrDomain
# (no-quadrature-error) reference within target ε? The minimum number of
# qubits then follows from the simulator's own truncation logic
# (`_truncate_energy_labels`, `_truncate_time_labels_for_oft`):
#
#     bits ≥ ⌈log₂(N_trunc)⌉,  N_trunc = length(truncated grid).
#
# This is the right cost question — increasing `num_energy_bits` at fixed `w0`
# does *nothing* because the simulator re-truncates. The bits-per-spacing
# tradeoff applies to *both* registers, with `w0 · t0 = 2π / N_full` tying
# them by Fourier on the un-truncated grid.
#
# Strategy
# --------
# 1. **Reference.** BohrDomain CKG for KMS — has no ω-quadrature error
#    (operates on Bohr frequencies directly via `pick_transition` + KMS shift).
#    Materialise via `build_dense_superoperator`.
# 2. **EnergyDomain sweep.** For each (β, s, w0), build the EnergyDomain L
#    via `apply_lindbladian!` + `build_dense_superoperator` and compute
#    ‖L_E(w0) − L_Bohr‖_op / ‖L_Bohr‖_op. Find the largest w0 such that this
#    relative error is ≤ ε. Record N_eng = length(truncated_energy_labels)
#    at that w0; bits_E = ⌈log₂(N_eng)⌉.
# 3. **TimeDomain sweep.** For each (β, s), fix `w0 = w0_max(ε)` from step 2,
#    then vary `t0`. Build the TimeDomain L and compute the same op-error
#    against L_Bohr. Find the largest t0 such that error ≤ ε. Record
#    N_time = length(oft_time_labels); bits_T = ⌈log₂(N_time)⌉.
# 4. **Total bits** = max(bits_E, bits_T) — both registers share storage in
#    the QPE-style algorithm.
#
# Sweep grid
# ----------
#   β ∈ {5, 10, 20}
#   s ∈ {0.0, 0.25}      (kinky vs the locked thesis default)
#   n = 3                  (small system; the cost question is about γ-shape,
#                           independent of Hamiltonian dim modulo Bohr spread)
#   w0 ∈ {0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625, 0.003125}
#   t0 ∈ {0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125}
#   target ε ∈ {1e-3, 1e-4, 1e-6}
#
# Output
# ------
# - drafts/figures/numerics/smooth_metro_qubit_cost.bson  (per-cell op-errors)
# - drafts/figures/numerics/smooth_metro_qubit_cost_w0.{png,pdf}  (energy)
# - drafts/figures/numerics/smooth_metro_qubit_cost_t0.{png,pdf}  (time)
# - drafts/figures/numerics/smooth_metro_qubit_cost_summary.{png,pdf}
# - Summary table to stdout: bits @ ε for each (β, s).

using Printf
using LinearAlgebra
using BSON
using QuantumFurnace
using QuantumFurnace: _load_hamiltonian_bson, _build_jump_set, _create_energy_labels,
                       _truncate_energy_labels, _truncate_time_labels_for_oft,
                       GaussianFilter
ENV["GKSwstype"] = "100"
using Plots

# ── Sweep grid ────────────────────────────────────────────────────────────────
const N_QUBITS    = 3
const BETA_VALUES = [5.0, 10.0, 20.0]
const S_VALUES    = [0.0, 0.25]
const W0_VALUES   = [0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625, 0.003125]
const T0_VALUES   = [0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125]
const EPS_TARGETS = [1e-3, 1e-4, 1e-6]
const NUM_ENERGY_BITS = 14   # generous; truncation does the actual cap

# ── Output paths ──────────────────────────────────────────────────────────────
const OUT_DIR  = joinpath(@__DIR__, "..", "drafts", "figures", "numerics")
mkpath(OUT_DIR)
const BSON_OUT = joinpath(OUT_DIR, "smooth_metro_qubit_cost.bson")
const FIG_W0_PNG = joinpath(OUT_DIR, "smooth_metro_qubit_cost_w0.png")
const FIG_W0_PDF = joinpath(OUT_DIR, "smooth_metro_qubit_cost_w0.pdf")
const FIG_T0_PNG = joinpath(OUT_DIR, "smooth_metro_qubit_cost_t0.png")
const FIG_T0_PDF = joinpath(OUT_DIR, "smooth_metro_qubit_cost_t0.pdf")
const FIG_SUM_PNG = joinpath(OUT_DIR, "smooth_metro_qubit_cost_summary.png")
const FIG_SUM_PDF = joinpath(OUT_DIR, "smooth_metro_qubit_cost_summary.pdf")

# ── Thesis colour palette ─────────────────────────────────────────────────────
const COLOR_S = Dict(0.0 => "#222222", 0.25 => "#5F8B8E")  # kinky black, smooth dustyteal
const STYLE_S = Dict(0.0 => :dash, 0.25 => :solid)
const COLOR_BETA = Dict(5.0 => "#5C7794", 10.0 => "#B5654A", 20.0 => "#7A2E39")

# ── Helpers ───────────────────────────────────────────────────────────────────
"""Materialise the dense d²×d² Lindbladian for a given (domain, β, s, w0, t0).
The η = 0.05 default matches the test-suite (`test_helpers.jl`); it only enters
the η-regularised b₊ used by TimeDomain (BohrDomain/EnergyDomain are unaffected
because they don't go through `_compute_b_plus_metro`)."""
function dense_lindbladian(domain, beta::Float64, s::Float64, w0::Float64, t0::Float64,
                            ham, jumps, n::Int)
    cfg = Config(;
        sim                       = Lindbladian(),
        domain                    = domain,
        construction              = KMS(),
        num_qubits                = n,
        with_linear_combination   = true,
        beta                      = beta,
        sigma                     = 1.0 / beta,
        a                         = 0.0,
        s                         = s,
        eta                       = 0.05,
        num_energy_bits           = NUM_ENERGY_BITS,
        w0                        = w0,
        t0                        = t0,
        num_trotter_steps_per_t0  = 10,
    )
    ws = Workspace(cfg, ham, jumps)
    d = size(ham.data, 1)
    L = build_dense_superoperator(d) do out, X
        apply_lindbladian!(ws, X, cfg, ham); copyto!(out, ws.scratch.rho_out); out
    end
    return L, cfg
end

"""Length of the truncated energy register at the given (β, s, w0)."""
function n_trunc_energy(beta::Float64, s::Float64, w0::Float64, ham, jumps, n::Int)
    cfg = Config(;
        sim=Lindbladian(), domain=EnergyDomain(), construction=KMS(),
        num_qubits=n, with_linear_combination=true,
        beta=beta, sigma=1.0/beta, a=0.0, s=s, eta=0.05,
        num_energy_bits=NUM_ENERGY_BITS, w0=w0, t0=2π/(2^NUM_ENERGY_BITS*w0),
        num_trotter_steps_per_t0=10,
    )
    return length(_truncate_energy_labels(_create_energy_labels(NUM_ENERGY_BITS, w0), cfg))
end

"""Length of the truncated OFT-time register at the given (β, s, w0, t0)."""
function n_trunc_time(beta::Float64, s::Float64, w0::Float64, t0::Float64,
                       ham, jumps, n::Int)
    cfg = Config(;
        sim=Lindbladian(), domain=TimeDomain(), construction=KMS(),
        num_qubits=n, with_linear_combination=true,
        beta=beta, sigma=1.0/beta, a=0.0, s=s, eta=0.05,
        num_energy_bits=NUM_ENERGY_BITS, w0=w0, t0=t0,
        num_trotter_steps_per_t0=10,
    )
    energy_labels = _create_energy_labels(NUM_ENERGY_BITS, w0)
    time_labels   = energy_labels .* (t0 / w0)
    sigma         = 1.0 / beta
    # CKG TimeDomain uses a GaussianFilter(sigma) for OFT-time truncation
    # (matches `_truncate_time_labels_for_oft` default arg).
    return length(_truncate_time_labels_for_oft(time_labels, sigma;
                                                 filter=GaussianFilter(sigma)))
end

# ── Driver ────────────────────────────────────────────────────────────────────
println("="^72)
println("Smooth-Metropolis qubit-cost analysis (qf-3il, corrected)")
println("="^72)
@printf("n           : %d\n", N_QUBITS)
@printf("β values    : %s\n", string(BETA_VALUES))
@printf("s values    : %s\n", string(S_VALUES))
@printf("w0 values   : %s\n", string(W0_VALUES))
@printf("t0 values   : %s\n", string(T0_VALUES))
@printf("ε targets   : %s\n", string(EPS_TARGETS))
@printf("num_energy_bits (max, before truncation) : %d\n", NUM_ENERGY_BITS)
println("="^72)
flush(stdout)

# Per-cell results: keyed by (beta, s) → vectors of NamedTuples for each (w0,) or (w0, t0)
results = Dict{Tuple{Float64, Float64}, Dict{Symbol, Any}}()
ham_cache = Dict{Float64, Any}()
jumps_cache = Dict{Float64, Any}()

ham_path = joinpath(dirname(@__DIR__), "hamiltonians",
                     "heis_xxx_zzdisordered_periodic_n$(N_QUBITS).bson")

for beta in BETA_VALUES
    if !haskey(ham_cache, beta)
        ham_cache[beta]   = _load_hamiltonian_bson(ham_path, beta)
        jumps_cache[beta] = _build_jump_set(ham_cache[beta], N_QUBITS)
    end
    ham = ham_cache[beta]
    jumps = jumps_cache[beta]
    for s in S_VALUES
        # Reference: BohrDomain (no ω-quadrature). Single materialisation per (β, s).
        L_bohr, _ = dense_lindbladian(BohrDomain(), beta, s, 0.05,
                                       2π/(2^NUM_ENERGY_BITS*0.05), ham, jumps, N_QUBITS)
        ref_norm = opnorm(L_bohr)

        # ── EnergyDomain sweep over w0 ──
        eng_rows = NamedTuple[]
        for w0 in W0_VALUES
            t0_paired = 2π / (2^NUM_ENERGY_BITS * w0)
            L_e, _ = dense_lindbladian(EnergyDomain(), beta, s, w0, t0_paired,
                                        ham, jumps, N_QUBITS)
            err_op = opnorm(L_e - L_bohr) / ref_norm
            N_e = n_trunc_energy(beta, s, w0, ham, jumps, N_QUBITS)
            push!(eng_rows, (w0=w0, op_err=err_op, N_trunc=N_e,
                             bits=ceil(Int, log2(max(N_e, 2)))))
        end

        # ── TimeDomain sweep over t0 (with each w0) ──
        time_rows = NamedTuple[]
        for w0 in W0_VALUES, t0 in T0_VALUES
            L_t, _ = dense_lindbladian(TimeDomain(), beta, s, w0, t0,
                                        ham, jumps, N_QUBITS)
            err_op = opnorm(L_t - L_bohr) / ref_norm
            N_t = n_trunc_time(beta, s, w0, t0, ham, jumps, N_QUBITS)
            push!(time_rows, (w0=w0, t0=t0, op_err=err_op, N_trunc=N_t,
                              bits=ceil(Int, log2(max(N_t, 2)))))
        end

        results[(beta, s)] = Dict(
            :ref_norm  => ref_norm,
            :eng_rows  => eng_rows,
            :time_rows => time_rows,
        )

        @printf("\nβ=%-5.1f s=%-4.2f   ‖L_Bohr‖_op = %.4f\n", beta, s, ref_norm)
        @printf("  EnergyDomain (w0 sweep):\n")
        for r in eng_rows
            @printf("    w0=%-8.5f  N_trunc=%-4d  bits=%-2d  rel_err=%.3e\n",
                    r.w0, r.N_trunc, r.bits, r.op_err)
        end
        flush(stdout)
    end
end

# ── Persist ───────────────────────────────────────────────────────────────────
bson_payload = Dict{Symbol, Any}(
    :results       => results,
    :n_qubits      => N_QUBITS,
    :beta_values   => BETA_VALUES,
    :s_values      => S_VALUES,
    :w0_values     => W0_VALUES,
    :t0_values     => T0_VALUES,
    :eps_targets   => EPS_TARGETS,
    :num_energy_bits_max => NUM_ENERGY_BITS,
    :git_sha       => try
        read(`git -C $(dirname(@__DIR__)) rev-parse HEAD`, String) |> strip
    catch
        "unknown"
    end,
)
BSON.bson(BSON_OUT, bson_payload)
@printf("\nWrote BSON: %s\n", BSON_OUT)

# ── Plot 1: rel_err vs w0 (EnergyDomain), one curve per (β, s) ───────────────
plots_w0 = Plots.Plot[]
for beta in BETA_VALUES
    p = plot(;
        xlabel = "w0",
        ylabel = "‖L_E - L_Bohr‖_op / ‖L_Bohr‖_op",
        title  = "EnergyDomain — β = $(beta)",
        xscale = :log10,
        yscale = :log10,
        legend = :bottomright,
        framestyle = :box,
        size = (520, 380),
    )
    for s in S_VALUES
        rec = results[(beta, s)][:eng_rows]
        ws = [r.w0 for r in rec]
        es = [max(r.op_err, 1e-17) for r in rec]
        plot!(p, ws, es;
            color = COLOR_S[s],
            linestyle = STYLE_S[s],
            linewidth = 2,
            marker = :circle, markersize = 4, markerstrokewidth = 0,
            label = "s = $(s)$(s == 0.0 ? " (kinky)" : "")",
        )
    end
    for eps in EPS_TARGETS
        hline!(p, [eps]; color=:black, linestyle=:dot, linewidth=0.7, alpha=0.5, label="")
    end
    push!(plots_w0, p)
end
fig_w0 = plot(plots_w0...; layout=(1, length(BETA_VALUES)),
              size=(520*length(BETA_VALUES), 380))
savefig(fig_w0, FIG_W0_PNG); savefig(fig_w0, FIG_W0_PDF)
@printf("Wrote: %s\n", FIG_W0_PNG)

# ── Plot 2: rel_err vs t0 at the converged w0_min, one curve per (β, s) ─────
# Pick "converged w0" = smallest w0 in our grid (~3e-3), so EnergyDomain error
# is in the 1e-15 range and the residual is the time-domain construction error.
const W0_REF = minimum(W0_VALUES)
plots_t0 = Plots.Plot[]
for beta in BETA_VALUES
    p = plot(;
        xlabel = "t0",
        ylabel = "‖L_T - L_Bohr‖_op / ‖L_Bohr‖_op",
        title  = "TimeDomain — β = $(beta), w0 = $(W0_REF)",
        xscale = :log10,
        yscale = :log10,
        legend = :bottomright,
        framestyle = :box,
        size = (520, 380),
    )
    for s in S_VALUES
        rec = results[(beta, s)][:time_rows]
        rec_ref = filter(r -> r.w0 == W0_REF, rec)
        sort!(rec_ref; by = r -> r.t0)
        ts = [r.t0 for r in rec_ref]
        es = [max(r.op_err, 1e-17) for r in rec_ref]
        plot!(p, ts, es;
            color = COLOR_S[s],
            linestyle = STYLE_S[s],
            linewidth = 2,
            marker = :diamond, markersize = 4, markerstrokewidth = 0,
            label = "s = $(s)$(s == 0.0 ? " (kinky)" : "")",
        )
    end
    for eps in EPS_TARGETS
        hline!(p, [eps]; color=:black, linestyle=:dot, linewidth=0.7, alpha=0.5, label="")
    end
    push!(plots_t0, p)
end
fig_t0 = plot(plots_t0...; layout=(1, length(BETA_VALUES)),
              size=(520*length(BETA_VALUES), 380))
savefig(fig_t0, FIG_T0_PNG); savefig(fig_t0, FIG_T0_PDF)
@printf("Wrote: %s\n", FIG_T0_PNG)

# ── Summary table: bits @ ε for each (β, s) ───────────────────────────────────
function bits_at_eps(rec_eng::Vector{<:NamedTuple}, eps::Float64)
    # Largest w0 such that op_err ≤ eps; bits = ceil(log2(N_trunc)).
    sort!(rec_eng; by = r -> -r.w0)   # descending in w0 (= coarsest first)
    for r in rec_eng
        r.op_err <= eps && return (w0=r.w0, bits=r.bits, N=r.N_trunc, op_err=r.op_err)
    end
    return (w0=NaN, bits=-1, N=0, op_err=NaN)
end

println("\n" * "="^88)
println("SUMMARY — energy-register qubits (bits_E) at each ε")
println("="^88)
@printf("%-6s %-6s %-15s %-15s %-15s\n", "β", "s", "ε=1e-3", "ε=1e-4", "ε=1e-6")
for beta in BETA_VALUES, s in S_VALUES
    res = results[(beta, s)][:eng_rows]
    out = String[]
    for eps in EPS_TARGETS
        b = bits_at_eps(copy(res), eps)
        push!(out, b.bits < 0 ? "(>$(NUM_ENERGY_BITS))" :
                   "$(b.bits) bits, w0=$(round(b.w0, sigdigits=3))")
    end
    @printf("%-6.1f %-6.2f %-15s %-15s %-15s\n", beta, s, out[1], out[2], out[3])
end

println("\n" * "="^88)
println("Bits saved by smooth (s=0.25) vs kinky (s=0)")
println("="^88)
@printf("%-6s %-15s %-15s %-15s\n", "β", "ε=1e-3", "ε=1e-4", "ε=1e-6")
for beta in BETA_VALUES
    res_k = results[(beta, 0.0)][:eng_rows]
    res_s = results[(beta, 0.25)][:eng_rows]
    out = String[]
    for eps in EPS_TARGETS
        bk = bits_at_eps(copy(res_k), eps)
        bs = bits_at_eps(copy(res_s), eps)
        if bk.bits < 0 || bs.bits < 0
            push!(out, "—")
        else
            push!(out, "$(bk.bits)→$(bs.bits)  ($(bk.bits - bs.bits) saved)")
        end
    end
    @printf("%-6.1f %-15s %-15s %-15s\n", beta, out[1], out[2], out[3])
end
println()
println("Done.")
