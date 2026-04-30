#!/usr/bin/env julia
#
# End-to-end resource estimation for Chen's quantum Gibbs sampling algorithm:
#   1. KMS path: truncated sim -> biexp fit -> mixing time -> Ham sim budget
#   2. GNS path: find sigma via Lindbladian eigendecomp (GNS approx gap ~ 1e-5),
#      then spectral gap -> analytical mixing time -> Ham sim budget
#   3. Side-by-side comparison table
#
# For GNS at gap ~ 1e-5, sigma is so small that simulation-based fitting is
# impractical (millions of steps). Instead, we extract the Lindbladian spectral
# gap from eigendecomposition (instant for 64x64) and compute T_mix analytically.
#
# Usage:  OPENBLAS_NUM_THREADS=4 julia -t4 --project scripts/simulation_time_budget.jl

using QuantumFurnace
using LinearAlgebra
using BSON
using Random
using Printf

BLAS.set_num_threads(4)

# ── 1. Load Hamiltonian ──────────────────────────────────────────────────────
num_qubits = 3
dim        = 2^num_qubits
beta       = 10.0

ham_path = joinpath(@__DIR__, "..", "hamiltonians",
                    "heis_disordered_periodic_n$(num_qubits).bson")

raw = open(ham_path) do io; BSON.parse(io) end
fields = raw[:hamiltonian][:data]

cache = IdDict()
data_matrix      = BSON.raise_recursive(fields[1],  cache, QuantumFurnace)::Matrix{ComplexF64}
base_terms       = Vector{Vector{Matrix{ComplexF64}}}(BSON.raise_recursive(fields[4], cache, QuantumFurnace))
base_coeffs      = BSON.raise_recursive(fields[5],  cache, QuantumFurnace)::Vector{Float64}
disordering_term = let dt = BSON.raise_recursive(fields[6], cache, QuantumFurnace)
    dt === nothing ? nothing : Vector{Matrix{ComplexF64}}(dt)
end
disordering_coeffs = let dc = BSON.raise_recursive(fields[7], cache, QuantumFurnace)
    dc === nothing ? nothing : Vector{Float64}(dc)
end
eigvals_vec      = BSON.raise_recursive(fields[8],  cache, QuantumFurnace)::Vector{Float64}
eigvecs_mat      = BSON.raise_recursive(fields[9],  cache, QuantumFurnace)::Matrix{ComplexF64}

raw_nt = (
    matrix             = data_matrix,
    terms              = base_terms,
    base_coeffs        = base_coeffs,
    disordering_term   = disordering_term,
    disordering_coeffs = disordering_coeffs,
    eigvals            = eigvals_vec,
    eigvecs            = eigvecs_mat,
    nu_min             = Float64(fields[10]),
    shift              = Float64(fields[11]),
    rescaling_factor   = Float64(fields[12]),
    periodic           = Bool(fields[13]),
)

hamiltonian = HamHam(raw_nt, beta)
gibbs       = hamiltonian.gibbs

# ── 2. Jump operators ────────────────────────────────────────────────────────
jump_paulis = [[X], [Y], [Z]]
n_jumps     = length(jump_paulis) * num_qubits
norm_factor = sqrt(n_jumps)

jumps = JumpOp[]
for pauli in jump_paulis
    for site in 1:num_qubits
        op       = Matrix(pad_term(pauli, num_qubits, site)) ./ norm_factor
        op_eigen = hamiltonian.eigvecs' * op * hamiltonian.eigvecs
        push!(jumps, JumpOp(op, op_eigen, op == transpose(op), op == op'))
    end
end

# ── 3. Shared parameters ────────────────────────────────────────────────────
delta      = 0.001
target     = 1e-4
r          = 12
w0         = 0.05
t0         = 2π / (2^r * w0)

println("=" ^ 72)
println("  System: $(num_qubits)-qubit disordered Heisenberg (periodic)")
@printf("  dim = %d, beta = %.1f, delta = %.4f, target eps = %.0e\n",
        dim, beta, delta, target)
println("  QPE grid: r = $r, N = $(2^r), w0 = $w0")
println("=" ^ 72)

# Helper: build a Thermalize Config (with optional grid override)
function make_therm_config(; construction, sigma, mixing_time, r_bits=r, w0_val=w0)
    t0_val = 2π / (2^r_bits * w0_val)
    Config(;
        sim                       = Thermalize(),
        domain                    = TimeDomain(),
        construction              = construction,
        num_qubits                = num_qubits,
        with_linear_combination   = true,
        beta                      = beta,
        sigma                     = sigma,
        a                         = beta / 30.0,
        s                         = 0.4,
        num_energy_bits           = r_bits,
        w0                        = w0_val,
        t0                        = t0_val,
        num_trotter_steps_per_t0  = 10,
        mixing_time               = mixing_time,
        delta                     = delta,
    )
end

# Helper: build a Lindbladian Config (with optional grid override)
function make_lind_config(; construction, sigma, r_bits=r, w0_val=w0)
    t0_val = 2π / (2^r_bits * w0_val)
    Config(;
        sim                       = Lindbladian(),
        domain                    = TimeDomain(),
        construction              = construction,
        num_qubits                = num_qubits,
        with_linear_combination   = true,
        beta                      = beta,
        sigma                     = sigma,
        a                         = beta / 30.0,
        s                         = 0.4,
        num_energy_bits           = r_bits,
        w0                        = w0_val,
        t0                        = t0_val,
        num_trotter_steps_per_t0  = 10,
    )
end

# Helper: extract Lindbladian spectral gap and fixed-point distance from eigendecomp
function lindbladian_analysis(; construction, sigma, r_bits=r, w0_val=w0)
    lind_config = make_lind_config(; construction, sigma, r_bits, w0_val)
    liouv = redirect_stdout(devnull) do
        construct_lindbladian(jumps, lind_config, hamiltonian)
    end
    eig = eigen(liouv)

    # Fixed point: eigenvalue closest to 0
    real_parts = real.(eig.values)
    ss_idx = argmin(abs.(real_parts))
    ss_vec = eig.vectors[:, ss_idx]
    ss_dm = reshape(ss_vec, dim, dim)
    ss_dm = (ss_dm + ss_dm') / 2
    ss_dm ./= tr(ss_dm)
    fp_dist = trace_distance_h(Hermitian(ss_dm), gibbs)

    # Spectral gap: second-smallest |Re(lambda)| (the gap to the steady state)
    sorted_re = sort(abs.(real_parts))
    spectral_gap = sorted_re[2]  # first is ~0 (steady state), second is the gap

    return (; fp_dist, spectral_gap, fixed_point=ss_dm)
end


# ══════════════════════════════════════════════════════════════════════════════
#  PART A: KMS — simulation + bi-exponential fit
# ══════════════════════════════════════════════════════════════════════════════
sigma_kms = 1.0 / beta

println()
println("-" ^ 72)
println("  KMS: Mixing time estimation  (sigma = $sigma_kms)")
println("-" ^ 72)

# Also get KMS Lindbladian spectral gap for reference (at r=12)
kms_analysis = lindbladian_analysis(; construction=KMS(), sigma=sigma_kms)
@printf("  Lindbladian spectral gap  = %.6f\n", kms_analysis.spectral_gap)
@printf("  Fixed-point dist to Gibbs = %.4e  (KMS exact DB -> ~0 in BohrDomain)\n",
        kms_analysis.fp_dist)

extrap_time   = 45.0
extrap_steps  = round(Int, extrap_time / delta)
@printf("  Running truncated simulation (T = %.1f, %d steps)...\n", extrap_time, extrap_steps)
flush(stdout)

config_kms_sim = make_therm_config(; construction=KMS(), sigma=sigma_kms, mixing_time=extrap_time)
result_kms = redirect_stdout(devnull) do
    run_thermalize(jumps, config_kms_sim, hamiltonian;
        rng        = MersenneTwister(42),
        save_every = 100,
    )
end

wall_kms = result_kms.metadata[:wall_time_seconds]
@printf("  Done in %.1fs  (%d data points)\n", wall_kms, length(result_kms.trace_distances))
@printf("  Final trace distance: %.4e\n\n", result_kms.trace_distances[end])

est_kms = estimate_mixing_time(result_kms;
    model          = :biexp,
    skip_initial   = 0.2,
    target_epsilon = target,
    extrapolate    = true,
)

bifit_kms = est_kms.biexp_fit_result
println("  Bi-exponential fit:")
@printf("    A_fast = %.4f   g_fast = %.4f\n", bifit_kms.amplitude_fast, bifit_kms.gap_fast)
@printf("    A_slow = %.4f   g_slow = %.4f  (spectral gap)\n", bifit_kms.amplitude, bifit_kms.gap)
@printf("    C      = %.4e  (asymptotic floor)\n", bifit_kms.offset)
@printf("    R^2    = %.6f,  converged: %s\n", est_kms.r_squared, est_kms.converged ? "yes" : "no")

if est_kms.mixing_time_extrapolated === nothing || isnan(est_kms.mixing_time_extrapolated)
    println("\n  ERROR: KMS extrapolation failed")
    exit(1)
end

T_mix_kms = est_kms.mixing_time_extrapolated
@printf("\n  --> T_mix (biexp) = %.2f  (%d steps)\n", T_mix_kms, ceil(Int, T_mix_kms / delta))


# ══════════════════════════════════════════════════════════════════════════════
#  PART B: GNS — scan (sigma, r) to find parameters that reach target
# ══════════════════════════════════════════════════════════════════════════════
println()
println("=" ^ 72)
println("  GNS: Scanning (sigma, r) for fixed-point gap < target")
println("=" ^ 72)

# Step 1: sigma scan at default r to see the non-monotonic floor
println()
println("  Step 1: sigma scan at r = $r (same grid as KMS)")
sigma_candidates = [0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001]
@printf("  %-12s  %-14s  %-14s\n", "sigma", "GNS gap", "spectral gap")
println("  " * "-" ^ 45)

for s in sigma_candidates
    info = lindbladian_analysis(; construction=GNS(), sigma=s)
    @printf("  %.4e    %.4e    %.6f\n", s, info.fp_dist, info.spectral_gap)
end

println()
println("  Floor is non-monotonic: large sigma -> large GNS approx error,")
println("  small sigma -> filter narrower than grid resolution.")
println("  Need to increase r (QPE bits) to resolve smaller sigma.")

# Step 2: for a range of sigma, adapt w0 to resolve the filter, find min r
# Rule: w0 = min(0.05, sigma/5) so the grid spacing resolves the Gaussian,
#        then increase r until 2^r * w0 covers the spectral range with margin.
println()
println("  Step 2: finding (sigma, w0, r) that achieve floor < $target")
println("  w0 adapted per sigma: w0 = min(0.05, sigma/5)")
@printf("  %-10s  %-10s  %-4s  %-8s  %-14s  %-14s\n",
        "sigma", "w0", "r", "N", "GNS gap", "spectral gap")
println("  " * "-" ^ 72)

gns_configs = []

for s in [0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001]
    w0_s = min(0.05, s / 5.0)
    found = false
    for r_try in 8:20
        N_try = 2^r_try
        # Ensure energy range covers spectrum: N*w0 > ~10 (spectrum is within a few units)
        if N_try * w0_s < 10.0
            continue
        end
        info = lindbladian_analysis(; construction=GNS(), sigma=s, r_bits=r_try, w0_val=w0_s)
        if info.fp_dist < target
            @printf("  %.4e  %.4e  %-4d  %-8d  %.4e    %.6f\n",
                    s, w0_s, r_try, N_try, info.fp_dist, info.spectral_gap)
            push!(gns_configs, (; sigma=s, w0_val=w0_s, r_bits=r_try,
                                  fp_dist=info.fp_dist, spectral_gap=info.spectral_gap, N=N_try))
            found = true
            break
        end
    end
    if !found
        @printf("  %.4e  %.4e  >20   (cannot reach %.0e)\n", s, w0_s, target)
    end
end

if isempty(gns_configs)
    println("\n  ERROR: No (sigma, w0, r) combination achieves GNS floor < $target")
    exit(1)
end

# Step 3: for each valid GNS config, compute total Ham sim time budget
println()
println("  Step 3: total Ham sim time for each valid GNS config")
@printf("  %-10s  %-10s  %-4s  %-10s  %-12s  %-12s  %-14s\n",
        "sigma", "w0", "r", "T_mix", "per_step", "n_steps", "TOTAL")
println("  " * "-" ^ 80)

rho0 = Matrix{ComplexF64}(I, dim, dim) / dim
d0 = trace_distance_h(Hermitian(rho0), gibbs)

best_gns = nothing
best_total = Inf

for gc in gns_configs
    eff_eps = target - gc.fp_dist
    t_mix = log(d0 / eff_eps) / gc.spectral_gap

    cfg = make_therm_config(; construction=GNS(), sigma=gc.sigma,
                              mixing_time=t_mix, r_bits=gc.r_bits, w0_val=gc.w0_val)
    budget = compute_simulation_time(cfg, hamiltonian, t_mix)

    @printf("  %.4e  %.4e  %-4d  %10.2f  %12.2f  %12d  %14.2f\n",
            gc.sigma, gc.w0_val, gc.r_bits, t_mix, budget.per_step_time, budget.n_steps, budget.total_time)

    if budget.total_time < best_total
        best_total = budget.total_time
        best_gns = (; gc..., t_mix, budget)
    end
end

sigma_gns    = best_gns.sigma
w0_gns       = best_gns.w0_val
r_gns        = best_gns.r_bits
gns_fp_dist  = best_gns.fp_dist
gns_spec_gap = best_gns.spectral_gap
T_mix_gns    = best_gns.t_mix
budget_gns   = best_gns.budget

println()
@printf("  Cheapest GNS: sigma = %.4e, w0 = %.4e, r = %d (N = %d)\n",
        sigma_gns, w0_gns, r_gns, 2^r_gns)
@printf("    Floor = %.4e, spectral gap = %.6f, T_mix = %.2f\n",
        gns_fp_dist, gns_spec_gap, T_mix_gns)

# KMS analytical mixing time at target (for comparison with GNS analytical)
effective_eps_kms = target - kms_analysis.fp_dist
T_mix_kms_analytical = log(d0 / effective_eps_kms) / kms_analysis.spectral_gap


# ══════════════════════════════════════════════════════════════════════════════
#  PART C: Hamiltonian simulation time budgets — side by side
# ══════════════════════════════════════════════════════════════════════════════
println()
println("=" ^ 72)
@printf("  Hamiltonian simulation time budgets  (target eps = %.0e)\n", target)
println("=" ^ 72)

# KMS budget (biexp-fitted T_mix)
config_kms = make_therm_config(; construction=KMS(), sigma=sigma_kms, mixing_time=T_mix_kms)
budget_kms = compute_simulation_time(config_kms, hamiltonian, T_mix_kms)

println()
@printf("  %-28s  %14s  %14s\n", "", "KMS", "GNS")
println("  " * "-" ^ 60)
@printf("  %-28s  %14.4f  %14.4e\n", "sigma", sigma_kms, sigma_gns)
@printf("  %-28s  %14.4f  %14.4e\n", "w0", w0, w0_gns)
@printf("  %-28s  %14d  %14d\n", "r (QPE bits)", r, r_gns)
@printf("  %-28s  %14d  %14d\n", "N (grid points)", 2^r, 2^r_gns)
@printf("  %-28s  %14.4e  %14.4e\n", "fixed-point dist to Gibbs", kms_analysis.fp_dist, gns_fp_dist)
@printf("  %-28s  %14.6f  %14.6f\n", "Lindbladian spectral gap", kms_analysis.spectral_gap, gns_spec_gap)
@printf("  %-28s  %14s  %14s\n", "T_mix method", "biexp fit", "analytical")
@printf("  %-28s  %14.2f  %14.2f\n", "T_mix", T_mix_kms, T_mix_gns)
@printf("  %-28s  %14d  %14d\n", "n_steps", budget_kms.n_steps, budget_gns.n_steps)
println("  " * "-" ^ 60)
@printf("  %-28s  %14.2f  %14.2f\n", "OFT time (per step)", budget_kms.oft_time, budget_gns.oft_time)
@printf("  %-28s  %14.2f  %14.2f\n", "B time (per step)", budget_kms.b_time, budget_gns.b_time)
@printf("  %-28s  %14.2f  %14.2f\n", "Per step (2*OFT + B)", budget_kms.per_step_time, budget_gns.per_step_time)
println("  " * "-" ^ 60)
@printf("  %-28s  %14.2f  %14.2f\n", "TOTAL Ham sim time", budget_kms.total_time, budget_gns.total_time)
println("  " * "-" ^ 60)

ratio = budget_gns.total_time / budget_kms.total_time
if ratio > 1
    @printf("\n  GNS is %.1fx MORE expensive than KMS\n", ratio)
else
    @printf("\n  GNS is %.1fx CHEAPER than KMS\n", 1.0 / ratio)
end

println()
println("  Key observations:")
println("  - KMS uses exact detailed balance (B term) at sigma = 1/beta = $(sigma_kms)")
println("  - GNS avoids the B term but needs much smaller sigma to control the")
println("    approximation gap, and more QPE bits (r = $r_gns vs $r) to resolve")
println("    the narrower filter.")
@printf("  - The larger grid (N = %d vs %d) inflates the OFT cost per step.\n",
        2^r_gns, 2^r)
if ratio > 1
    println("  - Net result: the QPE overhead outweighs saving the B term.")
end
println()
@printf("  Wall time: %.1fs\n", wall_kms)
println()
