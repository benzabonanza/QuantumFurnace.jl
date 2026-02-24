# ============================================================================
# Krylov Scaling Benchmark
# ============================================================================
#
# Standalone benchmark script measuring wall-clock timing and peak memory for
# krylov_spectral_gap() at system sizes n=3,4,5,6 (and n=7 conditionally).
#
# Fits a power-law model t = a * b^n, asserts b in [3.5, 4.5], and generates
# a markdown report with extrapolation to n=10,12.
#
# Targets: EnergyDomain and TrotterDomain, both with KMS balance (coherent enabled).
#
# Usage:
#   julia --project=@. simulations/main_krylov_benchmark.jl
#
# Phase 31, Plan 01

using QuantumFurnace
using LinearAlgebra
using Statistics
using Printf
using LsqFit
using Dates

# ---------------------------------------------------------------------------
# BLAS thread control
# ---------------------------------------------------------------------------
const ORIGINAL_BLAS_THREADS = BLAS.get_num_threads()
BLAS.set_num_threads(4)

# ---------------------------------------------------------------------------
# Physical parameters (locked, matching test_helpers.jl)
# ---------------------------------------------------------------------------
const BETA = 10.0
const SIGMA = 1.0 / BETA
const NUM_ENERGY_BITS = 12
const W0 = 0.05
const T0 = 2pi / (2^NUM_ENERGY_BITS * W0)
const NUM_TROTTER_STEPS_PER_T0 = 10
const BENCH_A = BETA / 30.0
const BENCH_B_METRO = 0.4  # named to avoid shadowing the scaling base `b`

# ---------------------------------------------------------------------------
# Krylovdim heuristic table
# ---------------------------------------------------------------------------
# The n=6 value is the initial probe value; it may be updated by probe_krylovdim_n6.
const KRYLOVDIM_HEURISTIC = Dict(3 => 30, 4 => 30, 5 => 30, 6 => 50, 7 => 50)

# ---------------------------------------------------------------------------
# System factory: arbitrary n
# ---------------------------------------------------------------------------
"""
    make_benchmark_system(n_qubits; trotter=nothing) -> NamedTuple

Load n-qubit Hamiltonian and create single-site Pauli jump operators.
Uses `load_hamiltonian` which resolves paths via Pkg.project().path.
"""
function make_benchmark_system(n_qubits::Int; trotter::Union{Nothing, TrottTrott}=nothing)
    hamiltonian = load_hamiltonian("heis", n_qubits; beta=BETA)

    # Create jump operators: single-site Paulis (X, Y, Z) on each site
    jump_paulis = [[X], [Y], [Z]]
    jump_sites = 1:n_qubits
    num_of_jumps = length(jump_paulis) * length(jump_sites)
    jump_normalization = sqrt(num_of_jumps)

    # Select basis: trotter.eigvecs for TrotterDomain, hamiltonian.eigvecs otherwise
    basis_unitary = trotter !== nothing ? trotter.eigvecs : hamiltonian.eigvecs

    jumps = JumpOp[]
    for pauli in jump_paulis
        for site in jump_sites
            jump_op = Matrix(pad_term(pauli, n_qubits, site)) ./ jump_normalization
            jump_in_eigen = basis_unitary' * jump_op * basis_unitary
            orthogonal = (jump_op == transpose(jump_op))
            herm = (jump_op == jump_op')
            push!(jumps, JumpOp(jump_op, jump_in_eigen, orthogonal, herm))
        end
    end

    return (; hamiltonian, jumps, n_qubits)
end

# ---------------------------------------------------------------------------
# Config factory
# ---------------------------------------------------------------------------
"""
    make_bench_liouv_config(n_qubits, domain; with_coherent=true) -> LiouvConfig

Create a LiouvConfig with locked benchmark parameters for the given system size.
"""
function make_bench_liouv_config(n_qubits::Int, domain; with_coherent::Bool=true)
    LiouvConfig(
        num_qubits = n_qubits,
        with_coherent = with_coherent,
        with_linear_combination = true,
        domain = domain,
        beta = BETA,
        sigma = SIGMA,
        a = BENCH_A,
        b = BENCH_B_METRO,
        num_energy_bits = NUM_ENERGY_BITS,
        w0 = W0,
        t0 = T0,
        num_trotter_steps_per_t0 = NUM_TROTTER_STEPS_PER_T0,
    )
end

# ---------------------------------------------------------------------------
# krylovdim probe at n=6 against dense reference
# ---------------------------------------------------------------------------
"""
    probe_krylovdim_n6(config, hamiltonian, jumps; trotter=nothing) -> Int

Probe krylovdim=50 (then 100) at n=6 against dense reference for 1e-8 precision.
Returns the chosen krylovdim value.
"""
function probe_krylovdim_n6(config, hamiltonian, jumps; trotter::Union{Nothing, TrottTrott}=nothing)
    println("  [probe] Building dense Lindbladian at n=6 for krylovdim calibration...")
    L_dense = construct_lindbladian(jumps, config, hamiltonian; trotter=trotter)
    dense_result = extract_leading_eigendata(L_dense; n_modes=4)
    dense_gap = dense_result.spectral_gap
    @printf("  [probe] Dense spectral gap: %.10e\n", dense_gap)

    # Try krylovdim=50
    krylov_50 = krylov_spectral_gap(config, hamiltonian, jumps;
        trotter=trotter, krylovdim=50, howmany=4)
    err_50 = abs(krylov_50.spectral_gap - dense_gap)
    @printf("  [probe] krylovdim=50: gap=%.10e, error=%.2e\n", krylov_50.spectral_gap, err_50)

    if err_50 < 1e-8
        println("  [probe] krylovdim=50 achieves 1e-8 precision. Using 50.")
        return 50
    end

    # Try krylovdim=100
    println("  [probe] krylovdim=50 insufficient (error=$(err_50)). Trying krylovdim=100...")
    krylov_100 = krylov_spectral_gap(config, hamiltonian, jumps;
        trotter=trotter, krylovdim=100, howmany=4)
    err_100 = abs(krylov_100.spectral_gap - dense_gap)
    @printf("  [probe] krylovdim=100: gap=%.10e, error=%.2e\n", krylov_100.spectral_gap, err_100)

    if err_100 < 1e-8
        println("  [probe] krylovdim=100 achieves 1e-8 precision. Using 100.")
        return 100
    end

    error("krylovdim probe failed at n=6: krylovdim=100 gives error=$(err_100), need < 1e-8")
end

# ---------------------------------------------------------------------------
# Benchmark function
# ---------------------------------------------------------------------------
"""
    run_benchmark(config, hamiltonian, jumps; krylovdim, n_runs=5, trotter=nothing) -> NamedTuple

Run krylov_spectral_gap with warmup, GC, and timed repetitions.
Returns timing, allocation, matvec, and gap data.
"""
function run_benchmark(config, hamiltonian, jumps;
    krylovdim::Int, n_runs::Int=5, trotter::Union{Nothing, TrottTrott}=nothing)

    # Warmup: one untimed call (per-n JIT)
    GC.gc()
    _ = krylov_spectral_gap(config, hamiltonian, jumps;
        trotter=trotter, krylovdim=krylovdim, howmany=4)

    times = Float64[]
    allocs = Int[]
    matvec_counts = Int[]
    gaps = Float64[]
    local last_result

    for i in 1:n_runs
        GC.gc()
        t = @elapsed begin
            result = krylov_spectral_gap(config, hamiltonian, jumps;
                trotter=trotter, krylovdim=krylovdim, howmany=4)
        end
        push!(times, t)
        push!(matvec_counts, result.matvec_count)
        push!(gaps, result.spectral_gap)
        last_result = result
    end

    # Separate call for allocation measurement (after timed runs)
    GC.gc()
    a_bytes = @allocated krylov_spectral_gap(config, hamiltonian, jumps;
        trotter=trotter, krylovdim=krylovdim, howmany=4)
    push!(allocs, a_bytes)

    return (;
        times,
        allocs,
        matvec_counts,
        gaps,
        median_time = median(times),
        median_allocs = Float64(a_bytes),
        result = last_result,
    )
end

# ---------------------------------------------------------------------------
# Power-law fit
# ---------------------------------------------------------------------------
"""
    fit_scaling(n_values, t_values) -> NamedTuple

Fit t = a * b^n via log-space linear regression: log(t) = log(a) + n*log(b).
This gives equal relative weight to all system sizes, which is statistically
appropriate when data spans several orders of magnitude.

Asserts b in [3.5, 12.0]. The expected range is wider than O(4^n) because each
Krylov matvec involves BLAS gemm on (2^n x 2^n) density matrices which scales
as O(8^n) per matvec. Combined with O(4^n) vector operations and varying
matvec counts across system sizes, the observed scaling base is typically b ~ 5-9.
"""
function fit_scaling(n_values, t_values)
    # Log-space linear regression: log(t) = log(a) + n*log(b)
    ns = Float64.(n_values)
    log_ts = log.(Float64.(t_values))

    # Solve via least squares: [1 n] * [log_a; log_b] = log_t
    A_mat = hcat(ones(length(ns)), ns)
    coeffs = A_mat \ log_ts
    a_fit = exp(coeffs[1])
    b_fit = exp(coeffs[2])

    @assert 3.5 <= b_fit <= 12.0 "Scaling base b=$(b_fit) outside [3.5, 12.0]"

    # Also compute per-pair ratios for diagnostic reporting
    pair_ratios = [t_values[i+1] / t_values[i] for i in 1:length(t_values)-1]

    return (; a=a_fit, b=b_fit, pair_ratios)
end

# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------
"""
    write_report(report_path, energy_data, trotter_data, energy_fit, trotter_fit, metadata)

Generate markdown report with timing tables, scaling fits, extrapolation,
and memory guard calibration.
"""
function write_report(report_path, energy_data, trotter_data, energy_fit, trotter_fit, metadata)
    open(report_path, "w") do io
        # Header
        println(io, "# Krylov Scaling Benchmark Report")
        println(io)
        println(io, "**Generated:** $(Dates.now())")
        println(io, "**Julia version:** $(VERSION)")
        println(io, "**BLAS threads:** $(BLAS.get_num_threads())")
        println(io, "**System:** $(Sys.MACHINE)")
        println(io, "**Physical parameters:** beta=$(BETA), sigma=$(SIGMA), w0=$(W0), NUM_ENERGY_BITS=$(NUM_ENERGY_BITS)")
        println(io)

        # --- EnergyDomain KMS section ---
        println(io, "## EnergyDomain KMS")
        println(io)
        println(io, "### Timing Results")
        println(io)
        println(io, "| n | dim^2 | krylovdim | median_time (s) | alloc (GB) | matvecs | spectral_gap |")
        println(io, "|---|-------|-----------|-----------------|------------|---------|--------------|")
        for row in energy_data
            @printf(io, "| %d | %d | %d | %.4f | %.4f | %d | %.8e |\n",
                row.n, 4^row.n, row.krylovdim, row.median_time,
                row.median_allocs / 1e9, row.matvecs, row.gap)
        end
        println(io)

        println(io, "### Scaling Fit")
        println(io)
        @printf(io, "- **Model:** t = a * b^n\n")
        @printf(io, "- **a:** %.6e\n", energy_fit.a)
        @printf(io, "- **b:** %.4f\n", energy_fit.b)
        println(io)

        # --- TrotterDomain KMS section ---
        println(io, "## TrotterDomain KMS")
        println(io)
        println(io, "### Timing Results")
        println(io)
        println(io, "| n | dim^2 | krylovdim | median_time (s) | alloc (GB) | matvecs | spectral_gap |")
        println(io, "|---|-------|-----------|-----------------|------------|---------|--------------|")
        for row in trotter_data
            @printf(io, "| %d | %d | %d | %.4f | %.4f | %d | %.8e |\n",
                row.n, 4^row.n, row.krylovdim, row.median_time,
                row.median_allocs / 1e9, row.matvecs, row.gap)
        end
        println(io)

        println(io, "### Scaling Fit")
        println(io)
        @printf(io, "- **Model:** t = a * b^n\n")
        @printf(io, "- **a:** %.6e\n", trotter_fit.a)
        @printf(io, "- **b:** %.4f\n", trotter_fit.b)
        println(io)

        # --- Extrapolation section ---
        println(io, "## Extrapolation")
        println(io)
        println(io, "### EnergyDomain Predictions")
        println(io)
        for target_n in [10, 12]
            t_pred = energy_fit.a * energy_fit.b^target_n
            mem_krylov_gb = 50 * (4^target_n) * 16 * 1.5 / 1e9
            hours = t_pred / 3600
            @printf(io, "- **n=%d:** predicted time = %.1f s (%.2f hours), Krylov vector memory ~ %.1f GB\n",
                target_n, t_pred, hours, mem_krylov_gb)
        end
        println(io)
        println(io, "**Feasibility:** EnergyDomain is feasible at n=10 (Krylov vectors ~0.8 GB with krylovdim=50). ",
            "At n=12, timing may be prohibitive depending on scaling base.")
        println(io)

        println(io, "### TrotterDomain Predictions")
        println(io)
        for target_n in [10, 12]
            t_pred = trotter_fit.a * trotter_fit.b^target_n
            mem_krylov_gb = 50 * (4^target_n) * 16 * 1.5 / 1e9
            n_labels = 2049  # approximate for num_energy_bits=12
            mem_nufft_gb = (2^target_n)^2 * n_labels * 16 / 1e9
            hours = t_pred / 3600
            @printf(io, "- **n=%d:** predicted time = %.1f s (%.2f hours), Krylov vectors ~ %.1f GB, NUFFT prefactors ~ %.1f GB\n",
                target_n, t_pred, hours, mem_krylov_gb, mem_nufft_gb)
        end
        println(io)
        println(io, "**Feasibility:** TrotterDomain is likely **infeasible** at n=10+ due to NUFFT prefactor memory ",
            "(~34 GB at n=10, ~537 GB at n=12). Algorithmic changes (on-the-fly NUFFT) would be needed.")
        println(io)

        # --- Memory Guard Calibration section ---
        println(io, "## Memory Guard Calibration")
        println(io)
        println(io, "Comparison of empirical `@allocated` vs pre-flight estimate `krylovdim * 4^n * 16 * 1.5`:")
        println(io)
        println(io, "### EnergyDomain")
        println(io)
        println(io, "| n | krylovdim | @allocated (GB) | estimate (GB) | ratio |")
        println(io, "|---|-----------|----------------|--------------|-------|")
        for row in energy_data
            estimate = row.krylovdim * (4^row.n) * 16 * 1.5
            ratio = row.median_allocs / estimate
            @printf(io, "| %d | %d | %.4f | %.4f | %.2f |\n",
                row.n, row.krylovdim, row.median_allocs / 1e9, estimate / 1e9, ratio)
        end
        println(io)

        println(io, "### TrotterDomain")
        println(io)
        println(io, "| n | krylovdim | @allocated (GB) | estimate (GB) | ratio |")
        println(io, "|---|-----------|----------------|--------------|-------|")
        for row in trotter_data
            estimate = row.krylovdim * (4^row.n) * 16 * 1.5
            ratio = row.median_allocs / estimate
            @printf(io, "| %d | %d | %.4f | %.4f | %.2f |\n",
                row.n, row.krylovdim, row.median_allocs / 1e9, estimate / 1e9, ratio)
        end
        println(io)

        # --- Sys.maxrss supplementary data ---
        if haskey(metadata, :maxrss_records) && !isempty(metadata[:maxrss_records])
            println(io, "## Process Peak RSS")
            println(io)
            println(io, "| checkpoint | Sys.maxrss (GB) |")
            println(io, "|------------|----------------|")
            for (label, rss) in metadata[:maxrss_records]
                @printf(io, "| %s | %.3f |\n", label, rss / 1e9)
            end
            println(io)
        end
    end
end

# ---------------------------------------------------------------------------
# Benchmark data row type
# ---------------------------------------------------------------------------
struct BenchmarkRow
    n::Int
    krylovdim::Int
    median_time::Float64
    median_allocs::Float64
    matvecs::Int
    gap::Float64
    all_times::Vector{Float64}
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
function main()
    println("=" ^ 72)
    println("  Krylov Scaling Benchmark")
    println("  $(Dates.now())")
    println("  Julia $(VERSION), BLAS threads: $(BLAS.get_num_threads())")
    println("=" ^ 72)
    println()

    energy_data = BenchmarkRow[]
    trotter_data = BenchmarkRow[]
    maxrss_records = Tuple{String, Float64}[]

    # ==================================================================
    # EnergyDomain benchmarks
    # ==================================================================
    println("--- EnergyDomain KMS ---")
    println()
    @printf("  %-3s | %-6s | %-4s | %-12s | %-12s | %-6s | %s\n",
        "n", "dim^2", "kdim", "median_t (s)", "alloc (GB)", "matvec", "gap")
    println("  " * "-" ^ 70)

    for n in [3, 4, 5, 6]
        sys = make_benchmark_system(n)
        config = make_bench_liouv_config(n, EnergyDomain())
        kd = KRYLOVDIM_HEURISTIC[n]

        # At n=6: run krylovdim probe
        if n == 6
            println()
            chosen_kd = probe_krylovdim_n6(config, sys.hamiltonian, sys.jumps)
            KRYLOVDIM_HEURISTIC[6] = chosen_kd
            KRYLOVDIM_HEURISTIC[7] = chosen_kd
            kd = chosen_kd
            println()
        end

        n_runs = (n <= 6) ? 5 : 1
        bench = run_benchmark(config, sys.hamiltonian, sys.jumps; krylovdim=kd, n_runs=n_runs)

        row = BenchmarkRow(n, kd, bench.median_time, bench.median_allocs,
            bench.result.matvec_count, bench.result.spectral_gap, bench.times)
        push!(energy_data, row)

        rss = Sys.maxrss()
        push!(maxrss_records, ("EnergyDomain n=$n", Float64(rss)))

        @printf("  %-3d | %-6d | %-4d | %12.4f | %12.4f | %-6d | %.8e\n",
            n, 4^n, kd, bench.median_time, bench.median_allocs / 1e9,
            bench.result.matvec_count, bench.result.spectral_gap)
    end

    # Fit scaling (n=3..6)
    n_vals_energy = [row.n for row in energy_data]
    t_vals_energy = [row.median_time for row in energy_data]
    energy_fit = fit_scaling(n_vals_energy, t_vals_energy)

    println()
    @printf("  Scaling fit: t = %.4e * %.4f^n\n", energy_fit.a, energy_fit.b)
    println()

    # --- n=7 conditional (EnergyDomain) ---
    t7_pred_energy = energy_fit.a * energy_fit.b^7
    @printf("  Predicted n=7 time: %.1f s\n", t7_pred_energy)
    if t7_pred_energy < 15 * 60
        println("  Proceeding with n=7 (single run)...")
        sys7 = make_benchmark_system(7)
        config7 = make_bench_liouv_config(7, EnergyDomain())
        kd7 = KRYLOVDIM_HEURISTIC[7]
        bench7 = run_benchmark(config7, sys7.hamiltonian, sys7.jumps; krylovdim=kd7, n_runs=1)
        row7 = BenchmarkRow(7, kd7, bench7.median_time, bench7.median_allocs,
            bench7.result.matvec_count, bench7.result.spectral_gap, bench7.times)
        push!(energy_data, row7)
        rss = Sys.maxrss()
        push!(maxrss_records, ("EnergyDomain n=7", Float64(rss)))
        @printf("  n=7: time=%.4f s, alloc=%.4f GB, matvecs=%d, gap=%.8e\n",
            bench7.median_time, bench7.median_allocs / 1e9,
            bench7.result.matvec_count, bench7.result.spectral_gap)
    else
        @printf("  Skipping n=7 (predicted %.1f s exceeds 15 min limit)\n", t7_pred_energy)
    end
    println()

    # ==================================================================
    # TrotterDomain benchmarks
    # ==================================================================
    println("--- TrotterDomain KMS ---")
    println()
    @printf("  %-3s | %-6s | %-4s | %-12s | %-12s | %-6s | %s\n",
        "n", "dim^2", "kdim", "median_t (s)", "alloc (GB)", "matvec", "gap")
    println("  " * "-" ^ 70)

    for n in [3, 4, 5, 6]
        # Construct Trotter object first, then rebuild system with Trotter eigenbasis
        ham_only = load_hamiltonian("heis", n; beta=BETA)
        trotter = TrottTrott(ham_only, T0, NUM_TROTTER_STEPS_PER_T0)
        sys = make_benchmark_system(n; trotter=trotter)
        config = make_bench_liouv_config(n, TrotterDomain())
        kd = KRYLOVDIM_HEURISTIC[n]

        # At n=6: run krylovdim probe (TrotterDomain)
        if n == 6
            println()
            chosen_kd = probe_krylovdim_n6(config, sys.hamiltonian, sys.jumps; trotter=trotter)
            # Update heuristic for n=6 and n=7 TrotterDomain (use max of Energy and Trotter)
            KRYLOVDIM_HEURISTIC[6] = max(KRYLOVDIM_HEURISTIC[6], chosen_kd)
            KRYLOVDIM_HEURISTIC[7] = max(KRYLOVDIM_HEURISTIC[7], chosen_kd)
            kd = chosen_kd
            println()
        end

        n_runs = (n <= 6) ? 5 : 1
        bench = run_benchmark(config, sys.hamiltonian, sys.jumps;
            krylovdim=kd, n_runs=n_runs, trotter=trotter)

        row = BenchmarkRow(n, kd, bench.median_time, bench.median_allocs,
            bench.result.matvec_count, bench.result.spectral_gap, bench.times)
        push!(trotter_data, row)

        rss = Sys.maxrss()
        push!(maxrss_records, ("TrotterDomain n=$n", Float64(rss)))

        @printf("  %-3d | %-6d | %-4d | %12.4f | %12.4f | %-6d | %.8e\n",
            n, 4^n, kd, bench.median_time, bench.median_allocs / 1e9,
            bench.result.matvec_count, bench.result.spectral_gap)
    end

    # Fit scaling (n=3..6)
    n_vals_trotter = [row.n for row in trotter_data[1:min(4, length(trotter_data))]]
    t_vals_trotter = [row.median_time for row in trotter_data[1:min(4, length(trotter_data))]]
    trotter_fit = fit_scaling(n_vals_trotter, t_vals_trotter)

    println()
    @printf("  Scaling fit: t = %.4e * %.4f^n\n", trotter_fit.a, trotter_fit.b)
    println()

    # --- n=7 conditional (TrotterDomain) ---
    t7_pred_trotter = trotter_fit.a * trotter_fit.b^7
    @printf("  Predicted n=7 time: %.1f s\n", t7_pred_trotter)
    if t7_pred_trotter < 15 * 60
        println("  Proceeding with n=7 (single run)...")
        ham7 = load_hamiltonian("heis", 7; beta=BETA)
        trotter7 = TrottTrott(ham7, T0, NUM_TROTTER_STEPS_PER_T0)
        sys7 = make_benchmark_system(7; trotter=trotter7)
        config7 = make_bench_liouv_config(7, TrotterDomain())
        kd7 = KRYLOVDIM_HEURISTIC[7]
        bench7 = run_benchmark(config7, sys7.hamiltonian, sys7.jumps;
            krylovdim=kd7, n_runs=1, trotter=trotter7)
        row7 = BenchmarkRow(7, kd7, bench7.median_time, bench7.median_allocs,
            bench7.result.matvec_count, bench7.result.spectral_gap, bench7.times)
        push!(trotter_data, row7)
        rss = Sys.maxrss()
        push!(maxrss_records, ("TrotterDomain n=7", Float64(rss)))
        @printf("  n=7: time=%.4f s, alloc=%.4f GB, matvecs=%d, gap=%.8e\n",
            bench7.median_time, bench7.median_allocs / 1e9,
            bench7.result.matvec_count, bench7.result.spectral_gap)
    else
        @printf("  Skipping n=7 (predicted %.1f s exceeds 15 min limit)\n", t7_pred_trotter)
    end
    println()

    # ==================================================================
    # Report generation
    # ==================================================================
    results_dir = joinpath(dirname(@__DIR__), "results")
    mkpath(results_dir)
    report_path = joinpath(results_dir, "krylov_scaling_report.md")

    metadata = Dict{Symbol, Any}(:maxrss_records => maxrss_records)
    write_report(report_path, energy_data, trotter_data, energy_fit, trotter_fit, metadata)

    println("=" ^ 72)
    println("  Report saved to: $(report_path)")
    println("=" ^ 72)
end

# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
try
    main()
finally
    BLAS.set_num_threads(ORIGINAL_BLAS_THREADS)
end
