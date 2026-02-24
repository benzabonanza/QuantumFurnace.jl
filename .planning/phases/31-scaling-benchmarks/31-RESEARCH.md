# Phase 31: Scaling Benchmarks - Research

**Researched:** 2026-02-24
**Domain:** Empirical timing and memory benchmarking of Krylov spectral gap estimation at n=3-7, scaling law fitting, and extrapolation to n=10,12
**Confidence:** HIGH

## Summary

Phase 31 creates a standalone benchmark script in `simulations/` that measures wall-clock timing and peak memory for `krylov_spectral_gap()` at system sizes n=3,4,5,6 (and n=7 conditionally), fits a power-law model `t = a * b^n` to confirm the expected 4^n scaling, and extrapolates resource requirements for n=10 and n=12 production runs. The script targets EnergyDomain and TrotterDomain with KMS balance (coherent term enabled), and produces a markdown report saved to `results/`.

The codebase already has all required infrastructure: `krylov_spectral_gap()` (Phase 29), cross-validated against dense at n=4 and n=6 (Phase 30), `LsqFit.jl` for curve fitting (already a dependency), pre-serialized Hamiltonians for n=3 through n=10, and factory functions for creating test systems at arbitrary n. The benchmark script needs to: (1) construct systems at each n, (2) run `krylov_spectral_gap()` with carefully chosen krylovdim heuristics, (3) measure wall-clock time and peak memory, (4) fit the scaling law, (5) assert the base is in [3.5, 4.5], and (6) generate the report.

Key design decisions from CONTEXT.md override the original phase spec: total eigsolve wall-clock time only (no per-matvec breakdown), and the krylovdim probe at n=6 against dense reference is the critical calibration step that determines timing data feeding into extrapolation.

**Primary recommendation:** Use `@elapsed` for wall-clock timing after JIT warmup, `@allocated` for total allocation measurement, and `Sys.maxrss()` for peak RSS tracking. LsqFit.jl's `curve_fit` handles the power-law regression. The script should follow the existing `simulations/main_liouv.jl` pattern as a standalone Julia file.

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

#### Benchmark scope
- Domains: EnergyDomain and TrotterDomain only, both with KMS balance (coherent term enabled)
- System sizes: n=3,4,5,6 guaranteed; n=7 conditional on extrapolated time fitting under 15 minutes
- Repetitions: 3-5 runs with median for n=3,4,5,6; single run for n=7 if attempted
- Execution: Standalone script in `simulations/` folder (fourth simulation kind alongside Lindbladian construction, DM, trajectory)

#### Output & reporting
- Print summary tables to stdout during execution
- Save markdown report to `results/` folder (top-level, same level as `simulations/`)
- Report includes extrapolation with time + memory + feasibility notes for n=10,12
- Separate sections for EnergyDomain and TrotterDomain (not side-by-side)

#### Breakdown granularity
- Total eigsolve wall-clock time only (no per-matvec breakdown)
- Peak memory allocation measured at each system size
- Memory guard thresholds to be calibrated from empirical data (current guards were set blindly)

#### Extrapolation approach
- Power-law fit: t = a * b^n, expect b ~ 4
- Hard assertion that b is in [3.5, 4.5] -- script errors if scaling is wrong
- Separate fits for EnergyDomain and TrotterDomain
- Report predicted time, memory, and feasibility notes for n=10 and n=12

#### Krylovdim selection
- Target Krylov precision: 1e-8 absolute (well below 1e-6 physical errors from Trotter/quadrature)
- No krylovdim sweep -- use heuristic values that increase with n
- At n=6: probe krylovdim=50 against dense reference for 1e-8 precision; if insufficient, try krylovdim=100 (memory permitting). Chosen krylovdim feeds into timing data and extrapolation
- Be memory-conscious at n=6+ (laptop RAM constraints)
- Only include krylovdim-vs-n recommendation table if research supports it with prior experience from the literature

### Claude's Discretion
- Exact krylovdim heuristic for n=3,4,5 (small systems, less sensitive)
- Memory measurement approach (peak allocation tracking method)
- Warmup strategy (JIT compilation handling)
- Report formatting and section structure
- How to handle n=7 time extrapolation logic

### Deferred Ideas (OUT OF SCOPE)
None -- discussion stayed within phase scope
</user_constraints>

## Note: BENCH-04 Override

The original phase spec (ROADMAP.md) includes Success Criterion #4: "Per-matvec timing breakdown separates BLAS matrix multiplication cost from precomputation lookup and Krylov iteration overhead" (BENCH-04). However, the CONTEXT.md locked decision explicitly states: "Total eigsolve wall-clock time only (no per-matvec breakdown)." **The CONTEXT.md decision overrides BENCH-04.** The planner should NOT include per-matvec timing instrumentation.

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| QuantumFurnace.jl | local | `krylov_spectral_gap()`, `extract_leading_eigendata()`, config/system factories | The code under benchmark |
| LsqFit.jl | 0.15 | Power-law curve fitting via Levenberg-Marquardt | Already a project dependency; used in v1.3 for exponential fitting |
| LinearAlgebra (stdlib) | -- | `BLAS.set_num_threads()` for controlled BLAS thread count | Already used throughout; benchmarks must control BLAS threads |
| Printf (stdlib) | -- | Formatted output tables | Already used in simulations and cross-validation tests |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| Base.Sys | stdlib | `Sys.free_memory()`, `Sys.maxrss()` | Memory measurement: free memory before runs, peak RSS after |
| Base | stdlib | `@elapsed`, `@allocated`, `GC.gc()` | Timing, allocation measurement, GC control between runs |
| Statistics (stdlib) | -- | `median()` for aggregating repeated measurements | Taking median of 3-5 timing runs |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| `@elapsed` | BenchmarkTools.jl `@belapsed` | BenchmarkTools runs many iterations to find minimum; unnecessary for benchmarks taking seconds to minutes. `@elapsed` after warmup is simpler and sufficient. |
| `Sys.maxrss()` | `/proc/self/status` parsing | `Sys.maxrss()` is Julia-native, cross-platform, and already available. Process-level RSS is cumulative (not per-call), but useful for peak tracking. |
| `@allocated` | Manual `Base.gc_live_bytes()` tracking | `@allocated` is the standard Julia idiom for total allocation measurement; `gc_live_bytes` has known accuracy issues. |

**No new dependencies needed.** All tools are either already in Project.toml or Julia stdlib.

## Architecture Patterns

### Recommended Script Structure
```
simulations/
    main_krylov_benchmark.jl    # [NEW] Standalone benchmark script
results/
    krylov_scaling_report.md    # [NEW] Generated markdown report
```

### Pattern 1: System Factory for Arbitrary n
**What:** Generalize the existing `make_n6_test_system()` pattern (from Phase 30) to arbitrary n, loading pre-serialized Hamiltonians and constructing jump operators.
**When to use:** Each benchmark size n=3,4,5,6,(7) needs a system.
**Example:**
```julia
function make_benchmark_system(n_qubits::Int; trotter::Union{Nothing, TrottTrott}=nothing)
    source_root = dirname(@__DIR__)  # or use Pkg.project() in standalone scripts
    ham_path = joinpath(source_root, "hamiltonians",
        "heis_disordered_periodic_n$(n_qubits).bson")
    hamiltonian = load_hamiltonian_bson(ham_path, BETA)

    jump_paulis = [[X], [Y], [Z]]
    n_jumps = 3 * n_qubits
    jump_norm = sqrt(n_jumps)
    basis = trotter !== nothing ? trotter.eigvecs : hamiltonian.eigvecs

    jumps = JumpOp[]
    for pauli in jump_paulis, site in 1:n_qubits
        op = Matrix(pad_term(pauli, n_qubits, site)) ./ jump_norm
        op_eigen = basis' * op * basis
        push!(jumps, JumpOp(op, op_eigen, op == transpose(op), op == op'))
    end
    return (; hamiltonian, jumps, n_qubits)
end
```
**Key note:** The script is standalone (not in test suite), so it must load Hamiltonians using the project-relative path, not `Pkg.project().path` which works differently outside `Pkg.test()`. Use `@__DIR__` to resolve paths relative to the script location.

### Pattern 2: Warmup + Timed Measurement
**What:** Run the benchmark once to trigger JIT compilation, then measure subsequent runs.
**When to use:** Every timing measurement.
**Example:**
```julia
function benchmark_krylov(config, ham, jumps; krylovdim, howmany=4, n_runs=5, trotter=nothing)
    # Warmup: trigger JIT compilation for this specific config/system combination
    GC.gc()
    _ = krylov_spectral_gap(config, ham, jumps;
        trotter=trotter, krylovdim=krylovdim, howmany=howmany)

    # Timed runs
    times = Float64[]
    allocs = Int[]
    for _ in 1:n_runs
        GC.gc()
        t = @elapsed begin
            result = krylov_spectral_gap(config, ham, jumps;
                trotter=trotter, krylovdim=krylovdim, howmany=howmany)
        end
        a = @allocated krylov_spectral_gap(config, ham, jumps;
            trotter=trotter, krylovdim=krylovdim, howmany=howmany)
        push!(times, t)
        push!(allocs, a)
    end
    return (; times, allocs, median_time=median(times), median_allocs=median(allocs))
end
```
**Critical note on warmup:** Julia's JIT compilation is per-method-specialization. A warmup at n=3 does NOT warm up n=4 (different matrix dimensions trigger different specializations of `mul!`, etc.). However, the workspace construction and closure setup share common code paths, so warmup at the first small size partially warms the compiler. The safest approach is: for each n, run one untimed call before the timed repetitions.

### Pattern 3: Power-Law Fitting with LsqFit
**What:** Fit `t = a * b^n` to timing data using nonlinear least squares.
**When to use:** After collecting timing data at all system sizes.
**Example:**
```julia
using LsqFit

model(n, p) = p[1] .* (p[2] .^ n)
n_values = Float64.([3, 4, 5, 6])
t_values = [t3, t4, t5, t6]  # median times

p0 = [1e-3, 4.0]  # initial guess: a ~ small, b ~ 4
fit = curve_fit(model, n_values, t_values, p0)
a, b = fit.param

# Hard assertion per locked decision
@assert 3.5 <= b <= 4.5 "Scaling base b=$b outside [3.5, 4.5]; expected ~4"

# Extrapolate
t10_predicted = a * b^10
t12_predicted = a * b^12
```
**Key considerations:**
- The model `t = a * b^n` is equivalent to `log(t) = log(a) + n*log(b)`, so one could also fit a linear model to `(n, log(t))`. However, LsqFit's nonlinear fit is more robust when timing data spans many orders of magnitude.
- The initial guess `p0 = [1e-3, 4.0]` should be reasonable for the expected scaling.
- LsqFit uses Levenberg-Marquardt and returns `fit.param`, `fit.resid`, and supports `confidence_interval(fit, alpha)` for uncertainty quantification.

### Pattern 4: krylovdim Probe at n=6
**What:** At n=6, compare Krylov gap with krylovdim=50 against dense reference to verify 1e-8 precision. If insufficient, try krylovdim=100.
**When to use:** The n=6 benchmark run, which is the critical calibration step.
**Example:**
```julia
# Dense reference at n=6
L_dense = construct_lindbladian(jumps, config, ham)
dense_result = extract_leading_eigendata(L_dense; n_modes=4)
dense_gap = dense_result.spectral_gap

# Probe krylovdim=50
krylov_50 = krylov_spectral_gap(config, ham, jumps; krylovdim=50, howmany=4)
err_50 = abs(krylov_50.spectral_gap - dense_gap)
if err_50 < 1e-8
    chosen_krylovdim = 50
else
    # Try krylovdim=100
    krylov_100 = krylov_spectral_gap(config, ham, jumps; krylovdim=100, howmany=4)
    err_100 = abs(krylov_100.spectral_gap - dense_gap)
    @assert err_100 < 1e-8 "krylovdim=100 insufficient at n=6: err=$err_100"
    chosen_krylovdim = 100
end
# Use chosen_krylovdim for the timed n=6 benchmark runs
```
**Important:** The dense reference `construct_lindbladian` at n=6 builds a 4096x4096 matrix. This is feasible (Phase 30 already does it in the cross-validation tests) but takes a few seconds. The dense computation is NOT timed as part of the Krylov benchmark.

### Pattern 5: n=7 Conditional Execution
**What:** Before running n=7, extrapolate timing from n=3-6 data to estimate n=7 time. Only proceed if the estimate is under 15 minutes.
**When to use:** After completing n=3-6 benchmarks.
**Example:**
```julia
# After fitting t = a * b^n to n=3..6 data
t7_predicted = a * b^7
if t7_predicted < 15 * 60  # 15 minutes
    @info "Predicted n=7 time: $(round(t7_predicted; digits=1))s -- proceeding with single run"
    # Run n=7 (single run, no repetitions per locked decision)
    GC.gc()
    t7 = @elapsed result_7 = krylov_spectral_gap(config, ham7, jumps7; krylovdim=krylovdim_7)
else
    @info "Predicted n=7 time: $(round(t7_predicted; digits=1))s -- skipping (exceeds 15min limit)"
    t7 = nothing
end
```

### Anti-Patterns to Avoid
- **Running benchmark without warmup:** First call includes JIT compilation time (can be 10-100x slower). Always run one untimed call per system size first.
- **Using BenchmarkTools for long-running benchmarks:** `@benchmark` tries to run multiple iterations to find a stable minimum. For benchmarks taking seconds to minutes, this wastes time. `@elapsed` with explicit repetitions and median is the right approach.
- **Not calling GC.gc() between runs:** Residual memory from previous runs can affect both timing (GC pauses) and memory measurements. Always GC between benchmark runs.
- **Fitting in log-space without weighting:** If fitting `log(t) = log(a) + n*log(b)` via linear regression, the fit is unweighted in log-space, meaning small-n data points (with smaller absolute errors) get disproportionate weight. The direct nonlinear fit `t = a * b^n` weights by absolute error, which is more appropriate when all data points are equally important.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Power-law curve fitting | Manual log-linear regression | `LsqFit.curve_fit(model, n_values, t_values, p0)` | Handles nonlinear models, provides residuals and confidence intervals |
| System construction at arbitrary n | Copy-paste from test_helpers.jl | Generalize `make_n6_test_system()` pattern | The test file already has the proven pattern; just parameterize on n |
| Hamiltonian loading in standalone script | Re-implement BSON parsing | Use `_load_hamiltonian_bson` from misc_tools.jl (via `load_hamiltonian`) or `_load_test_hamiltonian` pattern from test_helpers.jl | Both patterns are proven; the script needs to handle the path resolution carefully |
| Report markdown generation | String concatenation | Use `Printf.@sprintf` and `open(..., "w") do io; println(io, ...) end` | Standard Julia I/O pattern |

**Key insight:** This phase is primarily a data collection and analysis script, not new library code. The pattern is: call existing APIs, measure, fit, report. Almost no new "infrastructure" is needed.

## Common Pitfalls

### Pitfall 1: Hamiltonian Path Resolution in Standalone Scripts
**What goes wrong:** `load_hamiltonian("heis", n; beta=beta)` uses `Pkg.project().path` to find the hamiltonians directory. In standalone scripts invoked as `julia --project simulations/main_krylov_benchmark.jl`, `Pkg.project().path` points to the project root, which works correctly. But if run from a different directory, the path resolution can fail.
**Why it happens:** `Pkg.project().path` returns the `Project.toml` path, and `dirname` gives the project root. The hamiltonians directory is at `project_root/hamiltonians/`. This works when `--project` points to the right place.
**How to avoid:** Use `load_hamiltonian("heis", n; beta=beta)` directly, since the script will be run with `julia --project=@. simulations/main_krylov_benchmark.jl` from the project root. Alternatively, use `@__DIR__` and navigate up: `joinpath(dirname(@__DIR__), "hamiltonians", ...)`.
**Warning signs:** `SystemError: opening file "..." : No such file or directory` at Hamiltonian loading.

### Pitfall 2: JIT Warmup Only Covers Specific Type Specializations
**What goes wrong:** Warming up `krylov_spectral_gap` at n=3 (dim=8, dim^2=64) does NOT warm up the method for n=4 (dim=16, dim^2=256). The `mul!` specializations for different matrix sizes are compiled separately.
**Why it happens:** Julia's JIT compiler specializes on concrete types AND sizes of inputs. A `mul!(C::Matrix{ComplexF64}, A::Matrix{ComplexF64}, B::Matrix{ComplexF64})` call may not trigger recompilation for different sizes, but the internal BLAS dispatch and workspace allocation can differ.
**How to avoid:** For each n, run one untimed warmup call before the timed repetitions. This ensures all code paths for that specific system size are compiled. The warmup call is fast for small n (n=3,4) and non-trivial for large n (n=6,7), but the alternative (including compile time in measurements) is worse.
**Warning signs:** First run at each n is 2-10x slower than subsequent runs.

### Pitfall 3: Memory Measurement Includes Previous Allocations
**What goes wrong:** `Sys.maxrss()` returns the peak RSS for the entire Julia process, not for a single function call. If n=5 caused a peak RSS of 2 GB and n=6 needs 4 GB, the n=6 `Sys.maxrss()` reads correctly. But if n=6 needs only 1.5 GB, `Sys.maxrss()` still reads 2 GB (from the n=5 peak).
**Why it happens:** RSS is a monotonically non-decreasing process-level metric. It tracks the lifetime peak, not per-call peaks.
**How to avoid:** For each system size n, record `@allocated` which tracks total bytes allocated by a specific expression (heap allocations during that call). This is per-call, not cumulative. Also record `Sys.maxrss()` as a running baseline to confirm the process is not running out of memory. For memory extrapolation, use `@allocated` data (which is per-call and reproducible) rather than `Sys.maxrss()` (which is cumulative).
**Warning signs:** Memory measurements that are flat or non-monotonic across increasing n.

### Pitfall 4: NUFFTPrefactors Dominate Memory for TrotterDomain
**What goes wrong:** The NUFFT prefactor array `NUFFTPrefactors.data` has shape `(dim, dim, n_energy_labels)` where `dim = 2^n`. For n=6, dim=64, so the prefactors are `64 x 64 x n_labels * 16 bytes`. With num_energy_bits=12 and truncation, n_labels can be ~2049. That is 64*64*2049*16 = 134 MB for one domain. At n=7 (dim=128), it would be 128*128*2049*16 = 537 MB.
**Why it happens:** NUFFT precomputation stores a full dim x dim complex matrix for each energy label. This scales as O(4^n * n_energy_labels).
**How to avoid:** Be aware that TrotterDomain (and TimeDomain) benchmarks include NUFFT precomputation memory that EnergyDomain does not have. The benchmark should separately report workspace construction time and eigsolve time (but NOT per-matvec breakdown per CONTEXT.md). For memory, the key metric is total `@allocated` for `krylov_spectral_gap()`, which includes workspace construction.
**Warning signs:** TrotterDomain uses significantly more memory than EnergyDomain at the same n.

### Pitfall 5: BLAS Thread Count Affects Timing Reproducibility
**What goes wrong:** If BLAS thread count is not explicitly set, Julia defaults to the number of physical cores. Different machines produce different timing ratios (multi-threaded vs single-threaded BLAS scales differently with matrix size).
**Why it happens:** The locked decision specifies "4 BLAS threads" (from the original spec). BLAS parallelism affects the `mul!` calls in the matvec hot path.
**How to avoid:** Set `BLAS.set_num_threads(4)` at the start of the benchmark script. Record the BLAS thread count in the report. Restore the original thread count in a `try-finally` block if needed.
**Warning signs:** Timing numbers that vary wildly between runs due to background processes competing for BLAS threads.

### Pitfall 6: GC Pauses Contaminate Timing
**What goes wrong:** The Krylov eigsolve allocates one `dim^2`-length vector per matvec call (the `copy(vec(ws.rho_out))` in the closure). For n=6, dim^2 = 4096, so each vector is 32 KB. With krylovdim=50 and maxiter=100, that is up to 5000 allocations x 32 KB = 160 MB of garbage per eigsolve. GC can trigger during the timed section.
**Why it happens:** Julia's GC is stop-the-world. Large allocation rates trigger GC pauses that are included in `@elapsed` measurements.
**How to avoid:** Call `GC.gc()` before each timed run to start with a clean heap. Take the median of multiple runs to smooth out GC variability. Do NOT try to disable GC during benchmarks (it will cause OOM).
**Warning signs:** One outlier run that is 2-5x slower than the median.

## Code Examples

### LsqFit Power-Law Model
```julia
# Source: LsqFit.jl documentation (https://julianlsolvers.github.io/LsqFit.jl/latest/tutorial/)
using LsqFit

# Model: t = a * b^n
model(n, p) = p[1] .* (p[2] .^ n)

# Data
n_data = Float64.([3, 4, 5, 6])
t_data = [0.05, 0.2, 0.8, 3.2]  # hypothetical median times in seconds

# Fit
p0 = [1e-3, 4.0]
fit = curve_fit(model, n_data, t_data, p0)
a_fit, b_fit = fit.param

# Confidence intervals (95%)
ci = confidence_interval(fit, 0.05)
b_lo, b_hi = ci[2]  # confidence interval for b

# Extrapolate
t10 = a_fit * b_fit^10
t12 = a_fit * b_fit^12
```

### Julia Timing and Memory Measurement
```julia
# Source: Julia Base documentation

# Wall-clock timing (excludes compilation after warmup)
GC.gc()
t = @elapsed begin
    result = expensive_function(args...)
end

# Total heap allocations in bytes
a = @allocated expensive_function(args...)

# Peak process RSS (cumulative, not per-call)
peak_rss = Sys.maxrss()

# Free memory available
free_mem = Sys.free_memory()
```

### BLAS Thread Control
```julia
# Source: LinearAlgebra documentation
using LinearAlgebra

# Set BLAS threads for reproducible benchmarks
original_threads = BLAS.get_num_threads()
BLAS.set_num_threads(4)
try
    # ... run benchmarks ...
finally
    BLAS.set_num_threads(original_threads)
end
```

### Dense Reference at n=6 for krylovdim Probe
```julia
# Source: existing Phase 30 cross-validation pattern
using QuantumFurnace

# Build dense Lindbladian (4096x4096 at n=6, feasible in seconds)
L_dense = construct_lindbladian(jumps, config, ham)
dense_result = extract_leading_eigendata(L_dense; n_modes=4)
dense_gap = dense_result.spectral_gap

# Krylov probe
krylov_result = krylov_spectral_gap(config, ham, jumps; krylovdim=50, howmany=4)
precision = abs(krylov_result.spectral_gap - dense_gap)
```

### Markdown Report Generation
```julia
open(report_path, "w") do io
    println(io, "# Krylov Scaling Benchmark Report")
    println(io, "")
    println(io, "**Generated:** $(Dates.now())")
    println(io, "**BLAS threads:** 4")
    println(io, "")
    println(io, "## EnergyDomain KMS")
    println(io, "")
    println(io, "| n | dim^2 | krylovdim | median_time (s) | alloc (GB) | matvecs |")
    println(io, "|---|-------|-----------|-----------------|------------|---------|")
    for row in energy_rows
        @printf(io, "| %d | %d | %d | %.3f | %.3f | %d |\n",
            row.n, 4^row.n, row.krylovdim, row.median_time,
            row.median_alloc / 1e9, row.matvecs)
    end
end
```

## Discretion Area Resolutions

### Krylovdim Heuristic for n=3,4,5

**Recommendation:** Use krylovdim=30 for n=3,4,5 (matching the default in `krylov_spectral_gap`).

**Rationale:**
- At n=4 (dim^2=256), krylovdim=30 with tol=1e-10 converges reliably (verified in Phase 29 and Phase 30 tests). The cross-validation achieves atol < 1e-8 against dense.
- At n=3 (dim^2=64), krylovdim=30 is overkill (30/64 = 47% of the full space), but it is fast and consistent.
- At n=5 (dim^2=1024), krylovdim=30 should still converge to 1e-8 precision since the spectrum is well-separated (the spectral gap is O(1) while the smallest eigenvalue magnitude is O(100)).
- At n=6 (dim^2=4096), the krylovdim probe (50 then 100) determines the correct value.
- At n=7 (dim^2=16384), extrapolate from n=6: if n=6 needs krylovdim=50, use 50 for n=7; if n=6 needs 100, use 100 for n=7.

**Heuristic table:**
| n | dim^2 | krylovdim |
|---|-------|-----------|
| 3 | 64 | 30 |
| 4 | 256 | 30 |
| 5 | 1024 | 30 |
| 6 | 4096 | 50 (probe) or 100 |
| 7 | 16384 | same as n=6 |

### Memory Measurement Approach

**Recommendation:** Use `@allocated` for per-call total allocation and `Sys.maxrss()` for process-level peak RSS.

**Details:**
- `@allocated` measures total heap allocations during a specific expression. This is the primary metric for memory scaling analysis because it is per-call and reproducible.
- `Sys.maxrss()` provides peak resident set size for the Julia process. Useful as a sanity check (did we OOM?) but not reliable for per-call memory estimation because it is cumulative.
- For the memory guard calibration, compare `@allocated` from empirical runs against the pre-flight estimate `krylovdim * 4^n * 16 * 1.5`. If `@allocated` is significantly higher, the 1.5x factor needs adjustment.

**Strategy per n:**
1. Call `GC.gc()` before measurement
2. `alloc = @allocated krylov_spectral_gap(config, ham, jumps; krylovdim=kd)` -- this captures workspace construction + all Krylov iterations
3. Record `Sys.maxrss()` after the call as a running baseline

### Warmup Strategy

**Recommendation:** Per-n warmup with one untimed call at each system size.

**Details:**
1. At the very start of the script, run one `krylov_spectral_gap` call at n=3 to trigger initial package compilation (module loading, method compilation).
2. For each subsequent n, run one untimed call before the timed repetitions. This ensures type-specialized compilation (different dim values trigger different BLAS routines) is excluded from timing.
3. The warmup call's result is discarded.

**Timing impact:** The warmup at n=3 takes a few seconds (JIT). Subsequent warmups at n=4,5,6 take less time (most code paths are already compiled; only size-specific specializations need recompilation). The warmup at n=6 might take 10-30 seconds due to workspace construction + NUFFT precomputation.

### Report Formatting

**Recommendation:** Generate a markdown file with the following structure:

```
# Krylov Scaling Benchmark Report
## Metadata (date, BLAS threads, Julia version, system info)
## EnergyDomain KMS
### Timing Results (table: n, dim^2, krylovdim, median_time, alloc, matvecs)
### Scaling Fit (a, b, R^2)
## TrotterDomain KMS
### Timing Results (same table format)
### Scaling Fit (a, b, R^2)
## Extrapolation
### Predicted times for n=10 and n=12
### Memory extrapolation
### Feasibility notes
## Memory Guard Calibration
### Empirical vs pre-flight estimates
### Recommended threshold adjustment
```

### n=7 Conditional Logic

**Recommendation:** After fitting the n=3-6 data, compute `t7_predicted = a * b^7`. If `t7_predicted < 15 * 60` seconds (15 minutes), run a single timed iteration at n=7. Use the same krylovdim as n=6. Include the n=7 data point in the report (but NOT in the scaling fit, since it is a single run without repetitions).

If n=7 is attempted and succeeds, re-fit the scaling model including n=7 data and report both fits (with and without n=7) for comparison.

## Memory Budget Analysis

### Krylov Vector Storage (dominates at large n)
KrylovKit stores `krylovdim` vectors of length `dim^2 = 4^n`, each ComplexF64 (16 bytes):
- n=6, krylovdim=50: 50 * 4096 * 16 = 3.3 MB
- n=7, krylovdim=50: 50 * 16384 * 16 = 13.1 MB
- n=10, krylovdim=50: 50 * 1048576 * 16 = 838 MB
- n=12, krylovdim=50: 50 * 16777216 * 16 = 13.4 GB

### Workspace Scratch Matrices
`KrylovWorkspace` allocates 5 scratch matrices (jump_oft, tmp1, tmp2, LdagL, rho_out) each `dim x dim`:
- n=6: 5 * 64 * 64 * 16 = 327 KB
- n=7: 5 * 128 * 128 * 16 = 1.3 MB
- n=10: 5 * 1024 * 1024 * 16 = 83.9 MB

### NUFFT Prefactors (TrotterDomain / TimeDomain only)
`NUFFTPrefactors.data` has shape `(dim, dim, n_energy_labels)`:
- n=6, 2049 labels: 64 * 64 * 2049 * 16 = 134 MB
- n=7, 2049 labels: 128 * 128 * 2049 * 16 = 537 MB
- n=10, 2049 labels: 1024 * 1024 * 2049 * 16 = 34.4 GB

**Critical observation:** For TrotterDomain, the NUFFT prefactors dominate memory at all sizes and grow as O(4^n * n_labels). At n=10, this alone exceeds typical cluster node memory (34 GB). This means TrotterDomain production runs at n=10,12 are likely infeasible without algorithmic changes (e.g., on-the-fly NUFFT evaluation instead of precomputation). EnergyDomain does not have NUFFT prefactors and is dominated by Krylov vectors at large n.

The benchmark report should note this asymmetry: EnergyDomain is feasible at n=10 (Krylov vectors ~838 MB with krylovdim=50) while TrotterDomain is not (NUFFT prefactors alone ~34 GB).

### Memory Guard Calibration
The current pre-flight estimate is `krylovdim * 4^n * 16 * 1.5`. This captures Krylov vector storage but NOT:
- Workspace scratch matrices (~5 * dim^2 * 16)
- NUFFT prefactors (TrotterDomain/TimeDomain)
- Jump operator eigenbasis copies (n_jumps * dim^2 * 16)

The benchmark should compare `@allocated` against the formula and recommend an updated factor if the current 1.5x is too low.

## Open Questions

1. **Exact n_energy_labels for truncated labels at each n**
   - What we know: `num_energy_bits=12` gives 4096 raw labels, truncated by `_truncate_energy_labels`. The truncation depends on the Hamiltonian's energy spectrum.
   - What's unclear: How many labels survive truncation at n=5,6,7. This affects NUFFT memory for TrotterDomain.
   - Recommendation: The benchmark script should print `length(precomputed_data.energy_labels)` for each n as part of the report metadata. This is diagnostic, not a blocking question.

2. **Whether krylovdim=30 suffices at n=5 for 1e-8 precision**
   - What we know: n=4 achieves 1e-8 with krylovdim=30. n=6 may need 50 or 100.
   - What's unclear: n=5 is between the two. If krylovdim=30 is insufficient at n=5, the timing data would include retry overhead (1.5x increase from _eigsolve_with_retry).
   - Recommendation: Run the n=5 benchmark with krylovdim=30 and check if `result.converged >= howmany`. If not, increase krylovdim to 50 for n=5. This is handled automatically by the retry logic but should be noted in the report.

3. **Whether Sys.maxrss() is available on all target platforms**
   - What we know: `Sys.maxrss()` wraps `getrusage(2)` on Linux/macOS and `GetProcessMemoryInfo` on Windows. Julia >= 1.11 supports it.
   - What's unclear: Whether it returns accurate values in container environments (Docker/sandbox).
   - Recommendation: Use `Sys.maxrss()` as supplementary data. Primary memory metric is `@allocated`.

## Sources

### Primary (HIGH confidence)
- Codebase: `src/krylov_eigsolve.jl` -- `krylov_spectral_gap()` API, `_check_krylov_memory()`, `KrylovGapResult`
- Codebase: `src/krylov_workspace.jl` -- `KrylovWorkspace` constructor, memory allocation patterns
- Codebase: `src/nufft.jl` -- `NUFFTPrefactors` struct, `_prepare_oft_nufft_prefactors()` memory footprint
- Codebase: `test/test_krylov_crossvalidation.jl` -- `make_n6_test_system()`, `compare_krylov_dense()`, Phase 30 patterns
- Codebase: `test/test_helpers.jl` -- `make_test_system()`, config factories, physical parameters
- Codebase: `simulations/main_liouv.jl` -- existing standalone simulation script pattern
- Codebase: `hamiltonians/*.bson` -- pre-serialized Hamiltonians available for n=3 through n=10
- Codebase: `Project.toml` -- LsqFit, BenchmarkTools already in dependencies/compat
- [LsqFit.jl Tutorial](https://julianlsolvers.github.io/LsqFit.jl/latest/tutorial/) -- `curve_fit`, `confidence_interval` API
- [BenchmarkTools.jl Manual](https://juliaci.github.io/BenchmarkTools.jl/stable/manual/) -- `@belapsed`, warmup strategy, memory measurement

### Secondary (MEDIUM confidence)
- [Julia Discourse: Peak memory](https://discourse.julialang.org/t/peak-memory/9856) -- `Sys.maxrss()` discussion and limitations
- [Julia Discourse: Tracking allocations](https://m3g.github.io/JuliaNotes.jl/stable/memory/) -- `@allocated` usage patterns
- [Julia Performance Tips](https://docs.julialang.org/en/v1/manual/performance-tips/) -- warmup, allocation tracking, JIT compilation behavior
- [Julia Sys source](https://github.com/JuliaLang/julia/blob/master/base/sysinfo.jl) -- `Sys.maxrss()`, `Sys.free_memory()` implementation

### Tertiary (LOW confidence)
- None

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- all tools (LsqFit, @elapsed, @allocated, Sys.maxrss) are well-documented Julia stdlib or existing dependencies
- Architecture: HIGH -- script follows existing simulation patterns in the codebase; system factories generalize proven Phase 30 patterns
- Memory analysis: HIGH -- computed from known data structure sizes in nufft.jl, krylov_workspace.jl
- Scaling law fitting: MEDIUM -- LsqFit power-law model is straightforward, but the [3.5, 4.5] assertion range depends on the actual problem structure (verified in this specific codebase but theoretical justification is problem-dependent)
- Pitfalls: HIGH -- identified from direct code analysis (JIT warmup, NUFFT memory, Sys.maxrss limitations)

**Research date:** 2026-02-24
**Valid until:** 2026-03-24 (stable codebase, standard tools)
