#!/usr/bin/env julia
#
# Investigate why n=4 bi-exponential mixing time estimation has higher error.
# Compare the decay dynamics and fitting quality for n=3,4,5.
#
# Usage:  OPENBLAS_NUM_THREADS=4 julia -t4 --project scripts/investigate_n4_fitting.jl

using QuantumFurnace
using LinearAlgebra
using Random
using Printf
using Statistics: mean, std

BLAS.set_num_threads(4)

beta   = 10.0
target = 1e-4
coeffs = fill(1.0, 3)
delta  = 0.0005  # Same delta for all, to compare apples-to-apples

for num_qubits in [3, 4, 5]
    println("\n" * "="^70)
    @printf("  n = %d\n", num_qubits)
    println("="^70)

    # Generate Hamiltonian with Z+ZZ (qf-yi4: deterministic seed, no selector)
    raw = build_heis_1d(num_qubits, coeffs;
        seed=42, periodic=true,
        disordering_terms=Vector{Matrix{ComplexF64}}[[Z], [Z, Z]],
        disorder_strength=1.0)
    hamiltonian = HamHam(raw, beta)

    dim = 2^num_qubits
    @printf("  nu_min = %.6e\n", raw.nu_min)

    # Bohr frequency analysis
    bohr = hamiltonian.bohr_freqs
    bohr_upper = [bohr[i,j] for i in 1:dim for j in (i+1):dim]
    sort!(bohr_upper)
    @printf("  Min Bohr gap: %.4e\n", minimum(diff(bohr_upper)))

    # Look at the eigenvalue spectrum
    eigvals_sorted = sort(hamiltonian.eigvals)
    gaps = diff(eigvals_sorted)
    @printf("  Eigenvalue gaps (sorted): min=%.4e, max=%.4e\n", minimum(gaps), maximum(gaps))
    println("  All eigenvalue gaps:")
    for (i, g) in enumerate(gaps)
        @printf("    gap[%d] = %.6e\n", i, g)
    end

    # Gibbs state weights (in eigenbasis)
    gibbs_diag = real.(diag(hamiltonian.gibbs))
    println("  Gibbs weights (eigenbasis):")
    for (i, w) in enumerate(gibbs_diag)
        @printf("    w[%d] = %.6e\n", i, w)
    end
    @printf("  Max/min Gibbs ratio: %.2f\n", maximum(gibbs_diag) / minimum(gibbs_diag))

    # Build jumps
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

    sigma = 1.0 / beta
    w0    = 0.05
    num_energy_bits = 12
    t0    = 2π / (2^num_energy_bits * w0)

    # Run a LONG simulation to see the full decay curve and floor clearly
    long_time = 150.0
    save_every = 100  # more points

    config = Config(;
        sim=Thermalize(), domain=EnergyDomain(), construction=KMS(),
        num_qubits=num_qubits, with_linear_combination=true,
        beta=beta, sigma=sigma, a=beta/30.0, s=0.4,
        num_energy_bits=num_energy_bits, w0=w0, t0=t0,
        num_trotter_steps_per_t0=10,
        mixing_time=long_time, delta=delta,
    )

    println("\n  Running long simulation (T=$long_time, δ=$delta)...")
    flush(stdout)

    result = run_thermalize(jumps, config, hamiltonian;
        rng=MersenneTwister(42), save_every=save_every)

    wall = result.metadata[:wall_time_seconds]
    @printf("  Done in %.1fs, %d data points\n", wall, length(result.trace_distances))

    td = result.trace_distances
    ts = result.time_steps

    # Floor analysis — look at last quarter
    last_quarter = td[div(3*length(td),4):end]
    actual_floor = mean(last_quarter)
    @printf("  Actual floor (last quarter mean): %.6e\n", actual_floor)
    @printf("  Floor std dev: %.6e\n", std(last_quarter))
    @printf("  Final trace distance: %.6e\n", td[end])

    # Fit with different amounts of data (simulate different truncation points)
    # This reveals whether the fit quality depends on how much data we use
    println("\n  === Bi-exp fit sensitivity to truncation ===")
    for frac in [0.25, 0.33, 0.5, 0.67, 0.8, 1.0]
        # Simulate truncation by only using first frac of data
        n_pts = max(10, round(Int, frac * length(td)))
        # Create a "truncated" result by modifying the data we pass
        # We'll just fit directly on the data
        t_trunc = ts[1:n_pts]
        d_trunc = td[1:n_pts]
        T_trunc = t_trunc[end]

        est = estimate_mixing_time(result;
            model=:biexp, skip_initial=0.2,
            target_epsilon=target, extrapolate=true)

        # But to properly test truncation, we need to make a new result...
        # Instead, let's just report the full fit and then also see what happens
        # when we fit shorter data. Use the fit_biexponential_decay directly.
        skip_n = max(1, round(Int, 0.2 * n_pts))
        t_fit = t_trunc[skip_n:end]
        d_fit = d_trunc[skip_n:end]

        try
            bifit = fit_biexponential_decay(t_fit, d_fit)
            t_extrap = nothing
            if bifit.offset < target
                # Simple extrapolation: find where A_s*exp(-g_s*t) + A_f*exp(-g_f*t) + C = target
                # At large t, dominated by slow component: t ≈ -ln((ε-C)/A_s)/g_s
                t_approx = -log((target - bifit.offset) / bifit.amplitude) / bifit.gap
                t_extrap = t_approx
            end
            @printf("    T_trunc=%5.1f: A_f=%.3e g_f=%.3f  A_s=%.3e g_s=%.3f  C=%.3e",
                T_trunc, bifit.amplitude_fast, bifit.gap_fast,
                bifit.amplitude, bifit.gap, bifit.offset)
            if t_extrap !== nothing
                @printf("  t_cross≈%.1f", t_extrap)
            else
                print("  C≥ε!")
            end
            println()
        catch e
            @printf("    T_trunc=%5.1f: FIT FAILED (%s)\n", T_trunc, e)
        end
    end

    # Full fit with standard parameters
    est_full = estimate_mixing_time(result;
        model=:biexp, skip_initial=0.2,
        target_epsilon=target, extrapolate=true)
    bifit = est_full.biexp_fit_result
    @printf("\n  Full fit (T=%.0f, skip=0.2):\n", long_time)
    @printf("    A_fast=%.4e, gap_fast=%.4f\n", bifit.amplitude_fast, bifit.gap_fast)
    @printf("    A_slow=%.4e, gap_slow=%.4f\n", bifit.amplitude, bifit.gap)
    @printf("    Offset C=%.4e\n", bifit.offset)
    @printf("    R² = %.8f\n", est_full.r_squared)
    if est_full.mixing_time_extrapolated !== nothing && !isnan(est_full.mixing_time_extrapolated)
        @printf("    Extrapolated t = %.2f (steps = %d)\n",
            est_full.mixing_time_extrapolated, ceil(Int, est_full.mixing_time_extrapolated / delta))
    end

    # Find actual crossing time
    crossing_idx = findfirst(d -> d ≤ target, td)
    if crossing_idx !== nothing
        @printf("\n  ACTUAL crossing at t = %.4f (step %d)\n",
            ts[crossing_idx], round(Int, ts[crossing_idx] / delta))
    else
        @printf("\n  Target NOT reached. Min = %.4e at t = %.2f\n",
            minimum(td), ts[argmin(td)])
    end

    # Print trace distance at key time points
    println("\n  Decay curve samples:")
    sample_times = [5, 10, 20, 30, 40, 50, 60, 70, 80, 100, 120, 150]
    for st in sample_times
        idx = findfirst(t -> t >= st, ts)
        if idx !== nothing
            @printf("    t=%3d: d=%.4e\n", st, td[idx])
        end
    end

    # Residual analysis: how well does bi-exp fit the late-time data?
    println("\n  === Late-time residuals (bi-exp model vs data) ===")
    skip_n = max(1, round(Int, 0.2 * length(td)))
    t_fit = ts[skip_n:end]
    d_fit = td[skip_n:end]
    model_vals = bifit.amplitude_fast .* exp.(-bifit.gap_fast .* t_fit) .+
                 bifit.amplitude .* exp.(-bifit.gap .* t_fit) .+
                 bifit.offset
    residuals = d_fit .- model_vals

    # Split into time segments
    n_seg = length(t_fit)
    segments = [
        ("early (20-40%)", div(n_seg,5):div(2*n_seg,5)),
        ("mid   (40-60%)", div(2*n_seg,5):div(3*n_seg,5)),
        ("late  (60-80%)", div(3*n_seg,5):div(4*n_seg,5)),
        ("final (80-100%)", div(4*n_seg,5):n_seg),
    ]
    for (name, rng) in segments
        seg_res = residuals[rng]
        @printf("    %s: mean=%.3e, std=%.3e, max|r|=%.3e\n",
            name, mean(seg_res), std(seg_res), maximum(abs.(seg_res)))
    end
end
