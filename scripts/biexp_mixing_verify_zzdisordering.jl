#!/usr/bin/env julia
#
# Bi-exponential mixing time verification with Z+ZZ disordering.
#
# Generates fresh Hamiltonians with disordering_terms=[[Z], [Z,Z]] to break
# bipartite symmetry (eliminates Bohr frequency collisions in even-n chains).
# Tests both EnergyDomain and TrotterDomain for n = 3, 4, 5.
#
# Key insight: the bi-exp fit needs data well past the crossing point to properly
# separate the exponential decay from the floor offset. Use generous extrap times
# (~3× the expected mixing time).
#
# Usage:  OPENBLAS_NUM_THREADS=4 julia -t4 --project scripts/biexp_mixing_verify_zzdisordering.jl

using QuantumFurnace
using LinearAlgebra
using Random
using Printf

BLAS.set_num_threads(4)
println("BLAS threads: $(BLAS.get_num_threads()), Julia threads: $(Threads.nthreads())")

# ── Parameters ────────────────────────────────────────────────────────────────
beta       = 10.0
target     = 1e-4
coeffs     = fill(1.0, 3)   # XX + YY + ZZ Heisenberg
qubit_range = [3, 4, 5]
delta      = 0.0005          # Same for all n; floor ~ 3-6e-5, well below target

# Extrapolation time must extend well past the crossing (~3× expected mixing time).
# Crossing times (measured): n=3→~35, n=4→~51, n=5→~54.
extrap_time_map = Dict(3 => 100.0, 4 => 150.0, 5 => 150.0)

# ══════════════════════════════════════════════════════════════════════════════
# MAIN LOOP
# ══════════════════════════════════════════════════════════════════════════════

results_summary = []

for num_qubits in qubit_range
    println("\n" * "▓"^70)
    @printf("  n = %d qubits  (Z+ZZ disordering, δ = %.4f)\n", num_qubits, delta)
    println("▓"^70)

    # ── 1. Generate Hamiltonian with Z+ZZ disordering ─────────────────────
    println("\n  Generating Hamiltonian (batch_size=200)...")
    flush(stdout)

    raw = find_ideal_heisenberg(num_qubits, coeffs;
        batch_size     = 200,
        periodic       = true,
        disordering_terms = [[Z], [Z, Z]],
    )
    hamiltonian = HamHam(raw, beta)

    @printf("  dim = %d, nu_min = %.6e\n", 2^num_qubits, raw.nu_min)
    println("  disordering_terms: Z (single-site) + ZZ (two-site)")

    # Verify no Bohr frequency collisions
    bohr = hamiltonian.bohr_freqs
    dim = 2^num_qubits
    bohr_upper = [bohr[i,j] for i in 1:dim for j in (i+1):dim]
    sort!(bohr_upper)
    min_bohr_gap = minimum(diff(bohr_upper))
    @printf("  Min Bohr freq gap: %.4e  (collisions broken? %s)\n",
        min_bohr_gap, min_bohr_gap > 1e-10 ? "YES" : "NO")

    # ── 2. Build jump operators ───────────────────────────────────────────
    jump_paulis = [[X], [Y], [Z]]
    n_jumps     = length(jump_paulis) * num_qubits
    norm_factor = sqrt(n_jumps)

    jumps_energy = JumpOp[]
    for pauli in jump_paulis
        for site in 1:num_qubits
            op       = Matrix(pad_term(pauli, num_qubits, site)) ./ norm_factor
            op_eigen = hamiltonian.eigvecs' * op * hamiltonian.eigvecs
            push!(jumps_energy, JumpOp(op, op_eigen, op == transpose(op), op == op'))
        end
    end

    # ── 3. Shared config parameters ───────────────────────────────────────
    sigma  = 1.0 / beta
    w0     = 0.05
    num_energy_bits = 12
    t0     = 2π / (2^num_energy_bits * w0)
    num_trotter_steps_per_t0 = 10

    extrap_time = extrap_time_map[num_qubits]
    extrap_save_every = 100

    for domain_name in ["EnergyDomain", "TrotterDomain"]
        println("\n" * "="^70)
        @printf("  %s  (n=%d, δ=%.4f)\n", domain_name, num_qubits, delta)
        println("="^70)

        is_trotter = (domain_name == "TrotterDomain")
        domain = is_trotter ? TrotterDomain() : EnergyDomain()

        # Build Trotter if needed
        trotter = nothing
        jumps_domain = jumps_energy
        if is_trotter
            trotter = TrottTrott(hamiltonian, t0, num_trotter_steps_per_t0)
            trotter_error = compute_trotter_error(hamiltonian, trotter, 2^num_energy_bits * t0 / 2)
            @printf("  Trotter error over full T: %.2e\n", trotter_error)

            # Rebuild jumps in Trotter eigenbasis
            jumps_domain = JumpOp[]
            for pauli in jump_paulis
                for site in 1:num_qubits
                    op       = Matrix(pad_term(pauli, num_qubits, site)) ./ norm_factor
                    op_eigen = trotter.eigvecs' * op * trotter.eigvecs
                    push!(jumps_domain, JumpOp(op, op_eigen, op == transpose(op), op == op'))
                end
            end
        end

        # ── Phase 1: Truncated run + bi-exp extrapolation ─────────────────
        println("\n  PHASE 1: Truncated extrapolation run")
        extrap_steps = round(Int, extrap_time / delta)
        @printf("    T = %.1f → %d steps, save_every = %d\n", extrap_time, extrap_steps, extrap_save_every)
        println("    Running...")
        flush(stdout)

        config_extrap = Config(;
            sim                       = Thermalize(),
            domain                    = domain,
            construction              = KMS(),
            num_qubits                = num_qubits,
            with_linear_combination   = true,
            beta                      = beta,
            sigma                     = sigma,
            a                         = beta / 30.0,
            b                         = 0.4,
            num_energy_bits           = num_energy_bits,
            w0                        = w0,
            t0                        = t0,
            num_trotter_steps_per_t0  = num_trotter_steps_per_t0,
            mixing_time               = extrap_time,
            delta                     = delta,
        )

        result_extrap = run_thermalize(jumps_domain, config_extrap, hamiltonian, trotter;
            rng        = MersenneTwister(42),
            save_every = extrap_save_every,
        )

        wall1 = result_extrap.metadata[:wall_time_seconds]
        @printf("    Done in %.1fs\n", wall1)
        @printf("    Final trace distance: %.4e\n", result_extrap.trace_distances[end])

        # Bi-exponential fit
        est_bi = estimate_mixing_time(result_extrap;
            model          = :biexp,
            skip_initial   = 0.2,
            target_epsilon = target,
            extrapolate    = true,
        )

        bifit = est_bi.biexp_fit_result
        println("\n    --- Bi-exponential fit ---")
        @printf("    d(t) = %.4e × exp(-%.4f t) + %.4e × exp(-%.4f t) + %.4e\n",
            bifit.amplitude_fast, bifit.gap_fast, bifit.amplitude, bifit.gap, bifit.offset)
        @printf("    R²       = %.6f\n", est_bi.r_squared)
        @printf("    Slow gap = %.4f,  Fast gap = %.4f\n", bifit.gap, bifit.gap_fast)
        @printf("    Offset C = %.4e  (vs target ε = %.4e)\n", bifit.offset, target)
        @printf("    Converged = %s\n", est_bi.converged)

        if est_bi.mixing_time_extrapolated === nothing || isnan(est_bi.mixing_time_extrapolated)
            println("\n    ERROR: Bi-exponential extrapolation failed!")
            @printf("    Offset C = %.4e ≥ target ε = %.4e ?\n", bifit.offset, target)
            push!(results_summary, (n=num_qubits, domain=domain_name, status="FAIL (extrapolation)", error_pct=NaN))
            continue
        end

        predicted_time  = est_bi.mixing_time_extrapolated
        predicted_steps = ceil(Int, predicted_time / delta)
        @printf("\n    Extrapolated crossing: t = %.5f\n", predicted_time)
        @printf("    Predicted delta steps: %d\n", predicted_steps)

        # Single-exp comparison
        est_single = estimate_mixing_time(result_extrap;
            model          = :single,
            skip_initial   = 0.2,
            target_epsilon = target,
            extrapolate    = true,
        )
        if est_single.mixing_time_extrapolated !== nothing && !isnan(est_single.mixing_time_extrapolated)
            single_steps = ceil(Int, est_single.mixing_time_extrapolated / delta)
            @printf("    (Single-exp: %d steps, offset C = %.4e)\n", single_steps, est_single.offset)
        end

        # ── Phase 2: Verification ─────────────────────────────────────────
        println("\n  PHASE 2: Verification run")
        verify_mixing_time = predicted_steps * delta
        verify_save_every  = max(1, div(predicted_steps, 1000))

        @printf("    T = %.5f → %d steps, save_every = %d\n",
            verify_mixing_time, predicted_steps, verify_save_every)
        println("    Running...")
        flush(stdout)

        config_verify = Config(;
            sim                       = Thermalize(),
            domain                    = domain,
            construction              = KMS(),
            num_qubits                = num_qubits,
            with_linear_combination   = true,
            beta                      = beta,
            sigma                     = sigma,
            a                         = beta / 30.0,
            b                         = 0.4,
            num_energy_bits           = num_energy_bits,
            w0                        = w0,
            t0                        = t0,
            num_trotter_steps_per_t0  = num_trotter_steps_per_t0,
            mixing_time               = verify_mixing_time,
            delta                     = delta,
        )

        result_verify = run_thermalize(jumps_domain, config_verify, hamiltonian, trotter;
            rng        = MersenneTwister(42),
            save_every = verify_save_every,
        )

        wall2 = result_verify.metadata[:wall_time_seconds]
        final_dist = result_verify.trace_distances[end]
        min_dist   = minimum(result_verify.trace_distances)

        @printf("    Done in %.1fs\n", wall2)

        # ── Results ───────────────────────────────────────────────────────
        crossing_idx = findfirst(d -> d ≤ target, result_verify.trace_distances)

        println()
        if crossing_idx !== nothing
            crossing_time = result_verify.time_steps[crossing_idx]
            crossing_step = round(Int, crossing_time / delta)
            pct_err = abs(crossing_step - predicted_steps) / crossing_step * 100

            @printf("    ✓ Target REACHED at step %d  (t = %.5f)\n", crossing_step, crossing_time)
            @printf("      Bi-exp predicted %d steps → error: %.2f%%\n", predicted_steps, pct_err)

            if final_dist ≤ target
                @printf("      Final distance %.4e ≤ %.4e  ✓  PASS\n", final_dist, target)
                push!(results_summary, (n=num_qubits, domain=domain_name, status="PASS", error_pct=pct_err))
            else
                @printf("      Final distance %.4e > %.4e (rose above target by end)\n", final_dist, target)
                push!(results_summary, (n=num_qubits, domain=domain_name, status="PASS (transient)", error_pct=pct_err))
            end
        else
            @printf("    ✗ Target NOT reached. Final: %.4e, min: %.4e\n", final_dist, min_dist)
            push!(results_summary, (n=num_qubits, domain=domain_name, status="FAIL", error_pct=NaN))
        end

        @printf("    Wall time: %.1fs (extrap) + %.1fs (verify) = %.1fs\n", wall1, wall2, wall1 + wall2)
    end
end

# ══════════════════════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════════════════════
println("\n\n" * "▓"^70)
println("  SUMMARY  (Z+ZZ disordering, δ = $delta, target ε = $target)")
println("▓"^70)
println()
@printf("  %-4s  %-15s  %-18s  %s\n", "n", "Domain", "Status", "Error %")
println("  " * "-"^55)
for r in results_summary
    err_str = isnan(r.error_pct) ? "N/A" : @sprintf("%.2f%%", r.error_pct)
    @printf("  %-4d  %-15s  %-18s  %s\n", r.n, r.domain, r.status, err_str)
end
println()
