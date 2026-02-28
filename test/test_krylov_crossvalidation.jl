using Test
using LinearAlgebra
using Printf
using QuantumFurnace

# test_helpers.jl is already included by runtests.jl

# ============================================================================
# Cross-validation: Krylov spectral gap vs dense eigen() reference
# Phase 30: Establishes quantitative trust for n>6 production use
#
# XVAL-01: n=4 KMS cross-validation across all 4 domains (atol=1e-8)
# XVAL-03: L-vs-E convergence testing (O(delta^2) order >= 1.5)
# XVAL-04: n=4 GNS cross-validation across all 4 domains (atol=1e-8)
# ============================================================================

# ---------------------------------------------------------------------------
# Helper: compare Krylov and dense eigsolve on the same system
# ---------------------------------------------------------------------------

"""
    compare_krylov_dense(config_liouv, hamiltonian, jumps; kwargs...) -> NamedTuple

Run both Krylov matrix-free eigsolve and dense eigen() on the same Lindbladian,
returning structured comparison results.

Returns `(; krylov_result, dense_result, L_dense)`.
"""
function compare_krylov_dense(config_liouv, hamiltonian, jumps;
    trotter=nothing,
    krylovdim=30,
    howmany=4,
    tol=1e-10,
    n_dense_modes=4,
)
    # Dense reference: build full Lindbladian and extract leading eigendata
    L_dense = construct_lindbladian(jumps, config_liouv, hamiltonian; trotter=trotter)
    dense_result = extract_leading_eigendata(L_dense; n_modes=n_dense_modes)

    # Krylov result via matrix-free eigsolve
    krylov_result = krylov_spectral_gap(config_liouv, hamiltonian, jumps;
        trotter=trotter, krylovdim=krylovdim, howmany=howmany, tol=tol)

    return (; krylov_result, dense_result, L_dense)
end

# ---------------------------------------------------------------------------
# Diagnostic printing helpers
# ---------------------------------------------------------------------------

"""
    print_gap_summary(domain_name, balance_name, krylov_gap, dense_gap)

Print one-line summary: domain, balance, gap values, absolute and relative error.
Always printed on both success and failure for a record of achieved accuracy.
"""
function print_gap_summary(domain_name, balance_name, krylov_gap, dense_gap)
    err = abs(krylov_gap - dense_gap)
    rel_err = err / max(dense_gap, 1e-30)
    @printf("  %-15s %-5s | gap_krylov=%.8e  gap_dense=%.8e  err=%.2e  rel=%.2e\n",
            domain_name, balance_name, krylov_gap, dense_gap, err, rel_err)
end

"""
    print_eigenvalue_table(krylov_eigs, dense_eigs; top_k=6)

Print formatted comparison table of top-k eigenvalue real parts from Krylov vs dense.
"""
function print_eigenvalue_table(krylov_eigs, dense_eigs; top_k=6)
    n = min(top_k, length(krylov_eigs), length(dense_eigs))
    println("  Top-$n eigenvalue comparison:")
    println("  idx | Re(krylov)        | Re(dense)         | abs_err      | rel_err")
    println("  " * "-"^75)
    for k in 1:n
        kr = real(krylov_eigs[k])
        dr = real(dense_eigs[k])
        ae = abs(kr - dr)
        re = ae / max(abs(dr), 1e-30)
        @printf("  %3d | %+.10e | %+.10e | %.2e | %.2e\n", k, kr, dr, ae, re)
    end
end

"""
    on_failure_diagnostics(krylov_result, dense_result)

Print diagnostic information on test failure: top-k eigenvalue table and
Krylov convergence metadata.
"""
function on_failure_diagnostics(krylov_result, dense_result)
    println("  FAILURE DIAGNOSTICS:")
    print_eigenvalue_table(krylov_result.eigenvalues, dense_result.eigenvalues; top_k=6)
    @printf("  Krylov converged: %d, matvec_count: %d, restarts: %d\n",
            krylov_result.converged, krylov_result.matvec_count, krylov_result.num_restarts)
end

# ---------------------------------------------------------------------------
# L-vs-E convergence analysis helper
# ---------------------------------------------------------------------------

"""
    run_le_convergence(domain, hamiltonian, jumps; kwargs...) -> NamedTuple

Compute L-vs-E convergence: compare Lindbladian spectral gap (delta-independent)
against channel-derived gap at multiple delta values.

The faithful Chen channel has mu = exp(delta*lambda_L) + O(delta^2), so the
first-order conversion (mu-1)/delta introduces O(delta) error. Expected
convergence order is ~1.0; threshold is >= 0.9 with margin for sub-leading terms.

Returns `(; gap_L, rows, orders)` where rows is a vector of
`(; delta, gap_from_E, error)` named tuples and orders is a vector of
convergence orders from consecutive error pairs.
"""
function run_le_convergence(domain, hamiltonian, jumps;
    deltas=[0.1, 0.01, 0.001],
    trotter=nothing,
    krylovdim=30,
    tol=1e-10,
)
    # Lindbladian reference gap (delta-independent)
    config_liouv = make_config(Lindbladian(),domain; construction=KMS())
    gap_L = krylov_spectral_gap(config_liouv, hamiltonian, jumps;
        trotter=trotter, krylovdim=krylovdim, howmany=4, tol=tol).spectral_gap

    rows = NamedTuple{(:delta, :gap_from_E, :error), Tuple{Float64, Float64, Float64}}[]
    for delta in deltas
        config_therm = make_config(Thermalize(),domain; construction=KMS(), delta=delta)
        gap_from_E = krylov_spectral_gap(config_therm, hamiltonian, jumps;
            trotter=trotter, krylovdim=krylovdim, howmany=4, tol=tol).spectral_gap
        error_val = abs(gap_L - gap_from_E)
        push!(rows, (; delta, gap_from_E, error=error_val))
    end

    # Compute convergence orders from consecutive error pairs
    orders = Float64[]
    for i in 2:length(rows)
        if rows[i-1].error > 0 && rows[i].error > 0
            ratio = log(rows[i-1].error / rows[i].error) / log(rows[i-1].delta / rows[i].delta)
            push!(orders, ratio)
        end
    end

    # Print formatted convergence table
    domain_name = replace(string(typeof(domain)), "Domain" => "")
    println("  L-vs-E convergence ($domain_name, KMS):")
    @printf("  %-10s | %-18s | %-18s | %-12s | %s\n",
            "delta", "gap_L", "gap_from_E", "error", "order")
    println("  " * "-"^75)
    for (idx, row) in enumerate(rows)
        order_str = idx == 1 ? "  --" : @sprintf("%.2f", orders[idx-1])
        @printf("  %.2e | %.12e | %.12e | %.2e | %s\n",
                row.delta, gap_L, row.gap_from_E, row.error, order_str)
    end

    return (; gap_L, rows, orders)
end

# ============================================================================
# Test suites
# ============================================================================

@testset "Krylov Cross-Validation" begin

    # ========================================================================
    # XVAL-01: n=4 KMS cross-validation across all 4 domains
    # Tolerance: atol=1e-8 (KrylovKit tol=1e-10 provides margin)
    # ========================================================================
    @testset "n=4 KMS (all domains)" begin

        # Threshold rationale (atol=1e-8): KrylovKit tol=1e-10 bounds eigenvalue error.
        # atol=1e-8 gives 100x margin for iterative convergence variability on n=4 (DIM=16).
        @testset "EnergyDomain" begin
            comp = compare_krylov_dense(
                make_config(Lindbladian(),EnergyDomain(); construction=KMS()),
                TEST_HAM, TEST_JUMPS)
            print_gap_summary("EnergyDomain", "KMS",
                comp.krylov_result.spectral_gap, comp.dense_result.spectral_gap)
            gap_match = isapprox(comp.krylov_result.spectral_gap,
                comp.dense_result.spectral_gap; atol=1e-8)
            if !gap_match
                on_failure_diagnostics(comp.krylov_result, comp.dense_result)
            end
            @test gap_match
            @info "XVAL-01 gap (EnergyDomain KMS)" error=abs(comp.krylov_result.spectral_gap - comp.dense_result.spectral_gap) atol=1e-8
        end

        @testset "TimeDomain" begin
            comp = compare_krylov_dense(
                make_config(Lindbladian(),TimeDomain(); construction=KMS()),
                TEST_HAM, TEST_JUMPS)
            print_gap_summary("TimeDomain", "KMS",
                comp.krylov_result.spectral_gap, comp.dense_result.spectral_gap)
            gap_match = isapprox(comp.krylov_result.spectral_gap,
                comp.dense_result.spectral_gap; atol=1e-8)
            if !gap_match
                on_failure_diagnostics(comp.krylov_result, comp.dense_result)
            end
            @test gap_match
            @info "XVAL-01 gap (TimeDomain KMS)" error=abs(comp.krylov_result.spectral_gap - comp.dense_result.spectral_gap) atol=1e-8
        end

        @testset "TrotterDomain" begin
            comp = compare_krylov_dense(
                make_config(Lindbladian(),TrotterDomain(); construction=KMS()),
                TEST_HAM, TEST_TROTTER_JUMPS;
                trotter=TEST_TROTTER)
            print_gap_summary("TrotterDomain", "KMS",
                comp.krylov_result.spectral_gap, comp.dense_result.spectral_gap)
            gap_match = isapprox(comp.krylov_result.spectral_gap,
                comp.dense_result.spectral_gap; atol=1e-8)
            if !gap_match
                on_failure_diagnostics(comp.krylov_result, comp.dense_result)
            end
            @test gap_match
            @info "XVAL-01 gap (TrotterDomain KMS)" error=abs(comp.krylov_result.spectral_gap - comp.dense_result.spectral_gap) atol=1e-8
        end

        @testset "BohrDomain" begin
            comp = compare_krylov_dense(
                make_config(Lindbladian(),BohrDomain(); construction=KMS()),
                TEST_HAM, TEST_JUMPS)
            print_gap_summary("BohrDomain", "KMS",
                comp.krylov_result.spectral_gap, comp.dense_result.spectral_gap)
            gap_match = isapprox(comp.krylov_result.spectral_gap,
                comp.dense_result.spectral_gap; atol=1e-8)
            if !gap_match
                on_failure_diagnostics(comp.krylov_result, comp.dense_result)
            end
            @test gap_match
            @info "XVAL-01 gap (BohrDomain KMS)" error=abs(comp.krylov_result.spectral_gap - comp.dense_result.spectral_gap) atol=1e-8
        end

    end  # n=4 KMS

    # ========================================================================
    # XVAL-04: n=4 GNS cross-validation across all 4 domains
    # GNS uses make_config(Lindbladian(), ...; construction=GNS()) (with_coherent=false)
    # Tolerance: atol=1e-8
    # ========================================================================
    @testset "n=4 GNS (all domains)" begin

        # Threshold rationale (atol=1e-8): same as XVAL-01 KMS -- KrylovKit tol=1e-10, 100x margin.
        @testset "EnergyDomain" begin
            comp = compare_krylov_dense(
                make_config(Lindbladian(), EnergyDomain(); construction=GNS()),
                TEST_HAM, TEST_JUMPS)
            print_gap_summary("EnergyDomain", "GNS",
                comp.krylov_result.spectral_gap, comp.dense_result.spectral_gap)
            gap_match = isapprox(comp.krylov_result.spectral_gap,
                comp.dense_result.spectral_gap; atol=1e-8)
            if !gap_match
                on_failure_diagnostics(comp.krylov_result, comp.dense_result)
            end
            @test gap_match
            @info "XVAL-04 gap (EnergyDomain GNS)" error=abs(comp.krylov_result.spectral_gap - comp.dense_result.spectral_gap) atol=1e-8
        end

        @testset "TimeDomain" begin
            comp = compare_krylov_dense(
                make_config(Lindbladian(), TimeDomain(); construction=GNS()),
                TEST_HAM, TEST_JUMPS)
            print_gap_summary("TimeDomain", "GNS",
                comp.krylov_result.spectral_gap, comp.dense_result.spectral_gap)
            gap_match = isapprox(comp.krylov_result.spectral_gap,
                comp.dense_result.spectral_gap; atol=1e-8)
            if !gap_match
                on_failure_diagnostics(comp.krylov_result, comp.dense_result)
            end
            @test gap_match
            @info "XVAL-04 gap (TimeDomain GNS)" error=abs(comp.krylov_result.spectral_gap - comp.dense_result.spectral_gap) atol=1e-8
        end

        @testset "TrotterDomain" begin
            comp = compare_krylov_dense(
                make_config(Lindbladian(), TrotterDomain(); construction=GNS()),
                TEST_HAM, TEST_TROTTER_JUMPS;
                trotter=TEST_TROTTER)
            print_gap_summary("TrotterDomain", "GNS",
                comp.krylov_result.spectral_gap, comp.dense_result.spectral_gap)
            gap_match = isapprox(comp.krylov_result.spectral_gap,
                comp.dense_result.spectral_gap; atol=1e-8)
            if !gap_match
                on_failure_diagnostics(comp.krylov_result, comp.dense_result)
            end
            @test gap_match
            @info "XVAL-04 gap (TrotterDomain GNS)" error=abs(comp.krylov_result.spectral_gap - comp.dense_result.spectral_gap) atol=1e-8
        end

        @testset "BohrDomain" begin
            comp = compare_krylov_dense(
                make_config(Lindbladian(), BohrDomain(); construction=GNS()),
                TEST_HAM, TEST_JUMPS)
            print_gap_summary("BohrDomain", "GNS",
                comp.krylov_result.spectral_gap, comp.dense_result.spectral_gap)
            gap_match = isapprox(comp.krylov_result.spectral_gap,
                comp.dense_result.spectral_gap; atol=1e-8)
            if !gap_match
                on_failure_diagnostics(comp.krylov_result, comp.dense_result)
            end
            @test gap_match
            @info "XVAL-04 gap (BohrDomain GNS)" error=abs(comp.krylov_result.spectral_gap - comp.dense_result.spectral_gap) atol=1e-8
        end

    end  # n=4 GNS

    # ========================================================================
    # XVAL-03: L-vs-E convergence (KMS only, per locked decision)
    # Tests that channel-to-Lindbladian gap mapping converges with O(delta)
    # The faithful Chen channel gives mu = exp(delta*lambda_L) + O(delta^2),
    # so (mu-1)/delta has first-order error. Deltas: [0.1, 0.01, 0.001]
    # Hard assertion: convergence order >= 0.9 for each consecutive pair
    # ========================================================================
    # Threshold rationale (order >= 0.9): faithful Chen channel gives mu = exp(delta*lambda_L) + O(delta^2),
    # so first-order conversion (mu-1)/delta introduces O(delta) error. Expected convergence order ~1.0.
    # Threshold 0.9 gives margin for sub-leading terms and numerical noise.
    @testset "L-vs-E convergence (KMS)" begin

        @testset "EnergyDomain" begin
            result = run_le_convergence(EnergyDomain(), TEST_HAM, TEST_JUMPS)
            for (i, order) in enumerate(result.orders)
                @test order >= 0.9
            end
            @info "XVAL-03 convergence (EnergyDomain)" orders=result.orders threshold=0.9
        end

        @testset "TimeDomain" begin
            result = run_le_convergence(TimeDomain(), TEST_HAM, TEST_JUMPS)
            for (i, order) in enumerate(result.orders)
                @test order >= 0.9
            end
            @info "XVAL-03 convergence (TimeDomain)" orders=result.orders threshold=0.9
        end

        @testset "TrotterDomain" begin
            result = run_le_convergence(TrotterDomain(), TEST_HAM, TEST_TROTTER_JUMPS;
                trotter=TEST_TROTTER)
            for (i, order) in enumerate(result.orders)
                @test order >= 0.9
            end
            @info "XVAL-03 convergence (TrotterDomain)" orders=result.orders threshold=0.9
        end

        @testset "BohrDomain" begin
            result = run_le_convergence(BohrDomain(), TEST_HAM, TEST_JUMPS)
            for (i, order) in enumerate(result.orders)
                @test order >= 0.9
            end
            @info "XVAL-03 convergence (BohrDomain)" orders=result.orders threshold=0.9
        end

    end  # L-vs-E convergence

    # ========================================================================
    # XVAL-02: n=6 KMS cross-validation across all 4 domains
    # Gated behind QUANTUMFURNACE_FULL_TESTS=true (dense 4096x4096 eigen)
    # Tolerance: atol=1e-6 (looser than n=4's 1e-8 due to larger Krylov
    # subspace approximation error at dim^2=4096)
    # ========================================================================
    if get(ENV, "QUANTUMFURNACE_FULL_TESTS", "false") == "true"
        @testset "n=6 KMS (all domains)" begin

            # Construct n=6 test system once for non-Trotter domains
            n6_sys = make_test_system(; num_qubits=6,)
            n6_ham = n6_sys.hamiltonian

            # Construct n=6 Trotter system (separate eigenbasis for TrotterDomain)
            n6_trotter = TrottTrott(n6_ham, T0, NUM_TROTTER_STEPS_PER_T0)
            n6_trotter_sys = make_test_system(; num_qubits=6, trotter=n6_trotter)

            # Threshold rationale (atol=1e-6): n=6 system (DIM=64, DIM^2=4096) has larger Krylov
            # subspace approximation error than n=4. KrylovKit tol=1e-10 but larger system means
            # more FP accumulation. atol=1e-6 is 100x looser than n=4's 1e-8.
            @testset "EnergyDomain" begin
                config = make_config(Lindbladian(), EnergyDomain(); num_qubits=6, construction=KMS())
                comp = compare_krylov_dense(config, n6_ham, n6_sys.jumps)
                print_gap_summary("EnergyDomain", "KMS",
                    comp.krylov_result.spectral_gap, comp.dense_result.spectral_gap)
                gap_match = isapprox(comp.krylov_result.spectral_gap,
                    comp.dense_result.spectral_gap; atol=1e-6)
                if !gap_match
                    on_failure_diagnostics(comp.krylov_result, comp.dense_result)
                end
                @test gap_match
                @info "XVAL-02 gap (EnergyDomain KMS, n=6)" error=abs(comp.krylov_result.spectral_gap - comp.dense_result.spectral_gap) atol=1e-6
            end

            @testset "TimeDomain" begin
                config = make_config(Lindbladian(), TimeDomain(); num_qubits=6, construction=KMS())
                comp = compare_krylov_dense(config, n6_ham, n6_sys.jumps)
                print_gap_summary("TimeDomain", "KMS",
                    comp.krylov_result.spectral_gap, comp.dense_result.spectral_gap)
                gap_match = isapprox(comp.krylov_result.spectral_gap,
                    comp.dense_result.spectral_gap; atol=1e-6)
                if !gap_match
                    on_failure_diagnostics(comp.krylov_result, comp.dense_result)
                end
                @test gap_match
                @info "XVAL-02 gap (TimeDomain KMS, n=6)" error=abs(comp.krylov_result.spectral_gap - comp.dense_result.spectral_gap) atol=1e-6
            end

            @testset "TrotterDomain" begin
                config = make_config(Lindbladian(), TrotterDomain(); num_qubits=6, construction=KMS())
                comp = compare_krylov_dense(config, n6_ham, n6_trotter_sys.jumps;
                    trotter=n6_trotter)
                print_gap_summary("TrotterDomain", "KMS",
                    comp.krylov_result.spectral_gap, comp.dense_result.spectral_gap)
                gap_match = isapprox(comp.krylov_result.spectral_gap,
                    comp.dense_result.spectral_gap; atol=1e-6)
                if !gap_match
                    on_failure_diagnostics(comp.krylov_result, comp.dense_result)
                end
                @test gap_match
                @info "XVAL-02 gap (TrotterDomain KMS, n=6)" error=abs(comp.krylov_result.spectral_gap - comp.dense_result.spectral_gap) atol=1e-6
            end

            @testset "BohrDomain" begin
                config = make_config(Lindbladian(), BohrDomain(); num_qubits=6, construction=KMS())
                comp = compare_krylov_dense(config, n6_ham, n6_sys.jumps)
                print_gap_summary("BohrDomain", "KMS",
                    comp.krylov_result.spectral_gap, comp.dense_result.spectral_gap)
                gap_match = isapprox(comp.krylov_result.spectral_gap,
                    comp.dense_result.spectral_gap; atol=1e-6)
                if !gap_match
                    on_failure_diagnostics(comp.krylov_result, comp.dense_result)
                end
                @test gap_match
                @info "XVAL-02 gap (BohrDomain KMS, n=6)" error=abs(comp.krylov_result.spectral_gap - comp.dense_result.spectral_gap) atol=1e-6
            end

        end  # n=6 KMS
    else
        @info "Skipping n=6 cross-validation (set QUANTUMFURNACE_FULL_TESTS=true to run)"
    end

end  # @testset "Krylov Cross-Validation"
