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

The channel eigenvalue mapping error is O(delta^2), so the convergence order
should be >= 1.5 (with margin for sub-leading terms).

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
    config_liouv = make_liouv_config(domain; with_coherent=true)
    gap_L = krylov_spectral_gap(config_liouv, hamiltonian, jumps;
        trotter=trotter, krylovdim=krylovdim, howmany=4, tol=tol).spectral_gap

    rows = NamedTuple{(:delta, :gap_from_E, :error), Tuple{Float64, Float64, Float64}}[]
    for delta in deltas
        config_therm = make_thermalize_config(domain; with_coherent=true, delta=delta)
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

        @testset "EnergyDomain" begin
            comp = compare_krylov_dense(
                make_liouv_config(EnergyDomain(); with_coherent=true),
                TEST_HAM, TEST_JUMPS)
            print_gap_summary("EnergyDomain", "KMS",
                comp.krylov_result.spectral_gap, comp.dense_result.spectral_gap)
            gap_match = isapprox(comp.krylov_result.spectral_gap,
                comp.dense_result.spectral_gap; atol=1e-8)
            if !gap_match
                on_failure_diagnostics(comp.krylov_result, comp.dense_result)
            end
            @test gap_match
        end

        @testset "TimeDomain" begin
            comp = compare_krylov_dense(
                make_liouv_config(TimeDomain(); with_coherent=true),
                TEST_HAM, TEST_JUMPS)
            print_gap_summary("TimeDomain", "KMS",
                comp.krylov_result.spectral_gap, comp.dense_result.spectral_gap)
            gap_match = isapprox(comp.krylov_result.spectral_gap,
                comp.dense_result.spectral_gap; atol=1e-8)
            if !gap_match
                on_failure_diagnostics(comp.krylov_result, comp.dense_result)
            end
            @test gap_match
        end

        @testset "TrotterDomain" begin
            comp = compare_krylov_dense(
                make_liouv_config(TrotterDomain(); with_coherent=true),
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
        end

        @testset "BohrDomain" begin
            comp = compare_krylov_dense(
                make_liouv_config(BohrDomain(); with_coherent=true),
                TEST_HAM, TEST_JUMPS)
            print_gap_summary("BohrDomain", "KMS",
                comp.krylov_result.spectral_gap, comp.dense_result.spectral_gap)
            gap_match = isapprox(comp.krylov_result.spectral_gap,
                comp.dense_result.spectral_gap; atol=1e-8)
            if !gap_match
                on_failure_diagnostics(comp.krylov_result, comp.dense_result)
            end
            @test gap_match
        end

    end  # n=4 KMS

    # ========================================================================
    # XVAL-04: n=4 GNS cross-validation across all 4 domains
    # GNS uses make_liouv_config_gns (with_coherent=false)
    # Tolerance: atol=1e-8
    # ========================================================================
    @testset "n=4 GNS (all domains)" begin

        @testset "EnergyDomain" begin
            comp = compare_krylov_dense(
                make_liouv_config_gns(EnergyDomain()),
                TEST_HAM, TEST_JUMPS)
            print_gap_summary("EnergyDomain", "GNS",
                comp.krylov_result.spectral_gap, comp.dense_result.spectral_gap)
            gap_match = isapprox(comp.krylov_result.spectral_gap,
                comp.dense_result.spectral_gap; atol=1e-8)
            if !gap_match
                on_failure_diagnostics(comp.krylov_result, comp.dense_result)
            end
            @test gap_match
        end

        @testset "TimeDomain" begin
            comp = compare_krylov_dense(
                make_liouv_config_gns(TimeDomain()),
                TEST_HAM, TEST_JUMPS)
            print_gap_summary("TimeDomain", "GNS",
                comp.krylov_result.spectral_gap, comp.dense_result.spectral_gap)
            gap_match = isapprox(comp.krylov_result.spectral_gap,
                comp.dense_result.spectral_gap; atol=1e-8)
            if !gap_match
                on_failure_diagnostics(comp.krylov_result, comp.dense_result)
            end
            @test gap_match
        end

        @testset "TrotterDomain" begin
            comp = compare_krylov_dense(
                make_liouv_config_gns(TrotterDomain()),
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
        end

        @testset "BohrDomain" begin
            comp = compare_krylov_dense(
                make_liouv_config_gns(BohrDomain()),
                TEST_HAM, TEST_JUMPS)
            print_gap_summary("BohrDomain", "GNS",
                comp.krylov_result.spectral_gap, comp.dense_result.spectral_gap)
            gap_match = isapprox(comp.krylov_result.spectral_gap,
                comp.dense_result.spectral_gap; atol=1e-8)
            if !gap_match
                on_failure_diagnostics(comp.krylov_result, comp.dense_result)
            end
            @test gap_match
        end

    end  # n=4 GNS

    # ========================================================================
    # XVAL-03: L-vs-E convergence (KMS only, per locked decision)
    # Tests that channel-to-Lindbladian gap mapping converges with O(delta^2)
    # Deltas: [0.1, 0.01, 0.001]
    # Hard assertion: convergence order >= 1.5 for each consecutive pair
    # ========================================================================
    @testset "L-vs-E convergence (KMS)" begin

        @testset "EnergyDomain" begin
            result = run_le_convergence(EnergyDomain(), TEST_HAM, TEST_JUMPS)
            for (i, order) in enumerate(result.orders)
                @test order >= 1.5
            end
        end

        @testset "TimeDomain" begin
            result = run_le_convergence(TimeDomain(), TEST_HAM, TEST_JUMPS)
            for (i, order) in enumerate(result.orders)
                @test order >= 1.5
            end
        end

        @testset "TrotterDomain" begin
            result = run_le_convergence(TrotterDomain(), TEST_HAM, TEST_TROTTER_JUMPS;
                trotter=TEST_TROTTER)
            for (i, order) in enumerate(result.orders)
                @test order >= 1.5
            end
        end

        @testset "BohrDomain" begin
            result = run_le_convergence(BohrDomain(), TEST_HAM, TEST_JUMPS)
            for (i, order) in enumerate(result.orders)
                @test order >= 1.5
            end
        end

    end  # L-vs-E convergence

    # n=6 block will be added by Plan 30-02

end  # @testset "Krylov Cross-Validation"
