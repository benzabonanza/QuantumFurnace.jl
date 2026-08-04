"""
Allocation regression tests for optimized hot paths.

Verifies that allocation-reducing optimizations from Phase 11 (Plans 01 and 02)
remain effective. Each test calls the target function once for JIT warmup, then
measures allocations with @allocated. Thresholds are set to catch the eliminated
allocation patterns while allowing expected allocations (return values, scratch
buffers, and broadcasting overhead from closure-based element-wise operations).

Eliminated patterns that these tests guard against:
- B_bohr: per-frequency spzeros + sparse-dense multiply (O(num_freqs * dim^2))
- B_time/B_trotter: per-iteration Diagonal wrapper construction
- retained channel step: per-step R construction and eigendecomposition
"""

using QuantumFurnace: B_bohr, B_time, B_trotter,
                      _precompute_data, _apply_one_dm_substep!

@testset "Allocation Regression" begin

    @testset "B_bohr allocations" begin
        config = make_config(Lindbladian(), BohrDomain())
        jump = TEST_JUMPS[1]
        budget = 64 * 1024

        # Single-jump (wrapped as vector): warmup + measure
        B_ref = B_bohr(TEST_HAM, JumpOp[jump], config)
        allocs = @allocated B_bohr(TEST_HAM, JumpOp[jump], config)
        @test allocs <= budget
        @info "B_bohr allocations (single-jump)" allocs_bytes=allocs threshold=budget

        # Multi-jump: warmup + measure
        B_ref_multi = B_bohr(TEST_HAM, TEST_JUMPS, config)
        allocs_multi = @allocated B_bohr(TEST_HAM, TEST_JUMPS, config)
        # The row cache reuses the same fixed-size buffers for every jump.
        @test allocs_multi <= budget
        @info "B_bohr allocations (multi-jump)" allocs_bytes=allocs_multi threshold=budget
    end

    @testset "B_time allocations" begin
        config = make_config(Lindbladian(), TimeDomain())
        precomputed = _precompute_data(config, TEST_HAM)
        (; b_minus, b_plus) = precomputed
        jump = TEST_JUMPS[1]

        # Single-jump (wrapped as vector) warmup + measure
        B_ref = B_time(JumpOp[jump], TEST_HAM, b_minus, b_plus, T0, BETA, SIGMA)
        allocs = @allocated B_time(JumpOp[jump], TEST_HAM, b_minus, b_plus, T0, BETA, SIGMA)
        # Budget: pre-allocated buffers (diag_u, diag_u2 vectors; b_plus_summand, tmp, M, B matrices)
        # and lazy adjoint views from mul! calls. Must NOT include per-iteration Diagonal wrapper allocations.
        # qf-6af.5: when threading the inner τ × jumps loop and outer t-loop,
        # B_time additionally allocates per-task buffers (one set per chunk
        # plus a partial accumulator) and Task objects. Multiplied by
        # nthreads, the budget grows linearly. Even with the wider budget
        # the test still rejects the per-iteration Diagonal regression
        # (~num_b_plus × d² ≈ 183 × 256 × 16 ≈ 750 KB at NUM_QUBITS=4).
        d = DIM
        nt = max(Threads.nthreads(), 1)
        max_expected = (25 + 8 * (nt - 1)) * d^2 * sizeof(ComplexF64) + 4096 * nt
        @test allocs <= max_expected  # B_time single-jump: allow buffers, catch Diagonal wrapper reintroduction
        @info "B_time allocations (single-jump)" allocs_bytes=allocs threshold=max_expected

        # Multi-jump warmup + measure
        B_ref_m = B_time(TEST_JUMPS, TEST_HAM, b_minus, b_plus, T0, BETA, SIGMA)
        allocs_m = @allocated B_time(TEST_JUMPS, TEST_HAM, b_minus, b_plus, T0, BETA, SIGMA)
        # Multi-jump: additional mul! per jump per b_plus iteration
        num_jumps = length(TEST_JUMPS)
        num_b_plus = length(b_plus)
        max_expected_multi = max_expected + num_jumps * num_b_plus * d * sizeof(ComplexF64) + 4096  # empirical: per-jump adjoint overhead
        @test allocs_m <= max_expected_multi  # B_time multi-jump: linear scaling with n_jumps * n_b_plus
        @info "B_time allocations (multi-jump)" allocs_bytes=allocs_m threshold=max_expected_multi num_jumps=num_jumps
    end

    @testset "B_trotter allocations" begin
        config = make_config(Lindbladian(), TrotterDomain())
        precomputed = _precompute_data(config, TEST_TROTTER)
        (; b_minus, b_plus) = precomputed

        # Use pre-built Trotter-basis jumps from test_helpers.jl
        trotter_jumps = TEST_TROTTER_JUMPS
        jump = trotter_jumps[1]

        # Single-jump (wrapped as vector) warmup + measure
        B_ref = B_trotter(JumpOp[jump], TEST_TROTTER, b_minus, b_plus, BETA, SIGMA)
        allocs = @allocated B_trotter(JumpOp[jump], TEST_TROTTER, b_minus, b_plus, BETA, SIGMA)
        d = DIM
        # qf-6af.5: same widening rationale as `B_time allocations` above —
        # threading adds per-task buffers + Task overhead at construction.
        nt = max(Threads.nthreads(), 1)
        max_expected = (25 + 8 * (nt - 1)) * d^2 * sizeof(ComplexF64) + 4096 * nt
        @test allocs <= max_expected  # B_trotter single-jump: same budget rationale as B_time
        @info "B_trotter allocations (single-jump)" allocs_bytes=allocs threshold=max_expected

        # Multi-jump warmup + measure
        B_ref_m = B_trotter(trotter_jumps, TEST_TROTTER, b_minus, b_plus, BETA, SIGMA)
        allocs_m = @allocated B_trotter(trotter_jumps, TEST_TROTTER, b_minus, b_plus, BETA, SIGMA)
        num_jumps = length(trotter_jumps)
        num_b_plus = length(b_plus)
        max_expected_multi = max_expected + num_jumps * num_b_plus * d * sizeof(ComplexF64) + 4096  # empirical: per-jump adjoint overhead
        @test allocs_m <= max_expected_multi  # B_trotter multi-jump: same scaling as B_time
        @info "B_trotter allocations (multi-jump)" allocs_bytes=allocs_m threshold=max_expected_multi num_jumps=num_jumps
    end

    @testset "Retained channel hot-path allocations" begin
        config_therm = make_config(Thermalize(), TimeDomain(); delta=0.01)
        precomputed = _precompute_data(config_therm, TEST_HAM)
        ws = Workspace(config_therm, TEST_HAM, TEST_JUMPS)
        scratch = ws.scratch
        U_coherent = ws.U_coherents === nothing ? nothing : ws.U_coherents[1]
        evolving_dm = copy(Matrix{ComplexF64}(TEST_GIBBS))

        _apply_one_dm_substep!(
            evolving_dm, scratch, TEST_JUMPS[1], U_coherent,
            ws.K0s[1], ws.U_residuals[1], TEST_HAM, config_therm,
            precomputed, precomputed.gamma_norm_factor,
        )

        evolving_dm .= Matrix{ComplexF64}(TEST_GIBBS)
        substep_allocs = @allocated _apply_one_dm_substep!(
            evolving_dm, scratch, TEST_JUMPS[1], U_coherent,
            ws.K0s[1], ws.U_residuals[1], TEST_HAM, config_therm,
            precomputed, precomputed.gamma_norm_factor,
        )

        # Thread-task overhead is expected; rebuilding R/K0/U_residual or running
        # an eigendecomposition inside the step is not.
        substep_budget = 16_384 * max(Threads.nthreads(), 1)
        @test substep_allocs < substep_budget
        @info "_apply_one_dm_substep! allocations (TimeDomain)" allocs_bytes=substep_allocs threshold=substep_budget

        rho = copy(Matrix{ComplexF64}(TEST_GIBBS))
        apply_delta_channel!(ws, rho, config_therm, TEST_HAM)
        channel_allocs = @allocated apply_delta_channel!(ws, rho, config_therm, TEST_HAM)
        channel_budget = substep_budget * length(TEST_JUMPS)
        @test channel_allocs < channel_budget
        @info "apply_delta_channel! allocations (TimeDomain)" allocs_bytes=channel_allocs threshold=channel_budget
    end

end
