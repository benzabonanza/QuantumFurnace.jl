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
- _jump_contribution! (Time/Trotter): filter+abs vector intermediates
"""

using QuantumFurnace: B_bohr, B_time, B_trotter,
                      _precompute_data, _jump_contribution!,
                      ThermalizeScratch,
                      _build_trajectory_workspace
using Random

@testset "Allocation Regression" begin

    @testset "B_bohr allocations" begin
        config = make_config(Lindbladian(), BohrDomain())
        jump = TEST_JUMPS[1]
        num_freqs = length(keys(TEST_HAM.bohr_dict))

        # Single-jump (wrapped as vector): warmup + measure
        B_ref = B_bohr(TEST_HAM, JumpOp[jump], config)
        allocs = @allocated B_bohr(TEST_HAM, JumpOp[jump], config)
        # Budget: return matrix + scratch + per-frequency broadcasting overhead.
        # Must NOT include per-frequency sparse matrix allocations (the eliminated O(num_freqs * dim^2) pattern).
        max_expected = num_freqs * DIM^2 * sizeof(ComplexF64)  # generous upper bound for broadcasting
        @test allocs <= max_expected  # B_bohr single-jump: allow buffers + broadcasting, catch sparse matrix reintroduction
        @info "B_bohr allocations (single-jump)" allocs_bytes=allocs threshold=max_expected num_freqs=num_freqs

        # Multi-jump: warmup + measure
        B_ref_multi = B_bohr(TEST_HAM, TEST_JUMPS, config)
        allocs_multi = @allocated B_bohr(TEST_HAM, TEST_JUMPS, config)
        # Multi-jump: overhead scales linearly with number of jumps.
        num_jumps = length(TEST_JUMPS)
        max_expected_multi = num_jumps * max_expected
        @test allocs_multi <= max_expected_multi  # B_bohr multi-jump: linear scaling with n_jumps
        @info "B_bohr allocations (multi-jump)" allocs_bytes=allocs_multi threshold=max_expected_multi num_jumps=num_jumps
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
        d = DIM
        max_expected = 25 * d^2 * sizeof(ComplexF64) + 4096  # empirical: buffers + mul! adjoint wrappers + headroom
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
        max_expected = 25 * d^2 * sizeof(ComplexF64) + 4096  # empirical: same structure as B_time
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

    @testset "_jump_contribution! Time/Trotter filter allocation" begin
        # Verify that _jump_contribution! for Time/Trotter thermalize does not allocate
        # filter+abs vectors. This function operates in-place on scratch buffers,
        # so the only allocations should be from _finalize_kraus_step! (eigen decomposition).
        config_therm = make_config(Thermalize(), TimeDomain(); delta=0.01)
        precomputed = _precompute_data(config_therm, TEST_HAM)
        CT = ComplexF64
        d = DIM
        scratch = ThermalizeScratch(CT, d)
        jump = TEST_JUMPS[1]
        evolving_dm = copy(Matrix{CT}(TEST_GIBBS))

        # Warmup
        _jump_contribution!(evolving_dm, jump, TEST_HAM, config_therm, precomputed, scratch)

        # Measure
        evolving_dm .= Matrix{CT}(TEST_GIBBS)
        allocs = @allocated _jump_contribution!(evolving_dm, jump, TEST_HAM, config_therm, precomputed, scratch)

        # Budget: eigen() + mul! temporaries in _finalize_kraus_step!.
        # Must NOT include filter+abs vector overhead (2 * num_energies * sizeof(Float64) per iteration).
        max_expected = 50 * d^2 * sizeof(CT)  # generous for eigen decomposition + scratch
        @test allocs < max_expected  # _jump_contribution! in-place: allow eigen, catch filter vector reintroduction
        @info "_jump_contribution! allocations (TimeDomain)" allocs_bytes=allocs threshold=max_expected
    end

    @testset "step_along_trajectory! allocations" begin
        # Guards against GC pressure that would destroy parallel scaling.
        # Uses the 3-qubit SMALL system for fast execution.
        config = make_config(Thermalize(), TimeDomain();
            num_qubits=3, delta=0.01, mixing_time=1.0, construction=GNS())
        CT = ComplexF64
        dim = N3_DIM  # 8
        ws = QuantumFurnace._build_trajectory_workspace(config, N3_HAM, N3_JUMPS; delta=0.01)

        psi0 = zeros(CT, dim)
        psi0[1] = 1.0

        # Wrap in a function to avoid top-level scope allocation artifacts.
        # Julia's @allocated in global/testset scope can show spurious allocations
        # from boxing local variables; a function barrier ensures proper optimization.
        function _measure_step_allocs(ws, psi0)
            psi = copy(psi0)
            rng = Xoshiro(999)

            # Warmup: JIT compile all code paths (no-jump, residual, dissipative-jump)
            for _ in 1:100
                step_along_trajectory!(psi, ws, rng)
            end

            # Reset and measure
            copyto!(psi, psi0)
            rng2 = Xoshiro(999)
            return @allocated step_along_trajectory!(psi, ws, rng2)
        end

        allocs = _measure_step_allocs(ws, psi0)
        @test allocs == 0  # Hot path must be allocation-free for parallel scaling (Union{Nothing,Function} boxing handled by MATVEC_ALLOC_BUDGET elsewhere)
        @info "step_along_trajectory! allocations" allocs_bytes=allocs threshold=0
    end

end
