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
                      KrausScratch, transform_jumps_to_basis,
                      TrajectoryWorkspace, _evolve_along_trajectory!
using Random

@testset "Allocation Regression" begin

    @testset "B_bohr allocations" begin
        config = make_liouv_config(BohrDomain())
        jump = TEST_JUMPS[1]
        num_freqs = length(keys(TEST_HAM.bohr_dict))

        # Single-jump: warmup + measure
        B_ref = B_bohr(TEST_HAM, jump, config)
        allocs = @allocated B_bohr(TEST_HAM, jump, config)
        # B_bohr allocates: return matrix B (dim x dim), scratch f_A_nu_1 (dim x dim),
        # plus per-frequency broadcasting overhead from closure-based @. computation.
        # It should NOT allocate sparse matrices per frequency (the eliminated pattern).
        # The old spzeros+scatter pattern would add ~num_freqs * dim^2 * sizeof(ComplexF64) bytes
        # from the dense result of sparse*dense multiply per iteration.
        # Threshold: allow current baseline allocations (broadcasting + buffers) with headroom,
        # but catch reintroduction of per-frequency sparse matrix allocations.
        max_expected = num_freqs * DIM^2 * sizeof(ComplexF64)
        @test allocs <= max_expected

        # Multi-jump: warmup + measure
        B_ref_multi = B_bohr(TEST_HAM, TEST_JUMPS, config)
        allocs_multi = @allocated B_bohr(TEST_HAM, TEST_JUMPS, config)
        # Multi-jump: overhead scales linearly with number of jumps.
        num_jumps = length(TEST_JUMPS)
        max_expected_multi = num_jumps * max_expected
        @test allocs_multi <= max_expected_multi
    end

    @testset "B_time allocations" begin
        config = make_liouv_config(TimeDomain())
        precomputed = _precompute_data(config, TEST_HAM)
        (; b_minus, b_plus) = precomputed
        jump = TEST_JUMPS[1]

        # Single-jump warmup + measure
        B_ref = B_time(jump, TEST_HAM, b_minus, b_plus, T0, BETA, SIGMA)
        allocs = @allocated B_time(jump, TEST_HAM, b_minus, b_plus, T0, BETA, SIGMA)
        # Expected allocations: pre-allocated buffers (diag_u, diag_u2 vectors; b_plus_summand,
        # tmp, M, B matrices) and lazy adjoint views from mul! calls.
        # It should NOT include per-iteration Diagonal wrapper allocations (the eliminated pattern).
        # The old pattern would add ~(num_b_plus + num_b_minus) * (vector + Diagonal struct)
        # allocations, easily exceeding 100 * d^2 bytes.
        d = DIM
        max_expected = 25 * d^2 * sizeof(ComplexF64) + 4096
        @test allocs <= max_expected

        # Multi-jump warmup + measure
        B_ref_m = B_time(TEST_JUMPS, TEST_HAM, b_minus, b_plus, T0, BETA, SIGMA)
        allocs_m = @allocated B_time(TEST_JUMPS, TEST_HAM, b_minus, b_plus, T0, BETA, SIGMA)
        # Multi-jump has additional mul!(M, jump_eig', tmp) per jump per b_plus iteration.
        # Each adjoint-mul may allocate a lazy wrapper. Allow proportional budget.
        num_jumps = length(TEST_JUMPS)
        num_b_plus = length(b_plus)
        max_expected_multi = max_expected + num_jumps * num_b_plus * d * sizeof(ComplexF64) + 4096
        @test allocs_m <= max_expected_multi
    end

    @testset "B_trotter allocations" begin
        config = make_liouv_config(TrotterDomain())
        precomputed = _precompute_data(config, TEST_TROTTER)
        (; b_minus, b_plus) = precomputed

        # Need Trotter-basis jumps (callers now transform before calling B_trotter)
        trotter_jumps = transform_jumps_to_basis(TEST_JUMPS, TEST_TROTTER.eigvecs)
        jump = trotter_jumps[1]

        # Single-jump warmup + measure
        B_ref = B_trotter(jump, TEST_TROTTER, b_minus, b_plus, BETA, SIGMA)
        allocs = @allocated B_trotter(jump, TEST_TROTTER, b_minus, b_plus, BETA, SIGMA)
        d = DIM
        max_expected = 25 * d^2 * sizeof(ComplexF64) + 4096
        @test allocs <= max_expected

        # Multi-jump warmup + measure
        B_ref_m = B_trotter(trotter_jumps, TEST_TROTTER, b_minus, b_plus, BETA, SIGMA)
        allocs_m = @allocated B_trotter(trotter_jumps, TEST_TROTTER, b_minus, b_plus, BETA, SIGMA)
        num_jumps = length(trotter_jumps)
        num_b_plus = length(b_plus)
        max_expected_multi = max_expected + num_jumps * num_b_plus * d * sizeof(ComplexF64) + 4096
        @test allocs_m <= max_expected_multi
    end

    @testset "_jump_contribution! Time/Trotter filter allocation" begin
        # Verify that _jump_contribution! for Time/Trotter thermalize does not allocate
        # filter+abs vectors. This function operates in-place on scratch buffers,
        # so the only allocations should be from _finalize_kraus_step! (eigen decomposition).
        config_therm = make_thermalize_config(TimeDomain(); delta=0.01)
        precomputed = _precompute_data(config_therm, TEST_HAM)
        CT = ComplexF64
        d = DIM
        scratch = KrausScratch(CT, d)
        jump = TEST_JUMPS[1]
        evolving_dm = copy(Matrix{CT}(TEST_GIBBS))

        # Warmup
        _jump_contribution!(evolving_dm, jump, TEST_HAM, config_therm, precomputed, scratch)

        # Measure
        evolving_dm .= Matrix{CT}(TEST_GIBBS)
        allocs = @allocated _jump_contribution!(evolving_dm, jump, TEST_HAM, config_therm, precomputed, scratch)

        # _jump_contribution! should be mostly in-place. It still has mul! calls and
        # the eigen() in _finalize_kraus_step! which allocates. The key test:
        # allocations should NOT include the filter+abs vectors (2 * num_energies * sizeof(Float64)).
        # Since _finalize_kraus_step! inherently allocates (eigen decomposition), we test
        # that allocations are bounded and don't include the filter overhead.
        # Allow generous budget for eigen + mul! temporaries, but NOT for filter.
        @test allocs < 50 * d^2 * sizeof(CT)
    end

    @testset "step_along_trajectory! allocations" begin
        # Guards against GC pressure that would destroy parallel scaling.
        # Uses the 3-qubit SMALL system for fast execution.
        config = make_small_thermalize_config(TimeDomain();
            delta=0.01, mixing_time=1.0, with_coherent=false)
        precomputed = _precompute_data(config, SMALL_HAM)
        CT = ComplexF64
        dim = SMALL_DIM  # 8
        scratch = KrausScratch(CT, dim)
        fw = build_trajectoryframework(SMALL_JUMPS, SMALL_HAM, config, precomputed, scratch, 0.01)

        ws = TrajectoryWorkspace(CT, dim)

        psi0 = zeros(CT, dim)
        psi0[1] = 1.0

        # Wrap in a function to avoid top-level scope allocation artifacts.
        # Julia's @allocated in global/testset scope can show spurious allocations
        # from boxing local variables; a function barrier ensures proper optimization.
        function _measure_step_allocs(fw, ws, psi0)
            psi = copy(psi0)
            rng = Xoshiro(999)

            # Warmup: JIT compile all code paths (no-jump, residual, dissipative-jump)
            for _ in 1:100
                step_along_trajectory!(psi, fw, ws, rng)
            end

            # Reset and measure
            copyto!(psi, psi0)
            rng2 = Xoshiro(999)
            return @allocated step_along_trajectory!(psi, fw, ws, rng2)
        end

        allocs = _measure_step_allocs(fw, ws, psi0)
        @test allocs == 0
    end

end
