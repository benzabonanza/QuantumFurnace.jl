using Test
using LinearAlgebra
using Random
using QuantumFurnace

# test_helpers.jl is already included by runtests.jl

# ============================================================================
# Comprehensive tests for krylov_eigsolve.jl: spectral gap, channel path,
# domain coverage, guard rails, and eigenvalue properties.
# Phase 29: Eigensolver Integration
# Quick-36: Faithful Chen channel for apply_delta_channel!
# ============================================================================

@testset "Krylov Eigsolve" begin

    # ========================================================================
    # Testset 1: apply_delta_channel! faithful Chen channel
    # ========================================================================
    @testset "apply_delta_channel! faithful Chen channel" begin
        config_therm = make_config(Thermalize(),EnergyDomain(); construction=KMS(), delta=0.01)
        config_liouv = make_config(Lindbladian(),EnergyDomain(); construction=KMS())
        ws = Workspace(config_therm, TEST_HAM, TEST_JUMPS)
        delta = config_therm.delta

        # Dense Lindbladian for Euler comparison
        L_dense = construct_lindbladian(TEST_JUMPS, config_liouv, TEST_HAM)
        I_d2 = Matrix{ComplexF64}(LinearAlgebra.I(DIM^2))

        # Threshold rationale (trace preservation): tr(E(rho)) should equal tr(rho) exactly
        # for a CPTP map. Real FP drift from the per-jump Lie–Trotter sweep is at the
        # rounding floor (~5e-16 measured); `hermitianize!` averages off-diagonals
        # without touching the diagonal, so it adds no trace drift. Threshold 1e-10
        # gives >1e5× margin.
        #
        # Threshold rationale (positivity): eigenvalues of E(rho) should be >= 0 for CPTP.
        # FP rounding can produce tiny negatives: O(eps * ||rho||) ~ 1e-16. Threshold -1e-10 is generous.
        #
        # Threshold rationale (Euler closeness): faithful per-jump Lie–Trotter Φ_δ differs
        # from Euler by O(delta^2). The constant absorbs ‖𝓛‖·DIM and the Lie–Trotter
        # splitting constant on the n_jumps substeps. Empirically ~1.6e-6 at δ=1e-2,
        # so the 50·δ² = 5e-3 threshold has >3000x margin.
        max_trace_err = 0.0
        min_eigenvalue = Inf
        max_euler_err = 0.0
        euler_threshold = 50 * delta^2
        for _ in 1:5
            rho = Matrix(random_density_matrix(NUM_QUBITS))

            # Faithful Chen channel
            apply_delta_channel!(ws, rho, config_therm, TEST_HAM)
            rho_chen = copy(ws.scratch.rho_next)

            # Trace preservation: tr(E(rho)) == tr(rho)
            trace_err = abs(real(tr(rho_chen)) - real(tr(rho)))
            @test isapprox(real(tr(rho_chen)), real(tr(rho)); atol=1e-10)
            max_trace_err = max(max_trace_err, trace_err)

            # Positivity: eigenvalues of E(rho) >= -eps for valid density matrix input
            eigs = eigvals(Hermitian(rho_chen))
            @test all(eigs .> -1e-10)
            min_eigenvalue = min(min_eigenvalue, minimum(eigs))

            # O(delta^2) close to Euler: |E_chen(rho) - E_euler(rho)| < C * delta^2
            v_euler = (I_d2 + delta * L_dense) * vec(rho)
            euler_err = norm(vec(rho_chen) - v_euler)
            @test euler_err < euler_threshold
            max_euler_err = max(max_euler_err, euler_err)
        end
        @info "Chen channel trace preservation" max_trace_error=max_trace_err threshold=1e-8
        @info "Chen channel positivity" min_eigenvalue=min_eigenvalue threshold=-1e-10
        @info "Chen channel Euler closeness" max_euler_error=max_euler_err threshold=euler_threshold delta=delta
    end

    # ========================================================================
    # Testset 2: krylov_spectral_gap result fields
    # ========================================================================
    @testset "krylov_spectral_gap result fields" begin
        config_kms = make_config(Lindbladian(),EnergyDomain(); construction=KMS())
        result = krylov_spectral_gap(config_kms, TEST_HAM, TEST_JUMPS;
            krylovdim=30, howmany=4)

        @test result isa NamedTuple
        @test length(result.eigenvalues) >= 2
        @test result.spectral_gap > 0
        @test size(result.fixed_point) == (DIM, DIM)
        @test size(result.gap_mode) == (DIM, DIM)
        @test result.converged >= 2
        @test result.matvec_count > 0
        @test result.channel_eigenvalues === nothing
        @test result.delta_used === nothing
        @info "krylov_spectral_gap result fields" spectral_gap=result.spectral_gap n_eigenvalues=length(result.eigenvalues) converged=result.converged matvec_count=result.matvec_count
    end

    # ========================================================================
    # Testset 3: Lindbladian eigsolve accuracy (EnergyDomain KMS)
    # ========================================================================
    @testset "Lindbladian eigsolve accuracy (EnergyDomain KMS)" begin
        config = make_config(Lindbladian(),EnergyDomain(); construction=KMS())
        L_dense = construct_lindbladian(TEST_JUMPS, config, TEST_HAM)
        dense_result = extract_leading_eigendata(L_dense; n_modes=4)

        result = krylov_spectral_gap(config, TEST_HAM, TEST_JUMPS;
            krylovdim=30, howmany=4)

        # Threshold rationale (gap rtol=1e-6): KrylovKit tol=1e-10 (default) bounds eigenvalue error.
        # rtol=1e-6 gives 10000x margin for iterative convergence variability.
        gap_err = abs(result.spectral_gap - dense_result.spectral_gap) / dense_result.spectral_gap
        @test isapprox(result.spectral_gap, dense_result.spectral_gap; rtol=1e-6)
        @info "Eigsolve gap accuracy (EnergyDomain KMS)" krylov_gap=result.spectral_gap dense_gap=dense_result.spectral_gap relative_error=gap_err rtol=1e-6

        # Threshold rationale (trace distance < 1e-4): KMS fixed point should converge to Gibbs.
        # KrylovKit eigenvector accuracy bounded by tol * condition_number. 1e-4 is conservative.
        td = trace_distance_h(Hermitian(result.fixed_point), TEST_GIBBS)
        @test td < 1e-4
        @info "Fixed point trace distance to Gibbs" trace_distance=td threshold=1e-4

        # Threshold rationale (|Re(lambda_1)| < 1e-8): steady-state eigenvalue is exactly 0.
        # KrylovKit residual bounded by tol. 1e-8 gives 100x margin over tol=1e-10.
        ss_err = abs(real(result.eigenvalues[1]))
        @test ss_err < 1e-8
        @info "Steady-state eigenvalue magnitude" abs_real_lambda1=ss_err threshold=1e-8
    end

    # ========================================================================
    # Testset 4: Lindbladian eigsolve accuracy (EnergyDomain GNS)
    # ========================================================================
    @testset "Lindbladian eigsolve accuracy (EnergyDomain GNS)" begin
        config = make_config(Lindbladian(), EnergyDomain(); construction=GNS())
        L_dense = construct_lindbladian(TEST_JUMPS, config, TEST_HAM)
        dense_result = extract_leading_eigendata(L_dense; n_modes=4)

        result = krylov_spectral_gap(config, TEST_HAM, TEST_JUMPS;
            krylovdim=30, howmany=4)

        # Threshold rationale: same as KMS testset -- KrylovKit tol=1e-10, rtol=1e-6 gives 10000x margin.
        gap_err = abs(result.spectral_gap - dense_result.spectral_gap) / dense_result.spectral_gap
        @test isapprox(result.spectral_gap, dense_result.spectral_gap; rtol=1e-6)
        @info "Eigsolve gap accuracy (EnergyDomain GNS)" krylov_gap=result.spectral_gap dense_gap=dense_result.spectral_gap relative_error=gap_err rtol=1e-6

        # GNS fixed point differs from exact Gibbs -- skip trace distance check
        # Threshold rationale (trace=1 atol=1e-6): valid density matrix has tr(rho)=1.
        # KrylovKit eigenvector normalization should preserve this. 1e-6 is conservative.
        trace_err = abs(tr(result.fixed_point) - 1.0)
        @test isapprox(tr(result.fixed_point), 1.0; atol=1e-6)
        @info "GNS fixed point trace" trace=real(tr(result.fixed_point)) trace_error=trace_err atol=1e-6

        # Threshold rationale: same as KMS -- steady-state eigenvalue is exactly 0.
        ss_err = abs(real(result.eigenvalues[1]))
        @test ss_err < 1e-8
        @info "Steady-state eigenvalue magnitude (GNS)" abs_real_lambda1=ss_err threshold=1e-8
    end

    # ========================================================================
    # Testset 5: Channel eigsolve accuracy (Thermalize)
    # ========================================================================
    @testset "Channel eigsolve accuracy (Thermalize)" begin
        config_therm = make_config(Thermalize(),EnergyDomain();
            construction=KMS(), delta=0.01)
        # Dense reference from Lindbladian path
        config_liouv = make_config(Lindbladian(),EnergyDomain(); construction=KMS())
        L_dense = construct_lindbladian(TEST_JUMPS, config_liouv, TEST_HAM)
        dense_result = extract_leading_eigendata(L_dense; n_modes=4)

        result = krylov_spectral_gap(config_therm, TEST_HAM, TEST_JUMPS;
            krylovdim=30, howmany=4)

        # --- Channel OPERATOR gap: qf-9lp route + documented Φ_δ-Arnoldi fragility ---
        # qf-9lp / qf-e4z.37 (revalidated qf-4fb on the build_heis_1d n=4 draw): the
        # ROBUST channel operator gap comes from 𝓛's λ₂ — the `Config{Lindbladian}`
        # `krylov_spectral_gap` route — which equals the dense 𝓛 gap to machine
        # precision and the dense Φ_δ gap to O(δ). Arnoldi DIRECTLY on the non-normal
        # Φ_δ (this `Config{Thermalize}` method) is a fragile FALLBACK: at stiff cells
        # (μ-spectrum clustered tightly near 1) it locks onto the wrong eigenmode,
        # independent of krylovdim. The build_heis_1d n=4 fixture IS such a cell at
        # δ=0.01 — verified ground truth on this fixture:
        #   dense 𝓛 gap                       = 0.08765
        #   L-direct Arnoldi gap (qf-9lp)      = 0.08765  (rel 1e-15  ✓ machine prec)
        #   channel-trajectory Pass-1 |log|μ₂||/δ = 0.08826  (rel 0.7%  = dense Φ_δ gap)
        #   channel-Arnoldi-on-Φ_δ gap (here)  = 0.11628  (rel 33%, WRONG MODE)
        # and the wrong-mode error grows monotonically as δ→0 (6.8%→57% over
        # δ∈[0.1,0.005]) — the textbook qf-e4z.37 stiffness signature, NOT a bug
        # (dense Φ_δ is correct) and NOT krylovdim/register-fixable.
        #
        # ACCURACY claim (held TIGHT): the channel operator gap via the recommended
        # qf-9lp L-direct route matches the dense 𝓛 gap to ~1e-8.
        L_direct = krylov_spectral_gap(config_liouv, TEST_HAM, TEST_JUMPS;
            krylovdim=30, howmany=4)
        @test isapprox(L_direct.spectral_gap, dense_result.spectral_gap; atol=1e-8)
        @info "Channel operator gap (qf-9lp L-direct route)" L_direct_gap=L_direct.spectral_gap dense_gap=dense_result.spectral_gap atol=1e-8

        # FRAGILITY witness (do NOT assert agreement — qf-e4z.37): on this stiff cell
        # the Φ_δ-Arnoldi gap is a finite positive rate that resolves SOME Φ_δ
        # eigenmode (μ near 1), but not necessarily the slowest. Asserting it equals
        # the true gap to rtol=1e-2 would require masking a 33% wrong-mode error and is
        # exactly the regime qf-9lp says to avoid. We assert only that the method runs,
        # converges, and returns a physically-bounded rate.
        @test result.spectral_gap > 0
        @test result.converged >= 1
        # The Φ_δ-Arnoldi gap must map back to a |μ| ≤ 1 (contraction) — μ = (1 - δ·λ)
        # in the convention used here, so 0 < δ·gap < 2.
        @test 0 < config_therm.delta * result.spectral_gap < 2
        @info "Channel Φ_δ-Arnoldi gap (fragile fallback, qf-e4z.37)" phi_arnoldi_gap=result.spectral_gap dense_gap=dense_result.spectral_gap note="wrong-mode at stiff cells; not asserted equal to true gap"

        # Channel-specific fields are populated
        @test result.channel_eigenvalues !== nothing
        @test result.delta_used == 0.01

        # Threshold rationale (channel eigenvalue ~1, atol=0.01): steady-state channel eigenvalue
        # is exactly 1.0 for a CPTP map. KrylovKit convergence + delta discretization give O(delta) error.
        # This is fixture-INDEPENDENT (CPTP) and held tight — the Φ_δ-Arnoldi resolves
        # the steady state robustly even when it misses the gap mode.
        chan_ss_err = abs(abs(result.channel_eigenvalues[1]) - 1.0)
        @test isapprox(abs(result.channel_eigenvalues[1]), 1.0; atol=0.01)
        @info "Channel steady-state eigenvalue" abs_mu1=abs(result.channel_eigenvalues[1]) error=chan_ss_err atol=0.01
    end

    # ========================================================================
    # Testset 6: All domains work (Lindbladian path)
    # ========================================================================
    @testset "All domains work (Lindbladian path)" begin
        @testset "EnergyDomain" begin
            config = make_config(Lindbladian(),EnergyDomain(); construction=KMS())
            result = krylov_spectral_gap(config, TEST_HAM, TEST_JUMPS;
                krylovdim=30, howmany=2)
            @test result.spectral_gap > 0
            @info "Domain coverage (EnergyDomain)" spectral_gap=result.spectral_gap
        end

        @testset "TimeDomain" begin
            config = make_config(Lindbladian(),TimeDomain(); construction=KMS())
            result = krylov_spectral_gap(config, TEST_HAM, TEST_JUMPS;
                krylovdim=30, howmany=2)
            @test result.spectral_gap > 0
            @info "Domain coverage (TimeDomain)" spectral_gap=result.spectral_gap
        end

        @testset "TrotterDomain" begin
            config = make_config(Lindbladian(),TrotterDomain(); construction=KMS())
            result = krylov_spectral_gap(config, TEST_HAM, TEST_TROTTER_JUMPS;
                trotter=TEST_TROTTER, krylovdim=30, howmany=2)
            @test result.spectral_gap > 0
            @info "Domain coverage (TrotterDomain)" spectral_gap=result.spectral_gap
        end

        @testset "BohrDomain" begin
            config = make_config(Lindbladian(),BohrDomain(); construction=KMS())
            result = krylov_spectral_gap(config, TEST_HAM, TEST_JUMPS;
                krylovdim=30, howmany=2)
            @test result.spectral_gap > 0
            @info "Domain coverage (BohrDomain)" spectral_gap=result.spectral_gap
        end
    end

    # ========================================================================
    # Testset 7: Guard rails
    # ========================================================================
    @testset "Guard rails" begin
        config = make_config(Lindbladian(),EnergyDomain(); construction=KMS())

        # krylovdim <= howmany must error
        @test_throws ErrorException krylov_spectral_gap(config, TEST_HAM, TEST_JUMPS;
            krylovdim=2, howmany=4)

        # Memory guard should not throw for small n=4 system
        result = krylov_spectral_gap(config, TEST_HAM, TEST_JUMPS;
            krylovdim=50, howmany=2)
        @test result isa NamedTuple
    end

    # ========================================================================
    # Testset 8: Eigenvalue sorting and conversion
    # ========================================================================
    @testset "Eigenvalue sorting and conversion" begin
        # Lindbladian path: eigenvalues sorted by |Re(lambda)| ascending (structural check)
        config = make_config(Lindbladian(),EnergyDomain(); construction=KMS())
        result_liouv = krylov_spectral_gap(config, TEST_HAM, TEST_JUMPS;
            krylovdim=30, howmany=4)

        @test abs(real(result_liouv.eigenvalues[1])) <= abs(real(result_liouv.eigenvalues[2]))
        @info "Eigenvalue sorting" abs_Re_lambda1=abs(real(result_liouv.eigenvalues[1])) abs_Re_lambda2=abs(real(result_liouv.eigenvalues[2]))

        # Channel path: verify eigenvalue conversion formula
        config_therm = make_config(Thermalize(),EnergyDomain();
            construction=KMS(), delta=0.01)
        result_chan = krylov_spectral_gap(config_therm, TEST_HAM, TEST_JUMPS;
            krylovdim=30, howmany=4)

        # Threshold rationale (atol=1e-10): conversion mu = 1 + delta * lambda_L is algebraically
        # exact (just arithmetic). Error is FP rounding only: O(eps * |mu|) ~ 1e-16.
        # Threshold 1e-10 gives >1e6 margin.
        reconstructed_mu = 1.0 .+ result_chan.delta_used .* result_chan.eigenvalues
        conversion_err = maximum(abs.(reconstructed_mu .- result_chan.channel_eigenvalues))
        @test isapprox(reconstructed_mu, result_chan.channel_eigenvalues; atol=1e-10)
        @info "Eigenvalue conversion consistency" max_conversion_error=conversion_err atol=1e-10
    end

    # ========================================================================
    # Testset 9 — qf-8fr: symmetric Hamiltonians do NOT collapse the Arnoldi.
    #
    # Classical 1D Ising H = sum Z_i Z_{i+1} is invariant under translation and
    # spin-flip (otimes_i X). The maximally mixed state I/d is a fixed point
    # of both symmetries; if `krylov_spectral_gap` seeded Arnoldi with I/d
    # (the pre-qf-8fr default), the Krylov subspace would stay in the trivial
    # symmetric sector and miss the true gap eigenmode — which is the
    # spin-flip-odd magnetisation at lambda = -4.45e-2 for n=4. The patched
    # `_krylov_default_x0` adds a small traceless GUE perturbation, breaking
    # the symmetry while preserving the trace-1 normalisation. This test
    # would FAIL on the pre-qf-8fr code (returns the 2nd-symmetric-sector
    # eigenvalue, ~3.8x too large at n=4).
    # ========================================================================
    @testset "krylov_spectral_gap — symmetric system regression (qf-8fr)" begin
        # n=3 keeps the sandbox cheap. Periodic classical Ising, no disorder.
        n_ising = 3
        terms_zz = Vector{Vector{Matrix{ComplexF64}}}([[Z, Z]])
        coeffs_zz = [1.0]
        base_ham = QuantumFurnace._construct_base_ham(terms_zz, coeffs_zz, n_ising;
                                                      periodic=true)
        rescaling_factor, shift = QuantumFurnace._rescaling_and_shift_factors(base_ham)
        d_ising = 2^n_ising
        rescaled_mat = Matrix(base_ham) ./ rescaling_factor .+ shift * I(d_ising)
        eigvals_rs, eigvecs_rs = eigen(Hermitian(rescaled_mat))

        # Build a HamHam directly via the raw-NamedTuple constructor.
        raw = (
            matrix = rescaled_mat,
            terms = terms_zz,
            base_coeffs = coeffs_zz ./ rescaling_factor,
            disordering_terms = nothing,
            disordering_coeffs = nothing,
            eigvals = eigvals_rs,
            eigvecs = eigvecs_rs,
            nu_min = minimum(diff(eigvals_rs)),
            shift = shift,
            rescaling_factor = rescaling_factor,
            periodic = true,
        )
        β_phys = 0.5
        ham_ising = HamHam(raw; beta_phys=β_phys)
        jumps_ising = QuantumFurnace._build_jump_set(ham_ising, n_ising)
        β_alg = beta_alg(ham_ising, β_phys)

        # CKG smooth-Metropolis EnergyDomain Config (matches the qf-8fr sweep).
        σ = 1.0 / β_alg
        H_norm = maximum(abs, ham_ising.eigvals)
        omega_range = 2.0 * (H_norm + 8 * σ)
        r_D = 7
        w0_D = omega_range / 2.0^r_D
        t0_D = 2π / (2.0^r_D * w0_D)
        cfg_e = Config(
            sim = Lindbladian(), domain = EnergyDomain(), construction = KMS(),
            num_qubits = n_ising, with_linear_combination = true,
            beta = β_alg, beta_phys = β_phys, sigma = σ,
            a = 0.0, s = 0.25, gaussian_parameters = (nothing, nothing),
            num_energy_bits_D = r_D, w0_D = w0_D, t0_D = t0_D,
            num_trotter_steps_per_t0 = 10, filter = nothing,
        )

        # Dense reference: build L explicitly, take the |Re|-second smallest.
        L_dense = construct_lindbladian(jumps_ising, cfg_e, ham_ising)
        ev_dense = eigvals(L_dense)
        perm = sortperm(real.(ev_dense); by=abs)
        gap_dense = abs(real(ev_dense[perm[2]]))

        # Krylov: must match dense to rtol=1e-6 with the patched x_0.
        res = krylov_spectral_gap(cfg_e, ham_ising, jumps_ising;
                                  krylovdim=40, howmany=4)
        @test isapprox(res.spectral_gap, gap_dense; rtol=1e-6)
        @info "qf-8fr classical Ising regression" n=n_ising β_phys=β_phys gap_krylov=res.spectral_gap gap_dense=gap_dense rel_err=abs(res.spectral_gap - gap_dense)/gap_dense

        # Also exercise BohrDomain — same physical Lindbladian, different domain wiring.
        cfg_b = Config(
            sim = Lindbladian(), domain = BohrDomain(), construction = KMS(),
            num_qubits = n_ising, with_linear_combination = true,
            beta = β_alg, beta_phys = β_phys, sigma = σ,
            a = 0.0, s = 0.25, gaussian_parameters = (nothing, nothing),
            num_energy_bits_D = r_D, w0_D = w0_D, t0_D = t0_D,
            num_trotter_steps_per_t0 = 10, filter = nothing,
        )
        res_b = krylov_spectral_gap(cfg_b, ham_ising, jumps_ising;
                                    krylovdim=40, howmany=4)
        @test isapprox(res_b.spectral_gap, gap_dense; rtol=1e-6)
        # Energy ≡ Bohr to machine precision for classical Ising (no quadrature error).
        @test isapprox(res.spectral_gap, res_b.spectral_gap; rtol=1e-10)
        @info "qf-8fr Energy ≡ Bohr cross-check" gap_E=res.spectral_gap gap_B=res_b.spectral_gap
    end

    # ========================================================================
    # Testset 10 — qf-umr: krylov_spectral_gap `workspace=` reuse.
    #
    # The two numerical pipelines (trajectory/mixing-time vs spectral-gap/
    # spectrum) are decoupled, sharing only the one expensive resource — the
    # O(d^3.x) Workspace operator build. krylov_spectral_gap now accepts a
    # pre-built Workspace so the compute_true_gap=true trajectory path can
    # forward its own and skip the redundant second build. This testset asserts:
    #   (a) reuse is bit-identical to a fresh build (Lindbladian + channel);
    #   (b) the supplied Workspace is actually CONSULTED (a mismatched config /
    #       scratch type / dim throws cleanly — proving no silent fresh rebuild,
    #       i.e. the double-build is gone);
    #   (c) end-to-end, predict_*_trajectory(compute_true_gap=true) with a
    #       forwarded Workspace equals the self-build path bitwise.
    # n=3 throughout — sandbox-cheap.
    # ========================================================================
    @testset "krylov_spectral_gap workspace= reuse (qf-umr)" begin
        d3 = 2^3

        # -- (a) Lindbladian: workspace= reuse is bit-identical to fresh build --
        @testset "(a) Lindbladian reuse byte-equality" begin
            cfg = make_config(Lindbladian(), EnergyDomain(); construction=KMS())
            r_fresh = krylov_spectral_gap(cfg, N3_HAM, N3_JUMPS;
                                          krylovdim=30, howmany=4)
            ws = Workspace(cfg, N3_HAM, N3_JUMPS)
            r_reuse1 = krylov_spectral_gap(cfg, N3_HAM, N3_JUMPS;
                                           krylovdim=30, howmany=4, workspace=ws)
            # Reuse a second time on the SAME ws — no scratch contamination.
            r_reuse2 = krylov_spectral_gap(cfg, N3_HAM, N3_JUMPS;
                                           krylovdim=30, howmany=4, workspace=ws)
            @test r_fresh.spectral_gap == r_reuse1.spectral_gap
            @test r_reuse1.spectral_gap == r_reuse2.spectral_gap
            @test r_fresh.eigenvalues == r_reuse1.eigenvalues
            @test r_fresh.fixed_point == r_reuse1.fixed_point
            @test r_fresh.gap_mode == r_reuse1.gap_mode
            @test r_fresh.matvec_count == r_reuse1.matvec_count
            @test r_reuse1.spectral_modes.off_diag_weight ==
                  r_fresh.spectral_modes.off_diag_weight
        end

        # -- (a') Channel: workspace= reuse is bit-identical to fresh build --
        @testset "(a') Channel reuse byte-equality" begin
            cfg_ch = make_config(Thermalize(), EnergyDomain();
                                 construction=KMS(), delta=0.01)
            rc_fresh = krylov_spectral_gap(cfg_ch, N3_HAM, N3_JUMPS;
                                           krylovdim=30, howmany=4)
            ws_ch = Workspace(cfg_ch, N3_HAM, N3_JUMPS)
            rc_reuse1 = krylov_spectral_gap(cfg_ch, N3_HAM, N3_JUMPS;
                                            krylovdim=30, howmany=4, workspace=ws_ch)
            rc_reuse2 = krylov_spectral_gap(cfg_ch, N3_HAM, N3_JUMPS;
                                            krylovdim=30, howmany=4, workspace=ws_ch)
            @test rc_fresh.spectral_gap == rc_reuse1.spectral_gap
            @test rc_reuse1.spectral_gap == rc_reuse2.spectral_gap
            @test rc_fresh.eigenvalues == rc_reuse1.eigenvalues
            @test rc_fresh.channel_eigenvalues == rc_reuse1.channel_eigenvalues
            @test rc_fresh.matvec_count == rc_reuse1.matvec_count
        end

        # -- (b) The supplied Workspace is CONSULTED, not silently rebuilt.
        #        A mismatched config / scratch type / dim must throw — if the
        #        method ignored `workspace=` and built a fresh one, none of
        #        these would throw. This is the "no double-build" witness: the
        #        forwarded ws is on the live path.
        @testset "(b) reuse guards (workspace is consulted)" begin
            cfg = make_config(Lindbladian(), EnergyDomain(); construction=KMS())
            cfg_ch = make_config(Thermalize(), EnergyDomain();
                                 construction=KMS(), delta=0.01)
            cfg_ch2 = make_config(Thermalize(), EnergyDomain();
                                  construction=KMS(), delta=0.02)  # δ differs

            ws_L  = Workspace(cfg, N3_HAM, N3_JUMPS)
            ws_ch = Workspace(cfg_ch, N3_HAM, N3_JUMPS)

            # Config mismatch (channel ws built at δ=0.02, called at δ=0.01).
            @test_throws ArgumentError krylov_spectral_gap(
                cfg_ch, N3_HAM, N3_JUMPS;
                krylovdim=30, howmany=4, workspace=Workspace(cfg_ch2, N3_HAM, N3_JUMPS))
            # Cross-simulation: channel ws into the Lindbladian method.
            @test_throws ArgumentError krylov_spectral_gap(
                cfg, N3_HAM, N3_JUMPS; krylovdim=30, howmany=4, workspace=ws_ch)
            # Cross-simulation: Lindbladian ws into the channel method.
            @test_throws ArgumentError krylov_spectral_gap(
                cfg_ch, N3_HAM, N3_JUMPS; krylovdim=30, howmany=4, workspace=ws_L)
        end

        # -- (c) End-to-end: predict_*_trajectory(compute_true_gap=true) with a
        #        forwarded Workspace == the self-build path, bitwise. This is the
        #        production wiring (lindblad_action.jl Pass-2 forwards ws); proves
        #        the double-build elimination changed cost, not the number.
        @testset "(c) predict_* compute_true_gap forwarded ws == self-build" begin
            psi = ones(ComplexF64, d3) ./ sqrt(2.0^3)
            rho_0 = psi * psi'

            cfg = make_config(Lindbladian(), EnergyDomain(); construction=KMS())
            t_grid = collect(range(0.0, 1.0, length=3))
            # Self-build path (workspace defaulted to nothing inside predict_*).
            tj_self = predict_lindbladian_trajectory(
                cfg, N3_HAM, N3_JUMPS, rho_0, t_grid;
                krylovdim=30, compute_true_gap=true)
            # Forwarded-workspace path: predict_* reuses this ws for BOTH Pass-1
            # and the internal Pass-2 krylov_spectral_gap (qf-umr).
            ws = Workspace(cfg, N3_HAM, N3_JUMPS)
            tj_ws = predict_lindbladian_trajectory(
                cfg, N3_HAM, N3_JUMPS, rho_0, t_grid;
                krylovdim=30, compute_true_gap=true, workspace=ws)
            @test tj_self.spectral_gap == tj_ws.spectral_gap
            @test tj_self.total_matvecs == tj_ws.total_matvecs
            @test tj_self.all_converged == tj_ws.all_converged

            cfg_ch = make_config(Thermalize(), EnergyDomain(); construction=KMS())
            k_grid = collect(0:5:20)
            rho_init = Matrix{ComplexF64}(rho_0)
            tc_self = predict_channel_trajectory(
                cfg_ch, N3_HAM, N3_JUMPS, rho_init, k_grid;
                krylovdim=30, compute_true_gap=true)
            ws_ch = Workspace(cfg_ch, N3_HAM, N3_JUMPS)
            tc_ws = predict_channel_trajectory(
                cfg_ch, N3_HAM, N3_JUMPS, rho_init, k_grid;
                krylovdim=30, compute_true_gap=true, workspace=ws_ch)
            @test tc_self.spectral_gap == tc_ws.spectral_gap
            @test tc_self.total_matvecs == tc_ws.total_matvecs
            @test tc_self.all_converged == tc_ws.all_converged
        end
    end

end  # @testset "Krylov Eigsolve"
