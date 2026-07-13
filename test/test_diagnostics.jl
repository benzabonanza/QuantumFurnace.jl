@testset "Diagnostics (Phase 26)" begin

    # Build the 3-qubit Lindbladian (BohrDomain for exact Gibbs fixed point)
    config = make_config(Lindbladian(), BohrDomain(); num_qubits=3, construction=KMS())
    L_sparse = construct_lindbladian(N3_JUMPS, config, N3_HAM)
    L_dense = Matrix{ComplexF64}(L_sparse)

    # -----------------------------------------------------------------------
    # DIAG-01: Eigendata extraction
    # -----------------------------------------------------------------------
    @testset "DIAG-01: extract_leading_eigendata" begin
        result = extract_leading_eigendata(L_dense; n_modes=10)

        @test result isa EigenDecompositionResult

        # Correct number of modes
        @test length(result.eigenvalues) == 10

        # Eigenvalues sorted by |Re(lambda)| (ascending)
        @test issorted(abs.(real.(result.eigenvalues)))

        # First eigenvalue near zero (steady state)
        # Dense eigen() on DIM^2=64 matrix; steady-state eigenvalue is exactly 0 in theory,
        # error is O(DIM^2 * eps) ~ 64 * 2.2e-16 ~ 1.4e-14. Threshold 1e-10 gives ~7000x margin.
        @test abs(result.eigenvalues[1]) < 1e-10
        @info "DIAG-01: steady-state eigenvalue" abs_lambda1=abs(result.eigenvalues[1]) threshold=1e-10

        # Spectral gap is positive
        @test result.spectral_gap == abs(real(result.eigenvalues[2]))
        @test result.spectral_gap > 0.0
        @info "DIAG-01: spectral gap" gap=result.spectral_gap

        # Right/left eigenvector dimensions
        @test size(result.right_eigenvectors) == (N3_DIM^2, 10)
        @test size(result.left_eigenvectors) == (N3_DIM^2, 10)

        # Biorthonormality: V_left' * V_right approx I
        # Round-trip between two dense eigen decompositions on DIM^2=64 matrix.
        # Error accumulates as O(DIM^2 * eps) per entry ~ 1.4e-14; for the full 10x10 product
        # with 64-element dot products, expect ~ O(DIM^2 * eps * sqrt(10)) ~ 5e-14.
        # Threshold 1e-8 gives ~200,000x margin.
        biorth = result.left_eigenvectors' * result.right_eigenvectors
        biorth_err = maximum(abs.(biorth - I(10)))
        @test isapprox(biorth, I(10); atol=1e-8)
        @info "DIAG-01: biorthonormality" max_error=biorth_err threshold_atol=1e-8

        # Im/Re ratios
        @test length(result.im_re_ratios) == 10
        @test result.im_re_ratios[1] == 0.0  # Steady state
        @test all(r -> r >= 0.0, result.im_re_ratios)

        # Eigenvectors actually satisfy eigenvalue equation: L * r_k = lambda_k * r_k
        # Dense matrix-vector product on DIM^2=64 matrix. Error O(DIM^2 * eps) ~ 1.4e-14 per entry.
        # Threshold 1e-8 gives ~700,000x margin, appropriate for round-trip between two methods.
        max_eigvec_err = 0.0
        for k in 1:3  # Spot-check first 3
            r_k = result.right_eigenvectors[:, k]
            lam_k = result.eigenvalues[k]
            err_k = norm(L_dense * r_k - lam_k * r_k)
            max_eigvec_err = max(max_eigvec_err, err_k)
            @test isapprox(L_dense * r_k, lam_k * r_k; atol=1e-8)
        end
        @info "DIAG-01: eigenvector equation" max_residual=max_eigvec_err threshold_atol=1e-8 modes_checked=3
    end

    # Need eigendata for subsequent tests
    eigen_result = extract_leading_eigendata(L_dense; n_modes=10)

    # -----------------------------------------------------------------------
    # DIAG-02: Fixed point distance
    # -----------------------------------------------------------------------
    @testset "DIAG-02: compute_fixed_point_distance" begin
        fp = compute_fixed_point_distance(eigen_result, N3_GIBBS)

        @test fp isa FixedPointResult

        # Fixed point trace distance should be small for BohrDomain with KMS construction.
        # KMS gives exact detailed balance, so the Lindbladian fixed point equals the Gibbs state
        # up to Gaussian filter smoothing in the frequency domain. For 3-qubit, the smoothing
        # effect is small, giving trace_distance ~ O(1e-3). Threshold 0.01 gives ~10x margin.
        @test fp.trace_distance < 0.01
        @info "DIAG-02: Bohr fixed point trace distance" trace_distance=fp.trace_distance threshold=0.01

        # Fixed point is normalized: tr(rho) = 1.
        # Density matrix reconstruction from eigenvector; trace error O(DIM * eps) ~ 8 * 2.2e-16 ~ 1.8e-15.
        # Threshold 1e-12 gives ~550x margin.
        fp_tr = real(tr(fp.fixed_point))
        @test isapprox(tr(fp.fixed_point), 1.0; atol=1e-12)
        @info "DIAG-02: fixed point trace" trace=fp_tr deviation=abs(fp_tr - 1.0) threshold_atol=1e-12

        # Fixed point is Hermitian: ||rho - rho'|| < atol.
        # Eigenvector-based reconstruction preserves Hermiticity to machine precision.
        # Error O(DIM^2 * eps) ~ 64 * 2.2e-16 ~ 1.4e-14. Threshold 1e-12 gives ~70x margin.
        herm_err = maximum(abs.(fp.fixed_point - fp.fixed_point'))
        @test isapprox(fp.fixed_point, fp.fixed_point'; atol=1e-12)
        @info "DIAG-02: fixed point Hermiticity" max_error=herm_err threshold_atol=1e-12

        # Fixed point has correct dimension
        @test size(fp.fixed_point) == (N3_DIM, N3_DIM)

        # Fixed point eigenvalues are non-negative (it's a valid density matrix).
        # Eigenvalues of a Hermitian positive-semidefinite matrix; numerical error can make
        # smallest eigenvalue slightly negative. Threshold -1e-12 allows for FP rounding.
        fp_eigvals = eigvals(Hermitian(fp.fixed_point))
        min_eigval = minimum(fp_eigvals)
        @test all(v -> v >= -1e-12, fp_eigvals)
        @info "DIAG-02: fixed point positivity" min_eigenvalue=min_eigval threshold=-1e-12
    end

    # -----------------------------------------------------------------------
    # DIAG-03/04: Anti-Hermitian defect
    # -----------------------------------------------------------------------
    @testset "DIAG-03/04: compute_anti_hermitian_defect" begin
        defect = compute_anti_hermitian_defect(L_dense, N3_GIBBS)

        @test defect isa DefectResult

        # A_norm is non-negative (Frobenius norm of anti-Hermitian part)
        @test defect.A_norm >= 0.0
        @info "DIAG-03: anti-Hermitian norm" A_norm=defect.A_norm

        # H_gap is positive (non-trivial Hermitian part has a gap)
        @test defect.H_gap > 0.0
        @info "DIAG-04: Hermitian gap" H_gap=defect.H_gap

        # Consistency check: defect_ratio = A_norm / H_gap.
        # This is a pure arithmetic identity (division). Error is O(eps) relative.
        # Threshold 1e-14 is well above machine epsilon 2.2e-16.
        ratio_err = abs(defect.defect_ratio - defect.A_norm / defect.H_gap)
        @test isapprox(defect.defect_ratio, defect.A_norm / defect.H_gap; atol=1e-14)
        @info "DIAG-03/04: defect ratio consistency" defect_ratio=defect.defect_ratio error=ratio_err threshold_atol=1e-14

        # Threshold is 0.1
        @test defect.threshold == 0.1

        # Warning is a Bool
        @test defect.warning isa Bool

        # Warning matches threshold comparison
        @test defect.warning == (defect.defect_ratio > defect.threshold)
    end

    # -----------------------------------------------------------------------
    # DIAG-05: Overlap coefficients
    # -----------------------------------------------------------------------
    @testset "DIAG-05: compute_overlap_coefficients" begin
        # Build observables in eigenbasis
        V = N3_HAM.eigvecs
        n_qubits = 3

        # Z1 in eigenbasis
        Z1_comp = Matrix{ComplexF64}(pad_term([Z], n_qubits, 1))
        Z1_eigen = Matrix{ComplexF64}(V' * Z1_comp * V)

        # H diagonal in eigenbasis
        H_eigen = Matrix{ComplexF64}(diagm(ComplexF64.(N3_HAM.eigvals)))

        observables = Matrix{ComplexF64}[Z1_eigen, H_eigen]
        observable_names = String["Z1", "H"]

        # Initial state: maximally mixed (same in any basis)
        dim = N3_DIM
        rho0 = Matrix{ComplexF64}(I(dim) / dim)

        overlap = compute_overlap_coefficients(
            eigen_result, observables, observable_names, rho0, N3_GIBBS;
            n_modes=10, initial_state_name="maximally_mixed"
        )

        @test overlap isa OverlapResult

        # Coefficient matrix dimensions
        @test size(overlap.coefficients) == (2, 10)

        # Observable names preserved
        @test overlap.observable_names == ["Z1", "H"]
        @test overlap.initial_state_name == "maximally_mixed"

        # Steady-state coefficient c_1 should be near zero (rho_beta subtracted).
        # By construction, the overlap with mode k=1 (steady state) is Tr[O * (rho0 - rho_beta)],
        # which is numerically zero since rho_beta was subtracted. Error O(DIM * eps) ~ 1.8e-15.
        # Threshold 1e-8 gives ~5,000,000x margin.
        max_c1_mixed = 0.0
        for i in 1:2
            max_c1_mixed = max(max_c1_mixed, abs(overlap.coefficients[i, 1]))
            @test abs(overlap.coefficients[i, 1]) < 1e-8
        end
        @info "DIAG-05: steady-state overlap (maximally_mixed)" max_abs_c1=max_c1_mixed threshold=1e-8

        # Gap mode overlap vector
        @test length(overlap.gap_mode_overlap) == 2
        @test all(g -> g >= 0.0, overlap.gap_mode_overlap)

        # Test with |0>^n initial state (transform to eigenbasis)
        psi0_comp = zeros(ComplexF64, dim)
        psi0_comp[1] = 1.0
        psi0_eigen = V' * psi0_comp
        rho_up = psi0_eigen * psi0_eigen'

        overlap_up = compute_overlap_coefficients(
            eigen_result, observables, observable_names, rho_up, N3_GIBBS;
            n_modes=10, initial_state_name="all_up"
        )
        @test size(overlap_up.coefficients) == (2, 10)
        # c_1 still near zero (rho_beta subtracted from any initial state).
        # Same reasoning as above.
        max_c1_up = 0.0
        for i in 1:2
            max_c1_up = max(max_c1_up, abs(overlap_up.coefficients[i, 1]))
            @test abs(overlap_up.coefficients[i, 1]) < 1e-8
        end
        @info "DIAG-05: steady-state overlap (all_up)" max_abs_c1=max_c1_up threshold=1e-8
    end

    # -----------------------------------------------------------------------
    # DIAG-06: Symmetry labels
    # -----------------------------------------------------------------------
    @testset "DIAG-06: compute_sz_labels" begin
        labels = compute_sz_labels(eigen_result, N3_HAM; n_modes=10)

        @test length(labels) == 10

        # Purity bounds: 0 <= purity <= 1 (with FP tolerance).
        # Purity is sum of squared weights, bounded by construction.
        # Threshold 1e-12 on upper bound allows for FP accumulation in weight computation.
        max_purity = 0.0
        for label in labels
            @test label isa SzSectorLabel
            @test 0.0 <= label.purity <= 1.0 + 1e-12
            @test label.is_pure == (label.purity > 0.95)
            @test !isempty(label.sector_weights)
            max_purity = max(max_purity, label.purity)
        end
        @info "DIAG-06: Sz label purities" max_purity=max_purity n_labels=length(labels)

        # Steady-state mode (k=1) should have delta_sz = 0
        # (fixed point is diagonal => Sz(i) - Sz(j) = 0 for nonzero entries)
        @test labels[1].delta_sz == 0.0
        # Should be pure in delta_sz=0 sector (purity > 0.95 by definition of is_pure)
        @test labels[1].purity > 0.95
        @info "DIAG-06: steady-state Sz" delta_sz=labels[1].delta_sz purity=labels[1].purity
    end

    # -----------------------------------------------------------------------
    # Multiplet detection
    # -----------------------------------------------------------------------
    @testset "detect_multiplets" begin
        multiplets = detect_multiplets(eigen_result.eigenvalues)

        @test !isempty(multiplets)

        # All eigenvalue indices should be covered
        all_indices = vcat([g.eigenvalue_indices for g in multiplets]...)
        @test sort(all_indices) == 1:10

        # Within each multiplet, eigenvalues are close.
        # Relative tolerance 0.01 (1%) is the default grouping criterion.
        # Within a true multiplet, eigenvalues differ by O(DIM * eps) ~ 1e-14 relative,
        # so 1% is extremely generous. The threshold is set by the detection algorithm's default.
        max_rel_err = 0.0
        for group in multiplets
            if length(group.eigenvalue_indices) > 1
                vals = [eigen_result.eigenvalues[i] for i in group.eigenvalue_indices]
                for (a, b) in Iterators.product(vals, vals)
                    denom = max(abs(a), abs(b), 1e-10)
                    rel_err = abs(a - b) / denom
                    max_rel_err = max(max_rel_err, rel_err)
                    @test rel_err < 0.01
                end
            end
        end
        @info "Multiplet relative spread" max_relative_error=max_rel_err threshold=0.01 n_multiplets=length(multiplets)

        # Mean eigenvalue is reasonable.
        # Arithmetic mean of complex eigenvalues; error O(eps) relative.
        # Threshold 1e-12 is well above machine epsilon.
        max_mean_err = 0.0
        for group in multiplets
            vals = [eigen_result.eigenvalues[i] for i in group.eigenvalue_indices]
            err = abs(group.mean_eigenvalue - sum(vals) / length(vals))
            max_mean_err = max(max_mean_err, err)
            @test isapprox(group.mean_eigenvalue, sum(vals) / length(vals); atol=1e-12)
        end
        @info "Multiplet mean eigenvalue" max_error=max_mean_err threshold_atol=1e-12

        # Edge case: empty eigenvalues
        empty_result = detect_multiplets(ComplexF64[])
        @test isempty(empty_result)

        # Edge case: single eigenvalue
        single_result = detect_multiplets(ComplexF64[1.0 + 0.0im])
        @test length(single_result) == 1
        @test single_result[1].eigenvalue_indices == [1]
    end

    # -----------------------------------------------------------------------
    # Bundle: run_exact_diagnostics
    # -----------------------------------------------------------------------
    @testset "run_exact_diagnostics bundle" begin
        result = run_exact_diagnostics(
            L_dense, N3_HAM, N3_GIBBS;
            n_modes=10
        )

        @test result isa ExactDiagnosticsResult

        # Eigen sub-result
        @test result.eigen isa EigenDecompositionResult
        @test length(result.eigen.eigenvalues) == 10

        # Fixed point sub-result: trace distance < 0.01 (same rationale as DIAG-02 above)
        @test result.fixed_point isa FixedPointResult
        @test result.fixed_point.trace_distance < 0.01
        @info "Bundle: fixed point trace distance" trace_distance=result.fixed_point.trace_distance threshold=0.01

        # Defect sub-result
        @test result.defect isa DefectResult
        @test result.defect.H_gap > 0.0
        @info "Bundle: Hermitian gap" H_gap=result.defect.H_gap

        # Default 3 initial states: all_up, all_plus, maximally_mixed
        @test length(result.overlaps) == 3
        @test result.overlaps[1].initial_state_name == "all_up"
        @test result.overlaps[2].initial_state_name == "all_plus"
        @test result.overlaps[3].initial_state_name == "maximally_mixed"
        # Steady-state overlap coefficient c_1 near zero for all initial states and observables.
        # Same reasoning as DIAG-05: Tr[O * (rho0 - rho_beta)] ~ O(DIM * eps). Threshold 1e-8.
        max_c1_bundle = 0.0
        for ov in result.overlaps
            @test ov isa OverlapResult
            @test size(ov.coefficients, 2) == 10
            for i in 1:size(ov.coefficients, 1)
                max_c1_bundle = max(max_c1_bundle, abs(ov.coefficients[i, 1]))
                @test abs(ov.coefficients[i, 1]) < 1e-8
            end
        end
        @info "Bundle: steady-state overlap (all states)" max_abs_c1=max_c1_bundle threshold=1e-8

        # Sz labels
        @test length(result.sz_labels) == 10
        @test result.sz_labels[1].delta_sz == 0.0

        # Multiplets
        @test !isempty(result.multiplets)
        all_idx = sort(vcat([g.eigenvalue_indices for g in result.multiplets]...))
        @test all_idx == 1:10
    end

    # -----------------------------------------------------------------------
    # Bundle: custom observables and initial states
    # -----------------------------------------------------------------------
    @testset "run_exact_diagnostics custom inputs" begin
        V = N3_HAM.eigvecs
        n_qubits = 3
        dim = N3_DIM

        # Custom observable: just Z1
        Z1_comp = Matrix{ComplexF64}(pad_term([Z], n_qubits, 1))
        Z1_eigen = Matrix{ComplexF64}(V' * Z1_comp * V)

        # Custom initial state: maximally mixed
        rho_mixed = Matrix{ComplexF64}(I(dim) / dim)

        result = run_exact_diagnostics(
            L_dense, N3_HAM, N3_GIBBS;
            observables=Matrix{ComplexF64}[Z1_eigen],
            observable_names=String["Z1"],
            initial_states=Matrix{ComplexF64}[rho_mixed],
            initial_state_names=String["mixed"],
            n_modes=10
        )

        @test length(result.overlaps) == 1
        @test result.overlaps[1].initial_state_name == "mixed"
        @test size(result.overlaps[1].coefficients) == (1, 10)
    end

end

# ==========================================================================
# TrotterDomain Diagnostics
# ==========================================================================
@testset "Diagnostics TrotterDomain" begin

    # Build TrotterDomain Lindbladian (3-qubit)
    config_trott = make_config(Lindbladian(), TrotterDomain(); num_qubits=3, construction=KMS())
    L_trott_sparse = construct_lindbladian(N3_TROTTER_JUMPS, config_trott, N3_HAM; trotter=N3_TROTTER)
    L_trott = Matrix{ComplexF64}(L_trott_sparse)

    # Gibbs state in Trotter basis (matching furnace.jl pattern)
    gibbs_trott = Hermitian(N3_TROTTER.eigvecs' * N3_HAM.eigvecs *
                            Matrix(N3_GIBBS) * N3_HAM.eigvecs' * N3_TROTTER.eigvecs)

    n_qubits = 3
    dim = N3_DIM

    # -------------------------------------------------------------------
    # DIAG-01: Eigendata extraction on TrotterDomain L
    # -------------------------------------------------------------------
    @testset "DIAG-01: TrotterDomain eigendata" begin
        result = extract_leading_eigendata(L_trott; n_modes=10)

        @test result isa EigenDecompositionResult
        @test length(result.eigenvalues) == 10

        # Eigenvalues sorted by |Re(lambda)| (ascending)
        @test issorted(abs.(real.(result.eigenvalues)))

        # First eigenvalue near zero (steady state).
        # Same reasoning as BohrDomain DIAG-01: dense eigen on DIM^2=64 matrix,
        # error O(DIM^2 * eps) ~ 1.4e-14. Threshold 1e-10 gives ~7000x margin.
        @test abs(result.eigenvalues[1]) < 1e-10
        @info "DIAG-01 Trotter: steady-state eigenvalue" abs_lambda1=abs(result.eigenvalues[1]) threshold=1e-10

        # Spectral gap is positive
        @test result.spectral_gap > 0.0
        @info "DIAG-01 Trotter: spectral gap" gap=result.spectral_gap

        # Biorthonormality: V_left' * V_right approx I
        # Same error analysis as BohrDomain: O(DIM^2 * eps * sqrt(n_modes)) ~ 5e-14.
        # Threshold 1e-8 gives ~200,000x margin.
        biorth = result.left_eigenvectors' * result.right_eigenvectors
        biorth_err = maximum(abs.(biorth - I(10)))
        @test isapprox(biorth, I(10); atol=1e-8)
        @info "DIAG-01 Trotter: biorthonormality" max_error=biorth_err threshold_atol=1e-8

        # Eigenvectors satisfy eigenvalue equation.
        # Same error analysis as BohrDomain: O(DIM^2 * eps) per entry. Threshold 1e-8.
        max_eigvec_err = 0.0
        for k in 1:3
            r_k = result.right_eigenvectors[:, k]
            lam_k = result.eigenvalues[k]
            err_k = norm(L_trott * r_k - lam_k * r_k)
            max_eigvec_err = max(max_eigvec_err, err_k)
            @test isapprox(L_trott * r_k, lam_k * r_k; atol=1e-8)
        end
        @info "DIAG-01 Trotter: eigenvector equation" max_residual=max_eigvec_err threshold_atol=1e-8 modes_checked=3
    end

    eigen_trott = extract_leading_eigendata(L_trott; n_modes=10)

    # -------------------------------------------------------------------
    # DIAG-02: Fixed point distance with Trotter-basis Gibbs
    # -------------------------------------------------------------------
    @testset "DIAG-02: TrotterDomain fixed point distance" begin
        fp = compute_fixed_point_distance(eigen_trott, gibbs_trott)

        @test fp isa FixedPointResult

        # TrotterDomain fixed point distance should be NON-TRIVIALLY larger
        # than BohrDomain due to Trotter error shifting the fixed point.
        # Trotter error is O(dt^2) where dt = T0/n_steps. For 3-qubit with 10 Trotter steps,
        # the distance is small but nonzero. Threshold 0.5 is very conservative upper bound.
        @test fp.trace_distance > 0.0
        @test fp.trace_distance < 0.5
        @info "DIAG-02 Trotter: fixed point trace distance" trace_distance=fp.trace_distance upper_bound=0.5

        # Fixed point is normalized: same error analysis as BohrDomain DIAG-02.
        fp_tr = real(tr(fp.fixed_point))
        @test isapprox(tr(fp.fixed_point), 1.0; atol=1e-12)
        @info "DIAG-02 Trotter: fixed point trace" trace=fp_tr deviation=abs(fp_tr - 1.0) threshold_atol=1e-12

        # Fixed point is Hermitian: same error analysis as BohrDomain DIAG-02.
        herm_err = maximum(abs.(fp.fixed_point - fp.fixed_point'))
        @test isapprox(fp.fixed_point, fp.fixed_point'; atol=1e-12)
        @info "DIAG-02 Trotter: fixed point Hermiticity" max_error=herm_err threshold_atol=1e-12

        # Fixed point is a valid density matrix (non-negative eigenvalues)
        fp_eigvals = eigvals(Hermitian(fp.fixed_point))
        min_eigval = minimum(fp_eigvals)
        @test all(v -> v >= -1e-12, fp_eigvals)
        @info "DIAG-02 Trotter: fixed point positivity" min_eigenvalue=min_eigval threshold=-1e-12

        # Compare to BohrDomain fixed point distance -- both should be similar magnitude
        # (for 3-qubit with 10 Trotter steps, Trotter error is very small so distances
        # are nearly equal; we just verify they're in the same ballpark)
        config_bohr = make_config(Lindbladian(), BohrDomain(); num_qubits=3, construction=KMS())
        L_bohr_sparse = construct_lindbladian(N3_JUMPS, config_bohr, N3_HAM)
        L_bohr = Matrix{ComplexF64}(L_bohr_sparse)
        eigen_bohr = extract_leading_eigendata(L_bohr; n_modes=10)
        fp_bohr = compute_fixed_point_distance(eigen_bohr, N3_GIBBS)
        # Both distances are small with KMS (exact detailed balance).
        # Trotter has tiny residual Trotter error; Bohr is near machine precision.
        # Threshold 0.01 is conservative for both.
        @test fp.trace_distance < 0.01
        @test fp_bohr.trace_distance < 0.01
        @info "DIAG-02: Trotter vs Bohr comparison" trotter_dist=fp.trace_distance bohr_dist=fp_bohr.trace_distance threshold=0.01
    end

    # -------------------------------------------------------------------
    # DIAG-03/04: Anti-Hermitian defect with Trotter-basis Gibbs
    # -------------------------------------------------------------------
    @testset "DIAG-03/04: TrotterDomain anti-Hermitian defect" begin
        defect = compute_anti_hermitian_defect(L_trott, gibbs_trott)

        @test defect isa DefectResult

        # A_norm is non-negative (Frobenius norm)
        @test defect.A_norm >= 0.0
        @info "DIAG-03 Trotter: anti-Hermitian norm" A_norm=defect.A_norm

        # H_gap is positive
        @test defect.H_gap > 0.0
        @info "DIAG-04 Trotter: Hermitian gap" H_gap=defect.H_gap

        # Consistency: defect_ratio = A_norm / H_gap.
        # Pure arithmetic identity. Threshold 1e-14 well above eps.
        ratio_err = abs(defect.defect_ratio - defect.A_norm / defect.H_gap)
        @test isapprox(defect.defect_ratio, defect.A_norm / defect.H_gap; atol=1e-14)
        @info "DIAG-03/04 Trotter: defect ratio" defect_ratio=defect.defect_ratio error=ratio_err threshold_atol=1e-14

        # Threshold is 0.1
        @test defect.threshold == 0.1
        @test defect.warning == (defect.defect_ratio > defect.threshold)
    end

    # -------------------------------------------------------------------
    # DIAG-05: Overlap coefficients with Trotter-basis observables
    # -------------------------------------------------------------------
    @testset "DIAG-05: TrotterDomain overlap coefficients" begin
        Vt = N3_TROTTER.eigvecs

        # Z1 in Trotter basis
        Z1_comp = Matrix{ComplexF64}(pad_term([Z], n_qubits, 1))
        Z1_trott = Matrix{ComplexF64}(Vt' * Z1_comp * Vt)

        # I/dim as initial state (same in any basis)
        rho0 = Matrix{ComplexF64}(I(dim) / dim)

        overlap = compute_overlap_coefficients(
            eigen_trott, [Z1_trott], ["Z1"], rho0, gibbs_trott;
            n_modes=10, initial_state_name="maximally_mixed"
        )

        @test overlap isa OverlapResult
        @test size(overlap.coefficients) == (1, 10)
        @test overlap.observable_names == ["Z1"]

        # c_1 near zero (steady-state mode, rho_beta subtracted).
        # Same reasoning as BohrDomain DIAG-05. Threshold 1e-8.
        @test abs(overlap.coefficients[1, 1]) < 1e-8
        @info "DIAG-05 Trotter: steady-state overlap" abs_c1=abs(overlap.coefficients[1, 1]) threshold=1e-8

        # Gap mode overlap is non-negative
        @test length(overlap.gap_mode_overlap) == 1
        @test overlap.gap_mode_overlap[1] >= 0.0
    end

    # -------------------------------------------------------------------
    # DIAG-06: Sz labels with Trotter eigenvectors
    # -------------------------------------------------------------------
    @testset "DIAG-06: TrotterDomain Sz labels" begin
        labels = compute_sz_labels(eigen_trott, N3_TROTTER.eigvecs, n_qubits; n_modes=10)

        @test length(labels) == 10

        # Purity bounds with FP tolerance (same reasoning as BohrDomain DIAG-06)
        max_purity = 0.0
        for label in labels
            @test label isa SzSectorLabel
            @test 0.0 <= label.purity <= 1.0 + 1e-12
            @test label.is_pure == (label.purity > 0.95)
            @test !isempty(label.sector_weights)
            max_purity = max(max_purity, label.purity)
        end
        @info "DIAG-06 Trotter: Sz label purities" max_purity=max_purity n_labels=length(labels)

        # Steady-state mode (k=1) should have delta_sz ~ 0
        @test labels[1].delta_sz == 0.0
    end

    # -------------------------------------------------------------------
    # Bundle: run_exact_diagnostics with basis_eigvecs
    # -------------------------------------------------------------------
    @testset "run_exact_diagnostics TrotterDomain bundle" begin
        result = run_exact_diagnostics(L_trott, N3_HAM, gibbs_trott;
            n_modes=10, basis_eigvecs=N3_TROTTER.eigvecs)

        @test result isa ExactDiagnosticsResult

        # Eigen sub-result
        @test result.eigen isa EigenDecompositionResult
        @test length(result.eigen.eigenvalues) == 10

        # Fixed point sub-result: trace distance < 0.5 (Trotter has larger error than Bohr).
        # For 3-qubit with 10 Trotter steps, actual distance is typically < 0.01,
        # but 0.5 gives margin for different Trotter step counts.
        @test result.fixed_point isa FixedPointResult
        @test result.fixed_point.trace_distance < 0.5
        @info "Trotter bundle: fixed point trace distance" trace_distance=result.fixed_point.trace_distance threshold=0.5

        # Defect sub-result
        @test result.defect isa DefectResult
        @test result.defect.H_gap > 0.0
        @info "Trotter bundle: Hermitian gap" H_gap=result.defect.H_gap

        # Default 3 initial states
        @test length(result.overlaps) == 3
        @test result.overlaps[1].initial_state_name == "all_up"
        @test result.overlaps[2].initial_state_name == "all_plus"
        @test result.overlaps[3].initial_state_name == "maximally_mixed"
        # Steady-state overlap c_1 near zero for all states. Threshold 1e-8.
        max_c1_trott_bundle = 0.0
        for ov in result.overlaps
            @test ov isa OverlapResult
            @test size(ov.coefficients, 2) == 10
            for i in 1:size(ov.coefficients, 1)
                max_c1_trott_bundle = max(max_c1_trott_bundle, abs(ov.coefficients[i, 1]))
                @test abs(ov.coefficients[i, 1]) < 1e-8
            end
        end
        @info "Trotter bundle: steady-state overlap (all states)" max_abs_c1=max_c1_trott_bundle threshold=1e-8

        # Sz labels use Trotter eigenvectors
        @test length(result.sz_labels) == 10
        @test result.sz_labels[1].delta_sz == 0.0

        # Multiplets
        @test !isempty(result.multiplets)
        all_idx = sort(vcat([g.eigenvalue_indices for g in result.multiplets]...))
        @test all_idx == 1:10
    end

    # -------------------------------------------------------------------
    # Backward compatibility: BohrDomain without basis_eigvecs
    # -------------------------------------------------------------------
    @testset "backward compatibility: no basis_eigvecs" begin
        config_bohr = make_config(Lindbladian(), BohrDomain(); num_qubits=3, construction=KMS())
        L_bohr = Matrix{ComplexF64}(construct_lindbladian(N3_JUMPS, config_bohr, N3_HAM))

        result = run_exact_diagnostics(L_bohr, N3_HAM, N3_GIBBS; n_modes=10)

        @test result isa ExactDiagnosticsResult
        # Trace distance < 0.01 (BohrDomain + KMS = near-exact Gibbs fixed point)
        @test result.fixed_point.trace_distance < 0.01
        @info "Backward compat: Bohr trace distance" trace_distance=result.fixed_point.trace_distance threshold=0.01
        @test length(result.overlaps) == 3
        @test result.sz_labels[1].delta_sz == 0.0
    end

end
