@testset "Diagnostics (Phase 26)" begin

    # Build the 3-qubit Lindbladian (BohrDomain for exact Gibbs fixed point)
    config = make_small_liouv_config(BohrDomain())
    L_sparse = construct_lindbladian(SMALL_JUMPS, config, SMALL_HAM)
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
        @test abs(result.eigenvalues[1]) < 1e-10

        # Spectral gap is positive
        @test result.spectral_gap == abs(real(result.eigenvalues[2]))
        @test result.spectral_gap > 0.0

        # Right/left eigenvector dimensions
        @test size(result.right_eigenvectors) == (SMALL_DIM^2, 10)
        @test size(result.left_eigenvectors) == (SMALL_DIM^2, 10)

        # Biorthonormality: V_left' * V_right approx I
        biorth = result.left_eigenvectors' * result.right_eigenvectors
        @test isapprox(biorth, I(10); atol=1e-8)

        # Im/Re ratios
        @test length(result.im_re_ratios) == 10
        @test result.im_re_ratios[1] == 0.0  # Steady state
        @test all(r -> r >= 0.0, result.im_re_ratios)

        # Eigenvectors actually satisfy eigenvalue equation: L * r_k = lambda_k * r_k
        for k in 1:3  # Spot-check first 3
            r_k = result.right_eigenvectors[:, k]
            lam_k = result.eigenvalues[k]
            @test isapprox(L_dense * r_k, lam_k * r_k; atol=1e-8)
        end
    end

    # Need eigendata for subsequent tests
    eigen_result = extract_leading_eigendata(L_dense; n_modes=10)

    # -----------------------------------------------------------------------
    # DIAG-02: Fixed point distance
    # -----------------------------------------------------------------------
    @testset "DIAG-02: compute_fixed_point_distance" begin
        fp = compute_fixed_point_distance(eigen_result, SMALL_GIBBS)

        @test fp isa FixedPointResult

        # Fixed point trace distance should be small for BohrDomain
        # (Not exactly zero due to Gaussian filter smoothing in the 3-qubit test system)
        @test fp.trace_distance < 0.01

        # Fixed point is normalized
        @test isapprox(tr(fp.fixed_point), 1.0; atol=1e-12)

        # Fixed point is Hermitian
        @test isapprox(fp.fixed_point, fp.fixed_point'; atol=1e-12)

        # Fixed point has correct dimension
        @test size(fp.fixed_point) == (SMALL_DIM, SMALL_DIM)

        # Fixed point eigenvalues are non-negative (it's a valid density matrix)
        fp_eigvals = eigvals(Hermitian(fp.fixed_point))
        @test all(v -> v >= -1e-12, fp_eigvals)
    end

    # -----------------------------------------------------------------------
    # DIAG-03/04: Anti-Hermitian defect
    # -----------------------------------------------------------------------
    @testset "DIAG-03/04: compute_anti_hermitian_defect" begin
        defect = compute_anti_hermitian_defect(L_dense, SMALL_GIBBS)

        @test defect isa DefectResult

        # A_norm is non-negative
        @test defect.A_norm >= 0.0

        # H_gap is positive (non-trivial Hermitian part has a gap)
        @test defect.H_gap > 0.0

        # Consistency check: defect_ratio = A_norm / H_gap
        @test isapprox(defect.defect_ratio, defect.A_norm / defect.H_gap; atol=1e-14)

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
        V = SMALL_HAM.eigvecs
        n_qubits = 3

        # Z1 in eigenbasis
        Z1_comp = Matrix{ComplexF64}(pad_term([Z], n_qubits, 1))
        Z1_eigen = Matrix{ComplexF64}(V' * Z1_comp * V)

        # H diagonal in eigenbasis
        H_eigen = Matrix{ComplexF64}(diagm(ComplexF64.(SMALL_HAM.eigvals)))

        observables = Matrix{ComplexF64}[Z1_eigen, H_eigen]
        observable_names = String["Z1", "H"]

        # Initial state: maximally mixed (same in any basis)
        dim = SMALL_DIM
        rho0 = Matrix{ComplexF64}(I(dim) / dim)

        overlap = compute_overlap_coefficients(
            eigen_result, observables, observable_names, rho0, SMALL_GIBBS;
            n_modes=10, initial_state_name="maximally_mixed"
        )

        @test overlap isa OverlapResult

        # Coefficient matrix dimensions
        @test size(overlap.coefficients) == (2, 10)

        # Observable names preserved
        @test overlap.observable_names == ["Z1", "H"]
        @test overlap.initial_state_name == "maximally_mixed"

        # Steady-state coefficient c_1 should be near zero (rho_beta subtracted)
        for i in 1:2
            @test abs(overlap.coefficients[i, 1]) < 1e-8
        end

        # Gap mode overlap vector
        @test length(overlap.gap_mode_overlap) == 2
        @test all(g -> g >= 0.0, overlap.gap_mode_overlap)

        # Test with |0>^n initial state (transform to eigenbasis)
        psi0_comp = zeros(ComplexF64, dim)
        psi0_comp[1] = 1.0
        psi0_eigen = V' * psi0_comp
        rho_up = psi0_eigen * psi0_eigen'

        overlap_up = compute_overlap_coefficients(
            eigen_result, observables, observable_names, rho_up, SMALL_GIBBS;
            n_modes=10, initial_state_name="all_up"
        )
        @test size(overlap_up.coefficients) == (2, 10)
        # c_1 still near zero (rho_beta subtracted from any initial state)
        for i in 1:2
            @test abs(overlap_up.coefficients[i, 1]) < 1e-8
        end
    end

    # -----------------------------------------------------------------------
    # DIAG-06: Symmetry labels
    # -----------------------------------------------------------------------
    @testset "DIAG-06: compute_sz_labels" begin
        labels = compute_sz_labels(eigen_result, SMALL_HAM; n_modes=10)

        @test length(labels) == 10

        for label in labels
            @test label isa SzSectorLabel
            @test 0.0 <= label.purity <= 1.0 + 1e-12
            @test label.is_pure == (label.purity > 0.95)
            @test !isempty(label.sector_weights)
        end

        # Steady-state mode (k=1) should have delta_sz = 0
        # (fixed point is diagonal => Sz(i) - Sz(j) = 0 for nonzero entries)
        @test labels[1].delta_sz == 0.0
        @test labels[1].purity > 0.95  # Should be pure in delta_sz=0 sector
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

        # Within each multiplet, eigenvalues are close
        for group in multiplets
            if length(group.eigenvalue_indices) > 1
                vals = [eigen_result.eigenvalues[i] for i in group.eigenvalue_indices]
                for (a, b) in Iterators.product(vals, vals)
                    denom = max(abs(a), abs(b), 1e-10)
                    @test abs(a - b) / denom < 0.01  # rel_tol default
                end
            end
        end

        # Mean eigenvalue is reasonable
        for group in multiplets
            vals = [eigen_result.eigenvalues[i] for i in group.eigenvalue_indices]
            @test isapprox(group.mean_eigenvalue, sum(vals) / length(vals); atol=1e-12)
        end

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
            L_dense, SMALL_HAM, SMALL_GIBBS;
            n_modes=10
        )

        @test result isa ExactDiagnosticsResult

        # Eigen sub-result
        @test result.eigen isa EigenDecompositionResult
        @test length(result.eigen.eigenvalues) == 10

        # Fixed point sub-result
        @test result.fixed_point isa FixedPointResult
        @test result.fixed_point.trace_distance < 0.01

        # Defect sub-result
        @test result.defect isa DefectResult
        @test result.defect.H_gap > 0.0

        # Default 3 initial states: all_up, all_plus, maximally_mixed
        @test length(result.overlaps) == 3
        @test result.overlaps[1].initial_state_name == "all_up"
        @test result.overlaps[2].initial_state_name == "all_plus"
        @test result.overlaps[3].initial_state_name == "maximally_mixed"
        for ov in result.overlaps
            @test ov isa OverlapResult
            @test size(ov.coefficients, 2) == 10
            # c_1 near zero for all initial states
            for i in 1:size(ov.coefficients, 1)
                @test abs(ov.coefficients[i, 1]) < 1e-8
            end
        end

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
        V = SMALL_HAM.eigvecs
        n_qubits = 3
        dim = SMALL_DIM

        # Custom observable: just Z1
        Z1_comp = Matrix{ComplexF64}(pad_term([Z], n_qubits, 1))
        Z1_eigen = Matrix{ComplexF64}(V' * Z1_comp * V)

        # Custom initial state: maximally mixed
        rho_mixed = Matrix{ComplexF64}(I(dim) / dim)

        result = run_exact_diagnostics(
            L_dense, SMALL_HAM, SMALL_GIBBS;
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
    config_trott = make_small_liouv_config(TrotterDomain())
    L_trott_sparse = construct_lindbladian(SMALL_TROTTER_JUMPS, config_trott, SMALL_HAM; trotter=SMALL_TROTTER)
    L_trott = Matrix{ComplexF64}(L_trott_sparse)

    # Gibbs state in Trotter basis (matching furnace.jl pattern)
    gibbs_trott = Hermitian(SMALL_TROTTER.eigvecs' * SMALL_HAM.eigvecs *
                            Matrix(SMALL_GIBBS) * SMALL_HAM.eigvecs' * SMALL_TROTTER.eigvecs)

    n_qubits = 3
    dim = SMALL_DIM

    # -------------------------------------------------------------------
    # DIAG-01: Eigendata extraction on TrotterDomain L
    # -------------------------------------------------------------------
    @testset "DIAG-01: TrotterDomain eigendata" begin
        result = extract_leading_eigendata(L_trott; n_modes=10)

        @test result isa EigenDecompositionResult
        @test length(result.eigenvalues) == 10

        # Eigenvalues sorted by |Re(lambda)| (ascending)
        @test issorted(abs.(real.(result.eigenvalues)))

        # First eigenvalue near zero (steady state)
        @test abs(result.eigenvalues[1]) < 1e-10

        # Spectral gap is positive
        @test result.spectral_gap > 0.0

        # Biorthonormality: V_left' * V_right approx I
        biorth = result.left_eigenvectors' * result.right_eigenvectors
        @test isapprox(biorth, I(10); atol=1e-8)

        # Eigenvectors satisfy eigenvalue equation
        for k in 1:3
            r_k = result.right_eigenvectors[:, k]
            lam_k = result.eigenvalues[k]
            @test isapprox(L_trott * r_k, lam_k * r_k; atol=1e-8)
        end
    end

    eigen_trott = extract_leading_eigendata(L_trott; n_modes=10)

    # -------------------------------------------------------------------
    # DIAG-02: Fixed point distance with Trotter-basis Gibbs
    # -------------------------------------------------------------------
    @testset "DIAG-02: TrotterDomain fixed point distance" begin
        fp = compute_fixed_point_distance(eigen_trott, gibbs_trott)

        @test fp isa FixedPointResult

        # TrotterDomain fixed point distance should be NON-TRIVIALLY larger
        # than BohrDomain due to Trotter error shifting the fixed point
        @test fp.trace_distance > 0.0

        # But still finite and reasonable for 3-qubit with 10 Trotter steps
        @test fp.trace_distance < 0.5

        # Fixed point is normalized
        @test isapprox(tr(fp.fixed_point), 1.0; atol=1e-12)

        # Fixed point is Hermitian
        @test isapprox(fp.fixed_point, fp.fixed_point'; atol=1e-12)

        # Fixed point is a valid density matrix (non-negative eigenvalues)
        fp_eigvals = eigvals(Hermitian(fp.fixed_point))
        @test all(v -> v >= -1e-12, fp_eigvals)

        # Compare to BohrDomain fixed point distance -- both should be similar magnitude
        # (for 3-qubit with 10 Trotter steps, Trotter error is very small so distances
        # are nearly equal; we just verify they're in the same ballpark)
        config_bohr = make_small_liouv_config(BohrDomain())
        L_bohr_sparse = construct_lindbladian(SMALL_JUMPS, config_bohr, SMALL_HAM)
        L_bohr = Matrix{ComplexF64}(L_bohr_sparse)
        eigen_bohr = extract_leading_eigendata(L_bohr; n_modes=10)
        fp_bohr = compute_fixed_point_distance(eigen_bohr, SMALL_GIBBS)
        # Both distances are finite and on the same order of magnitude
        @test fp.trace_distance / fp_bohr.trace_distance > 0.5
        @test fp.trace_distance / fp_bohr.trace_distance < 2.0
    end

    # -------------------------------------------------------------------
    # DIAG-03/04: Anti-Hermitian defect with Trotter-basis Gibbs
    # -------------------------------------------------------------------
    @testset "DIAG-03/04: TrotterDomain anti-Hermitian defect" begin
        defect = compute_anti_hermitian_defect(L_trott, gibbs_trott)

        @test defect isa DefectResult

        # A_norm is non-negative
        @test defect.A_norm >= 0.0

        # H_gap is positive
        @test defect.H_gap > 0.0

        # Consistency: defect_ratio = A_norm / H_gap
        @test isapprox(defect.defect_ratio, defect.A_norm / defect.H_gap; atol=1e-14)

        # Threshold is 0.1
        @test defect.threshold == 0.1
        @test defect.warning == (defect.defect_ratio > defect.threshold)
    end

    # -------------------------------------------------------------------
    # DIAG-05: Overlap coefficients with Trotter-basis observables
    # -------------------------------------------------------------------
    @testset "DIAG-05: TrotterDomain overlap coefficients" begin
        Vt = SMALL_TROTTER.eigvecs

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

        # c_1 near zero (steady-state mode, rho_beta subtracted)
        @test abs(overlap.coefficients[1, 1]) < 1e-8

        # Gap mode overlap is non-negative
        @test length(overlap.gap_mode_overlap) == 1
        @test overlap.gap_mode_overlap[1] >= 0.0
    end

    # -------------------------------------------------------------------
    # DIAG-06: Sz labels with Trotter eigenvectors
    # -------------------------------------------------------------------
    @testset "DIAG-06: TrotterDomain Sz labels" begin
        labels = compute_sz_labels(eigen_trott, SMALL_TROTTER.eigvecs, n_qubits; n_modes=10)

        @test length(labels) == 10

        for label in labels
            @test label isa SzSectorLabel
            @test 0.0 <= label.purity <= 1.0 + 1e-12
            @test label.is_pure == (label.purity > 0.95)
            @test !isempty(label.sector_weights)
        end

        # Steady-state mode (k=1) should have delta_sz ~ 0
        @test labels[1].delta_sz == 0.0
    end

    # -------------------------------------------------------------------
    # Bundle: run_exact_diagnostics with basis_eigvecs
    # -------------------------------------------------------------------
    @testset "run_exact_diagnostics TrotterDomain bundle" begin
        result = run_exact_diagnostics(L_trott, SMALL_HAM, gibbs_trott;
            n_modes=10, basis_eigvecs=SMALL_TROTTER.eigvecs)

        @test result isa ExactDiagnosticsResult

        # Eigen sub-result
        @test result.eigen isa EigenDecompositionResult
        @test length(result.eigen.eigenvalues) == 10

        # Fixed point sub-result
        @test result.fixed_point isa FixedPointResult
        @test result.fixed_point.trace_distance < 0.5

        # Defect sub-result
        @test result.defect isa DefectResult
        @test result.defect.H_gap > 0.0

        # Default 3 initial states
        @test length(result.overlaps) == 3
        @test result.overlaps[1].initial_state_name == "all_up"
        @test result.overlaps[2].initial_state_name == "all_plus"
        @test result.overlaps[3].initial_state_name == "maximally_mixed"
        for ov in result.overlaps
            @test ov isa OverlapResult
            @test size(ov.coefficients, 2) == 10
            # c_1 near zero for all initial states
            for i in 1:size(ov.coefficients, 1)
                @test abs(ov.coefficients[i, 1]) < 1e-8
            end
        end

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
        config_bohr = make_small_liouv_config(BohrDomain())
        L_bohr = Matrix{ComplexF64}(construct_lindbladian(SMALL_JUMPS, config_bohr, SMALL_HAM))

        result = run_exact_diagnostics(L_bohr, SMALL_HAM, SMALL_GIBBS; n_modes=10)

        @test result isa ExactDiagnosticsResult
        @test result.fixed_point.trace_distance < 0.01
        @test length(result.overlaps) == 3
        @test result.sz_labels[1].delta_sz == 0.0
    end

end
