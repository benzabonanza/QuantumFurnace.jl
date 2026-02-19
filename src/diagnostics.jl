# ============================================================================
# Exact Diagnostics: ground-truth spectral data and structural analysis
# ============================================================================
#
# Provides six diagnostic functions (DIAG-01 through DIAG-06) plus a
# `run_exact_diagnostics()` bundle for validating trajectory-based gap
# estimates at n=4,6.  All computations use dense eigen() which is
# feasible for dim^2 <= 4096 (n <= 6 qubits).
#
# DIAG-01: extract_leading_eigendata      -- leading eigenvalues + left/right eigenvectors
# DIAG-02: compute_fixed_point_distance   -- fixed point vs Gibbs trace distance
# DIAG-03/04: compute_anti_hermitian_defect -- KMS similarity transform defect ratio
# DIAG-05: compute_overlap_coefficients   -- observable overlap with eigenmodes
# DIAG-06: compute_sz_labels              -- Delta_Sz symmetry sector labeling
#          detect_multiplets              -- near-degenerate eigenvalue grouping

# ---------------------------------------------------------------------------
# Result structs
# ---------------------------------------------------------------------------

"""
    EigenDecompositionResult

Result of DIAG-01: leading eigenvalue extraction with left and right eigenvectors.

# Fields
- `eigenvalues`: Leading n_modes eigenvalues sorted by |Re(lambda)|.
- `right_eigenvectors`: dim^2 x n_modes matrix (columns = right eigenvectors).
- `left_eigenvectors`: dim^2 x n_modes matrix (columns = left eigenvectors from inv(V)').
- `spectral_gap`: abs(real(eigenvalues[2])).
- `im_re_ratios`: |Im(lambda_k)/Re(lambda_k)| for each leading mode.
"""
struct EigenDecompositionResult
    eigenvalues::Vector{ComplexF64}
    right_eigenvectors::Matrix{ComplexF64}
    left_eigenvectors::Matrix{ComplexF64}
    spectral_gap::Float64
    im_re_ratios::Vector{Float64}
end

"""
    FixedPointResult

Result of DIAG-02: Lindbladian fixed point compared to Gibbs state.

# Fields
- `fixed_point`: Normalized density matrix from lambda_1 eigenvector.
- `trace_distance`: Trace distance to Gibbs state.
"""
struct FixedPointResult
    fixed_point::Matrix{ComplexF64}
    trace_distance::Float64
end

"""
    DefectResult

Result of DIAG-03/04: anti-Hermitian defect of the KMS-similarity-transformed Lindbladian.

# Fields
- `A_norm`: Operator 2-norm of anti-Hermitian part.
- `H_gap`: Gap of Hermitian part of similarity-transformed Lindbladian.
- `defect_ratio`: A_norm / H_gap.
- `warning`: true if defect_ratio > threshold (advisory only, does NOT gate anything).
- `threshold`: The warning threshold used.
"""
struct DefectResult
    A_norm::Float64
    H_gap::Float64
    defect_ratio::Float64
    warning::Bool
    threshold::Float64
end

"""
    OverlapResult

Result of DIAG-05: observable overlap coefficients with Lindbladian eigenmodes.

# Fields
- `coefficients`: n_obs x n_modes overlap coefficients.
- `observable_names`: Names of observables.
- `initial_state_name`: Name of initial state used.
- `gap_mode_overlap`: |c_2| per observable.
"""
struct OverlapResult
    coefficients::Matrix{ComplexF64}
    observable_names::Vector{String}
    initial_state_name::String
    gap_mode_overlap::Vector{Float64}
end

"""
    SzSectorLabel

Result of DIAG-06: Delta_Sz symmetry sector label for a single Lindbladian eigenvector.

# Fields
- `delta_sz`: Dominant Delta_Sz quantum number.
- `purity`: Fraction of weight in dominant sector.
- `is_pure`: purity > 0.95.
- `sector_weights`: All Delta_Sz weights.
"""
struct SzSectorLabel
    delta_sz::Float64
    purity::Float64
    is_pure::Bool
    sector_weights::Dict{Float64, Float64}
end

"""
    MultipletGroup

A group of near-degenerate Lindbladian eigenvalues forming a multiplet.

# Fields
- `eigenvalue_indices`: Indices into eigenvalue array.
- `mean_eigenvalue`: Mean eigenvalue of the group.
- `sz_labels`: Per-eigenvector SzSectorLabel (may be empty if not yet computed).
"""
struct MultipletGroup
    eigenvalue_indices::Vector{Int}
    mean_eigenvalue::ComplexF64
    sz_labels::Vector{SzSectorLabel}
end

"""
    ExactDiagnosticsResult

Bundle result from `run_exact_diagnostics()` containing all six DIAG outputs.

# Fields
- `eigen`: EigenDecompositionResult from DIAG-01.
- `fixed_point`: FixedPointResult from DIAG-02.
- `defect`: DefectResult from DIAG-03/04.
- `overlaps`: Vector{OverlapResult}, one per initial state, from DIAG-05.
- `sz_labels`: Vector{SzSectorLabel}, one per mode, from DIAG-06.
- `multiplets`: Vector{MultipletGroup}, grouped near-degenerate modes.
"""
struct ExactDiagnosticsResult
    eigen::EigenDecompositionResult
    fixed_point::FixedPointResult
    defect::DefectResult
    overlaps::Vector{OverlapResult}
    sz_labels::Vector{SzSectorLabel}
    multiplets::Vector{MultipletGroup}
end

# ---------------------------------------------------------------------------
# DIAG-01: Leading eigenvalue extraction
# ---------------------------------------------------------------------------

"""
    extract_leading_eigendata(L::Matrix{ComplexF64}; n_modes::Int=20) -> EigenDecompositionResult

Extract leading eigenvalues and both left and right eigenvectors from a dense
Lindbladian matrix via full eigendecomposition.

Uses dense `eigen(L)` (feasible for n<=6, i.e. 4096x4096) and extracts left
eigenvectors via `inv(V)` to ensure biorthonormality: `V_left' * V_right = I`.

Eigenvalues are sorted by |Re(lambda)| (proximity to zero), with the steady
state (lambda ~ 0) at index 1.
"""
function extract_leading_eigendata(L::Matrix{ComplexF64}; n_modes::Int=20)
    d2 = size(L, 1)
    n_modes = min(n_modes, d2)

    # Dense eigendecomposition (feasible for n<=6: 4096x4096)
    F = eigen(L)

    # Sort by proximity to zero: |Re(lambda)|
    perm = sortperm(abs.(real.(F.values)))
    eigenvalues = F.values[perm[1:n_modes]]

    # Right eigenvectors (columns) for selected modes
    V_full = F.vectors[:, perm]
    V_right = V_full[:, 1:n_modes]

    # Left eigenvectors via inv(V_full): rows of V^{-1} are left eigenvectors.
    # Transpose to get them as columns: V_left = V_inv[1:n_modes, :]'
    # Then V_left' * V_right = I (biorthonormality).
    V_inv = inv(V_full)
    V_left = V_inv[1:n_modes, :]'

    spectral_gap = abs(real(eigenvalues[2]))

    # Im/Re ratios: for mode 1 (steady state) set to 0.0
    im_re_ratios = Vector{Float64}(undef, n_modes)
    im_re_ratios[1] = 0.0
    for k in 2:n_modes
        im_re_ratios[k] = abs(imag(eigenvalues[k])) / max(abs(real(eigenvalues[k])), 1e-30)
    end

    return EigenDecompositionResult(eigenvalues, V_right, V_left, spectral_gap, im_re_ratios)
end

# ---------------------------------------------------------------------------
# DIAG-02: Fixed point distance
# ---------------------------------------------------------------------------

"""
    compute_fixed_point_distance(eigen_result::EigenDecompositionResult, gibbs::Hermitian) -> FixedPointResult

Compute trace distance between the Lindbladian fixed point (lambda_1 eigenvector)
and the Gibbs state. Near-zero for well-constructed Lindbladians (BohrDomain).
"""
function compute_fixed_point_distance(eigen_result::EigenDecompositionResult, gibbs::Hermitian)
    fp_vec = eigen_result.right_eigenvectors[:, 1]
    dim = isqrt(length(fp_vec))
    fp_dm = reshape(copy(fp_vec), dim, dim)

    # Normalize: make Hermitian, trace = 1
    hermitianize!(fp_dm)
    fp_dm ./= tr(fp_dm)

    dist = trace_distance_h(Hermitian(fp_dm), gibbs)
    return FixedPointResult(fp_dm, dist)
end

# ---------------------------------------------------------------------------
# DIAG-03/04: Anti-Hermitian defect via KMS similarity transform
# ---------------------------------------------------------------------------

"""
    compute_anti_hermitian_defect(L::Matrix{ComplexF64}, gibbs::Hermitian; eps_trunc::Float64=1e-12) -> DefectResult

Compute the anti-Hermitian defect ratio of the KMS-similarity-transformed Lindbladian.

The KMS transform uses the diagonal Gibbs state in the eigenbasis:
`D = diag(rho^{-1/4} kron rho^{-1/4}) * L * diag(rho^{1/4} kron rho^{1/4})`

Returns the defect ratio `||A|| / lambda_gap(H_D)` where H_D is the Hermitian
part of D. An advisory warning is emitted when the ratio exceeds 0.1.
"""
function compute_anti_hermitian_defect(L::Matrix{ComplexF64}, gibbs::Hermitian;
                                        eps_trunc::Float64=1e-12)
    # Extract Gibbs diagonal (Gibbs is diagonal in Hamiltonian eigenbasis)
    gibbs_diag = real.(diag(Matrix(gibbs)))
    gibbs_diag_safe = max.(gibbs_diag, eps_trunc)

    rho_quarter = gibbs_diag_safe .^ 0.25
    rho_inv_quarter = gibbs_diag_safe .^ (-0.25)

    # Diagonal KMS transform: D_ij = d_left[i] * L_ij * d_right[j]
    d_right = kron(rho_quarter, rho_quarter)
    d_left = kron(rho_inv_quarter, rho_inv_quarter)

    D_matrix = d_left .* L .* d_right'

    # Hermitian / anti-Hermitian split
    H_part = (D_matrix + D_matrix') / 2
    A_part = (D_matrix - D_matrix') / 2

    A_norm = opnorm(A_part)  # Operator 2-norm (largest singular value)

    # Gap of Hermitian part: second smallest absolute eigenvalue
    H_eigenvalues = eigvals(Hermitian(H_part))
    sorted_H_abs = sort(abs.(H_eigenvalues))
    H_gap = sorted_H_abs[2]

    defect_ratio = A_norm / max(H_gap, 1e-30)

    threshold = 0.1
    warning = defect_ratio > threshold
    if warning
        @warn "Anti-Hermitian defect ratio $(round(defect_ratio; digits=4)) > $(threshold) threshold. " *
              "Non-normality effects may cause oscillatory transients -- consider oscillatory fit model."
    end

    return DefectResult(A_norm, H_gap, defect_ratio, warning, threshold)
end
