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
    # Materialise the discriminant and split into Hermitian / anti-Hermitian
    # parts.  See src/discriminant.jl for the column-stacking convention.
    D_matrix = materialize_discriminant(L, gibbs; eps_trunc=eps_trunc)
    H_part, A_part = hermitian_antihermitian_split(D_matrix)

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

# ---------------------------------------------------------------------------
# DIAG-05: Observable overlap coefficients
# ---------------------------------------------------------------------------

"""
    compute_overlap_coefficients(eigen_result, observables, observable_names, rho0, rho_beta;
                                  n_modes=20, initial_state_name="custom") -> OverlapResult

Compute overlap coefficients between observables and Lindbladian eigenmodes.

Uses the formula: `c_k = Tr[O R_k] * Tr[L_k^dagger (rho_0 - rho_beta)]`
where R_k are right eigenvectors and L_k are left eigenvectors (biorthonormal).

The coefficient c_1 (steady-state mode) should be near zero due to the
explicit subtraction of rho_beta.
"""
function compute_overlap_coefficients(
    eigen_result::EigenDecompositionResult,
    observables::Vector{<:Matrix{<:Complex}},
    observable_names::Vector{String},
    rho0::Matrix{<:Complex},
    rho_beta::Hermitian;
    n_modes::Int=20,
    initial_state_name::String="custom",
)
    dim = size(rho0, 1)
    n_modes_actual = min(n_modes, length(eigen_result.eigenvalues))
    n_obs = length(observables)
    rho_diff = rho0 - Matrix(rho_beta)

    coeffs = zeros(ComplexF64, n_obs, n_modes_actual)

    for k in 1:n_modes_actual
        R_k = reshape(eigen_result.right_eigenvectors[:, k], dim, dim)
        L_k = reshape(eigen_result.left_eigenvectors[:, k], dim, dim)

        # Tr[L_k^dagger (rho_0 - rho_beta)] = dot(vec(L_k), vec(rho_diff))
        # Julia's dot conjugates the first argument: dot(a,b) = sum(conj(a).*b)
        lk_factor = dot(vec(L_k), vec(rho_diff))

        for (i, O) in enumerate(observables)
            ok_factor = tr(O * R_k)
            coeffs[i, k] = ok_factor * lk_factor
        end
    end

    gap_mode_overlap = Float64[abs(coeffs[i, 2]) for i in 1:n_obs]

    return OverlapResult(coeffs, observable_names, initial_state_name, gap_mode_overlap)
end

# ---------------------------------------------------------------------------
# DIAG-06: Delta_Sz symmetry sector labeling
# ---------------------------------------------------------------------------

"""
    compute_sz_labels(eigen_result, eigvecs, n_qubits; n_modes=20) -> Vector{SzSectorLabel}

Assign Delta_Sz quantum numbers to each Lindbladian eigenvector based on the
density matrix support structure.

For each eigenvector R_k (reshaped as a dim x dim matrix), computes the weight
in each (i,j) element and groups by Delta_Sz = Sz(E_i) - Sz(E_j). Reports the
dominant sector and purity fraction.

# Arguments
- `eigen_result`: EigenDecompositionResult from DIAG-01.
- `eigvecs`: Unitary matrix whose columns define the working basis (e.g.
  `hamiltonian.eigvecs` for BohrDomain, `trotter.eigvecs` for TrotterDomain).
- `n_qubits`: Number of qubits in the system.
- `n_modes`: Number of leading modes to label (default: 20).
"""
function compute_sz_labels(eigen_result::EigenDecompositionResult, eigvecs::Matrix{<:Complex},
                            n_qubits::Int; n_modes::Int=20)
    dim = size(eigvecs, 1)
    n_modes_actual = min(n_modes, length(eigen_result.eigenvalues))

    # Build total Sz operator in computational basis: sum(Z_i)/2
    Sz_comp = zeros(ComplexF64, dim, dim)
    for site in 1:n_qubits
        Sz_comp .+= Matrix{ComplexF64}(pad_term([Z], n_qubits, site))
    end
    Sz_comp ./= 2

    # Transform to working eigenbasis
    V = eigvecs
    Sz_eigen = V' * Sz_comp * V
    sz_vals = real.(diag(Sz_eigen))

    labels = Vector{SzSectorLabel}(undef, n_modes_actual)

    for k in 1:n_modes_actual
        M_k = reshape(eigen_result.right_eigenvectors[:, k], dim, dim)
        weights = abs2.(M_k)

        # Compute Delta_Sz weight map
        delta_sz_map = Dict{Float64, Float64}()
        for j in 1:dim, i in 1:dim
            w = weights[i, j]
            w < 1e-14 && continue
            dsz = round(sz_vals[i] - sz_vals[j]; digits=6)
            delta_sz_map[dsz] = get(delta_sz_map, dsz, 0.0) + w
        end

        total_weight = sum(values(delta_sz_map))

        # Find dominant sector
        dominant_dsz = 0.0
        dominant_weight = 0.0
        for (dsz, wt) in delta_sz_map
            if wt > dominant_weight
                dominant_weight = wt
                dominant_dsz = dsz
            end
        end

        purity = dominant_weight / max(total_weight, 1e-30)
        is_pure = purity > 0.95

        labels[k] = SzSectorLabel(dominant_dsz, purity, is_pure, delta_sz_map)
    end

    return labels
end

"""
    compute_sz_labels(eigen_result, hamiltonian::HamHam; n_modes=20) -> Vector{SzSectorLabel}

Convenience method: delegates to the eigvecs-based method using `hamiltonian.eigvecs`.
"""
function compute_sz_labels(eigen_result::EigenDecompositionResult, hamiltonian::HamHam;
                            n_modes::Int=20)
    n_qubits = Int(log2(size(hamiltonian.data, 1)))
    return compute_sz_labels(eigen_result, hamiltonian.eigvecs, n_qubits; n_modes=n_modes)
end

# ---------------------------------------------------------------------------
# Multiplet detection: near-degenerate eigenvalue grouping
# ---------------------------------------------------------------------------

"""
    detect_multiplets(eigenvalues::Vector{ComplexF64}; rel_tol=0.01) -> Vector{MultipletGroup}

Group near-degenerate eigenvalues into multiplets.

Two eigenvalues are in the same multiplet if
`|lambda_i - lambda_j| / max(|lambda_i|, |lambda_j|, 1e-10) < rel_tol`.

Uses sequential grouping on eigenvalues sorted by absolute value.
"""
function detect_multiplets(eigenvalues::Vector{ComplexF64}; rel_tol::Float64=0.01)
    n = length(eigenvalues)
    n == 0 && return MultipletGroup[]

    # Sort by absolute value for sequential grouping
    sorted_perm = sortperm(abs.(eigenvalues))
    sorted_vals = eigenvalues[sorted_perm]

    groups = Vector{MultipletGroup}()
    current_indices = [sorted_perm[1]]
    current_sum = sorted_vals[1]

    for i in 2:n
        idx = sorted_perm[i]
        val = sorted_vals[i]
        prev_val = sorted_vals[i-1]

        denom = max(abs(val), abs(prev_val), 1e-10)
        if abs(val - prev_val) / denom < rel_tol
            push!(current_indices, idx)
            current_sum += val
        else
            # Finalize current group
            mean_val = current_sum / length(current_indices)
            push!(groups, MultipletGroup(copy(current_indices), mean_val, SzSectorLabel[]))
            current_indices = [idx]
            current_sum = val
        end
    end

    # Finalize last group
    mean_val = current_sum / length(current_indices)
    push!(groups, MultipletGroup(copy(current_indices), mean_val, SzSectorLabel[]))

    return groups
end

# ---------------------------------------------------------------------------
# Bundle: run_exact_diagnostics
# ---------------------------------------------------------------------------

"""
    run_exact_diagnostics(L, hamiltonian, gibbs; kwargs...) -> ExactDiagnosticsResult

Run all six DIAG diagnostics in a single call, returning a bundled result.

# Arguments
- `L::Matrix{ComplexF64}`: Full dense Lindbladian superoperator.
- `hamiltonian::HamHam`: Hamiltonian with eigenbasis data.
- `gibbs::Hermitian`: Gibbs state (in the working basis -- Hamiltonian eigenbasis
  for BohrDomain, Trotter eigenbasis for TrotterDomain).

# Keyword Arguments
- `basis_eigvecs`: Unitary matrix defining the working basis. When `nothing` (default),
  uses `hamiltonian.eigvecs` (backward compatible with BohrDomain). For TrotterDomain,
  pass `trotter.eigvecs` so that default observables, initial states, and Sz labels
  are all constructed in the Trotter eigenbasis.
- `observables`: Observable matrices (in working basis). Default: [Z1, H].
- `observable_names`: Observable names. Default: ["Z1", "H"].
- `initial_states`: Initial density matrices. Default: [|0>^n, |+>^n, I/dim].
- `initial_state_names`: Names for initial states. Default: ["all_up", "all_plus", "maximally_mixed"].
- `n_modes`: Number of leading modes to extract (default: 20).
- `eps_trunc`: Truncation threshold for KMS transform (default: 1e-12).
"""
function run_exact_diagnostics(
    L::Matrix{ComplexF64},
    hamiltonian::HamHam,
    gibbs::Hermitian;
    basis_eigvecs::Union{Nothing, Matrix{<:Complex}}=nothing,
    observables::Union{Nothing, Vector{<:Matrix{<:Complex}}}=nothing,
    observable_names::Union{Nothing, Vector{String}}=nothing,
    initial_states::Union{Nothing, Vector{<:Matrix{<:Complex}}}=nothing,
    initial_state_names::Union{Nothing, Vector{String}}=nothing,
    n_modes::Int=20,
    eps_trunc::Float64=1e-12,
)
    dim = size(hamiltonian.data, 1)
    n = Int(log2(dim))

    # Working basis: trotter.eigvecs for TrotterDomain, hamiltonian.eigvecs otherwise
    V = basis_eigvecs === nothing ? hamiltonian.eigvecs : Matrix{ComplexF64}(basis_eigvecs)

    # DIAG-01: eigendata extraction
    eigen_result = extract_leading_eigendata(L; n_modes=n_modes)

    # DIAG-02: fixed point distance
    fp_result = compute_fixed_point_distance(eigen_result, gibbs)

    # DIAG-03/04: anti-Hermitian defect
    defect_result = compute_anti_hermitian_defect(L, gibbs; eps_trunc=eps_trunc)

    # Build default observables if not provided
    if observables === nothing
        # Z1 in working basis
        Z1_comp = Matrix{ComplexF64}(pad_term([Z], n, 1))
        Z1_eigen = Matrix{ComplexF64}(V' * Z1_comp * V)
        # H in working basis (diagonal when V = hamiltonian.eigvecs, non-diagonal otherwise)
        H_eigen = Matrix{ComplexF64}(V' * hamiltonian.data * V)
        observables = Matrix{ComplexF64}[Z1_eigen, H_eigen]
        observable_names = String["Z1", "H"]
    end

    # Build default initial states if not provided
    if initial_states === nothing
        # |0>^n (all spins up) -- transform to working basis
        psi0_comp = zeros(ComplexF64, dim)
        psi0_comp[1] = 1.0
        psi0_eigen = V' * psi0_comp
        rho_up = psi0_eigen * psi0_eigen'

        # |+>^n (all X-plus) -- transform to working basis
        psi_plus_comp = fill(ComplexF64(1 / sqrt(2^n)), 2^n)
        psi_plus_eigen = V' * psi_plus_comp
        rho_plus = psi_plus_eigen * psi_plus_eigen'

        # I/dim (maximally mixed) -- same in any basis
        rho_mixed = Matrix{ComplexF64}(I(dim) / dim)

        initial_states = Matrix{ComplexF64}[rho_up, rho_plus, rho_mixed]
        initial_state_names = String["all_up", "all_plus", "maximally_mixed"]
    end

    # DIAG-05: overlap coefficients for each initial state
    overlaps_vec = OverlapResult[]
    for (rho0, name) in zip(initial_states, initial_state_names)
        overlap = compute_overlap_coefficients(
            eigen_result, observables, observable_names, rho0, gibbs;
            n_modes=n_modes, initial_state_name=name,
        )
        push!(overlaps_vec, overlap)
    end

    # DIAG-06: symmetry sector labels (use working basis eigvecs)
    sz_labels = compute_sz_labels(eigen_result, V, n; n_modes=n_modes)

    # Multiplet detection
    multiplets = detect_multiplets(eigen_result.eigenvalues)

    # Fill multiplet sz_labels from computed labels
    for group in multiplets
        for idx in group.eigenvalue_indices
            if idx <= length(sz_labels)
                push!(group.sz_labels, sz_labels[idx])
            end
        end
    end

    return ExactDiagnosticsResult(eigen_result, fp_result, defect_result,
                                   overlaps_vec, sz_labels, multiplets)
end
