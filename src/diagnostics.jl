# Dense exact diagnostics for small systems.

"""
    EigenDecompositionResult

Leading eigenvalues and biorthogonal left/right eigenvectors.

# Fields
- `eigenvalues`: Modes sorted by `abs(real(lambda))`.
- `right_eigenvectors`, `left_eigenvectors`: Biorthogonal mode columns.
- `spectral_gap`: `abs(real(eigenvalues[2]))`.
- `im_re_ratios`: Oscillation-to-decay ratio per mode.
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

Normalised Lindbladian fixed point and its Gibbs-state trace distance.

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

Anti-Hermitian defect of the KMS quantum discriminant.

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

Observable overlaps with Lindbladian eigenmodes.

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

Dominant `Delta_Sz` sector and its weight for one eigenmode.

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

Group of near-degenerate Lindbladian eigenvalues.

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
    SpectralModeDiagnostics

Per-mode properties of a biorthogonal Krylov decomposition.

# Fields
- `off_diag_weight`: Hilbert–Schmidt fraction outside the diagonal.
- `c_abs2`: Squared initial-state modal coefficient, or `NaN`.
- `modal_hs_weight`: Scale-invariant value `abs2(c[k]) * norm(R[k])^2`, or `NaN`.
- `mode_spacing`: Complex-plane distance to the next captured eigenvalue.

Channel predictors report spacing in channel-eigenvalue units; gap solvers
first convert channel eigenvalues to generator rates.
"""
struct SpectralModeDiagnostics
    off_diag_weight::Vector{Float64}
    c_abs2::Vector{Float64}
    modal_hs_weight::Vector{Float64}
    mode_spacing::Vector{Float64}
end

"""
    spectral_mode_diagnostics(eigenvalues, R_modes, c=nothing)

Compute coherence, amplitude, and spacing diagnostics for captured modes.

# Arguments
- `eigenvalues`: Captured generator or channel eigenvalues.
- `R_modes`: Corresponding right modes as matrices.
- `c`: Optional initial-state coefficients.

# Returns
A `SpectralModeDiagnostics`; coefficient-dependent fields are `NaN` when
`c === nothing`.
"""
function spectral_mode_diagnostics(
    eigenvalues::AbstractVector{<:Complex},
    R_modes::AbstractVector{<:AbstractMatrix},
    c::Union{Nothing, AbstractVector{<:Complex}} = nothing,
)
    m = length(R_modes)
    length(eigenvalues) == m || throw(ArgumentError(
        "eigenvalues and R_modes must have equal length (got $(length(eigenvalues)), $m)"))
    c === nothing || length(c) == m || throw(ArgumentError(
        "c must match R_modes length when provided (got $(length(c)), $m)"))
    odw = Vector{Float64}(undef, m)
    ca2 = Vector{Float64}(undef, m)
    mhw = Vector{Float64}(undef, m)
    spc = Vector{Float64}(undef, m)
    @inbounds for k in 1:m
        R = R_modes[k]
        hs2 = sum(abs2, R)
        diag2 = 0.0
        for i in 1:min(size(R, 1), size(R, 2))
            diag2 += abs2(R[i, i])
        end
        odw[k] = hs2 > 0 ? clamp(1.0 - diag2 / hs2, 0.0, 1.0) : 0.0
        if c === nothing
            ca2[k] = NaN
            mhw[k] = NaN
        else
            ca2[k] = abs2(c[k])
            mhw[k] = ca2[k] * hs2
        end
        spc[k] = k < m ? abs(eigenvalues[k] - eigenvalues[k + 1]) : Inf
    end
    return SpectralModeDiagnostics(odw, ca2, mhw, spc)
end

"""
    ExactDiagnosticsResult

Combined result from `run_exact_diagnostics`.

# Fields
- `eigen`: Dense eigendecomposition.
- `fixed_point`: Fixed-point comparison.
- `defect`: Detailed-balance defect.
- `overlaps`: One observable-overlap result per initial state.
- `sz_labels`: One symmetry label per mode.
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

"""
    extract_leading_eigendata(L::Matrix{ComplexF64}; n_modes::Int=20) -> EigenDecompositionResult

Extract leading left and right modes from a dense Lindbladian.

# Keywords
- `n_modes`: Maximum number of modes to retain.

# Returns
An `EigenDecompositionResult` sorted by `abs(real(lambda))`, with left modes
constructed so `\$V_L^dagger V_R = I\$`.
"""
function extract_leading_eigendata(L::Matrix{ComplexF64}; n_modes::Int=20)
    d2 = size(L, 1)
    n_modes = min(n_modes, d2)

    F = eigen(L)

    # Sort modes by decay rate, placing the stationary mode first.
    perm = sortperm(abs.(real.(F.values)))
    eigenvalues = F.values[perm[1:n_modes]]

    V_full = F.vectors[:, perm]
    V_right = V_full[:, 1:n_modes]

    # Rows of `inv(V_full)` are the biorthogonal left eigenvectors.
    V_inv = inv(V_full)
    V_left = V_inv[1:n_modes, :]'

    spectral_gap = abs(real(eigenvalues[2]))

    im_re_ratios = Vector{Float64}(undef, n_modes)
    im_re_ratios[1] = 0.0
    for k in 2:n_modes
        im_re_ratios[k] = abs(imag(eigenvalues[k])) / max(abs(real(eigenvalues[k])), 1e-30)
    end

    return EigenDecompositionResult(eigenvalues, V_right, V_left, spectral_gap, im_re_ratios)
end

"""
    compute_fixed_point_distance(eigen_result::EigenDecompositionResult, gibbs::Hermitian) -> FixedPointResult

Normalise the stationary right mode and compare it with the Gibbs state.

# Returns
A `FixedPointResult` containing the density matrix and trace distance.
"""
function compute_fixed_point_distance(eigen_result::EigenDecompositionResult, gibbs::Hermitian)
    fp_vec = eigen_result.right_eigenvectors[:, 1]
    dim = isqrt(length(fp_vec))
    fp_dm = reshape(copy(fp_vec), dim, dim)

    # Math: $rho_infinity <- herm(R_1) / tr(herm(R_1))$.
    hermitianize!(fp_dm)
    fp_dm ./= tr(fp_dm)

    dist = trace_distance_h(Hermitian(fp_dm), gibbs)
    return FixedPointResult(fp_dm, dist)
end

"""
    compute_anti_hermitian_defect(L::Matrix{ComplexF64}, gibbs::Hermitian; eps_trunc::Float64=1e-12) -> DefectResult

Compute the anti-Hermitian defect of the KMS quantum discriminant.

# Keywords
- `eps_trunc`: Eigenvalue floor in Gibbs fractional powers.

# Returns
A `DefectResult` with ratio `norm(A) / gap(H)` and an advisory warning flag.
"""
function compute_anti_hermitian_defect(L::Matrix{ComplexF64}, gibbs::Hermitian;
                                        eps_trunc::Float64=1e-12)
    # Math: $D = H + A$, where $H = (D + D^dagger)/2$ and $A = (D-D^dagger)/2$.
    D_matrix = materialize_discriminant(L, gibbs; eps_trunc=eps_trunc)
    H_part, A_part = hermitian_antihermitian_split(D_matrix)

    A_norm = opnorm(A_part)

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

"""
    compute_overlap_coefficients(eigen_result, observables, observable_names, rho0, rho_beta;
                                  n_modes=20, initial_state_name="custom") -> OverlapResult

Compute observable coefficients in a biorthogonal Lindbladian expansion.

# Returns
An `OverlapResult` with one row per observable and one column per retained mode.
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

        # Math: $c_k = tr(O R_k) tr(L_k^dagger (rho_0-rho_beta))$.
        lk_factor = dot(vec(L_k), vec(rho_diff))

        for (i, O) in enumerate(observables)
            ok_factor = tr(O * R_k)
            coeffs[i, k] = ok_factor * lk_factor
        end
    end

    gap_mode_overlap = Float64[abs(coeffs[i, 2]) for i in 1:n_obs]

    return OverlapResult(coeffs, observable_names, initial_state_name, gap_mode_overlap)
end

"""
    compute_sz_labels(eigen_result, eigvecs, n_qubits; n_modes=20) -> Vector{SzSectorLabel}

Assign a dominant `Delta_Sz` sector to each Lindbladian right mode.

# Arguments
- `eigen_result`: Captured right modes.
- `eigvecs`: Columns defining the working basis.
- `n_qubits`: System size.
- `n_modes`: Maximum number of modes to label.

# Returns
A vector of `SzSectorLabel` values.
"""
function compute_sz_labels(eigen_result::EigenDecompositionResult, eigvecs::Matrix{<:Complex},
                            n_qubits::Int; n_modes::Int=20)
    dim = size(eigvecs, 1)
    n_modes_actual = min(n_modes, length(eigen_result.eigenvalues))

    # Math: $S_z = 1/2 sum_i Z_i$.
    Sz_comp = zeros(ComplexF64, dim, dim)
    for site in 1:n_qubits
        Sz_comp .+= Matrix{ComplexF64}(pad_term([Z], n_qubits, site))
    end
    Sz_comp ./= 2

    V = eigvecs
    Sz_eigen = V' * Sz_comp * V
    sz_vals = real.(diag(Sz_eigen))

    labels = Vector{SzSectorLabel}(undef, n_modes_actual)

    for k in 1:n_modes_actual
        M_k = reshape(eigen_result.right_eigenvectors[:, k], dim, dim)
        weights = abs2.(M_k)

        delta_sz_map = Dict{Float64, Float64}()
        for j in 1:dim, i in 1:dim
            w = weights[i, j]
            w < 1e-14 && continue
            dsz = round(sz_vals[i] - sz_vals[j]; digits=6)
            delta_sz_map[dsz] = get(delta_sz_map, dsz, 0.0) + w
        end

        total_weight = sum(values(delta_sz_map))

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

Label modes in `hamiltonian.eigvecs`.
"""
function compute_sz_labels(eigen_result::EigenDecompositionResult, hamiltonian::HamHam;
                            n_modes::Int=20)
    n_qubits = Int(log2(size(hamiltonian.data, 1)))
    return compute_sz_labels(eigen_result, hamiltonian.eigvecs, n_qubits; n_modes=n_modes)
end

"""
    detect_multiplets(eigenvalues::Vector{ComplexF64}; rel_tol=0.01) -> Vector{MultipletGroup}

Group adjacent eigenvalues by relative complex-plane spacing.

# Keywords
- `rel_tol`: Maximum relative spacing within a group.

# Returns
Multiplets ordered by eigenvalue magnitude.
"""
function detect_multiplets(eigenvalues::Vector{ComplexF64}; rel_tol::Float64=0.01)
    n = length(eigenvalues)
    n == 0 && return MultipletGroup[]

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
            mean_val = current_sum / length(current_indices)
            push!(groups, MultipletGroup(copy(current_indices), mean_val, SzSectorLabel[]))
            current_indices = [idx]
            current_sum = val
        end
    end

    mean_val = current_sum / length(current_indices)
    push!(groups, MultipletGroup(copy(current_indices), mean_val, SzSectorLabel[]))

    return groups
end

"""
    run_exact_diagnostics(L, hamiltonian, gibbs; kwargs...) -> ExactDiagnosticsResult

Run the complete dense diagnostic bundle.

# Arguments
- `L`: Dense Lindbladian superoperator.
- `hamiltonian`: Hamiltonian and basis data.
- `gibbs`: Gibbs state in the selected working basis.

# Keywords
- `basis_eigvecs`: Working basis; defaults to the Hamiltonian eigenbasis.
- `observables`, `observable_names`: Optional observable set and labels.
- `initial_states`, `initial_state_names`: Optional initial states and labels.
- `n_modes`: Maximum number of modes.
- `eps_trunc`: Gibbs fractional-power floor.

# Returns
An `ExactDiagnosticsResult`.
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

    V = basis_eigvecs === nothing ? hamiltonian.eigvecs : Matrix{ComplexF64}(basis_eigvecs)

    eigen_result = extract_leading_eigendata(L; n_modes=n_modes)

    fp_result = compute_fixed_point_distance(eigen_result, gibbs)

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

    overlaps_vec = OverlapResult[]
    for (rho0, name) in zip(initial_states, initial_state_names)
        overlap = compute_overlap_coefficients(
            eigen_result, observables, observable_names, rho0, gibbs;
            n_modes=n_modes, initial_state_name=name,
        )
        push!(overlaps_vec, overlap)
    end

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
