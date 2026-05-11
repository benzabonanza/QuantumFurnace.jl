"""
    compute_oft_trotter_error(jump_data, hamiltonian, trotter; kwargs...)

Compare Operator Fourier Transforms computed with exact time evolution (TimeDomain)
vs Trotter time evolution (TrotterDomain).

For each sampled energy label ω, computes:
  - A_exact(ω) via `time_oft!` (uses exact eigenvalues exp(iEt))
  - A_trotter(ω) via `trotter_oft!` (uses Trotter eigenvalues λ^n)

Both are transformed to computational basis and their operator-norm difference is recorded.

Returns a NamedTuple with fields:
  - `max_error`: Maximum ‖A_exact(ω) - A_trotter(ω)‖ across sampled ω
  - `mean_error`: Mean error across sampled ω
  - `errors`: Vector of per-ω errors
  - `energy_labels_sampled`: The ω values that were sampled
"""
function compute_oft_trotter_error(
    jump_data::Matrix{ComplexF64},
    hamiltonian::HamHam,
    trotter::AbstractTrotter;
    num_energy_bits::Int = 12,
    w0::Float64 = 0.05,
    sigma::Float64 = 0.08,
    energy_sample_stride::Int = 1,
)
    dim = size(jump_data, 1)
    t0 = Float64(trotter.t0)

    # Create labels (same as simulation pipeline)
    energy_labels = _create_energy_labels(num_energy_bits, w0)
    time_labels = energy_labels .* (t0 / w0)

    # JumpOps in both eigenbases
    jump_ham_eigbasis = hamiltonian.eigvecs' * jump_data * hamiltonian.eigvecs
    jump_trott_eigbasis = trotter.eigvecs' * jump_data * trotter.eigvecs

    is_orth = isapprox(jump_data, transpose(jump_data))
    is_herm = isapprox(jump_data, jump_data')

    jump_ham = JumpOp(jump_data, jump_ham_eigbasis, is_orth, is_herm)
    jump_trott = JumpOp(jump_data, jump_trott_eigbasis, is_orth, is_herm)

    # Allocate caches
    caches_time = OFTCaches{Float64}(dim)
    caches_trott = OFTCaches{Float64}(dim)
    out_time = zeros(ComplexF64, dim, dim)
    out_trotter = zeros(ComplexF64, dim, dim)

    # Scratch for basis transformation
    tmp = zeros(ComplexF64, dim, dim)

    # Sample energy labels
    sample_indices = 1:energy_sample_stride:length(energy_labels)
    errors = Vector{Float64}(undef, length(sample_indices))

    for (k, idx) in enumerate(sample_indices)
        energy = energy_labels[idx]

        time_oft!(out_time, caches_time, jump_ham, energy, hamiltonian, time_labels, sigma)
        trotter_oft!(out_trotter, caches_trott, jump_trott, energy, trotter, time_labels, sigma)

        # Transform to computational basis
        # A_comp = V * A_eig * V'
        mul!(tmp, hamiltonian.eigvecs, out_time)
        mul!(out_time, tmp, hamiltonian.eigvecs')  # reuse out_time as A_exact_comp

        mul!(tmp, trotter.eigvecs, out_trotter)
        mul!(out_trotter, tmp, trotter.eigvecs')  # reuse out_trotter as A_trott_comp

        out_time .-= out_trotter
        errors[k] = norm(out_time)
    end

    return (
        max_error = maximum(errors),
        mean_error = sum(errors) / length(errors),
        errors = errors,
        energy_labels_sampled = energy_labels[sample_indices],
    )
end

"""
    compute_oft_trotter_error_all_jumps(jump_paulis, num_qubits, hamiltonian, trotter; kwargs...)

Compute the OFT Trotter error aggregated over all jump operators (Pauli on each site).
Returns per-jump errors plus aggregate statistics.
"""
function compute_oft_trotter_error_all_jumps(
    jump_paulis::Vector{Vector{Matrix{ComplexF64}}},
    num_qubits::Int,
    hamiltonian::HamHam,
    trotter::AbstractTrotter;
    kwargs...
)
    norm_factor = sqrt(length(jump_paulis) * num_qubits)
    results = []

    for pauli in jump_paulis
        for site in 1:num_qubits
            op = Matrix(pad_term(pauli, num_qubits, site)) ./ norm_factor
            r = compute_oft_trotter_error(op, hamiltonian, trotter; kwargs...)
            push!(results, (pauli=pauli, site=site, max_error=r.max_error, mean_error=r.mean_error))
        end
    end

    max_errors = [r.max_error for r in results]
    mean_errors = [r.mean_error for r in results]

    return (
        per_jump = results,
        aggregate_max = maximum(max_errors),
        aggregate_mean = sum(mean_errors) / length(mean_errors),
    )
end
