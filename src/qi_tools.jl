using LinearAlgebra
using SparseArrays
using Random
using Printf

"""
    hermitianize!(A::AbstractMatrix) -> A

Replace `A` by its Hermitian part and return it.
"""
function hermitianize!(A::AbstractMatrix)
    # Math: $A <- (A + A^dagger) / 2$.
    A .= 0.5 .* (A .+ A')
    return A
end

"""
    _kron!(C, A, B, alpha) -> C

Accumulate `alpha * kron(A, B)` into `C` without materialising the product.
"""
function _kron!(
    C::AbstractMatrix,
    A::AbstractMatrix,
    B::AbstractMatrix,
    alpha::Number
)
    m_a, n_a = size(A)
    m_b, n_b = size(B)
    
    for j in 1:n_a
        for i in 1:m_a
            a_ij = A[i, j]
            iszero(a_ij) && continue

            c_row_offset = (i - 1) * m_b
            c_col_offset = (j - 1) * n_b

            val = alpha * a_ij
            for l in 1:n_b
                for k in 1:m_b
                    @inbounds C[c_row_offset + k, c_col_offset + l] += val * B[k, l]
                end
            end
        end
    end
    return C
end

"""
    _vectorize_liouv_diss_and_add!(L_target, jump, scalar, ws) -> L_target

Add a single-jump dissipator to a column-stacked Liouvillian.
"""
function _vectorize_liouv_diss_and_add!(
    L_target::AbstractMatrix{<:Complex},
    jump::AbstractMatrix{<:Complex},
    scalar::Number,
    ws::DenseLindbladianWorkspace,
)
    Id = ws.Id

    jump_conj = ws.scratch.jump_conj
    jump_dag_jump = ws.scratch.jump_dag_jump

    # Math: $vec(J rho J^dagger) = (conj(J) tensor J) vec(rho)$.
    @. jump_conj = conj(jump)
    _kron!(L_target, jump_conj, jump, scalar)

    mul!(jump_dag_jump, jump', jump)
    _kron!(L_target, Id, jump_dag_jump, -0.5 * scalar)
    _kron!(L_target, transpose(jump_dag_jump), Id, -0.5 * scalar)

    return L_target
end

"""
    _vectorize_liouv_diss_and_add!(L_target, jump_1, jump_2, scalar, ws)

Add a two-jump dissipator to a column-stacked Liouvillian.
"""
function _vectorize_liouv_diss_and_add!(
    L_target::AbstractMatrix{<:Complex},
    jump_1::AbstractMatrix{<:Complex},
    jump_2::AbstractMatrix{<:Complex},
    scalar::Number,
    ws::DenseLindbladianWorkspace,
)
    Id = ws.Id
    jump2_jump1 = ws.scratch.jump2_jump1

    # Math: $vec(J_1 rho J_2) = (J_2^T tensor J_1) vec(rho)$.
    _kron!(L_target, transpose(jump_2), jump_1, scalar)

    mul!(jump2_jump1, jump_2, jump_1)
    _kron!(L_target, Id, jump2_jump1, -0.5 * scalar)
    _kron!(L_target, transpose(jump2_jump1), Id, -0.5 * scalar)

    return L_target
end

function _vectorize_liouvillian_coherent!(
    L_target::AbstractMatrix{<:Complex},
    coherent_term::AbstractMatrix{<:Complex},
    ws::DenseLindbladianWorkspace,
)
    Id = ws.Id
    _kron!(L_target, coherent_term, Id, -1im)
    _kron!(L_target, Id, transpose(coherent_term), +1im)
    return L_target
end

"""Return the trace distance between Hermitian matrices."""
function trace_distance_h(rho::Union{Hermitian{<:Real}, Hermitian{<:Complex}}, 
    sigma::Union{Hermitian{<:Real}, Hermitian{<:Complex}})
    return sum(abs.(eigvals(rho - sigma))) / 2
end

"""Return the trace distance between general matrices using singular values."""
function trace_distance_nh(rho::Union{Matrix{<:Real}, Matrix{<:Complex}}, 
    sigma::Union{Matrix{<:Real}, Matrix{<:Complex}})
    return sum(svdvals(rho - sigma)) / 2
end

"""Return the trace norm of a Hermitian matrix from its eigenvalues."""
function trace_norm_h(rho::Union{Hermitian{<:Real}, Hermitian{<:Complex}})
    return sum(abs.(eigvals(rho)))
end

"""Return the trace norm of a general matrix from its singular values."""
function trace_norm_nh(rho::Union{Matrix{<:Real}, Matrix{<:Complex}})
    return sum(svdvals(rho))
end

"""
    fidelity(rho, sigma; validate=true) -> Real

Return the squared quantum-state fidelity, optionally validating both inputs.
"""
function fidelity(rho::Union{Hermitian{<:Real}, Hermitian{<:Complex}}, 
    sigma::Union{Hermitian{<:Real}, Hermitian{<:Complex}}; validate::Bool = true)

    if validate && (!is_density_matrix(rho) || !is_density_matrix(sigma))
        throw(ArgumentError("Input matrices are not density matrices"))
    end

    eig_vals = real(eigvals(rho * sigma))
    return real(sum(sqrt.(eig_vals[eig_vals.>0])))^2
end

"""
    is_density_matrix(rho) -> true

Validate Hermiticity, positive semidefiniteness, and unit trace.

Throws `ArgumentError` when an invariant fails.
"""
function is_density_matrix(rho::Union{Hermitian{<:Real}, Hermitian{<:Complex}})
    if !isapprox(rho, rho')
        throw(ArgumentError("Input matrix is not Hermitian"))
    end

    eig_vals = real(round.(eigvals(rho), digits=15))
    if any(eig_vals .< 0)
        throw(ArgumentError("Input matrix has negative eigenvalues"))
    end

    if !isapprox(sum(eig_vals), 1.0)
        throw(ArgumentError("Input matrix has got trace different from 1"))
    end

    return true
end

function is_density_matrix(rho::Hermitian{Complex{T}, Matrix{Complex{T}}}) where {T<:AbstractFloat}
    if !isapprox(rho, rho')
        throw(ArgumentError("Input matrix is not Hermitian"))
    end

    eig_vals = real(round.(eigvals(rho), digits=13))
    if any(eig_vals .< 0)
        throw(ArgumentError("Input matrix has negative eigenvalues"))
    end

    if !isapprox(sum(eig_vals), 1.0)
        throw(ArgumentError("Input matrix has got trace different from 1"))
    end

    return true
end

"""
    gibbs_state(hamiltonian, beta) -> Matrix

Return `\$rho_beta = exp(-beta H) / Z\$` in the computational basis.
"""
function gibbs_state(hamiltonian::HamHam{T}, beta::Real) where {T<:AbstractFloat}
    CT = Complex{T}
    Z = sum(exp.(-beta * hamiltonian.eigvals))
    rho = sum([exp(-beta * hamiltonian.eigvals[i]) * hamiltonian.eigvecs[:, i] * hamiltonian.eigvecs[:, i]'
                                                                                    for i in 1:length(hamiltonian.eigvals)])
    return Matrix{CT}(rho / Z)
end

"""
    gibbs_state_in_eigen(hamiltonian, beta) -> Matrix

Return `\$rho_beta = exp(-beta H) / Z\$` in the Hamiltonian eigenbasis.
"""
function gibbs_state_in_eigen(hamiltonian::HamHam{T}, beta::Real) where {T<:AbstractFloat}
    CT = Complex{T}
    eigvecs_in_eigen = I(size(hamiltonian.data)[1])
    Z = sum(exp.(-beta * hamiltonian.eigvals))
    rho = sum([exp(-beta * hamiltonian.eigvals[i]) * eigvecs_in_eigen[:, i] * eigvecs_in_eigen[:, i]'
                                                                                    for i in 1:length(hamiltonian.eigvals)])
    return Matrix{CT}(rho / Z)
end

"""Return a random `num_qubits`-qubit density matrix from a Ginibre draw."""
function random_density_matrix(num_qubits::Int)
    A = randn(ComplexF64, 2^num_qubits, 2^num_qubits)
    ρ = A * A'
    ρ /= tr(ρ)

    return Hermitian(ρ)
end

"""
    validate_jump_pairing(jumps; allow_unpaired_nonhermitian=false, atol=1e-12)

Validate that stored jump operators form an adjoint-closed multiset.

# Keywords
- `allow_unpaired_nonhermitian`: Skip closure checks for non-physical diagnostics.
- `atol`: Absolute tolerance in both stored bases.

# Returns
`nothing`; throws `ArgumentError` on invalid Hermitian flags or missing adjoints.

KMS detailed balance requires
`\$alpha(omega_1, omega_2) = alpha(-omega_2, -omega_1) exp(-beta(omega_1 + omega_2)/2)\$`.
"""
function validate_jump_pairing(jumps::AbstractVector{<:JumpOp};
                                allow_unpaired_nonhermitian::Bool = false,
                                atol::Real = 1e-12)
    invalid_hermitian_indices = Int[]
    for k in eachindex(jumps)
        jumps[k].hermitian || continue
        data_ok = isapprox(jumps[k].data, jumps[k].data'; atol=atol)
        eigenbasis_ok = isapprox(
            jumps[k].in_eigenbasis, jumps[k].in_eigenbasis'; atol=atol)
        (data_ok && eigenbasis_ok) || push!(invalid_hermitian_indices, k)
    end
    isempty(invalid_hermitian_indices) || throw(ArgumentError(
        "Jump(s) marked hermitian=true are not Hermitian in both stored bases " *
        "at index/indices $(invalid_hermitian_indices)."))

    allow_unpaired_nonhermitian && return nothing

    matched = falses(length(jumps))
    unpaired_indices = Int[]
    for k in eachindex(jumps)
        (jumps[k].hermitian || matched[k]) && continue
        data_adjoint = jumps[k].data'
        eigenbasis_adjoint = jumps[k].in_eigenbasis'
        partner = nothing
        for j in eachindex(jumps)
            (j == k || matched[j] || jumps[j].hermitian) && continue
            data_ok = size(jumps[j].data) == size(data_adjoint) &&
                isapprox(jumps[j].data, data_adjoint; atol=atol)
            eigenbasis_ok = size(jumps[j].in_eigenbasis) == size(eigenbasis_adjoint) &&
                isapprox(jumps[j].in_eigenbasis, eigenbasis_adjoint; atol=atol)
            if data_ok && eigenbasis_ok
                partner = j
                break
            end
        end
        if partner === nothing
            push!(unpaired_indices, k)
        else
            matched[k] = true
            matched[partner] = true
        end
    end

    isempty(unpaired_indices) && return nothing

    throw(ArgumentError(
        "KMS detailed balance requires (A, A†) pairs for non-Hermitian jumps. " *
        "Found $(length(unpaired_indices)) unpaired non-Hermitian jump(s) at " *
        "index/indices $(unpaired_indices). Pass " *
        "`allow_unpaired_nonhermitian=true` for unit-test fixtures only."))
end
