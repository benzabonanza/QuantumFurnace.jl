using LinearAlgebra
using SparseArrays
using Random
using Printf

"""
    hermitianize!(A::AbstractMatrix) -> A

In-place Hermitianization: A .= (A + A') / 2.
Used to enforce Hermiticity after numerical accumulation.
"""
function hermitianize!(A::AbstractMatrix)
    A .= 0.5 .* (A .+ A')
    return A
end

"""
    Computes C .+= alpha .* kron(A, B) completely in-place, without allocating
    the result of the Kronecker product. Speed.
"""
function _kron!(
    C::AbstractMatrix,
    A::AbstractMatrix,
    B::AbstractMatrix,
    alpha::Number
)
    m_a, n_a = size(A)
    m_b, n_b = size(B)
    
    # This check is good for debugging but can be removed for performance
    # @assert size(C) == (m_a * m_b, n_a * n_b) "Output matrix C has incorrect dimensions."

    for j in 1:n_a
        for i in 1:m_a
            a_ij = A[i, j]
            # If the element in A is zero, the whole block is zero.
            iszero(a_ij) && continue
            
            # Calculate the top-left corner of the block in C
            c_row_offset = (i - 1) * m_b
            c_col_offset = (j - 1) * n_b
            
            # Iterate over the B matrix
            val = alpha * a_ij
            for l in 1:n_b
                for k in 1:m_b
                    # C[block_row, block_col] += alpha * A[i,j] * B[k,l]
                    @inbounds C[c_row_offset + k, c_col_offset + l] += val * B[k, l]
                end
            end
        end
    end
    return C
end

"""
    Vectorize a single-operator Lindblad dissipator and add it to the target Liouvillian.

    Dissipator: J * rho * J' - 0.5 * {J'J, rho}

    Vectorization (Watrous/column-stacking convention):
      kron(conj(J), J) * vec(rho)         = vec(J * rho * J')       [sandwich]
      kron(I, J'J) * vec(rho)             = vec(J'J * rho)          [left anticommutator]
      kron((J'J)^T, I) * vec(rho)         = vec(rho * J'J)          [right anticommutator]

    Coded to minimise allocations: the Liouvillian is allocated once, and partial
    Kronecker products are accumulated in-place one by one.
"""
function _vectorize_liouv_diss_and_add!(
    L_target::AbstractMatrix{<:Complex},
    jump::AbstractMatrix{<:Complex},
    scalar::Number,
    ws::Workspace{Lindbladian},
)
    Id = ws.Id

    # scratch buffers (nested in LiouvillianScratch)
    jump_conj = ws.scratch.jump_conj
    jump_dag_jump = ws.scratch.jump_dag_jump

    @. jump_conj = conj(jump)
    _kron!(L_target, jump_conj, jump, scalar)                        # kron(conj(J), J) => J*rho*J'

    mul!(jump_dag_jump, jump', jump)                                 # J'J (Hermitian, convention-independent)
    _kron!(L_target, Id, jump_dag_jump, -0.5 * scalar)              # kron(I, J'J) => J'J*rho
    _kron!(L_target, transpose(jump_dag_jump), Id, -0.5 * scalar)   # kron((J'J)^T, I) => rho*J'J

    return L_target
end

"""
    Vectorize a two-operator Lindblad dissipator and add it to the target Liouvillian.

    Dissipator: J1 * X * J2 - 0.5 * (J2 * J1 * X + X * J2 * J1)

    Vectorization:
      kron(J2^T, J1) * vec(rho) = vec(J1 * rho * J2)       [sandwich]
      kron(I, J2*J1)            = vec(J2*J1 * rho)          [left anticommutator]
      kron((J2*J1)^T, I)        = vec(rho * J2*J1)          [right anticommutator]
"""
function _vectorize_liouv_diss_and_add!(
    L_target::AbstractMatrix{<:Complex},
    jump_1::AbstractMatrix{<:Complex},
    jump_2::AbstractMatrix{<:Complex},
    scalar::Number,
    ws::Workspace{Lindbladian},
)
    Id = ws.Id
    jump2_jump1 = ws.scratch.jump2_jump1

    _kron!(L_target, transpose(jump_2), jump_1, scalar)              # kron(J2^T, J1) => J1*rho*J2

    mul!(jump2_jump1, jump_2, jump_1)                                # J2*J1
    _kron!(L_target, Id, jump2_jump1, -0.5 * scalar)                # kron(I, J2*J1) => J2*J1*rho
    _kron!(L_target, transpose(jump2_jump1), Id, -0.5 * scalar)     # kron((J2*J1)^T, I) => rho*J2*J1

    return L_target
end

function _vectorize_liouvillian_coherent!(
    L_target::AbstractMatrix{<:Complex},
    coherent_term::AbstractMatrix{<:Complex},
    ws::Workspace{Lindbladian},
)
    Id = ws.Id
    _kron!(L_target, coherent_term, Id, -1im)
    _kron!(L_target, Id, transpose(coherent_term), +1im)
    return L_target
end

### ----------------------------- 
function trace_distance_h(rho::Union{Hermitian{<:Real}, Hermitian{<:Complex}}, 
    sigma::Union{Hermitian{<:Real}, Hermitian{<:Complex}})
    """Qutip apparently uses some sparse eigval solver, but let's go with the dense one for now."""
    return sum(abs.(eigvals(rho - sigma))) / 2
end

function trace_distance_nh(rho::Union{Matrix{<:Real}, Matrix{<:Complex}}, 
    sigma::Union{Matrix{<:Real}, Matrix{<:Complex}})
    return sum(svdvals(rho - sigma)) / 2
end

function trace_norm_h(rho::Union{Hermitian{<:Real}, Hermitian{<:Complex}})
    return sum(abs.(eigvals(rho)))
end

function trace_norm_nh(rho::Union{Matrix{<:Real}, Matrix{<:Complex}})
    return sum(svdvals(rho))
end

function fidelity(rho::Union{Hermitian{<:Real}, Hermitian{<:Complex}}, 
    sigma::Union{Hermitian{<:Real}, Hermitian{<:Complex}}; validate::Bool = true)

    if validate && (!is_density_matrix(rho) || !is_density_matrix(sigma))
        throw(ArgumentError("Input matrices are not density matrices"))
    end

    eig_vals = real(eigvals(rho * sigma))
    return real(sum(sqrt.(eig_vals[eig_vals.>0])))^2
end

function is_density_matrix(rho::Union{Hermitian{<:Real}, Hermitian{<:Complex}})
    if !isapprox(rho, rho')
        throw(ArgumentError("Input matrix is not Hermitian"))
    end

    eig_vals = real(round.(eigvals(rho), digits=15))
    # check if eigenvalues are approximately nonnegative
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
    # check if eigenvalues are approximately nonnegative
    if any(eig_vals .< 0)
        throw(ArgumentError("Input matrix has negative eigenvalues"))
    end

    if !isapprox(sum(eig_vals), 1.0)
        throw(ArgumentError("Input matrix has got trace different from 1"))
    end

    return true
end

function gibbs_state(hamiltonian::HamHam{T}, beta::Real) where {T<:AbstractFloat}
    """Computes Gibbs state in computational basis!"""
    CT = Complex{T}
    Z = sum(exp.(-beta * hamiltonian.eigvals))
    rho = sum([exp(-beta * hamiltonian.eigvals[i]) * hamiltonian.eigvecs[:, i] * hamiltonian.eigvecs[:, i]'
                                                                                    for i in 1:length(hamiltonian.eigvals)])
    return Matrix{CT}(rho / Z)
end

function gibbs_state_in_eigen(hamiltonian::HamHam{T}, beta::Real) where {T<:AbstractFloat}
    """Computes Gibbs state in eigenbasis"""
    CT = Complex{T}
    eigvecs_in_eigen = I(size(hamiltonian.data)[1])
    Z = sum(exp.(-beta * hamiltonian.eigvals))
    rho = sum([exp(-beta * hamiltonian.eigvals[i]) * eigvecs_in_eigen[:, i] * eigvecs_in_eigen[:, i]'
                                                                                    for i in 1:length(hamiltonian.eigvals)])
    return Matrix{CT}(rho / Z)
end

function random_density_matrix(num_qubits::Int)
    # Generate a random complex matrix
    A = randn(ComplexF64, 2^num_qubits, 2^num_qubits)

    # Compute A * A^†
    ρ = A * A'

    # Normalize the matrix to make the trace equal to 1
    ρ /= tr(ρ)

    return Hermitian(ρ)
end

"""
    validate_jump_pairing(jumps; allow_unpaired_nonhermitian=false, atol=1e-12)

Check that every non-Hermitian jump in `jumps` has a partner whose `data ≈
A†` (within `atol`) elsewhere in the set. KMS detailed balance for the
KMS-CKG construction requires this pairing — without it, the Kossakowski
α-skew-symmetry `α(ω₁,ω₂) = α(-ω₂,-ω₁) e^{-β(ω₁+ω₂)/2}` cannot be
satisfied (Chen et al. 2025, qf-bm1 Q1).

`allow_unpaired_nonhermitian=true` skips the check; intended only for unit
tests of internal code paths (e.g., serial-vs-threaded equivalence) that
compare two evaluations of the same physics on a non-physical fixture.
"""
function validate_jump_pairing(jumps::AbstractVector{<:JumpOp};
                                allow_unpaired_nonhermitian::Bool = false,
                                atol::Real = 1e-12)
    allow_unpaired_nonhermitian && return nothing

    unpaired_indices = Int[]
    for k in eachindex(jumps)
        jumps[k].hermitian && continue
        Adag = jumps[k].data'
        # Search for any other jump whose data ≈ A†
        found = false
        for j in eachindex(jumps)
            j == k && continue
            if size(jumps[j].data) == size(Adag) && isapprox(jumps[j].data, Adag; atol=atol)
                found = true
                break
            end
        end
        found || push!(unpaired_indices, k)
    end

    isempty(unpaired_indices) && return nothing

    throw(ArgumentError(
        "KMS detailed balance requires (A, A†) pairs for non-Hermitian jumps. " *
        "Found $(length(unpaired_indices)) unpaired non-Hermitian jump(s) at " *
        "index/indices $(unpaired_indices). Pass " *
        "`allow_unpaired_nonhermitian=true` for unit-test fixtures only."))
end

