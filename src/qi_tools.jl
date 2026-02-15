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
    transform_jumps_to_basis(jumps, eigvecs) -> Vector{JumpOp}

Transform jump operators from computational basis to a new eigenbasis.
Used for Trotter basis transforms in TrotterDomain simulations.
"""
function transform_jumps_to_basis(jumps::AbstractVector{<:JumpOp}, eigvecs::AbstractMatrix)
    return JumpOp[JumpOp(j.data, eigvecs' * j.data * eigvecs, j.orthogonal, j.hermitian) for j in jumps]
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
    Takes in a jump operator, vectorizes it and adds it to the target liouvillian (Watrous convention).
    L = J1 * X * J2 - 0.5 * (J2 * J1 * X + X * J2 * J1)
    The way it is coded, makes it sure, that we allocate the least amount of times, i.e. the liouvillian is only allocated 
    once, there is no unnecessary copies during the calculations. 
    Note, that since it adds the liouvillian parts one by one to the liouvillian, on large scales
    there is some arithmetic error to an implementation that adds the parts together and then to the liouvillian.
"""
function _vectorize_liouv_diss_and_add!(
    L_target::AbstractMatrix{<:Complex},
    jump::AbstractMatrix{<:Complex},
    scalar::Number,
    ws::LindbladianWorkspace
)

    # scratch buffers
    jump_conj = ws.jump_conj
    jump_dag_jump = ws.jump_dag_jump
    Id = ws.Id

    @. jump_conj = conj(jump)
    _kron!(L_target, jump, jump_conj, scalar)
    
    mul!(jump_dag_jump, jump', jump)
    _kron!(L_target, jump_dag_jump, Id, -0.5 * scalar)
    _kron!(L_target, Id, transpose(jump_dag_jump), -0.5 * scalar)
    
    return L_target
end

function _vectorize_liouv_diss_and_add!(
    L_target::AbstractMatrix{<:Complex},
    jump_1::AbstractMatrix{<:Complex},
    jump_2::AbstractMatrix{<:Complex},
    scalar::Number,
    ws::LindbladianWorkspace
)
    Id = ws.Id
    jump2_jump1 = ws.jump2_jump1

    _kron!(L_target, jump_1, transpose(jump_2), scalar)

    mul!(jump2_jump1, jump_2, jump_1)
    _kron!(L_target, jump2_jump1, Id, -0.5 * scalar)
    _kron!(L_target, Id, transpose(jump2_jump1), -0.5 * scalar)
    
    return L_target
end

function _vectorize_liouvillian_coherent!(
    L_target::AbstractMatrix{<:Complex},
    coherent_term::AbstractMatrix{<:Complex},
    ws::LindbladianWorkspace)

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

function frobenius_norm(A::Matrix{<:Complex})
    eig_vals = eigvals(A)
    return sqrt(sum(abs.(eig_vals).^2))
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

