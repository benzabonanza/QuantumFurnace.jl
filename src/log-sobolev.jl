using LinearAlgebra
using Optim
using Printf

struct LSIFramework{T}
    dim::Int

    A::Matrix{T}            # Parameter matrix
    AdagA::Matrix{T}            # B = A'A
    sigma_AdagA::Matrix{T}  # sig^1/4 A'A sig^1/4
    gradient::Matrix{T}     # Gradient accumulator

    temp1::Matrix{T}       
    temp2::Matrix{T}        
    temp3::Matrix{T}        

    sigma_quarter::Matrix{T}   # Sigma^1/4
    sigma_half::Matrix{T}      # Sigma^1/2
    sigma_log::Matrix{T}       # log(Sigma)

    AdagA_vec::Vector{T}    # vec(A'A)
    L_AdagA_vec::Vector{T}  # vec(L(A'A))
end

function LSIFramework(dim::Int)
    T = ComplexF64
    dim2 = dim^2
    return LSIFramework{T}(
        dim,
        Matrix{T}(undef, dim, dim), Matrix{T}(undef, dim, dim), Matrix{T}(undef, dim, dim), Matrix{T}(undef, dim, dim),
        Matrix{T}(undef, dim, dim), Matrix{T}(undef, dim, dim), Matrix{T}(undef, dim, dim),
        Matrix{T}(undef, dim, dim), Matrix{T}(undef, dim, dim), Matrix{T}(undef, dim, dim),
        Vector{T}(undef, dim2), Vector{T}(undef, dim2)
    )
end

function compute_LSI_alpha2(result::HotSpectralResults; n_restarts=3, g_tol=1e-5)

    # In respective eigenspace either (H or Trotter)
    L_mat = result.data
    sigma = Hermitian(result.fixed_point)
    sigma_eigen = eigen(sigma)
    sigma_eigen.values .= max.(sigma_eigen.values, 1e-14) # To ensure full-rank numerically

    dim = size(sigma, 1)
    dim2 = dim^2
    fw = LSIFramework(dim)

    # sigma^1/4 in sigma eigenspace
    @. fw.temp1 = sigma_eigen.vectors * (sigma_eigen.values' ^ 0.25)
    mul!(fw.sigma_quarter, fw.temp1, sigma_eigen.vectors')

    # sigma^1/2 in sigma eigenspace
    @. fw.temp1 = sigma_eigen.vectors * (sigma_eigen.values' ^ 0.5)
    mul!(fw.sigma_half, fw.temp1, sigma_eigen.vectors')

    # log(sigma) in sigma eigenspace
    @. fw.temp1 = sigma_eigen.vectors * log(sigma_eigen)
    mul!(fw.sigma_log, fw.temp1, sigma_eigen.vectors')

    function fg!(cost, gradient, A_flat)
        n_params = length(A_flat) / 2
        A_real = unsafe_wrap(Array, pointer(A_flat), (d, d))
        A_imag = unsafe_wrap(Array, pointer(A_flat, n_params * sizeof(Float64) + 1), (d, d))

        @. fw.A = Complex(A_real, A_imag)
        mul!(fw.AdagA, fw.A', fw.A)

        # for log stability
        @inbounds for i in i:dim
            fw.AdagA[i, i] += 1e-12
        end

        # Y = sigma^1/4 * A' A  * sigma^1/4
        mul!(fw.temp1, fw.AdagA, fw.sigma_quarter)
        mul!(fw.sigma_AdagA, fw.sigma_quarter, fw.temp1)

        # Normalize (wlog)
        norm_val = real(dot(fw.sigma_AdagA, fw.sigma_AdagA))
        scale = 1.0 / sqrt(norm_val)
        rmul!(fw.AdagA, scale)
        rmul!(fw.sigma_AdagA, scale)

        # +++ Cost function: dirichlet / entropy

        # --- Relative entropy D(rho || sigma)
        # Faster to work in eigenspace of rho:
        copyto!(fw.temp2, fw.sigma_AdagA)
        eigen_sigma_AdagA = eigen!(Hermitian(fw.temp2))
        
        # Term 1: Tr(rho ln rho) : rho = sigma_AdagA^2
        entropy_term1 = 0.0
        @inbounds for eigval in eigen_sigma_AdagA.values
            if eigval > 1e-12
                eigval_sq = eigval^2
                term_1 += eigval_sq * log(eigval_sq)
            end
        end

        mul!(fw.temp3, fw.sigma_AdagA, fw.sigma_AdagA)
        entropy_term2 = real(dot(fw.temp3, fw.sigma_log))

        relative_entropy = entropy_term1 - entropy_term2

        # --- Dirichlet form
        copyto!(fw.AdagA_vec, fw.AdagA)
        mul!(fw.L_AdagA_vec, L_mat, fw.AdagA_vec)  # L_AdagA_vec = L(A'A)

        # <X, L(X)>_sigma := Tr(sigma^1/2 X sigma^1/2 L(X)) : X = A'A
        mul!(fw.temp2, fw.AdagA, fw.sigma_half)
        mul!(fw.temp3, fw.sigma_half, fw.temp2)  # temp3 = sigma^1/2 A'A sigma^1/2

        dirichlet = -real(dot(fw.temp3, fw.L_AdagA_vec))

        # +++
        if cost !== nothing
            if abs(relative_entropy) < 1e-9  # Avoid entropy 0
                return 1000.0
            end
            lsi_value =  energy / entropy  # LSI constant
        end

        # +++ Gradient = (grad dirichlet - (dirichlet / entropy^2) *  grad entropy) / entropy
        if gradient !== nothing
            inverse_entropy = 1.0 / relative_entropy
            ratio = dirichlet / relative_entropy^2
        
            # --- Entropy gradient = sigma^1/4 [2Y(ln(Y^2)+1) - {Y, ln sigma}] sigma^1/4 : Y = sigma^1/4 A'A sigma^1/4
            U = eigen_sigma_AdagA.vectors  # To work in eigenbasis of Y

            fill!(fw.temp1, 0.0)  # temp1 to store: U * (diagonal core)
            @inbounds for j in 1:dim  # U * (diagonal core)
                eignval_j = eigen_sigma_AdagA.values[j]
                value = (eignval_j > 1e-12) ? (2 * eignval_j * (2 * log(eignval_j) + 1)) : 0.0
                for i in 1:dim
                    fw.temp1[i, j] = U[i, j] * value
                end
            end

            mul!(fw.temp2, fw.temp1, U')  # temp2 = U (diagonal core) U' = part 1 of entropy gradient

            mul!(fw.temp1, fw.AdagA, fw.sigma_log)  # temp1 = Y ln sigma; temp1' = ln sigma Y (both Hermitian)
            @. fw.temp2 = fw.temp2 - (fw.temp1 + fw.temp1')  # temp2 = Entropy gradient without sigma sandwich

            # Sandwich in sigmas
            mul!(fw.temp1, fw.temp2, fw.sigma_quarter)
            mul!(fw.gradient, fw.sigma_quarter, fw.temp1)  # fw.gradient = Entropy gradient

            # --- Dirichlet gradient = -0.5 * [sigma^1/2 L(X) sigma^1/2 + L'(sigma^1/2 X sigma^1/2)]
            L_AdagA_mat = reshape(fw.L_AdagA_vec, dim2, dim2)
            mul!(fw.temp1, L_AdagA_mat, fw.sigma_half)
            mul!(fw.temp2, fw.sigma_half, fw.temp1)  # temp2 = sigma^1/2 L(X) sigma^1/2

            mul!(fw.temp3, fw.AdagA, fw.sigma_half)
            mul!(fw.temp1, fw.sigma_half, fw.temp3)  # temp1 = sigma^1/2 X sigma^1/2  #TODO: Can't we write a Gamma function to sandwich nicer?

            mul!(fw.L_AdagA_vec, L_mat', fw.temp1)  # = L^\dagger(sigma^1/2 A'A sigma^1/2)
            Ldag_AdagA_mat = reshape(fw.L_AdagA_vec, dim2, dim2)

            # Combine #?
            @inbounds for i in 1:length(fw.gradient)
            end

        end
    end
end