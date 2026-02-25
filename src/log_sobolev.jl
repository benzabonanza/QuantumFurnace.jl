# using LinearAlgebra
# using Optim
# using Printf

# """
#     _sandwich!(target, middle, bread, cache)

#     Computes `target = bread * middle * bread` in-place.
#     Helps with all the KMS sigma conjugations.
# """
# function _sandwich!(target::AbstractMatrix, middle::AbstractMatrix, bread::AbstractMatrix, cache::AbstractMatrix)
#     mul!(cache, middle, bread)  # Cache = middle * bread
#     mul!(target, bread, cache)  # Target = bread * middle * bread
#     return target
# end

# #TODO: Rewrite this with apply_lindbladian!()
# function compute_LSI_alpha2(result::LindbladianResult; n_restarts=3, g_tol=1e-5)

#     # In respective eigenspace either (H or Trotter)
#     L_mat = result.liouvillian
#     sigma = Hermitian(result.fixed_point)
#     sigma_eigen = eigen(sigma)
#     sigma_eigen.values .= max.(sigma_eigen.values, 1e-14) # To ensure full-rank numerically

#     dim = size(sigma, 1)
#     fw = LSIFramework(dim)

#     # sigma^1/4 in sigma eigenspace
#     @. fw.temp1 = sigma_eigen.vectors * (sigma_eigen.values' ^ 0.25)
#     mul!(fw.sigma_quarter, fw.temp1, sigma_eigen.vectors')

#     # sigma^1/2 in sigma eigenspace
#     @. fw.temp1 = sigma_eigen.vectors * (sigma_eigen.values' ^ 0.5)
#     mul!(fw.sigma_half, fw.temp1, sigma_eigen.vectors')

#     # log(sigma) in sigma eigenspace
#     @. fw.temp1 = sigma_eigen.vectors * log(sigma_eigen)
#     mul!(fw.sigma_log, fw.temp1, sigma_eigen.vectors')

#     function fg!(cost, gradient, A_flat)
#         n_params = length(A_flat) / 2
#         # The next lines are like reshape(view(A_flat, ...)) but they return Array type that BLAS operations like.
#         A_real = unsafe_wrap(Array, pointer(A_flat), (d, d))
#         A_imag = unsafe_wrap(Array, pointer(A_flat, n_params * sizeof(Float64) + 1), (d, d))

#         @. fw.A = Complex(A_real, A_imag)
#         mul!(fw.AdagA, fw.A', fw.A)

#         # for log stability
#         @inbounds for i in 1:dim
#             fw.AdagA[i, i] += 1e-12
#         end

#         # Gamma2_AdagA = sigma^1/4 * A' A  * sigma^1/4
#         _sandwich!(fw.Gamma2_AdagA, fw.AdagA, fw.sigma_quarter, fw.temp1)

#         # Normalize (wlog)
#         norm_val = real(dot(fw.Gamma2_AdagA, fw.Gamma2_AdagA))
#         scale = 1.0 / sqrt(norm_val)
#         rmul!(fw.AdagA, scale)
#         rmul!(fw.Gamma2_AdagA, scale)

#         # +++ Cost function: dirichlet / entropy

#         # --- Relative entropy D(rho || sigma)
#         # Faster to work in eigenspace of rho:
#         copyto!(fw.temp2, fw.Gamma2_AdagA)
#         eigen_Gamma2_AdagA = eigen!(Hermitian(fw.temp2))
        
#         # Term 1: Tr(rho ln rho) : rho = Gamma2_AdagA^2
#         entropy_term1 = 0.0
#         @inbounds for eigval in eigen_Gamma2_AdagA.values
#             if eigval > 1e-12
#                 eigval_sq = eigval^2
#                 entropy_term1 += eigval_sq * log(eigval_sq)
#             end
#         end

#         mul!(fw.temp3, fw.Gamma2_AdagA, fw.Gamma2_AdagA)
#         entropy_term2 = real(dot(fw.temp3, fw.sigma_log))

#         relative_entropy = entropy_term1 - entropy_term2

#         # --- Dirichlet form
#         copyto!(fw.AdagA_vec, fw.AdagA)
#         mul!(fw.L_AdagA_vec, L_mat, fw.AdagA_vec)  # L_AdagA_vec = L(A'A)

#         # <X, L(X)>_sigma := Tr(sigma^1/2 X sigma^1/2 L(X)) : X = A'A
#         _sandwich!(fw.temp3, fw.AdagA, fw.sigma_half, fw.temp2)  # temp3 = \sigma^1/2 A'A \sigma^1/2

#         dirichlet = -real(dot(fw.temp3, fw.L_AdagA_vec))

#         # +++
#         if cost !== nothing
#             if abs(relative_entropy) < 1e-9  # Avoid entropy 0
#                 return 1000.0
#             end
#             return dirichlet / relative_entropy  # LSI constant
#         end

#         # +++ Gradient = (grad dirichlet - (dirichlet / entropy^2) *  grad entropy) / entropy
#         if gradient !== nothing
#             inverse_entropy = 1.0 / relative_entropy
#             ratio = dirichlet / (relative_entropy^2)
        
#             # --- Entropy gradient = sigma^1/4 [2Y(ln(Y^2)+1) - {Y, ln sigma}] sigma^1/4 : Y = sigma^1/4 A'A sigma^1/4
#             U = eigen_Gamma2_AdagA.vectors  # To work in eigenbasis of Y

#             fill!(fw.temp1, 0.0)  # temp1 to store: U * (diagonal core)
#             @inbounds for j in 1:dim  # U * (diagonal core)
#                 eignval_j = eigen_Gamma2_AdagA.values[j]
#                 value = (eignval_j > 1e-12) ? (2 * eignval_j * (2 * log(eignval_j) + 1)) : 0.0
#                 for i in 1:dim
#                     fw.temp1[i, j] = U[i, j] * value
#                 end
#             end

#             mul!(fw.temp2, fw.temp1, U')  # temp2 = U (diagonal core) U' = part 1 of entropy gradient

#             mul!(fw.temp1, fw.AdagA, fw.sigma_log)  # temp1 = Y ln sigma; temp1' = ln sigma Y (both Hermitian)
#             @. fw.temp2 = fw.temp2 - (fw.temp1 + fw.temp1')  # temp2 = Entropy gradient without sigma sandwich

#             # Sandwich in sigmas
#             _sandwich!(fw.gradient, fw.temp2, fw.sigma_quarter, fw.temp1)  # fw.gradient = Entropy gradient

#             # --- Dirichlet gradient = -0.5 * [sigma^1/2 L(X) sigma^1/2 + L'(sigma^1/2 X sigma^1/2)]
#             L_AdagA_mat = reshape(fw.L_AdagA_vec, dim, dim)
#             _sandwich!(fw.temp2, L_AdagA_mat, fw.sigma_half, fw.temp1)  # temp2 = sigma^1/2 L(X) sigma^1/2

#             mul!(fw.L_AdagA_vec, L_mat', fw.temp3)  # = L^\dagger(sigma^1/2 A'A sigma^1/2)
#             Ldag_AdagA_mat = reshape(fw.L_AdagA_vec, dim, dim)

#             # Combine
#             factor_dirichlet = -0.5 * inverse_entropy
#             factor_entropy = -ratio

#             @inbounds for i in 1:length(fw.gradient)
#                 gradient_dirichlet = real(fw.temp2[i] + Ldag_AdagA_mat[i])
#                 gradient_entropy = real(fw.gradient[i])
#                 fw.gradient[i] = factor_dirichlet * gradient_dirichlet + factor_entropy * gradient_entropy  # Grad_(A'A)
#             end

#             # Chain rule for A (the optimized over variable): Grad_A = 2 * A * Grad_(A'A)
#             mul!(fw.temp1, fw.A, fw.gradient)

#             # Write back to flat gradient
#             gradient_real = unsafe_wrap(Array, pointer(gradient), (dim, dim))
#             gradient_imag = unsafe_wrap(Array, pointer(gradient, n_params * sizeof(Float64) + 1), (dim, dim))

#             @inbounds for i in 1:dim^2
#                 val = 2 * fw.temp1[i]
#                 gradient_real[i] = real(val)
#                 gradient_imag[i] = imag(val)
#             end
#         end
#     end

#     # +++ Optimization
#     optim_options = Optim.Options(
#         g_tol = g_tol,
#         iterations = 1000,
#         store_trace = false,
#         show_trace = false
#     )

#     x0 = Vector{Float64}(undef, 2 * dim^2)
#     best_LSI_alpha2 = Inf

#     for i in 0:n_restarts

#         # Start in the spectral gap valley once, to see where that leads to.
#         # Every random initial state in this convex valley would lead to the same LSI constant.
#         if i == 0  
#             v2_mat = result.gap_mode
#             hermitianize!(v2_mat)
#             v2_hermitian = Hermitian(v2_mat)
#             v2_normalized = v2_hermitian / opnorm(v2_hermitian)

#             # X -> A : X = A'A (Cholesky parametrization to avoid constrained optim. with just X)
#             X_initial = Matrix{ComplexF64}(I, dim, dim) + 0.1 * v2_normalized  # close to equilibrium init.
#             X_eigen = eigen(Hermitian(X_initial))
#             X_eigen.values .= max.(X_eigen.values, 1e-14)  # Numerical stability

#             A_initial = X_eigen.vectors * Diagonal(sqrt.(X_eigen.values)) * X_eigen.vectors'

#             # Flatten
#             x0[1:dim^2] .= real(vec(A_initial))
#             x0[dim^2+1:end] .= imag(vec(A_initial))
#         else  # Random initialization
#             rand!(x0)  # Gaussian random matrix -> makes X = A'A a Wishart matrix.
#             rmul!(x0, 0.5)  # Scaled down a bit.
#         end

#         # --- Run optimizer
#         optim_result = optimize(Optim.only_fg!(fg!), x0, LBFGS(), optim_options)
#         LSI_alpha2 = Optim.minimum(optim_result)

#         if LSI_alpha2 < best_LSI_alpha2
#             best_LSI_alpha2 = LSI_alpha2
#         end

#         @printf "Run %d: alpha = %.4f (Converged: %s)\n" i LSI_alpha2 Optim.converged(optim_result)
#     end

#     return best_LSI_alpha2
# end