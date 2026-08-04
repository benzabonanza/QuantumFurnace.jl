# =============================================================================
# Staged code (qf-6z9.4): deprecated direct OFT routines + their cache struct.
#
# Why staged: `time_oft!` and `trotter_oft!` were the original direct-summation
# Operator Fourier Transform routines for TimeDomain and TrotterDomain. Mainline
# code paths replaced them with the NUFFT-based prefactor pipeline
# (`_prepare_oft_nufft_prefactors`), which is faster and produces the same
# values (cross-validated to ~1e-10 op-norm in the regression suite, since
# removed). They were kept for personal reference and future reproducibility.
#
# Why commented out: the file is excluded from the active `include()` chain.
# Reattaching requires (i) uncommenting, (ii) re-exporting from
# `src/QuantumFurnace.jl` (`time_oft!`, `trotter_oft!`,
# `compute_oft_trotter_error`, `compute_oft_trotter_error_all_jumps`),
# (iii) restoring the `OFTCaches` struct definition, and (iv) re-introducing
# the corresponding regression tests in `test/test_dm_scaling.jl` (DMTST-06
# and the ofts-dependent half of DMTST-06b).
# =============================================================================

# # Became obsolete with NUFFTCaches. But used for debugging.
# struct OFTCaches{T<:AbstractFloat}
#     prefactors::Vector{Complex{T}}
#     U::Diagonal{Complex{T}, Vector{Complex{T}}}
#     temp_op::Matrix{Complex{T}}
#
#     function OFTCaches{T}(dim::Int) where {T<:AbstractFloat}
#         CT = Complex{T}
#         prefactors = zeros(CT, 0) # Will be resized later
#         U = Diagonal(zeros(CT, dim))
#         temp_op = zeros(CT, dim, dim)
#         new{T}(prefactors, U, temp_op)
#     end
# end
#
# # Depricated but used for tests. We use precomputed NUFFT prefactors now in jump_contributions!()
# function time_oft!(
#     out_matrix::Matrix{<:Complex},
#     caches::OFTCaches,
#     jump::JumpOp,
#     energy::Real,
#     hamiltonian::HamHam,
#     time_labels::AbstractVector{<:Real},
#     sigma::Real;
#     filter::AbstractFilter = GaussianFilter(sigma),
#     )
#
#     # Ensure the prefactor cache is the right size
#     if length(caches.prefactors) != length(time_labels)
#         resize!(caches.prefactors, length(time_labels))
#     end
#
#     # In-place calculation of prefactors via the filter abstraction.
#     # Explicit loop (vs `@fastmath @.`) so DLL's complex time_kernel works.
#     @inbounds for k in eachindex(time_labels)
#         caches.prefactors[k] = time_kernel(filter, time_labels[k]) * cis(-energy * time_labels[k])
#     end
#
#     zero_index = findfirst(t -> t >= -1e-12, time_labels)
#
#     # Zero out the output matrix before we start accumulating
#     fill!(out_matrix, 0.0)
#
#     # --- Re-use the cache matrices U and temp_op inside the loops ---
#     if jump.orthogonal  # Orthogonal (X, Z)
#         # t = 0.0 case: U = I
#         @. out_matrix += caches.prefactors[zero_index] * jump.in_eigenbasis
#
#         for i in (zero_index + 1):length(time_labels)
#             t = time_labels[i]
#             @fastmath caches.U.diag .= exp.(1im .* hamiltonian.eigvals .* t)
#
#             copyto!(caches.temp_op, jump.in_eigenbasis)
#             # temp_op = U*jump*U', for diagonal U's:
#             caches.temp_op .*= (caches.U.diag * caches.U.diag')
#
#             @. out_matrix += caches.prefactors[i] * caches.temp_op
#             @. out_matrix += conj(caches.prefactors[i]) * $(transpose(caches.temp_op))  # We learnt: @. makes transpose.()
#         end
#     else  # Non-orthogonal (Y)
#         for i in eachindex(time_labels)
#             t = time_labels[i]
#             @fastmath caches.U.diag .= exp.(1im .* hamiltonian.eigvals .* t)
#
#             mul!(caches.temp_op, caches.U, jump.in_eigenbasis)  # temp = U * jump
#             mul!(out_matrix, caches.temp_op, caches.U', caches.prefactors[i], 1.0) # out += prefactor * U*jump*U'
#         end
#     end
#
#     return out_matrix
# end
#
# function trotter_oft!(
#     out_matrix::Matrix{<:Complex},
#     caches::OFTCaches,
#     jump::JumpOp,
#     energy::Real,
#     trotter::TrottTrott,
#     time_labels::AbstractVector{<:Real},
#     sigma::Real;
#     filter::AbstractFilter = GaussianFilter(sigma),
#     )
#
#     if length(caches.prefactors) != length(time_labels)
#         resize!(caches.prefactors, length(time_labels))
#     end
#
#     # In-place calculation of prefactors via the filter abstraction.
#     # Explicit loop (vs `@fastmath @.`) so DLL's complex time_kernel works.
#     @inbounds for k in eachindex(time_labels)
#         caches.prefactors[k] = time_kernel(filter, time_labels[k]) * cis(-energy * time_labels[k])
#     end
#
#     zero_index = findfirst(t -> t >= -1e-12, time_labels)
#
#     fill!(out_matrix, 0.0)
#
#     if jump.orthogonal
#         # t = 0.0 case: U = I
#         @fastmath out_matrix .+= caches.prefactors[zero_index] .* jump.in_eigenbasis
#
#         for i in (zero_index + 1):length(time_labels)
#             num_t0_steps = i - zero_index
#             @fastmath caches.U.diag .= trotter.eigvals_t0 .^ num_t0_steps
#
#             copyto!(caches.temp_op, jump.in_eigenbasis)
#             # temp_op = U*jump*U', for diagonal U's:
#             caches.temp_op .*= (caches.U.diag * caches.U.diag')
#
#             # Accumulate both terms in-place
#             LinearAlgebra.axpby!(caches.prefactors[i], caches.temp_op, 1.0, out_matrix)
#             LinearAlgebra.axpby!(conj(caches.prefactors[i]), transpose(caches.temp_op), 1.0, out_matrix)
#         end
#     else # Non-orthogonal jumps
#         for i in eachindex(time_labels)
#             num_t0_steps = i - zero_index
#             @fastmath caches.U.diag .= trotter.eigvals_t0 .^ num_t0_steps
#
#             # temp_op = U * jump.in_eigenbasis
#             mul!(caches.temp_op, caches.U, jump.in_eigenbasis)
#
#             # out_matrix += prefactor * (temp_op * U')
#             mul!(out_matrix, caches.temp_op, caches.U', caches.prefactors[i], 1.0)
#         end
#     end
#
#     return out_matrix
# end
