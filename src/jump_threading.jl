# Shared threading policy and ordered frequency work-list construction for jump
# operators. This file is included before every construction and runtime caller.

# Minimum number of frequency labels to enable omega-loop parallelism.
# Benchmarks on four Julia threads place the task-overhead crossover at ten.
const OMEGA_THREAD_THRESHOLD = 10

"""
    _partition_range(range, n_chunks) -> Vector{UnitRange{Int}}

Partition a range into approximately equal chunks for parallel execution.
"""
function _partition_range(range::UnitRange{Int}, n_chunks::Int)
    len = length(range)
    n_chunks = min(n_chunks, len)
    base = div(len, n_chunks)
    remainder = rem(len, n_chunks)
    chunks = Vector{UnitRange{Int}}(undef, n_chunks)
    start = first(range)
    for i in 1:n_chunks
        chunk_size = base + (i <= remainder ? 1 : 0)
        chunks[i] = start:(start + chunk_size - 1)
        start += chunk_size
    end
    return chunks
end

# Build flat `(jump_idx, label_idx)` work in jump-major, label-minor order.
# Hermitian jumps use the non-positive half-grid and reconstruct the explicit
# negative-frequency partner in the caller; non-Hermitian jumps use every
# signed label.
function _populate_jump_frequency_work_list!(
    work::Vector{Tuple{Int, Int}},
    jump_hermitian::Vector{Bool},
    energy_labels::AbstractVector{<:Real},
)
    n_jumps = length(jump_hermitian)
    n_labels = length(energy_labels)
    empty!(work)
    sizehint!(work, n_jumps * n_labels)
    @inbounds for k in 1:n_jumps
        is_herm = jump_hermitian[k]
        for li in 1:n_labels
            if is_herm && energy_labels[li] > 1e-12
                continue
            end
            push!(work, (k, li))
        end
    end
    return work
end
