"""
    KrylovWorkspace{T, PD}

Pre-allocated workspace for matrix-free Lindbladian action (matvec).

Stores precomputed data (transition function, energy labels, gamma normalization),
the coherent correction matrix B (if applicable), jump operator references,
and scratch matrices for zero-allocation dissipator accumulation.

Tied to a specific (config, hamiltonian) pair at construction time.

# Type Parameters
- `T <: Complex`: element type of all matrices (e.g. `ComplexF64`)
- `PD <: NamedTuple`: concrete type of the precomputed data tuple (varies by domain)

# Fields
- `precomputed_data::PD`: from `_precompute_data(config, ham_or_trott)`
- `B_total::Union{Nothing, Matrix{T}}`: precomputed coherent B (nothing for GNS / with_coherent=false)
- `jumps::Vector{JumpOp}`: reference to jump operators

Scratch matrices (all dim x dim, zeroed or overwritten each matvec call):
- `jump_oft::Matrix{T}`: A(omega) buffer (written by `oft!`)
- `tmp1::Matrix{T}`: scratch for `mul!` results
- `tmp2::Matrix{T}`: scratch for `mul!` results
- `LdagL::Matrix{T}`: scratch for L'*L product
- `rho_out::Matrix{T}`: output accumulator (zeroed at start of each matvec)
"""
struct KrylovWorkspace{T<:Complex, PD<:NamedTuple}
    # Precomputed data (immutable after construction)
    precomputed_data::PD
    B_total::Union{Nothing, Matrix{T}}
    jumps::Vector{JumpOp}

    # Scratch matrices for dissipator accumulation (dim x dim)
    jump_oft::Matrix{T}
    tmp1::Matrix{T}
    tmp2::Matrix{T}
    LdagL::Matrix{T}
    rho_out::Matrix{T}
end

"""
    KrylovWorkspace(config, hamiltonian, jumps; trotter=nothing)

Construct a `KrylovWorkspace` pre-allocating all scratch matrices for the given
(config, hamiltonian) pair. Mirrors `construct_lindbladian` setup in `furnace.jl`.

# Arguments
- `config::AbstractLiouvConfig`: Lindbladian configuration (EnergyDomain, TimeDomain, etc.)
- `hamiltonian::HamHam`: Hamiltonian with eigenbasis data
- `jumps::Vector{JumpOp}`: Jump operators (stored by reference)
- `trotter::Union{TrottTrott, Nothing}=nothing`: Trotter object (required for TrotterDomain)
"""
function KrylovWorkspace(
    config::AbstractLiouvConfig,
    hamiltonian::HamHam,
    jumps::Vector{JumpOp};
    trotter::Union{TrottTrott, Nothing}=nothing,
)
    # Determine ham_or_trott (mirrors construct_lindbladian in furnace.jl)
    ham_or_trott = if config.domain isa TrotterDomain
        trotter === nothing && error("A Trotter object must be provided for the TrotterDomain")
        trotter
    else
        hamiltonian
    end

    # Precompute domain-specific data (transition, energy_labels, gamma_norm_factor, ...)
    precomputed_data = _precompute_data(config, ham_or_trott)

    # Precompute coherent B_total (returns nothing for GNS / with_coherent=false)
    B_total = _precompute_coherent_total_B(jumps, ham_or_trott, config, precomputed_data)

    # Determine dimensions and element type
    dim = size(hamiltonian.data, 1)
    CT = Complex{eltype(hamiltonian.eigvals)}

    # Allocate scratch matrices
    jump_oft = zeros(CT, dim, dim)
    tmp1     = zeros(CT, dim, dim)
    tmp2     = zeros(CT, dim, dim)
    LdagL    = zeros(CT, dim, dim)
    rho_out  = zeros(CT, dim, dim)

    return KrylovWorkspace{CT, typeof(precomputed_data)}(
        precomputed_data, B_total, jumps,
        jump_oft, tmp1, tmp2, LdagL, rho_out,
    )
end
