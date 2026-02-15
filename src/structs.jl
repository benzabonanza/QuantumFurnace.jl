# Domains
abstract type AbstractDomain end

struct BohrDomain <: AbstractDomain end
struct EnergyDomain <: AbstractDomain end
struct TimeDomain <: AbstractDomain end
struct TrotterDomain <: AbstractDomain end

# Became obsolete with NUFFTCaches. But used for debugging.
struct OFTCaches
    prefactors::Vector{ComplexF64}
    U::Diagonal{ComplexF64, Vector{ComplexF64}}
    temp_op::Matrix{ComplexF64}
    
    function OFTCaches(dim::Int)
        prefactors = zeros(ComplexF64, 0) # Will be resized later
        U = Diagonal(zeros(ComplexF64, dim))
        temp_op = zeros(ComplexF64, dim, dim)
        new(prefactors, U, temp_op)
    end
end

"""Workspace for building a dense Liouvillian matrix with minimal allocations.

Used by `construct_liouvillian` when accumulating the full vectorized Lindbladian
(`dim^2 × dim^2`).

Buffers:
  - `Id`: identity on the system Hilbert space.
  - `jump_tmp`: generic scratch (e.g. OFT output, alpha-weighted jump).
  - `jump_conj`: scratch for elementwise conjugate of `jump_tmp`.
  - `jump_dag_jump`: scratch for `jump_tmp' * jump_tmp`.
  - `jump2_jump1`: scratch for `jump_2 * jump_1` in mixed (Bohr) note.
"""
struct LindbladianWorkspace
    Id::Matrix{ComplexF64}
    jump_tmp::Matrix{ComplexF64}
    jump_conj::Matrix{ComplexF64}
    jump_dag_jump::Matrix{ComplexF64}
    jump2_jump1::Matrix{ComplexF64}

    function LindbladianWorkspace(dim::Int)
        Id = Matrix{ComplexF64}(I, dim, dim)
        jump_tmp = zeros(ComplexF64, dim, dim)
        jump_conj = zeros(ComplexF64, dim, dim)
        jump_dag_jump = zeros(ComplexF64, dim, dim)
        jump2_jump1 = zeros(ComplexF64, dim, dim)
        new(Id, jump_tmp, jump_conj, jump_dag_jump, jump2_jump1)
    end
end

abstract type AbstractConfig{D<:AbstractDomain} end
abstract type AbstractLiouvConfig{D<:AbstractDomain} <: AbstractConfig{D} end
abstract type AbstractThermalizeConfig{D<:AbstractDomain} <: AbstractConfig{D} end

"""
    _build_common_fields(; kwargs...) -> NamedTuple

Shared constructor helper for all 4 config types. Validates and returns
the common fields as a NamedTuple. Centralizes default logic so that
individual struct definitions stay DRY in construction.
"""
function _build_common_fields(;
    num_qubits::Int64,
    with_coherent::Bool,
    with_linear_combination::Bool,
    beta::Float64,
    sigma::Float64,
    gaussian_parameters::Union{Tuple{Float64, Float64}, Tuple{Nothing, Nothing}} = (nothing, nothing),
    a::Union{Float64, Nothing} = nothing,
    b::Union{Float64, Nothing} = nothing,
    num_energy_bits::Union{Int, Nothing} = nothing,
    t0::Union{Float64, Nothing} = nothing,
    w0::Union{Float64, Nothing} = nothing,
    eta::Union{Float64, Nothing} = nothing,
    num_trotter_steps_per_t0::Union{Int, Nothing} = nothing,
)
    return (; num_qubits, with_coherent, with_linear_combination,
            beta, sigma, gaussian_parameters, a, b,
            num_energy_bits, t0, w0, eta, num_trotter_steps_per_t0)
end

# Let's keep this structure, and have the "give w0 for desired energy integral error" type of config optimization
# before the construct_liouvillian function
"""
    LiouvConfig

    A configuration object that holds all the parameters for the core function: `run_liouvillian`, which constructs the Lindbladian of the thermalizing system.

    # Fields
    - `num_qubits`: The number of system qubits.
    - `with_coherent`: The option to add (=true) or omit (=false) the coherent term in the Lindbladian.\nIf added, the target state of the evolution will be the exactly the Gibbs state, otherwise only approximately.
    - `with_linear_combination`: The option to choose if we want to apply a convex combination of Lindbladians for a faster mixing. Could add extra complexities if the resulting transition function is not smooth. (See more in `Theory`).
    - `a` and `b`: The parameters that specify the type of linear combination.
    - `eta`: in the case of the Metropolis linear combination, η is an additional coefficient that determines the accuracy of the time domain approximation.
    - `domain`: The domain the simulation runs in (`BOHR`, `ENERGY`, `TIME`, `TROTTER`). The choice of the domain represents the levels of approximations we need to get form theory down to quantum circuitry.
    - `num_energy_bits`: Determines the how coarse the energy and time grid is and thus how accurate the approximations between each domain are.
    - `t0` and `w0`: are the time and energy units we are working with in the Riemann summed integrals. Of course, the smaller the better but also the costlier, and the two are intertwined due to Fourier: ω₀t₀ = 2π / N.
    - `num_trotter_steps_per_t0`: The number of Trotter steps used for a unit of time t₀.

    ## Currently possible linear combinations:
    (a, b) =
    - (0, 0) - no linear combination, simple Gaussian
    - (>0, 0) - linear combination that results in Metropolis-like transition
    - (>0, >0) - linear combination that results in Glauber transition (smoother)

    ## Available domains:
    The `domain` field can be set to one of the following options:
    - **`BohrDomain()`**: The highest level domain where the jump operators and thus the Lindbladian are written in a decomposition of Bohr frequencies.
    - **`EnergyDomain()`**: A level lower, in which the operators are approximated by energy integrals.
    - **`TimeDomain()`**: Another level lower, in which the energy approximates are written up as Fourier's of the temporal equals.
    - **`TrotterDomain()`**: The lowest level, thus also the only one implementable on a quantum computer, in which all time evolutions are replaced via their Trotter series.
"""
@kwdef struct LiouvConfig{D <: AbstractDomain} <: AbstractLiouvConfig{D}
    num_qubits::Int64
    with_coherent::Bool
    with_linear_combination::Bool
    domain::D
    beta::Float64
    sigma::Float64
    gaussian_parameters::Union{Tuple{Float64, Float64}, Tuple{Nothing, Nothing}} = (nothing, nothing)  # (ω_γ, σ_γ)
    a::Union{Float64, Nothing} = nothing
    b::Union{Float64, Nothing} = nothing
    num_energy_bits::Union{Int, Nothing} = nothing
    t0::Union{Float64, Nothing} = nothing
    w0::Union{Float64, Nothing} = nothing
    eta::Union{Float64, Nothing} = nothing
    num_trotter_steps_per_t0::Union{Int, Nothing} = nothing
end
"""
    LiouvConfigGNS

    Configuration for Liouvillian construction for Chen's **approx. GNS-detailed-balance** Lindbladian.

    This is the "GNS-DB" line: it uses the **unshifted** transition weight \tilde{γ}(ω) (KMS-conditioned),
    and (by design) **omits** the coherent correction term `B`.

    Fields are shared with `LiouvConfig`.
"""
struct LiouvConfigGNS{D <: AbstractDomain} <: AbstractLiouvConfig{D}
    num_qubits::Int64
    with_coherent::Bool
    with_linear_combination::Bool
    domain::D
    beta::Float64
    sigma::Float64
    gaussian_parameters::Union{Tuple{Float64, Float64}, Tuple{Nothing, Nothing}}
    a::Union{Float64, Nothing}
    b::Union{Float64, Nothing}
    num_energy_bits::Union{Int, Nothing}
    t0::Union{Float64, Nothing}
    w0::Union{Float64, Nothing}
    eta::Union{Float64, Nothing}
    num_trotter_steps_per_t0::Union{Int, Nothing}

    function LiouvConfigGNS{D}(
        num_qubits::Int64, with_coherent::Bool, with_linear_combination::Bool, domain::D,
        beta::Float64, sigma::Float64,
        gaussian_parameters::Union{Tuple{Float64, Float64}, Tuple{Nothing, Nothing}},
        a::Union{Float64, Nothing}, b::Union{Float64, Nothing},
        num_energy_bits::Union{Int, Nothing}, t0::Union{Float64, Nothing},
        w0::Union{Float64, Nothing}, eta::Union{Float64, Nothing},
        num_trotter_steps_per_t0::Union{Int, Nothing}
    ) where {D}
        with_coherent && error("GNS configs must have with_coherent=false")
        new{D}(num_qubits, with_coherent, with_linear_combination, domain,
               beta, sigma, gaussian_parameters, a, b,
               num_energy_bits, t0, w0, eta, num_trotter_steps_per_t0)
    end
end

# Keyword constructor for LiouvConfigGNS (replaces @kwdef)
function LiouvConfigGNS(;
    num_qubits::Int64,
    with_coherent::Bool = false,
    with_linear_combination::Bool,
    domain::D,
    beta::Float64,
    sigma::Float64,
    gaussian_parameters::Union{Tuple{Float64, Float64}, Tuple{Nothing, Nothing}} = (nothing, nothing),
    a::Union{Float64, Nothing} = nothing,
    b::Union{Float64, Nothing} = nothing,
    num_energy_bits::Union{Int, Nothing} = nothing,
    t0::Union{Float64, Nothing} = nothing,
    w0::Union{Float64, Nothing} = nothing,
    eta::Union{Float64, Nothing} = nothing,
    num_trotter_steps_per_t0::Union{Int, Nothing} = nothing,
) where {D <: AbstractDomain}
    LiouvConfigGNS{D}(num_qubits, with_coherent, with_linear_combination, domain,
                      beta, sigma, gaussian_parameters, a, b,
                      num_energy_bits, t0, w0, eta, num_trotter_steps_per_t0)
end

"""
    ThermalizeConfig

    Configuration for the thermalization process, that emulates the quantum algorithm step-by-step.

    Inherits core physical parameters from the logic in [`LiouvConfig`](@ref), but includes
    simulation-specific settings, e.g. `mixing_time` and `delta`

    # Specific Fields
    - `mixing_time`: Total duration of the time evolution.
    - `delta`: Time step size for the weak-measurement emulation.
    """
@kwdef struct ThermalizeConfig{D <: AbstractDomain}  <: AbstractThermalizeConfig{D}
    num_qubits::Int64
    with_coherent::Bool
    with_linear_combination::Bool
    domain::D
    beta::Float64
    sigma::Float64
    gaussian_parameters::Union{Tuple{Float64, Float64}, Tuple{Nothing, Nothing}} = (nothing, nothing)  # (ω_γ, σ_γ)
    a::Union{Float64, Nothing} = nothing
    b::Union{Float64, Nothing} = nothing
    num_energy_bits::Union{Int, Nothing} = nothing
    t0::Union{Float64, Nothing} = nothing
    w0::Union{Float64, Nothing} = nothing
    eta::Union{Float64, Nothing} = nothing
    num_trotter_steps_per_t0::Union{Int, Nothing} = nothing

    # For thermalization the configs:
    mixing_time::Float64
    delta::Float64
end

"""
    ThermalizeConfigGNS

    Configuration for the step-by-step weak-measurement thermalization emulation for
    Chen's **approx. GNS-detailed-balance** Lindbladian.

    This line uses the unshifted transition weight and (by design) omits the coherent term `B`.

    Fields are shared with `ThermalizeConfig`.
"""
struct ThermalizeConfigGNS{D <: AbstractDomain} <: AbstractThermalizeConfig{D}
    num_qubits::Int64
    with_coherent::Bool
    with_linear_combination::Bool
    domain::D
    beta::Float64
    sigma::Float64
    gaussian_parameters::Union{Tuple{Float64, Float64}, Tuple{Nothing, Nothing}}
    a::Union{Float64, Nothing}
    b::Union{Float64, Nothing}
    num_energy_bits::Union{Int, Nothing}
    t0::Union{Float64, Nothing}
    w0::Union{Float64, Nothing}
    eta::Union{Float64, Nothing}
    num_trotter_steps_per_t0::Union{Int, Nothing}

    mixing_time::Float64
    delta::Float64

    function ThermalizeConfigGNS{D}(
        num_qubits::Int64, with_coherent::Bool, with_linear_combination::Bool, domain::D,
        beta::Float64, sigma::Float64,
        gaussian_parameters::Union{Tuple{Float64, Float64}, Tuple{Nothing, Nothing}},
        a::Union{Float64, Nothing}, b::Union{Float64, Nothing},
        num_energy_bits::Union{Int, Nothing}, t0::Union{Float64, Nothing},
        w0::Union{Float64, Nothing}, eta::Union{Float64, Nothing},
        num_trotter_steps_per_t0::Union{Int, Nothing},
        mixing_time::Float64, delta::Float64
    ) where {D}
        with_coherent && error("GNS configs must have with_coherent=false")
        new{D}(num_qubits, with_coherent, with_linear_combination, domain,
               beta, sigma, gaussian_parameters, a, b,
               num_energy_bits, t0, w0, eta, num_trotter_steps_per_t0,
               mixing_time, delta)
    end
end

# Keyword constructor for ThermalizeConfigGNS (replaces @kwdef)
function ThermalizeConfigGNS(;
    num_qubits::Int64,
    with_coherent::Bool = false,
    with_linear_combination::Bool,
    domain::D,
    beta::Float64,
    sigma::Float64,
    gaussian_parameters::Union{Tuple{Float64, Float64}, Tuple{Nothing, Nothing}} = (nothing, nothing),
    a::Union{Float64, Nothing} = nothing,
    b::Union{Float64, Nothing} = nothing,
    num_energy_bits::Union{Int, Nothing} = nothing,
    t0::Union{Float64, Nothing} = nothing,
    w0::Union{Float64, Nothing} = nothing,
    eta::Union{Float64, Nothing} = nothing,
    num_trotter_steps_per_t0::Union{Int, Nothing} = nothing,
    mixing_time::Float64,
    delta::Float64,
) where {D <: AbstractDomain}
    ThermalizeConfigGNS{D}(num_qubits, with_coherent, with_linear_combination, domain,
                           beta, sigma, gaussian_parameters, a, b,
                           num_energy_bits, t0, w0, eta, num_trotter_steps_per_t0,
                           mixing_time, delta)
end

"""
    JumpOp

    Represents an operator from which we can build the Lindbladian jump operators later.

    # Fields
    - `data`: The operator in the computational basis.
    - `in_eigenbasis`: The operator transformed into the Hamiltonian's eigenbasis (or Trotter basis).
    - `orthogonal`: Boolean flag indicating if this operator is self-orthogonal. If yes, the algorithm simplifies a bit.
"""
struct JumpOp{T <: AbstractMatrix{ComplexF64}}
    data::T
    in_eigenbasis::Matrix{ComplexF64}
    orthogonal::Bool
    hermitian::Bool
end

"""
        HotAlgorithmResults{D}

    Results from the step-by-step quantum algorithm emulation on thermalization.

    # Fields
    - `evolved_dm`: The final density matrix after evolution.
    - `distances_to_gibbs`: Trace distances to the target Gibbs state at each time step.
    - `time_steps`: Vector of time points where data was recorded.
    - `hamiltonian`: The [`HamHam`](@ref) data used.
    - `trotter`: The [`TrottTrott`](@ref) data used, in case of a TrotterDomain simulation.
    - `config`: The given configuration used.
"""
@kwdef struct HotAlgorithmResults{D}
    evolved_dm::Matrix{ComplexF64}
    distances_to_gibbs::Vector{Float64}
    time_steps::Vector{Float64}
    hamiltonian::HamHam
    trotter::Union{TrottTrott,Nothing} = nothing
    config::AbstractThermalizeConfig{D}
end

"""
        HotSpectralResults{D}

    Results from the spectral analysis of the Liouvillian.

    # Fields
    - `data`: The Liouvillian matrix.
    - `fixed_point`: The steady state found via spectral analysis.
    - `gap_mode`: The next eigenmode after the steady state.
    - `spectral_gap`: The first non-zero eigenvalue (gap).
    - `hamiltonian`: The [`HamHam`](@ref) data used.
    - `trotter`: The [`TrottTrott`](@ref) data used, in case of a TrotterDomain simulation.
    - `config`: The given configuration used.
"""
@kwdef struct HotSpectralResults{D}
    data::Matrix{ComplexF64}  #! Remove when space matters
    fixed_point::Matrix{ComplexF64}
    gap_mode::Matrix{ComplexF64}
    spectral_gap::ComplexF64
    hamiltonian::HamHam
    trotter::Union{TrottTrott,Nothing} = nothing
    config::AbstractLiouvConfig{D}
end

struct LSIFramework{T}
    dim::Int

    A::Matrix{T}            # Parameter matrix
    AdagA::Matrix{T}            # B = A'A
    Gamma2_AdagA::Matrix{T}  # Γ_2(A'A) = sig^1/4 A'A sig^1/4
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