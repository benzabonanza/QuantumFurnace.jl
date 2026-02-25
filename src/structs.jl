# Domains
abstract type AbstractDomain end

struct BohrDomain <: AbstractDomain end
struct EnergyDomain <: AbstractDomain end
struct TimeDomain <: AbstractDomain end
struct TrotterDomain <: AbstractDomain end

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
struct LindbladianWorkspace{T<:AbstractFloat}
    Id::Matrix{Complex{T}}
    jump_tmp::Matrix{Complex{T}}
    jump_conj::Matrix{Complex{T}}
    jump_dag_jump::Matrix{Complex{T}}
    jump2_jump1::Matrix{Complex{T}}

    function LindbladianWorkspace{T}(dim::Int) where {T<:AbstractFloat}
        CT = Complex{T}
        Id = Matrix{CT}(I, dim, dim)
        jump_tmp = zeros(CT, dim, dim)
        jump_conj = zeros(CT, dim, dim)
        jump_dag_jump = zeros(CT, dim, dim)
        jump2_jump1 = zeros(CT, dim, dim)
        new{T}(Id, jump_tmp, jump_conj, jump_dag_jump, jump2_jump1)
    end
end
# Simulation types
abstract type AbstractSimulation end
struct Lindbladian    <: AbstractSimulation end
struct Thermalize     <: AbstractSimulation end
struct KrylovSpectrum <: AbstractSimulation end
struct Trajectory     <: AbstractSimulation end

# Construction types (detailed balance)
abstract type AbstractConstruction end
struct KMS <: AbstractConstruction end
struct GNS <: AbstractConstruction end
struct DLL <: AbstractConstruction end

# Trait: coherent term presence (derived from construction type)
with_coherent(::KMS) = true
with_coherent(::GNS) = false
with_coherent(::DLL) = true  # placeholder for Ding et al.

"""
    Config{S, D, C, T}

A unified configuration object holding all parameters for quantum Gibbs sampler simulations.

Type parameters encode the three dispatch axes:
- `S <: AbstractSimulation`: simulation kind (`Lindbladian`, `Thermalize`, `KrylovSpectrum`, `Trajectory`)
- `D <: AbstractDomain`: domain level (`BohrDomain`, `EnergyDomain`, `TimeDomain`, `TrotterDomain`)
- `C <: AbstractConstruction`: detailed-balance construction (`KMS`, `GNS`, `DLL`)
- `T <: AbstractFloat`: numeric precision

Whether the coherent correction term is included is derived from the construction type
via the trait function `with_coherent(construction)`, not stored as a field.

# Fields
## Type-encoding singletons
- `sim`: The simulation singleton (e.g. `Lindbladian()`).
- `domain`: The domain singleton (e.g. `EnergyDomain()`).
- `construction`: The construction singleton (e.g. `KMS()`).

## System parameters
- `num_qubits`: The number of system qubits.
- `with_linear_combination`: Whether to apply a convex combination of Lindbladians for faster mixing.

## Physics parameters
- `beta`: Inverse temperature.
- `sigma`: Gaussian width parameter.
- `gaussian_parameters`: Optional `(omega_gamma, sigma_gamma)` tuple for secondary Gaussian.
- `a` and `b`: Parameters for the linear combination type.

## Grid parameters
- `num_energy_bits`: Coarseness of energy/time grid.
- `t0` and `w0`: Time and energy units for Riemann-summed integrals (related by Fourier: w0*t0 = 2pi/N).
- `eta`: Accuracy coefficient for Metropolis linear combination in time domain.
- `num_trotter_steps_per_t0`: Trotter steps per unit time t0.

## Thermalize-specific
- `mixing_time`: Total duration of time evolution (only for `Thermalize` simulations).
- `delta`: Time step size for weak-measurement emulation (only for `Thermalize` simulations).

## Currently possible linear combinations:
(a, b) =
- (0, 0) - no linear combination, simple Gaussian
- (>0, 0) - linear combination that results in Metropolis-like transition
- (>0, >0) - linear combination that results in Glauber transition (smoother)

## Available domains:
- **`BohrDomain()`**: Highest level -- Lindbladian in Bohr frequency decomposition.
- **`EnergyDomain()`**: Operators approximated by energy integrals.
- **`TimeDomain()`**: Energy approximations as Fourier transforms of temporal equivalents.
- **`TrotterDomain()`**: Lowest level -- all time evolutions replaced by Trotter series.
"""
@kwdef struct Config{S <: AbstractSimulation, D <: AbstractDomain, C <: AbstractConstruction, T <: AbstractFloat}
    # Type-encoding singletons
    sim::S
    domain::D
    construction::C

    # System parameters
    num_qubits::Int
    with_linear_combination::Bool

    # Physics parameters
    beta::T
    sigma::T
    gaussian_parameters::Union{Tuple{T, T}, Tuple{Nothing, Nothing}} = (nothing, nothing)
    a::Union{T, Nothing} = nothing
    b::Union{T, Nothing} = nothing

    # Grid parameters
    num_energy_bits::Union{Int, Nothing} = nothing
    t0::Union{T, Nothing} = nothing
    w0::Union{T, Nothing} = nothing
    eta::Union{T, Nothing} = nothing
    num_trotter_steps_per_t0::Union{Int, Nothing} = nothing

    # Thermalize-specific
    mixing_time::Union{T, Nothing} = nothing
    delta::Union{T, Nothing} = nothing
end

# Outer constructor: infer S, D, C from singletons and T from beta.
# Required because @kwdef with 4 type parameters needs help with type inference.
function Config(;
    sim::S, domain::D, construction::C,
    beta::T,
    kwargs...
) where {S <: AbstractSimulation, D <: AbstractDomain, C <: AbstractConstruction, T <: AbstractFloat}
    Config{S, D, C, T}(; sim=sim, domain=domain, construction=construction, beta=beta, kwargs...)
end

"""
    JumpOp

    Represents an operator from which we can build the Lindbladian jump operators later.

    # Fields
    - `data`: The operator in the computational basis.
    - `in_eigenbasis`: The operator transformed into the Hamiltonian's eigenbasis (or Trotter basis).
    - `orthogonal`: Boolean flag indicating if this operator is self-orthogonal. If yes, the algorithm simplifies a bit.
"""
struct JumpOp{T <: AbstractMatrix{<:Complex}}
    data::T
    in_eigenbasis::Matrix{<:Complex}
    orthogonal::Bool
    hermitian::Bool
end

"""
        DMSimulationResult{T}

    Results from the step-by-step quantum algorithm emulation on thermalization.

    Slimmed version of the former HotAlgorithmResults: carries only the simulation
    output, not the hamiltonian/trotter/config (callers already have those at the call site).

    # Fields
    - `final_dm`: The final density matrix after evolution.
    - `trace_distances`: Trace distances to the target Gibbs state at each time step.
    - `time_steps`: Vector of time points where data was recorded.
"""
@kwdef struct DMSimulationResult{T<:AbstractFloat}
    final_dm::Matrix{Complex{T}}
    trace_distances::Vector{T}
    time_steps::Vector{T}
end

"""
        LindbladianResult{T}

    Results from the spectral analysis of the Liouvillian.

    Slimmed version of the former HotSpectralResults: carries only the spectral
    output, not the hamiltonian/trotter/config (callers already have those at the call site).

    # Fields
    - `liouvillian`: The Liouvillian matrix.
    - `fixed_point`: The steady state found via spectral analysis.
    - `gap_mode`: The next eigenmode after the steady state.
    - `spectral_gap`: The first non-zero eigenvalue (gap).
"""
@kwdef struct LindbladianResult{T<:AbstractFloat}
    liouvillian::Matrix{Complex{T}}
    fixed_point::Matrix{Complex{T}}
    gap_mode::Matrix{Complex{T}}
    spectral_gap::Complex{T}
end

"""
    ConvergenceData

Stores convergence metrics at batch checkpoints during trajectory sampling.
Scalars only (no density matrix snapshots) to keep memory O(n_batches).

# Fields (Phase 16 -- core convergence tracking)
- `batch_sizes`: Number of trajectories in each batch.
- `cumulative_n_traj`: Running total of trajectories after each batch.
- `trace_distances`: Trace distance to Gibbs state at each checkpoint.
- `observable_names`: Names of the tracked observables (e.g. "ZZ_12", "H").
- `observable_values`: Observable expectation values, n_obs x n_checkpoints.
- `observable_gibbs_values`: Reference Gibbs expectation values for each observable.

# Fields (Phase 17 -- adaptive diagnostics)
- `converged`: Did adaptive stopping trigger? (false for fixed-count runs)
- `final_relative_change`: Windowed relative change at termination (NaN for fixed-count runs).
- `consecutive_stable_batches`: How many consecutive stable checks achieved at termination.
- `total_batches`: Number of batches actually run.
"""
struct ConvergenceData
    # Phase 16: core convergence tracking
    batch_sizes::Vector{Int}
    cumulative_n_traj::Vector{Int}
    trace_distances::Vector{Float64}
    observable_names::Vector{String}
    observable_values::Matrix{Float64}      # n_obs x n_checkpoints
    observable_gibbs_values::Vector{Float64} # <O_i>_gibbs reference values
    # Phase 17: adaptive diagnostics
    converged::Bool
    final_relative_change::Float64
    consecutive_stable_batches::Int
    total_batches::Int
end

# Backward-compatible 6-argument outer constructor (Phase 16 callers pass 6 args).
# Uses broad types to accept BSON-deserialized data (e.g. Vector{Any} for strings).
function ConvergenceData(
    batch_sizes, cumulative_n_traj, trace_distances,
    observable_names, observable_values, observable_gibbs_values,
)
    ConvergenceData(
        batch_sizes, cumulative_n_traj, trace_distances,
        observable_names, observable_values, observable_gibbs_values,
        false, NaN, 0, length(batch_sizes),
    )
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

# Became obsolete with NUFFTCaches. But used for debugging.
struct OFTCaches{T<:AbstractFloat}
    prefactors::Vector{Complex{T}}
    U::Diagonal{Complex{T}, Vector{Complex{T}}}
    temp_op::Matrix{Complex{T}}

    function OFTCaches{T}(dim::Int) where {T<:AbstractFloat}
        CT = Complex{T}
        prefactors = zeros(CT, 0) # Will be resized later
        U = Diagonal(zeros(CT, dim))
        temp_op = zeros(CT, dim, dim)
        new{T}(prefactors, U, temp_op)
    end
end