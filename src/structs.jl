# Domains
abstract type AbstractDomain end

struct BohrDomain <: AbstractDomain end
struct EnergyDomain <: AbstractDomain end
struct TimeDomain <: AbstractDomain end
struct TrotterDomain <: AbstractDomain end

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

# ---------------------------------------------------------------------------
# New typed Result structs (Phase 36)
# ---------------------------------------------------------------------------

abstract type AbstractResults end

"""
    LindbladResults{T<:AbstractFloat} <: AbstractResults

Results from dense Liouvillian spectral analysis (`run_lindblad`).
Stores leading eigenvalues, fixed point, gap mode, and spectral gap.
Does NOT store the full Liouvillian matrix (prohibitively large at scale).
"""
struct LindbladResults{T<:AbstractFloat} <: AbstractResults
    config::Config
    eigenvalues::Vector{Complex{T}}
    fixed_point::Matrix{Complex{T}}
    gap_mode::Matrix{Complex{T}}
    spectral_gap::Complex{T}
    metadata::Dict{Symbol, Any}
end

"""
    ThermalizeResults{T<:AbstractFloat} <: AbstractResults

Results from density-matrix Kraus evolution (`run_thermalize`).
Includes trace distances to the Gibbs state over time for convergence plotting.
"""
struct ThermalizeResults{T<:AbstractFloat} <: AbstractResults
    config::Config
    final_dm::Matrix{Complex{T}}
    trace_distances::Vector{T}
    time_steps::Vector{T}
    metadata::Dict{Symbol, Any}
end

"""
    KrylovSpectrumResults{T<:AbstractFloat} <: AbstractResults

Results from Krylov-based spectral gap estimation (`run_krylov_spectrum`).
Stores eigenvalues, gap, fixed point, convergence info, and optional channel data.
"""
struct KrylovSpectrumResults{T<:AbstractFloat} <: AbstractResults
    config::Config
    eigenvalues::Vector{Complex{T}}
    spectral_gap::T
    fixed_point::Matrix{Complex{T}}
    gap_mode::Matrix{Complex{T}}
    converged::Int
    matvec_count::Int
    num_restarts::Int
    normres::Vector{T}
    channel_eigenvalues::Union{Nothing, Vector{Complex{T}}}
    delta_used::Union{Nothing, T}
    metadata::Dict{Symbol, Any}
end

"""
    TrajectoryResults{T<:AbstractFloat} <: AbstractResults

Results from trajectory-based quantum simulation (`run_trajectory`).
Observable and convergence data are `Union{Nothing, ...}` -- only populated
when those features were used.
"""
struct TrajectoryResults{T<:AbstractFloat} <: AbstractResults
    config::Config
    rho_mean::Matrix{Complex{T}}
    n_trajectories::Int
    seed::Int
    times::Union{Nothing, Vector{Float64}}
    measurements_mean::Union{Nothing, Matrix{Float64}}
    convergence::Union{Nothing, ConvergenceData}
    metadata::Dict{Symbol, Any}
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

# ---------------------------------------------------------------------------
# Scratch sub-structs for Workspace (Phase 35)
# ---------------------------------------------------------------------------

"""
    LiouvillianScratch{T<:Complex}

Scratch buffers for dense Liouvillian construction (`construct_lindbladian`).
Replaces the old `LindbladianWorkspace` (Id is now computed inline at call sites).
"""
struct LiouvillianScratch{T<:Complex}
    jump_tmp::Matrix{T}
    jump_conj::Matrix{T}
    jump_dag_jump::Matrix{T}
    jump2_jump1::Matrix{T}
end

function LiouvillianScratch(::Type{CT}, dim::Int) where {CT<:Complex}
    Zm() = zeros(CT, dim, dim)
    return LiouvillianScratch{CT}(Zm(), Zm(), Zm(), Zm())
end

"""
    ThermalizeScratch{T<:Complex}

Scratch buffers for DM Kraus evolution (`run_thermalization`).
Replaces the old `KrausScratch` with physics-descriptive names and dead K0 removed.
"""
struct ThermalizeScratch{T<:Complex}
    jump_oft::Matrix{T}
    LdagL::Matrix{T}
    R::Matrix{T}
    rho_jump::Matrix{T}
    sandwich_tmp::Matrix{T}    # was tmp1
    rho_work::Matrix{T}        # was tmp2
    rho_next::Matrix{T}
end

function ThermalizeScratch(::Type{CT}, dim::Int) where {CT<:Complex}
    Zm() = zeros(CT, dim, dim)
    return ThermalizeScratch{CT}(Zm(), Zm(), Zm(), Zm(), Zm(), Zm(), Zm())
end

"""
    KrylovScratch{T<:Complex}

Scratch buffers for Krylov matvec and eigsolve hot paths.
Fields use physics-descriptive names (sandwich_tmp replaces tmp1, sandwich_out replaces LdagL).
"""
struct KrylovScratch{T<:Complex}
    jump_oft::Matrix{T}
    sandwich_tmp::Matrix{T}    # was tmp1 (BLAS gemm scratch)
    sandwich_out::Matrix{T}    # was LdagL (sandwich result)
    rho_out::Matrix{T}
    channel_rho_jump::Union{Nothing, Matrix{T}}  # Thermalize-channel only
end

function KrylovScratch(::Type{CT}, dim::Int; with_channel_rho_jump::Bool=false) where {CT<:Complex}
    Zm() = zeros(CT, dim, dim)
    crj = with_channel_rho_jump ? Zm() : nothing
    return KrylovScratch{CT}(Zm(), Zm(), Zm(), Zm(), crj)
end

"""
    TrajectoryScratch{T<:Complex}

Scratch buffers for trajectory simulation hot paths (`step_along_trajectory!`).
All fields are mutable per-trajectory working memory. Each thread needs its own
TrajectoryScratch to avoid shared mutable state.
"""
struct TrajectoryScratch{T<:Complex}
    jump_oft::Matrix{T}
    psi_tmp::Vector{T}
    Rpsi::Vector{T}
    rho_acc::Matrix{T}
end

function TrajectoryScratch(::Type{CT}, dim::Int) where {CT<:Complex}
    TrajectoryScratch{CT}(
        zeros(CT, dim, dim),  # jump_oft
        zeros(CT, dim),       # psi_tmp
        zeros(CT, dim),       # Rpsi
        zeros(CT, dim, dim),  # rho_acc
    )
end

"""
    Workspace{S, D, C, T}

Unified parametric workspace for all simulation paths (KrylovSpectrum, Lindbladian,
Thermalize, Trajectory).

Type parameters:
- `S <: AbstractSimulation`: simulation kind (KrylovSpectrum, Lindbladian, Thermalize, Trajectory)
- `D <: AbstractDomain`: domain (BohrDomain, EnergyDomain, TimeDomain, TrotterDomain)
- `C <: AbstractConstruction`: detailed-balance construction (KMS, GNS, DLL)
- `T <: AbstractFloat`: numeric precision

Dispatch signatures use partial parameterization:
`ws::Workspace{KrylovSpectrum}`, `ws::Workspace{Lindbladian}`, `ws::Workspace{Trajectory}`, etc.

Scratch sub-structs are accessed via type assertions at function entry points
(e.g. `sc = ws.scratch::KrylovScratch{T}`) for type-stable hot-path access.
"""
struct Workspace{S<:AbstractSimulation, D<:AbstractDomain, C<:AbstractConstruction, T<:AbstractFloat}
    # Physics data (Krylov/Thermalize)
    jump_eigenbases::Union{Nothing, Vector{Matrix{Complex{T}}}}
    jump_hermitian::Union{Nothing, Vector{Bool}}
    jumps::Union{Nothing, Vector{JumpOp}}
    B_total::Union{Nothing, Matrix{Complex{T}}}

    # Krylov effective Hamiltonian (Lindbladian mode)
    G_left::Union{Nothing, Matrix{Complex{T}}}
    G_right::Union{Nothing, Matrix{Complex{T}}}
    G_left_adj::Union{Nothing, Matrix{Complex{T}}}
    G_right_adj::Union{Nothing, Matrix{Complex{T}}}

    # CPTP channel (Krylov Thermalize mode, and Thermalize DM; also per-operator alpha/delta for Trajectory)
    K0::Union{Nothing, Matrix{Complex{T}}}
    U_residual::Union{Nothing, Matrix{Complex{T}}}
    U_coherent::Union{Nothing, Matrix{Complex{T}}}
    alpha::Union{Nothing, Float64}
    delta::Union{Nothing, Float64}

    # Domain-specific precomputed data (absorbed from NamedTuple)
    transition::Union{Nothing, Function}
    gamma_norm_factor::Union{Nothing, Float64}
    energy_labels::Union{Nothing, Vector{Float64}}
    oft_domain_prefactor::Union{Nothing, Float64}
    oft_nufft_prefactors::Any  # NUFFTPrefactors or Nothing
    bohr_alpha::Any            # BohrDomain alpha function (renamed to avoid clash with CPTP alpha)
    bohr_keys::Any             # BohrDomain keys
    bohr_is::Any               # BohrDomain row indices
    bohr_js::Any               # BohrDomain column indices
    b_minus::Any               # Time/TrotterDomain coherent
    b_plus::Any                # Time/TrotterDomain coherent

    # Thermalize DM-specific (coherent unitaries)
    coherent_unitaries::Union{Nothing, Vector{Matrix{Complex{T}}}}

    # Trajectory-specific fields (per-operator Lie-Trotter splitting)
    ham_or_trott::Any          # HamHam or TrottTrott (needed for EnergyDomain oft!())
    n_jumps::Union{Nothing, Int}
    scaled_prefactor::Union{Nothing, Float64}
    sigma::Union{Nothing, Float64}
    Rs::Union{Nothing, Vector{Matrix{Complex{T}}}}         # per-jump R^a
    K0s::Union{Nothing, Vector{Matrix{Complex{T}}}}        # per-jump K0^a
    U_residuals::Union{Nothing, Vector{Matrix{Complex{T}}}}  # per-jump U_residual^a
    U_Bs::Union{Nothing, Vector{Union{Nothing, Matrix{Complex{T}}}}}  # per-jump coherent unitary

    # Identity matrix (Lindbladian construction path)
    Id::Union{Nothing, Matrix{Complex{T}}}

    # Scratch buffers (nested, simulation-path-specific)
    scratch
end