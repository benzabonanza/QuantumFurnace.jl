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
- `a` and `s`: Parameters for the linear combination type.

## Grid parameters

The simulator uses up to **three independent register triples** `(num_energy_bits_X, t0_X, w0_X)`,
one per term that lives on its own QPE register on the quantum side. Each triple obeys
its own Fourier relation `w0_X * t0_X = 2π / 2^{num_energy_bits_X}` (with `t0_X` not
required for EnergyDomain since the dissipator uses analytical `A(ω)` there):

- `num_energy_bits_D`, `t0_D`, `w0_D` — **dissipative** OFT register (used by the
  CKG/GNS/DLL dissipator: OFT integral `Σ_t̄ b̄(t̄) e^{-iωt̄} A(t̄)`).
- `num_energy_bits_b_minus`, `t0_b_minus`, `w0_b_minus` — **outer** coherent
  integration register (`b_-(t)` loop in `B`); KMS-only.
- `num_energy_bits_b_plus`, `t0_b_plus`, `w0_b_plus` — **inner** coherent
  integration register (`b_+(τ)` loop in `B`); KMS-only.

For backward compatibility, the **legacy** kwargs `num_energy_bits`, `t0`, `w0` still
work and auto-promote to the three triples (`X_D = X_b_minus = X_b_plus = legacy_X`).
Mixing legacy and new on the same field is rejected at validation.

- `eta`: Accuracy coefficient for Metropolis linear combination in time domain.
- `num_trotter_steps_per_t0`: Trotter steps per unit time t0_D.

## Thermalize-specific
- `mixing_time`: Total duration of time evolution (only for `Thermalize` simulations).
- `delta`: Time step size for weak-measurement emulation (only for `Thermalize` simulations).

## GQSP-specific (Thermalize/Trajectory coherent step)
- `with_gqsp`: If `true`, use the GQSP polynomial approximation of `exp(-iδ B_a)` for the
  coherent step instead of the exact matrix exponential. Requires
  `with_coherent(construction)=true` and `domain isa Union{TimeDomain, TrotterDomain}`.
- `gqsp_degree`: Truncation degree `d ≥ 1` of the Jacobi-Anger polynomial. Default `1`
  is faithful to `O((δα)²)` and matches the splitting error.

## Generator splitting (Thermalize/Trajectory dissipative step)
- `jump_selection`: `:sweep` (default, thesis-preferred) deterministically cycles through
  the jump set per outer δ-step `Φ_𝓐 = e^{δ𝓛_S} ∘ ⋯ ∘ e^{δ𝓛_1} ≈ e^{δ𝓛}` with bare-δ
  rates per substep; `:random` picks one jump uniformly per outer step with rates
  rescaled by `S = |𝓐|` so that `E[step] ≈ e^{δ𝓛}`. The deterministic sweep has a
  strictly smaller leading-order spectral-gap perturbation (Methods §generator-splitting).

## OFT filter (DLL-1)
- `filter`: Optional `AbstractFilter` for the Operator Fourier Transform. `nothing`
  (default) selects the CKG Gaussian path with width `sigma` — byte-identical to the
  pre-filter codebase. A `DLLGaussianFilter(beta)` selects the Ding–Li–Lin Gevrey
  filter (Time/TrotterDomain only, currently the NUFFT prefactor path).

## Currently possible linear combinations:
(a, s) =
- (0, 0)   - plain Metropolis (kinky, eta-regularized in time domain)
- (0, >0)  - smooth Metropolis (eta-regularized, kink-smoothed by s; thesis-main case)
- (>0, 0)  - a-regularized smooth Metro (alternative regularization, no kink-smoothing)
- (>0, >0) - a-regularized Glauberish (smooth in both senses)

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
    s::Union{T, Nothing} = nothing

    # Grid parameters — per-register triples (qf-9z0). Each X-register obeys its
    # own Fourier relation `w0_X * t0_X = 2π / 2^{num_energy_bits_X}`. EnergyDomain
    # never needs `t0_X`. DLL TimeDomain never needs `w0_D`. Coherent registers
    # (`b_minus`, `b_plus`) are KMS-only.
    num_energy_bits_D::Union{Int, Nothing} = nothing
    t0_D::Union{T, Nothing} = nothing
    w0_D::Union{T, Nothing} = nothing
    num_energy_bits_b_minus::Union{Int, Nothing} = nothing
    t0_b_minus::Union{T, Nothing} = nothing
    w0_b_minus::Union{T, Nothing} = nothing
    num_energy_bits_b_plus::Union{Int, Nothing} = nothing
    t0_b_plus::Union{T, Nothing} = nothing
    w0_b_plus::Union{T, Nothing} = nothing
    # Legacy single-register kwargs — auto-promote to all three triples in
    # `validate_config!`. Test fixtures and scripts may still set them directly.
    num_energy_bits::Union{Int, Nothing} = nothing
    t0::Union{T, Nothing} = nothing
    w0::Union{T, Nothing} = nothing
    eta::Union{T, Nothing} = nothing
    num_trotter_steps_per_t0::Union{Int, Nothing} = nothing

    # Thermalize-specific
    mixing_time::Union{T, Nothing} = nothing
    delta::Union{T, Nothing} = nothing

    # GQSP-specific (Thermalize/Trajectory coherent step)
    with_gqsp::Bool = false
    gqsp_degree::Int = 1

    # Dissipative jump-selection rule (Thermalize/Trajectory): :sweep | :random.
    # :sweep is the thesis-preferred deterministic Lie-Trotter sweep over {A^a};
    # :random keeps the legacy uniform-random sampling with 1/p_a rate rescaling.
    jump_selection::Symbol = :sweep

    # OFT filter (DLL-1): nothing -> CKG Gaussian via config.sigma
    filter::Union{Nothing, AbstractFilter} = nothing
end

# ---------------------------------------------------------------------------
# Per-register accessors (qf-9z0).
#
# Each helper resolves to the explicit per-term field if set, else falls back
# to the legacy single-register field. Use these throughout `src/` instead of
# `cfg.t0` / `cfg.w0` / `cfg.num_energy_bits` so that downstream code can
# transparently consume either the legacy or the per-term API.
# ---------------------------------------------------------------------------

"""
    register_t0_D(cfg)        register_w0_D(cfg)        register_r_D(cfg)

Resolve the **dissipative** time/energy/bit triple `(t0_D, w0_D, r_D)`. If the
explicit per-term field is `nothing`, fall back to the legacy
`cfg.t0` / `cfg.w0` / `cfg.num_energy_bits`.
"""
@inline register_t0_D(cfg::Config) = cfg.t0_D !== nothing ? cfg.t0_D : cfg.t0
@inline register_w0_D(cfg::Config) = cfg.w0_D !== nothing ? cfg.w0_D : cfg.w0
@inline register_r_D(cfg::Config)  = cfg.num_energy_bits_D !== nothing ? cfg.num_energy_bits_D : cfg.num_energy_bits

"""
    register_t0_b_minus(cfg)  register_w0_b_minus(cfg)  register_r_b_minus(cfg)

Resolve the **outer coherent** triple `(t0_b_minus, w0_b_minus, r_b_minus)` —
the spacing of the `b_-(t)` outer Riemann sum in `B`. KMS-only. Falls back to
the legacy field when the explicit per-term field is `nothing`.
"""
@inline register_t0_b_minus(cfg::Config) = cfg.t0_b_minus !== nothing ? cfg.t0_b_minus : cfg.t0
@inline register_w0_b_minus(cfg::Config) = cfg.w0_b_minus !== nothing ? cfg.w0_b_minus : cfg.w0
@inline register_r_b_minus(cfg::Config)  = cfg.num_energy_bits_b_minus !== nothing ? cfg.num_energy_bits_b_minus : cfg.num_energy_bits

"""
    register_t0_b_plus(cfg)   register_w0_b_plus(cfg)   register_r_b_plus(cfg)

Resolve the **inner coherent** triple `(t0_b_plus, w0_b_plus, r_b_plus)` — the
spacing of the `b_+(τ)` inner Riemann sum in `B`. KMS-only. Falls back to the
legacy field when the explicit per-term field is `nothing`.
"""
@inline register_t0_b_plus(cfg::Config) = cfg.t0_b_plus !== nothing ? cfg.t0_b_plus : cfg.t0
@inline register_w0_b_plus(cfg::Config) = cfg.w0_b_plus !== nothing ? cfg.w0_b_plus : cfg.w0
@inline register_r_b_plus(cfg::Config)  = cfg.num_energy_bits_b_plus !== nothing ? cfg.num_energy_bits_b_plus : cfg.num_energy_bits

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

# Note: `OFTCaches` (cache struct used by the deprecated `time_oft!` /
# `trotter_oft!` direct-summation OFT routines) has been retired alongside
# those functions. The original definition is preserved for reference at
# `src/staging/ofts.jl`. Mainline code uses `NUFFTCaches` (above).

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

Scratch buffers for DM Kraus evolution (`run_thermalize`).
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

`task_scratches` is a pre-allocated pool of per-thread scratches used by the
qf-in3 ω-loop threading dispatch in `apply_lindbladian!`; `work_list` is a
pre-built `(jump_idx, label_idx)` schedule populated by the Workspace
constructor. Together they keep the per-matvec allocation budget down to the
intrinsic `Threads.@spawn` Task overhead (~1 kB) even when threading
dispatches.
"""
struct KrylovScratch{T<:Complex}
    jump_oft::Matrix{T}
    sandwich_tmp::Matrix{T}    # was tmp1 (BLAS gemm scratch)
    sandwich_out::Matrix{T}    # was LdagL (sandwich result)
    rho_out::Matrix{T}
    channel_rho_jump::Union{Nothing, Matrix{T}}  # Thermalize-channel only
    task_scratches::Vector{KrylovScratch{T}}     # per-thread pool (qf-in3.4)
    work_list::Vector{Tuple{Int, Int}}           # pre-built (k, li) schedule (qf-in3.4)
end

function KrylovScratch(::Type{CT}, dim::Int;
                       with_channel_rho_jump::Bool=false,
                       num_threads::Int=Threads.nthreads()) where {CT<:Complex}
    Zm() = zeros(CT, dim, dim)
    crj = with_channel_rho_jump ? Zm() : nothing

    # Per-thread scratch pool — only needed when threading. Each task scratch
    # itself carries its own empty `task_scratches` and `work_list` vectors
    # (terminates recursion / reduces footprint, and avoids any future
    # aliasing surprise if the chunk functions ever start reading those
    # fields).
    task_pool = if num_threads > 1
        [KrylovScratch{CT}(Zm(), Zm(), Zm(), Zm(),
                           with_channel_rho_jump ? Zm() : nothing,
                           KrylovScratch{CT}[],
                           Tuple{Int, Int}[]) for _ in 1:num_threads]
    else
        KrylovScratch{CT}[]
    end

    return KrylovScratch{CT}(Zm(), Zm(), Zm(), Zm(), crj, task_pool, Tuple{Int, Int}[])
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

    # DLL per-jump Bohr-domain Lindblad operators (Ding–Li–Lin 2024 Eq. 3.4 first form)
    dll_lindblads::Union{Nothing, Vector{Matrix{Complex{T}}}}

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
    oft_prefactors_energy::Any # EnergyDomainPrefactors or Nothing (qf-e60.1)
    bohr_alpha::Any            # BohrDomain alpha function (renamed to avoid clash with CPTP alpha)
    bohr_keys::Any             # BohrDomain keys
    bohr_is::Any               # BohrDomain row indices
    bohr_js::Any               # BohrDomain column indices
    b_minus::Any               # Time/TrotterDomain coherent
    b_plus::Any                # Time/TrotterDomain coherent

    # Thermalize DM-specific (coherent unitaries)
    coherent_unitaries::Union{Nothing, Vector{Matrix{Complex{T}}}}

    # Trajectory-specific fields (per-operator Lie-Trotter splitting)
    ham_or_trott::Any          # HamHam or TrottTrott (carried alongside cached prefactors)
    n_jumps::Union{Nothing, Int}
    scaled_prefactor::Union{Nothing, Float64}
    sigma::Union{Nothing, Float64}
    Rs::Union{Nothing, Vector{Matrix{Complex{T}}}}         # per-jump R^a
    K0s::Union{Nothing, Vector{Matrix{Complex{T}}}}        # per-jump K0^a
    U_residuals::Union{Nothing, Vector{Matrix{Complex{T}}}}  # per-jump U_residual^a
    U_Bs::Union{Nothing, Vector{Union{Nothing, Matrix{Complex{T}}}}}  # per-jump coherent unitary
    jump_selection::Union{Nothing, Symbol}  # :sweep | :random | nothing (non-Trajectory)

    # Identity matrix (Lindbladian construction path)
    Id::Union{Nothing, Matrix{Complex{T}}}

    # Scratch buffers (nested, simulation-path-specific)
    scratch
end