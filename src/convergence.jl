# ============================================================================
# Convergence Tracking: batch-level monitoring for trajectory sampling
# ============================================================================
#
# ConvergenceData struct is defined in structs.jl (must be available before
# trajectories.jl which embeds it in TrajectoryResult).
#

# ---------------------------------------------------------------------------
# Observable builders
# ---------------------------------------------------------------------------

"""
    build_convergence_observables(hamiltonian::HamHam, num_qubits::Int)

Build nearest-neighbor Z_iZ_{i+1} correlation matrices and the energy observable <H>,
all in the Hamiltonian eigenbasis.

Returns `(observables::Vector{Matrix{ComplexF64}}, names::Vector{String})`.

CRITICAL: All observables are in the Hamiltonian eigenbasis because `run_trajectories`
returns `rho_mean` in the eigenbasis for EnergyDomain. The Gibbs state `hamiltonian.gibbs`
is also in eigenbasis. Using computational basis observables would give wrong `tr(rho * O)` values.
"""
function build_convergence_observables(hamiltonian::HamHam, num_qubits::Int)
    V = hamiltonian.eigvecs
    observables = Matrix{ComplexF64}[]
    names = String[]

    # Nearest-neighbor Z_iZ_{i+1} correlations (periodic)
    for i in 1:num_qubits
        ZZ_comp = Matrix{ComplexF64}(pad_term([Z, Z], num_qubits, i; periodic=true))
        # Transform to eigenbasis: O_eigen = V' * O_comp * V
        ZZ_eigen = V' * ZZ_comp * V
        push!(observables, ZZ_eigen)
        j = (i % num_qubits) + 1
        push!(names, "ZZ_$(i)$(j)")
    end

    # Energy observable: diagonal matrix of eigenvalues in eigenbasis
    H_eigen = Matrix{ComplexF64}(diagm(ComplexF64.(hamiltonian.eigvals)))
    push!(observables, H_eigen)
    push!(names, "H")

    return observables, names
end

"""
    build_convergence_observables_trotter(hamiltonian::HamHam, trotter::TrottTrott, num_qubits::Int)

Build nearest-neighbor Z_iZ_{i+1} correlation matrices and the energy observable <H>,
all in the Trotter eigenbasis (for TrotterDomain simulations).

Returns `(observables::Vector{Matrix{ComplexF64}}, names::Vector{String})`.
"""
function build_convergence_observables_trotter(hamiltonian::HamHam, trotter::TrottTrott, num_qubits::Int)
    V_T = trotter.eigvecs
    observables = Matrix{ComplexF64}[]
    names = String[]

    # Nearest-neighbor Z_iZ_{i+1} correlations (periodic)
    for i in 1:num_qubits
        ZZ_comp = Matrix{ComplexF64}(pad_term([Z, Z], num_qubits, i; periodic=true))
        # Transform to Trotter eigenbasis: O_trotter = V_T' * O_comp * V_T
        ZZ_trotter = V_T' * ZZ_comp * V_T
        push!(observables, ZZ_trotter)
        j = (i % num_qubits) + 1
        push!(names, "ZZ_$(i)$(j)")
    end

    # Energy observable in Trotter basis: V_T' * H_comp * V_T (NOT diagonal)
    H_trotter = V_T' * hamiltonian.data * V_T
    push!(observables, Matrix{ComplexF64}(H_trotter))
    push!(names, "H")

    return observables, names
end

"""
    build_total_magnetization(hamiltonian::HamHam, num_qubits::Int; trotter=nothing)

Build the per-site total magnetization observable M_z = sum(Z_i) / n,
transformed into the Hamiltonian eigenbasis (default) or the Trotter eigenbasis
(when `trotter` is supplied).

Per-site normalization (division by `num_qubits`) is applied so that the
amplitude is system-size-independent.  The decay rate (spectral gap) is
unaffected by normalization.

Returns `(observables::Vector{Matrix{ComplexF64}}, names::Vector{String})`
matching the existing observable builder convention.
"""
function build_total_magnetization(hamiltonian::HamHam, num_qubits::Int;
                                    trotter::Union{TrottTrott, Nothing}=nothing)
    dim = size(hamiltonian.data, 1)

    # Build M_z in computational basis: sum(Z_i) / n
    Mz_comp = zeros(ComplexF64, dim, dim)
    for i in 1:num_qubits
        Mz_comp .+= Matrix{ComplexF64}(pad_term([Z], num_qubits, i))
    end
    Mz_comp ./= num_qubits  # Per-site normalization (locked decision)

    # Select basis transformation matrix
    V = trotter !== nothing ? trotter.eigvecs : hamiltonian.eigvecs

    # Transform to eigenbasis: O_eigen = V' * O_comp * V
    Mz_eigen = Matrix{ComplexF64}(V' * Mz_comp * V)

    return [Mz_eigen], ["Mz"]
end

"""
    build_gap_estimation_observables(hamiltonian::HamHam, num_qubits::Int; trotter=nothing)

Build the observable bundle for spectral gap estimation, in the Hamiltonian
eigenbasis (default) or the Trotter eigenbasis (when `trotter` is supplied).

Returns `(observables::Vector{Matrix{ComplexF64}}, names::Vector{String})`
with 5 observables: `["H", "Mz", "XX_avg", "YY_avg", "ZZ_avg"]`.

- `H`: Energy observable.
- `Mz`: Per-site total magnetization (sum of Z_i / n).
- `XX_avg`, `YY_avg`, `ZZ_avg`: Per-bond averaged nearest-neighbor two-site
  correlations (sum over all periodic bonds, divided by n).
"""
function build_gap_estimation_observables(hamiltonian::HamHam, num_qubits::Int;
                                           trotter::Union{TrottTrott, Nothing}=nothing)
    # Get M_z from the dedicated builder (DRY)
    mz_obs, mz_names = build_total_magnetization(hamiltonian, num_qubits; trotter=trotter)

    # Build H in the correct basis
    if trotter !== nothing
        # H is NOT diagonal in Trotter basis — full basis transform required
        H = Matrix{ComplexF64}(trotter.eigvecs' * hamiltonian.data * trotter.eigvecs)
    else
        # H IS diagonal in its own eigenbasis
        H = Matrix{ComplexF64}(diagm(ComplexF64.(hamiltonian.eigvals)))
    end

    # Mutable collections for push! in the correlation loop below
    observables = Matrix{ComplexF64}[H]
    names = String["H"]
    append!(observables, mz_obs)
    append!(names, mz_names)

    # Averaged two-site correlations (nearest-neighbor, periodic)
    V = trotter !== nothing ? trotter.eigvecs : hamiltonian.eigvecs
    for (pauli_pair, pair_name) in [([X, X], "XX_avg"), ([Y, Y], "YY_avg"), ([Z, Z], "ZZ_avg")]
        dim = size(hamiltonian.data, 1)
        PP_sum = zeros(ComplexF64, dim, dim)
        for i in 1:num_qubits
            PP_sum .+= Matrix{ComplexF64}(pad_term(pauli_pair, num_qubits, i; periodic=true))
        end
        PP_sum ./= num_qubits  # Per-bond average
        PP_eigen = Matrix{ComplexF64}(V' * PP_sum * V)
        push!(observables, PP_eigen)
        push!(names, pair_name)
    end

    return observables, names
end

# ---------------------------------------------------------------------------
# Gibbs state helpers
# ---------------------------------------------------------------------------

"""
    _gibbs_in_trotter_basis(hamiltonian::HamHam, trotter::TrottTrott) -> Hermitian

Compute the Gibbs state in the Trotter eigenbasis.
Follows the same transform as furnace.jl run_thermalization (line 116-118).
"""
function _gibbs_in_trotter_basis(hamiltonian::HamHam, trotter::TrottTrott)
    return Hermitian(trotter.eigvecs' * hamiltonian.eigvecs * hamiltonian.gibbs *
                     hamiltonian.eigvecs' * trotter.eigvecs)
end

"""
    _compute_gibbs_observable_values(gibbs::Hermitian, observables::Vector{<:Matrix{<:Complex}})

Compute the Gibbs expectation value for each observable: tr(gibbs * O_i).
"""
function _compute_gibbs_observable_values(gibbs::Hermitian, observables::Vector{<:Matrix{<:Complex}})
    return [real(tr(Matrix(gibbs) * obs)) for obs in observables]
end

# ---------------------------------------------------------------------------
# Windowed relative change (adaptive convergence detection)
# ---------------------------------------------------------------------------

"""
    _windowed_relative_change(trace_dists::Vector{Float64}, window_size::Int) -> Float64

Compute the relative change between two consecutive non-overlapping windows
of trace distance values. Compares the mean of the most recent `window_size`
entries against the mean of the previous `window_size` entries.

Returns `Inf` if fewer than `2 * window_size` data points are available
(not enough data for two full windows).
"""
function _windowed_relative_change(trace_dists::Vector{Float64}, window_size::Int)
    n = length(trace_dists)
    n < 2 * window_size && return Inf
    mean_recent = sum(@view trace_dists[n - window_size + 1 : n]) / window_size
    mean_previous = sum(@view trace_dists[n - 2*window_size + 1 : n - window_size]) / window_size
    return abs(mean_recent - mean_previous) / max(abs(mean_previous), eps(Float64))
end

# ---------------------------------------------------------------------------
# Main convergence runner
# ---------------------------------------------------------------------------

"""
    run_trajectories_convergence(jumps, config, psi0, hamiltonian; kwargs...)

Run trajectory simulations in batches, measuring convergence metrics after each batch.

Returns a `TrajectoryResult` with the `convergence` field populated.

# Arguments
- `jumps`: Jump operators.
- `config`: Thermalization configuration.
- `psi0`: Initial state vector.
- `hamiltonian`: Hamiltonian data.

# Keyword Arguments
- `gibbs`: Reference Gibbs state (Hermitian, in the correct basis).
- `observables`: Observable matrices (in the same basis as rho).
- `observable_names`: Names for each observable.
- `batch_size`: Number of trajectories per batch (default: 1000).
- `n_batches`: Number of batches to run (default: 10).
- `seed`: Master RNG seed (default: random from system entropy).
- `trotter`: TrottTrott object (required for TrotterDomain).
- `total_time`: Evolution time per trajectory (default: config.mixing_time).
- `delta`: Time step (default: config.delta).
"""
function run_trajectories_convergence(
    jumps::Vector{JumpOp},
    config::AbstractThermalizeConfig,
    psi0::Vector{<:Complex},
    hamiltonian::HamHam;
    gibbs::Hermitian,
    observables::Vector{<:Matrix{<:Complex}},
    observable_names::Vector{String},
    batch_size::Int = 1000,
    n_batches::Int = 10,
    seed::Union{Int,Nothing} = nothing,
    trotter::Union{TrottTrott,Nothing} = nothing,
    total_time::Real = config.mixing_time,
    delta::Real = config.delta,
)
    # Build framework ONCE (not per batch)
    fw, actual_seed = _build_framework_and_seed(
        jumps, config, psi0, hamiltonian;
        trotter=trotter, delta=delta, seed=seed,
    )

    CT = eltype(psi0)
    dim = length(psi0)
    n_obs = length(observables)

    # Pre-allocate running sum accumulator
    rho_acc = zeros(CT, dim, dim)
    n_total = 0

    # Pre-compute Gibbs observable values
    obs_gibbs = _compute_gibbs_observable_values(gibbs, observables)

    # Pre-allocate convergence data storage
    trace_dists = Vector{Float64}(undef, n_batches)
    obs_values = Matrix{Float64}(undef, n_obs, n_batches)
    cummulative_n_traj = Vector{Int}(undef, n_batches)
    batch_sizes_vec = Vector{Int}(undef, n_batches)

    for batch_idx in 1:n_batches
        # Seed offset ensures non-overlapping trajectory seed sequences
        batch_seed = actual_seed + n_total

        # Run batch using pre-built framework (no rebuild per batch)
        rho_batch = _run_batch_no_obs!(fw, psi0, batch_size, batch_seed, total_time)

        # Accumulate running sum (not average)
        rho_acc .+= rho_batch .* batch_size
        n_total += batch_size

        # Compute running average
        rho_running = rho_acc ./ n_total
        hermitianize!(rho_running)

        # Measure trace distance to Gibbs
        trace_dists[batch_idx] = trace_distance_h(Hermitian(rho_running), gibbs)

        # Measure observable values: tr(rho * O)
        for i in 1:n_obs
            obs_values[i, batch_idx] = real(tr(rho_running * observables[i]))
        end

        cummulative_n_traj[batch_idx] = n_total
        batch_sizes_vec[batch_idx] = batch_size
    end

    # Compute final averaged density matrix
    rho_final = rho_acc ./ n_total
    hermitianize!(rho_final)

    conv_data = ConvergenceData(
        batch_sizes_vec,
        cummulative_n_traj,
        trace_dists,
        observable_names,
        obs_values,
        obs_gibbs,
    )

    return TrajectoryResult(rho_final, n_total, actual_seed, nothing, nothing, conv_data)
end

# ---------------------------------------------------------------------------
# Adaptive convergence runner
# ---------------------------------------------------------------------------

"""
    run_trajectories_adaptive(jumps, config, psi0, hamiltonian; kwargs...)

Run trajectory simulations in adaptive batches, stopping automatically when
trace distance convergence is detected or a hard trajectory cap is reached.

Returns a `TrajectoryResult` with the `convergence` field populated
(same as `run_trajectories_convergence`, but with adaptive diagnostic fields
in the `ConvergenceData`).

# Convergence criterion
Convergence is declared when the windowed relative change of trace distances
stays below `convergence_threshold` for `patience` consecutive batch checkpoints.
The windowed relative change compares the mean trace distance of the last
`window_size` batches against the mean of the previous `window_size` batches.
Only trace distance gates stopping; observables are tracked but do not affect
the stopping decision.

# Arguments
- `jumps`: Jump operators.
- `config`: Thermalization configuration.
- `psi0`: Initial state vector.
- `hamiltonian`: Hamiltonian data.

# Keyword Arguments
- `gibbs`: Reference Gibbs state (Hermitian, in the correct basis).
- `observables`: Observable matrices (in the same basis as rho).
- `observable_names`: Names for each observable.
- `batch_size`: Number of trajectories per batch (default: 200). Fixed throughout.
- `n_max`: Hard cap on total trajectories (default: 20_000). Ceiling division
  is used, so the actual maximum may slightly exceed `n_max` if it is not a
  multiple of `batch_size`.
- `convergence_threshold`: Relative change threshold for convergence (default: 0.01).
- `patience`: Number of consecutive stable checks required (default: 3).
- `min_batches`: Minimum batches before convergence can trigger (default: 5).
  Automatically clamped to `2 * window_size` if smaller.
- `window_size`: Window size for windowed average comparison (default: 3).
- `seed`: Master RNG seed (default: random from system entropy).
- `trotter`: TrottTrott object (required for TrotterDomain).
- `total_time`: Evolution time per trajectory (default: config.mixing_time).
- `delta`: Time step (default: config.delta).
"""
function run_trajectories_adaptive(
    jumps::Vector{JumpOp},
    config::AbstractThermalizeConfig,
    psi0::Vector{<:Complex},
    hamiltonian::HamHam;
    gibbs::Hermitian,
    observables::Vector{<:Matrix{<:Complex}},
    observable_names::Vector{String},
    batch_size::Int = 200,
    n_max::Int = 20_000,
    convergence_threshold::Float64 = 0.01,
    patience::Int = 3,
    min_batches::Int = 5,
    window_size::Int = 3,
    seed::Union{Int,Nothing} = nothing,
    trotter::Union{TrottTrott,Nothing} = nothing,
    total_time::Real = config.mixing_time,
    delta::Real = config.delta,
)
    # Build framework ONCE (not per batch)
    fw, actual_seed = _build_framework_and_seed(
        jumps, config, psi0, hamiltonian;
        trotter=trotter, delta=delta, seed=seed,
    )

    max_batches = cld(n_max, batch_size)

    # Ensure enough data for two full windows before convergence checking
    effective_min = max(min_batches, 2 * window_size)

    CT = eltype(psi0)
    dim = length(psi0)
    n_obs = length(observables)

    # Running sum accumulator
    rho_acc = zeros(CT, dim, dim)
    n_total = 0

    # Pre-compute Gibbs observable values
    obs_gibbs = _compute_gibbs_observable_values(gibbs, observables)

    # Dynamic storage (unknown final batch count)
    trace_dists = Float64[]
    obs_values_list = Vector{Float64}[]
    cum_n_traj = Int[]
    batch_sizes_vec = Int[]

    converged = false
    consecutive_stable = 0
    last_relative_change = NaN

    for batch_idx in 1:max_batches
        # Seed offset ensures non-overlapping trajectory seed sequences
        batch_seed = actual_seed + n_total

        # Run batch using pre-built framework (no rebuild per batch)
        rho_batch = _run_batch_no_obs!(fw, psi0, batch_size, batch_seed, total_time)

        # Accumulate running sum (not average)
        rho_acc .+= rho_batch .* batch_size
        n_total += batch_size

        # Compute running average
        rho_running = rho_acc ./ n_total
        hermitianize!(rho_running)

        # Measure trace distance to Gibbs
        push!(trace_dists, trace_distance_h(Hermitian(rho_running), gibbs))

        # Measure observable values: tr(rho * O)
        push!(obs_values_list, [real(tr(rho_running * observables[i])) for i in 1:n_obs])
        push!(cum_n_traj, n_total)
        push!(batch_sizes_vec, batch_size)

        # Convergence check (only after burn-in with enough data for windowed comparison)
        if batch_idx >= effective_min
            rel_change = _windowed_relative_change(trace_dists, window_size)
            last_relative_change = rel_change

            if rel_change < convergence_threshold
                consecutive_stable += 1
            else
                consecutive_stable = 0
            end

            if consecutive_stable >= patience
                converged = true
                break
            end
        end
    end

    # Compute final averaged density matrix
    rho_final = rho_acc ./ n_total
    hermitianize!(rho_final)

    # Build observable values matrix from collected vectors
    obs_values = reduce(hcat, obs_values_list)

    conv_data = ConvergenceData(
        batch_sizes_vec, cum_n_traj, trace_dists,
        observable_names, obs_values, obs_gibbs,
        converged, last_relative_change, consecutive_stable, length(trace_dists),
    )

    return TrajectoryResult(rho_final, n_total, actual_seed, nothing, nothing, conv_data)
end
