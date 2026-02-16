# ============================================================================
# Convergence Tracking: batch-level monitoring for trajectory sampling
# ============================================================================

"""
    ConvergenceData

Stores convergence metrics at batch checkpoints during trajectory sampling.
Scalars only (no density matrix snapshots) to keep memory O(n_batches).

# Fields
- `batch_sizes`: Number of trajectories in each batch.
- `cumulative_n_traj`: Running total of trajectories after each batch.
- `trace_distances`: Trace distance to Gibbs state at each checkpoint.
- `observable_names`: Names of the tracked observables (e.g. "ZZ_12", "H").
- `observable_values`: Observable expectation values, n_obs x n_checkpoints.
- `observable_gibbs_values`: Reference Gibbs expectation values for each observable.
"""
struct ConvergenceData
    batch_sizes::Vector{Int}
    cumulative_n_traj::Vector{Int}
    trace_distances::Vector{Float64}
    observable_names::Vector{String}
    observable_values::Matrix{Float64}      # n_obs x n_checkpoints
    observable_gibbs_values::Vector{Float64} # <O_i>_gibbs reference values
end

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
# Main convergence runner
# ---------------------------------------------------------------------------

"""
    run_trajectories_convergence(jumps, config, psi0, hamiltonian; kwargs...)

Run trajectory simulations in batches, measuring convergence metrics after each batch.

Returns `(TrajectoryResult, ConvergenceData)`.

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
    actual_seed = seed === nothing ? Int(rand(Random.RandomDevice(), UInt64) >> 1) : seed

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
    cum_n_traj = Vector{Int}(undef, n_batches)
    batch_sizes_vec = Vector{Int}(undef, n_batches)

    for batch_idx in 1:n_batches
        # Seed offset ensures non-overlapping trajectory seed sequences
        batch_seed = actual_seed + n_total

        # Run batch using existing run_trajectories (observables=nothing to avoid
        # time-resolved measurement tracking, which is different from batch-level convergence)
        result = run_trajectories(
            jumps, config, psi0, hamiltonian;
            trotter = trotter,
            total_time = total_time,
            delta = delta,
            ntraj = batch_size,
            seed = batch_seed,
        )

        # Accumulate running sum (not average)
        rho_acc .+= result.rho_mean .* batch_size
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

        cum_n_traj[batch_idx] = n_total
        batch_sizes_vec[batch_idx] = batch_size
    end

    # Compute final averaged density matrix
    rho_final = rho_acc ./ n_total
    hermitianize!(rho_final)

    conv_data = ConvergenceData(
        batch_sizes_vec,
        cum_n_traj,
        trace_dists,
        observable_names,
        obs_values,
        obs_gibbs,
    )

    traj_result = TrajectoryResult(rho_final, n_total, actual_seed, nothing, nothing)

    return traj_result, conv_data
end
