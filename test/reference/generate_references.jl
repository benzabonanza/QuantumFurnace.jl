"""
One-time script to generate frozen BSON reference data for regression tests (TINF-02).

Run via: julia --project test/reference/generate_references.jl

Generates 4 reference files:
  - energy_dm_reference.bson       (EnergyDomain, DM evolution)
  - energy_traj_reference.bson     (EnergyDomain, trajectory average)
  - trotter_coherent_dm_reference.bson   (TrotterDomain+coherent, DM evolution)
  - trotter_coherent_traj_reference.bson (TrotterDomain+coherent, trajectory average)

Each BSON file stores a plain Matrix{ComplexF64} density matrix plus metadata.
"""

using QuantumFurnace, LinearAlgebra, Random, BSON

# Load shared test fixtures (SMALL system: 3-qubit)
include(joinpath(@__DIR__, "..", "test_helpers.jl"))

# ---------------------------------------------------------------------------
# Constants for reference generation
# ---------------------------------------------------------------------------
const REF_DIR = @__DIR__
const REF_DELTA = 0.1
const REF_SEED = 12345
const REF_NTRAJ = 1000

# ---------------------------------------------------------------------------
# DM reference generator
# ---------------------------------------------------------------------------
"""
    generate_dm_reference(domain; with_coherent, filename)

Compute a single-step density matrix evolution via exp(delta*L) and save as BSON.
"""
function generate_dm_reference(domain; with_coherent::Bool, filename::String)
    liouv_config = make_small_liouv_config(domain; with_coherent=with_coherent)
    trotter_kw = domain isa TrotterDomain ? (; trotter=SMALL_TROTTER) : (;)
    L = construct_lindbladian(SMALL_JUMPS, liouv_config, SMALL_HAM; trotter_kw...)

    psi0 = fill(ComplexF64(1.0), SMALL_DIM) / sqrt(SMALL_DIM)
    rho0 = psi0 * psi0'
    rho_dm = reshape(exp(REF_DELTA * L) * vec(rho0), SMALL_DIM, SMALL_DIM)
    rho_dm = (rho_dm + rho_dm') / 2  # Hermitianize

    BSON.bson(joinpath(REF_DIR, filename), Dict(
        :rho => Matrix(rho_dm),
        :delta => REF_DELTA,
        :domain => string(typeof(domain)),
        :with_coherent => with_coherent,
        :num_qubits => 3,
    ))
    println("Saved: $filename")
end

# ---------------------------------------------------------------------------
# Trajectory reference generator
# ---------------------------------------------------------------------------
"""
    generate_traj_reference(domain; with_coherent, filename)

Run REF_NTRAJ single-step trajectories with fixed seed and save averaged rho as BSON.
"""
function generate_traj_reference(domain; with_coherent::Bool, filename::String)
    therm_config = make_small_thermalize_config(domain;
        with_coherent=with_coherent, delta=REF_DELTA, mixing_time=Float64(REF_DELTA))
    ham_or_trott = domain isa TrotterDomain ? SMALL_TROTTER : SMALL_HAM
    precomputed = precompute_data(domain, therm_config, ham_or_trott)
    scratch = KrausScratch(ComplexF64, SMALL_DIM)
    fw = build_trajectoryframework(SMALL_JUMPS, ham_or_trott, therm_config,
        precomputed, scratch, REF_DELTA)

    psi0 = fill(ComplexF64(1.0), SMALL_DIM) / sqrt(SMALL_DIM)

    Random.seed!(REF_SEED)
    rho_traj = zeros(ComplexF64, SMALL_DIM, SMALL_DIM)
    for _ in 1:REF_NTRAJ
        psi = copy(psi0)
        step_along_trajectory!(psi, fw)
        rho_traj .+= psi * psi'
    end
    rho_traj ./= REF_NTRAJ
    rho_traj = (rho_traj + rho_traj') / 2  # Hermitianize

    BSON.bson(joinpath(REF_DIR, filename), Dict(
        :rho => Matrix(rho_traj),
        :delta => REF_DELTA,
        :domain => string(typeof(domain)),
        :with_coherent => with_coherent,
        :num_qubits => 3,
        :seed => REF_SEED,
        :ntraj => REF_NTRAJ,
    ))
    println("Saved: $filename")
end

# ---------------------------------------------------------------------------
# Generate all 4 reference files
# ---------------------------------------------------------------------------
println("Generating frozen BSON reference data...")
println()

generate_dm_reference(EnergyDomain(); with_coherent=false, filename="energy_dm_reference.bson")
generate_traj_reference(EnergyDomain(); with_coherent=false, filename="energy_traj_reference.bson")
generate_dm_reference(TrotterDomain(); with_coherent=true, filename="trotter_coherent_dm_reference.bson")
generate_traj_reference(TrotterDomain(); with_coherent=true, filename="trotter_coherent_traj_reference.bson")

println()
println("All reference files generated successfully.")
