module QuantumFurnace

using Pkg
using Base
using Printf
using BSON
using Arpack
using LinearAlgebra
using SparseArrays
using Random
using ProgressMeter
using Distributed
using Roots
using DataStructures
using SpecialFunctions: erfc
using QuadGK
using Optim
using Base.Threads
using LinearMaps
using SharedArrays
using FINUFFT

# --- Public API ---
export AbstractConfig, AbstractLiouvConfig, AbstractThermalizeConfig,
       LiouvConfig, LiouvConfigGNS, ThermalizeConfig, ThermalizeConfigGNS, HamHam, TrottTrott, HotAlgorithmResults, HotSpectralResults, JumpOp,
       BohrDomain, EnergyDomain, TimeDomain, TrotterDomain, LindbladWorkspace, LSIFramework,
       OFTCaches, NUFFTPrefactors
export run_lindbladian, run_thermalization, construct_lindbladian, B_time, B_trotter, coherent_bohr
export generate_filename, validate_config!, create_trotter, compute_trotter_error, gibbs_state, gibbs_state_in_eigen,
       create_bohr_dict, pad_term, pick_transition, create_hamham, find_ideal_heisenberg, create_alpha, 
       expm_pauli_padded, finalize_hamham, load_hamiltonian, oft!, prepare_oft_nufft_prefactors, prefactor_view,
       precompute_coherent_terms, precompute_coherent_total_B
export create_alpha_gns
# Quantum Trajectory
export TrajectoryFramework, TrajectoryWorkspace, build_trajectoryframework, step_along_trajectory!, evolve_along_trajectory,
       evolve_and_measure_along_trajectory, run_trajectories, precompute_R,
       precompute_data
export KrausScratch, jump_contribution!
# Log Sobolev bound
export compute_LSI_alpha2
export X, Y, Z, Had

# --- Internal Implementation ---
include("constants.jl")
include("hamiltonian.jl")
include("trotter_domain.jl")
include("structs.jl")
include("qi_tools.jl")
include("misc_tools.jl")
include("time_domain.jl")
include("nufft.jl")
include("ofts.jl")
include("errors.jl")
include("kraus.jl")
include("energy_domain.jl")
include("bohr_domain.jl")
include("coherent.jl")
include("jump_workers.jl")
include("trajectories.jl")
include("furnace_utensils.jl")
include("furnace.jl")
include("log_sobolev.jl")
include("log_sobolev_manopt.jl")
include("linearmaps_liouv.jl")

end