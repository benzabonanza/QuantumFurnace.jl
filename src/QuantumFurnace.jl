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
using ClusterManagers
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
export LiouvConfig, ThermalizeConfig, HamHam, TrottTrott, HotAlgorithmResults, HotSpectralResults, JumpOp,
       BohrDomain, EnergyDomain, TimeDomain, TrotterDomain, LindbladWorkspace, KrausFramework, LSIFramework,
       OFTCaches, NUFFTPrefactors
export run_liouvillian, run_thermalization, construct_liouvillian, B_time, B_trotter, coherent_bohr
export generate_filename, validate_config!, create_trotter, compute_trotter_error, gibbs_state, gibbs_state_in_eigen,
       create_bohr_dict, pad_term, pick_transition, create_hamham, find_ideal_heisenberg, create_alpha, expm_pauli_padded, 
       finalize_hamham, load_hamiltonian, oft!, prepare_oft_nufft_prefactors, prefactor_view
# Quantum Trajectory
export krausframework, step_along_the_trajectory!, evolve_along_trajectory, construct_gksl_lindbladian,
       apply_jump_contribution!, apply_lindbladian_dagger!, apply_lindbladian!, precompute_B, precompute_R, 
       precompute_kraus_jumps, precompute_data, verify_completeness
# Log Sobolev bound
export compute_LSI_alpha2
export X, Y, Z, Had

# --- Internal Implementation ---
include("constants.jl")
include("structs.jl")
include("hamiltonian.jl")
include("qi_tools.jl")
include("misc_tools.jl")
include("nufft.jl")
include("ofts.jl")
include("errors.jl")
include("jump_workers.jl")
include("coherent.jl")
include("bohr_domain.jl")
include("energy_domain.jl")
include("time_domain.jl")
include("timelike_tools.jl")
include("trotter_domain.jl")
include("trajectories.jl")
include("furnace_utensils.jl")
include("furnace.jl")
include("log_sobolev.jl")
include("log_sobolev_manopt.jl")
include("linearmaps_liouv.jl")

end