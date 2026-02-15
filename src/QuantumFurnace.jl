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

# Types: Simulation
export AbstractConfig, AbstractLiouvConfig, AbstractThermalizeConfig,
       LiouvConfig, LiouvConfigGNS, ThermalizeConfig, ThermalizeConfigGNS,
       HamHam, TrottTrott,
       HotAlgorithmResults, HotSpectralResults,
       JumpOp

# Types: Domains
export BohrDomain, EnergyDomain, TimeDomain, TrotterDomain

# Types: Trajectory
export TrajectoryFramework

# Types: Log-Sobolev
export LSIFramework

# Simulation
export run_lindbladian, run_thermalization, construct_lindbladian,
       run_trajectories, build_trajectoryframework, step_along_trajectory!

# QI Tools
export trace_distance_h, trace_distance_nh, trace_norm_h, trace_norm_nh,
       fidelity, frobenius_norm, is_density_matrix, random_density_matrix,
       hermitianize!, transform_jumps_to_basis

# Gibbs & Hamiltonian
export gibbs_state, gibbs_state_in_eigen,
       find_ideal_heisenberg, load_hamiltonian,
       create_bohr_dict, compute_trotter_error

# Transition functions & Kossakowski matrix
export pick_transition, create_alpha, create_alpha_gns, create_alpha_gauss,
       create_f, create_f_gauss, check_alpha_skew_symmetry

# Coherent terms (B operators)
export B_time, B_trotter, coherent_bohr

# Pauli & Trotter building blocks
export X, Y, Z, Had,
       pad_term, expm_pauli_padded, pauli_string_to_matrix,
       trotterize, group_hamiltonian_terms

# Config validation
export validate_config!

# Log-Sobolev bound
export compute_LSI_alpha2

# OFT (kept for debugging / pedagogy)
export oft!

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