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
using Roots
using DataStructures
using SpecialFunctions: erfc
using QuadGK
using Base.Threads
using FINUFFT
using KrylovKit
using LibGit2
using LsqFit
using Dates

# --- Public API ---

# --- Lindbladian ---
export run_lindblad, construct_lindbladian
export LindbladResults
export apply_lindbladian!, apply_adjoint_lindbladian!
export krylov_spectral_gap, apply_delta_channel!

# --- Thermalize ---
export run_thermalize
export ThermalizeResults

# --- Krylov ---
export run_krylov_spectrum
export KrylovSpectrumResults

# --- Trajectory ---
export run_trajectory
export run_trajectories, run_observable_trajectories
export TrajectoryResults
export TrajectoryResult, ObservableTrajectoryResult
export step_along_trajectory!
export ConvergenceData, run_trajectories_convergence, run_trajectories_adaptive
export build_preset_trajectory_observables

# --- Diagnostics ---
export EigenDecompositionResult, FixedPointResult, DefectResult, OverlapResult,
       SzSectorLabel, MultipletGroup, ExactDiagnosticsResult,
       extract_leading_eigendata, compute_fixed_point_distance,
       compute_anti_hermitian_defect, compute_overlap_coefficients,
       compute_sz_labels, detect_multiplets, run_exact_diagnostics

# --- Common ---
export Config, AbstractSimulation, Lindbladian, Thermalize, KrylovSpectrum, Trajectory
export AbstractConstruction, KMS, GNS, DLL, with_coherent
export Workspace, LiouvillianScratch, ThermalizeScratch, KrylovScratch, TrajectoryScratch
export AbstractResults, save_result, load_result
export BohrDomain, EnergyDomain, TimeDomain, TrotterDomain
export HamHam, TrottTrott, JumpOp
export trace_distance_h, trace_distance_nh, trace_norm_h, trace_norm_nh,
       fidelity, frobenius_norm, is_density_matrix, random_density_matrix,
       hermitianize!
export gibbs_state, gibbs_state_in_eigen,
       find_ideal_heisenberg, load_hamiltonian,
       create_bohr_dict, compute_trotter_error
export pick_transition, create_alpha, create_alpha_gns, create_alpha_gauss,
       create_f, create_f_gauss, check_alpha_skew_symmetry
export B_time, B_trotter, B_bohr
export X, Y, Z, Had,
       pad_term, expm_pauli_padded, pauli_string_to_matrix,
       trotterize, group_hamiltonian_terms
export validate_config!
export oft!, time_oft!, trotter_oft!

# STAGING: estimate_spectral_gap, OverlapAnalysisResult, eigenbasis_overlap_analysis
export fit_exponential_decay, FitResult
export estimate_mixing_time, MixingTimeEstimate
# STAGING: LSIFramework, compute_LSI_alpha2

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
include("energy_domain.jl")
include("bohr_domain.jl")
include("coherent.jl")
include("jump_workers.jl")
include("trajectories.jl")
include("furnace_utensils.jl")
include("furnace.jl")
include("krylov_workspace.jl")
include("krylov_matvec.jl")
include("krylov_eigsolve.jl")
include("convergence.jl")
include("diagnostics.jl")
include("results.jl")
include("fitting.jl")
include("mixing.jl")

end