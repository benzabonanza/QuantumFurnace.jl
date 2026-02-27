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
using SharedArrays
using FINUFFT
using LsqFit
using KrylovKit
using LibGit2
using Dates

# --- Public API ---

# Types: Simulation
export Config,
       AbstractSimulation, Lindbladian, Thermalize, KrylovSpectrum, Trajectory,
       AbstractConstruction, KMS, GNS, DLL,
       with_coherent,
       HamHam, TrottTrott,
       LindbladianResult, DMSimulationResult,
       JumpOp

# Workspace (Phase 35)
export Workspace, LiouvillianScratch, ThermalizeScratch, KrylovScratch, TrajectoryScratch

# Types: Domains
export BohrDomain, EnergyDomain, TimeDomain, TrotterDomain

# Types: Log-Sobolev
# export LSIFramework, compute_LSI_alpha2

# Trajectory
export TrajectoryResult, ObservableTrajectoryResult, step_along_trajectory!, run_observable_trajectories

# Results / Data persistence
export ExperimentResult, save_experiment, load_experiment

# New typed Results (Phase 36)
export AbstractResults, LindbladResults, ThermalizeResults, KrylovSpectrumResults, TrajectoryResults
export save_result, load_result

# Convergence tracking
export ConvergenceData, run_trajectories_convergence, run_trajectories_adaptive, build_preset_trajectory_observables

# Fitting
export fit_exponential_decay, FitResult

# Gap estimation
export SpectralGapResult, estimate_spectral_gap, OverlapAnalysisResult, eigenbasis_overlap_analysis

# Diagnostics (Phase 26)
export EigenDecompositionResult, FixedPointResult, DefectResult, OverlapResult,
       SzSectorLabel, MultipletGroup, ExactDiagnosticsResult,
       extract_leading_eigendata, compute_fixed_point_distance,
       compute_anti_hermitian_defect, compute_overlap_coefficients,
       compute_sz_labels, detect_multiplets, run_exact_diagnostics

# Simulation
export run_lindbladian, run_thermalization, construct_lindbladian,
       run_trajectories

# New entry points (Phase 36)
export run_lindblad, run_thermalize, run_krylov_spectrum, run_trajectory

# Krylov matvec (Phase 27)
export apply_lindbladian!, apply_adjoint_lindbladian!

# Krylov eigsolve (Phase 29)
export KrylovGapResult, krylov_spectral_gap, apply_delta_channel!

# QI Tools
export trace_distance_h, trace_distance_nh, trace_norm_h, trace_norm_nh,
       fidelity, frobenius_norm, is_density_matrix, random_density_matrix,
       hermitianize!

# Gibbs & Hamiltonian & Trotter
export gibbs_state, gibbs_state_in_eigen,
       find_ideal_heisenberg, load_hamiltonian,
       create_bohr_dict, compute_trotter_error

# Transition functions & Kossakowski matrices
export pick_transition, create_alpha, create_alpha_gns, create_alpha_gauss,
       create_f, create_f_gauss, check_alpha_skew_symmetry

# Coherent terms (B operators)
export B_time, B_trotter, B_bohr

# Pauli & Trotter building blocks
export X, Y, Z, Had,
       pad_term, expm_pauli_padded, pauli_string_to_matrix,
       trotterize, group_hamiltonian_terms

# Config validation
export validate_config!

# OFT (kept for debugging / pedagogy)
export oft!, time_oft!, trotter_oft!

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
include("krylov_workspace.jl")
include("krylov_matvec.jl")
include("krylov_eigsolve.jl")
include("log_sobolev.jl")
include("convergence.jl")
include("fitting.jl")
include("gap_estimation.jl")
include("diagnostics.jl")
include("results.jl")

end