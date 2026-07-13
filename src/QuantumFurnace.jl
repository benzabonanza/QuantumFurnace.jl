module QuantumFurnace

using Pkg
using Base
using Printf
using BSON
using Arpack
using LinearAlgebra
using SparseArrays
using Random
using Statistics: median
using ProgressMeter
using Roots
using DataStructures
using SpecialFunctions: erfc, besselj
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
export krylov_spectral_gap, apply_delta_channel!, apply_adjoint_delta_channel!

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
export SpectralModeDiagnostics, spectral_mode_diagnostics

# --- Discriminant ---
export DiscriminantBuffers, gibbs_fractional_powers, apply_discriminant!
export materialize_discriminant, materialize_discriminant!,
       hermitian_antihermitian_split, hermitian_antihermitian_split!
export DiscriminantSpectrum, discriminant_spectrum
export DBVerificationResult, verify_detailed_balance

# --- KMS geometry (qf-mto.{1,2,3}, parked diagnostic — see src/kms_geometry.jl) ---
# These compare Lindbladian magnitudes (CKG vs DLL); they are not used to
# rescale generators in mainline simulations.
export kms_inner_product, kms_norm, kms_variance, kms_dirichlet_form
export build_dense_superoperator
export spectral_gap_kms, max_dirichlet_rate_kms, intrinsic_mixing_ratio
export dissipator_one_to_one_norm_bound, dissipator_trace_alpha, hs_operator_norm
export hs_operator_norm_krylov

# --- Common ---
export Config, AbstractSimulation, Lindbladian, Thermalize, KrylovSpectrum, Trajectory
export AbstractConstruction, KMS, GNS, DLL, with_coherent
export Workspace, LiouvillianScratch, ThermalizeScratch, KrylovScratch, TrajectoryScratch
export AbstractResults, save_result, load_result
export BohrDomain, EnergyDomain, TimeDomain, TrotterDomain
export HamHam, AbstractTrotter, TrottTrott, TrotterTriple, JumpOp
export trace_distance_h, trace_distance_nh, trace_norm_h, trace_norm_nh,
       fidelity, is_density_matrix, random_density_matrix,
       hermitianize!, validate_jump_pairing
export gibbs_state, gibbs_state_in_eigen,
       build_heis_1d, build_tfim_2d, load_hamiltonian,
       create_bohr_dict, compute_trotter_error, make_trotter_for_config
# β_phys ↔ β_alg helpers (qf-6vr) — convert physical and algorithm-side
# inverse temperatures through `ham.rescaling_factor`. See docstrings in
# `src/hamiltonian.jl` for the convention.
export beta_alg, beta_phys
export pick_transition, pick_gamma_sup, create_alpha, create_alpha_gns, create_alpha_gauss,
       create_f, create_f_gauss, check_alpha_skew_symmetry
# `default_smooth_s` (src/bohr_domain.jl) intentionally not exported — internal
# helper used by production β-sweep scripts. Call as `QuantumFurnace.default_smooth_s`.
export B_time, B_trotter, B_bohr
export X, Y, Z, Had,
       pad_term, expm_pauli_padded, pauli_string_to_matrix,
       trotterize, group_hamiltonian_terms
export validate_config!
export register_t0_D, register_w0_D, register_r_D,
       register_t0_b_minus, register_w0_b_minus, register_r_b_minus,
       register_t0_b_plus, register_w0_b_plus, register_r_b_plus,
       register_M_D, register_M_b_minus, register_M_b_plus
export oft!

# --- Filters (DLL-1) ---
export AbstractFilter, GaussianFilter, DLLGaussianFilter, DLLMetropolisFilter
export time_kernel, freq_kernel, filter_time_cutoff

# --- Multi-channel DLL (qf-7go epic, parked diagnostic — see src/dll_multichannel.jl) ---
export DLLMultiChannelFilter, ShiftedSymmetricFilter, dll_multichannel_translates

# --- DLL dissipator helpers (DLL-2) ---
export dll_lindblad_op_bohr, dll_lindblad_op_time

# --- DLL coherent helpers (DLL-3) ---
export dll_coherent_op_bohr, dll_coherent_op_time
# dll_coherent_kernel_bohr is intentionally not exported — internal kernel
# (Eq. 3.5) used by the test-only legacy reference path. Reach via
# QuantumFurnace.dll_coherent_kernel_bohr if needed in scripts.

# --- DLL Kossakowski (DLL-4) ---
export dll_kossakowski_bohr

# --- Lindbladian-action ODE integrator (qf-lkb.1) ---
export lindblad_action_integrate, discriminant_action_integrate, integrate_to_gibbs, sweep_mixing_times

# --- Krylov spectral-expansion trajectory (qf-ev5) ---
export predict_lindbladian_trajectory, predict_channel_trajectory

# --- Matrix-free superoperator trace-norm distance (qf-72g) ---
export PropagatorArm, propagator_trace_distance, propagator_fixed_point_distance
export lindbladian_arm, channel_arm
# --- Slow-subspace generator-mismatch distance (qf-e4z.45) ---
export slow_subspace_generator_distance
# --- Robust fixed-point extraction + fixed-point-vs-Gibbs distance (qf-e4z.48) ---
export arm_fixed_point, fixed_point_gibbs_distance
# --- Anti-Hermitian quantum-discriminant norm: channel KMS-DB violation (qf-e4z.50) ---
export discriminant_antiherm_norm, channel_discriminant_antiherm_norm,
       lindbladian_discriminant_antiherm_norm

# --- Channel sweep harness (qf-e4z.2 / P0b) ---
export sweep_channel_mixing

# STAGING: estimate_spectral_gap, OverlapAnalysisResult, eigenbasis_overlap_analysis
# STAGING (qf-6z9.4): compute_oft_trotter_error, compute_oft_trotter_error_all_jumps moved to src/staging/errors.jl
export fit_exponential_decay, FitResult
export fit_biexponential_decay, BiexpFitResult
export estimate_mixing_time, MixingTimeEstimate
export eigenmode_mixing_time
export SimulationTimeBudget, compute_simulation_time
export TrotterStepBudget, count_trotter_steps
export RxxBudget, estimate_rxx_count, load_rxx_table

# --- Empirical scaling-law extraction (qf-now) ---
export ScalingFit, fit_scaling, predict_scaling, aicc_weights, compare_models,
       formula_string, scaling_fit_grid
# STAGING: LSIFramework, compute_LSI_alpha2

# --- Internal Implementation ---
include("constants.jl")
include("hamiltonian.jl")
include("trotter_domain.jl")
include("filters.jl")
include("structs.jl")
include("qi_tools.jl")
include("misc_tools.jl")
include("time_domain.jl")
include("nufft.jl")
include("ofts.jl")
include("energy_domain.jl")
include("bohr_domain.jl")
include("coherent.jl")
include("dll.jl")
include("jump_workers.jl")
include("trajectories.jl")
include("furnace_utensils.jl")
include("dll_multichannel.jl")
include("furnace.jl")
include("krylov_workspace.jl")
include("krylov_matvec.jl")
include("krylov_eigsolve.jl")
include("convergence.jl")
include("diagnostics.jl")
include("discriminant.jl")
include("kms_geometry.jl")
include("lindblad_action.jl")
include("superop_distance.jl")
include("results.jl")
include("fitting.jl")
include("mixing.jl")
include("scaling_fit.jl")
# `errors.jl` retired in qf-6z9.4 — moved to `src/staging/errors.jl` (commented out)
include("simulation_time.jl")

# Classical Gibbs-sampling baseline (sign-free SSE QMC) — the directly-comparable classical
# competitor to the KMS Lindbladian sampler (qf-h23). Self-contained sub-module; re-export its
# public API into the QuantumFurnace namespace.
include("classical_qmc.jl")
using .ClassicalQMC
export SSEResult, run_sse, build_sse_heis_model, build_sse_tfim_model,
    sse_exact_reference, sse_reconstruction_error

end