using Test
using QuantumFurnace
using LinearAlgebra

# Shared fixtures and constants
include("test_helpers.jl")

@testset "QuantumFurnace.jl" begin
    include("test_aqua.jl")
    include("test_register_validation.jl")
    include("test_register_independence.jl")
    include("test_hamiltonian.jl")
    include("test_compilation.jl")
    include("test_trotter_caches.jl")
    include("test_trajectory_fixes.jl")
    include("test_cptp.jl")
    include("test_dm_detailed_balance.jl")
    include("test_dm_scaling.jl")
    include("test_regression.jl")
    include("test_allocation.jl")
    include("test_workspace_independence.jl")
    include("test_threading.jl")
    include("test_qf_6af_construction_threading.jl")
    include("test_gns_trajectory.jl")
    include("test_results.jl")
    include("test_convergence.jl")
    include("test_save_every.jl")
    include("test_observable_trajectories.jl")
    include("test_diagnostics.jl")
    include("test_dll_filter.jl")
    include("test_dll_multichannel_filter.jl")
    include("test_dll_multichannel_bohr.jl")
    include("test_dll_multichannel_time.jl")
    include("test_dll_multichannel_simulator.jl")
    include("test_dll_dissipator.jl")
    include("test_dll_coherent.jl")
    include("test_dll_kossakowski.jl")
    include("test_dll_kms_db.jl")
    include("test_discriminant.jl")
    include("test_kms_geometry.jl")
    include("test_lindblad_action.jl")
    include("test_predict_lindbladian.jl")
    include("test_predict_channel.jl")
    include("test_sweep_channel_mixing.jl")
    include("test_krylov_matvec.jl")
    include("test_krylov_eigsolve.jl")
    include("test_krylov_crossvalidation.jl")
    include("test_fitting.jl")
    include("test_mixing.jl")
    include("test_simulation_time.jl")
    include("test_smooth_metro_eta.jl")
    include("test_gqsp_config.jl")
    include("test_gqsp_polynomial.jl")
    include("test_gqsp_thermalize.jl")
    include("test_jump_selection.jl")
    include("test_gamma_norm_invariance.jl")
    include("test_non_hermitian_jumps.jl")
    include("test_validate_jump_pairing.jl")
    include("test_qf_bm1_kwarg_threading.jl")
    include("test_oft_prefactors_energy.jl")

    # Trajectory validation (slow: ~minutes, gated behind env variable)
    if get(ENV, "QUANTUMFURNACE_FULL_TESTS", "false") == "true"
        include("trajectory_validation/run_trajectory_validation.jl")
        include("trajectory_validation/run_convergence_tests.jl")
    else
        @info "Skipping trajectory validation tests (set QUANTUMFURNACE_FULL_TESTS=true to run)"
    end
end
