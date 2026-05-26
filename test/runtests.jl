using Test
using QuantumFurnace
using LinearAlgebra

# Shared fixtures and constants
include("test_helpers.jl")

# ----------------------------------------------------------------------------
# SANDBOX vs NO_SANDBOX split (qf-5nz).
#
# Test files in `SANDBOX_FILES` must run within the 3.5 GB / few-minute
# sandbox container (per the project's sandbox profile, see
# `.claude/rules/julia-code.md`). They form the default `Pkg.test()` suite.
#
# `NO_SANDBOX_FILES` are heavier tests intentionally kept out of the default
# run because their physics-meaningful assertions cannot be tightened
# without compromising correctness within the sandbox's resource envelope.
# Trajectory-validation tests live in `trajectory_validation/` and are
# gated by the same env switch.
#
# To run the full suite (NO_SANDBOX + trajectory_validation), set:
#   QUANTUMFURNACE_FULL_TESTS=true julia --project -e 'using Pkg; Pkg.test()'
#
# Some individual subtests inside SANDBOX files are also gated by
# `QUANTUMFURNACE_FULL_TESTS` — see `test_dll_kms_db.jl::(j)` and
# `test_lindblad_action.jl::(i)` for the canonical examples.
# ----------------------------------------------------------------------------

const SANDBOX_FILES = String[
    "test_aqua.jl",
    "test_register_validation.jl",
    "test_register_independence.jl",
    "test_hamiltonian.jl",
    "test_build_hamiltonian_extra.jl",
    "test_boundary_conditions.jl",
    "test_compilation.jl",
    "test_trotter_caches.jl",
    "test_trajectory_fixes.jl",
    "test_cptp.jl",
    "test_dm_detailed_balance.jl",
    "test_dm_scaling.jl",
    "test_regression.jl",
    "test_allocation.jl",
    "test_workspace_independence.jl",
    "test_threading.jl",
    "test_qf_6af_construction_threading.jl",
    "test_gns_trajectory.jl",
    "test_results.jl",
    "test_convergence.jl",
    "test_save_every.jl",
    "test_observable_trajectories.jl",
    "test_diagnostics.jl",
    "test_spectral_mode_diagnostics.jl",
    "test_dll_filter.jl",
    "test_dll_multichannel_filter.jl",
    "test_dll_multichannel_bohr.jl",
    "test_dll_multichannel_time.jl",
    "test_dll_multichannel_simulator.jl",
    "test_dll_dissipator.jl",
    "test_dll_coherent.jl",
    "test_dll_kossakowski.jl",
    "test_dll_kms_db.jl",
    "test_dll_kms_db_sandbox.jl",
    "test_discriminant.jl",
    "test_predict_workspace_reuse.jl",
    "test_sweep_channel_mixing.jl",
    "test_krylov_matvec.jl",
    "test_krylov_eigsolve.jl",
    "test_krylov_crossvalidation.jl",
    "test_fitting.jl",
    "test_mixing.jl",
    "test_scaling_fit.jl",
    "test_simulation_time.jl",
    "test_smooth_metro_eta.jl",
    "test_gqsp_config.jl",
    "test_gqsp_polynomial.jl",
    "test_gqsp_thermalize.jl",
    "test_jump_selection.jl",
    "test_gamma_norm_invariance.jl",
    "test_non_hermitian_jumps.jl",
    "test_validate_jump_pairing.jl",
    "test_qf_bm1_kwarg_threading.jl",
    "test_qf_sta_b_bohr_cache.jl",
    "test_beta_phys_conversion.jl",
    "test_beta_phys_sweep.jl",
    # qf-x56: sandbox shadows of the NO_SANDBOX files.
    "test_kms_geometry_sandbox.jl",
    "test_lindblad_action_sandbox.jl",
    "test_predict_sandbox.jl",
    "test_faithful_apply_delta_channel_sandbox.jl",
    "test_trajectory_validation_sandbox.jl",
    # qf-0fv: adversarial gating test for predict_*_trajectory's
    # compute_true_gap kwarg. n=3 only, ~12s, ~50 MB.
    "test_qf_0fv_verifier.jl",
]

# Whole-file NO_SANDBOX test list. Heavy files that don't fit the cumulative
# sandbox envelope even after the per-file `GC.gc(true)` between includes.
# Run them with `QUANTUMFURNACE_FULL_TESTS=true julia --project -e 'using Pkg; Pkg.test()'`
# outside the sandbox.
const NO_SANDBOX_FILES = String[
    # Each of these builds large dense Lindbladians, large NUFFT-prefactor
    # workspaces (num_energy_bits=12 ⇒ N=4096), or both, at multiple (n, β)
    # cells. In isolation they fit (~1 GB peak), but combined with the ~30
    # preceding files they push cumulative RSS past the 3.5 GB sandbox cap.
    # The β_phys / β_alg correctness check that exercises these paths is
    # instead kept in `test_beta_phys_conversion.jl::(f)`, sized for sandbox.
    "test_kms_geometry.jl",
    "test_lindblad_action.jl",
    "test_predict_lindbladian.jl",
    "test_predict_channel.jl",
    "test_faithful_apply_delta_channel.jl",
]

const RUN_FULL = get(ENV, "QUANTUMFURNACE_FULL_TESTS", "false") == "true"

@testset "QuantumFurnace.jl" begin
    for (i, f) in enumerate(SANDBOX_FILES)
        rss_before = round(Int, Sys.maxrss()/1024^2)
        t = @elapsed include(f)
        # Force a full GC between test files. Each file accumulates large
        # NUFFT working buffers (FINUFFT plans, Lindbladian dense matrices)
        # that Julia's pool would keep around indefinitely. Without this
        # cumulative pressure pushes the suite over the sandbox cap, even
        # though every individual file fits comfortably (qf-5nz).
        GC.gc(true)
        rss_after = round(Int, Sys.maxrss()/1024^2)
        println(stderr, "[", lpad(i, 2), "/", length(SANDBOX_FILES), "] ",
                rpad(f, 45), " rss=", lpad(rss_after, 5), "MB ",
                "Δ=", lpad(rss_after - rss_before, 5), "MB  ",
                round(t, digits=1), "s")
    end

    if RUN_FULL
        for f in NO_SANDBOX_FILES
            include(f)
            GC.gc(true)
        end
        include("trajectory_validation/run_trajectory_validation.jl")
        include("trajectory_validation/run_convergence_tests.jl")
    else
        if !isempty(NO_SANDBOX_FILES)
            @info "Skipping NO_SANDBOX tests (set QUANTUMFURNACE_FULL_TESTS=true to run)" n=length(NO_SANDBOX_FILES) files=NO_SANDBOX_FILES
        end
        @info "Skipping trajectory validation tests (set QUANTUMFURNACE_FULL_TESTS=true to run)"
    end
end
