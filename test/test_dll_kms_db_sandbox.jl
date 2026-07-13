# test/test_dll_kms_db_sandbox.jl
#
# Sandbox shadow of test_dll_kms_db.jl::(j) (qf-x56.5). The heavy (j)
# testset is gated NO_SANDBOX and runs a multi-t0 convergence sweep at
# Nt = 4096 (r_D = 12) that does not fit the sandbox RSS envelope even
# when isolated in its own file. This shadow keeps a single regression-
# guard cell at the same β = 10 fixture but at r_D = 11 (Nt = 2048,
# 4× less NUFFT memory than r_D = 12), which still extends coverage one
# bit above the existing inline (i) testset's r_D = 10 default.
#
# Why a separate file: inlined inside test_dll_kms_db.jl the (a)..(i)
# testsets accumulate Nt = 1024 NUFFT workspaces; layering a fresh
# Nt = 2048 pair on top crosses the 1.5 GB heap-size-hint, so we let
# the inter-file `GC.gc(true)` in runtests.jl clear the prior state.
#
# Threshold: 1e-3. PHYSICS CHECK: at fixed w0 = 0.05 the time-grid extent
# `t_max = π / w0 ≈ 62.8` is Nt-independent — moving r_D = 12 → r_D = 11
# halves Nt but does NOT shrink t_max (only t0 = 2π/(Nt·w0) doubles, with
# the same number of grid points spanning the same [0, t_max] window).
# So the Metropolis time_kernel tail |f(t_max)| is unchanged; only the
# Gaussian-quadrature precision shifts (Nt halved ⇒ NUFFT noise floor up
# by ≈ 4×). At r_D = 12 the heavy (j) measured `err10_coarse ≈ 2e-5`;
# at r_D = 11 the same fixture sits in the same regime (empirical: ~3e-5).
# Threshold 1e-3 leaves a ~30× margin against any regression that loosens
# the tail or breaks the half-grid fold.

using LinearAlgebra: opnorm
using Test
using QuantumFurnace


function _sandbox_dll_meta_cfg_t(beta::Real, num_energy_bits::Int)
    Config(;
        sim = Lindbladian(), domain = TimeDomain(), construction = DLL(),
        num_qubits = 3, with_linear_combination = true,
        beta = beta, sigma = 1.0 / beta, a = beta / 30, s = 0.4,
        num_energy_bits = num_energy_bits,
        t0 = 2π / (2^num_energy_bits * 0.05),
        num_trotter_steps_per_t0 = 10,
        filter = DLLMetropolisFilter(beta; S = 2.0),
    )
end

function _sandbox_dll_meta_cfg_b(beta::Real, num_energy_bits::Int)
    Config(;
        sim = Lindbladian(), domain = BohrDomain(), construction = DLL(),
        num_qubits = 3, with_linear_combination = true,
        beta = beta, sigma = 1.0 / beta, a = beta / 30, s = 0.4,
        num_energy_bits = num_energy_bits,
        t0 = 2π / (2^num_energy_bits * 0.05),
        num_trotter_steps_per_t0 = 10,
        filter = DLLMetropolisFilter(beta; S = 2.0),
    )
end


@testset "DLL KMS-DB (j-sb) sandbox shadow (qf-x56.5)" begin
    @testset "(j-sb) DLL Metropolis Bohr ↔ Time @ β=10, r_D=11" begin
        beta = 10.0
        r_D = 11
        sys = make_dll_n3_system(beta)
        L_b = Matrix(construct_lindbladian(sys.jumps,
                                            _sandbox_dll_meta_cfg_b(beta, r_D),
                                            sys.ham))
        L_t = Matrix(construct_lindbladian(sys.jumps,
                                            _sandbox_dll_meta_cfg_t(beta, r_D),
                                            sys.ham))
        err = opnorm(L_b - L_t)
        # Drop the two dense 64×64 matrices + their NUFFT workspaces
        # immediately so the per-file Δ RSS reported by runtests.jl is small.
        L_b = nothing; L_t = nothing; sys = nothing
        GC.gc()

        @test err < 1e-3
        @info "(j-sb) DLL Metropolis Bohr ↔ Time" β=beta r_D err threshold=1e-3
    end
end
