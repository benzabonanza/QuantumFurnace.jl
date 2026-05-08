# scripts/scratch_krylov_trajectory_demo.jl
#
# qf-ev5 demo: Krylov spectral-expansion trajectory predictors.
#
# Demonstrates two new public functions in src/lindblad_action.jl:
#
#   * predict_lindbladian_trajectory(cfg::Config{Lindbladian}, ham, jumps, rho_0, t_grid; ...)
#       Reconstructs e^{tL} rho_0 on a t-grid from a single Arnoldi
#       factorization of the matrix-free apply_lindbladian!. Same
#       NamedTuple shape as `lindblad_action_integrate`, so the bi-exp
#       fitter / `estimate_mixing_time` plugs in unchanged.
#
#   * predict_channel_trajectory(cfg::Config{Thermalize}, ham, jumps, rho_0, k_grid; ...)
#       Reconstructs (Phi_delta)^k rho_0 from a single Arnoldi of the
#       FAITHFUL implemented channel. The forward matvec is byte-for-byte
#       the same per-jump :sweep substep run_thermalize uses, so the
#       reconstruction matches a stepped trajectory to ~ 1e-12 abs.
#
# Both bypass the per-step iteration cost of `lindblad_action_integrate` /
# `run_thermalize`, replacing length(t_grid) * krylov_subspace matvecs (or
# T/delta substeps) with one Arnoldi factorisation of size krylovdim.
#
# Run:  julia --project scripts/scratch_krylov_trajectory_demo.jl

using QuantumFurnace
using LinearAlgebra
using Printf

include(joinpath(dirname(@__DIR__), "test", "test_helpers.jl"))


function demo_lindbladian(n::Int, beta::Real)
    println("="^72)
    println("Lindbladian (ideal e^{tL}) at n=$n, beta=$beta CKG smooth-Metro")
    println("="^72)

    fixture_loader = if n == 3
        β -> begin
            f = make_dll_n3_system(β)
            (f.ham, f.jumps)
        end
    else
        β -> begin
            ham_path = joinpath(dirname(@__DIR__), "hamiltonians",
                                "heis_disordered_periodic_n$(n).bson")
            ham = QuantumFurnace._load_hamiltonian_bson(ham_path, β)
            jp = ([X], [Y], [Z])
            num_jumps = length(jp) * n
            jn = sqrt(num_jumps)
            jumps = JumpOp[]
            for pauli in jp, site in 1:n
                op = Matrix(pad_term(pauli, n, site)) ./ jn
                op_eb = ham.eigvecs' * op * ham.eigvecs
                push!(jumps, JumpOp(op, op_eb, op == transpose(op), op == op'))
            end
            (ham, jumps)
        end
    end
    ham, jumps = fixture_loader(beta)
    d = size(ham.data, 1)

    cfg = Config(
        sim = Lindbladian(),
        domain = EnergyDomain(),
        construction = KMS(),
        num_qubits = n,
        with_linear_combination = true,
        beta = beta, sigma = 1.0 / beta,
        a = 0.0, s = 0.25,
        num_energy_bits = 12, w0 = 0.05,
        t0 = 2π / (2^12 * 0.05),
        num_trotter_steps_per_t0 = 10,
    )

    rho_0 = Matrix{ComplexF64}(I(d) / d)
    sigma_beta = Matrix{ComplexF64}(ham.gibbs)

    gap_res = krylov_spectral_gap(cfg, ham, jumps; krylovdim=30, howmany=4, tol=1e-10)
    gap = gap_res.spectral_gap
    println(@sprintf("spectral gap = %.6e", gap))
    t_grid = collect(range(0.0, 5.0 / gap, length=41))

    # ODE-based reference
    t0_ode = time()
    res_ode = integrate_to_gibbs(cfg, ham, jumps, rho_0, t_grid;
                                  mode=:L, krylovdim=30, tol=1e-10)
    t_ode = time() - t0_ode

    # Krylov spectral
    t0_kr = time()
    res_kr = predict_lindbladian_trajectory(cfg, ham, jumps, rho_0, t_grid;
                                              krylovdim=40, tol=1e-10)
    t_kr = time() - t0_kr

    println(@sprintf("ODE     %5d matvecs  wall %.3fs", res_ode.total_matvecs, t_ode))
    println(@sprintf("Krylov  %5d matvecs  wall %.3fs   |   speedup %.1fx",
                     res_kr.total_matvecs, t_kr, t_ode / t_kr))
    max_diff = maximum(abs.(res_kr.distances .- res_ode.distances))
    println(@sprintf("max |Krylov - ODE| over trajectory: %.3e", max_diff))

    # τ_mix from each
    est_ode = estimate_mixing_time(res_ode; model=:biexp,
                                    target_epsilon=1e-3, extrapolate=true)
    est_kr = estimate_mixing_time(res_kr; model=:biexp,
                                   target_epsilon=1e-3, extrapolate=true)
    if isfinite(est_ode.mixing_time) && isfinite(est_kr.mixing_time)
        rel = abs(est_ode.mixing_time - est_kr.mixing_time) / est_ode.mixing_time
        println(@sprintf("τ_mix (ODE)    = %.4f", est_ode.mixing_time))
        println(@sprintf("τ_mix (Krylov) = %.4f   rel diff = %.3e",
                         est_kr.mixing_time, rel))
    end
    println()
end


function demo_channel(n::Int, beta::Real, delta::Real, num_steps::Int = 30000)
    println("="^72)
    println("Channel (implemented Phi_delta) at n=$n, beta=$beta, delta=$delta")
    println("="^72)

    fixture_loader = if n == 3
        β -> begin
            f = make_dll_n3_system(β)
            (f.ham, f.jumps)
        end
    else
        β -> begin
            ham_path = joinpath(dirname(@__DIR__), "hamiltonians",
                                "heis_disordered_periodic_n$(n).bson")
            ham = QuantumFurnace._load_hamiltonian_bson(ham_path, β)
            jp = ([X], [Y], [Z])
            num_jumps = length(jp) * n
            jn = sqrt(num_jumps)
            jumps = JumpOp[]
            for pauli in jp, site in 1:n
                op = Matrix(pad_term(pauli, n, site)) ./ jn
                op_eb = ham.eigvecs' * op * ham.eigvecs
                push!(jumps, JumpOp(op, op_eb, op == transpose(op), op == op'))
            end
            (ham, jumps)
        end
    end
    ham, jumps = fixture_loader(beta)
    d = size(ham.data, 1)

    cfg = Config(
        sim = Thermalize(),
        domain = BohrDomain(),
        construction = KMS(),
        num_qubits = n,
        with_linear_combination = true,
        beta = beta, sigma = 1.0 / beta,
        a = 0.0, s = 0.25,
        num_energy_bits = 12, w0 = 0.05,
        t0 = 2π / (2^12 * 0.05),
        num_trotter_steps_per_t0 = 10,
        delta = delta,
        mixing_time = num_steps * delta,
        jump_selection = :sweep,
    )

    rho_0 = Matrix{ComplexF64}(I(d) / d)

    save_every = max(1, num_steps ÷ 30)
    k_grid = collect(0:save_every:num_steps)

    t0_kr = time()
    res_kr = predict_channel_trajectory(cfg, ham, jumps, rho_0, k_grid;
                                         krylovdim=40, tol=1e-10)
    t_kr = time() - t0_kr

    t0_th = time()
    res_th = run_thermalize(jumps, cfg, ham; initial_dm=copy(rho_0),
                             save_every=save_every)
    t_th = time() - t0_th

    n_compare = min(length(res_kr.distances), length(res_th.trace_distances))
    max_abs = maximum(abs.(res_kr.distances[1:n_compare] .-
                            res_th.trace_distances[1:n_compare]))

    println(@sprintf("Krylov channel:  %3d matvecs  wall %.3fs", res_kr.total_matvecs, t_kr))
    println(@sprintf("run_thermalize:  %d steps     wall %.3fs", num_steps, t_th))
    println(@sprintf("speedup: %.0fx wall  /   max abs diff = %.3e", t_th / t_kr, max_abs))
    println(@sprintf("μ_1 = %.10f  (CPTP fixed point should be 1.0)", abs(res_kr.eigenvalues[1])))
    println(@sprintf("μ_2 = %.10f  ⇒  spectral gap (in 1/δ) = %.4f",
                     abs(res_kr.eigenvalues[2]), res_kr.spectral_gap))
    println()
end


function main()
    println()
    println("===  qf-ev5 demo: Krylov spectral-expansion trajectories  ===")
    println()

    # Lindbladian-side: scales with n.
    for n in [3, 4]
        demo_lindbladian(n, 10.0)
    end

    # Channel-side: matches run_thermalize byte-for-byte at n=3.
    demo_channel(3, 10.0, 1e-3, 30000)

    println("="^72)
    println("Done.")
end

main()
