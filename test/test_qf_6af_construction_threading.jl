# test_qf_6af_construction_threading.jl
#
# Regression tests for construction-time threading added under qf-6af. Each
# threaded helper is invoked directly and its result compared against a hand-
# rolled serial reference *within the same Julia process*. The same-process
# comparison isolates the threading correctness question (chunk-reduction
# accumulation order) from the orthogonal cross-process eigendecomposition
# phase ambiguity issue (`eigen()` of a hermitian matrix returns eigvecs with
# arbitrary sign/phase).
#
# Coverage:
#   • _accumulate_R_total_threaded_energy!   (Workspace EnergyDomain)
#   • _accumulate_R_total_threaded_timetrot! (Workspace Time + TrotterDomain)
#   • _accumulate_R_total_threaded_bohr!     (Workspace BohrDomain)
#   • _accumulate_R_total_dll_chunk!         (Workspace DLL BohrDomain)
#   • _b_time_inner_threaded / _b_time_outer_threaded   (B_time)
#   • _b_trotter_inner_threaded / _b_trotter_outer_threaded (B_trotter)
#   • _B_bohr_threaded                       (BohrDomain coherent precompute)
#
# All threaded paths are entered when `Threads.nthreads() > 1` and `n_work >=
# OMEGA_THREAD_THRESHOLD = 10`; the test fixtures comfortably meet both.

using LinearAlgebra
using Random

@testset "qf-6af construction-time threading" begin
    # ------------------------------------------------------------
    # (a) R_total accumulation — EnergyDomain (Hermitian + non-Hermitian)
    # ------------------------------------------------------------
    @testset "(a) R_total threaded vs serial: EnergyDomain" begin
        if Threads.nthreads() > 1
            cfg = make_config(Lindbladian(), EnergyDomain(); construction=KMS())
            precomp = QuantumFurnace._precompute_data(cfg, TEST_HAM)

            for (jumps, name) in [
                (TEST_JUMPS, "Hermitian"),
                (begin
                    rng = MersenneTwister(123)
                    raw = randn(rng, ComplexF64, DIM, DIM) ./ sqrt(DIM)
                    JumpOp[JumpOp(raw, TEST_HAM.eigvecs' * raw * TEST_HAM.eigvecs,
                        false, false)]
                end, "non-Hermitian"),
            ]
                jump_eigenbases = [Matrix{ComplexF64}(j.in_eigenbasis) for j in jumps]
                jump_hermitian = [j.hermitian for j in jumps]

                R_threaded = zeros(ComplexF64, DIM, DIM)
                QuantumFurnace._accumulate_R_total!(R_threaded, jump_eigenbases,
                    jump_hermitian, precomp, cfg, TEST_HAM)

                # Hand-rolled serial reference (replicates the inline serial body
                # of `_accumulate_R_total!` for EnergyDomain).
                R_serial = zeros(ComplexF64, DIM, DIM)
                jump_oft = zeros(ComplexF64, DIM, DIM)
                LdagL = zeros(ComplexF64, DIM, DIM)
                inv_4sigma2 = 1.0 / (4 * cfg.sigma^2)
                pf = precomp.oft_domain_prefactor * precomp.gamma_norm_factor
                for (k, eigenbasis) in enumerate(jump_eigenbases)
                    is_herm = jump_hermitian[k]
                    if is_herm
                        for w_raw in precomp.energy_labels
                            w_raw > 1e-12 && continue
                            w = abs(w_raw)
                            QuantumFurnace.oft!(jump_oft, eigenbasis,
                                TEST_HAM.bohr_freqs, w, inv_4sigma2)
                            r2 = pf * precomp.transition(w)
                            mul!(LdagL, jump_oft', jump_oft)
                            R_serial .+= r2 .* LdagL
                            if w > 1e-12
                                r2n = pf * precomp.transition(-w)
                                mul!(LdagL, jump_oft, jump_oft')
                                R_serial .+= r2n .* LdagL
                            end
                        end
                    else
                        for w in precomp.energy_labels
                            QuantumFurnace.oft!(jump_oft, eigenbasis,
                                TEST_HAM.bohr_freqs, w, inv_4sigma2)
                            r2 = pf * precomp.transition(w)
                            mul!(LdagL, jump_oft', jump_oft)
                            R_serial .+= r2 .* LdagL
                        end
                    end
                end

                err = norm(R_threaded .- R_serial)
                rel = err / max(norm(R_serial), 1.0)
                @test rel < 1e-12
                @info "qf-6af R_total energy" path=name err=err rel=rel
            end
        else
            @info "Skipping qf-6af R_total energy test (nthreads=$(Threads.nthreads()))"
            @test true
        end
    end

    # ------------------------------------------------------------
    # (b) R_total — TimeDomain (NUFFT lookup path)
    # ------------------------------------------------------------
    @testset "(b) R_total threaded vs serial: TimeDomain" begin
        if Threads.nthreads() > 1
            cfg = make_config(Lindbladian(), TimeDomain(); construction=KMS())
            precomp = QuantumFurnace._precompute_data(cfg, TEST_HAM)
            jump_eigenbases = [Matrix{ComplexF64}(j.in_eigenbasis) for j in TEST_JUMPS]
            jump_hermitian = [j.hermitian for j in TEST_JUMPS]

            R_threaded = zeros(ComplexF64, DIM, DIM)
            QuantumFurnace._accumulate_R_total!(R_threaded, jump_eigenbases,
                jump_hermitian, precomp, cfg, TEST_HAM)

            R_serial = zeros(ComplexF64, DIM, DIM)
            jump_oft = zeros(ComplexF64, DIM, DIM)
            LdagL = zeros(ComplexF64, DIM, DIM)
            pf = precomp.oft_domain_prefactor * precomp.gamma_norm_factor
            for (k, eigenbasis) in enumerate(jump_eigenbases)
                is_herm = jump_hermitian[k]
                if is_herm
                    for w_raw in precomp.energy_labels
                        w_raw > 1e-12 && continue
                        w = abs(w_raw)
                        nufft_pf = QuantumFurnace._prefactor_view(
                            precomp.oft_nufft_prefactors, w)
                        @. jump_oft = eigenbasis * nufft_pf
                        r2 = pf * precomp.transition(w)
                        mul!(LdagL, jump_oft', jump_oft)
                        R_serial .+= r2 .* LdagL
                        if w > 1e-12
                            r2n = pf * precomp.transition(-w)
                            mul!(LdagL, jump_oft, jump_oft')
                            R_serial .+= r2n .* LdagL
                        end
                    end
                else
                    for (li, w) in enumerate(precomp.energy_labels)
                        nufft_pf = @view precomp.oft_nufft_prefactors.data[:, :, li]
                        @. jump_oft = eigenbasis * nufft_pf
                        r2 = pf * precomp.transition(w)
                        mul!(LdagL, jump_oft', jump_oft)
                        R_serial .+= r2 .* LdagL
                    end
                end
            end

            err = norm(R_threaded .- R_serial)
            rel = err / max(norm(R_serial), 1.0)
            @test rel < 1e-12
            @info "qf-6af R_total time" err=err rel=rel
        else
            @test true
        end
    end

    # ------------------------------------------------------------
    # (c) R_total — TrotterDomain (NUFFT path on Trotter eigenbasis)
    # ------------------------------------------------------------
    @testset "(c) R_total threaded vs serial: TrotterDomain" begin
        if Threads.nthreads() > 1
            cfg = make_config(Lindbladian(), TrotterDomain(); construction=KMS())
            precomp = QuantumFurnace._precompute_data(cfg, TEST_TROTTER)
            jump_eigenbases = [Matrix{ComplexF64}(j.in_eigenbasis)
                               for j in TEST_TROTTER_JUMPS]
            jump_hermitian = [j.hermitian for j in TEST_TROTTER_JUMPS]

            R_threaded = zeros(ComplexF64, DIM, DIM)
            QuantumFurnace._accumulate_R_total!(R_threaded, jump_eigenbases,
                jump_hermitian, precomp, cfg, TEST_TROTTER)

            R_serial = zeros(ComplexF64, DIM, DIM)
            jump_oft = zeros(ComplexF64, DIM, DIM)
            LdagL = zeros(ComplexF64, DIM, DIM)
            pf = precomp.oft_domain_prefactor * precomp.gamma_norm_factor
            for (k, eigenbasis) in enumerate(jump_eigenbases)
                is_herm = jump_hermitian[k]
                if is_herm
                    for w_raw in precomp.energy_labels
                        w_raw > 1e-12 && continue
                        w = abs(w_raw)
                        nufft_pf = QuantumFurnace._prefactor_view(
                            precomp.oft_nufft_prefactors, w)
                        @. jump_oft = eigenbasis * nufft_pf
                        r2 = pf * precomp.transition(w)
                        mul!(LdagL, jump_oft', jump_oft)
                        R_serial .+= r2 .* LdagL
                        if w > 1e-12
                            r2n = pf * precomp.transition(-w)
                            mul!(LdagL, jump_oft, jump_oft')
                            R_serial .+= r2n .* LdagL
                        end
                    end
                else
                    for (li, w) in enumerate(precomp.energy_labels)
                        nufft_pf = @view precomp.oft_nufft_prefactors.data[:, :, li]
                        @. jump_oft = eigenbasis * nufft_pf
                        r2 = pf * precomp.transition(w)
                        mul!(LdagL, jump_oft', jump_oft)
                        R_serial .+= r2 .* LdagL
                    end
                end
            end

            err = norm(R_threaded .- R_serial)
            rel = err / max(norm(R_serial), 1.0)
            @test rel < 1e-12
            @info "qf-6af R_total trot" err=err rel=rel
        else
            @test true
        end
    end

    # ------------------------------------------------------------
    # (d) R_total — BohrDomain
    # ------------------------------------------------------------
    @testset "(d) R_total threaded vs serial: BohrDomain" begin
        if Threads.nthreads() > 1
            cfg = make_config(Lindbladian(), BohrDomain(); construction=KMS())
            precomp = QuantumFurnace._precompute_data(cfg, TEST_HAM)
            jump_eigenbases = [Matrix{ComplexF64}(j.in_eigenbasis) for j in TEST_JUMPS]
            jump_hermitian = [j.hermitian for j in TEST_JUMPS]

            R_threaded = zeros(ComplexF64, DIM, DIM)
            QuantumFurnace._accumulate_R_total!(R_threaded, jump_eigenbases,
                jump_hermitian, precomp, cfg, TEST_HAM)

            R_serial = zeros(ComplexF64, DIM, DIM)
            jump_oft = zeros(ComplexF64, DIM, DIM)
            A_nu2_dag = zeros(ComplexF64, DIM, DIM)
            bohr_keys = collect(keys(TEST_HAM.bohr_dict))
            for (k, eigenbasis) in enumerate(jump_eigenbases)
                for nu_2 in bohr_keys
                    @. jump_oft = precomp.alpha(TEST_HAM.bohr_freqs, nu_2) * eigenbasis
                    fill!(A_nu2_dag, 0)
                    indices = TEST_HAM.bohr_dict[nu_2]
                    @inbounds for idx in indices
                        i = idx[1]; j = idx[2]
                        A_nu2_dag[j, i] = conj(eigenbasis[i, j])
                    end
                    mul!(R_serial, A_nu2_dag, jump_oft, precomp.gamma_norm_factor, 1.0)
                end
            end

            err = norm(R_threaded .- R_serial)
            rel = err / max(norm(R_serial), 1.0)
            @test rel < 1e-12
            @info "qf-6af R_total bohr" err=err rel=rel
        else
            @test true
        end
    end

    # ------------------------------------------------------------
    # (e) DLL R_total + per-jump operators — BohrDomain
    # ------------------------------------------------------------
    @testset "(e) DLL R_total + per-jump ordering" begin
        if Threads.nthreads() > 1
            beta = 10.0
            cfg = Config(;
                sim = Lindbladian(),
                domain = BohrDomain(),
                construction = DLL(),
                num_qubits = 3,
                with_linear_combination = true,
                beta = beta,
                sigma = 1.0/beta,
                a = beta/30.0,
                s = 0.4,
                num_energy_bits = 12,
                t0 = 2pi / (2^12 * 0.05),
                num_trotter_steps_per_t0 = 10,
                filter = DLLGaussianFilter(beta),
            )
            sys = make_dll_n3_system(beta)
            precomp = QuantumFurnace._precompute_data(cfg, sys.ham)

            R_threaded = zeros(ComplexF64, N3_DIM, N3_DIM)
            dll_threaded = Vector{Matrix{ComplexF64}}()
            QuantumFurnace._accumulate_R_total_dll!(
                R_threaded, dll_threaded, sys.jumps,
                precomp, cfg, sys.ham)

            R_serial = zeros(ComplexF64, N3_DIM, N3_DIM)
            dll_serial = Vector{Matrix{ComplexF64}}()
            for jump in sys.jumps
                L_or_Ls = QuantumFurnace.dll_lindblad_op_bohr(jump, sys.ham, precomp.filter)
                if L_or_Ls isa AbstractMatrix
                    L_a = Matrix{ComplexF64}(L_or_Ls)
                    push!(dll_serial, L_a)
                    mul!(R_serial, L_a', L_a, 1.0, 1.0)
                else
                    for L_one in L_or_Ls
                        L_a = Matrix{ComplexF64}(L_one)
                        push!(dll_serial, L_a)
                        mul!(R_serial, L_a', L_a, 1.0, 1.0)
                    end
                end
            end

            @test length(dll_threaded) == length(dll_serial)
            for k in eachindex(dll_serial)
                @test isapprox(dll_threaded[k], dll_serial[k]; atol=1e-15)
            end

            err = norm(R_threaded .- R_serial)
            rel = err / max(norm(R_serial), 1.0)
            @test rel < 1e-12
            @info "qf-6af R_total dll" err=err rel=rel n_ops=length(dll_serial)
        else
            @test true
        end
    end

    # ------------------------------------------------------------
    # (f) B_time threaded (inner τ × jumps + outer t)
    # ------------------------------------------------------------
    @testset "(f) B_time threaded vs serial" begin
        if Threads.nthreads() > 1
            cfg = make_config(Lindbladian(), TimeDomain(); construction=KMS())
            precomp = QuantumFurnace._precompute_data(cfg, TEST_HAM)
            B_threaded = QuantumFurnace.B_time(TEST_JUMPS, TEST_HAM,
                precomp.b_minus, precomp.b_plus,
                QuantumFurnace.register_t0_b_minus(cfg),
                QuantumFurnace.register_t0_b_plus(cfg),
                cfg.beta, cfg.sigma)

            # Serial reference: explicit nested loop matching the inline body.
            d = DIM
            CT = ComplexF64
            eigvals = TEST_HAM.eigvals

            b_plus_summand = zeros(CT, d, d)
            diag_u  = Vector{CT}(undef, d)
            diag_u2 = Vector{CT}(undef, d)
            tmp     = Matrix{CT}(undef, d, d)
            M       = Matrix{CT}(undef, d, d)
            for tau in keys(precomp.b_plus)
                t_tau = tau * cfg.beta
                @. diag_u  = exp(1im * eigvals * t_tau)
                @. diag_u2 = exp(-2im * eigvals * t_tau)
                diag_u_row = transpose(diag_u)
                for jump_a in TEST_JUMPS
                    jump_eig = jump_a.in_eigenbasis
                    @. tmp = diag_u2 * jump_eig
                    mul!(M, jump_eig', tmp)
                    b_plus_summand .+= precomp.b_plus[tau] .* diag_u .* M .* diag_u_row
                end
            end
            B_serial = zeros(CT, d, d)
            for t in keys(precomp.b_minus)
                @. diag_u = exp(1im * eigvals * (t / cfg.sigma))
                diag_u_row = transpose(diag_u)
                B_serial .+= precomp.b_minus[t] .* conj.(diag_u) .* b_plus_summand .* diag_u_row
            end
            t0o = QuantumFurnace.register_t0_b_minus(cfg)
            t0i = QuantumFurnace.register_t0_b_plus(cfg)
            B_serial .*= t0o * t0i

            err = norm(B_threaded .- B_serial)
            rel = err / max(norm(B_serial), 1.0)
            @test rel < 1e-12
            @info "qf-6af B_time threaded" err=err rel=rel
        else
            @test true
        end
    end

    # ------------------------------------------------------------
    # (g) B_trotter threaded
    # ------------------------------------------------------------
    @testset "(g) B_trotter threaded vs serial" begin
        if Threads.nthreads() > 1
            cfg = make_config(Lindbladian(), TrotterDomain(); construction=KMS())
            precomp = QuantumFurnace._precompute_data(cfg, TEST_TROTTER)
            t0o = QuantumFurnace.register_t0_b_minus(cfg)
            t0i = QuantumFurnace.register_t0_b_plus(cfg)
            B_threaded = QuantumFurnace.B_trotter(TEST_TROTTER_JUMPS, TEST_TROTTER,
                precomp.b_minus, precomp.b_plus, t0o, t0i, cfg.beta, cfg.sigma)

            d = DIM
            CT = ComplexF64
            eigvals_outer = TEST_TROTTER.eigvals_t0_b_minus !== nothing ?
                TEST_TROTTER.eigvals_t0_b_minus : TEST_TROTTER.eigvals_t0
            eigvals_inner = TEST_TROTTER.eigvals_t0_b_plus !== nothing ?
                TEST_TROTTER.eigvals_t0_b_plus : TEST_TROTTER.eigvals_t0
            t0_step_outer = TEST_TROTTER.t0_b_minus !== nothing ?
                TEST_TROTTER.t0_b_minus : TEST_TROTTER.t0
            t0_step_inner = TEST_TROTTER.t0_b_plus !== nothing ?
                TEST_TROTTER.t0_b_plus : TEST_TROTTER.t0

            b_plus_summand = zeros(CT, d, d)
            diag_u  = Vector{CT}(undef, d)
            diag_u2 = Vector{CT}(undef, d)
            tmp     = Matrix{CT}(undef, d, d)
            M       = Matrix{CT}(undef, d, d)
            for (tau, b_tau) in precomp.b_plus
                num_t0_steps = Int(round(tau * cfg.beta / t0_step_inner))
                @. diag_u  = eigvals_inner ^ num_t0_steps
                @. diag_u2 = eigvals_inner ^ (-2 * num_t0_steps)
                diag_u_row = transpose(diag_u)
                for jump_a in TEST_TROTTER_JUMPS
                    jump_a_eig = jump_a.in_eigenbasis
                    @. tmp = diag_u2 * jump_a_eig
                    mul!(M, jump_a_eig', tmp)
                    b_plus_summand .+= b_tau .* diag_u .* M .* diag_u_row
                end
            end
            B_serial = zeros(CT, d, d)
            for (t, b_t) in precomp.b_minus
                num_t0_steps = Int(round(t / (cfg.sigma * t0_step_outer)))
                @. diag_u = eigvals_outer ^ num_t0_steps
                diag_u_row = transpose(diag_u)
                B_serial .+= b_t .* conj.(diag_u) .* b_plus_summand .* diag_u_row
            end
            B_serial .*= t0o * t0i

            err = norm(B_threaded .- B_serial)
            rel = err / max(norm(B_serial), 1.0)
            @test rel < 1e-12
            @info "qf-6af B_trotter threaded" err=err rel=rel
        else
            @test true
        end
    end

    # ------------------------------------------------------------
    # (h) B_bohr threaded (jump × freq)
    # ------------------------------------------------------------
    @testset "(h) B_bohr threaded vs serial" begin
        if Threads.nthreads() > 1
            cfg = make_config(Lindbladian(), BohrDomain(); construction=KMS())
            B_threaded = QuantumFurnace.B_bohr(TEST_HAM, TEST_JUMPS, cfg)

            f = QuantumFurnace._pick_f(cfg)
            unique_freqs = collect(keys(TEST_HAM.bohr_dict))
            d = DIM
            B_serial = zeros(ComplexF64, d, d)
            f_A_nu_1 = zeros(ComplexF64, d, d)
            for jump in TEST_JUMPS
                for nu_2 in unique_freqs
                    indices = TEST_HAM.bohr_dict[nu_2]
                    @. f_A_nu_1 = f(TEST_HAM.bohr_freqs, nu_2) * jump.in_eigenbasis
                    @inbounds for idx in indices
                        i, j = idx[1], idx[2]
                        val = conj(jump.in_eigenbasis[i, j])
                        @inbounds for col in 1:d
                            B_serial[j, col] += val * f_A_nu_1[i, col]
                        end
                    end
                end
            end

            err = norm(B_threaded .- B_serial)
            rel = err / max(norm(B_serial), 1.0)
            @test rel < 1e-12
            @info "qf-6af B_bohr threaded" err=err rel=rel
        else
            @test true
        end
    end

    # ------------------------------------------------------------
    # (i) Channel jump-sandwich threaded — BohrDomain
    # ------------------------------------------------------------
    @testset "(i) Channel sandwich threaded vs serial: BohrDomain" begin
        if Threads.nthreads() > 1
            cfg = make_config(Thermalize(), BohrDomain(); construction=KMS())
            ws = Workspace(cfg, TEST_HAM, TEST_JUMPS)
            Random.seed!(456)
            rho = Matrix(random_density_matrix(NUM_QUBITS))

            sc = ws.scratch
            fill!(sc.channel_rho_jump, 0)
            QuantumFurnace._accumulate_jump_sandwich!(sc.channel_rho_jump,
                ws, rho, ws.delta, cfg, TEST_HAM)
            out_threaded = copy(sc.channel_rho_jump)

            d = DIM
            bohr_alpha_fn = ws.bohr_alpha
            gnf = ws.gamma_norm_factor
            bohr_keys = ws.bohr_keys === nothing ?
                collect(keys(TEST_HAM.bohr_dict)) : ws.bohr_keys
            out_serial = zeros(ComplexF64, d, d)
            for (k, eigenbasis) in enumerate(ws.jump_eigenbases)
                for nu_2 in bohr_keys
                    f_A = similar(eigenbasis)
                    @. f_A = bohr_alpha_fn(TEST_HAM.bohr_freqs, nu_2) * eigenbasis
                    rho_A = zeros(ComplexF64, d, d)
                    indices = TEST_HAM.bohr_dict[nu_2]
                    for idx in indices
                        i, j = idx[1], idx[2]
                        v = conj(eigenbasis[i, j])
                        for p in 1:d
                            rho_A[p, i] += rho[p, j] * v
                        end
                    end
                    out_serial .+= ws.delta * gnf * (f_A * rho_A)
                end
            end

            err = norm(out_threaded .- out_serial)
            rel = err / max(norm(out_serial), 1.0)
            @test rel < 1e-12
            @info "qf-6af channel sandwich bohr" err=err rel=rel
        else
            @test true
        end
    end

    # ------------------------------------------------------------
    # (j) BLAS thread restoration on every threaded entry
    # ------------------------------------------------------------
    @testset "(j) BLAS thread restoration after threaded entry" begin
        if Threads.nthreads() > 1
            old_blas = LinearAlgebra.BLAS.get_num_threads()

            cfg_e = make_config(Lindbladian(), EnergyDomain(); construction=KMS())
            ws_e = Workspace(cfg_e, TEST_HAM, TEST_JUMPS)
            @test LinearAlgebra.BLAS.get_num_threads() == old_blas

            cfg_t = make_config(Lindbladian(), TrotterDomain(); construction=KMS())
            ws_t = Workspace(cfg_t, TEST_HAM, TEST_TROTTER_JUMPS;
                trotter=TEST_TROTTER)
            @test LinearAlgebra.BLAS.get_num_threads() == old_blas

            cfg_b = make_config(Lindbladian(), BohrDomain(); construction=KMS())
            ws_b = Workspace(cfg_b, TEST_HAM, TEST_JUMPS)
            @test LinearAlgebra.BLAS.get_num_threads() == old_blas

            cfg_c = make_config(Thermalize(), BohrDomain(); construction=KMS())
            ws_c = Workspace(cfg_c, TEST_HAM, TEST_JUMPS)
            Random.seed!(0)
            rho = Matrix(random_density_matrix(NUM_QUBITS))
            apply_delta_channel!(ws_c, rho, cfg_c, TEST_HAM)
            @test LinearAlgebra.BLAS.get_num_threads() == old_blas

            @info "qf-6af BLAS restoration" old_blas=old_blas
        else
            @test true
        end
    end
end
