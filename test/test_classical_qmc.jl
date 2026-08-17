# Sandbox tests for the SSE quantum Monte Carlo implementation, used as a
# classical-computation baseline for quantum Gibbs sampling.
#
# Correctness gate = (1) the reconstructed physical Hamiltonian matches
# the fixture's un-rescaled Hamiltonian to machine precision, and (2) MC thermal observables
# reproduce the exact dense ρ_β within a few Monte-Carlo error bars.
#
# Non-flaky design: FIXED RNG seeds + enough sweeps that the MC error is small, then assert
# |MC − exact| < max(5·MC_err, abs_tol). The 5σ band (not 1σ) makes the test deterministic
# despite the ~1/150 single-σ outliers expected from honest statistics. Small n + modest
# sweeps keep it to a few seconds.

using QuantumFurnace
using Random
using Test

@testset "classical_qmc module (SSE quantum Monte Carlo)" begin

    @testset "validated SSE regime and statistics" begin
        @test_throws ArgumentError build_sse_heis_model(
            1; seed=1, periodic=false, disorder_strength=0.0)
        @test_throws ArgumentError build_sse_heis_model(
            4; seed=1, periodic=false, disorder_strength=0.0, J=0.0)
        @test_throws ArgumentError build_sse_heis_model(
            4; seed=1, periodic=false, disorder_strength=NaN)
        @test_throws ArgumentError build_sse_tfim_model(
            0, 2; seed=1, h=1.0, disorder_strength=0.0)
        @test_throws ArgumentError build_sse_tfim_model(
            2, 2; seed=1, h=0.0, disorder_strength=0.0)
        @test_throws ArgumentError build_sse_tfim_model(
            2, 2; seed=1, h=1.0, disorder_strength=-0.1)

        even_model = build_sse_heis_model(
            4; seed=1, periodic=true, disorder_strength=0.0)
        odd_model = build_sse_heis_model(
            3; seed=1, periodic=true, disorder_strength=0.0)
        @test_throws ArgumentError run_sse(
            odd_model, 0.5; pairs=[(1, 2)], nsweeps=100, nwarm=10)
        @test_throws ArgumentError run_sse(
            even_model, 0.0; pairs=[(1, 2)], nsweeps=100, nwarm=10)
        @test_throws ArgumentError run_sse(
            even_model, 0.5; pairs=[(1, 2)], nsweeps=1, nwarm=10)
        @test_throws ArgumentError run_sse(
            even_model, 0.5; pairs=[(1, 2)], nsweeps=100, nwarm=-1)
        @test_throws ArgumentError run_sse(
            even_model, 0.5; pairs=[(0, 2)], nsweeps=100, nwarm=10)
        @test_throws ArgumentError run_sse(
            even_model, 0.5; pairs=[(1, 1)], nsweeps=100, nwarm=10)
        @test_throws ArgumentError sse_exact_reference(
            even_model, Inf; pairs=[(1, 2)])
        @test_throws ArgumentError sse_exact_reference(
            even_model, 0.5; pairs=[(1, 5)])

        qmc = QuantumFurnace.ClassicalQMC
        @test qmc.count_signflips([1.0, 0.0, -1.0, 0.0, 1.0]) == 2
        @test qmc.count_signflips([0.0, 1.0, 0.0, 1.0, 0.0]) == 0

        rng = MersenneTwister(20260816)
        correlated = zeros(50_000)
        for i in 2:length(correlated)
            correlated[i] = 0.8 * correlated[i - 1] + randn(rng)
        end
        tau, tau_err, window = qmc.tau_int(correlated; max_pairs=5_000_000)
        @test 3.5 < tau < 5.5
        @test tau_err > 0
        @test window >= 6 * tau

        unresolved = collect(1.0:10_000.0)
        @test_logs (:warn, r"No self-consistent autocorrelation window") begin
            truncated_tau, _, truncated_window = qmc.tau_int(
                unresolved; max_pairs=10_000)
            @test truncated_tau >= 0.5
            @test truncated_window == 1
        end
        @test_logs (:warn, r"measurement run is too short") begin
            @test qmc.blocking_nblocks(100, 20.0) == 2
        end
    end

    @testset "seeded TFIM update traces remain exact" begin
        qmc = QuantumFurnace.ClassicalQMC

        function opstring_signature(ops)
            signature = UInt64(0xcbf29ce484222325)
            @inbounds for op in ops
                for value in (UInt64(op.kind), UInt64(op.idx), UInt64(op.offdiag))
                    signature = (signature ⊻ value) * UInt64(0x00000100000001b3)
                end
            end
            return signature
        end

        function trace_updates(update!)
            model = build_sse_tfim_model(
                2, 2; seed=17, h=1.0, disorder_strength=0.2, J=1.0,
                periodic_x=true, periodic_y=true)
            rng = MersenneTwister(314159)
            state = qmc.SSEState(
                model.n,
                rand(rng, (-1, 1), model.n),
                qmc.Op[qmc.IDOP for _ in 1:20],
                20,
                0,
            )
            menu = qmc.tfim_menu(model)
            trace = NamedTuple[]
            for _ in 1:6
                qmc.tfim_diagonal_update!(state, model, menu, 0.7, rng)
                counts = update!(state, model, rng)
                push!(trace, (
                    spins=Tuple(state.spins),
                    n_op=state.n_op,
                    counts=counts,
                    signature=opstring_signature(state.opstring),
                ))
            end
            return trace
        end

        expected_cluster = [
            (spins=(-1, -1, -1, -1), n_op=12, counts=(3, 1), signature=0x2fc7e6c087551901),
            (spins=(-1, -1, -1, -1), n_op=14, counts=(3, 1), signature=0x8c4593ec9e3a03e9),
            (spins=(-1, -1, -1, -1), n_op=14, counts=(2, 1), signature=0xf5c33c54df8af1af),
            (spins=(-1, -1, -1, -1), n_op=14, counts=(2, 0), signature=0x7f4e0fac874d009a),
            (spins=(-1, -1, -1, -1), n_op=12, counts=(2, 1), signature=0x8d443d5f87911d78),
            (spins=(1, 1, 1, 1), n_op=8, counts=(1, 1), signature=0x721238b319cd3d88),
        ]
        expected_local = [
            (spins=(-1, -1, -1, -1), n_op=12, counts=(5, 1), signature=0x2fc7e6c087551901),
            (spins=(-1, -1, -1, -1), n_op=14, counts=(6, 1), signature=0x7414c1b143af1510),
            (spins=(-1, -1, -1, -1), n_op=12, counts=(5, 0), signature=0x4ba94f165c487d24),
            (spins=(-1, -1, -1, -1), n_op=11, counts=(5, 1), signature=0x24db771a5be21be2),
            (spins=(-1, -1, -1, -1), n_op=9, counts=(5, 1), signature=0x0540f654b7ff8389),
            (spins=(-1, -1, -1, -1), n_op=10, counts=(4, 0), signature=0x13579ae0e16226d3),
        ]

        @test trace_updates(qmc.tfim_cluster_update!) == expected_cluster
        @test trace_updates(qmc.tfim_local_update!) == expected_local
    end

    # --- Reconstruction: H_phys matches the fixture un-rescaled Hamiltonian to ~1e-10 -----------
    @testset "model reconstruction vs fixture" begin
        for n in (3, 4), ds in (0.0, 0.1)
            raw = build_heis_1d(n, [1.0, 1.0, 1.0]; seed = 46, periodic = true, disorder_strength = ds)
            m   = build_sse_heis_model(n; seed = 46, periodic = true, disorder_strength = ds)
            @test sse_reconstruction_error(m, raw) ≤ 1e-10
        end
        for (Lx, Ly) in ((2, 2),), ds in (0.0, 1e-3), h in (3.5, 1.0)
            raw = build_tfim_2d(Lx, Ly; J = 1.0, h = h, seed = 46, disorder_strength = ds)
            m   = build_sse_tfim_model(Lx, Ly; seed = 46, h = h, disorder_strength = ds)
            @test sse_reconstruction_error(m, raw) ≤ 1e-10
        end
        for (px, py) in ((true, false), (false, true), (false, false))
            raw = build_tfim_2d(
                2, 3; J=1.0, h=1.0, seed=46, disorder_strength=1e-3,
                periodic_x=px, periodic_y=py)
            m = build_sse_tfim_model(
                2, 3; J=1.0, h=1.0, seed=46, disorder_strength=1e-3,
                periodic_x=px, periodic_y=py)
            @test m.periodic_x === px
            @test m.periodic_y === py
            @test sse_reconstruction_error(m, raw) <= 1e-10
        end
    end

    # Helper: assert an MC estimate sits within a generous, deterministic band of the exact value.
    function check(mc, err, exact; abs_tol = 5e-3, k = 5.0)
        @test isfinite(mc) && isfinite(err) && err ≥ 0
        @test abs(mc - exact) ≤ max(k * err, abs_tol)
    end

    # --- Heisenberg deterministic operator-loop: clean + disordered, even-n PBC ------------------
    @testset "Heisenberg loop vs exact ρ_β" begin
        pairs = [(1, 2), (1, 3)]
        for ds in (0.0, 0.1)
            m  = build_sse_heis_model(4; seed = 46, periodic = true, disorder_strength = ds)
            β  = 0.5
            ex = sse_exact_reference(m, β; pairs = pairs)
            res = run_sse(m, β; pairs = pairs, nsweeps = 120_000, nwarm = 20_000, seed = 909)
            check(res.energy, res.energy_err, ex.energy)
            check(res.mz2, res.mz2_err, ex.mz2)
            for pr in pairs
                (mc, err) = res.zz[pr]
                check(mc, err, ex.zz[pr])
            end
            @test res.avg_sign == 1.0            # even-n PBC is bipartite ⇒ sign-free
            @test 0.3 ≤ res.loop_accept ≤ 0.7    # deterministic loop flips ~½
            @test res.tau_E[1] ≥ 0.5             # τ_int ≥ ½ by definition
        end
    end

    # --- TFIM Swendsen–Wang cluster: ordered + disordered phase ----------------------------------
    @testset "TFIM cluster vs exact ρ_β" begin
        pairs = [(1, 2)]
        for (h, ds) in ((1.0, 0.0), (3.5, 1e-3))
            m  = build_sse_tfim_model(2, 2; seed = 46, h = h, disorder_strength = ds)
            β  = 1.0
            ex = sse_exact_reference(m, β; pairs = pairs)
            res = run_sse(m, β; pairs = pairs, nsweeps = 120_000, nwarm = 20_000, seed = 707)
            check(res.energy, res.energy_err, ex.energy)
            check(res.mz2, res.mz2_err, ex.mz2)
            (mc, err) = res.zz[(1, 2)]
            check(mc, err, ex.zz[(1, 2)])
            @test res.avg_sign == 1.0            # TFIM is stoquastic ⇒ sign-free
        end
    end

    # --- Local-update baseline: correct where it mixes; metastable in deep order -----------------
    @testset "TFIM local update — correctness + freezing" begin
        # Where it mixes (paramagnet / mild order) the local sampler must reproduce exact ρ_β.
        for (h, β) in ((3.5, 0.25), (1.0, 1.0))
            m  = build_sse_tfim_model(2, 2; seed = 46, h = h, disorder_strength = 1e-3)
            ex = sse_exact_reference(m, β; pairs = [(1, 2)])
            res = run_sse(m, β; pairs = [(1, 2)], nsweeps = 120_000, nwarm = 20_000,
                          seed = 505, update = :local)
            check(res.energy, res.energy_err, ex.energy)
            check(res.mz2, res.mz2_err, ex.mz2)
            h == 3.5 && @test res.n_signflips > 100   # paramagnet: tunnels freely (ergodic)
        end
        # On this finite run in the deep ordered phase, the local update remains in one
        # sector (⟨m_z⟩ ↛ 0, few/no flips), while the cluster update tunnels
        # (⟨m_z⟩ ≈ 0, many flips). Sector-symmetric ⟨m_z²⟩ is
        # blind to this — the freezing shows only in the sector-odd ⟨m_z⟩ / the flip count.
        m  = build_sse_tfim_model(2, 4; seed = 46, h = 1.0, disorder_strength = 1e-3)
        rc = run_sse(m, 2.0; pairs = [(1, 2)], nsweeps = 60_000, nwarm = 20_000, seed = 7, update = :cluster)
        rl = run_sse(m, 2.0; pairs = [(1, 2)], nsweeps = 60_000, nwarm = 20_000, seed = 7, update = :local)
        @test abs(rc.mz_mean) < 0.1              # cluster: exact ⟨m_z⟩ = 0 reproduced
        @test rc.n_signflips > 1000              # cluster: many tunnelling events
        @test abs(rl.mz_mean) > 0.5              # local: trapped near ±m_typ (≈0.97)
        @test rl.n_signflips < 100               # local: tunnelling suppressed (frozen)
        @test_throws ArgumentError run_sse(m, 2.0; pairs = [(1, 2)], nsweeps = 10, update = :bogus)
    end

    # --- Sign bookkeeping: frustrated bond only on odd-n PBC -------------------------------------
    @testset "frustrated-bond / sign sectors" begin
        m4pbc = build_sse_heis_model(4; seed = 46, periodic = true, disorder_strength = 0.0)
        m3obc = build_sse_heis_model(3; seed = 46, periodic = false, disorder_strength = 0.0)
        m3pbc = build_sse_heis_model(3; seed = 46, periodic = true, disorder_strength = 0.0)
        @test count(m4pbc.frustrated) == 0       # even ring bipartite
        @test count(m3obc.frustrated) == 0       # open chain (path) bipartite
        @test count(m3pbc.frustrated) == 1       # odd ring: one frustrated wrap bond

        # Odd-n OPEN chain is bipartite ⇒ sign-free ⇒ validates against exact.
        pairs = [(1, 2)]
        β  = 0.5
        ex = sse_exact_reference(m3obc, β; pairs = pairs)
        res = run_sse(m3obc, β; pairs = pairs, nsweeps = 120_000, nwarm = 20_000, seed = 313)
        @test res.avg_sign == 1.0
        check(res.energy, res.energy_err, ex.energy)
        (mc, err) = res.zz[(1, 2)]
        check(mc, err, ex.zz[(1, 2)])
    end
end
