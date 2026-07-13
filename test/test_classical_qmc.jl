# Sandbox tests for the classical SSE QMC sampler (src/classical_qmc.jl, beads qf-h23.3).
#
# The sampler is the sign-free classical Gibbs-sampling baseline competing with the KMS
# Lindbladian sampler. Correctness gate = (1) the reconstructed physical Hamiltonian matches
# the fixture's un-rescaled Hamiltonian to machine precision, and (2) MC thermal observables
# reproduce the exact dense ρ_β within a few Monte-Carlo error bars.
#
# Non-flaky design: FIXED RNG seeds + enough sweeps that the MC error is small, then assert
# |MC − exact| < max(5·MC_err, abs_tol). The 5σ band (not 1σ) makes the test deterministic
# despite the ~1/150 single-σ outliers expected from honest statistics. Small n + modest
# sweeps keep it to a few seconds.

using QuantumFurnace
using Test

@testset "classical_qmc (SSE QMC) — qf-h23.3" begin

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

    # --- Local-update baseline: correct where it mixes; freezes (ergodicity breaks) in deep order --
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
        # Deep ordered phase: the local update freezes into one sector (⟨m_z⟩ ↛ 0, few/no flips),
        # while the cluster stays ergodic (⟨m_z⟩ ≈ 0, many flips). Sector-symmetric ⟨m_z²⟩ is
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
