"""
qf-e4z.28 — Empirical test: does X+ZZ disorder remove the parity-trap so
that single-pass Krylov suffices for both the trajectory AND the
spectral gap of `predict_lindbladian_trajectory`?

Background (see [[predict-trajectory-two-pass-qf-e4z-27]]):
  `predict_lindbladian_trajectory` currently runs TWO Krylov passes —
    Pass 1: single-seed Arnoldi from `vec(rho_0)` (this is the
            `_krylov_spectral_decomposition` call). Its eigenvalues are
            returned in `traj.eigenvalues`; `traj.eigenvalues[2]` is
            the SINGLE-PASS gap.
    Pass 2: a separate `krylov_spectral_gap` call seeded with
            `_krylov_default_x0 = vec(I/d + 1e-10·H_GUE_traceless)`,
            using KrylovKit's thick-restart Krylov-Schur. Its
            `spectral_gap` is returned in `traj.spectral_gap` and is
            the TRUE Lindbladian gap on parity-symmetric fixtures.
  Pass 2 costs an extra ~50-100 matvecs that we'd like to drop on
  fixtures where the bug cannot manifest.

The bug is a symmetry trap: when L and rho_0 share an exact symmetry
P̂ (e.g. Z^⊗N on legacy 1D Heisenberg + Z/ZZ disorder fixtures where
[[Z],[Z,Z]] disorder commutes with P = Z^⊗N), Arnoldi from
`vec(rho_0=I/d)` stays inside the parity-even sector. The true gap
lives in the parity-odd sector and is missed. The X+ZZ fixtures
(qf-e4z.27) have X-disorder that anti-commutes with P = Z^⊗N, so
[P, H] ≠ 0 and there is no parity sector to be trapped in.

This script tests:
  (A) On X+ZZ at n=5, n=7, β_phys=2.5, seed=42: does single-pass with
      rho_0 = I/d give a gap that matches the 2-pass reference?
      Accept if rel_err < 1e-6.
  (B) Cross-check on legacy Z+ZZ at the same cells: does single-pass
      FAIL as expected (rel_err ≈ 0.17-0.25)?
  (C) Optional alternative initial state `rho_0 = |+⟩⟨+|^⊗N` —
      irrelevant on X+ZZ (single-pass already works) but tested for
      curiosity on legacy Z+ZZ: does the |+⟩^⊗N state break the
      parity trap there too?

Run:
    julia --project scripts/scratch_qf_e4z_28_single_pass.jl

Coherent term is ON throughout (default `include_coherent=true` in
all Lindbladian calls — physical KMS Lindbladian, as the parity bug
is itself a property of the full L).
"""

using LinearAlgebra
using Printf
using QuantumFurnace
using QuantumFurnace: X, Y, Z, _parse_hamiltonian_bson, HamHam, beta_alg,
                     _build_jump_set, _krylov_spectral_decomposition,
                     _krylov_default_x0

const FIXTURE_DIR = joinpath(@__DIR__, "output", "multiseed_fixtures")

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

"""
    build_cfg(ham, n, β_phys) -> Config

Identical to `scratch_qf_e4z_27_verification.jl`: canonical CKG smooth-
Metropolis Lindbladian/EnergyDomain Config (r_D = 7, s = 0.25).
"""
function build_cfg(ham, n, β_phys)
    β_alg = beta_alg(ham, β_phys)
    σ = 1.0 / β_alg
    H_norm = maximum(abs, ham.eigvals)
    omega_range = 2.0 * (H_norm + 8 * σ)
    r_D = 7
    w0_D = omega_range / 2.0^r_D
    t0_D = 2π / (2.0^r_D * w0_D)
    return Config(
        sim = Lindbladian(), domain = EnergyDomain(), construction = KMS(),
        num_qubits = n, with_linear_combination = true,
        beta = β_alg, beta_phys = β_phys, sigma = σ,
        a = 0.0, s = 0.25, gaussian_parameters = (nothing, nothing),
        num_energy_bits_D = r_D, w0_D = w0_D, t0_D = t0_D,
        num_trotter_steps_per_t0 = 10, filter = nothing,
    )
end

"""
    rho_plus_tensor(n) -> Matrix{ComplexF64}

Construct `|+⟩⟨+|^⊗n` where |+⟩ = (|0⟩ + |1⟩)/√2. Equals
`((I + X)/2)^⊗n` — one Hadamard layer applied to |0⟩^⊗N.

Trace 1, rank 1, no Lindbladian-spectrum knowledge. On any Lindbladian
whose only conserved generator is P = Z^⊗N, this state has support on
both P-even and P-odd sectors (since H|0⟩ = |+⟩ has equal amplitude on
|0⟩ (P-even when n even) and |1⟩ (P-odd when n even); for n odd the
single-site parity flips and the analysis is sector-by-sector).
"""
function rho_plus_tensor(n)
    plus = ComplexF64[0.5 0.5; 0.5 0.5]   # |+⟩⟨+|
    rho = plus
    for _ in 2:n
        rho = kron(rho, plus)
    end
    return rho
end

"""
    parity_commutator_norm(ham, n) -> Float64

‖[P, H]‖_∞ with P = Z^⊗N, H the un-rescaled physical Hamiltonian.
0 iff H has exact Z^⊗N symmetry.
"""
function parity_commutator_norm(ham, n)
    Z_pauli = ComplexF64[1.0 0.0; 0.0 -1.0]
    P = Z_pauli
    for _ in 2:n
        P = kron(P, Z_pauli)
    end
    H_phys = Matrix(ham.data .* ham.rescaling_factor .- ham.shift .* I(2^n))
    comm = P * H_phys - H_phys * P
    return opnorm(comm)
end

"""
    run_one(fixture_stem, n, β_phys, seed; rho_0_kind) -> NamedTuple

Run a single (fixture, n, β_phys, seed, rho_0_kind) cell. Returns a
NamedTuple containing:
  - `gap_single_pass`  : `abs(real(traj.eigenvalues[2]))` — from
                         `_krylov_spectral_decomposition` seeded with
                         `vec(rho_0)`. This is what we'd get from a
                         single-pass predictor.
  - `gap_reference`    : `traj.spectral_gap` — 2-pass result from
                         `krylov_spectral_gap` with `_krylov_default_x0`
                         + KrylovKit thick-restart. The reference.
  - `rel_err`          : |gap_single_pass - gap_reference| / gap_reference
  - `traj_matvecs`     : total matvecs used (Pass 1 + Pass 2).
"""
function run_one(fixture_stem, n, β_phys, seed;
                 rho_0_kind::Symbol = :identity,
                 krylovdim_traj::Int = 40)
    ham_path = joinpath(FIXTURE_DIR,
        @sprintf("%s_periodic_n%d_seed%d.bson", fixture_stem, n, seed))
    isfile(ham_path) || error("missing fixture: $ham_path")
    raw = _parse_hamiltonian_bson(ham_path)
    ham = HamHam(raw; beta_phys=β_phys)
    jumps = _build_jump_set(ham, n)
    cfg = build_cfg(ham, n, β_phys)

    d = 2^n
    rho_0 = if rho_0_kind === :identity
        Matrix{ComplexF64}(I(d) / d)
    elseif rho_0_kind === :plus_tensor
        rho_plus_tensor(n)
    elseif rho_0_kind === :gue_perturbed
        # rho_0 = I/d + 1e-10·H_GUE_traceless, the seed `_krylov_default_x0`
        # uses for `krylov_spectral_gap`. This is the smallest perturbation
        # that survives KrylovKit `tol=1e-10` orthogonalisation while keeping
        # I/d numerically dominant. NB: trace(GUE-traceless) = 0 so trace
        # remains 1 exactly.
        Matrix{ComplexF64}(reshape(_krylov_default_x0(d), d, d))
    else
        throw(ArgumentError("unknown rho_0_kind = $rho_0_kind"))
    end

    # Reference: stand-alone krylov_spectral_gap call. This is what
    # `predict_lindbladian_trajectory`'s Pass 2 does internally.
    gap_res = krylov_spectral_gap(cfg, ham, jumps;
                                  krylovdim=30, howmany=4, tol=1e-10)
    gap_reference = gap_res.spectral_gap

    # Single-pass result: drive `predict_lindbladian_trajectory` (which
    # uses `rho_0` as the Arnoldi seed for Pass 1) and read out its
    # `eigenvalues` field. `traj.eigenvalues[2]` is exactly the
    # single-pass gap (the same number Pass 1 would return if we
    # dropped Pass 2). `traj.spectral_gap` is the 2-pass patched value.
    t_grid = collect(range(0.0, 10.0 / gap_reference, length=11))
    traj = predict_lindbladian_trajectory(cfg, ham, jumps, rho_0, t_grid;
                                          krylovdim=krylovdim_traj, tol=1e-10)
    gap_single_pass = abs(real(traj.eigenvalues[2]))
    rel_err = abs(gap_single_pass - gap_reference) / gap_reference

    return (
        gap_single_pass = gap_single_pass,
        gap_reference   = gap_reference,
        gap_2pass       = traj.spectral_gap,
        rel_err         = rel_err,
        traj_matvecs    = traj.total_matvecs,
    )
end

# -------------------------------------------------------------------------
# Driver
# -------------------------------------------------------------------------

function main()
    println("=" ^ 78)
    println("qf-e4z.28 — single-pass Krylov on X+ZZ disorder")
    println("=" ^ 78)
    summary = String[]

    cells = [(5, 2.5, 42), (7, 2.5, 42)]

    # ---------- Sanity: parity commutator norms ----------
    println("\n--- ‖[P, H]‖ sanity check (P = Z^⊗N) ---")
    println("X+ZZ should have ‖[P,H]‖ > 0 (parity broken).")
    println("Z+ZZ should have ‖[P,H]‖ ≈ 0 (parity preserved).")
    @printf("%-6s %-30s %-30s\n", "n", "X+ZZ ‖[P,H]‖", "Z+ZZ ‖[P,H]‖")
    for (n, _, seed) in cells
        path_xz = joinpath(FIXTURE_DIR,
            @sprintf("heis_xxx_XZdisordered_periodic_n%d_seed%d.bson", n, seed))
        path_zz = joinpath(FIXTURE_DIR,
            @sprintf("heis_xxx_disordered_periodic_n%d_seed%d.bson", n, seed))
        ham_xz = HamHam(_parse_hamiltonian_bson(path_xz); beta_phys=1.0)
        ham_zz = HamHam(_parse_hamiltonian_bson(path_zz); beta_phys=1.0)
        c_xz = parity_commutator_norm(ham_xz, n)
        c_zz = parity_commutator_norm(ham_zz, n)
        @printf("n=%d   %-30.4e %-30.4e\n", n, c_xz, c_zz)
    end

    # ---------- Test A: X+ZZ + rho_0 = I/d (the main question) ----------
    println("\n" * "=" ^ 78)
    println("Test A: X+ZZ fixtures, rho_0 = I/d (current default)")
    println("Acceptance: rel_err < 1e-6 ⇒ single-pass is correct on X+ZZ.")
    println("=" ^ 78)
    @printf("%-6s %-12s %-12s %-10s %-8s\n",
            "cell", "gap_1pass", "gap_ref", "rel_err", "matvecs")
    println("-" ^ 60)
    test_A_pass = true
    for (n, β_phys, seed) in cells
        r = run_one("heis_xxx_XZdisordered", n, β_phys, seed;
                    rho_0_kind = :identity)
        ok = r.rel_err < 1e-6
        test_A_pass &= ok
        push!(summary,
            @sprintf("A [n=%d β=%.1f X+ZZ I/d]: rel_err=%.2e  %s",
                n, β_phys, r.rel_err, ok ? "PASS" : "FAIL"))
        @printf("n=%d   %.4e  %.4e  %.2e  %d %s\n",
            n, r.gap_single_pass, r.gap_reference, r.rel_err, r.traj_matvecs,
            ok ? "" : "<-- FAIL")
    end

    # ---------- Test B: Z+ZZ + rho_0 = I/d (the bug control) ----------
    println("\n" * "=" ^ 78)
    println("Test B: Z+ZZ legacy fixtures, rho_0 = I/d (control)")
    println("Expectation: rel_err >> 1e-6 (the parity-trap bug should manifest).")
    println("=" ^ 78)
    @printf("%-6s %-12s %-12s %-10s %-8s\n",
            "cell", "gap_1pass", "gap_ref", "rel_err", "matvecs")
    println("-" ^ 60)
    test_B_bug_present = true
    for (n, β_phys, seed) in cells
        r = run_one("heis_xxx_disordered", n, β_phys, seed;
                    rho_0_kind = :identity)
        bug_present = r.rel_err > 0.01    # qf-e4z.27 saw ~17-25%
        test_B_bug_present &= bug_present
        push!(summary,
            @sprintf("B [n=%d β=%.1f Z+ZZ I/d]: rel_err=%.2e  bug %s",
                n, β_phys, r.rel_err,
                bug_present ? "PRESENT (expected)" : "MISSING (anomaly!)"))
        @printf("n=%d   %.4e  %.4e  %.2e  %d %s\n",
            n, r.gap_single_pass, r.gap_reference, r.rel_err, r.traj_matvecs,
            bug_present ? "(parity trap as expected)" : "<-- bug missing?")
    end

    # ---------- Test C: Z+ZZ + rho_0 = |+⟩⟨+|^⊗N (optional) ----------
    println("\n" * "=" ^ 78)
    println("Test C: Z+ZZ legacy fixtures, rho_0 = |+⟩⟨+|^⊗N")
    println("Does a |+⟩^⊗N starting state break the parity trap on Z+ZZ?")
    println("|+⟩ = H|0⟩, single Hadamard layer — easy to prepare. If yes,")
    println("we could swap the default rho_0 instead of running Pass 2.")
    println("=" ^ 78)
    @printf("%-6s %-12s %-12s %-10s %-8s\n",
            "cell", "gap_1pass", "gap_ref", "rel_err", "matvecs")
    println("-" ^ 60)
    for (n, β_phys, seed) in cells
        r = run_one("heis_xxx_disordered", n, β_phys, seed;
                    rho_0_kind = :plus_tensor)
        ok = r.rel_err < 1e-6
        push!(summary,
            @sprintf("C [n=%d β=%.1f Z+ZZ |+⟩^⊗N]: rel_err=%.2e  %s",
                n, β_phys, r.rel_err, ok ? "trap broken" : "trap holds"))
        @printf("n=%d   %.4e  %.4e  %.2e  %d %s\n",
            n, r.gap_single_pass, r.gap_reference, r.rel_err, r.traj_matvecs,
            ok ? "(trap broken)" : "(trap holds)")
    end

    # ---------- Test D: X+ZZ + rho_0 = |+⟩⟨+|^⊗N (consistency) ----------
    println("\n" * "=" ^ 78)
    println("Test D: X+ZZ fixtures, rho_0 = |+⟩⟨+|^⊗N (consistency)")
    println("If Test A passed, this should also pass — no symmetry left to")
    println("trap any initial state.")
    println("=" ^ 78)
    @printf("%-6s %-12s %-12s %-10s %-8s\n",
            "cell", "gap_1pass", "gap_ref", "rel_err", "matvecs")
    println("-" ^ 60)
    for (n, β_phys, seed) in cells
        r = run_one("heis_xxx_XZdisordered", n, β_phys, seed;
                    rho_0_kind = :plus_tensor)
        ok = r.rel_err < 1e-6
        push!(summary,
            @sprintf("D [n=%d β=%.1f X+ZZ |+⟩^⊗N]: rel_err=%.2e  %s",
                n, β_phys, r.rel_err, ok ? "PASS" : "FAIL"))
        @printf("n=%d   %.4e  %.4e  %.2e  %d %s\n",
            n, r.gap_single_pass, r.gap_reference, r.rel_err, r.traj_matvecs,
            ok ? "" : "<-- FAIL")
    end

    # ---------- Test E: GUE-perturbed seed on BOTH fixtures ----------
    println("\n" * "=" ^ 78)
    println("Test E: rho_0 = I/d + 1e-10·H_GUE_traceless (Pass 2's seed)")
    println("If single-pass with THIS seed gives the true gap on X+ZZ AND")
    println("Z+ZZ, we could merge Pass 1 with this seed and drop Pass 2.")
    println("Cost vs current Pass 1 seed: t=0 reconstruction lost to ~1e-10")
    println("on `rho_0 = I/d` trajectories (unmeasurable for our use).")
    println("=" ^ 78)
    @printf("%-18s %-12s %-12s %-10s %-8s\n",
            "cell", "gap_1pass", "gap_ref", "rel_err", "matvecs")
    println("-" ^ 65)
    test_E_pass_xz = true
    test_E_pass_zz = true
    for fix_stem in ("heis_xxx_XZdisordered", "heis_xxx_disordered")
        for (n, β_phys, seed) in cells
            r = run_one(fix_stem, n, β_phys, seed; rho_0_kind = :gue_perturbed)
            ok = r.rel_err < 1e-6
            tag = fix_stem == "heis_xxx_XZdisordered" ? "X+ZZ" : "Z+ZZ"
            tag == "X+ZZ" && (test_E_pass_xz &= ok)
            tag == "Z+ZZ" && (test_E_pass_zz &= ok)
            push!(summary,
                @sprintf("E [n=%d β=%.1f %s I/d+GUE]: rel_err=%.2e  %s",
                    n, β_phys, tag, r.rel_err, ok ? "PASS" : "FAIL"))
            @printf("n=%d %s   %.4e  %.4e  %.2e  %d %s\n",
                n, tag, r.gap_single_pass, r.gap_reference, r.rel_err,
                r.traj_matvecs, ok ? "" : "<-- FAIL")
        end
    end

    # ---------- Summary ----------
    println("\n" * "=" ^ 78)
    println("Summary:")
    println("=" ^ 78)
    for line in summary
        println("  ", line)
    end
    println()
    if test_A_pass && test_B_bug_present
        println("✓ qf-e4z.28 main result: X+ZZ structurally fixes the parity trap.")
        println("  Single-pass from rho_0 = I/d matches the 2-pass reference to")
        println("  < 1e-6 on X+ZZ at n=5,7. qf-e4z.27 Pass 2 unnecessary on X+ZZ.")
    elseif !test_A_pass
        println("✗ qf-e4z.28 main result: X+ZZ structural parity-break is NOT")
        println("  enough — single-pass from vec(I/d) still misses the gap on")
        println("  X+ZZ at n=5,7 with the same rel_err as Z+ZZ. The X-disorder")
        println("  strength (0.1) is too weak to enter the Arnoldi subspace.")
        if test_E_pass_xz && test_E_pass_zz
            println("✓ But Test E PASSED on both fixtures: a 1e-10 GUE-perturbed")
            println("  seed alone suffices for single-pass on X+ZZ AND Z+ZZ.")
            println("  Recommendation: change `_krylov_spectral_decomposition`'s")
            println("  default seed for `predict_lindbladian_trajectory` to")
            println("  `vec(rho_0 + 1e-10·H_GUE_traceless)` and drop Pass 2.")
        else
            println("✗ Test E also FAILED — keep the qf-e4z.27 two-pass.")
        end
    end
end

main()
