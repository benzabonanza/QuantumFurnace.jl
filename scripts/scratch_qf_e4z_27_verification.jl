"""
qf-e4z.27 verification (Part D of the issue plan).

Runs the targeted spot cells specified in `bd show qf-e4z.27` to confirm:

1. Bug regression — old [[Z],[Z,Z]] fixtures at (n=5, β_phys=2.5, seed=42)
   and (n=7, β_phys=2.5, seed=42): post-refactor `predict_lindbladian_trajectory`
   gap matches `krylov_spectral_gap` reference (= true Lindbladian gap).
2. No-regression on parity-trivial cells: n=3, 4, 6 at β_phys=2.5 on the
   old fixtures — post-refactor gap matches the pre-refactor behaviour to
   high tolerance (no degradation on cases without symmetry traps).
3. New X-disorder fixtures — bug doesn't manifest: regenerate
   (n=5, n=7, β_phys=2.5, seed=42) with [[X],[Z,Z]] disorder; verify that
   `predict_lindbladian_trajectory` and `krylov_spectral_gap` AGREE (no
   parity-trap to begin with).
4. Cross-check Hamiltonian symmetry: verify `[P, H] ≠ 0` for the new
   [[X],[Z,Z]] fixtures at all n ∈ {3..8} (P = Z^⊗N).

Run:
    julia --project scripts/scratch_qf_e4z_27_verification.jl

Output (stdout):
    Per-cell rel_err table for all four parts, plus a PASS/FAIL summary.
"""

using LinearAlgebra
using Printf
using BSON
using QuantumFurnace
using QuantumFurnace: X, Y, Z, _parse_hamiltonian_bson, HamHam, beta_alg,
                     _build_jump_set, pad_term, _construct_disordering_terms

const FIXTURE_DIR = joinpath(@__DIR__, "output", "multiseed_fixtures")

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

"""
    build_cfg(ham, n, β_phys) -> Config

Standard CKG smooth-Metropolis Lindbladian/EnergyDomain Config for the
qf-e4z.27 spot cells. r_D = 7, s = 0.25 (canonical per `.claude/CLAUDE.md`).
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
    run_cell(fixture_name, n, β_phys, seed) -> (gap_traj, gap_ref, rel_err, matvecs)

Load the named fixture, build a CKG smooth-Metropolis Config, and compare
`predict_lindbladian_trajectory`'s `spectral_gap` against
`krylov_spectral_gap` (the reference).
"""
function run_cell(fixture_name, n, β_phys, seed; krylovdim_traj=40)
    ham_path = joinpath(FIXTURE_DIR,
        @sprintf("%s_periodic_n%d_seed%d.bson", fixture_name, n, seed))
    isfile(ham_path) || error("missing fixture: $ham_path")
    raw = _parse_hamiltonian_bson(ham_path)
    ham = HamHam(raw; beta_phys=β_phys)
    jumps = _build_jump_set(ham, n)
    cfg = build_cfg(ham, n, β_phys)

    gap_res = krylov_spectral_gap(cfg, ham, jumps;
                                   krylovdim=30, howmany=4, tol=1e-10)
    d = 2^n
    rho_0 = Matrix{ComplexF64}(I(d) / d)
    t_grid = collect(range(0.0, 10.0 / gap_res.spectral_gap, length=11))
    traj = predict_lindbladian_trajectory(cfg, ham, jumps, rho_0, t_grid;
                                            krylovdim=krylovdim_traj, tol=1e-10)
    rel_err = abs(traj.spectral_gap - gap_res.spectral_gap) /
              gap_res.spectral_gap
    return (gap_traj = traj.spectral_gap, gap_ref = gap_res.spectral_gap,
            rel_err = rel_err, matvecs = traj.total_matvecs)
end

"""
    parity_commutator_norm(ham, n) -> Float64

Compute `‖[P, H]‖_∞` where P = Z^⊗N is the global Z-parity operator and
H is the (unrescaled) base Hamiltonian. Returns 0 iff the Hamiltonian
has exact Z^⊗N symmetry.
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

# -------------------------------------------------------------------------
# Driver
# -------------------------------------------------------------------------

function main()
    println("=" ^ 78)
    println("qf-e4z.27 verification (Part D of the issue plan)")
    println("=" ^ 78)
    summary = String[]

    # ---------- Part 1: Bug-regression cells (old Z+ZZ fixtures) ----------
    println("\n--- Part 1: Bug regression on old [[Z],[Z,Z]] fixtures ---")
    println("Expect: predict_lindbladian_trajectory.spectral_gap matches")
    println("        krylov_spectral_gap reference (which is qf-8fr-fixed).")
    @printf("%-6s %-12s %-12s %-10s %-8s\n", "cell", "gap_traj", "gap_ref", "rel_err", "matvecs")
    println("-" ^ 60)
    for (n, β_phys) in ((5, 2.5), (7, 2.5))
        r = run_cell("heis_xxx_disordered", n, β_phys, 42)
        ok = r.rel_err < 1e-8
        push!(summary,
            @sprintf("Part 1 [n=%d β_phys=%.1f Z+ZZ]: rel_err=%.2e  %s",
                n, β_phys, r.rel_err, ok ? "PASS" : "FAIL"))
        @printf("n=%d   %.4e  %.4e  %.2e  %d %s\n",
            n, r.gap_traj, r.gap_ref, r.rel_err, r.matvecs,
            ok ? "" : "<-- FAIL")
    end

    # ---------- Part 2: No-regression on parity-trivial cells ----------
    println("\n--- Part 2: No-regression on parity-trivial cells (old Z+ZZ) ---")
    println("These cells (n=3,4,6) at β_phys=2.5 may still have parity")
    println("symmetry, but small enough that pre-fix worked. New code should")
    println("also work, with high agreement to the krylov_spectral_gap reference.")
    @printf("%-6s %-12s %-12s %-10s %-8s\n", "cell", "gap_traj", "gap_ref", "rel_err", "matvecs")
    println("-" ^ 60)
    for n in (3, 4, 6)
        r = run_cell("heis_xxx_disordered", n, 2.5, 42)
        ok = r.rel_err < 1e-6
        push!(summary,
            @sprintf("Part 2 [n=%d β_phys=2.5 Z+ZZ]: rel_err=%.2e  %s",
                n, r.rel_err, ok ? "PASS" : "FAIL"))
        @printf("n=%d   %.4e  %.4e  %.2e  %d %s\n",
            n, r.gap_traj, r.gap_ref, r.rel_err, r.matvecs,
            ok ? "" : "<-- FAIL")
    end

    # ---------- Part 3: New X-disorder fixtures ----------
    println("\n--- Part 3: New [[X],[Z,Z]] disorder fixtures (no parity trap) ---")
    println("Z^⊗N parity is broken at the Hamiltonian level; the trajectory")
    println("predictor and krylov_spectral_gap should AGREE without any")
    println("seed-perturbation trick (which has nothing to fight against).")
    @printf("%-6s %-12s %-12s %-10s %-8s\n", "cell", "gap_traj", "gap_ref", "rel_err", "matvecs")
    println("-" ^ 60)
    for n in (5, 7)
        r = run_cell("heis_xxx_XZdisordered", n, 2.5, 42)
        ok = r.rel_err < 1e-6
        push!(summary,
            @sprintf("Part 3 [n=%d β_phys=2.5 X+ZZ]: rel_err=%.2e  %s",
                n, r.rel_err, ok ? "PASS" : "FAIL"))
        @printf("n=%d   %.4e  %.4e  %.2e  %d %s\n",
            n, r.gap_traj, r.gap_ref, r.rel_err, r.matvecs,
            ok ? "" : "<-- FAIL")
    end

    # ---------- Part 4: Symmetry cross-check ----------
    println("\n--- Part 4: [P, H] for new X+ZZ disorder (P = Z^⊗N) ---")
    println("Expect: ‖[P, H]‖ > 0 — i.e., X-disorder genuinely breaks Z^⊗N.")
    println("Compare to old Z+ZZ fixtures where [P, H] = 0 exactly.")
    @printf("%-6s %-30s %-30s\n", "n", "‖[P, H]‖ (X+ZZ)", "‖[P, H]‖ (Z+ZZ legacy)")
    println("-" ^ 60)
    for n in 3:8
        path_xz = joinpath(FIXTURE_DIR,
            @sprintf("heis_xxx_XZdisordered_periodic_n%d_seed%d.bson", n, 42))
        path_zz = joinpath(FIXTURE_DIR,
            @sprintf("heis_xxx_disordered_periodic_n%d_seed%d.bson", n, 42))
        # Skip n values not generated yet for X+ZZ
        if !isfile(path_xz)
            @printf("n=%d   [skip — fixture not generated]\n", n)
            continue
        end
        if !isfile(path_zz)
            @printf("n=%d   ‖[P, H]‖ (X+ZZ) only [missing Z+ZZ]\n", n)
            continue
        end
        raw_xz = _parse_hamiltonian_bson(path_xz)
        raw_zz = _parse_hamiltonian_bson(path_zz)
        ham_xz = HamHam(raw_xz; beta_phys=1.0)
        ham_zz = HamHam(raw_zz; beta_phys=1.0)
        comm_xz = parity_commutator_norm(ham_xz, n)
        comm_zz = parity_commutator_norm(ham_zz, n)
        ok = comm_xz > 1e-10
        push!(summary,
            @sprintf("Part 4 [n=%d  ‖[P,H]‖_X+ZZ=%.2e  ‖[P,H]‖_Z+ZZ=%.2e]: %s",
                n, comm_xz, comm_zz, ok ? "PASS" : "FAIL"))
        @printf("n=%d   %-30.4e %-30.4e %s\n",
            n, comm_xz, comm_zz,
            ok ? "" : "<-- FAIL: X disorder did NOT break parity")
    end

    # ---------- Summary ----------
    println("\n" * "=" ^ 78)
    println("Summary:")
    println("=" ^ 78)
    for line in summary
        println("  ", line)
    end
    n_pass = count(s -> contains(s, "PASS"), summary)
    n_fail = count(s -> contains(s, "FAIL"), summary)
    @printf("\nTotal: %d PASS, %d FAIL\n", n_pass, n_fail)
    if n_fail == 0
        println("\n✓ All qf-e4z.27 spot cells PASS.")
    else
        println("\n✗ Some qf-e4z.27 spot cells FAIL — investigate.")
    end
end

main()
