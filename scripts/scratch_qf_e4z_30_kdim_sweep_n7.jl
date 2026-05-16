"""
qf-e4z.30 — Krylov-dim convergence of single-pass gap on
predict_lindbladian_trajectory at n=7 with rho_0 = |+⟩⟨+|^⊗N.

qf-e4z.28 measured rel_err 2.6e-5 (X+ZZ) / 9.2e-5 (Z+ZZ) at n=7
krylovdim=40 with the |+⟩⟨+|^⊗N seed — already 3-4 orders of
magnitude better than the parity trap from `vec(I/d)`, but above the
1e-6 threshold for thesis-grade gap reporting.

qf-e4z.29 verified at n=5 that |+⟩⟨+|^⊗N seed gives machine-precision
agreement with the dense reference. The n=7 residual is plain Krylov
truncation (40 eigenmodes captured out of d² = 16384), not symmetry.

This script measures whether increasing krylovdim brings the gap
agreement under 1e-7 at n=7.

For each (fixture, krylovdim):
  • Compute `gap_ref = krylov_spectral_gap(...).spectral_gap` ONCE per
    fixture (independent of krylovdim) — KrylovKit thick-restart
    Krylov-Schur, the reference.
  • Run `_krylov_spectral_decomposition` directly (= Pass 1 only)
    with `vec(|+⟩⟨+|^⊗N)` as the seed, at the swept krylovdim.
  • Compare `abs(real(decomp.eigenvalues[2]))` to `gap_ref`.

Coherent term ON throughout. Same Config as qf-e4z.28: r_D = 7,
s = 0.25, β_phys = 2.5, seed = 42. n = 7.

Run:
    julia --project scripts/scratch_qf_e4z_30_kdim_sweep_n7.jl
"""

using LinearAlgebra
using Printf
using QuantumFurnace
using QuantumFurnace: X, Y, Z, _parse_hamiltonian_bson, HamHam, beta_alg,
                     _build_jump_set, _krylov_spectral_decomposition,
                     Workspace, apply_lindbladian!

const FIXTURE_DIR = joinpath(@__DIR__, "output", "multiseed_fixtures")

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

function rho_plus_tensor(n)
    plus = ComplexF64[0.5 0.5; 0.5 0.5]
    rho = plus
    for _ in 2:n
        rho = kron(rho, plus)
    end
    return rho
end

"""
    single_pass_gap(ws, cfg, ham, rho_0, krylovdim) -> Tuple{Float64, Int}

Pass 1 only: run `_krylov_spectral_decomposition` from `vec(rho_0)`
and return (|Re(eigenvalues[2])|, matvec_count).
"""
function single_pass_gap(ws, cfg, ham, rho_0, krylovdim)
    d = size(rho_0, 1)
    fwd! = let ws = ws, cfg = cfg, ham = ham
        (out::AbstractMatrix, x::AbstractMatrix) -> begin
            apply_lindbladian!(ws, x, cfg, ham)
            copyto!(out, ws.scratch.rho_out)
            return out
        end
    end
    decomp = _krylov_spectral_decomposition(fwd!, rho_0, d;
                                            krylovdim=krylovdim, tol=1e-10)
    return abs(real(decomp.eigenvalues[2])), decomp.matvec_count
end

function main()
    n = 7
    β_phys = 2.5
    seed = 42
    kdims = (40, 60, 80, 100, 120)
    fixtures = ("heis_xxx_XZdisordered" => "X+ZZ",
                "heis_xxx_disordered" => "Z+ZZ")

    println("=" ^ 78)
    @printf("qf-e4z.30 — krylovdim sweep at n=%d, β_phys=%.2f, seed=%d, rho_0 = |+⟩⟨+|^⊗N\n",
            n, β_phys, seed)
    println("=" ^ 78)
    println("Acceptance threshold: rel_err < 1e-7")
    println()

    rho_0 = rho_plus_tensor(n)

    for (stem, label) in fixtures
        ham_path = joinpath(FIXTURE_DIR,
            @sprintf("%s_periodic_n%d_seed%d.bson", stem, n, seed))
        raw = _parse_hamiltonian_bson(ham_path)
        ham = HamHam(raw; beta_phys=β_phys)
        jumps = _build_jump_set(ham, n)
        cfg = build_cfg(ham, n, β_phys)

        # Build workspace ONCE per fixture — reused across all krylovdims.
        ws = Workspace(cfg, ham, jumps)

        # Reference gap: KrylovKit thick-restart Krylov-Schur, computed
        # once per fixture (it's krylovdim-independent for our purposes;
        # uses kdim=30 internally which is enough for KrylovKit to
        # converge under restart).
        gap_res = krylov_spectral_gap(cfg, ham, jumps;
                                       krylovdim=30, howmany=4, tol=1e-10)
        gap_ref = gap_res.spectral_gap

        @printf("\n--- Fixture: %s ---\n", label)
        @printf("gap_ref (krylov_spectral_gap, kdim=30 + restart): %.10e (%d matvecs)\n",
                gap_ref, gap_res.matvec_count)
        @printf("%-8s %-18s %-14s %-10s %-8s\n",
                "kdim", "gap_1pass", "rel_err", "verdict", "matvecs")
        println("-" ^ 65)

        for kdim in kdims
            gap_1, mv = single_pass_gap(ws, cfg, ham, rho_0, kdim)
            rel = abs(gap_1 - gap_ref) / gap_ref
            ok = rel < 1e-7
            verdict = ok ? "PASS" : "FAIL"
            @printf("%-8d %-18.10e %-14.4e %-10s %-8d\n",
                    kdim, gap_1, rel, verdict, mv)
        end
    end

    println("\n" * "=" ^ 78)
    println("Note: 'matvecs' is Pass 1 only (= krylovdim). The reference's")
    println("matvecs cost is amortised — one call per fixture, reused.")
    println()
    println("For comparison, qf-e4z.28 reported rel_err at kdim=40 (X+ZZ) = 2.59e-5,")
    println("(Z+ZZ) = 9.17e-5. We expect those numbers to reproduce here at kdim=40")
    println("and shrink with kdim.")
end

main()
