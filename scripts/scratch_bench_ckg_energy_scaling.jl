#!/usr/bin/env julia
#
# Scratch: CKG EnergyDomain scaling benchmark (qf-lkb.11.5).
#
# Strategy:
#   For n ∈ {3, 4, 5} (and optionally n=6 if a fixture is available),
#   β = 10, smooth Metropolis (a=0, s=0.25), measure:
#     (a) wall time per `apply_lindbladian!` matvec for EnergyDomain CKG,
#     (b) wall time of an end-to-end `integrate_to_gibbs` ODE solve over a
#         representative t_grid.
#
# Compare against the (b) pre-existing BohrDomain timing at n=3,4 from the
# qf-lkb.6 baseline (n=3 ≈ 4.7s, n=4 ≈ 100s).
#
# PHYSICS CHECK: a=0, s=0.25 (thesis convention).
# PHYSICS CHECK: krylovdim=20 matches the production sweep horizon. t_grid
# chosen to span ~5/gap_est on a small but non-trivial integration window.

using Printf
using LinearAlgebra
using Random
using QuantumFurnace
const QF = QuantumFurnace

const HAM_DIR = joinpath(@__DIR__, "..", "hamiltonians")

function build_jumps(ham::HamHam, n::Integer)
    paulis = ([X], [Y], [Z])
    num_jumps = length(paulis) * n
    jump_norm = sqrt(num_jumps)
    jumps = JumpOp[]
    for pauli in paulis
        for site in 1:n
            op = Matrix(pad_term(pauli, n, site)) ./ jump_norm
            op_eb = ham.eigvecs' * op * ham.eigvecs
            push!(jumps, JumpOp(op, op_eb,
                                op == transpose(op), op == op'))
        end
    end
    return jumps
end

function build_config(domain, n::Integer, beta::Real)
    Config(;
        sim = Lindbladian(),
        domain = domain,
        construction = KMS(),
        num_qubits = n,
        with_linear_combination = true,
        beta = beta,
        sigma = 1.0 / beta,
        a = 0.0,
        s = 0.25,
        num_energy_bits = 12,
        w0 = 0.05,
        t0 = 2π / (2^12 * 0.05),
        num_trotter_steps_per_t0 = 10,
    )
end

function bench_matvec(config, ham, jumps; n_warmup=3, n_samples=20)
    ws = Workspace(config, ham, jumps)
    d = size(ham.data, 1)
    rho = Matrix(QF.random_density_matrix(config.num_qubits))
    # warmup
    for _ in 1:n_warmup
        apply_lindbladian!(ws, rho, config, ham)
    end
    # measurement
    t0 = time()
    for _ in 1:n_samples
        apply_lindbladian!(ws, rho, config, ham)
    end
    elapsed = time() - t0
    return elapsed / n_samples
end

function bench_integrate(config, ham, jumps; t_max=2.0, t_steps=21, krylovdim=20, tol=1e-10)
    d = size(ham.data, 1)
    rho_0 = Matrix{ComplexF64}(I(d) / d)
    t_grid = collect(range(0.0, t_max, length=t_steps))
    t0 = time()
    res = integrate_to_gibbs(config, ham, jumps, rho_0, t_grid;
                              mode=:L, krylovdim=krylovdim, tol=tol)
    elapsed = time() - t0
    return elapsed, res.total_matvecs
end

println("="^80)
println("CKG EnergyDomain scaling benchmark (qf-lkb.11.5)")
println("smooth Metropolis a=0, s=0.25; β=10")
println("="^80)
@printf("Julia threads: %d, BLAS threads: %d\n",
        Threads.nthreads(), BLAS.get_num_threads())
println()
@printf("%-3s %-5s %-15s %-15s %-15s %-15s %-12s %-10s\n",
        "n", "dim", "matvec (Bohr)", "matvec (Energy)", "t_int (Bohr)",
        "t_int (Energy)", "matvecs", "speedup")
println("-"^120)

beta = 10.0
n_values = [3, 4, 5, 6, 7]

# Bohr timing inflates fast — limit to n≤4 to keep benchmark wall time reasonable
# At n=4, Bohr matvec is ~21x slower than Energy (40 ms vs 1.85 ms); n=5 would
# take ~5 minutes per measurement and is the entire pain point qf-lkb.11 fixes.
bohr_n_max = 4

# Pick the fixture family. Legacy heis_xxx_zzdisordered_periodic_n* covers n=3..5;
# the newer heis_xxx_zzdisordered_periodic_n* family extends to n=10 and
# eliminates the bipartite-pairing collisions for even n (see memory:
# Bohr Frequency Collision Root Cause). For the cross-n benchmark we want
# uniform physics, so use the newer family throughout.
function fixture_path(n::Integer)
    legacy = joinpath(HAM_DIR, "heis_xxx_zzdisordered_periodic_n$(n).bson")
    new_fam = joinpath(HAM_DIR, "heis_xxx_zzdisordered_periodic_n$(n).bson")
    return n <= 5 && isfile(legacy) ? legacy : new_fam
end

for n in n_values
    ham_path = fixture_path(n)
    if !isfile(ham_path)
        @printf("%-3d  fixture missing — skipping\n", n)
        continue
    end
    ham = QF._load_hamiltonian_bson(ham_path, beta)
    jumps = build_jumps(ham, n)
    d = size(ham.data, 1)

    # EnergyDomain
    cfg_e = build_config(EnergyDomain(), n, beta)
    t_mv_e = bench_matvec(cfg_e, ham, jumps)
    t_int_e, mvs_e = bench_integrate(cfg_e, ham, jumps;
                                       t_max = 5.0 / (1.0 / beta),  # ~5/gap_min heuristic
                                       t_steps = 21)

    # BohrDomain (only for small n)
    if n <= bohr_n_max
        cfg_b = build_config(BohrDomain(), n, beta)
        t_mv_b = bench_matvec(cfg_b, ham, jumps)
        t_int_b, mvs_b = bench_integrate(cfg_b, ham, jumps;
                                           t_max = 5.0 / (1.0 / beta),
                                           t_steps = 21)
        speedup_mv = t_mv_b / t_mv_e
        @printf("%-3d %-5d %-15.4e %-15.4e %-15.2f %-15.2f %-12d %-10.2fx\n",
                n, d, t_mv_b, t_mv_e, t_int_b, t_int_e, mvs_e, speedup_mv)
    else
        @printf("%-3d %-5d %-15s %-15.4e %-15s %-15.2f %-12d %-10s\n",
                n, d, "(skip)", t_mv_e, "(skip)", t_int_e, mvs_e, "—")
    end
    flush(stdout)
end

println("\nDone.")
