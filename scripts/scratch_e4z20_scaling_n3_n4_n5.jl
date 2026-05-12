using QuantumFurnace
using LinearAlgebra
using Printf
using BenchmarkTools

const QF = QuantumFurnace
const OMEGA_RANGE = 2.5
const T_MINUS = 18.0
const T_PLUS = 12.0
const BETA = 10.0
const SIGMA = 1.0 / BETA

function jumps_in_basis(V, n)
    jumps = JumpOp[]
    for pauli in [[X], [Y], [Z]], site in 1:n
        op = Matrix(pad_term(pauli, n, site)) ./ sqrt(3 * n)
        push!(jumps, JumpOp(op, V' * op * V, op == transpose(op), op == op'))
    end
    return jumps
end

function make_cfg(n; r_D, r_bm, r_bp, M_D, M_bm, M_bp)
    w0_D = OMEGA_RANGE / 2^r_D
    t0_D = 2pi / (2^r_D * w0_D)
    t0_bm = 2 * T_MINUS / 2^r_bm
    w0_bm = 2pi / (2^r_bm * t0_bm)
    t0_bp = 2 * T_PLUS / 2^r_bp
    w0_bp = 2pi / (2^r_bp * t0_bp)
    Config(
        sim = Lindbladian(), domain = TrotterDomain(), construction = KMS(),
        num_qubits = n, with_linear_combination = true,
        beta = BETA, sigma = SIGMA, s = 0.25, a = BETA/30, eta = 1e-3,
        num_energy_bits_D = r_D, t0_D = t0_D, w0_D = w0_D,
        num_energy_bits_b_minus = r_bm, t0_b_minus = t0_bm, w0_b_minus = w0_bm,
        num_energy_bits_b_plus = r_bp,  t0_b_plus  = t0_bp,  w0_b_plus = w0_bp,
        num_trotter_steps_per_t0 = M_D,
        num_trotter_steps_per_t0_b_minus = M_bm,
        num_trotter_steps_per_t0_b_plus  = M_bp,
    )
end

function probe(n, ham, label; r_D, r_bm, r_bp, M_D, M_bm, M_bp)
    cfg = make_cfg(n; r_D, r_bm, r_bp, M_D, M_bm, M_bp)
    validate_config!(cfg)
    trotter = make_trotter_for_config(ham, cfg)
    jumps = jumps_in_basis(trotter.eigvecs, n)
    d = size(ham.data, 1)
    t_L = @elapsed L = construct_lindbladian(jumps, cfg, ham; trotter=trotter)
    sigma_beta = Hermitian(trotter.eigvecs' * ham.eigvecs * ham.gibbs *
                            ham.eigvecs' * trotter.eigvecs)
    residue = norm(L * vec(Matrix(sigma_beta)))
    triple_bytes = (sizeof(trotter.D.eigvecs) + sizeof(trotter.D.eigvals_t0) +
                    sizeof(trotter.D.bohr_freqs) +
                    sizeof(trotter.b_minus.eigvecs) + sizeof(trotter.b_minus.eigvals_t0) +
                    sizeof(trotter.b_minus.bohr_freqs) +
                    sizeof(trotter.b_plus.eigvecs) + sizeof(trotter.b_plus.eigvals_t0) +
                    sizeof(trotter.b_plus.bohr_freqs) +
                    sizeof(trotter.R_bm_in_D) + sizeof(trotter.R_bp_in_D) +
                    sizeof(trotter.R_bm_in_bp))
    @printf("    %-50s  residue=%.3e  cache=%4d KB  L=%6.2f s\n",
            label, residue, triple_bytes ÷ 1024, t_L)
    flush(stdout)
    return (cfg, trotter, jumps, residue, triple_bytes, t_L)
end

# ---------------------------------------------------------------------------
# Pass 1 — residue at each n
# ---------------------------------------------------------------------------
for n in (3, 4, 5)
    println("="^80)
    @printf("n = %d  (d = %d)\n", n, 2^n)
    println("="^80)
    flush(stdout)
    ham = QF._load_hamiltonian_bson(
        joinpath(pwd(), "hamiltonians", "heis_xxx_zzdisordered_periodic_n$n.bson"), BETA)
    probe(n, ham, "qf-7xt:      (r_D=7, r_bm=6, r_bp=14, M=64,64,1)";
          r_D=7, r_bm=6, r_bp=14, M_D=64, M_bm=64, M_bp=1)
    probe(n, ham, "+1 bit r_D:  (r_D=8, r_bm=6, r_bp=14, M=128,128,1)";
          r_D=8, r_bm=6, r_bp=14, M_D=128, M_bm=128, M_bp=1)
    probe(n, ham, "+2 bit r_D:  (r_D=9, r_bm=6, r_bp=14, M=128,128,1)";
          r_D=9, r_bm=6, r_bp=14, M_D=128, M_bm=128, M_bp=1)
    println()
    flush(stdout)
end

# ---------------------------------------------------------------------------
# Pass 2 — cost breakdown at the chosen recipe per n.
# ---------------------------------------------------------------------------
println("="^80)
println("Cost breakdown: rotation vs B_trotter total at recipe achieving ≤1e-6")
println("="^80)
flush(stdout)
@printf("  %-3s  %-5s  %-12s  %-12s  %-12s  %-12s  %-12s\n",
        "n", "d", "B_trotter", "rotation", "ratio", "triple MB", "NUFFT KB")
flush(stdout)
for n in (3, 4, 5)
    ham = QF._load_hamiltonian_bson(
        joinpath(pwd(), "hamiltonians", "heis_xxx_zzdisordered_periodic_n$n.bson"), BETA)
    r_D = n == 3 ? 7 : (n == 4 ? 8 : 9)
    cfg, triple, jumps, residue, triple_bytes, _ =
        probe(n, ham, "(r_D=$r_D, r_bm=6, r_bp=14, M=128,128,1)";
              r_D=r_D, r_bm=6, r_bp=14, M_D=128, M_bm=128, M_bp=1)

    pd = QF._precompute_data(cfg, triple)
    b_minus = pd.b_minus
    b_plus = pd.b_plus
    t0_outer = QF.register_t0_b_minus(cfg)
    t0_inner = QF.register_t0_b_plus(cfg)
    nufft_bytes = sizeof(pd.oft_nufft_prefactors.data)

    # Warm up
    _ = B_trotter(jumps, triple, b_minus, b_plus, t0_outer, t0_inner, BETA, SIGMA)
    bench_total = @belapsed B_trotter($jumps, $triple, $b_minus, $b_plus,
                                       $t0_outer, $t0_inner, $BETA, $SIGMA) samples=3

    R_bp_in_D = triple.R_bp_in_D
    bench_rot = @belapsed begin
        @inbounds for j in $jumps
            _ = Matrix($R_bp_in_D * j.in_eigenbasis * $R_bp_in_D')
        end
    end samples=5

    @printf("  %-3d  %-5d  %8.2f ms  %7.2f μs   %.4f%%        %5.2f      %d\n",
            n, 2^n, 1000 * bench_total, 1e6 * bench_rot,
            100 * bench_rot / bench_total, triple_bytes / (1024^2), nufft_bytes ÷ 1024)
    flush(stdout)
end
