#!/usr/bin/env julia
# ============================================================================
# Momentum Sector Diagnostic Script
# ============================================================================
#
# Diagnoses why n=6 periodic Heisenberg chain has zero overlap between all
# observables and the Lindbladian gap mode, while n=4 works fine.
#
# Hypothesis: The Lindbladian inherits translational symmetry from the
# periodic chain. Eigenmodes have definite crystal momentum k = 2*pi*m/n.
# All 5 observables (H, Mz, XX_avg, YY_avg, ZZ_avg) are translationally
# invariant (k=0 sector). If the gap mode (second eigenvalue) lives in a
# nonzero momentum sector, all k=0 observables get exactly zero overlap.
#
# Usage:
#   cd QuantumFurnace.jl && julia --project experiments/diagnose_gap_momentum.jl
# ============================================================================

using QuantumFurnace
using LinearAlgebra
using Printf

# ============================================================================
# Section 0: Constants and Setup (from validate_spectral_gap.jl)
# ============================================================================

const BETA = 10.0
const DELTA = 0.01

# Grid parameters matching test suite conventions
const NUM_ENERGY_BITS = 12
const W0 = 0.05
const T0 = 2pi / (2^NUM_ENERGY_BITS * W0)
const NUM_TROTTER_STEPS_PER_T0 = 10
const SIGMA = 1.0 / BETA  # 0.1

# ---------------------------------------------------------------------------
# Helpers (copied from validate_spectral_gap.jl)
# ---------------------------------------------------------------------------

function make_system(n, beta)
    ham = HamHam([[X,X],[Y,Y],[Z,Z]], [1.0,1.0,1.0], n, beta; periodic=true)
    dim = 2^n

    jump_paulis = [[X], [Y], [Z]]
    num_of_jumps = 3 * n
    jump_normalization = sqrt(num_of_jumps)
    V = ham.eigvecs

    jumps = JumpOp[]
    for pauli in jump_paulis
        for site in 1:n
            jump_op = Matrix(pad_term(pauli, n, site)) ./ jump_normalization
            jump_in_eigen = V' * jump_op * V
            orthogonal = (jump_op == transpose(jump_op))
            herm = (jump_op == jump_op')
            push!(jumps, JumpOp(jump_op, jump_in_eigen, orthogonal, herm))
        end
    end

    return ham, jumps, dim
end

function make_liouv_config(n)
    LiouvConfig(;
        num_qubits = n,
        with_coherent = false,
        with_linear_combination = true,
        domain = TimeDomain(),
        beta = BETA,
        sigma = SIGMA,
        a = BETA / 30.0,
        b = 0.4,
        num_energy_bits = NUM_ENERGY_BITS,
        w0 = W0,
        t0 = T0,
        num_trotter_steps_per_t0 = NUM_TROTTER_STEPS_PER_T0,
    )
end

# ============================================================================
# Section 1: Build translation operator in Hilbert space
# ============================================================================

"""
    build_translation_operator(n) -> Matrix{ComplexF64}

Build the cyclic translation operator T for n qubits.
T|q_1, q_2, ..., q_n> = |q_2, q_3, ..., q_n, q_1>

For each computational basis state j (0-indexed), the permuted index is:
    new_index = ((j << 1) | (j >> (n-1))) & ((1 << n) - 1)

Julia uses 1-based indexing, so we add 1 to both indices.
"""
function build_translation_operator(n::Int)
    dim = 2^n
    mask = (1 << n) - 1
    T = zeros(ComplexF64, dim, dim)
    for j in 0:(dim - 1)
        # Cyclic left shift of n-bit integer
        new_j = ((j << 1) | (j >> (n - 1))) & mask
        # T maps |j> to |new_j>: T[new_j, j] = 1
        T[new_j + 1, j + 1] = 1.0
    end
    return T
end

# ============================================================================
# Main analysis
# ============================================================================

println("=" ^ 70)
println("  MOMENTUM SECTOR DIAGNOSTIC")
println("  Hypothesis: gap mode lives in non-zero momentum sector")
println("  beta=$BETA")
println("=" ^ 70)
println()

for n in [4, 6]
    dim = 2^n
    d2 = dim^2  # Liouville space dimension

    println("=" ^ 70)
    @printf("  System: n=%d periodic Heisenberg chain (dim=%d, Liouville %dx%d)\n",
            n, dim, d2, d2)
    println("=" ^ 70)
    println()

    # --- Step 1: Build system and Lindbladian ---
    @printf("Building system and Lindbladian for n=%d...\n", n)
    ham, jumps, _ = make_system(n, BETA)
    config_l = make_liouv_config(n)
    liouv_result = run_lindbladian(jumps, config_l, ham)
    L = Matrix(liouv_result.liouvillian)
    exact_gap = abs(real(liouv_result.spectral_gap))
    @printf("Exact gap (ARPACK): %.10f\n\n", exact_gap)

    # --- Step 2: Build translation operator T in Hilbert space ---
    @printf("--- Step 2: Translation operator T (Hilbert space) ---\n")
    T_hilbert = build_translation_operator(n)

    # Verify T^n = I
    T_power = T_hilbert^n
    tn_error = norm(T_power - I(dim))
    @printf("||T^%d - I|| = %.3e  (%s)\n", n, tn_error,
            tn_error < 1e-12 ? "PASS" : "FAIL")

    # Verify T commutes with H_comp (computational basis Hamiltonian)
    H_comp = ham.data
    TH_comm = norm(T_hilbert * H_comp - H_comp * T_hilbert) / norm(H_comp)
    @printf("||[T, H]|| / ||H|| = %.3e  (%s)\n", TH_comm,
            TH_comm < 1e-12 ? "PASS: T commutes with H" : "FAIL")
    println()

    # --- Step 3: Build T_L in Liouville space (eigenbasis) ---
    @printf("--- Step 3: Translation operator T_L (Liouville/eigenbasis) ---\n")

    # Transform T to Hamiltonian eigenbasis: T_eigen = V' * T * V
    V = ham.eigvecs
    T_eigen = V' * T_hilbert * V

    # T_L in eigenbasis: T_L_eigen = kron(T_eigen, conj(T_eigen))
    # Watrous convention: vec(rho) column-stacking, T_L * vec(rho) = vec(T * rho * T')
    T_L = kron(T_eigen, conj(T_eigen))

    # Verify T_L commutes with L
    comm_norm = norm(T_L * L - L * T_L) / norm(L)
    @printf("||[T_L, L]|| / ||L|| = %.3e  (%s)\n", comm_norm,
            comm_norm < 1e-10 ? "PASS: T_L commutes with L" : "FAIL")
    println()

    # --- Step 4: Eigendecompose L and measure gap mode momentum ---
    @printf("--- Step 4: Eigendecomposition and momentum analysis ---\n")
    F = eigen(L)
    sorted_idx = sortperm(abs.(real.(F.values)))

    # Report first 10 eigenvalues and their momentum sectors
    n_report = min(10, length(sorted_idx))
    @printf("\nFirst %d eigenmodes (sorted by |Re(lambda)|):\n", n_report)
    @printf("%-4s  %16s  %16s  %14s  %8s  %6s\n",
            "Mode", "Re(lambda)", "Im(lambda)", "T_L eigenval", "|T_L ev|", "k/2pi*n")
    @printf("%-4s  %16s  %16s  %14s  %8s  %6s\n",
            "-"^4, "-"^16, "-"^16, "-"^14, "-"^8, "-"^6)

    for i in 1:n_report
        idx = sorted_idx[i]
        lambda = F.values[idx]
        v = F.vectors[:, idx]

        # Compute T_L eigenvalue on this eigenvector via Rayleigh quotient
        Tv = T_L * v
        exp_ik = dot(v, Tv) / dot(v, v)
        abs_exp_ik = abs(exp_ik)

        # Extract momentum
        k = angle(exp_ik)
        m = round(Int, k * n / (2 * pi))
        # Normalize m to [0, n-1]
        m = mod(m, n)

        @printf("%-4d  %16.10f  %16.10f  %7.4f%+7.4fi  %8.6f  m=%d\n",
                i, real(lambda), imag(lambda),
                real(exp_ik), imag(exp_ik), abs_exp_ik, m)
    end
    println()

    # --- Detailed gap mode analysis ---
    gap_idx = sorted_idx[2]
    v_gap = F.vectors[:, gap_idx]
    lambda_gap = F.values[gap_idx]

    Tv_gap = T_L * v_gap
    exp_ik_gap = dot(v_gap, Tv_gap) / dot(v_gap, v_gap)
    k_gap = angle(exp_ik_gap)
    m_gap = mod(round(Int, k_gap * n / (2 * pi)), n)

    @printf("Gap mode (mode 2) detailed analysis:\n")
    @printf("  Eigenvalue lambda_2 = %.10f %+.10fi\n", real(lambda_gap), imag(lambda_gap))
    @printf("  |Re(lambda_2)| = %.10f (exact gap)\n", abs(real(lambda_gap)))
    @printf("  T_L eigenvalue = %.6f %+.6fi\n", real(exp_ik_gap), imag(exp_ik_gap))
    @printf("  |T_L eigenvalue| = %.10f (should be ~1.0)\n", abs(exp_ik_gap))
    @printf("  Momentum k = %.6f rad\n", k_gap)
    @printf("  Momentum sector m = %d (out of n=%d)\n", m_gap, n)
    @printf("  Gap mode in k=0 sector: %s\n", m_gap == 0 ? "YES" : "NO")
    println()

    # --- Check for near-degenerate modes near the gap ---
    gap_eigenval = abs(real(lambda_gap))
    @printf("--- Near-degenerate modes (within 10%% of gap eigenvalue) ---\n")
    for i in 2:min(20, length(sorted_idx))
        idx_i = sorted_idx[i]
        lam_i = F.values[idx_i]
        if abs(abs(real(lam_i)) - gap_eigenval) / gap_eigenval < 0.10
            v_i = F.vectors[:, idx_i]
            Tv_i = T_L * v_i
            exp_ik_i = dot(v_i, Tv_i) / dot(v_i, v_i)
            m_i = mod(round(Int, angle(exp_ik_i) * n / (2 * pi)), n)
            @printf("  Mode %d: Re(lam)=%.10f, Im(lam)=%.10f, sector m=%d\n",
                    i, real(lam_i), imag(lam_i), m_i)
        end
    end
    println()

    # --- Step 5: Verify observable momentum sectors ---
    @printf("--- Step 5: Observable momentum sector analysis ---\n")
    obs, obs_names = build_preset_trajectory_observables(ham, n)

    @printf("%-12s  %14s  %8s  %6s  %s\n",
            "Observable", "T_L Rayleigh", "|value|", "Sector", "Is k=0?")
    @printf("%-12s  %14s  %8s  %6s  %s\n",
            "-"^12, "-"^14, "-"^8, "-"^6, "-"^7)

    all_k0 = true
    for (i, name) in enumerate(obs_names)
        O = obs[i]
        O_vec = vec(O)  # Vectorize observable (column-stacking)

        # Rayleigh quotient: <O_vec | T_L_eigen | O_vec> / <O_vec | O_vec>
        rq = dot(O_vec, T_L * O_vec) / dot(O_vec, O_vec)
        abs_rq = abs(rq)

        # Check if k=0 (eigenvalue should be ~1.0)
        is_k0 = abs(rq - 1.0) < 1e-8
        if !is_k0
            all_k0 = false
        end

        @printf("%-12s  %7.4f%+7.4fi  %8.6f  %6s  %s\n",
                name, real(rq), imag(rq), abs_rq,
                is_k0 ? "k=0" : "k!=0",
                is_k0 ? "YES" : "NO")
    end
    println()
    @printf("All observables in k=0 sector: %s\n\n", all_k0 ? "YES" : "NO")

    # --- Step 6: Summary and conclusion ---
    println("=" ^ 50)
    @printf("  SUMMARY for n=%d\n", n)
    println("=" ^ 50)
    @printf("  ||[T_L, L]|| / ||L||  = %.3e  (%s)\n", comm_norm,
            comm_norm < 1e-10 ? "translational symmetry confirmed" : "WARNING")
    @printf("  Gap eigenvalue        = %.10f\n", abs(real(lambda_gap)))
    @printf("  Gap mode T_L eigenval = %.4f%+.4fi\n", real(exp_ik_gap), imag(exp_ik_gap))
    @printf("  Gap mode momentum     = k = 2*pi*%d/%d\n", m_gap, n)
    @printf("  Gap mode in k=0       = %s\n", m_gap == 0 ? "YES" : "NO")
    @printf("  All obs in k=0        = %s\n", all_k0 ? "YES" : "NO")
    println()

    if m_gap != 0 && all_k0
        @printf("  CONFIRMED: gap mode is in k=%d sector (k=2*pi*%d/%d), not k=0\n",
                m_gap, m_gap, n)
        @printf("  -- explains zero overlap: k=0 observables have exactly zero\n")
        @printf("     projection onto k!=%d eigenmodes by orthogonality of\n", 0)
        @printf("     translational symmetry sectors.\n")
    elseif m_gap == 0
        @printf("  REFUTED: gap mode IS in k=0 sector -- momentum is NOT the explanation\n")
        @printf("  for zero overlap. Another mechanism must be responsible.\n")
    else
        @printf("  PARTIAL: gap mode in k=%d, but not all obs are k=0.\n", m_gap)
        @printf("  Situation is more complex than the simple hypothesis.\n")
    end
    println()
end

println("=" ^ 70)
println("  DIAGNOSTIC COMPLETE")
println("=" ^ 70)
