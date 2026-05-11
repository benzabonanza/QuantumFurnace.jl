#!/usr/bin/env python
"""Debug the GQSP Qiskit circuit by comparing the realised unitary directly."""
import sys
from pathlib import Path

import numpy as np
from qiskit.quantum_info import Operator

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from gqsp.circuit import build_block_encoding_one_anc, build_walk, gqsp_circuit
from gqsp.jacobi_anger import jacobi_anger_coeffs, evaluate_polynomial_on_circle
from gqsp.berntson_sunderhauf import complementary_polynomial_bs
from gqsp.motlagh_wiebe import gqsp_extract_angles


def all_w_circuit_matrix_compute(W, lam, thetas, phis):
    """Manually compute the all-W GQSP unitary as a numpy matrix following MW's factorization.

    MW factorization (BS Eq. 1.2):  M(z) = L_0 · ∏_{j=1}^N D(z) · R_j  with D = diag(z, 1).
    In matrix product convention (left = applied last in time), this translates to the
    circuit time order: R_N first, then ctrl-W, ..., R_1, then ctrl-W, then L_0 last.

    Open controls (apply W when QSP=|0⟩) realise the substitution z → W in D(z).
    """
    n = W.shape[0]
    In = np.eye(n, dtype=complex)

    c_th0, s_th0 = np.cos(thetas[0]), np.sin(thetas[0])
    L0 = np.array([
        [np.exp(1j * (lam + phis[0])) * c_th0, np.exp(1j * lam) * s_th0],
        [np.exp(1j * phis[0]) * s_th0, -c_th0],
    ], dtype=complex)

    P0 = np.array([[1, 0], [0, 0]], dtype=complex)
    P1 = np.array([[0, 0], [0, 1]], dtype=complex)
    ctrl_W = np.kron(P0, W) + np.kron(P1, In)

    N = len(thetas) - 1
    # Time order: R_N first, ctrl-W, R_{N-1}, ..., R_1, ctrl-W, L_0 last.
    # Matrix product (= product in reverse time): L_0 · ctrl-W · R_1 · ctrl-W · R_2 · ... · ctrl-W · R_N
    U = np.eye(2 * n, dtype=complex)
    for j in range(N, 0, -1):  # j = N, N-1, …, 1: R_j applied first (rightmost in matrix product order)
        c_th, s_th = np.cos(thetas[j]), np.sin(thetas[j])
        eφ = np.exp(1j * phis[j])
        Rj = np.array([
            [eφ * c_th, s_th],
            [eφ * s_th, -c_th],
        ], dtype=complex)
        Rj_full = np.kron(Rj, In)
        U = Rj_full @ U
        U = ctrl_W @ U
    # Finally L_0 (applied last in time = leftmost in matrix product)
    U = np.kron(L0, In) @ U
    return U


def main():
    np.random.seed(42)
    n_sys = 1
    d = 1
    delta = 0.1

    # Random Hermitian H
    rng = np.random.default_rng(seed=43)
    G = rng.standard_normal((2 ** n_sys, 2 ** n_sys)) + 1j * rng.standard_normal((2 ** n_sys, 2 ** n_sys))
    H = (G + G.conj().T) / 2
    H = 0.9 * H / np.linalg.norm(H, ord=2)
    print("H:", H)

    UH = build_block_encoding_one_anc(H)
    W = build_walk(UH)
    print(f"\nW (shape {W.shape}):")
    print(np.round(W, 4))

    # Get angles from Jacobi-Anger + BS + MW
    p = jacobi_anger_coeffs(d, delta)
    theta_grid = np.linspace(0, 2 * np.pi, 4096, endpoint=False)
    Pmax = float(np.max(np.abs(evaluate_polynomial_on_circle(p, theta_grid))))
    rescale = 0.99 / max(Pmax, 0.99)
    p_rs = rescale * p
    q = complementary_polynomial_bs(p_rs, max(256, 32 * (d + 1)))
    lam, thetas, phis = gqsp_extract_angles(p_rs, q)
    print(f"\nlam={lam:.4f}, thetas={thetas}, phis={phis}")
    print(f"rescale={rescale:.4f}")

    # Manually compute the all-W GQSP unitary
    U_manual = all_w_circuit_matrix_compute(W, lam, thetas, phis)
    print(f"\nManual U_manual (shape {U_manual.shape}):")
    n = W.shape[0]
    # Top-left QSP block: ⟨0_Q| U |0_Q⟩ on the target = U[:n, :n]
    P_W_manual = U_manual[:n, :n]
    print("P(W) from manual circuit (top-left QSP block):")
    print(np.round(P_W_manual, 4))

    # Direct P(W) from polynomial coefficients (rescaled)
    P_W_direct = sum(p_rs[m] * np.linalg.matrix_power(W, m) for m in range(len(p_rs)))
    print("\nP(W) from direct polynomial:")
    print(np.round(P_W_direct, 4))
    print(f"\n||P(W)_manual - P(W)_direct||₂ = {np.linalg.norm(P_W_manual - P_W_direct, ord=2):.3e}")

    # L_d(W) = P(W) * W^{-d}
    Wmd = np.linalg.matrix_power(W.conj().T, d)
    Ld_manual = P_W_manual @ Wmd
    print(f"\nL_d(W) from manual circuit (top-left QSP block @ W^-d):")
    print(np.round(Ld_manual / rescale, 4))

    Ld_direct = P_W_direct @ Wmd
    print(f"\nL_d(W) direct:")
    print(np.round(Ld_direct / rescale, 4))
    print(f"||L_d_manual - L_d_direct||₂ = {np.linalg.norm(Ld_manual - Ld_direct, ord=2):.3e}")

    # exp(-iδH) on system from L_d(W) projected to anc=|0⟩
    n_sys_dim = 2 ** n_sys
    # In W's basis (anc MSB), anc=|0⟩ block is the top n_sys_dim x n_sys_dim
    Ld_anc0_manual = Ld_manual[:n_sys_dim, :n_sys_dim] / rescale
    Ld_anc0_direct = Ld_direct[:n_sys_dim, :n_sys_dim] / rescale
    print(f"\nL_d anc=|0⟩ block (manual):")
    print(np.round(Ld_anc0_manual, 4))
    print(f"L_d anc=|0⟩ block (direct):")
    print(np.round(Ld_anc0_direct, 4))

    eigvals, eigvecs = np.linalg.eigh(H)
    target = eigvecs @ np.diag(np.exp(-1j * delta * eigvals)) @ eigvecs.conj().T
    print(f"exp(-iδH) target:")
    print(np.round(target, 4))
    print(f"||L_d_anc0_manual - target||₂ = {np.linalg.norm(Ld_anc0_manual - target, ord=2):.3e}")

    # Now compare with the Qiskit-built circuit
    qc = gqsp_circuit(UH, lam, thetas, phis, d, realise="L_d")
    U_qiskit = Operator(qc).data
    print(f"\nU_qiskit shape: {U_qiskit.shape}")

    # Extract ⟨0_qsp, 0_anc | U_qiskit | 0_qsp, 0_anc⟩ on sys
    # Qiskit basis: index = sys * 4 + anc * 2 + qsp (for n_sys=1)
    block_qiskit = np.zeros((n_sys_dim, n_sys_dim), dtype=complex)
    for i in range(n_sys_dim):
        for j in range(n_sys_dim):
            i_full = i * 4 + 0  # qsp=0, anc=0, sys=i
            j_full = j * 4 + 0
            block_qiskit[i, j] = U_qiskit[i_full, j_full]
    block_qiskit = block_qiskit / rescale
    print(f"Qiskit ⟨0_qsp 0_anc| U |0_qsp 0_anc⟩ on sys:")
    print(np.round(block_qiskit, 4))
    print(f"||block_qiskit - target||₂ = {np.linalg.norm(block_qiskit - target, ord=2):.3e}")
    print(f"||block_qiskit - L_d_anc0_manual||₂ = {np.linalg.norm(block_qiskit - Ld_anc0_manual, ord=2):.3e}")


if __name__ == "__main__":
    main()
