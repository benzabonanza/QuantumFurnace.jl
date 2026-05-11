#!/usr/bin/env python
"""Integration tests for the Python GQSP package.

Run:
    .venv-uv/bin/python src/python/tests/test_gqsp.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from qiskit.quantum_info import Operator

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from gqsp.berntson_sunderhauf import complementary_polynomial_bs
from gqsp.circuit import (
    _gqsp_circuit_bs_legacy,
    build_block_encoding_one_anc,
    build_walk,
    gqsp_circuit,
)
from gqsp.jacobi_anger import (
    evaluate_polynomial_on_circle,
    jacobi_anger_coeffs,
    laurent_target_cos,
    sup_norm_error,
)
from gqsp.motlagh_wiebe import (
    _gqsp_extract_angles_bs_legacy,
    _reconstruct_PQ_bs_legacy,
    gqsp_extract_angles,
    reconstruct_PQ,
)


def test_jacobi_anger_d1_closed_form():
    """L_1(e^{iθ}) = J_0(δα) - 2i J_1(δα) cos θ."""
    from scipy.special import jv

    delta_alpha = 0.7
    theta = np.linspace(0, 2 * np.pi, 257)[:-1]
    c = jacobi_anger_coeffs(1, delta_alpha)
    P = evaluate_polynomial_on_circle(c, theta)
    L1 = np.exp(-1j * theta) * P
    L1_closed = jv(0, delta_alpha) - 2j * jv(1, delta_alpha) * np.cos(theta)
    assert np.max(np.abs(L1 - L1_closed)) < 1e-13
    print("  ✓ d=1 closed form")


def test_jacobi_anger_quadratic_error():
    """At d=1, ‖L_1 - target‖_∞ = O((δα)²) for small δα."""
    for delta_alpha in (0.1, 0.05, 0.01):
        err, _, _ = sup_norm_error(1, delta_alpha)
        assert err < 0.5 * delta_alpha ** 2, f"d=1 quadratic error at δα={delta_alpha}"
    print("  ✓ d=1 quadratic error")


def test_bs_complementarity():
    """|P|² + |Q|² = 1 on the unit circle after BS Algorithm 1."""
    for delta_alpha in (0.1, 0.5, 1.0):
        for d in (1, 2, 3):
            p = jacobi_anger_coeffs(d, delta_alpha)
            theta = np.linspace(0, 2 * np.pi, 4096, endpoint=False)
            Pmax = float(np.max(np.abs(evaluate_polynomial_on_circle(p, theta))))
            rescale = 0.99 / max(Pmax, 0.99)
            p_rs = rescale * p
            N = max(256, 32 * (d + 1))
            while N & (N - 1) != 0:
                N += 1
            q = complementary_polynomial_bs(p_rs, N)
            P_vals = evaluate_polynomial_on_circle(p_rs, theta)
            Q_vals = evaluate_polynomial_on_circle(q, theta)
            err = float(np.max(np.abs(np.abs(P_vals) ** 2 + np.abs(Q_vals) ** 2 - 1)))
            assert err < 1e-10, f"BS complementarity at δα={delta_alpha}, d={d}: err={err}"
    print("  ✓ BS complementarity holds to ≤1e-10")


def test_mw_round_trip():
    """MW Algorithm 1 angles reconstruct (P, Q) via the matrix product."""
    test_cases = [(0.1, 1), (0.5, 1), (1.0, 1), (0.5, 2), (1.0, 2), (1.0, 3)]
    for delta_alpha, d in test_cases:
        p = jacobi_anger_coeffs(d, delta_alpha)
        theta = np.linspace(0, 2 * np.pi, 4096, endpoint=False)
        Pmax = float(np.max(np.abs(evaluate_polynomial_on_circle(p, theta))))
        rescale = 0.99 / max(Pmax, 0.99)
        p_rs = rescale * p
        N = max(256, 32 * (d + 1))
        while N & (N - 1) != 0:
            N += 1
        q = complementary_polynomial_bs(p_rs, N)
        lam, thetas, phis = gqsp_extract_angles(p_rs, q)
        P_rec, Q_rec = reconstruct_PQ(lam, thetas, phis)
        n = max(len(p_rs), len(P_rec))
        P_pad = np.concatenate([p_rs, np.zeros(n - len(p_rs), dtype=complex)])
        Prec_pad = np.concatenate([P_rec, np.zeros(n - len(P_rec), dtype=complex)])
        err = float(np.max(np.abs(P_pad - Prec_pad)))
        assert err < 1e-10, f"MW round-trip at δα={delta_alpha}, d={d}: err={err}"
    print("  ✓ MW round-trip ≤1e-10")


def test_block_encoding_unitarity():
    """1-anc qubitization U_H is unitary and has H in top-left block."""
    rng = np.random.default_rng(seed=42)
    for n_sys in (1, 2):
        dim = 2 ** n_sys
        G = rng.standard_normal((dim, dim)) + 1j * rng.standard_normal((dim, dim))
        H = (G + G.conj().T) / 2
        H = 0.9 * H / np.linalg.norm(H, ord=2)
        UH = build_block_encoding_one_anc(H)
        # Unitarity
        I_full = np.eye(2 * dim, dtype=complex)
        assert np.linalg.norm(UH @ UH.conj().T - I_full, ord=2) < 1e-12
        # Top-left block = H
        assert np.linalg.norm(UH[:dim, :dim] - H, ord=2) < 1e-12
    print("  ✓ block encoding U_H is unitary with H in top-left")


def test_walk_eigenstructure():
    """Walk W = (Z_anc) U_H has eigenvalues e^{±iθ} with cos θ = λ on qubitization subspaces."""
    rng = np.random.default_rng(seed=42)
    H = rng.standard_normal((2, 2)) + 1j * rng.standard_normal((2, 2))
    H = (H + H.conj().T) / 2
    H = 0.9 * H / np.linalg.norm(H, ord=2)
    eigvals, _ = np.linalg.eigh(H)
    UH = build_block_encoding_one_anc(H)
    W = build_walk(UH)
    Wevals = np.linalg.eigvals(W)
    # All eigenvalues of W on the unit circle:
    assert np.allclose(np.abs(Wevals), 1)
    # cos θ for each eigenvalue should match an H eigenvalue
    cos_thetas = np.real(Wevals)
    for lam in eigvals:
        assert any(abs(c - lam) < 1e-10 for c in cos_thetas), f"eigenvalue {lam} not found"
    print("  ✓ walk eigenstructure cos θ = λ")


def _build_qiskit_block(H, delta, d, *, realise):
    """Helper: build the GQSP circuit and extract the QSP=|0>, anc=|0> block."""
    UH = build_block_encoding_one_anc(H)
    p = jacobi_anger_coeffs(d, delta)
    theta = np.linspace(0, 2 * np.pi, 4096, endpoint=False)
    Pmax = float(np.max(np.abs(evaluate_polynomial_on_circle(p, theta))))
    rescale = 0.99 / max(Pmax, 0.99)
    p_rs = rescale * p
    N_bs = max(256, 32 * (d + 1))
    while N_bs & (N_bs - 1) != 0:
        N_bs += 1
    q = complementary_polynomial_bs(p_rs, N_bs)
    lam, thetas, phis = gqsp_extract_angles(p_rs, q)
    qc = gqsp_circuit(UH, lam, thetas, phis, d, realise=realise)
    U = Operator(qc).data
    dim = H.shape[0]
    block = np.zeros((dim, dim), dtype=complex)
    for i in range(dim):
        for j in range(dim):
            block[i, j] = U[i * 4, j * 4]
    return block / rescale


def _slope2_check(realise, label):
    rng = np.random.default_rng(seed=43)
    for n_sys in (1, 2):
        dim = 2 ** n_sys
        G = rng.standard_normal((dim, dim)) + 1j * rng.standard_normal((dim, dim))
        H = (G + G.conj().T) / 2
        H = 0.9 * H / np.linalg.norm(H, ord=2)
        eigvals, eigvecs = np.linalg.eigh(H)

        prev_err = None
        prev_delta = None
        for delta in (1e-1, 1e-2, 1e-3):
            d = 1
            block = _build_qiskit_block(H, delta, d, realise=realise)
            target = eigvecs @ np.diag(np.exp(-1j * delta * eigvals)) @ eigvecs.conj().T
            err = float(np.linalg.norm(block - target, ord=2))
            if prev_err is not None:
                slope = np.log(err / prev_err) / np.log(delta / prev_delta)
                assert abs(slope - 2.0) < 0.01, (
                    f"slope-2 violated for {label} at n_sys={n_sys}, δ={delta}: slope={slope}"
                )
            prev_err = err
            prev_delta = delta
    print(f"  ✓ Qiskit GQSP slope-2 ({label}, n_sys ∈ {{1, 2}}, δ ∈ [1e-3, 1e-1])")


def test_qiskit_slope2():
    """Form C (3d queries, MW Eq. 52): slope-2 in δ."""
    _slope2_check("L_d", "Form C: 2d ctrl-W + W^-d tail")


def test_qiskit_slope2_form_b():
    """Form B (2d queries, MW Eq. 46): slope-2 in δ."""
    _slope2_check("L_d_form_B", "Form B: d ctrl-W + d closed-ctrl-W†")


def _build_qiskit_block_legacy_bs(H, delta, d, *, realise):
    """BS-legacy variant of _build_qiskit_block using BS angles + BS circuit."""
    UH = build_block_encoding_one_anc(H)
    p = jacobi_anger_coeffs(d, delta)
    theta = np.linspace(0, 2 * np.pi, 4096, endpoint=False)
    Pmax = float(np.max(np.abs(evaluate_polynomial_on_circle(p, theta))))
    rescale = 0.99 / max(Pmax, 0.99)
    p_rs = rescale * p
    N_bs = max(256, 32 * (d + 1))
    while N_bs & (N_bs - 1) != 0:
        N_bs += 1
    q = complementary_polynomial_bs(p_rs, N_bs)
    lam, thetas, phis = _gqsp_extract_angles_bs_legacy(p_rs, q)
    qc = _gqsp_circuit_bs_legacy(UH, lam, thetas, phis, d, realise=realise)
    from qiskit.quantum_info import Operator
    U = Operator(qc).data
    dim = H.shape[0]
    block = np.zeros((dim, dim), dtype=complex)
    for i in range(dim):
        for j in range(dim):
            block[i, j] = U[i * 4, j * 4]
    return block / rescale


def test_mw_native_equals_bs_legacy_polynomial():
    """MW-native and BS-legacy circuits must realise the SAME polynomial P(W).

    The two forms differ in matrix-product convention (M_MW = M_BS^T, with
    R_MW = R_BS^T and L_0_MW = L_0_BS^T), gate matrices applied at each layer,
    and the time position of L_0 (first vs last). Despite these differences,
    the post-selected QSP=|0>, anc=|0> block must coincide because both forms
    realise the same Laurent polynomial L_d(W) on that block.
    """
    rng = np.random.default_rng(seed=45)
    for n_sys in (1, 2):
        dim = 2 ** n_sys
        G = rng.standard_normal((dim, dim)) + 1j * rng.standard_normal((dim, dim))
        H = (G + G.conj().T) / 2
        H = 0.9 * H / np.linalg.norm(H, ord=2)
        for delta in (1e-1, 1e-2):
            for d in (1, 2):
                for realise in ("L_d", "L_d_form_B"):
                    block_mw = _build_qiskit_block(H, delta, d, realise=realise)
                    block_bs = _build_qiskit_block_legacy_bs(H, delta, d, realise=realise)
                    err = float(np.linalg.norm(block_mw - block_bs, ord=2))
                    assert err < 1e-12, (
                        f"MW-native ≢ BS-legacy at n_sys={n_sys}, δ={delta}, d={d}, "
                        f"realise={realise}: err={err:.3e}"
                    )
    print("  ✓ MW-native ≡ BS-legacy realised polynomial (n_sys ∈ {1, 2}, δ ∈ [1e-2, 1e-1], "
          "d ∈ {1, 2}, both forms) ≤ 1e-12")


def test_form_b_equivalent_to_form_c():
    """Form B ≡ Form C: same MW+BS angles produce the same QSP=|0> block.

    Algebraic identity from MW2024 Eqs. 49--53. Equivalence is on the post-
    selected QSP=|0>, anc=|0> block (the two forms differ off the post-
    selected branch).
    """
    rng = np.random.default_rng(seed=44)
    for n_sys in (1, 2):
        dim = 2 ** n_sys
        G = rng.standard_normal((dim, dim)) + 1j * rng.standard_normal((dim, dim))
        H = (G + G.conj().T) / 2
        H = 0.9 * H / np.linalg.norm(H, ord=2)
        for delta in (1e-1, 1e-2):
            for d in (1, 2):
                block_C = _build_qiskit_block(H, delta, d, realise="L_d")
                block_B = _build_qiskit_block(H, delta, d, realise="L_d_form_B")
                err = float(np.linalg.norm(block_C - block_B, ord=2))
                assert err < 1e-12, (
                    f"Form B ≢ Form C at n_sys={n_sys}, δ={delta}, d={d}: err={err:.3e}"
                )
    print("  ✓ Form B ≡ Form C (n_sys ∈ {1, 2}, δ ∈ [1e-2, 1e-1], d ∈ {1, 2}) ≤ 1e-12")


def main():
    print("== Python GQSP integration tests ==")
    test_jacobi_anger_d1_closed_form()
    test_jacobi_anger_quadratic_error()
    test_bs_complementarity()
    test_mw_round_trip()
    test_block_encoding_unitarity()
    test_walk_eigenstructure()
    test_qiskit_slope2()
    test_qiskit_slope2_form_b()
    test_mw_native_equals_bs_legacy_polynomial()
    test_form_b_equivalent_to_form_c()
    print("\nAll tests passed.")


if __name__ == "__main__":
    main()
