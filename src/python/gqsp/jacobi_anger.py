"""Jacobi-Anger polynomial coefficients for GQSP.

Convention (from 1_preliminaries.tex):
    target  f(theta) = exp(-i delta alpha cos theta)
    Laurent L_d(z)   = sum_{k=-d}^{d} (-i)^k J_k(delta alpha) z^k
    ordinary P(z)    = z^d L_d(z) = sum_{m=0}^{2d} c_m z^m
                       c_m = (-i)^(m-d) J_{m-d}(delta alpha)
"""
from __future__ import annotations

import numpy as np
from scipy.special import jv

D_MAX = 100


def jacobi_anger_coeffs(d: int, delta_alpha: float) -> np.ndarray:
    """Return c[0..2d] with c[m] = (-i)^(m-d) J_{m-d}(delta*alpha).

    The ordinary polynomial P(z) = sum_m c[m] z^m has degree 2d and satisfies
    P(z) = z^d L_d(z) where L_d(e^{i theta}) approximates exp(-i delta alpha cos theta).
    """
    if not (0 <= d <= D_MAX):
        raise ValueError(f"d must be in [0, {D_MAX}], got {d}")
    c = np.zeros(2 * d + 1, dtype=complex)
    for m in range(2 * d + 1):
        k = m - d
        # (-i)^k = exp(-i pi k / 2)
        c[m] = np.exp(-1j * np.pi * k / 2) * jv(k, delta_alpha)
    return c


def laurent_target_cos(delta_alpha: float, theta: np.ndarray) -> np.ndarray:
    """Evaluate the target f(theta) = exp(-i delta alpha cos theta) on a theta grid."""
    return np.exp(-1j * delta_alpha * np.cos(theta))


def evaluate_polynomial_on_circle(c: np.ndarray, theta: np.ndarray) -> np.ndarray:
    """Evaluate P(z) = sum_m c[m] z^m at z = exp(i theta)."""
    deg = len(c) - 1
    z = np.exp(1j * theta)
    out = np.zeros_like(z, dtype=complex)
    zk = np.ones_like(z, dtype=complex)
    for m in range(deg + 1):
        out = out + c[m] * zk
        zk = zk * z
    return out


def sup_norm_error(d: int, delta_alpha: float, ngrid: int = 4 * 4096) -> tuple[float, float, float]:
    """Return (sup_norm_err, P_max, tail_bound) for the Jacobi-Anger Laurent truncation."""
    theta = np.linspace(0.0, 2 * np.pi, ngrid, endpoint=False)
    c = jacobi_anger_coeffs(d, delta_alpha)
    P_vals = evaluate_polynomial_on_circle(c, theta)
    L_vals = np.exp(-1j * d * theta) * P_vals  # |z|=1, so |L_d| = |P|
    target = laurent_target_cos(delta_alpha, theta)
    err = float(np.max(np.abs(L_vals - target)))
    Pmax = float(np.max(np.abs(P_vals)))
    tail = 2.0 * (np.e * delta_alpha / (2 * (d + 1))) ** (d + 1)
    return err, Pmax, tail
