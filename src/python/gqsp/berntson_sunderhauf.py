"""Berntson-Sünderhauf canonical complementary polynomial via FFT.

Implements Algorithm 1 from
  Berntson & Sünderhauf, "Complementary Polynomials in Quantum Signal Processing",
  Comm. Math. Phys. (2025) 406:161, arXiv:2406.04246.

Given P in C[z] with deg P = d and ||P||_{infty, T} <= 1 - delta for some
delta > 0, returns Q in C[z] with deg Q = d such that |P|^2 + |Q|^2 = 1
on the unit circle T, with Q outer (no zeros in the open unit disk).
"""
from __future__ import annotations

import numpy as np


def complementary_polynomial_bs(p: np.ndarray, N: int) -> np.ndarray:
    """Return monomial coefficients q[0..d] of canonical complementary polynomial Q.

    Args:
        p: monomial coefficients of P, length d+1.
        N: FFT length, must be a power of 2 with N > d. For accuracy epsilon,
           N ~ (d/epsilon) log(d/epsilon) is sufficient (BS Theorem 3).

    Returns:
        q: monomial coefficients of Q, length d+1.
    """
    d = len(p) - 1
    if not (N > d):
        raise ValueError(f"N must exceed degree d (got N={N}, d={d})")
    if N & (N - 1) != 0:
        raise ValueError(f"N must be a power of 2 (got {N})")

    # Step 1: evaluate P at N-th roots of unity by IFFT of (p, 0, ..., 0).
    p_pad = np.zeros(N, dtype=complex)
    p_pad[: d + 1] = p
    # numpy.fft.ifft has the 1/N normalization; scale by N to get sum P(omega_N^n).
    P_at_roots = N * np.fft.ifft(p_pad)  # P_at_roots[n] = P(omega_N^n), n = 0..N-1

    # Step 2: Fourier coefficients of log(1 - |P|^2) via FFT.
    f_at_roots = np.log(np.maximum(1.0 - np.abs(P_at_roots) ** 2, np.finfo(float).eps))
    a_hat = np.fft.fft(f_at_roots) / N
    # a_hat[k] is tilde a_n for n = k mod N (mapping n in -N/2+1..N/2 to k in 0..N-1
    # via k = n if n >= 0 else N + n).

    # Step 3: Apply Fourier multiplier Pi: keep n=0 halved, n in 1..N/2 unchanged,
    # zero out n in -N/2+1..-1 (mapped to N/2+1..N-1).
    g = np.zeros(N, dtype=complex)
    g[0] = a_hat[0] / 2  # n = 0
    g[1 : N // 2 + 1] = a_hat[1 : N // 2 + 1]  # n = 1..N/2
    # n = N/2+1..N-1 set to 0 (already by initialization)

    # Step 4: log Q at roots = N * IFFT(g).
    G_at_roots = N * np.fft.ifft(g)
    Q_at_roots = np.exp(G_at_roots)

    # Step 5: q[m] = (1/N) sum_n Q(omega_N^n) omega_N^{-mn} = FFT(Q_at_roots)[m] / N.
    q_full = np.fft.fft(Q_at_roots) / N

    # Truncate to degree d.
    q = q_full[: d + 1]
    return q
