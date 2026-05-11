"""Motlagh-Wiebe Algorithm 1: extract GQSP angles from (P, Q).

Implements Theorem 1 / Algorithm 1 of
  Motlagh & Wiebe, "Generalized Quantum Signal Processing", arXiv:2308.01501.

We use the MW-native matrix form (their Eq. 42):

    [ P(z)   *
      Q(z)   * ]  =  ( prod_{j=N..1} R_MW(theta_j, phi_j) D(z) )
                       L_0(lambda, theta_0, phi_0)

where
    D(z)              = diag(z, 1)
    R_MW(theta, phi)  = [ e^{i phi} cos theta   e^{i phi} sin theta
                          sin theta            -cos theta            ]
    L_0(λ, θ_0, φ_0)  = R_MW(theta_0, phi_0) · diag(e^{i lambda}, 1)
                      = [ e^{i (lambda + phi_0)} cos theta_0   e^{i phi_0} sin theta_0
                          e^{i lambda} sin theta_0             -cos theta_0           ]

R_MW is the transpose of the rotation matrix R_BS used in the
Berntson--Sünderhauf paper (BS Eq. 1.2), and L_0 is the transpose of L_0_BS.
In MW-native form, L_0 is the *first* gate applied in time (rightmost in the
matrix product). The layer indices increase in time: layer j = 1 is applied
just after L_0, and layer j = N is applied last (= leftmost factor in the
matrix product).

The post-selected polynomial lives in column 0 of M (M[0,0] = P(z),
M[1,0] = Q(z)). Because M_MW = M_BS^T for matched angles, column 0 of
M_MW equals row 0 of M_BS, and the same polynomial pair (p, q) feeds both
recursions. The numerical angle values produced are identical to the BS-form
output --- what differs is the gate matrix (R_MW = R_BS^T) and the time
ordering of L_0 in the realised circuit.

Recursion: peel V_j = R_MW(theta_j, phi_j) D(z) from the LEFT, iterating
j = N, N-1, ..., 1 so that iteration index = angle index. At each step,
(theta_j, phi_j) is determined from the leading coefficients of the
column-0 polynomials (P_{j-1}, Q_{j-1}); after peeling, P_j loses its
constant term (divided by z), and Q_j loses its leading term.

The BS-form recursion is preserved as ``_extract_angles_bs_legacy`` /
``_reconstruct_PQ_bs_legacy`` for the equivalence regression test only.
"""
from __future__ import annotations

import numpy as np


def gqsp_extract_angles(p: np.ndarray, q: np.ndarray, tol: float = 1e-8):
    """Extract MW-native angles. Returns (lambda, thetas, phis) with len(thetas) = N+1, len(phis) = N+1.

    thetas[0], phis[0] are the L_0 angles (with global phase lambda).
    thetas[j], phis[j] for j >= 1 are the layer-j single-qubit rotation angles
    (R_MW(theta_j, phi_j) sandwiched with D(z) in the matrix product).

    Inputs (p, q) are the Laurent coefficients of column 0 of M:
        p = M[0, 0],  q = M[1, 0],
    stored low-to-high (p[0] = constant, p[N] = leading z^N coefficient).
    """
    N = len(p) - 1
    if len(q) != len(p):
        raise ValueError("P and Q must have the same number of coefficients")

    P = p.astype(complex).copy()
    Q = q.astype(complex).copy()
    thetas = np.zeros(N + 1)
    phis = np.zeros(N + 1)

    # Iterate j = N, N-1, ..., 1: peel the LEFTMOST factor of the matrix product
    # (= R_MW(theta_j, phi_j) D for the current j) at each step. Iteration index
    # equals angle index, so thetas[j] is the angle of the j-th gate in matrix-
    # product order (and the angle of the j-th gate from the right in time).
    for j in range(N, 0, -1):
        last = len(P) - 1  # leading coefficient index = j
        a = P[last]
        b = Q[last]
        ra = abs(a)
        rb = abs(b)
        rho = np.hypot(ra, rb)
        if rho < tol:
            # Effective degree below current length: layer is trivial.
            thetas[j] = 0.0
            phis[j] = 0.0
            P = P[:-1]
            Q = Q[:-1]
            continue
        c_theta = ra / rho
        s_theta = rb / rho
        theta = np.arctan2(rb, ra)
        phi = np.angle(a) - np.angle(b)
        thetas[j] = theta
        phis[j] = phi

        # Inverse layer (peel R_MW(theta_j, phi_j) D from the LEFT):
        #   z P_j(z) = e^{-i phi} cos theta * P_{j-1}(z) + sin theta * Q_{j-1}(z)
        #   Q_j(z)   = e^{-i phi} sin theta * P_{j-1}(z) - cos theta * Q_{j-1}(z)
        # P_j has degree one lower than P_{j-1} (we divide by z, dropping the
        # constant term); Q_j has degree one lower (we drop the leading term).
        e_mphi = np.exp(-1j * phi)
        Pnew = e_mphi * c_theta * P + s_theta * Q
        Qnew = e_mphi * s_theta * P - c_theta * Q
        if abs(Pnew[0]) > tol:
            raise RuntimeError(
                f"Constant term of z*P_{{j}} should vanish at j={j}, got {abs(Pnew[0])}"
            )
        if abs(Qnew[last]) > tol:
            raise RuntimeError(
                f"Leading coef of Q_{{j}} should vanish at j={j}, got {abs(Qnew[last])}"
            )
        P = Pnew[1:]   # divide by z (drop constant, shift indices down)
        Q = Qnew[:-1]  # truncate to drop leading

    # Base case: P, Q are constants = column 0 of M_N = L_0
    #   M_N[0, 0] = e^{i (lambda + phi_0)} cos theta_0
    #   M_N[1, 0] = e^{i lambda} sin theta_0
    P0 = P[0]
    Q0 = Q[0]
    rP = abs(P0)
    rQ = abs(Q0)
    rho = np.hypot(rP, rQ)
    if rho < tol:
        raise RuntimeError("Base-case constants vanish")
    thetas[0] = np.arctan2(rQ, rP)
    lam = np.angle(Q0)              # arg(Q0) = lambda
    phis[0] = np.angle(P0) - lam    # arg(P0) - lambda = phi_0

    return lam, thetas, phis


def reconstruct_PQ(lam: float, thetas: np.ndarray, phis: np.ndarray):
    """Reconstruct (P, Q) = column 0 of M by direct matrix product.

    Builds M(z) = ( prod_{j=1..N} R_MW(theta_j, phi_j) D(z) ) * L_0
    polynomial-entry by polynomial-entry, then returns column 0 = (M[0,0], M[1,0]).
    Used for round-trip checks of `gqsp_extract_angles`.
    """
    N = len(thetas) - 1
    if len(phis) != N + 1:
        raise ValueError("thetas and phis must have the same length")

    # M_N = L_0(lambda, theta_0, phi_0) — constant 2x2 polynomial matrix:
    #   [ e^{i(lam+phi_0)} c_0,  e^{i phi_0} s_0;
    #     e^{i lam} s_0,         -c_0                 ]
    c_th0, s_th0 = np.cos(thetas[0]), np.sin(thetas[0])
    A = np.array([np.exp(1j * (lam + phis[0])) * c_th0], dtype=complex)
    B = np.array([np.exp(1j * phis[0]) * s_th0], dtype=complex)
    C = np.array([np.exp(1j * lam) * s_th0], dtype=complex)
    D = np.array([-c_th0], dtype=complex)

    def poly_shift(a, k):
        return np.concatenate([np.zeros(k, dtype=complex), a])

    def poly_add(a, b):
        n = max(len(a), len(b))
        out = np.zeros(n, dtype=complex)
        out[: len(a)] += a
        out[: len(b)] += b
        return out

    # Forward: build M = (prod_{j=N..1} R_MW(theta_j, phi_j) D(z)) * L_0 by
    # left-multiplying R_MW(theta_j, phi_j) D(z) onto the running product,
    # iterating j = 1, 2, ..., N so that the j=N factor ends up leftmost
    # (= last applied in time) and j=1 ends up just left of L_0
    # (= second applied in time, immediately after L_0).
    # Acting on a 2x2 polynomial matrix [[A, B], [C, D]]:
    #   D-step (left-multiply by diag(z, 1)): row 0 -> z*row 0.
    #   R_MW-step (left-multiply by [[e^{iphi}c, e^{iphi}s], [s, -c]]):
    #     A_new = e^{iphi} c * A_top + e^{iphi} s * C_bot     (where A_top = z*A, after D)
    #     B_new = e^{iphi} c * B_top + e^{iphi} s * D_bot
    #     C_new = s * A_top - c * C_bot
    #     D_new = s * B_top - c * D_bot
    for j in range(1, N + 1):
        c_th, s_th = np.cos(thetas[j]), np.sin(thetas[j])
        e_phi = np.exp(1j * phis[j])
        A_top = poly_shift(A, 1)  # z * A
        B_top = poly_shift(B, 1)
        C_bot = C
        D_bot = D
        A_new = poly_add(e_phi * c_th * A_top, e_phi * s_th * C_bot)
        B_new = poly_add(e_phi * c_th * B_top, e_phi * s_th * D_bot)
        C_new = poly_add(s_th * A_top, -c_th * C_bot)
        D_new = poly_add(s_th * B_top, -c_th * D_bot)
        A, B, C, D = A_new, B_new, C_new, D_new
    return A, C  # column 0 gives (P, Q)


# ---------------------------------------------------------------------------
# BS-form (legacy) — kept as private helpers so the regression test can verify
# that BS-form angles and MW-native angles realise the same polynomial P(W).
# Do not use these in production thesis paths; ``gqsp_extract_angles`` /
# ``reconstruct_PQ`` above are MW-native and are the canonical interface.
# ---------------------------------------------------------------------------
def _gqsp_extract_angles_bs_legacy(p: np.ndarray, q: np.ndarray, tol: float = 1e-8):
    """BS-form angle extraction. Matrix form (BS Eq. 1.2 = transpose of MW Eq. 7):

        [ P(z)   Q(z)
          *      *  ]  =  L_0(lambda, theta_0, phi_0)
                           prod_{j=1..N} D(z) R_BS(theta_j, phi_j)

    with R_BS(theta, phi) = [[e^{iphi} c, s], [e^{iphi} s, -c]] and
    L_0 = diag(e^{i lambda}, 1) * R_BS(theta_0, phi_0). Tracks the top row.
    Recursion peels D R_BS from the RIGHT, j = N, N-1, ..., 1.
    """
    N = len(p) - 1
    if len(q) != len(p):
        raise ValueError("P and Q must have the same number of coefficients")

    P = p.astype(complex).copy()
    Q = q.astype(complex).copy()
    thetas = np.zeros(N + 1)
    phis = np.zeros(N + 1)

    for j in range(N, 0, -1):
        a = P[j]
        b = Q[j]
        ra = abs(a)
        rb = abs(b)
        rho = np.hypot(ra, rb)
        if rho < tol:
            thetas[j] = 0.0
            phis[j] = 0.0
            P = P[:j]
            Q = Q[:j]
            continue
        c_theta = ra / rho
        s_theta = rb / rho
        theta = np.arctan2(rb, ra)
        phi = np.angle(a) - np.angle(b)
        thetas[j] = theta
        phis[j] = phi

        e_mphi = np.exp(-1j * phi)
        Pnew_top = e_mphi * c_theta * P + s_theta * Q
        Qnew = e_mphi * s_theta * P - c_theta * Q
        if abs(Pnew_top[0]) > tol:
            raise RuntimeError(
                f"BS legacy: constant of z*P_{{j-1}} should vanish at j={j}, got {abs(Pnew_top[0])}"
            )
        if abs(Qnew[j]) > tol:
            raise RuntimeError(
                f"BS legacy: leading of Q_{{j-1}} should vanish at j={j}, got {abs(Qnew[j])}"
            )
        P = Pnew_top[1:]
        Q = Qnew[:j]

    P0 = P[0]
    Q0 = Q[0]
    rP = abs(P0)
    rQ = abs(Q0)
    rho = np.hypot(rP, rQ)
    if rho < tol:
        raise RuntimeError("BS legacy: base-case constants vanish")
    thetas[0] = np.arctan2(rQ, rP)
    lam = np.angle(Q0)
    phis[0] = np.angle(P0) - np.angle(Q0)

    return lam, thetas, phis


def _reconstruct_PQ_bs_legacy(lam: float, thetas: np.ndarray, phis: np.ndarray):
    """BS-form reconstruction. Returns (P, Q) = top row of M = L_0 prod D R_BS."""
    N = len(thetas) - 1
    if len(phis) != N + 1:
        raise ValueError("thetas and phis must have the same length")

    c_th0, s_th0 = np.cos(thetas[0]), np.sin(thetas[0])
    A = np.array([np.exp(1j * (lam + phis[0])) * c_th0], dtype=complex)
    B = np.array([np.exp(1j * lam) * s_th0], dtype=complex)
    C = np.array([np.exp(1j * phis[0]) * s_th0], dtype=complex)
    D = np.array([-c_th0], dtype=complex)

    def poly_add(a, b):
        n = max(len(a), len(b))
        out = np.zeros(n, dtype=complex)
        out[: len(a)] += a
        out[: len(b)] += b
        return out

    def poly_shift(a, k):
        return np.concatenate([np.zeros(k, dtype=complex), a])

    def poly_scale(a, c):
        return c * a

    for j in range(1, N + 1):
        c_th, s_th = np.cos(thetas[j]), np.sin(thetas[j])
        eφ = np.exp(1j * phis[j])
        A_new = poly_add(poly_shift(poly_scale(A, eφ * c_th), 1), poly_scale(B, eφ * s_th))
        B_new = poly_add(poly_shift(poly_scale(A, s_th), 1), poly_scale(B, -c_th))
        C_new = poly_add(poly_shift(poly_scale(C, eφ * c_th), 1), poly_scale(D, eφ * s_th))
        D_new = poly_add(poly_shift(poly_scale(C, s_th), 1), poly_scale(D, -c_th))
        A, B, C, D = A_new, B_new, C_new, D_new
    return A, B  # top row gives (P, Q)
