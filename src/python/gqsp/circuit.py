"""Qiskit GQSP circuit construction (MW-native form, Motlagh-Wiebe Eq. 42).

Two equivalent realisations of the Laurent target L_d(z) = z^{-d} P(z) are supported:

* `realise="L_d"` (Form C, MW2024 Eq. 52): all 2d slots are open-controlled W
  (the A = |0><0| (x) U + |1><1| (x) I of MW Eq. 5, fires on QSP=|0>),
  followed by an uncontrolled W^{-d} = (W^dagger)^d tail. Cost: 3d block-
  encoding queries.
* `realise="L_d_form_B"` (Form B, MW2024 Eq. 46): the slot pattern uses
  open-controlled W on one half of the slots and closed-controlled
  W^dagger (A' = |0><0| (x) I + |1><1| (x) U^dagger of MW Eq. 45 — fires on
  QSP=|1>) on the other half, with no uncontrolled tail. Cost: 2d block-
  encoding queries (1.5x cheaper than Form C). Same MW+BS angles.

MW-native time ordering: L_0 is applied FIRST in time. The 2d layers
follow with angle index j = 1..2d INCREASING in time (each layer is a
controlled signal slot followed by R_MW(theta_j, phi_j) on the QSP qubit).
For Form C, the uncontrolled W^{-d} tail is applied LAST in time.

Matrix product convention for the unitary realised in the QSP=|0>, anc=|0>
post-selected block:

    M(z) = ( prod_{j=N..1} R_MW(theta_j, phi_j) D(z) ) * L_0(lambda, theta_0, phi_0)

with D(z) = diag(z, 1), R_MW(theta, phi) = [[e^{iphi} c, e^{iphi} s], [s, -c]]
(= R_BS^T) and L_0 = R_MW(theta_0, phi_0) * diag(e^{i lambda}, 1) (= L_0_BS^T).
The matrix product reads left-to-right with the j=N factor leftmost; in time,
the j=N gate is applied last and the j=1 gate is applied first after L_0.
The post-selected polynomial lives in column 0 of M (M[0,0] = P(z)).

For ``_gqsp_circuit_bs_legacy`` we keep the BS-form circuit (L_0 last in time,
R_BS gates, j = N..1 left-to-right) for the equivalence regression test.

We use a 1-ancilla qubitization block encoding for testing:
    U_H = [[ H,             sqrt(I - H^2) ],
           [ sqrt(I - H^2), -H            ]]
with alpha = 1, ||H|| <= 1. The walk is W = (Z_anc otimes I_sys) * U_H.
"""
from __future__ import annotations

import numpy as np
from qiskit import QuantumCircuit, QuantumRegister
from qiskit.circuit.library import UnitaryGate


def build_block_encoding_one_anc(H: np.ndarray) -> np.ndarray:
    """1-ancilla qubitization block encoding U_H of a Hermitian H with ||H|| <= 1.

    Returns the (2*n_sys) x (2*n_sys) unitary matrix in the basis
    {|0>_anc, |1>_anc} otimes {|0>, ..., |n_sys-1>}.
    """
    H = np.asarray(H, dtype=complex)
    if not np.allclose(H, H.conj().T):
        raise ValueError("H must be Hermitian")
    if np.linalg.norm(H, ord=2) > 1 + 1e-12:
        raise ValueError("||H|| must be <= 1 for this construction")
    eigvals, eigvecs = np.linalg.eigh(H)
    sqrtIH = eigvecs @ np.diag(np.sqrt(np.maximum(1 - eigvals ** 2, 0))) @ eigvecs.conj().T
    n = H.shape[0]
    UH = np.zeros((2 * n, 2 * n), dtype=complex)
    UH[:n, :n] = H
    UH[:n, n:] = sqrtIH
    UH[n:, :n] = sqrtIH
    UH[n:, n:] = -H
    return UH


def build_walk(UH: np.ndarray) -> np.ndarray:
    """Walk operator W = (Z_anc tensor I_sys) * U_H."""
    n2 = UH.shape[0]
    n = n2 // 2
    Z_anc = np.zeros((n2, n2), dtype=complex)
    Z_anc[:n, :n] = np.eye(n)
    Z_anc[n:, n:] = -np.eye(n)
    return Z_anc @ UH


def _R_mw(theta: float, phi: float) -> np.ndarray:
    """MW-native single-qubit rotation R_MW(theta, phi) = R_BS(theta, phi)^T."""
    c = np.cos(theta)
    s = np.sin(theta)
    eφ = np.exp(1j * phi)
    return np.array([[eφ * c, eφ * s],
                     [s,      -c    ]], dtype=complex)


def _L0_mw(lam: float, theta0: float, phi0: float) -> np.ndarray:
    """MW-native L_0 = R_MW(theta_0, phi_0) * diag(e^{i lambda}, 1) = L_0_BS^T."""
    c = np.cos(theta0)
    s = np.sin(theta0)
    return np.array([[np.exp(1j * (lam + phi0)) * c, np.exp(1j * phi0) * s],
                     [np.exp(1j * lam) * s,          -c                  ]], dtype=complex)


def gqsp_circuit(
    UH: np.ndarray,
    lam: float,
    thetas: np.ndarray,
    phis: np.ndarray,
    d: int,
    *,
    realise: str = "L_d_form_B",
) -> QuantumCircuit:
    """Build the Qiskit GQSP circuit on (QSP, anc, sys) in MW-native ordering.

    Args:
        UH: 1-ancilla block encoding unitary (2*n_sys x 2*n_sys numpy array).
        lam, thetas, phis: angles from `gqsp_extract_angles` (MW-native), where
            thetas[0], phis[0] parameterise L_0 with global phase lambda, and
            thetas[j], phis[j] for j = 1..2d parameterise the j-th MW layer.
            Numerically these match BS-form angles for the same target P, but
            here they parameterise R_MW = R_BS^T gates and a transposed L_0.
        d: half the polynomial degree (so 2d slots).
        realise:
          * 'P' — realise the ordinary polynomial P(W) on the QSP=|0> block (no tail).
          * 'L_d' — Form C realisation of L_d(W) = W^{-d} P(W) (MW2024 Eq. 52):
              all 2d slots are open-controlled W (A = |0><0| (x) U + |1><1| (x) I,
              fires on QSP=|0>), followed by an uncontrolled W^{-d} = (W^dagger)^d
              tail on (anc, sys), applied LAST in time. Cost: 2d controlled-W +
              d uncontrolled-W^dagger = 3d block-encoding queries.
          * 'L_d_form_B' (default) — Form B realisation of L_d(W) (MW2024 Eq. 46):
              slots j = 1..d (FIRST in time, immediately after L_0) are open-
              controlled W, slots j = d+1..2d (LAST in time) are closed-
              controlled W^dagger. No uncontrolled tail. Cost: 2d block-
              encoding queries (1.5x fewer than Form C).

    Returns:
        QuantumCircuit acting on (1 QSP qubit, 1 LCU/block-encoding ancilla qubit,
        n_sys system qubits) where n_sys = log2(UH.shape[0] / 2).

    Qubit order: [QSP, anc, sys_0, sys_1, ...] (Qiskit little-endian within sys).
    """
    if realise not in ("P", "L_d", "L_d_form_B"):
        raise ValueError(
            f"realise must be 'P', 'L_d', or 'L_d_form_B', got {realise!r}"
        )

    n2 = UH.shape[0]
    n_sys_dim = n2 // 2
    n_sys_qubits = int(np.log2(n_sys_dim))
    if 2 ** n_sys_qubits != n_sys_dim:
        raise ValueError(f"UH dimension {n2} is not 2*2^n_sys")

    N = len(thetas) - 1
    if N != 2 * d:
        raise ValueError(f"len(thetas) - 1 = {N} should equal 2*d = {2*d}")
    if len(phis) != N + 1:
        raise ValueError("len(phis) must equal len(thetas)")

    qsp = QuantumRegister(1, "qsp")
    anc = QuantumRegister(1, "anc")
    sys = QuantumRegister(n_sys_qubits, "sys")
    qc = QuantumCircuit(qsp, anc, sys)

    W = build_walk(UH)
    W_dag = W.conj().T

    # MW-native time order:
    #   1. L_0 first (rightmost factor in matrix product M = (prod R_MW D) L_0).
    #   2. For j = 1, 2, ..., N: signal slot, then R_MW(theta_j, phi_j) on QSP.
    #      (Each matrix factor (R_MW(theta_j) D) reads in time as: D first
    #      (= controlled signal slot), then R_MW second.)
    #   3. Form C: append uncontrolled W^{-d} tail LAST in time
    #      (= leftmost factor of total product W^{-d} * M_MW).
    L0 = _L0_mw(lam, thetas[0], phis[0])
    qc.append(UnitaryGate(L0, label="L_0"), [qsp[0]])

    for j in range(1, N + 1):
        # Controlled signal slot (D = ctrl-W) first
        if realise == "L_d_form_B" and j > d:
            # Slots j = d+1..2d (last d in time): closed-controlled W^dagger.
            # A' = MW Eq. 45, fires on QSP=|1>. No uncontrolled tail.
            ctrl_Wdag = UnitaryGate(W_dag, label="W†").control(
                num_ctrl_qubits=1, ctrl_state=1
            )
            qc.append(ctrl_Wdag, [qsp[0], *sys, anc[0]])
        else:
            # Open-controlled W (A = MW Eq. 5, fires on QSP=|0>).
            ctrl_W = UnitaryGate(W, label="W").control(
                num_ctrl_qubits=1, ctrl_state=0
            )
            qc.append(ctrl_W, [qsp[0], *sys, anc[0]])

        # R_MW(theta_j, phi_j) on QSP after the signal slot
        Rj = _R_mw(thetas[j], phis[j])
        qc.append(UnitaryGate(Rj, label=f"R_{j}"), [qsp[0]])

    if realise == "L_d":
        # Form C tail: uncontrolled W^{-d} applied LAST in time on (anc, sys).
        # Total unitary = W^{-d} * M_MW gives L_d(W) on the QSP=|0> block.
        Wmd = np.linalg.matrix_power(W_dag, d)
        qc.append(UnitaryGate(Wmd, label=f"W^-{d}"), [*sys, anc[0]])

    return qc


# ---------------------------------------------------------------------------
# BS-form (legacy) — kept so the regression test can build the BS-form circuit
# from BS-form angles and verify that it realises the same polynomial as the
# MW-native circuit. Do not use in production thesis paths.
# ---------------------------------------------------------------------------
def _gqsp_circuit_bs_legacy(
    UH: np.ndarray,
    lam: float,
    thetas: np.ndarray,
    phis: np.ndarray,
    d: int,
    *,
    realise: str = "L_d_form_B",
) -> QuantumCircuit:
    """BS-form circuit: M_BS = L_0_BS prod_{j=1..N} D R_BS(theta_j, phi_j),
    L_0 last in time, R_BS gates, slot order matching BS Eq. 1.2.

    Same Form B / Form C semantics as the MW-native gqsp_circuit; only the
    matrix-product convention and gate matrices differ.
    """
    if realise not in ("P", "L_d", "L_d_form_B"):
        raise ValueError(
            f"realise must be 'P', 'L_d', or 'L_d_form_B', got {realise!r}"
        )

    n2 = UH.shape[0]
    n_sys_dim = n2 // 2
    n_sys_qubits = int(np.log2(n_sys_dim))
    if 2 ** n_sys_qubits != n_sys_dim:
        raise ValueError(f"UH dimension {n2} is not 2*2^n_sys")

    N = len(thetas) - 1
    if N != 2 * d:
        raise ValueError(f"len(thetas) - 1 = {N} should equal 2*d = {2*d}")
    if len(phis) != N + 1:
        raise ValueError("len(phis) must equal len(thetas)")

    qsp = QuantumRegister(1, "qsp")
    anc = QuantumRegister(1, "anc")
    sys = QuantumRegister(n_sys_qubits, "sys")
    qc = QuantumCircuit(qsp, anc, sys)

    W = build_walk(UH)
    W_dag = W.conj().T

    for j in range(N, 0, -1):
        # R_BS(theta_j, phi_j) on QSP first (= rightmost in matrix factor D R_BS,
        # which is leftmost-in-time within this iteration of the matrix product
        # since matrix product reads right-to-left in time).
        c_th, s_th = np.cos(thetas[j]), np.sin(thetas[j])
        eφ = np.exp(1j * phis[j])
        Rj = np.array([[eφ * c_th, s_th],
                       [eφ * s_th, -c_th]], dtype=complex)
        qc.append(UnitaryGate(Rj, label=f"R_{j}"), [qsp[0]])

        if realise == "L_d_form_B" and j > d:
            ctrl_Wdag = UnitaryGate(W_dag, label="W†").control(
                num_ctrl_qubits=1, ctrl_state=1
            )
            qc.append(ctrl_Wdag, [qsp[0], *sys, anc[0]])
        else:
            ctrl_W = UnitaryGate(W, label="W").control(
                num_ctrl_qubits=1, ctrl_state=0
            )
            qc.append(ctrl_W, [qsp[0], *sys, anc[0]])

    # L_0_BS = diag(e^{i lambda}, 1) * R_BS(theta_0, phi_0) on QSP last
    c_th0, s_th0 = np.cos(thetas[0]), np.sin(thetas[0])
    L0 = np.array([
        [np.exp(1j * (lam + phis[0])) * c_th0, np.exp(1j * lam) * s_th0],
        [np.exp(1j * phis[0]) * s_th0, -c_th0],
    ], dtype=complex)
    qc.append(UnitaryGate(L0, label="L_0"), [qsp[0]])

    if realise == "L_d":
        Wmd = np.linalg.matrix_power(W_dag, d)
        qc.append(UnitaryGate(Wmd, label=f"W^-{d}"), [*sys, anc[0]])

    return qc
