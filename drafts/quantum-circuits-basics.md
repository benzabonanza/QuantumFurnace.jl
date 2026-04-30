# Quantum Circuits — Basics

> **Insertion target:** `supplementary-informations/1_preliminaries.tex`, body of `\subsection{Basics.}` inside `\chapter{Quantum Circuits}` (run-in style — title flows into the first paragraph on the same line).
> **Existing labels referenced:** `eq:HS`, `eq:post_meas_state`, `eq:kraus`, `sec:prelim-weak-meas`, `circ:trotter-strang`.
> **New labels introduced:** `eq:paulis`, `eq:Ry`.
> **New citation key:** `pulidoMateo2024arbitrary` (bibtex stub at the bottom).

## Qubits from bits

The smallest non-trivial system in $\mathcal{D}(\mathcal{H})$ is the qubit, the quantum analogue of a classical bit. Its Hilbert space is $\mathcal{H} = \mathbb{C}^2$ and a pure state takes the form $\ket{\psi} = \alpha\ket{0} + \beta\ket{1}$ with $|\alpha|^2 + |\beta|^2 = 1$. An $n$-qubit register lives in the tensor product $\mathcal{H}^{\otimes n} = (\mathbb{C}^2)^{\otimes n}$ of dimension $2^n$, with computational basis $\{\ket{x}\}_{x\in\{0,1\}^n}$ where $\ket{x} = \ket{x_1}\otimes\cdots\otimes\ket{x_n}$. Mixed states $\rho$ on this register are density operators on $\mathcal{H}^{\otimes n}$ in the sense of (§Quantum~states).

## Quantum gates and Paulis

A *quantum gate* acting on $k$ qubits is a unitary $U \in \mathcal{B}((\mathbb{C}^2)^{\otimes k})$. Any unitary on $n$ qubits can be approximated to arbitrary accuracy by composing gates from the universal set consisting of all single-qubit unitaries together with the two-qubit `\textsc{CNOT}` [Nielsen–Chuang] [CITE: bible]. The single-qubit Pauli matrices

$$
X = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}, \quad
Y = \begin{pmatrix} 0 & -\ii \\ \ii & 0 \end{pmatrix}, \quad
Z = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix},
$$ <!-- \label{eq:paulis} -->

are simultaneously Hermitian and unitary and satisfy $X^2 = Y^2 = Z^2 = \idm$. Labelling them $\sigma_i$ with $i\in\{x,y,z\}$, their commutator and anti-commutator read

$$
[\sigma_i, \sigma_j] = 2\ii\,\epsilon_{ijk}\,\sigma_k, \qquad \{\sigma_i, \sigma_j\} = 2\,\delta_{ij}\,\idm,
$$

with $\epsilon_{ijk}$ the totally antisymmetric Levi-Civita symbol and summation over the repeated index $k$ implied; equivalently, $\sigma_i\sigma_j = \delta_{ij}\,\idm + \ii\,\epsilon_{ijk}\,\sigma_k$. In particular distinct Paulis anti-commute, $\{X,Y\} = \{Y,Z\} = \{X,Z\} = 0$. Together with $\idm$ they span $\mathcal{B}(\mathbb{C}^2)$, and tensor products of Paulis form an orthogonal basis of $\mathcal{B}((\mathbb{C}^2)^{\otimes n})$ under the Hilbert–Schmidt inner product (Eq. eq:HS).

## Pauli rotations and the $Y_\theta$ shorthand

Each Pauli generates a one-parameter family of *Pauli rotations* $R_\sigma(\theta) := \ee^{-\ii\theta\sigma/2}$ for $\sigma\in\{X,Y,Z\}$; the same convention extends to multi-qubit Pauli strings (e.g. $R_{\sigma_j\sigma_k}$ in Fig. circ:trotter-strang). Of particular interest is

$$
R_Y(\theta) = \cos(\tfrac{\theta}{2})\,\idm - \ii\sin(\tfrac{\theta}{2})\,Y,
$$ <!-- \label{eq:Ry} -->

the workhorse of *amplitude encoding*: $R_Y(2\arcsin\sqrt{p})\ket{0} = \sqrt{1-p}\ket{0} + \sqrt{p}\ket{1}$ deterministically loads a probability $p\in[0,1]$ into the amplitude of $\ket{1}$. We will use the shorthand $Y_\theta := R_Y(2\arcsin\sqrt{\theta})$ throughout, particularly for the Boltzmann and weak-measurement ancillas in (§sec:prelim-weak-meas).

## Quantum circuits

A *quantum circuit* is a horizontal-time diagram in which each wire carries one qubit (or a bundled register), boxes denote gates applied at that time slice, and a filled dot $\bullet$ on a wire indicates a control: the gate on the target wire fires conditionally on the control being $\ket{1}$ (an open dot $\circ$ controls on $\ket{0}$). A half-filled dot \tikz[baseline=-0.5ex]{\node[halfctrl]{};} on a bundled register denotes a *uniform* (SELECT-style) control, under which the target gate is applied in superposition to every computational-basis state of the control register, each contributing its own branch to the state; this notation appears in the block-encoding and dissipation circuits of later chapters (Fig. fig:block-encoding-ba-halfctrl, Fig. fig:u-diss). All circuit figures in this chapter follow this convention via the `quantikz2` package.

## Measurement

Read-out is performed by *computational-basis measurement*, i.e. the projective measurement (Eq. eq:post_meas_state) associated with the spectral decomposition of the $Z$-operators on each qubit. On a pure state $\ket{\psi}\in(\mathbb{C}^2)^{\otimes n}$ this returns the bit-string $x\in\{0,1\}^n$ with probability $|\langle x|\psi\rangle|^2$ and collapses the register onto $\ket{x}$. *Mid-circuit* measurements are likewise allowed: outcomes can be discarded (yielding a non-unitary channel of Kraus form (Eq. eq:kraus)) or post-selected on a desired result --- both modes appear in the weak-measurement scheme of (§sec:prelim-weak-meas).

## From gates to quantum channels

The three primitives introduced so far --- a universal gate set, fresh ancilla qubits prepared in $\ket{0}$, and computational-basis measurement --- already suffice to realise *any* CPTP map $\Phi: \mathcal{D}(\mathcal{H}^{\otimes n}) \to \mathcal{D}(\mathcal{H}^{\otimes n})$. By Stinespring's dilation theorem (see §OQS for the channel formalism and [CITE: bible]), every quantum channel can be lifted to a unitary on a larger Hilbert space: there exist an $m$-qubit ancilla $A$ and a joint unitary $U$ on $n+m$ qubits such that

$$
\Phi(\rho) = \tr_A\!\bigl[\,U\,\bigl(\rho \otimes \ket{0}\bra{0}^{\otimes m}\bigr)\,U^\dagger\,\bigr],
$$

with the partial trace implemented either by physically discarding the ancilla or, equivalently, by measuring it in the computational basis and forgetting the outcome --- the resulting average over unread outcomes reproduces the Kraus form (Eq. eq:kraus). Since single-qubit unitaries together with `\textsc{CNOT}` decompose $U$ exactly, a circuit implementing $\Phi$ to any desired accuracy requires only the three primitives above; this is the principle we invoke implicitly whenever a Lindbladian-generated channel (§sec:OQS) is compiled to gates in later chapters.

## Hardware timescales (practical aside)

A practical aside on *timescales*. The gates above are physical processes on real hardware, and their wall-clock duration is a hard constant in any complexity statement. On state-of-the-art *superconducting* processors, single- and two-qubit gates run in tens to hundreds of nanoseconds: a recent survey of IBM's fleet reports a median CZ two-qubit gate time of $68\,\text{ns}$ on the 156-qubit Heron r2 (`ibm_fez`) and $84\,\text{ns}$ on Heron r1 (`ibm_torino`), against $460$–$665\,\text{ns}$ for the echoed cross-resonance (ECR) two-qubit gate on the preceding Eagle r3 generation [AbuGhanem 2025] [CITE: abughanem2025ibm]; the latest Heron r3 retains the same tunable-coupler CZ architecture, pushing median two-qubit fidelity above $99.9\%$. On the *trapped-ion* platform pursued at Leibniz Universität Hannover and PTB Braunschweig within Quantum Valley Lower Saxony, the most recent fully chip-integrated two-qubit register on $^{9}\text{Be}^{+}$ runs microwave-driven Mølmer–Sørensen gates of duration $\tau \approx 1.15\,\text{ms}$ with composite process fidelity $96.6(4)\%$ under cycle benchmarking [Pulido-Mateo et al. 2024] [CITE: pulidoMateo2024arbitrary]. By contrast, a logic gate on classical CMOS hardware switches in $\sim 10\,\text{ps}$. The quantum–classical wall-clock gap is therefore three to eight orders of magnitude per gate, a constant overhead that any claim of *quantum advantage* for Gibbs sampling must amortise over the gate count.

---

## Citations

- **`bible`** — Nielsen, M.~A. & Chuang, I.~L. *Quantum Computation and Quantum Information*. Cambridge University Press (already in `references.bib`).
- **`pulidoMateo2024arbitrary`** — Pulido-Mateo, Mendpara, Duwe, Dubielzig, Zarantonello, Krinner, Ospelkaus, *"Arbitrary quantum circuits on a fully integrated two-qubit computation register for a trapped-ion quantum processor"*, **Phys. Rev. Research 6, L022067 (2024)**, arXiv:2403.19809.
- **`abughanem2025ibm`** — AbuGhanem, M., *"IBM quantum computers: evolution, performance, and future directions"*, **The Journal of Supercomputing 81, 687 (2025)**, arXiv:2410.00916. Reports per-backend calibration data (gate durations, fidelities, coherence times) across IBM's Canary → Eagle → Heron lineage; the numbers cited here are from its Heron / Eagle summary tables.

```bibtex
@article{pulidoMateo2024arbitrary,
  author        = {Pulido-Mateo, N. and Mendpara, H. and Duwe, M. and Dubielzig, T. and Zarantonello, G. and Krinner, L. and Ospelkaus, C.},
  title         = {Arbitrary quantum circuits on a fully integrated two-qubit computation register for a trapped-ion quantum processor},
  journal       = {Phys. Rev. Research},
  volume        = {6},
  pages         = {L022067},
  year          = {2024},
  doi           = {10.1103/PhysRevResearch.6.L022067},
  eprint        = {2403.19809},
  archivePrefix = {arXiv},
}

@article{abughanem2025ibm,
  author        = {AbuGhanem, M.},
  title         = {{IBM} quantum computers: evolution, performance, and future directions},
  journal       = {The Journal of Supercomputing},
  volume        = {81},
  number        = {6},
  pages         = {687},
  year          = {2025},
  doi           = {10.1007/s11227-025-07047-7},
  eprint        = {2410.00916},
  archivePrefix = {arXiv},
}
```

---

## Writing notes

- **Coverage**: 7 paragraphs (qubits-from-bits → gates + universal set + Paulis → Pauli rotations + $Y_\theta$ shorthand → circuit notation → measurement → gates-to-channels/Stinespring → hardware-timescale aside) + 3 displays.
- **Notation conventions used**: thesis macros `\ii`, `\ee`, `\idm` (defined in `main.tex`). Paulis are uppercase Roman `X, Y, Z`; `\sigma` appears only as a Pauli metavariable inside `R_\sigma(\theta)`, matching the Trotter caption convention `R_{\sigma_j\sigma_k}(\theta) \equiv \ee^{-\ii\theta\,\sigma_j\sigma_k/2}` (Fig. circ:trotter-strang).
- **Length**: comfortably within one printed page in the 11pt / a4paper / `\setstretch{1.1}` layout.
- **Inline half-filled control symbol**: rendered as `\tikz[baseline=-0.5ex]{\node[halfctrl]{};}`, which reuses the `halfctrl` tikz style already defined in `main.tex` (lines 28–42). No new preamble needed. If the inline form clutters body text, the user may define `\newcommand{\halfctrlsym}{\tikz[baseline=-0.5ex]{\node[halfctrl]{};}}` in the preamble and write `\halfctrlsym` instead.
- **Two draft-checker fixes already incorporated**:
  - "Mixed states $\rho$ … in the sense of (§Quantum~states)" replaces an earlier "as introduced in the preceding subsection" (which was wrong after the chapter break).
  - "post-selected on a desired result" replaces "fed back as classical controls" (the thesis never actually uses classical-feedback control).
- **Open minor items left to user judgement** (none mathematically wrong):
  - amplitude letters $\alpha,\beta$ briefly clash with $\beta$ for inverse temperature in later chapters;
  - $|\langle x|\psi\rangle|^2$ uses the textbook bracket form rather than the `\bra{}\ket{}` style used elsewhere;
  - `\tfrac` inside the $R_Y$ display vs. `\frac` elsewhere — purely stylistic.
