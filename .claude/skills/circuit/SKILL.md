---
name: circuit
description: Generate publication-quality quantum circuit diagrams in LaTeX using quantikz2. Use when the thesis needs a circuit figure.
argument-hint: "[circuit description, e.g. 'DissipativeStep from Algorithm 1b' or 'basic QSP structure']"
model: opus
effort: max
allowed-tools: Read Write Glob Grep Bash WebSearch WebFetch
---

# /circuit --- Draw a Quantum Circuit Figure

Produce a publication-quality quantum circuit diagram in LaTeX using the **quantikz2** TikZ library.

**Argument**: `$ARGUMENTS` --- description of what circuit to draw (e.g., "DissipativeStep weak measurement circuit", "QSP signal processing operator", "Strang splitting Trotterization for n=4").

## Process

1. **Parse input**: Determine what circuit is requested from `$ARGUMENTS`.

2. **Derive slug**: Create a short kebab-case name (e.g., `dissipative-step`, `qsp-basic`, `trotter-strang`).

3. **Research before drawing** --- you MUST read relevant context:
   - Read `.claude/skills/circuit/quantikz2-reference.md` for the full quantikz2 syntax reference
   - Read relevant pages of `supplementary-informations/thesis.pdf` to match notation and style
   - Read relevant drafts from `drafts/` that describe the circuit algorithmically (especially `drafts/dissipative-step.md`, `drafts/coherent-step.md`, `drafts/algorithm1-outer.md`, `drafts/algorithm-subsection.md`)
   - Read relevant reference papers from `supplementary-informations/` if the circuit comes from or extends a published construction (Chen et al., Ding et al., Low & Chuang, etc.)
   - Check existing circuits in `drafts/circuits/` to avoid duplicates and maintain visual consistency

4. **Write the circuit** to `drafts/circuits/<slug>.tex`. The file must be:
   - A **standalone compilable** LaTeX document (with `\documentclass`, preamble, `\begin{document}`)
   - Also contain the circuit as a **reusable fragment** wrapped in a `\newcommand` or `\begin{figure}` that can be `\input{}`-ed into the thesis
   - Well-commented explaining what each section of the circuit does

5. **Skip compilation**: Do NOT attempt to compile the LaTeX. The user will compile it themselves locally. Focus on writing correct quantikz2 syntax.

6. **Beads integration**: If a beads issue is associated (user provides ID or `bd search` finds one):
   - Mark in progress: `bd update <id> -s in_progress`
   - On completion: `bd note <id> "circuit drawn: drafts/circuits/<slug>.tex"` and `bd close <id> -r "circuit complete"`

7. **Report**: Path to the .tex file, brief description of the circuit, and any design choices made.

## Output Format

Each `.tex` file in `drafts/circuits/` must follow this structure:

```latex
% Circuit: <Human-readable title>
% Generated for: QuantumFurnace.jl PhD thesis
% Usage: \input{circuits/<slug>.tex} in main thesis, or compile standalone
%
\documentclass{article}
\usepackage[margin=1cm]{geometry}
\usepackage{tikz}
\usetikzlibrary{quantikz2}
\usepackage{amsmath,amssymb}
\usepackage{braket}

% Rounded corners on all gates globally
\tikzset{operator/.append style={rounded corners}}

\begin{document}
\pagestyle{empty}

\begin{figure}[ht]
  \centering
  \begin{quantikz}
    % ... circuit here ...
  \end{quantikz}
  \caption{<2--4 sentence caption explaining the circuit.>}
  \label{fig:<slug>}
\end{figure}

\end{document}
```

**Figure rule**: Every circuit MUST be wrapped in a `\begin{figure}[ht]` environment with `\centering`, a `\caption{...}`, and a `\label{fig:<slug>}`. This gives proper "Figure N: ..." numbering and allows `\ref{fig:<slug>}` cross-references. The standalone file uses `article` class with `\pagestyle{empty}` and tight margins so it compiles cleanly on its own. When embedding in the thesis, copy the entire `\begin{figure}...\end{figure}` block (omit the preamble). Keep captions to 2--4 sentences.

## Design Principles

### Visual Style
- **Rounded corners on all gates**: always include `\tikzset{operator/.append style={rounded corners}}` in the preamble. This applies rounded corners globally to every gate box. Do NOT add `rounded corners` to individual gate styles --- the global setting handles it.
- **Clean and readable**: prefer whitespace over density. Use `\qw` liberally for wire segments.
- **Consistent register labeling**: always use `\lstick{...}` for input labels with standard ket notation
- **Group related operations**: use `\gategroup` or `\slice{}` to visually separate stages of the circuit
- **Bundled registers**: use `\qwbundle{r}` for multi-qubit registers (e.g., frequency register $\Omega$ with $r$ qubits) rather than drawing every wire
- **Color sparingly**: use light fills (`fill=blue!10`) to highlight key gates or distinguish circuit stages, but don't overdo it

### Notation Conventions (match the thesis)
- **Upright math constants**: use `\ee` for Euler's number and `\ii` for the imaginary unit (defined as `\newcommand{\ee}{\mathrm{e}}` and `\newcommand{\ii}{\mathrm{i}}`). In standalone files, define these in the preamble. Never use bare `e` or `i` for these constants.
- System register: $S$ or $\ket{\psi}_S$
- Frequency/time registers: $\Omega$, $T_-$, $T_+$ with qubit counts $r$, $r_-$, $r_+$
- Boltzmann qubit: $q_\gamma$
- Weak-measurement qubit: $q_\delta$
- Jump operators: $A_a$ (Pauli operators in our setting)
- Trotterized evolution: $S_p^{(M)}(\bar t)$ or simplified as $\ee^{-\ii H\bar t}$
- Controlled operations: filled dot for control, gate box for target
- State preparations: labeled as $\text{Prep}$ or $\ket{b_\pm}$ or $\ket{f}$
- QFT: labeled box
- Measurements: `\meter{}` with optional basis label

### Circuit-Specific Guidance

**DissipativeStep (weak measurement)**:
- Registers: $\Omega$ (bundled, $r$ qubits), $S$ (system), $q_\gamma$, $q_\delta$
- Flow: StatePrep -> Ctrl-Trott$^-$ -> $A_a$ -> Ctrl-Trott$^+$ -> QFT -> Boltzmann rotation -> Weak measurement -> [branch on $q_\delta$: no-jump uncomputes, jump doesn't] -> Measure ancillas
- The controlled-on-$q_\delta$ branching can be shown with controlled gates or a dashed box
- Reference: \cite[Figure 3]{chen2023quantum}, but simplified for unitary Pauli jumps (no block-encoding ancillas)

**CoherentStep (block encoding via nested LCU)**:
- Registers: $T_-$ (bundled, $r_-$ qubits), $T_+$ (bundled, $r_+$ qubits), $S$ (system)
- Flow: Prep$|b_-\rangle$ -> Ctrl-Trott outer$^-$ -> Prep$|b_+\rangle$ -> [inner triple: Trott$^+$ -> $A_a^\dagger$ -> Trott$^{-2}$ -> $A_a$ -> Trott$^+$] -> Unprep$|b_+\rangle$ -> Ctrl-Trott outer$^+$ -> Unprep$|b_-\rangle$
- The inner Heisenberg evolution ($A^\dagger(\beta t') A(-\beta t')$) is the nested core
- Use `\gategroup` to highlight the inner LCU sandwich

**GQSP (Generalized Quantum Signal Processing)**:
- Show the signal-processing structure: alternating signal unitaries $U$ and phase rotations $e^{i\phi_k Z}$
- One ancilla qubit for the signal processing
- Polynomial degree $d$ determines number of repetitions
- Reference: Motlagh & Wiebe 2024

**Trotterization**:
- Show Strang splitting structure: $S_2(t) = e^{-iH_1 t/2} e^{-iH_2 t} e^{-iH_1 t/2}$
- For Heisenberg chain: show bond-by-bond structure with even/odd grouping
- Controlled version: controlled on time register qubits

**QPE / QFT**:
- Standard textbook circuits with thesis notation
- QPE uses controlled unitaries + inverse QFT

**LCU (Linear Combination of Unitaries)**:
- Generic LCU template: Prep -> Ctrl-$U_j$ -> Unprep -> measure ancilla
- Show how block encoding emerges from the sandwich structure

### Notes
- Compilation is done by the user locally --- do NOT run pdflatex/lualatex
- quantikz2 requires a recent TikZ (>= 3.1.10). If unsure about compatibility, note in comments that the user can fall back to `\usepackage{quantikz}` (v1 syntax)
- For very wide circuits, use `[row sep=0.3cm, column sep=0.4cm]` options on the quantikz environment to compress
- For circuits that don't fit on one line, split into sub-circuits with labeled wire connections ($\cdots$)
