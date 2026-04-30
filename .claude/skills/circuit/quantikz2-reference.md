# quantikz2 Reference Guide

Comprehensive syntax reference for the quantikz2 TikZ library (v5+), the de facto standard for quantum circuit diagrams in LaTeX. Based on the official tutorial by Alastair Kay (arXiv:1809.03842).

## Preamble

```latex
\usepackage{tikz}
\usetikzlibrary{quantikz2}
```

For arXiv submissions, bundle `tikzlibraryquantikz2.code.tex` alongside your `.tex` source.

If quantikz2 is unavailable (older TeX distributions), fall back to:
```latex
\usepackage{quantikz}
```

## Basic Structure

Circuits are laid out as a matrix: rows = wires, columns = time steps. Ampersands `&` separate columns, `\\` separates rows.

```latex
\begin{quantikz}
  & \gate{H} & \ctrl{1} & \qw \\
  & \qw      & \targ{}  & \qw
\end{quantikz}
```

Every cell that carries a wire but has no gate must contain `\qw` (quantum wire).

### Environment Options

```latex
\begin{quantikz}[row sep={0.6cm,between origins}, column sep=0.5cm]
  ...
\end{quantikz}
```

Key options:
- `row sep=0.5cm` --- vertical spacing between wires
- `column sep=0.4cm` --- horizontal spacing between time steps
- `thin lines` --- thinner default wires
- `transparent` --- makes gate backgrounds transparent
- `wire types={q,q,c,b}` --- set wire types per row: `q`=quantum, `c`=classical, `b`=bundled, `n`=none
- `slice all` --- add vertical slice lines at every column

## Wire Types and Labels

### Input/Output Labels

```latex
\lstick{\ket{0}}     % left label (input)
\rstick{\ket{\psi}}  % right label (output)
\lstick[3]{$S$}      % label spanning 3 wires
\midstick{=}         % mid-circuit label
```

### Wire Types

```latex
\qw          % quantum wire (default)
\cw          % classical wire
\qwbundle{r} % bundled wire (shows "/" with label r)
\setwiretype{c}  % change wire type mid-circuit (placed after a gate)
\setwiretype{n}  % no wire (invisible)
```

### Bundled (Multi-Qubit) Registers

For registers like $\Omega$ with $r$ qubits, use bundled wires:

```latex
\begin{quantikz}[wire types={b,q}]
  \lstick{$\Omega$\;($r$)} & \gate{\text{Prep}} & \qwbundle{r} & \gate[2]{U} & \qw \\
  \lstick{$\ket{\psi}_S$}  & \qw                & \qw          &             & \qw
\end{quantikz}
```

Or set it inline:
```latex
\lstick{$T_-$} & \qwbundle{r_-} & \gate{\text{Prep}} & ...
```

## Gates

### Single-Qubit Gates

```latex
\gate{H}              % named gate (Hadamard)
\gate{X}              % Pauli-X
\gate{R_Z(\theta)}    % parameterized gate
\gate[style={fill=blue!10}]{U}  % styled gate
\phase{\alpha}         % phase gate (small dot with label)
\gate{R_Y(\theta)}    % rotation gate
```

### Multi-Qubit Gates (Spanning Wires)

```latex
% Gate spanning 2 adjacent wires
\begin{quantikz}
  & \gate[2]{U_{B_a}} & \qw \\
  &                    & \qw
\end{quantikz}

% Gate spanning 3 wires
\begin{quantikz}
  & \gate[3]{S_2^{(M)}(t)} & \qw \\
  &                         & \qw \\
  &                         & \qw
\end{quantikz}

% Gate spanning non-adjacent wires (skipping wire 2)
\begin{quantikz}
  & \gate[3,disable auto height]{U} & \qw \\
  & \qw                              & \qw \\  % this wire passes through
  &                                  & \qw
\end{quantikz}
```

The number in `\gate[N]{...}` is the number of wires the gate spans. The gate box is drawn from the current row downward for N rows. Leave the spanned rows' gate cells empty.

### Gate Styling

```latex
\gate[style={fill=yellow!20, draw=blue}]{U}        % colored fill and border
\gate[style={dashed}]{U}                            % dashed border
\gate[label style={font=\scriptsize}]{R_Y(\theta)}  % smaller label
\gate[style={fill=green!10}, label style={blue}]{A_a} % both
```

## Controlled Operations

### Controlled Gates (Filled Dot)

```latex
\ctrl{1}     % control dot, target 1 row below
\ctrl{-1}    % control dot, target 1 row above
\ctrl{2}     % control dot, target 2 rows below
```

### Targets

```latex
\targ{}      % circled-plus (CNOT target)
\control{}   % filled dot (for CZ: control on both ends)
\targX{}     % cross (SWAP target)
```

### Open Controls (White Dot)

```latex
\octrl{1}    % open control, target 1 row below
```

### Controlled-U

```latex
\begin{quantikz}
  & \ctrl{1}           & \qw \\
  & \gate{R_Y(\theta)} & \qw
\end{quantikz}
```

### Multi-Controlled Gates

```latex
% Toffoli (CCX)
\begin{quantikz}
  & \ctrl{1} & \qw \\
  & \ctrl{1} & \qw \\
  & \targ{}  & \qw
\end{quantikz}
```

### Controlled on Bundled Register

For controlled operations where the control is a multi-qubit register (e.g., controlled Hamiltonian evolution controlled on time register), draw the control from the bundled wire:

```latex
\begin{quantikz}[wire types={b,q}]
  \lstick{$T$} & \qwbundle{r} & \ctrl{1}                   & \qw \\
  \lstick{$S$} & \qw          & \gate{S_2^{(M)}(\bar{t})}  & \qw
\end{quantikz}
```

## SWAP

```latex
\begin{quantikz}
  & \swap{1}  & \qw \\
  & \targX{}  & \qw
\end{quantikz}
```

## Measurements

```latex
\meter{}       % standard measurement (meter symbol)
\meter{Z}      % measurement with basis label
\meter[2]{B}   % two-wire measurement (Bell)
\meterD{0}     % measurement with outcome label
```

### Measurement with Classical Output

```latex
\begin{quantikz}
  & \gate{H} & \meter{} \\
\end{quantikz}
```

To continue with a classical wire after measurement:
```latex
& \meter{} \setwiretype{c} & \ctrl{1} & \cw \\
```

## Slices (Vertical Separators)

Useful for marking stages of an algorithm:

```latex
\begin{quantikz}
  & \gate{\text{Prep}} & \slice{$U_\text{diss}$} & \gate{\text{QFT}} & \slice{filter} & \gate{R_Y} & \qw \\
\end{quantikz}
```

Options:
```latex
\slice[style={dashed, blue}]{label}  % styled slice
```

## Gate Groups

Highlight a region of the circuit:

```latex
\begin{quantikz}
  & \gate{H} & \ctrl{1}                                  & \gate{H} & \qw \\
  & \qw      & \targ{} \gategroup[2,steps=1]{inner LCU}  & \qw      & \qw \\
  & \qw      & \gate{A_a}                                & \qw      & \qw
\end{quantikz}
```

Syntax: `\gategroup[num_qubits, steps=num_cols, style={...}, label style={...}]{label}`

```latex
\gategroup[3, steps=5, style={dashed, rounded corners, fill=blue!5, inner sep=4pt},
           background, label style={label position=below, yshift=-0.4cm}]
           {\scriptsize Inner Heisenberg evolution}
```

Key options:
- `steps=N` --- how many columns the group spans
- `style={...}` --- TikZ style for the box
- `background` --- draw behind the gates (not on top)
- `label style={label position=below}` --- where to put the label

## Vertical Wires and Classical Control

For feed-forward based on measurement results:

```latex
\begin{quantikz}
  & \meter{} & \ctrl[vertical wire=c]{1} & \cw \\
  & \qw      & \gate{X}                  & \qw
\end{quantikz}
```

`\ctrl[vertical wire=c]{N}` draws a classical (double) vertical wire from the control to the target.

## Ground / Trash

```latex
\ground{}   % ground symbol (trace out)
\trash{}    % trash symbol (discard)
```

## Ghost (Invisible Gate for Alignment)

When you need to align gates across wires without drawing a visible gate:

```latex
\ghost{U}   % invisible placeholder matching the size of gate U
```

## Phase / Rotation Notation

```latex
\phase{\phi}        % small filled dot with phase label
\phase[style={fill=white}]{\phi}  % open phase dot
```

## Barriers

Draw a barrier (vertical line) across selected wires without a label:

```latex
& \slice{} &   % empty slice = just the line
```

## Advanced: Custom TikZ Overlays

Since quantikz2 is built on TikZ, you can add arbitrary TikZ annotations:

```latex
\begin{quantikz}
  & \gate{U} & \qw & \\
  & \gate{V} & \qw &
\end{quantikz}
% Can add TikZ overlays using the node names generated by quantikz
```

## Advanced: Midcircuit Reset

```latex
& \gate[style={starburst, draw=red, fill=yellow}]{\text{reset}} &
```

Or simply:
```latex
& \rstick{\ket{0}} \setwiretype{q} &
```

## Advanced: Parallel / Grouped Wires with Labels

For showing register structure:

```latex
\begin{quantikz}[wire types={b,b,q,q,q}]
  \lstick{$\Omega$\;($r$)}     & \qwbundle{r}   & \gate{\text{Prep}} & ...   \\
  \lstick{$T_+$\;($r_+$)}      & \qwbundle{r_+}  & \qw               & ...   \\
  \lstick{$\ket{\psi}_S$}       & \qw             & \qw               & ...   \\
  \lstick{$q_\gamma$}           & \qw             & \qw               & ...   \\
  \lstick{$q_\delta$}           & \qw             & \qw               & ...
\end{quantikz}
```

## Complete Examples

### Example 1: Simple Block Encoding (LCU Sandwich)

```latex
\begin{quantikz}
  \lstick{$\ket{0}_a$} & \gate{\text{Prep}} & \ctrl{1}  & \gate{\text{Prep}^\dagger} & \meter{} \\
  \lstick{$\ket{\psi}$} & \qw               & \gate{U_j} & \qw                       & \qw
\end{quantikz}
```

### Example 2: Controlled Hamiltonian Simulation on Bundled Register

```latex
\begin{quantikz}[wire types={b,q}, row sep=0.5cm]
  \lstick{$T$\;($r$ qubits)} & \qwbundle{r} & \gate{\text{Prep}\;\ket{b}} & \ctrl{1} & \gate{\text{Prep}^\dagger} & \qw \\
  \lstick{$\ket{\psi}_S$}     & \qw          & \qw                         & \gate{S_2^{(M)}(\bar{t})} & \qw   & \qw
\end{quantikz}
```

### Example 3: Weak Measurement Pattern

```latex
\begin{quantikz}[row sep=0.4cm]
  \lstick{$\ket{0}_\Omega$} & \gate{\text{Prep}\;\ket{f}}
    & \ctrl{1} & \qw & \ctrl{1} & \gate{\text{QFT}}
    & \ctrl{2} & \qw & \qw \\
  \lstick{$\ket{\psi}_S$} & \qw
    & \gate{e^{+iH\bar{t}}} & \gate{A_a} & \gate{e^{-iH\bar{t}}} & \qw
    & \qw & \qw & \qw \\
  \lstick{$\ket{0}_{q_\gamma}$} & \qw
    & \qw & \qw & \qw & \qw
    & \gate{R_Y} & \ctrl{1} & \qw \\
  \lstick{$\ket{0}_{q_\delta}$} & \qw
    & \qw & \qw & \qw & \qw
    & \qw & \gate{R_Y(\delta)} & \meter{}
\end{quantikz}
```

### Example 4: QSP / GQSP Structure

```latex
\begin{quantikz}[column sep=0.4cm]
  \lstick{$\ket{+}$} & \gate{e^{i\phi_0 Z}} & \ctrl{1} & \gate{e^{i\phi_1 Z}}
    & \ctrl{1} & \ \ldots\ \qw & \ctrl{1} & \gate{e^{i\phi_d Z}} & \meter{} \\
  \lstick{$\ket{\psi}$} & \qw & \gate{U} & \qw
    & \gate{U} & \ \ldots\ \qw & \gate{U} & \qw & \qw
\end{quantikz}
```

### Example 5: Strang Splitting (Second-Order Trotter)

```latex
\begin{quantikz}[column sep=0.35cm]
  \lstick{$\ket{\psi}$}
    & \gate{e^{-iH_1 t/2}} & \gate{e^{-iH_2 t}} & \gate{e^{-iH_1 t/2}}
    & \qw
\end{quantikz}
```

### Example 6: Inner Heisenberg Evolution (Triple Product)

```latex
\begin{quantikz}[wire types={b,q}, row sep=0.5cm, column sep=0.35cm]
  \lstick{$T_+$} & \qwbundle{r_+} & \ctrl{1} & \qw & \ctrl{1} & \qw & \ctrl{1} & \qw \\
  \lstick{$S$} & \qw
    & \gate{S_2^{(M_+)}\!(+\beta\bar{t}')} & \gate{A_a^\dagger}
    & \gate{S_2^{(M_+)}\!(-2\beta\bar{t}')} & \gate{A_a}
    & \gate{S_2^{(M_+)}\!(+\beta\bar{t}')} & \qw
\end{quantikz}
```

### Example 7: Nested LCU (CoherentStep Skeleton)

```latex
\begin{quantikz}[wire types={b,b,q}, row sep=0.5cm, column sep=0.3cm]
  \lstick{$T_-$} & \qwbundle{r_-}
    & \gate{\text{Prep}\,\ket{b_-}} & \ctrl{2} & \qw & \qw & \qw & \qw & \qw
    & \ctrl{2} & \gate{\text{Prep}^\dagger} & \qw \\
  \lstick{$T_+$} & \qwbundle{r_+}
    & \qw & \qw & \gate{\text{Prep}\,\ket{b_+}} & \ctrl{1} & \gate{\text{Prep}^\dagger} & \qw & \qw
    & \qw & \qw & \qw \\
  \lstick{$S$} & \qw
    & \qw & \gate{e^{-iH\bar{t}/\sigma}} & \qw
    & \gate[style={fill=blue!8}]{\;A_a^\dagger(\beta\bar{t}')\,A_a(-\beta\bar{t}')\;}
    & \qw & \qw
    & \gate{e^{+iH\bar{t}/\sigma}} & \qw & \qw & \qw
\end{quantikz}
```

## Tips for Thesis-Quality Circuits

1. **Use `\text{}` inside math gates** for roman-font labels: `\gate{\text{QFT}}`, `\gate{\text{Prep}}`
2. **Subscript register names** in `\lstick`: `\lstick{$\ket{0}_{q_\gamma}$}`
3. **Group related stages** with `\gategroup` or `\slice{}` to guide the reader's eye
4. **Keep circuits horizontal** --- if too wide, break into sub-figures or use `\scalebox{0.85}{...}`
5. **Match thesis notation exactly**: $S_p^{(M)}$, $\tilde{A}_a(\bar\omega)$, $\mathcal{L}$, etc.
6. **For controlled-on-register operations**: draw control from the bundled wire to emphasize that all qubits in the register participate
7. **Ellipsis for repeated structure**: use `\ \ldots\ \qw` between repeated gates (with explicit spaces)
8. **Consistent sizing**: use the same `row sep` and `column sep` across all thesis figures
9. **Gate fill colors** for visual distinction: use very light fills (10-15% saturation) like `fill=blue!10`, `fill=red!8`, `fill=green!8`
10. **For the thesis preamble**, add once:
    ```latex
    \usepackage{tikz}
    \usetikzlibrary{quantikz2}
    ```
    Then each circuit is just `\begin{quantikz}...\end{quantikz}` inside a `figure` environment.
