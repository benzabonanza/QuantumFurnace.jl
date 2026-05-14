---
name: feedback-thesis-scope-finite-dim-only
description: "Thesis scope is finite-dim systems only. Infinite-dim / bosonic / continuous-variable Gibbs samplers should be mentioned but never prioritised in syntheses, frontier rankings, or chapter focus."
metadata: 
  node_type: memory
  type: feedback
  originSessionId: ecc5ecf9-b993-4095-b2ea-a932eae84f21
---

The thesis works exclusively with **finite-dimensional Hilbert spaces** — spin systems, fermions on a lattice (truncated to bounded local Fock dim), qubits / qudits. Infinite-dim / continuous-variable / bosonic Gibbs samplers (e.g. Becker 2026 single-mode photon-number, Bose–Hubbard companion, oscillator chains) are outside scope.

**Why:** the whole machinery (KMS-DB Lindbladian construction, the energy / Bohr / Trotter quadrature ladders, Krylov spectral analysis, the algorithm correctness proofs, the numerics chapter) is built for and tested on finite-dim. Infinite-dim requires a different functional-analytic framework (domain conditions for unbounded generators, photon-number truncation error, etc.) that the thesis does not develop.

**How to apply:**
- In any literature synthesis or frontier ranking, infinite-dim cells (Becker 2026, bosonic, continuous-variable) get a **single honest mention** that the cell exists, but are **never placed in a top-N "where quantum could first beat classical" ranking**. The ranking should focus on cells where the thesis machinery can actually compete.
- When the chapter writer surfaces "open frontiers", infinite-dim should appear in a separate aside / footnote, not in the main spine.
- If the user asks for cell-by-cell results, infinite-dim can stay in the cell table (with `[OUT-OF-SCOPE: infinite-dim]` tag) so the chapter explicitly explains why it is not pursued — silently dropping it would look like an oversight.
- Same rule applies to the `[OPEN-FRONTIER]` and `[VERDICT-UNCLEAR]` categories: an infinite-dim contact gets tagged but does not count toward "the frontier the thesis is pushing".

User feedback delivered 2026-05-14 after seeing Becker 2026 (bosonic / infinite-dim) ranked at position #4 in the [[lit-synthesis-all]] top-5 frontier — quote: "infinite dim is of less importance here, the whole thesis only worked and meant to work with finite dim systems. However of course infinite should be mentioned as its also out there."
