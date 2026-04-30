---
name: TODO — bring R_b reflection into 2_methods.tex GQSP section
description: Methods chapter's GQSP section mentions U_B alone as the walk; the prelim figure correctly shows W = R_b U_B. Bring the reflection over (circuit + text) to make Chapter 5 consistent with the preliminaries.
type: project
---

# TODO: Bring R_b (qubitization reflection) over to `2_methods.tex`

## What

The preliminaries subsection `sec:prelim-QSP` commits (correctly) to GQSP applied to the qubitization walk $W = R_b U_B$ with $R_b = 2\ket{0^b}\bra{0^b} - I$ — this is what the figure `circ:gqsp` at `1_preliminaries.tex:1068-1113` draws, with $\Lambda(U_B)$ and $\Lambda(R_b)$ as two separate open-controlled gates per layer.

But `2_methods.tex:1731` currently says "Since $U_{B_a}$ is a unitary that block-encodes a Hermitian operator … its eigenvalues come in conjugate pairs $\ee^{\pm\ii\theta}$ where $\sin\theta = \lambda$". This phrasing treats $U_{B_a}$ itself as the walk, without mention of the $R_b$ reflection. The two are only consistent if either (a) $U_{B_a}$'s construction already folds in a reflection, or (b) the sentence is using loose language — in which case the explicit walk $W = R_b U_{B_a}$ should be stated.

**Why:** the preliminaries and methods chapter must agree on which unitary GQSP is being applied to, otherwise readers get a convention mismatch. The reflection also shows up in the cost accounting (one $R_b$ per layer = $d$ extra reflections per coherent step).

**How to apply:**
- When polishing Chapter 5's QSP section (around `2_methods.tex:1730-1738`), bring the $R_b$ reflection into both the text (walk operator $W = R_b U_{B_a}$) and the still-to-be-drawn Algorithm 1a circuit (`\todo{add circuit figures}` at `2_methods.tex:1935`).
- Also resolve the $\sin\theta = \lambda$ vs $\cos\theta = \lambda$ convention — standard qubitization gives $\cos\theta = \lambda$; the thesis currently says $\sin\theta = \lambda$. Either swap to $\cos$ or document the specific convention that yields $\sin$ (there are variants using a different reflection).
- Keep the preliminaries figure `circ:gqsp` as the canonical GQSP circuit; Chapter 5 should explicitly reference it rather than redefining.
