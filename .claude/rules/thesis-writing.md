---
paths:
  - "drafts/**"
  - "proofs/**"
  - "supplementary-informations/**"
---

# Thesis Writing Rules

## Context

PhD thesis on quantum Gibbs sampling via KMS detailed-balance Lindbladians.

## Reference Papers (in `supplementary-informations/`)

All knowledge-heavy work (proofs, drafts, reviews) must consult these as needed:

- Chen et al. 2023 — Quantum Thermal State Preparation
- Chen et al. 2025 — Efficient noncommutative quantum Gibbs sampler
- Ding et al. 2024 — Efficient quantum Gibbs samplers with KMS
- Ding et al. 2024 — Polynomial-Time Preparation of Low-Temperature Gibbs states
- Ding et al. 2025 — End-to-End Efficient Quantum Thermal and Ground State Preparation
- Li & Wang — Simulating Markovian open quantum systems (higher order series expansion)
- Lin 2025 — Dissipative Preparation of Many-Body Quantum States
- Ramkumar & Soleimanifar 2024 — Mixing time for random sparse Hamiltonians
- Scandi & Alhambra 2025 — Thermalization and KMS detailed balance

The thesis PDF is at `supplementary-informations/thesis.pdf`.

## Output discipline — .tex files are read-only, drafts are pure Markdown

Files under `supplementary-informations/` (`1_preliminaries.tex`, `2_methods.tex`, `main.tex`, …) are **copies** the user manually syncs from their canonical thesis project. **Never** call `Edit` or `Write` on them. Read them freely for context, but all generated content goes to `drafts/<slug>.md`. The same applies to review/fix patches: emit them as `.md`, do not patch the `.tex` in place. Compiling `main.tex` to sanity-check is fine; modifying its sources is not.

**Format the draft as readable Markdown, not as a LaTeX dump in a code block.** The user reads the `.md` side-by-side with the `.tex` they type into manually, so the draft must render legibly in a Markdown viewer:

- `#` / `##` / `###` Markdown headers for structure (not `\section{}` / `\subsection{}`).
- Prose in plain Markdown with `*italic*` and `**bold**` (not `\textit{}` / `\textbf{}`).
- Inline math with `$...$` and display math with `$$...$$` — both render in standard Markdown viewers (GitHub, VS Code preview, MacDown, etc.). Keep the math source identical to what should appear in the `.tex` (so the user can copy a formula verbatim).
- Display equations that need a LaTeX label: write the label as an HTML comment after the block, e.g. `$$ X = \cdots $$ <!-- \label{eq:paulis} -->`.
- Citations: write inline as `[Author Year]` followed by a footnote-style `[CITE: bibkey]` annotation, and list the bibtex keys with full author/title/journal in a "Citations" section at the bottom of the draft. New bibtex entries go in a fenced ```bibtex block in that section.
- Cross-references: write inline as `(see §Quantum~states)`, `(Eq. eq:HS)`, `(Fig. circ:trotter-strang)` — leave the bare label so the user knows exactly which `\ref{}` / `\eqref{}` / `\cref{}` to type.
- LaTeX-only constructs that have no clean Markdown analogue (`\textsc{CNOT}`, `\smallskip`, custom environments, complex matrix layouts) stay inline as raw LaTeX — the user reads them fluently.
- Add an opening "Insertion target" line naming the target `.tex` file and the slot to replace, but do not bracket the entire draft in a `latex` fence.

The goal: the user opens the rendered `.md` in a viewer, sees real headings and real math, and types the corresponding LaTeX into their canonical thesis with as little visual translation as possible.

## Style Requirements

- Match existing thesis notation: $\mathcal{L}$ for Lindbladian, $\beta$ for inverse temperature, etc.
- Use LaTeX notation (`$...$` inline, `$$...$$` display)
- Rigorous claims need proofs or citations — never assert a mathematical fact without backing
- Flag uncertainties with `[CHECK]` rather than asserting
- Every paragraph should advance the argument; cut padding
- Cross-reference style: `\ref{...}`, citation style: `\cite{key}`

## Proofs

- Minimize informal reasoning; every step must be justifiable
- State assumptions explicitly — do not hide hypotheses
- Be precise about quantifiers: $\forall$, $\exists$, domains, ranges
- If you cannot prove a statement: write the strongest partial result, clearly mark the gap

## Reviews

- Be specific: "notation is inconsistent" is useless; "$\mathcal{L}$ is used for both X and Y in paragraphs 2 and 5" is useful
- Verify before flagging — don't flag false positives
- Separate severity levels: typo != logical error != fatal flaw
- Every criticism must come with a concrete fix suggestion
