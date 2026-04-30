---
name: Thesis .tex files are read-only; drafts must be readable Markdown (not LaTeX dumps)
description: Never Edit/Write .tex files under supplementary-informations/. Drafts go to drafts/<slug>.md as readable Markdown with inline $math$ — not as raw LaTeX in a code fence — so the user can read the .md while typing into the canonical .tex.
type: feedback
originSessionId: 5deed69d-1d32-4683-b7d3-f5c9df2eb1a5
---
The `.tex` files in `supplementary-informations/` (e.g. `1_preliminaries.tex`, `2_methods.tex`, `main.tex`) are **copies** that the user manually syncs from their real thesis project. Treat them as read-only references — never call `Edit` or `Write` on them. All thesis content goes to `drafts/<slug>.md`, formatted as **readable Markdown**.

**Why:** The user maintains the canonical thesis in a separate folder and types every change into the real `.tex` by hand. They keep the rendered `.md` open side-by-side while typing, so a `.md` that is just a `\subsection{...}` LaTeX block stuffed into a code fence is unreadable for that workflow — it forces them to re-read raw LaTeX they could have copy-pasted from the `.tex` itself. Real Markdown (headers, prose, inline `$math$`, display `$$math$$`) renders in any viewer (GitHub, VS Code preview, MacDown) and lets them read the *content* while typing the LaTeX *form* in the other window.

**How to apply** — for every thesis draft or fix patch:

- Output file: `drafts/<slug>.md`. Never `Edit`/`Write` on any `.tex` file in `supplementary-informations/`.
- Headers: use Markdown `#` / `##` / `###`, not `\section{}` / `\subsection{}` (those belong only in the user's `.tex`).
- Prose: plain Markdown with `*italic*` / `**bold**`, not `\textit{}` / `\textbf{}`.
- Math: inline `$...$` and display `$$...$$` (both render in standard Markdown viewers). Keep the math source identical to what should land in `.tex` so the user copies formulas verbatim.
- Display-equation labels: append `<!-- \label{eq:foo} -->` after the closing `$$` so the user sees which label to type.
- Cross-references: parenthetical with bare label, e.g. `(Eq. eq:HS)`, `(§Quantum~states)`, `(Fig. circ:trotter-strang)` — tells the user exactly which `\ref{}` / `\eqref{}` / `\cref{}` to write.
- Citations: inline `[Author Year]` + `[CITE: bibkey]` annotation, with a "Citations" section at the bottom of the draft listing each new bibkey with full bibliographic data plus a fenced ```bibtex block.
- LaTeX-only constructs without a clean Markdown analogue (`\textsc{CNOT}`, `\smallskip`, custom environments, complicated matrix layouts) stay inline as raw LaTeX — the user reads them fluently.
- Open with one short "Insertion target" callout naming the target `.tex` file and the slot to replace; do **not** wrap the entire draft in a single ```latex fence.
- Reading the `.tex` files (with `Read`, `Grep`, `Glob`) to gather context, verify cross-references, count line numbers, etc. is fine and encouraged.
- Compiling `main.tex` with `pdflatex` to sanity-check the user's own changes is fine; just don't modify the source.
- The same rule applies to draft-checker fix patches: emit them as Markdown patches inside the review (e.g. "replace 'foo' with 'bar'"), never as `Edit` calls on the `.tex`.
- Canonical example: `drafts/quantum-circuits-basics.md` is kept in this format as the reference template.
