---
name: lit-review
description: Spawn lit-reviewer agent to distill a paper into a concise, regime-tagged review for the thesis Review chapter on quantum vs classical Gibbs sampling. Use when adding a paper to the literature corpus.
argument-hint: "[path to PDF in supplementary-informations/, arXiv ID/URL, or paper title]"
model: opus
effort: high
context: fork
agent: lit-reviewer
allowed-tools: Read Write Glob Grep Bash WebFetch WebSearch
---

# /lit-review — Literature Review for Thesis Review Chapter

Distill one paper into a structured, regime-tagged Markdown review focused on the **quantum-vs-classical Gibbs sampling comparison**. Reviews are short by design and accumulate under `drafts/literature-review/` so a later synthesis pass can group them by regime (high-T, low-T, commuting, noncommuting, fermionic, bosonic, …).

**Argument**: `$ARGUMENTS` — one of:
- A **path** to a paper PDF (e.g., `supplementary-informations/Chen et al. - 2023 - Quantum Thermal State Preparation.pdf`)
- An **arXiv ID or URL** (e.g., `2303.18224` or `https://arxiv.org/abs/2303.18224`)
- A **paper title or DOI** — agent will locate it via web search
- Any of the above plus a **focus hint** (e.g., `"…  focus: mixing-time bound only"`)

## Process

1. **Parse input**: Determine whether `$ARGUMENTS` is a local PDF path, an arXiv reference, or a free-text title. Strip any trailing focus hint and pass it through to the agent.

2. **Derive slug**: kebab-case `lastname-year-shortkey` (e.g., `chen-2023-thermal-prep`, `ramkumar-2024-sparse-mixing`, `ding-2024-low-temp`). The slug determines the output filename.

3. **Check for duplicates**: `ls drafts/literature-review/` (create the directory if missing). If a review with the same slug already exists, ask the user whether to replace, version (`-v2`), or abort.

4. **Spawn the agent** with the paper reference and focus hint. The agent will:
   - Read the paper (skim, not transcribe)
   - Read minimal thesis/memory context for regime alignment
   - Skim adjacent existing reviews for cross-links
   - Write the review to `drafts/literature-review/<slug>.md` with mandatory regime frontmatter

5. **Beads integration**: If a beads issue is associated with this review (user provides an issue ID, or `bd search "literature review"` finds an open epic):
   - Note the result: `bd note <id> "lit review: <slug> — <regime cell>. drafts/literature-review/<slug>.md"`

6. **Report**:
   - Path to the review file
   - The frontmatter regime cell (temperature × commutativity × locality × particle-statistics)
   - The one-sentence takeaway
   - `key-scaling` line from frontmatter
   - Any `[CHECK]` flags
   - Closely related existing reviews (so the user knows what to cross-read)
