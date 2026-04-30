---
name: draft-check
description: Review a thesis draft or section for logical errors, notation inconsistencies, missing citations, and exposition quality. Use after writing a draft or to review parts of the thesis.
argument-hint: "[path to .md draft, thesis page range like 'pages 12-18', or free-text excerpt]"
model: opus
effort: max
context: fork
agent: draft-checker
allowed-tools: Read Write Glob Grep Bash
---

# /draft-check — Thesis Draft Review

Rigorously review a thesis draft or existing thesis section.

**Argument**: `$ARGUMENTS` — one of:
- A **file path** to a draft `.md` file (e.g., `drafts/mixing-time-scaling.md`)
- A **thesis page range** (e.g., `pages 12-18` or `Chapter 3`)
- A **free-text description** of what to review, possibly with inline content

## Process

1. **Parse input**: Determine the review target from `$ARGUMENTS`:
   - If it looks like a file path -> read that file
   - If it mentions pages/chapters of the thesis -> read those pages from `supplementary-informations/thesis.pdf`
   - If it's free text -> treat it as the material to review

2. **Derive slug**: Create a short kebab-case name for the review output (e.g., `mixing-time-scaling-review`, `ch3-kms-review`)

3. Before reviewing, you MUST:
   1. Read enough of `supplementary-informations/thesis.pdf` to understand notation conventions, thesis structure, and adjacent sections (at minimum: table of contents, preliminaries/notation section, and sections immediately before/after the reviewed material)
   2. Read relevant **reference papers** from `supplementary-informations/` — verify every citation, check for missing references
   3. Read `.claude-memory/MEMORY.md` and referenced memory files for thesis context and known findings
   4. If the material discusses code or numerics: read relevant files in `src/` and `scripts/`
   5. If there are related drafts in `drafts/`, skim them for consistency

4. Write your review to `drafts/<slug>-review.md`.

5. **Beads integration**: If a beads issue is associated with this review (user provides an issue ID, or `bd search "<topic keywords>"` finds one):
   - Note the review result: `bd note <id> "review: <OVERALL_RATING> — <critical count> critical, <notation count> notation, <citation count> citation issues. Review at drafts/<slug>-review.md"`

6. **Report**:
   - Path to the review file
   - The overall assessment (1-2 sentences)
   - Count of critical issues, notation issues, missing citations, and suggestions
   - The summary scorecard table
   - Any items flagged as particularly well-done (so the user knows what to preserve)
