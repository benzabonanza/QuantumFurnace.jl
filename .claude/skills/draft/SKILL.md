---
name: draft
description: Spawn thesis-writer agent to produce a polished thesis section draft. Use when writing or rewriting any part of the thesis.
argument-hint: "[section description or path to source .md file]"
model: opus
effort: max
context: fork
agent: thesis-writer
allowed-tools: Read Write Glob Grep Bash
---

# /draft — Write a Thesis Section Draft

Produce a polished draft for a given thesis section or topic.

**Argument**: `$ARGUMENTS` — description of what to draft (e.g., "Ding & Chen comparison intro", "numerics chapter: mixing time scaling plots") or a path to a finalized proof file to convert into thesis prose.

## Process

1. **Parse input**: Determine whether `$ARGUMENTS` is:
   - A **free-text prompt** describing what section to write
   - A **file path** to a proof `.md` or other source material to convert into thesis prose

2. **Derive slug**: Create a short kebab-case name for the output file (e.g., `ding-chen-intro`, `mixing-time-scaling`)

3. Before writing, you MUST:
   1. Read relevant pages of `supplementary-informations/thesis.pdf` to match style and notation
   2. Read relevant **reference papers** from `supplementary-informations/` (Chen, Ding, Li & Wang, Lin, Ramkumar & Soleimanifar, Scandi & Alhambra — whichever relate to the topic). Cite them accurately.
   3. Read `.claude-memory/MEMORY.md` and referenced memory files for thesis structure context
   4. Read any source files in `src/` that are relevant to the topic
   5. **Check for existing reviews**: Glob for `drafts/*-review.md` and read any review file related to this section. Treat critical issues and notation fixes from reviews as hard requirements.

4. Write your draft to `drafts/<slug>.md`.

5. **Beads integration**: If a beads issue is associated with this draft (user provides an issue ID, or `bd search "<topic keywords>"` finds one):
   - Mark it in progress at start: `bd update <id> -s in_progress`
   - On completion: `bd note <id> "draft written: drafts/<slug>.md"` and `bd close <id> -r "draft complete"`

6. **Report**: path to the draft file, brief summary of what was written, and any `[CHECK]` flags or open questions.
