# QuantumFurnace.jl

Julia package for quantum thermal simulation — the computational backbone of a PhD thesis on Gibbs sampling via KMS detailed-balance Lindbladians.

## Quick Reference

```
julia --project -e 'using Pkg; Pkg.test()'   # run tests
julia --project scripts/<name>.jl             # run a script
julia --project                               # REPL, then: using QuantumFurnace
```

This is a **Julia** project, not Python. Ignore any parent-level CLAUDE.md references to Python/pip/pytest.

## Project Structure

```
src/                           Julia source (QuantumFurnace module)
  QuantumFurnace.jl              Module definition, includes, exports
  fitting.jl                     Curve fitting (single-exp, bi-exp via LsqFit.jl)
  mixing.jl                      Mixing time estimation (post-processing)
  errors.jl                      Error analysis, floor models
test/                          Test suite
  runtests.jl                    Entry point (includes test_*.jl files)
scripts/                       Standalone analysis and diagnostic scripts
supplementary-informations/    Thesis PDF and reference materials
proofs/                        Formal proof artifacts (statement, proof, review, ledger per topic)
drafts/                        Thesis section drafts
  circuits/                      Quantum circuit figures (quantikz2 LaTeX)
.claude/agents/                Agent definitions (8 agents)
.claude/skills/                Skill definitions (each skill is a directory with SKILL.md)
.claude/rules/                 Path-scoped rules (load based on active file paths)
.claude-memory/                Persistent memory (symlinked into Claude's memory system)
```

## Thesis Context

PhD thesis on quantum Gibbs sampling. Remaining goals:

1. **Numerics chapter** — plots cross-checking all analytical claims (see memory: `thesis_numerics_plan.md`)
2. **Ding & Chen comparison** — implement their KMS DB Lindbladian variant, compare mixing times and Kossakowski matrices, write 5–10 page chapter
3. **Quantum advantage literature review** — what's efficient classically vs quantum for Gibbs sampling, identify parameter regimes where quantum advantage may emerge
4. **Polish** — rigorous proofs, clean writing, complete all sections

## Path-Scoped Rules

Context-specific rules load automatically based on which files you're working with:

- `.claude/rules/julia-code.md` — coding conventions, architecture, atomic commits, physics checks (loads for `src/`, `test/`)
- `.claude/rules/scripts.md` — standalone script conventions (loads for `scripts/`)
- `.claude/rules/thesis-writing.md` — notation, reference papers, proof/review style (loads for `drafts/`, `proofs/`, `supplementary-informations/`)

## Agents

All agents run with `model: opus` and `effort: max`.

**Producers** (have skills preloaded for planning/orchestration):

| Agent | Purpose | Produces | Skills |
|-------|---------|----------|--------|
| `prover` | Formal mathematical/physics proofs | `proofs/<topic>/proof-vN.md` | — |
| `thesis-writer` | Thesis section drafts in existing style | `drafts/<section>.md` | `plan` |
| `sci-coder` | Standalone Julia implementations | `scripts/scratch_<name>.jl` | `plan` |
| `code-integrator` | Merge new code into codebase | Edits to `src/`, `test/` | `plan` |

**Reviewers** (no skills — stay independent and unbiased):

| Agent | Purpose | Produces |
|-------|---------|----------|
| `proof-reviewer` | Skeptical proof verification | `proofs/<topic>/review-vN.md` |
| `draft-checker` | Rigorous thesis draft/section review | `drafts/<slug>-review.md` |
| `code-verifier` | Test and verify code | New tests, test reports |
| `physics-checker` | Validate physical/mathematical sense | Analysis reports |

## Skills

All skills run with `model: opus` and `effort: max`.

| Skill | Usage | What it does | Execution |
|-------|-------|--------------|-----------|
| `/prove` | `/prove "statement"` | Prover↔reviewer loop (max 10 iterations) with ledger | inline |
| `/draft` | `/draft "section description"` | Produce a polished thesis section draft | forked (`thesis-writer`) |
| `/draft-check` | `/draft-check "path or pages"` | Review draft for errors, notation, citations, exposition | forked (`draft-checker`) |
| `/plan` | `/plan <description>` | Research and plan an implementation phase | inline |
| `/execute` | `/execute <plan>` | Execute a plan with atomic commits | inline |
| `/circuit` | `/circuit "description"` | Generate quantum circuit diagram (quantikz2 LaTeX) | inline |

## Beads Issue Tracking

This project uses [beads](https://github.com/steveyegge/beads) (`bd`) for issue tracking, backed by Dolt. Dolt is the only supported backend in `bd` v0.62 (the `no-db: true` config option is silently ignored).

### Sandbox Dolt fix

This Docker sandbox has `sleep infinity` as PID 1 with no `tini`/`dumb-init` subreaper, so any orphaned dolt subprocess becomes an unreapable zombie. The fix:

1. **Vendored static `tini` binary** at `.beads/bin/tini` (~550KB, in the persistent workspace volume so it survives container rebuilds).
2. **Wrapper script** replaces the real `bd` binary at `/usr/local/share/npm-global/lib/node_modules/@beads/bd/bin/bd` and invokes `tini -s -- bd.real "$@"`. The original binary is preserved as `bd.real` alongside.
3. **Idempotent installer** at `.beads/bin/install-bd-wrapper.sh` runs from a `SessionStart` hook in `.claude/settings.json`, so the wrapper is reinstalled automatically after any container rebuild.
4. `.beads/config.yaml` has `dolt.auto-start: true`, so the per-project dolt sql-server comes up on the first `bd` call in a fresh session.

If `bd` ever stops working: `bash .beads/bin/install-bd-wrapper.sh` re-applies the wrapper, and `bd init --prefix qf --from-jsonl --force` rebuilds the dolt database from `.beads/issues.jsonl` (the git-tracked source of truth).

A small trickle of dolt zombies (~0.5/call) still appears because tini exits with `bd.real` and can't catch grandchildren that double-fork into PID 1's tree. Cosmetic, not functional — clears on container restart, well below PID exhaustion limits for any realistic session.

### Key Commands

```
bd list                                    # show all issues (tree view)
bd show <id>                               # issue details
bd create "title" -t task -p 2             # create task (P0=critical, P4=low)
bd create "title" -t epic                  # create epic
bd create "title" --parent <epic-id>       # create child of epic
bd update <id> -s in_progress              # mark in progress
bd update <id> -s open                     # revert to open
bd close <id> -r "reason"                  # close with reason
bd note <id> "text"                        # append a note
bd search "query"                          # text search
bd ready                                   # show unblocked work
```

### Integration Rules

- **Automatic hooks**: `TaskCreated` and `TaskCompleted` hooks in `settings.json` auto-sync Claude Code tasks to beads issues (create on task creation, close on completion)
- **`/plan`** creates a beads epic with child issues for each planned task
- **`/execute`** marks issues `in_progress` when starting each task and closes them on completion
- **Code agents** (`sci-coder`, `code-integrator`, `code-verifier`) update beads issues directly when given an issue ID
- **Orchestrating skills** (`/prove`, `/draft`, `/draft-check`) note results on associated beads issues after their agents finish
- When starting any work associated with a beads issue, mark it `in_progress`
- When completing work, close the issue with a brief reason
- Never create duplicate issues — check `bd list` or `bd search` first

## Docker Persistence

This workspace runs inside a Docker container that may be rebuilt. Everything under the workspace directory persists (mounted volume). The memory system at `.claude-memory/` is symlinked into Claude's internal memory path via a `SessionStart` hook in `.claude/settings.json`. No action needed — it reconnects automatically on each session.
