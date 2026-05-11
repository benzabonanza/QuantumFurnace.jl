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

## Conventions: β_phys vs β_alg (qf-6vr / Phase qf-bphys)

The Hamiltonian fixtures store a **rescaled** spectrum in `H.data` / `H.eigvals`
that lives in `[0, 0.45]` for *every* `n` (see `_rescaling_and_shift_factors`
in `src/hamiltonian.jl`). The single user-facing scale is therefore split:

- **`β_phys`** — physical inverse temperature against the un-rescaled
  Hamiltonian `H_phys`. The thing a physicist would type.
- **`β_alg = β_phys · ham.rescaling_factor`** — algorithm-side inverse
  temperature, against the rescaled spectrum that the simulator actually
  sees. Equal to the legacy `cfg.beta`.

For an extensive Hamiltonian (Heisenberg, TFIM) `ham.rescaling_factor`
grows roughly linearly with `n`, so a sweep "at fixed `β_alg` across n"
silently varies `β_phys` by the same factor — and a "fixed `β_phys` across
n" sweep silently varies `β_alg`. Always be explicit about which one you
mean. The qf-6vr layer keeps `cfg.beta` = β_alg for back-compat; the new
fields and helpers are:

| Object | Field / function | Convention |
|---|---|---|
| `HamHam` | `HamHam(raw, β)` (positional) | `β` = β_alg (legacy) |
| `HamHam` | `HamHam(raw; beta_phys=…)` (keyword) | derives `β_alg = β_phys · rescale` |
| `HamHam` | `beta_alg(ham, β_phys)` / `beta_phys(ham, β_alg)` | scalar conversion helpers |
| `Config` | `cfg.beta` (= `beta_alg(cfg)`) | β_alg, required |
| `Config` | `cfg.beta_phys` (= `beta_phys(cfg)`) | β_phys, optional (Union{T, Nothing}) |
| `validate_config!(cfg, ham)` | 2-arg method | enforces `cfg.beta ≈ cfg.beta_phys · ham.rescaling_factor` when β_phys is set |
| `sweep_mixing_times`, `sweep_channel_mixing` | positional `beta_values` | β_alg (legacy) |
| `sweep_mixing_times`, `sweep_channel_mixing` | kwarg `beta_phys_values` | β_phys (qf-6vr); mutually exclusive with the positional |
| Sidecar BSON | `:beta_phys`, `:beta_alg`, `:rescaling_factor` | always emitted |
| Sidecar filename | `beta<β_alg>` vs `betaphys<β_phys>` | swap on mode (no collisions) |
| `fit_scaling(::Vector{<:NamedTuple})` | `beta_kind = :auto` (default) | prefers `:beta_phys`, falls back to `:beta_alg`/`:beta` |
| `ScalingFit` | `beta_kind ∈ {:phys, :alg}` | sets the formula label (`β_phys^y` vs `β_alg^y`) |
| Test constants | `BETA` (= `BETA_ALG`) | β_alg, legacy semantics preserved |
| Test constants | `BETA_PHYS`, `N3_BETA_PHYS` | β_phys derived from each fixture's rescaling_factor |

The canonical β_phys sweep grid is **`{0.25, 0.5, 1.0}`** (decided
2026-05-11, replaces the legacy β_alg grid `{5.0, 10.0, 20.0}`).
Rationale: 0.25 is the smallest β_phys with meaningful thermal contrast
(`S(ρ_β)/log(d)` ≈ 0.80 uniformly across n=3..10 — below that the Gibbs
state is essentially uniform); 1.0 is the practical upper bound (σ =
1/β_alg ≤ 0.04 at n ≥ 8, and the `default_smooth_s(β, σ) = (0.05/σ)²`
rule would force `s` ≈ O(10) at n=11 to preserve the absolute kink width
`σ·√s = 0.05` — a smooth-Metro kernel regime we have not characterised).

All new numerics drivers write β_phys into the sidecar; the harness
derives β_alg per cell. `migrate_bson_beta_phys.jl` annotates legacy
β_alg-keyed sidecars with the new triple so plot scripts can read them
under the same `fit_scaling` contract.

**Caveat on the smooth-Metropolis `s` at large β_alg.** The legacy fixed
`s = 0.25` is held in the current drivers; if a β_phys=1 cell at large n
mixes terribly because the kink-width is now σ·√s ≈ 0.0025 (much narrower
than the legacy 0.05), the choice is between (a) re-enabling
`default_smooth_s` and verifying that s ≈ 25 still gives a physically
sensible γ-rate, (b) holding σ = 0.1 fixed (decoupling the OFT filter
width from β), or (c) reinstating the fixed s = 0.25 value. The choice
hasn't been made yet — make the call after the first full sweep.

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


<!-- BEGIN BEADS INTEGRATION v:1 profile:minimal hash:ca08a54f -->
## Beads Issue Tracker

This project uses **bd (beads)** for issue tracking. Run `bd prime` to see full workflow context and commands.

### Quick Reference

```bash
bd ready              # Find available work
bd show <id>          # View issue details
bd update <id> --claim  # Claim work
bd close <id>         # Complete work
```

### Rules

- Use `bd` for ALL task tracking — do NOT use TodoWrite, TaskCreate, or markdown TODO lists
- Run `bd prime` for detailed command reference and session close protocol
- Use `bd remember` for persistent knowledge — do NOT use MEMORY.md files

## Session Completion

**When ending a work session**, you MUST complete ALL steps below. Work is NOT complete until `git push` succeeds.

**MANDATORY WORKFLOW:**

1. **File issues for remaining work** - Create issues for anything that needs follow-up
2. **Run quality gates** (if code changed) - Tests, linters, builds
3. **Update issue status** - Close finished work, update in-progress items
4. **PUSH TO REMOTE** - This is MANDATORY:
   ```bash
   git pull --rebase
   bd dolt push
   git push
   git status  # MUST show "up to date with origin"
   ```
5. **Clean up** - Clear stashes, prune remote branches
6. **Verify** - All changes committed AND pushed
7. **Hand off** - Provide context for next session

**CRITICAL RULES:**
- Work is NOT complete until `git push` succeeds
- NEVER stop before pushing - that leaves work stranded locally
- NEVER say "ready to push when you are" - YOU must push
- If push fails, resolve and retry until it succeeds
<!-- END BEADS INTEGRATION -->
