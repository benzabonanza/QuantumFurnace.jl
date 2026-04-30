---
name: prove
description: Orchestrate a rigorous proof through iterative prover/reviewer cycles with ledger tracking. Use when a mathematical or physics claim needs formal proof.
argument-hint: "[statement to prove or path to statement.md]"
model: opus
effort: max
allowed-tools: Agent Read Write Edit Glob Grep Bash
---

# /prove — Formal Proof Orchestration

Orchestrate a rigorous proof of a mathematical/physics statement through iterative prover/reviewer cycles with ledger tracking.

**Argument**: `$ARGUMENTS` — either the statement to prove (free text) or a path to an existing `statement.md` file.

## Setup

1. **Parse input**: If `$ARGUMENTS` is a file path, read it. Otherwise, treat it as the statement text.
2. **Create slug**: derive a short kebab-case topic slug from the statement (e.g., `kms-mixing-bound`, `spectral-gap-scaling`)
3. **Create directory**: `proofs/<slug>/`
4. **Write statement file**: `proofs/<slug>/statement.md` with the precise mathematical statement. Include:
   - The claim in formal mathematical language
   - Context: what objects are involved, what space they live in
   - Any known constraints or hints from the user
5. **Initialize ledger**: `proofs/<slug>/ledger.md` with header only

## Iteration Loop

Run up to **10 iterations**. For iteration N:

### Step 1: Spawn Prover Agent

Spawn the `prover` agent with this prompt:

> Read the statement in `proofs/<slug>/statement.md`.
> [If N > 1]: Read the ledger in `proofs/<slug>/ledger.md` to learn from previous attempts — do NOT repeat failed approaches.
> Write your proof to `proofs/<slug>/proof-v<N>.md`.

Wait for the prover to finish. Read its output file.

### Step 2: Spawn Proof Reviewer Agent

Spawn the `proof-reviewer` agent with this prompt:

> Read the statement in `proofs/<slug>/statement.md`.
> Read the proof attempt in `proofs/<slug>/proof-v<N>.md`.
> Write your review to `proofs/<slug>/review-v<N>.md`.

Wait for the reviewer to finish. Read its output file.

### Step 3: Update Ledger

Read `review-v<N>.md`. Append to `proofs/<slug>/ledger.md`:

```markdown
## Iteration N
**Approach**: [1-sentence summary of proof strategy used]
**Verdict**: [PASS / PARTIAL / FAIL]
**Greenlit steps**: [list of step numbers confirmed correct]
**Failed steps**: [list with 1-line explanation each]
**Key insight for next attempt**: [what the next prover should do differently]
```

Keep entries concise — the ledger is a compressed signal for the next prover, NOT a copy of the proof or review.

### Step 4: Check Termination

- **PASS**: Report success, break loop
- **Iteration = 10**: Report best result, break loop
- **PARTIAL with progress**: Continue — the proof is improving
- **FAIL with no progress for 2 consecutive iterations**: Ask the user whether to continue, reformulate the statement, or abandon

## Completion

Report to user:
- **Final verdict**: PASS / best PARTIAL / FAIL
- **Proof location**: path to the final proof file
- **Journey summary**: approaches tried, key difficulties encountered, iterations used
- **If PARTIAL/FAIL**: what remains unresolved, suggested reformulations

## Beads Integration

If a beads issue is associated with this proof (user provides an issue ID, or search `bd search "<statement keywords>"` finds one):
1. Mark it in progress at start: `bd update <id> -s in_progress`
2. After each iteration, note the verdict: `bd note <id> "iteration N: VERDICT — <1-line summary>"`
3. On completion:
   - **PASS**: `bd close <id> -r "proved in N iterations: proofs/<slug>/proof-vN.md"`
   - **FAIL**: `bd note <id> "FAIL after N iterations — <what remains unresolved>"`

## Edge Cases

- If the prover says the statement is **false**: report this immediately with the counterexample/argument. Ask user whether to prove a corrected version.
- If the prover and reviewer **disagree on basic definitions**: flag the ambiguity and ask the user to clarify the statement.
- If the statement is **trivial** (proved in iteration 1 with PASS): just report success, no need to celebrate.
