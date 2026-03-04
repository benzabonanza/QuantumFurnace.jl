# /plan — Research and Plan a Phase

Given a phase number or description, produce a thorough PLAN.md by first researching the codebase.

**Argument**: `$ARGUMENTS` (phase number like "45" or a free-text description)

## Phase 1: Research

Before writing any plan, you MUST read and understand:

1. **Project context**:
   - `.planning/ROADMAP.md` — current milestone structure
   - `.planning/STATE.md` — what's done, what's in progress
   - `.planning/PROJECT.md` — project-level goals and constraints
   - The relevant phase directory in `.planning/phases/` if it exists

2. **Codebase files that will be touched**:
   - Read the actual source files, not just filenames
   - Trace call paths: who calls what, what depends on what
   - Check struct definitions and their constructors
   - Note existing patterns (how similar features were implemented before)

3. **Prior work**:
   - Read completed phase summaries in `.planning/phases/` for similar work
   - Check if functions or patterns already exist that can be reused

4. **Pitfall identification**:
   - Numerical stability: condition numbers, floating-point edge cases
   - Type stability: will new code introduce `Any` in hot paths?
   - Breaking changes: does this modify public API or struct layouts?
   - Test gaps: what existing tests might break?

## Phase 2: Plan Creation

Write the plan to `.planning/phases/<phase-dir>/PLAN.md` where `<phase-dir>` matches the existing directory pattern (e.g., `45-description/`).

### Plan Structure

```markdown
# Phase [N]: [Title]

## Goal
One sentence: what is true when this phase is done.

## Research Summary
Key findings from the research phase. What exists, what needs to change, what patterns to follow.

## Tasks

### Task 1: [Title]
**Files**: list of files to create or modify
**Action**: what specifically to do
**Dependencies**: which tasks must complete first (if any)

### Task 2: ...

## Must-Haves
Observable truths that must hold when done. Each one is verifiable by reading code or running tests.
- [ ] ...

## Risk Flags
What could go wrong and how to mitigate:
- **[Risk]**: [Mitigation]

## Test Strategy
What tests to add or modify. What edge cases to cover.
```

### Quality Bar

- Every task must specify exact files
- Must-haves must be objectively verifiable (not "code is clean" but "all exported functions have tests")
- Tasks should be ordered so each can be atomically committed
- If a task is too large to commit atomically, split it
