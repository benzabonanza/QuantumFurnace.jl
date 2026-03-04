# /execute — Execute a Plan

Execute a PLAN.md with atomic commits and continuous testing.

**Argument**: `$ARGUMENTS` (path to PLAN.md, phase number, or description to find the plan)

## Setup

1. **Find the plan**: Locate the PLAN.md in `.planning/phases/`
2. **Read it fully**: Understand all tasks, dependencies, must-haves, and risks
3. **Read STATE.md**: Know current project state
4. **Verify preconditions**: Check that dependency tasks are already done

## Execution Loop

For each task in order:

### 1. Implement
- Read the files listed in the task before modifying them
- Follow existing code patterns (check nearby functions for style)
- Write the minimum code needed — no over-engineering

### 2. Test
- Run `julia --project -e 'using Pkg; Pkg.test()'` after each task
- If tests fail: fix the issue before moving on
- If fix requires architectural change: STOP and report to user

### 3. Commit
- Stage only the specific files changed for this task (`git add <file1> <file2> ...`)
- **Never use `git add .` or `git add -A`**
- Commit message format: `<type>(<phase>): <description>`
  - Types: `feat`, `fix`, `refactor`, `test`, `docs`
  - Example: `feat(45): add BiexpFitResult struct`
- Include `Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>`

### 4. Proceed or Halt

**Continue** if:
- Tests pass
- Implementation matches the task description
- No unintended side effects

**Halt and ask user** if:
- Tests fail after reasonable fix attempt
- Task requires changing the plan's approach
- Discovered a risk flag from the plan materializing
- Unsure about a design decision not covered in the plan

## Completion

After all tasks:

1. Run full test suite one final time
2. Update `.planning/STATE.md` with completion status
3. Report summary:
   - Tasks completed: N/N
   - Commits made: list with hashes
   - Tests: pass/fail count
   - Any deviations from plan
