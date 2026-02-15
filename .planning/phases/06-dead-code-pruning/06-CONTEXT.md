# Phase 6: Dead Code Pruning - Context

**Gathered:** 2026-02-15
**Status:** Ready for planning

<domain>
## Phase Boundary

Remove all dead code from the codebase: ~930 lines of commented-out code blocks across 9 files, ~35 unused functions, and dead structs (LindbladianJumpCaches, LiouvLiouv, non-mutating `oft` wrapper). All 224 existing tests must pass after pruning.

</domain>

<decisions>
## Implementation Decisions

### Preservation policy
- Git history is sufficient -- no need to document removed code elsewhere before deletion
- TODO/FIXME/HACK comments are preserved -- only dead code blocks are removed
- Explanatory comments that describe *why* something was disabled or a decision was made are kept (these are documentation, not dead code)
- Dead structs (LindbladianJumpCaches, LiouvLiouv) are deleted cleanly with no tombstone comments

### Borderline functions
- Claude decides per-function based on reachability and v1.1 relevance
- Functions relevant to upcoming phases 7-11 are preserved even if currently unused
- **Explicit keep list:**
  - `time_oft!()` and `trotter_oft!()` -- only the non-mutating `oft` wrapper is removed
  - `qi_tools` functions -- keep useful quantum information utilities even if not called in simulation paths (pedagogical/user-facing value)
  - `linearmaps` functions -- future project dependency, preserve
  - Log Sobolev functions -- future project dependency, preserve
- Everything else: prune if unreachable from public API or test suite AND not relevant to v1.1 phases

### Claude's Discretion
- Removal ordering strategy (all-at-once vs. file-by-file vs. category-by-category)
- Judgment calls on individual borderline functions not covered by the explicit keep list
- How to verify reachability (static analysis approach)

</decisions>

<specifics>
## Specific Ideas

No specific requirements -- open to standard approaches for dead code identification and removal.

</specifics>

<deferred>
## Deferred Ideas

None -- discussion stayed within phase scope.

</deferred>

---

*Phase: 06-dead-code-pruning*
*Context gathered: 2026-02-15*
