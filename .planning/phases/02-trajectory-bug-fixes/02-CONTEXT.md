# Phase 2: Trajectory Bug Fixes - Context

**Gathered:** 2026-02-14
**Status:** Ready for planning

<domain>
## Phase Boundary

Fix trajectory simulation bugs (U_B ordering, normalization, PSD guard, jump sampling faithfulness) and verify CPTP channel correctness. This phase makes the trajectory code correct — cross-validation against DM results is Phase 4.

</domain>

<decisions>
## Implementation Decisions

### Jump sampling fix approach
- Claude determines the correct jump sampling scheme by cross-checking THREE sources: trajectory code, DM code (`jump_contribution!`), and the paper equations (Chen 2023 Theorem III.1, eq. 3.2-3.3)
- Chen 2023 describes the weak measurement scheme explicitly with statevector equations; Chen 2025 adds the coherent term (deterministic coherent evolution via U_B), which is likely already correct in code
- If the jump sampling is structurally wrong (different algorithm, not just a bug), rewrite it to match the paper's construction
- The coherent evolution (U_B) is probably correct; the jump sampling is the likely problem

### U_B ordering
- U_B (coherent unitary) must be applied BEFORE branch selection, matching the DM simulator's ordering
- Claude should verify the DM simulator's ordering and make the trajectory simulator match

### Normalization guard
- Warning + continue (not assertion): emit a brief message like "Normalization violation: sum = X.XXX at step N"
- No full diagnostic dump — keep output minimal
- Check p_nojump + p_res + p_jump_total ≈ 1.0 at each trajectory step

### PSD fallback
- When S matrix is not positive semi-definite, clamp negative eigenvalues to zero (project onto PSD cone)
- Silent fallback — no warning or error, just fix it numerically

### CPTP verification
- Tolerance: 1e-10 (practical, allows small numerical accumulation)
- Single delta value (TEST_DELTA from Phase 1 fixtures)
- Test all three domains: Energy, Time, Trotter
- Claude verifies the correct completeness relation formula from the paper (don't assume K0'K0 + delta*R + U_res'U_res = I is exactly right)
- Separate test file from the bug fix tests

### Bug fix validation
- One dedicated test per TFIX requirement (TFIX-02, TFIX-03, TFIX-04, TFIX-05) for clear traceability
- Single-step tests (one call to `step_along_trajectory!`) — sufficient to isolate each fix
- CPTP verification (TVAL-01) in its own separate test file

### Claude's Discretion
- Internal code organization and refactoring approach for the fixes
- Test helper functions and fixtures beyond what Phase 1 established
- Exact warning message formatting for normalization check
- How to structure the paper cross-check (inline comments, separate analysis doc, etc.)

</decisions>

<specifics>
## Specific Ideas

- Chen 2023 Theorem III.1, equations 3.2 and 3.3 are the primary reference for the weak measurement scheme
- Chen 2025 adds the coherent term B which produces U_B — this is deterministic evolution, separate from the stochastic jump sampling
- The DM simulator's `jump_contribution!` function is the reference implementation for correct channel structure
- Cross-check must verify trajectory code, DM code, AND paper all agree — don't assume any single source is correct

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 02-trajectory-bug-fixes*
*Context gathered: 2026-02-14*
