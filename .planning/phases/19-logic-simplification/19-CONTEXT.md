# Phase 19: Logic Simplification - Context

**Gathered:** 2026-02-16
**Status:** Ready for planning

<domain>
## Phase Boundary

Simplify overly complex logic accumulated during v1.2 development. Three targets: flatten the deep trajectory call chain, eliminate redundant jump basis transforms, and simplify result struct hierarchy. No new capabilities — pure simplification of existing code.

</domain>

<decisions>
## Implementation Decisions

### Call chain flattening
- Current chain is too deep: `run_experiment() -> run_trajectories_adaptive() -> run_trajectories() -> _evolve_along_trajectory!() -> step_along_trajectory!()` — 5 levels before actual computation
- Reduce indirection so the path from experiment entry point to trajectory stepping is shorter and more readable
- API signatures, structs, and internal organization can all change freely — no backward compatibility constraint

### Jump basis construction
- Eliminate the `jumps_for_diss_raw` / `transform_jumps_to_basis` / `convert` code block that appears in multiple places
- Fix at the source: when constructing `JumpOp`, set `in_eigenbasis` to the correct basis immediately
- For TrotterDomain: use `trotter.eigvecs` for the basis transform when building the JumpOp
- For other domains: use `hamiltonian.eigvecs` as before
- This is safe because TrotterDomain computation is entirely done in Trotter basis
- The downstream `transform_jumps_to_basis` call and `convert(Vector{JumpOp{Matrix{CT}}}, ...)` become unnecessary

### Result struct simplification
- Current structs (ExperimentResult, TrajectoryResult, ConvergenceData) feel overly complex
- Simplify to distinct result types per simulation method: Lindbladian results, DM simulator results, Trajectory results
- Each type should be clean and self-contained rather than one flexible type with optional fields

### Claude's Discretion
- How to flatten the call chain — which layers to merge, inline, or keep separate
- Internal architecture of simplified result structs
- Whether `step_along_trajectory!` remains its own function or gets merged
- How to handle ConvergenceData in the new struct hierarchy

</decisions>

<specifics>
## Specific Ideas

- User provided exact code for the jump basis fix:
  ```julia
  basis_unitary = (domain isa TrotterDomain) ? trotter.eigvecs : hamiltonian.eigvecs
  jump_op_in_eigenbasis = basis_unitary' * jump_op * basis_unitary
  ```
  This replaces the current pattern of building in energy basis then transforming later.

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 19-logic-simplification*
*Context gathered: 2026-02-16*
