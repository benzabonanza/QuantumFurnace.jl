# Phase 11: Allocation Optimization - Context

**Gathered:** 2026-02-15
**Status:** Ready for planning

<domain>
## Phase Boundary

Eliminate unnecessary heap allocations in core simulation hot paths — sparse matrix construction in coherent_bohr, Diagonal wrapper creation in B_time/B_trotter, filter intermediates in jump_workers, and redundant basis transforms in B_trotter multi-jump. All optimizations target existing functions; no new simulation capabilities.

</domain>

<decisions>
## Implementation Decisions

### Sparse matrix strategy (coherent_bohr)
- Replace per-iteration `spzeros` + scatter with index-based accumulation: loop over `bohr_dict` indices directly, accumulate contributions element-wise into B without building any intermediate matrix
- The `f_A_nu_1` full-matrix broadcast (`f(bohr_freqs, nu_2) * jump.in_eigenbasis`) is acceptable — element-wise broadcast into a pre-allocated dense matrix is cheap
- Apply the same index-based approach consistently to both single-jump and multi-jump `B_bohr` variants
- Do NOT precompute/cache per-frequency A_nu matrices — Bohr frequency count explodes at larger system sizes (12 qubits for DM/trajectory)

### Diagonal elimination (B_time / B_trotter)
- Replace `Diagonal(exp.(...))` and `Diagonal(eigvals .^ n)` wrappers with pre-allocated vector buffers + element-wise broadcasting (`.* `)
- No `Diagonal` objects should be created in any loop iteration
- Apply to all 4 variants: B_time single-jump, B_time multi-jump, B_trotter single-jump, B_trotter multi-jump

### Redundant basis transforms (B_trotter)
- Both single-jump and multi-jump `B_trotter` variants currently compute `trotter.eigvecs' * jump.data * trotter.eigvecs` — this is redundant because `JumpOp.in_eigenbasis` already stores the jump in Trotter eigenbasis when using TrotterDomain
- Fix: use `jump.in_eigenbasis` directly instead of recomputing the basis transform
- This applies to both single-jump (`jump_in_trotter = ...`) and multi-jump (`jump_a_trotter = ...` inside the inner loop)

### Filter intermediate (jump_workers)
- The `abs.(filter(w -> w < 1e-12, energy_labels))` allocation may already be eliminated from earlier refactoring
- Researcher should verify whether this hotspot still exists in current code
- If still present: inline the condition in the loop body (skip + abs), keeping current Hermitian branch logic
- Note: `abs()` for energies may be obsolete since energy labels are now returned as a symmetrized grid around 0

### Struct mutation scope
- Prefer function-local pre-allocated buffers over adding new struct fields
- Slightly open to struct changes ONLY if truly needed AND doesn't overcomplicate existing struct design
- No new structs
- Do NOT add domain-specific fields to cross-domain structs (structs have good structure across domains now)
- Keep in mind scale targets: DM/trajectory should run for 12 qubits, Lindbladian constructor for ~8 qubits

### Verification
- Function-level `@allocated` tests for each optimized hot path: B_bohr, B_time, B_trotter, _jump_contribution! (if filter fix needed)
- Existing 224 tests continue to verify correctness
- No entry-point level allocation tests needed

### Claude's Discretion
- Exact implementation of index-based accumulation loop structure
- Choice of broadcasting vs manual loops for element-wise diagonal scaling
- How to structure @allocated test assertions (exact zero vs threshold)
- Whether f_A_nu_1 buffer needs pre-allocation or can rely on in-place broadcast

</decisions>

<specifics>
## Specific Ideas

- "JumpOp.in_eigenbasis should save the jump operator in the Trotter eigenbasis for TrotterDomain, since we should never need a jump in the Hamiltonian eigenbasis in that config" — this is the key insight for the B_trotter fix
- "Number of Bohr frequencies will explode for larger systems" — rules out per-frequency caching strategies
- abs() on energy labels may be obsolete now that truncation returns symmetrized grids

</specifics>

<deferred>
## Deferred Ideas

- **Peak memory estimation / allocation transparency** — ability to see how many large objects (Lindbladians, DM-sized matrices) are simultaneously allocated during simulation. Should become its own phase for memory profiling/budgeting, enabling users to estimate RAM requirements for cluster simulations.

</deferred>

---

*Phase: 11-allocation-optimization*
*Context gathered: 2026-02-15*
