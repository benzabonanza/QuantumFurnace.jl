# Phase 10: API Surface Cleanup - Context

**Gathered:** 2026-02-15
**Status:** Ready for planning

<domain>
## Phase Boundary

Curate the public API so it exposes exactly what users and researchers need as building blocks for pedagogy and research workflows. Internalize implementation details. No new functionality — only reorganize what's exported vs internal.

</domain>

<decisions>
## Implementation Decisions

### Pedagogical boundary
- Export all functions with **physics meaning** — anything a researcher would call to study or construct quantum objects
- This includes: Gibbs states, Pauli gates, Bohr/alpha tools, B-terms, transition functions, Kraus-related physics, QI distance measures, Hamiltonian construction helpers
- Users should be able to create their own Hamiltonians, Kraus jumps, transition functions, and B-terms from exported building blocks
- Internalize **infrastructure/plumbing**: workspaces (OFTCaches, NUFFTPrefactors, KrausScratch, TrajectoryWorkspace), precompute helpers, filename generation, config validation

### Breaking change strategy
- Clean break — remove from export list, no deprecation warnings
- Package is pre-1.0 (0.x), breaking changes in minor versions are expected
- No special versioning or migration path needed

### Internal function naming
- Prefix all internal functions with `_` to signal "don't use this" (but power users can still call via `QuantumFurnace._func_name`)
- Apply `_` prefix to **all** internal functions across the codebase, not just de-exported ones
- Types keep their original names (no `_` prefix on types) — types appear in method signatures and dispatch

### Additional exports (newly exported)
- **All qi_tools.jl functions**: trace_distance_h, trace_distance_nh, trace_norm_h, trace_norm_nh, fidelity, frobenius_norm, is_density_matrix, random_density_matrix, hermitianize!, transform_jumps_to_basis
- **Physics helpers from other files**: create_f, create_f_gauss, create_alpha_gauss, check_alpha_skew_symmetry, trotterize, group_hamiltonian_terms, pauli_string_to_matrix, and other functions with physics meaning
- **load_hamiltonian** stays exported — useful for quick experimentation with pre-built Hamiltonians

### Export block organization
- Reorganize the export block in QuantumFurnace.jl into logical groups with comments (e.g., `# Simulation`, `# QI Tools`, `# Hamiltonian building blocks`)

### errors.jl
- Do NOT export anything from errors.jl for now
- errors.jl is reserved for future physical error rate computation (quadrature errors, Trotter errors across domains), not exception types

### Claude's Discretion
- Exact categorization of each function as physics-meaningful vs internal
- Grouping labels for the organized export block
- Order of exports within each group

</decisions>

<specifics>
## Specific Ideas

- The user wants researchers to be able to build their own quantum objects from exported building blocks — Hamiltonians, Kraus jumps, transition functions, B-terms
- All QI tools (trace distance, fidelity, etc.) should be readily available without qualification
- Internal functions get `_` prefix as a convention signal, but remain accessible via qualified name for power users

</specifics>

<deferred>
## Deferred Ideas

- Physical error rate computation functions in errors.jl — future phase
- Documentation/docstrings for exported API — not in scope for this phase

</deferred>

---

*Phase: 10-api-surface-cleanup*
*Context gathered: 2026-02-15*
