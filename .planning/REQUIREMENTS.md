# Requirements: QuantumFurnace.jl v1.1 Reduce

**Defined:** 2026-02-14
**Core Value:** Correct and efficient classical simulation of Lindbladian-based quantum Gibbs samplers

## v1.1 Requirements

Requirements for the Reduce milestone. Each maps to roadmap phases.

### Dead Code Pruning

- [ ] **PRUNE-01**: Remove all commented-out code blocks (~930 lines across 9 files: qi_tools, energy_domain, trotter_domain, coherent, trajectories, jump_workers, errors, misc_tools, ofts)
- [ ] **PRUNE-02**: Remove all unused active functions (~35 functions across bohr_domain, coherent, trotter_domain, qi_tools, misc_tools, ofts, errors)
- [ ] **PRUNE-03**: Remove dead structs LindbladianJumpCaches, LiouvLiouv, and non-mutating `oft` allocating wrapper. Keep `time_oft!`, `trotter_oft!`, and `OFTCaches` for NUFFT validation

### Struct Simplification

- [ ] **STRUCT-01**: Consolidate 4 config structs (LiouvConfig/GNS, ThermalizeConfig/GNS) into ≤2 structs with a db_type field
- [ ] **STRUCT-02**: Refactor HamHam to eliminate two-step initialization pattern (remove Nothing-typed fields, compute bohr_freqs/gibbs in constructor)
- [ ] **STRUCT-03**: Make TrottTrott immutable and fix num_trotter_steps_per_t0 type from Float64 to Int

### Type Parameterization

- [ ] **TYPE-01**: Parameterize HamHam on element type `{T<:AbstractFloat}`
- [ ] **TYPE-02**: Parameterize LindbladianWorkspace on element type
- [ ] **TYPE-03**: Parameterize config structs on float type

### API Surface Cleanup

- [ ] **API-01**: Remove dead/deprecated exports (non-mutating oft, dead structs)
- [ ] **API-02**: Internalize implementation-detail exports (~18 items: workspaces, precompute helpers, internal dispatch functions)
- [ ] **API-03**: Export trace_distance_h for convergence analysis

### DRY Refactoring

- [ ] **DRY-01**: Extract `hermitianize!` helper to replace 8+ repeated Hermitianization patterns
- [ ] **DRY-02**: Extract shared CPTP channel application function from 3 identical 30-line blocks in jump_workers.jl (K0/residual/Cholesky sequence)
- [ ] **DRY-03**: Extract shared coherent unitary application helper from 3 identical blocks in jump_workers.jl
- [ ] **DRY-04**: Deduplicate Trotter basis transform of jumps (3 identical locations in furnace.jl + trajectories.jl)

### Allocation Optimization

- [ ] **ALLOC-01**: Eliminate sparse matrix allocation in coherent_bohr inner loop (pre-allocate or precompute A_nu matrices)
- [ ] **ALLOC-02**: Eliminate Diagonal allocation in B_time/B_trotter closures (in-place cis computation)
- [ ] **ALLOC-03**: Fix abs.(filter(...)) allocation in Time/Trotter thermalize hot path in jump_workers.jl
- [ ] **ALLOC-04**: Precompute Trotter basis transforms in B_trotter multi-jump variant (avoid recomputation in inner loop)

## Future Requirements

Deferred from v1.1 — tracked for future milestones.

- **PERF-01**: Multi-threaded trajectory sampling with shared precomputed data
- **DOC-01**: API docs via Documenter.jl
- **DOC-02**: Theory tutorials via Literate.jl

## Out of Scope

| Feature | Reason |
|---------|--------|
| New simulation capabilities | v1.1 is cleanup only — no new physics |
| linearmaps/log-sobolev file deletion | Kept intentionally for future milestones |
| Float32 simulation testing | Type params enable it but testing F32 paths is future work |
| Performance benchmarking suite | Allocation fixes are correctness-driven, formal benchmarks are future |
| Refactoring test code | Tests are the validation gate, not the target |
| errors.jl rewrite | Keep stubs for future; broken functions pruned but file structure preserved |

## Traceability

Which phases cover which requirements. Updated during roadmap creation.

| Requirement | Phase | Status |
|-------------|-------|--------|
| PRUNE-01 | — | Pending |
| PRUNE-02 | — | Pending |
| PRUNE-03 | — | Pending |
| STRUCT-01 | — | Pending |
| STRUCT-02 | — | Pending |
| STRUCT-03 | — | Pending |
| TYPE-01 | — | Pending |
| TYPE-02 | — | Pending |
| TYPE-03 | — | Pending |
| API-01 | — | Pending |
| API-02 | — | Pending |
| API-03 | — | Pending |
| DRY-01 | — | Pending |
| DRY-02 | — | Pending |
| DRY-03 | — | Pending |
| DRY-04 | — | Pending |
| ALLOC-01 | — | Pending |
| ALLOC-02 | — | Pending |
| ALLOC-03 | — | Pending |
| ALLOC-04 | — | Pending |

**Coverage:**
- v1.1 requirements: 20 total
- Mapped to phases: 0
- Unmapped: 20 ⚠️

---
*Requirements defined: 2026-02-14*
*Last updated: 2026-02-14 after scoping with user*
