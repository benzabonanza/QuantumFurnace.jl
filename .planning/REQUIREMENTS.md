# Requirements: QuantumFurnace.jl v1.1 Reduce

**Defined:** 2026-02-14
**Core Value:** Correct and efficient classical simulation of Lindbladian-based quantum Gibbs samplers

## v1.1 Requirements

### Dead Code Pruning

- [ ] **PRUNE-01**: Dead functions in `bohr_domain.jl` removed (coherent_bohr_gauss, transition_bohr_gauss_vectorized, transition_bohr_gauss_gibbsed_vectorized, thermalize_bohr_gauss_vectorized, B_nu_gauss, R_nu_gauss, check_alpha_skew_symmetry, find_all_nu1s_to_nu2)
- [ ] **PRUNE-02**: Dead functions in `coherent.jl` removed (coherent_term_time, coherent_term_trotter, coherent_term_time_metro_exact, coherent_term_timedomain_integrated_gauss, coherent_term_time_integrated_metro, coherent_term_time_integrated_eh, coherent_term_time_integrated_eh_b, check_B_gauss, check_B_metro)
- [ ] **PRUNE-03**: Dead functions in `trotter_domain.jl` removed (trotterize, trotter_diag, trotter2_diag, trotter2_t0_multiple)
- [ ] **PRUNE-04**: Broken/dead functions in `qi_tools.jl` removed (are_we_tp, frobenius_norm if unused)
- [ ] **PRUNE-05**: Deprecated OFT functions removed (time_oft!, trotter_oft! in `ofts.jl`)
- [ ] **PRUNE-06**: Dead code in `errors.jl` removed (compute_errors, stubs that would crash)
- [ ] **PRUNE-07**: Entirely-dead files cleaned up (linearmaps_liouv.jl commented-out content, log_sobolev_manopt.jl empty file)
- [ ] **PRUNE-08**: Dead structs removed (LiouvLiouv, LindbladianJumpCaches, OFTCaches if unused)
- [ ] **PRUNE-09**: All 224 existing tests still pass after pruning

### Struct Simplification

- [ ] **STRUCT-01**: Unused struct definitions removed from `structs.jl`
- [ ] **STRUCT-02**: Remaining structs reviewed for unnecessary fields or complexity
- [ ] **STRUCT-03**: TrottTrott mutability assessed — made immutable if no fields are mutated after construction
- [ ] **STRUCT-04**: All 224 tests still pass after struct changes

### Type Parameterization

- [ ] **TYPE-01**: HamHam parameterized on float precision (ComplexF64 → Complex{T} where T<:AbstractFloat)
- [ ] **TYPE-02**: Result structs (HotAlgorithmResults, HotSpectralResults) parameterized on precision
- [ ] **TYPE-03**: Workspace structs (LindbladWorkspace, LindbladianWorkspace) parameterized on precision where not already
- [ ] **TYPE-04**: Type consistency enforced — no mixed Float64/ComplexF64 hardcoding in parameterized structs
- [ ] **TYPE-05**: Default behavior unchanged (Float64 precision by default, no user-visible API change)
- [ ] **TYPE-06**: All 224 tests still pass with type changes

### API Surface

- [ ] **API-01**: Exports audited — internal helpers that shouldn't be public are unexported
- [ ] **API-02**: Useful building blocks for researchers are exported (trace_distance, fidelity, quantum info tools)
- [ ] **API-03**: Exports list organized by category (types, simulation, building blocks, utilities)
- [ ] **API-04**: No exported symbol references dead code
- [ ] **API-05**: All 224 tests still pass after export changes (test imports updated if needed)

### Redundancy Removal

- [ ] **REDUN-01**: Redundant normalization checks across simulation layers identified and removed
- [ ] **REDUN-02**: Duplicate validation logic consolidated (e.g., config validation not repeated unnecessarily)
- [ ] **REDUN-03**: All 224 tests still pass after redundancy removal

### Allocation Optimization

- [ ] **ALLOC-01**: Per-step allocation in jump_contribution! for Time/TrotterDomain thermalization eliminated (abs.(filter(...)) in jump_workers.jl)
- [ ] **ALLOC-02**: spzeros allocation inside Bohr frequency loops replaced with pre-allocated or accumulated approach
- [ ] **ALLOC-03**: Any other per-step allocations in core simulation paths identified and eliminated
- [ ] **ALLOC-04**: All 224 tests still pass after allocation changes

## Future Requirements

Deferred from v1.1 — tracked for future milestones.

- **PERF-01**: Multi-threaded trajectory sampling with shared precomputed data
- **DOC-01**: API docs via Documenter.jl
- **DOC-02**: Theory tutorials via Literate.jl

## Out of Scope

| Feature | Reason |
|---------|--------|
| New simulation capabilities | v1.1 is cleanup only — no new physics |
| linearmaps/log-sobolev/errors file deletion | Kept intentionally for future milestones |
| Float32 simulation testing | Type params enable it but testing F32 paths is future work |
| Performance benchmarking suite | Allocation fixes are correctness-driven, formal benchmarks are future |
| Refactoring test code | Tests are the validation gate, not the target |

## Traceability

| Requirement | Phase | Status |
|-------------|-------|--------|
| PRUNE-01 | — | Pending |
| PRUNE-02 | — | Pending |
| PRUNE-03 | — | Pending |
| PRUNE-04 | — | Pending |
| PRUNE-05 | — | Pending |
| PRUNE-06 | — | Pending |
| PRUNE-07 | — | Pending |
| PRUNE-08 | — | Pending |
| PRUNE-09 | — | Pending |
| STRUCT-01 | — | Pending |
| STRUCT-02 | — | Pending |
| STRUCT-03 | — | Pending |
| STRUCT-04 | — | Pending |
| TYPE-01 | — | Pending |
| TYPE-02 | — | Pending |
| TYPE-03 | — | Pending |
| TYPE-04 | — | Pending |
| TYPE-05 | — | Pending |
| TYPE-06 | — | Pending |
| API-01 | — | Pending |
| API-02 | — | Pending |
| API-03 | — | Pending |
| API-04 | — | Pending |
| API-05 | — | Pending |
| REDUN-01 | — | Pending |
| REDUN-02 | — | Pending |
| REDUN-03 | — | Pending |
| ALLOC-01 | — | Pending |
| ALLOC-02 | — | Pending |
| ALLOC-03 | — | Pending |
| ALLOC-04 | — | Pending |

**Coverage:**
- v1.1 requirements: 31 total
- Mapped to phases: 0
- Unmapped: 31

---
*Requirements defined: 2026-02-14*
*Last updated: 2026-02-14 after initial definition*
