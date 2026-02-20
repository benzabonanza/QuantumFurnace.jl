# Roadmap: QuantumFurnace.jl

## Milestones

- ✅ **v1.0 Trajectories** -- Phases 1-5 (shipped 2026-02-14)
- ✅ **v1.1 Reduce** -- Phases 6-11 (shipped 2026-02-15)
- ✅ **v1.2 Multi-threading** -- Phases 12-19 (shipped 2026-02-16)
- ✅ **v1.3 Mixing Time Estimation** -- Phases 20-25 (shipped 2026-02-19)
- ✅ **v1.4 Spectral Gap Refinement** -- Phase 26 (partial, shipped 2026-02-20)

## Phases

<details>
<summary>✅ v1.0 Trajectories (Phases 1-5) -- SHIPPED 2026-02-14</summary>

- [x] Phase 1: Foundation and Compilation (1/1 plans) -- completed 2026-02-13
- [x] Phase 2: Trajectory Bug Fixes (2/2 plans) -- completed 2026-02-14
- [x] Phase 3: DM Reference Test Suite (3/3 plans) -- completed 2026-02-14
- [x] Phase 4: Trajectory Cross-Validation (2/2 plans) -- completed 2026-02-14
- [x] Phase 5: Statistical Validation and Regression (2/2 plans) -- completed 2026-02-14

Full details: [milestones/v1.0-ROADMAP.md](milestones/v1.0-ROADMAP.md)

</details>

<details>
<summary>✅ v1.1 Reduce (Phases 6-11) -- SHIPPED 2026-02-15</summary>

- [x] Phase 6: Dead Code Pruning (2/2 plans) -- completed 2026-02-15
- [x] Phase 7: DRY Refactoring (2/2 plans) -- completed 2026-02-15
- [x] Phase 8: Struct Simplification (3/3 plans) -- completed 2026-02-15
- [x] Phase 9: Type Parameterization (3/3 plans) -- completed 2026-02-15
- [x] Phase 10: API Surface Cleanup (3/3 plans) -- completed 2026-02-15
- [x] Phase 11: Allocation Optimization (3/3 plans) -- completed 2026-02-15

Full details: [milestones/v1.1-ROADMAP.md](milestones/v1.1-ROADMAP.md)

</details>

<details>
<summary>✅ v1.2 Multi-threading (Phases 12-19) -- SHIPPED 2026-02-16</summary>

- [x] Phase 12: Workspace Refactor (2/2 plans) -- completed 2026-02-15
- [x] Phase 13: Multi-Threaded Trajectory Engine (2/2 plans) -- completed 2026-02-16
- [x] Phase 14: GNS Trajectory Path (1/1 plans) -- completed 2026-02-16
- [x] Phase 15: Data Architecture (2/2 plans) -- completed 2026-02-16
- [x] Phase 16: Convergence Tracking (2/2 plans) -- completed 2026-02-16
- [x] Phase 17: Adaptive Sampling (2/2 plans) -- completed 2026-02-16
- [x] Phase 18: KMS-vs-GNS Experiments (1/1 plans) -- completed 2026-02-16
- [x] Phase 19: Logic Simplification (3/3 plans) -- completed 2026-02-16

Full details: [milestones/v1.2-ROADMAP.md](milestones/v1.2-ROADMAP.md)

</details>

<details>
<summary>✅ v1.3 Mixing Time Estimation (Phases 20-25) -- SHIPPED 2026-02-19</summary>

- [x] Phase 20: Observable Infrastructure (1/1 plans) -- completed 2026-02-17
- [x] Phase 21: Exponential Fitting (1/1 plans) -- completed 2026-02-17
- [x] Phase 22: Observable-Only Trajectory Runner (1/1 plans) -- completed 2026-02-17
- [x] Phase 23: Gap Estimation API (1/1 plans) -- completed 2026-02-17
- [x] Phase 24: Cross-Validation (3/3 plans) -- completed 2026-02-17
- [x] Phase 25: Spectral Gap Validation Overhaul (3/3 plans) -- completed 2026-02-18

Full details: [milestones/v1.3-ROADMAP.md](milestones/v1.3-ROADMAP.md)

</details>

<details>
<summary>✅ v1.4 Spectral Gap Refinement (Phase 26) -- PARTIAL, SHIPPED 2026-02-20</summary>

- [x] Phase 26: Exact Reference and Structural Diagnostics (2/2 plans) -- completed 2026-02-19

Phases 27-30 deferred (two-exponential fitting, effective rates, bootstrap, dashboard).

Full details: [milestones/v1.4-ROADMAP.md](milestones/v1.4-ROADMAP.md)

</details>

## Progress

| Phase | Milestone | Plans Complete | Status | Completed |
|-------|-----------|----------------|--------|-----------|
| 1. Foundation and Compilation | v1.0 | 1/1 | Complete | 2026-02-13 |
| 2. Trajectory Bug Fixes | v1.0 | 2/2 | Complete | 2026-02-14 |
| 3. DM Reference Test Suite | v1.0 | 3/3 | Complete | 2026-02-14 |
| 4. Trajectory Cross-Validation | v1.0 | 2/2 | Complete | 2026-02-14 |
| 5. Statistical Validation and Regression | v1.0 | 2/2 | Complete | 2026-02-14 |
| 6. Dead Code Pruning | v1.1 | 2/2 | Complete | 2026-02-15 |
| 7. DRY Refactoring | v1.1 | 2/2 | Complete | 2026-02-15 |
| 8. Struct Simplification | v1.1 | 3/3 | Complete | 2026-02-15 |
| 9. Type Parameterization | v1.1 | 3/3 | Complete | 2026-02-15 |
| 10. API Surface Cleanup | v1.1 | 3/3 | Complete | 2026-02-15 |
| 11. Allocation Optimization | v1.1 | 3/3 | Complete | 2026-02-15 |
| 12. Workspace Refactor | v1.2 | 2/2 | Complete | 2026-02-15 |
| 13. Multi-Threaded Trajectory Engine | v1.2 | 2/2 | Complete | 2026-02-16 |
| 14. GNS Trajectory Path | v1.2 | 1/1 | Complete | 2026-02-16 |
| 15. Data Architecture | v1.2 | 2/2 | Complete | 2026-02-16 |
| 16. Convergence Tracking | v1.2 | 2/2 | Complete | 2026-02-16 |
| 17. Adaptive Sampling | v1.2 | 2/2 | Complete | 2026-02-16 |
| 18. KMS-vs-GNS Experiments | v1.2 | 1/1 | Complete | 2026-02-16 |
| 19. Logic Simplification | v1.2 | 3/3 | Complete | 2026-02-16 |
| 20. Observable Infrastructure | v1.3 | 1/1 | Complete | 2026-02-17 |
| 21. Exponential Fitting | v1.3 | 1/1 | Complete | 2026-02-17 |
| 22. Observable-Only Trajectory Runner | v1.3 | 1/1 | Complete | 2026-02-17 |
| 23. Gap Estimation API | v1.3 | 1/1 | Complete | 2026-02-17 |
| 24. Cross-Validation | v1.3 | 3/3 | Complete | 2026-02-17 |
| 25. Spectral Gap Validation Overhaul | v1.3 | 3/3 | Complete | 2026-02-18 |
| 26. Exact Reference and Structural Diagnostics | v1.4 | 2/2 | Complete | 2026-02-19 |
