# Milestones: QuantumFurnace.jl

## v1.0 Trajectories (Shipped: 2026-02-14)

**Started:** 2026-02-13 | **Shipped:** 2026-02-14
**Phases:** 1-5 (10 plans + 12 quick tasks) | **Tests:** 0 → 224 | **Commits:** 85
**Julia LOC:** 13,082 total (2,034 test) | **Files changed:** 18 (+1,793 / -243)
**Git range:** ecd17ae..5dc2f45

**Delivered:** Trajectory simulation fixed, validated against density matrix evolution, and locked down with a comprehensive 224-test correctness suite covering CPTP channels, detailed balance, error scaling, cross-validation, and regression.

**Key accomplishments:**
1. Fixed trajectory compilation (5 bugs) and the critical U_B ordering bug in EnergyDomain
2. CPTP channel verification: K0†K0 + δR + U_res†U_res = I to 1e-10 for all domains
3. DM ground truth: BohrDomain detailed balance at 1.6e-15, domain error hierarchy verified, O(δ²)/O(δ) scaling confirmed
4. Trajectory-vs-DM cross-validation: single-step O(δ²) scaling and multi-step Gibbs convergence for Energy, Time, Trotter
5. Per-operator Lie-Trotter splitting refactor: TrotterDomain Gibbs fixed point distance 0.004 → 9e-9
6. Portable regression framework: DM-based comparison replacing frozen BSON for cross-platform/version stability

**Archives:** [v1.0-ROADMAP.md](milestones/v1.0-ROADMAP.md) | [v1.0-REQUIREMENTS.md](milestones/v1.0-REQUIREMENTS.md)

---

