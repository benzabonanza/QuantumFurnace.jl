# QuantumFurnace historical archive

This tree contains code retired from the active package on 2026-08-04. It is
kept as historical reference only: archived files are frozen, unmaintained,
not included by `QuantumFurnace`, and not run by `Pkg.test()`.

- `legacy_stochastic_trajectories/` preserves the obsolete stochastic
  pure-state/observable simulator and its tests.
- `inactive_staging/` preserves non-included staging implementations and
  staging tests that were removed from the active source and test trees.

The private repository has a complementary
`private/archive/legacy_stochastic_trajectories/` tree containing historical
experiment and simulation drivers. Keeping those drivers private avoids
publishing thesis-internal workflows in the package repository.

The last green baseline before archival was the Task 2 default suite:
6,867/6,867 tests passed. Archived material depends on the source layout and
internal APIs as they existed at that baseline and is not expected to run
against later package versions without reconstruction work.
