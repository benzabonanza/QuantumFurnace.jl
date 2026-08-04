# Legacy stochastic trajectories

Archived on 2026-08-04 after the Task 2 default suite passed 6,867/6,867 tests.
This subsystem is frozen and unmaintained. It is not loaded by the active
`QuantumFurnace` module and its tests are not part of `Pkg.test()`.

The archival boundary is the stochastic pure-state/observable simulator. The
following former public entry points and types were deliberately removed:

- `Trajectory`, `TrajectoryScratch`, `TrajectoryResult`,
  `ObservableTrajectoryResult`, `TrajectoryResults`, and `ConvergenceData`
- `run_trajectory`, `run_trajectories`, `run_observable_trajectories`, and
  `step_along_trajectory!`
- `run_trajectories_convergence`, `run_trajectories_adaptive`, and
  `build_preset_trajectory_observables`

The matrix-free `predict_lindbladian_trajectory` and
`predict_channel_trajectory` functions are unrelated Krylov spectral
predictors. They remain active, exported, and tested.

Layout:

- `src/` contains the retired simulator and convergence implementation.
- `src/support/` contains exact pre-archival snapshots of the shared struct and
  result-serialization files so removed type definitions and branches remain
  recoverable. These snapshots also contain retained declarations and must not
  be included wholesale in a current checkout.
- `test/` contains whole-file stochastic tests and the former full-tier
  validation runners. Stochastic sections removed from mixed active tests can
  be recovered from repository history or the support snapshots where
  applicable.
- Private experiment and simulation callers live in the companion archive at
  `private/archive/legacy_stochastic_trajectories/`.

These files depend on the package layout and internal helpers as they existed
at archival time. They are historical source material, not a standalone or
supported compatibility package.
