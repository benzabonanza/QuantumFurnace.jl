---
paths:
  - "src/**/*.jl"
  - "test/**/*.jl"
  - "experiments/**/*.jl"
  - "simulations/**/*.jl"
  - "hamiltonians/**/*.jl"
---

# Julia Code Rules

## Coding Conventions

- Parametric types: `T` consistently, no hardcoded `Float64` in generic code
- Immutable structs for result types
- Numerical tests use `atol`/`rtol`, never exact equality for floats
- Library code lives in `src/`; read nearby code before writing — match existing patterns
- Use existing deps: `LinearAlgebra`, `SparseArrays`, `LsqFit`, `Roots` — don't add new deps without justification
- Performance: preallocate, `@views` for slices, avoid allocations in hot loops
- No `Any` in struct fields or hot-path return types

## Reuse Over New Surface Area (codebase only)

**Before writing any new function or struct in `src/` or `test/`, search the codebase for existing functionality that already does the job or could be extended.**

- Grep for the operation, the relevant types, and the surrounding domain (e.g. `construct_lindbladian`, `apply_lindbladian!`, `materialize_discriminant`, `kms_*`, `oft!`, `lindblad_action_integrate`, `predict_lindbladian_trajectory`) before introducing a parallel path.
- If existing code is 80% of what you need, **extend it** (extra kwarg, generalize a type parameter, factor out a helper) rather than writing a near-duplicate.
- No "hacky" parallel functions — quick standalone helpers that re-implement existing behavior just to avoid touching the original are the wrong move; they add layers, drift from the canonical path, and rot tests.
- Prefer extending an existing struct field / Config kwarg over inventing a new struct that wraps the old one.
- This applies to **codebase work only**. Standalone scripts in `scripts/` (e.g. `scratch_*.jl`) can do whatever is fastest to understand a problem before promoting the result into `src/`.

## Threading Defaults for Lindbladian Work

For any Lindbladian simulation / matvec / Krylov work:

- `JULIA_NUM_THREADS = max` (all available cores)
- `BLAS.set_num_threads(1)`

This is the empirically optimal regime for our workloads (per qf-e60 / Krylov-trajectory benchmarks). Julia-thread parallelism over ω-grids / per-jump assembly outweighs nested BLAS threading on our matrix sizes. Don't override this in code without explicit measurement.

## Krylov First

When solving for fixed points, spectral gaps, or trajectory simulations, prefer Krylov methods over dense eigendecompositions whenever possible:

- Fixed point + gap: `krylov_spectral_gap(config, ham, jumps; ...)` from `src/krylov_eigsolve.jl`
- Trajectory simulation: `predict_lindbladian_trajectory` (Lindbladian) / `predict_channel_trajectory` (Channel) from `src/lindblad_action.jl`
- Lindbladian action: `apply_lindbladian!` (matrix-free) and `lindblad_action_integrate`
- Only build dense `d² × d²` superoperators when genuinely needed (norm comparisons, KMS geometry, dense cross-checks at small n)
- `eigen()` on a dense superoperator for fixed-point / gap analysis is forbidden in production code paths

## Register Sizing from the Sweep Table

For QPE / quadrature register sizes (`num_energy_bits_X`, `t0_X`, `w0_X` for X ∈ {D, b_minus, b_plus}), use the values from the canonical sweep table — do not invent or guess.

- Reference: [Quadrature register recipe qf-7xt](`.claude-memory/quadrature_register_recipe_qf_7xt.md`) and the companion summary at `drafts/error-analysis/quadrature-convergence-summary.md`.
- Use Config's per-term register triples (`register_t0_X / register_w0_X / register_r_X`) rather than the legacy global `t0/w0/num_energy_bits` kwargs in new code.
- If a fixture lands outside the measured `(n, β)` range and you need to extend it, run a sweep script in `scripts/` first; don't hand-pick numbers in `src/`.

## Where Plans Live

- **`drafts/` is for thesis drafts only.** Do not write coding plans, refactor notes, or implementation TODOs there unless explicitly told to.
- **Use `.planning/` for implementation plans**, but skip the old "Phases" plugin convention (`phases/NN-name/...`) — that pattern was retired for being too verbose. Write a single concise plan file or short topic directory; whatever fits the task.
- Per-task tracking belongs in `bd` (beads), not in markdown files under either `drafts/` or `.planning/`.

## Architecture

- **Domain hierarchy**: `BohrDomain` -> `EnergyDomain` -> `TimeDomain` -> `TrotterDomain`
- **Core structs**: `HamHam` (Hamiltonian + disorder), `TrottTrott` (Trotterized), `Workspace`
- **Fitting**: `FitResult`, `BiexpFitResult`, `MixingTimeEstimate` — all immutable, via LsqFit.jl
- **Simulation**: Krylov (sparse, large n) vs dense (exact, small n)
- **Disordering**: `find_ideal_heisenberg` with `disordering_terms` kwarg; Z-only or Z+ZZ for symmetry breaking

## Test Suite: SANDBOX vs NO_SANDBOX (qf-5nz)

The test suite is split into two tiers in `test/runtests.jl`:

- **`SANDBOX_FILES`** (default): every test file that runs inside the 3.5 GB / few-minute sandbox container. This is the entire suite the assistant runs by default with `Pkg.test()`.
- **`NO_SANDBOX_FILES`** + `test/trajectory_validation/`: heavier tests gated behind `QUANTUMFURNACE_FULL_TESTS=true`. Some individual subtests inside SANDBOX files are also gated inline by the same env var (e.g. `test_dll_kms_db.jl::(j)`, `test_lindblad_action.jl::(i)`); they are tagged `[NO_SANDBOX]` in the testset name and explained in a comment above.

Default test command (sandbox-safe — this is what the assistant must use):

```
JULIA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 \
  julia --project -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=1500M"])'
```

The `--heap-size-hint=1500M` flag is **required** inside the sandbox: it caps
Julia's heap target so cumulative buffers from the ~50 test files do not push
RSS past the ~1.8 GB available to the test process. Empirically the suite
peaks at ~1.3 GB RSS at this setting, leaving comfortable headroom.
`runtests.jl` prints a per-file RSS trace on stderr so OOM regressions are
immediately visible.

Full suite — only safe **outside the sandbox** (developer machine with ≥ 16 GB RAM):

```
QUANTUMFURNACE_FULL_TESTS=true julia --project -e 'using Pkg; Pkg.test()'
```

**Inside the sandbox, only run the SANDBOX tier.** Do not set `QUANTUMFURNACE_FULL_TESTS=true` from a sandbox session — it will OOM or time out, exactly as before this split was introduced.

When adding a new test:
- Default it to `SANDBOX_FILES`.
- If a subtest cannot fit the sandbox envelope without losing physics validity, gate just that subtest inline with `if get(ENV, "QUANTUMFURNACE_FULL_TESTS", "false") == "true"`. Tag the testset name with `[NO_SANDBOX]` and document the reason in a comment above the gate.
- Only add a whole file to `NO_SANDBOX_FILES` when most of it is unavoidably heavy. Prefer inline gating.

When tuning a sandbox-bound test, prefer reducing register sizes (`num_energy_bits`) over loosening tolerances. Empirically (qf-5nz), DLL TimeDomain Bohr↔Time agreement on n=3 is FINUFFT-floor-limited (~3e-9) by `Nt = 256`; `Nt = 4096` (N=12) buys nothing but uses 16× more memory. The DLL-test fixtures now run at `N = 10`.

## Atomic Commits

Every code change MUST be an atomic git commit. One logical change per commit.

- Each commit: single, self-contained change (one new function, one bug fix, one struct, one test file)
- Tests must pass after each commit
- Commit message format: `<type>(<scope>): <description>`
- Stage only the specific files — never `git add .` or `git add -A`
- Include `Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>`
- `src/` + `test/` in one commit is fine if the test is specifically for that change

## Physics Sanity Rule

Whenever code makes a decision based on physical or mathematical reasoning — parameter choices, expected scalings, algorithm selection, convergence criteria — spawn the `physics-checker` agent to verify. Flag any `# PHYSICS CHECK:` comments for review.
