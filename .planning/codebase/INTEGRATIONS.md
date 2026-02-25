# External Integrations

**Analysis Date:** 2026-02-25

## APIs & External Services

None. QuantumFurnace.jl is a self-contained scientific computing library with no external API calls at runtime. All computation is local.

## Data Storage

**Databases:**
- None. No database connections.

**File Storage (Local Filesystem):**
- Experiment results saved as BSON + companion `.txt` files
- Default output path: `results/{subdir}/` relative to package root (resolved via `Pkg.project().path` in `src/results.jl`)
- Subdirectory convention: `results/approx_gns/` for GNS configs, `results/kms/` for KMS configs
- Filename convention: `{db}_{n}_{beta}_{domain}_{date}.bson` (e.g., `gns_n4_beta10_energy_20260225.bson`)
- Companion `.txt` human-readable summaries written alongside each `.bson` file
- Archive of prior experiment data present in `experiments/` directory (`.bson` + `.txt` pairs)

**Serialization Format:**
- BSON via `BSON` package (`BSON.bson(path, dict)` / `BSON.load(path)`)
- All structs converted to plain `Dict{Symbol, Any}` before serialization (see `src/results.jl`)
- No direct struct serialization — avoids BSON type instability

**Caching:**
- None.

## Authentication & Identity

- None required. No network services accessed at runtime.

## Monitoring & Observability

**Error Tracking:**
- None. No external error tracking service.

**Logs:**
- Julia standard `@warn` macros for advisory warnings (e.g., KrylovKit partial convergence in `src/krylov_eigsolve.jl`, memory estimates)
- `ProgressMeter` progress bars for trajectory loops (`src/trajectories.jl`, `src/convergence.jl`)
- Companion `.txt` files written alongside BSON results for human-readable run summaries

**Experiment Provenance (captured automatically in `src/results.jl`):**
- Julia version (`string(VERSION)`)
- Timestamp (`Dates.now()`)
- Git commit hash (via `LibGit2.GitRepo` + `LibGit2.head_oid`; falls back to `"unknown"`)
- Thread count (`Threads.nthreads()`)
- Wall time in seconds

## CI/CD & Deployment

**Documentation Hosting:**
- GitHub Actions workflow at `.github/workflows/Documentation.yml`
- Triggers: push to `main`, any tag, any pull request
- Runner: `ubuntu-latest`
- Deploys Documenter.jl-generated docs using `DOCUMENTER_KEY` secret
- Secrets required: `GITHUB_TOKEN` (automatic), `DOCUMENTER_KEY` (repository secret)

**No test CI/CD pipeline detected.** The GitHub Actions workflow only builds documentation. Tests are run manually with `julia --project=@. -e 'using Pkg; Pkg.test()'`.

## HPC / Cluster Deployment

**Distributed Execution:**
- `ClusterManagers` package in `[extras]` enables Slurm cluster integration
- Simulation scripts (`simulations/main_liouv.jl`, `simulations/main_thermalize.jl`) use Julia `Distributed` with `@everywhere` and optional `addprocs(N, exeflags="--project=@.")`
- Pattern: `addprocs` commented out by default; uncommented for cluster runs
- Slurm batch submission: `simulations/run_julia.sbatch` (referenced path; no direct `.sbatch` found in repo root)

**BLAS Thread Control:**
- Benchmark scripts explicitly set `BLAS.set_num_threads(4)` before heavy computation
- No global BLAS thread configuration in library code itself

## Webhooks & Callbacks

**Incoming:** None.

**Outgoing:** None. All CI/CD callbacks are handled by GitHub Actions infrastructure automatically via `GITHUB_TOKEN`.

## Environment Configuration

**Required environment variables at runtime:** None.

**CI/CD secrets (GitHub Actions only):**
- `GITHUB_TOKEN` - Auto-provided by GitHub for docs deployment authentication
- `DOCUMENTER_KEY` - Repository secret for Documenter.jl SSH deployment key

**No `.env` files present.** No secrets in codebase. No credentials stored.

---

*Integration audit: 2026-02-25*
