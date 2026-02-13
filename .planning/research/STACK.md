# Stack Research

**Domain:** Julia quantum Lindbladian simulation package (open quantum systems, Gibbs sampling)
**Researched:** 2026-02-13
**Confidence:** MEDIUM-HIGH (most recommendations verified via official docs/repos; some version numbers from training data marked accordingly)

## Current State of QuantumFurnace.jl

Before recommending the stack, here is what the project already uses (from `Project.toml`):

| Dependency | Role | Status |
|------------|------|--------|
| LinearAlgebra (stdlib) | Dense matrix ops, eigen decomposition | Keep |
| SparseArrays (stdlib) | Sparse Hamiltonians, Liouvillians | Keep |
| Arpack | Sparse eigensolve for Liouvillian steady state | **Replace** |
| FINUFFT | Non-uniform FFT for time/energy domain | Keep |
| Distributed / SharedArrays / ClusterManagers | Parallel jump operator construction | Keep (with caveats) |
| LinearMaps | Matrix-free Liouvillian representation | Keep |
| QuadGK | Numerical integration | Keep |
| Roots | Root finding | Keep |
| Optim | Optimization (log-Sobolev bounds) | Keep |
| SpecialFunctions | erfc for filter functions | Keep |
| DataStructures | OrderedDict for Bohr frequencies | Keep |
| BSON | Serialization | Keep |
| ProgressMeter | Progress bars | Keep |
| Plots | Visualization | **Reconsider** |
| Revise / Debugger / BenchmarkTools | Dev tooling | Keep (dev-only) |
| Documenter / Literate / DocumenterTools | Documentation | Keep, upgrade |

Julia version: **1.12.4** (Manifest confirms; 1.12.5 available as of Feb 2026).

## Recommended Stack

### Core Technologies (Keep / Already Present)

| Technology | Version | Purpose | Why Recommended |
|------------|---------|---------|-----------------|
| Julia | >= 1.12 | Language runtime | Already in use. Julia 1.12 brings redefinable types (better dev workflow) and code trimming. LTS is 1.10 but 1.12 is fine for a research package. |
| LinearAlgebra (stdlib) | 1.12.x | Dense linear algebra, eigen, matrix exp | Standard, zero-dep. Handles all dense ops for <= 12 qubits (4096x4096 matrices). No reason to change. |
| SparseArrays (stdlib) | 1.12.x | Sparse matrix storage and ops | CSC format. Standard for quantum Hamiltonians. Keep. |
| FINUFFT.jl | >= 3.4.2 | Non-uniform FFT | Core to the time-domain and energy-domain Lindbladian construction. No Julia-native competitor of comparable quality. Keep. |
| LinearMaps.jl | latest | Matrix-free linear operators | Used for Liouvillian-as-linear-map (avoids materializing d^4 superoperator). Essential for scaling. Keep. |
| QuadGK.jl | latest | Adaptive Gauss-Kronrod quadrature | Used for numerical integrals in filter functions. Lightweight, standard. Keep. |
| Roots.jl | latest | Root finding | Lightweight. Keep. |
| Optim.jl | >= 1.13 | Optimization | Used for log-Sobolev bound computation. Mature, well-maintained. Keep. |
| SpecialFunctions.jl | latest | erfc and related | Only importing erfc. Lightweight dep. Keep. |
| DataStructures.jl | latest | OrderedDict for Bohr frequency bookkeeping | Lightweight. Keep. |
| BSON.jl | latest | Data serialization | Used for saving/loading results. Keep. |

**Confidence: HIGH** -- These are all verified in the existing Project.toml and are standard Julia scientific computing choices.

### New Dependencies to Add

| Technology | Version | Purpose | Why Recommended |
|------------|---------|---------|-----------------|
| KrylovKit.jl | >= 0.8 | Sparse eigensolve (replace Arpack) | **Replace Arpack.jl.** Arpack wraps Fortran ARPACK which is (a) not thread-safe (single-threaded only), (b) lacks specialized Hermitian routines, (c) less stable than KrylovKit for complex Hermitian problems. KrylovKit provides native Julia Krylov-Schur/Lanczos with Hermitian specialization, thread-safety, and works with any `AbstractMatrix` or callable. Used by ITensors, QuantumOptics ecosystem. |
| ExponentialUtilities.jl | >= 1.27 | Matrix exponential, Krylov expmv | For trajectory simulation: `exp(t*A)*v` without forming full matrix exponential. Krylov-based `expv!` is O(m*n) vs O(n^3) for dense expm. At 12 qubits the Liouvillian is 16M entries dense -- Krylov expmv avoids this. Also useful for Trotterized time evolution validation. Part of SciML ecosystem. |
| PythonCall.jl | >= 0.9.23 | Julia-Python interop for Qiskit | **Use PythonCall, not PyCall.** PythonCall is the modern replacement: type-stable returns, no lossy conversions, no numpy dependency, uses PATH Python by default (or CondaPkg). Needed for Qiskit circuit generation/resource estimation. PyCall is legacy. |
| CondaPkg.jl | latest | Manage Python deps for PythonCall | Companion to PythonCall. Declares Python dependencies (qiskit, qiskit-aer) in a CondaPkg.toml so they install automatically. Reproducible Python env. |
| OrdinaryDiffEqTsit5.jl | >= 1.x | Lightweight ODE solver | For trajectory validation: compare trajectory-averaged density matrix evolution against direct Lindblad master equation integration. Use the split sub-package (not full DifferentialEquations.jl which pulls 100+ deps). Tsit5 is the recommended default non-stiff solver. |
| Aqua.jl | >= 0.8 | Automated package quality assurance | Test for method ambiguities, stale deps, missing compat entries, undefined exports. Standard for Julia packages aiming for registry. Add to test deps. |
| JuliaFormatter.jl | >= 1.0 | Code formatting | Enforces consistent style. Use `.JuliaFormatter.toml` with SciML style (4-space indent, consistent with scientific Julia convention). Add as dev tool, optionally enforce in CI. |

**Confidence: HIGH** for KrylovKit (verified: widely used, official docs confirm Hermitian specialization), PythonCall (verified: official comparison page confirms advantages over PyCall), Aqua (verified: standard practice).
**Confidence: MEDIUM** for ExponentialUtilities (version number from training data; functionality verified via SciML docs), OrdinaryDiffEqTsit5 (split package approach verified; exact version from training data).

### Testing Stack

| Tool | Purpose | When to Use |
|------|---------|-------------|
| Test (stdlib) | Core test macros (@test, @testset, @test_throws) | Always. Foundation of all testing. |
| Aqua.jl | Package quality checks | Run in CI. Catches ambiguities, stale deps, missing compat before they become bugs. |
| SafeTestsets.jl | Isolated test sets | Use for test files that might leak state. Each @safetestset runs in its own module. Prevents test interference. |
| BenchmarkTools.jl | Performance regression testing | Already a dep. Use @btime/@benchmark in dedicated performance tests, not in CI (too noisy). |

**Confidence: HIGH** -- Test stdlib is standard. SafeTestsets verified via JuliaTesting org. Aqua verified.

### Documentation Stack

| Tool | Version | Purpose | Notes |
|------|---------|---------|-------|
| Documenter.jl | >= 1.16 | Generate HTML docs from docstrings + markdown | Already in extras. v1.16.1 is latest (Nov 2025). Julia 1.12 workspace support simplifies docs/Project.toml setup. |
| Literate.jl | >= 2.20 | Literate programming for tutorials | Already in extras. Generates Documenter pages AND Jupyter notebooks from same .jl source. Essential for physics packages -- users expect runnable examples. |
| DocumenterTools.jl | >= 0.1.20 | Doc generation helpers | Already present. Keep. |

**Confidence: HIGH** -- Already in use, versions verified from Project.toml compat entries and Documenter.jl releases page.

### CI / Infrastructure

| Tool | Purpose | Notes |
|------|---------|-------|
| GitHub Actions | CI/CD | Standard for Julia packages. Use julia-actions/setup-julia, julia-actions/julia-runtest, julia-actions/julia-processcoverage. |
| julia-actions/setup-julia | Install Julia in CI | Use version matrix: ['1.12', '1'] to test current and latest. |
| julia-actions/julia-runtest | Run Pkg.test() in CI | Set annotate: true for PR annotations on test failures. |
| julia-actions/julia-processcoverage | Generate coverage reports | Outputs lcov.info for Codecov. |
| codecov/codecov-action | Upload coverage | Works without token for public repos. |
| JuliaFormatter check | Formatting CI | Run `using JuliaFormatter; format("src", verbose=true, overwrite=false)` and fail if changes detected. |

**Confidence: HIGH** -- julia-actions verified via GitHub repos. Standard Julia package CI pattern.

### Development Tools

| Tool | Purpose | Notes |
|------|---------|-------|
| Revise.jl | Hot-reload during development | Already a dep. Essential for Julia dev workflow. Keep as dev-only. |
| Debugger.jl | Step debugging | Already a dep. Keep as dev-only. Not needed at runtime. |
| BenchmarkTools.jl | Microbenchmarking | Already a dep. Keep. |
| PrecompileTools.jl | Reduce TTFX | Add @compile_workload blocks for common entry points (run_lindbladian, run_thermalization). Julia 1.12 benefits from native code caching. Worth adding once API stabilizes. |

**Confidence: HIGH** for existing tools. **MEDIUM** for PrecompileTools (defer until API is stable; premature optimization otherwise).

## Alternatives Considered

| Recommended | Alternative | When to Use Alternative |
|-------------|-------------|-------------------------|
| KrylovKit.jl | Arpack.jl | Never for new code. Arpack is not thread-safe and lacks Hermitian specialization. Only keep if you need exact backward compatibility with existing saved eigendecompositions. |
| PythonCall.jl | PyCall.jl | Only if you must interop with an existing Julia package that hard-depends on PyCall (e.g., older qiskit-alt). For new Qiskit interop, PythonCall is strictly better. |
| OrdinaryDiffEqTsit5 | DifferentialEquations.jl | Use full DifferentialEquations.jl only if you need stiff solvers (Rosenbrock, BDF) or SDE solvers. For non-stiff trajectory validation ODEs, the split sub-package keeps deps minimal. |
| OrdinaryDiffEqTsit5 | Hand-rolled RK4 | Only for pedagogical purposes or if the ODE is trivial. SciML solvers handle adaptive stepping, error control, and dense output -- don't reinvent this. |
| ExponentialUtilities.jl | Manual dense expm (LinearAlgebra.exp) | At 12 qubits the Hilbert space is 4096, Liouvillian is 16M. Dense expm is fine for Hilbert-space ops but not for superoperators. For trajectory simulation (Hilbert-space-sized), dense expm from LinearAlgebra may suffice. ExponentialUtilities matters for larger problems or superoperator expm. |
| Test + SafeTestsets | TestItemRunner.jl | TestItemRunner is VS Code-centric (run individual test items from editor). Good for DX but adds coupling to IDE. SafeTestsets is lighter and CI-friendly. Consider TestItemRunner only if heavy VS Code usage. |
| Plots.jl | Makie.jl (CairoMakie) | Makie is faster, more modern, GPU-accelerated. But Plots.jl is already a dep and has simpler API for quick research plots. Switch to CairoMakie only if you need publication-quality figures or interactive 3D visualization. |
| BSON.jl | JLD2.jl | JLD2 is HDF5-based, more robust for large arrays and Julia types. BSON can silently lose type information. Consider migrating to JLD2 for result serialization -- but not urgent. |

## What NOT to Use

| Avoid | Why | Use Instead |
|-------|-----|-------------|
| PyCall.jl (for new code) | Legacy package. Lossy type conversions, requires numpy, not type-stable. Single-threaded Python GIL interaction is poorly handled. | PythonCall.jl |
| Arpack.jl (for new eigensolve code) | Fortran ARPACK is not re-entrant (single-thread only). No Hermitian specialization. Less stable for complex eigenvalue problems. | KrylovKit.jl |
| Full DifferentialEquations.jl | Pulls ~100 transitive dependencies. Massive compilation overhead. You only need one ODE solver for validation. | OrdinaryDiffEqTsit5 (split sub-package) |
| QuantumOptics.jl / QuantumToolbox.jl (as dependency) | These are full frameworks with their own type systems (Ket, Bra, Operator with basis tracking). QuantumFurnace has its own types (HamHam, TrottTrott, JumpOp, LiouvConfig). Importing a framework would force type conversion layers everywhere and add massive dep trees. | Keep custom types. Optionally implement conversion functions for interop. |
| Yao.jl (as dependency) | Circuit-centric framework. QuantumFurnace does continuous-time Lindbladian evolution, not gate-based. Yao adds no value here. For circuit resource estimation, talk to Qiskit via PythonCall. | PythonCall + Qiskit |
| CUDA.jl / GPU computing | At <= 12 qubits, Hilbert space dim is <= 4096. Dense matrices are 4096x4096 = 128 MB. GPU overhead (data transfer, kernel launch) dominates at this scale. GPU only helps at >= 16 qubits. | CPU with BLAS threading |
| Distributed.jl for < 8 qubits | Interprocess communication overhead exceeds computation time for small systems. | Base.Threads (already using via Base.Threads) for < 8 qubits. Distributed for >= 8 qubits where jump operator construction is expensive. |

## Stack Patterns by Variant

**If adding Qiskit circuit generation:**
- Add PythonCall.jl + CondaPkg.jl
- Create `CondaPkg.toml` with `qiskit >= 1.0` and `qiskit-aer`
- Wrap Qiskit calls in a dedicated `src/qiskit_interop.jl` behind an extension or optional import
- Because: Not all users need Python. Use Julia package extensions (Julia 1.9+) to make PythonCall a weak dependency

**If validating trajectories against master equation:**
- Add OrdinaryDiffEqTsit5 for reference Lindblad integration
- Use ExponentialUtilities.jl for matrix-exponential-based time stepping
- Because: Trajectory simulation correctness needs an independent reference. ODE integration of the master equation provides this.

**If targeting Julia General Registry:**
- Add Aqua.jl tests, compat bounds for all deps, CI with julia-runtest
- Remove dev-only deps from [deps] (Revise, Debugger, Plots, BenchmarkTools should be in [extras] only)
- Because: Registry requires compat entries. Aqua catches missing ones. Dev tools in [deps] force all users to install them.

**If supporting > 12 qubits in the future:**
- Replace dense matrix operations with sparse-only paths
- Consider KrylovKit + ExponentialUtilities for all eigensolve/expmv
- Because: At 14 qubits, Hilbert dim = 16384, Liouvillian dim = 268M entries. Dense is infeasible.

## Version Compatibility

| Package A | Compatible With | Notes |
|-----------|-----------------|-------|
| Julia 1.12 | All recommended packages | 1.12 is recent; all active packages support it. FINUFFT.jl is tested up to 1.11 per docs -- verify 1.12 compat (likely works, flag for testing). |
| PythonCall >= 0.9 | Julia >= 1.10, Python >= 3.10 | Requires Python in PATH or CondaPkg. Does NOT work with PyCall simultaneously -- they conflict on Python process ownership. |
| OrdinaryDiffEqTsit5 | Julia >= 1.10 | Part of SciML ecosystem. Depends on SciMLBase, which may pull more deps than expected. Test actual dep tree. |
| KrylovKit >= 0.8 | Julia >= 1.6 | Pure Julia, minimal deps. Drop-in replacement for Arpack eigs() calls. |
| ExponentialUtilities >= 1.27 | Julia >= 1.10 | Part of SciML. Check transitive deps. |
| Documenter >= 1.16 | Julia >= 1.12 workspace feature | Use docs/ workspace in Project.toml for clean separation. |

## Cleanup: Dependencies to Move or Remove

The current `Project.toml` has dev-only tools in `[deps]` that should move:

| Package | Current Location | Should Be | Why |
|---------|-----------------|-----------|-----|
| Revise | [deps] | [extras] only | Dev tool. Users don't need it. |
| Debugger | [deps] | [extras] only | Dev tool. |
| BenchmarkTools | [deps] | [extras] only | Already in [extras], remove from [deps]. |
| Plots | [deps] | [extras] or separate | Visualization. Heavy dep. Not needed by core simulation. |
| Pkg | [deps] | Remove | stdlib, auto-available. No need to declare. |
| BSON | [deps] | Keep, but consider JLD2 long-term | BSON can lose Julia type info on roundtrip. |
| DocumenterTools | [deps] | [extras] with docs target | Only needed for doc generation. |

## Installation

```julia
# In Project.toml [deps], add:
# KrylovKit = "0f1e7b(...)  (use Pkg.add to get UUID)
# ExponentialUtilities = "d4d017(...)"

# From Julia REPL:
using Pkg
Pkg.add("KrylovKit")
Pkg.add("ExponentialUtilities")

# For Qiskit interop (when ready):
Pkg.add("PythonCall")
Pkg.add("CondaPkg")

# For trajectory validation:
Pkg.add("OrdinaryDiffEqTsit5")

# For test quality:
Pkg.add("Aqua")  # add to test deps

# For formatting:
Pkg.add("JuliaFormatter")  # dev tool only
```

## Sources

- [QuantumToolbox.jl paper (Quantum journal, Sept 2025)](https://quantum-journal.org/papers/q-2025-09-29-1866/) -- MEDIUM confidence, benchmarks and architecture patterns
- [KrylovKit.jl official docs (Oct 2025)](https://jutho.github.io/KrylovKit.jl/stable/) -- HIGH confidence, Hermitian Lanczos specialization verified
- [PythonCall.jl vs PyCall comparison (official docs)](https://juliapy.github.io/PythonCall.jl/v0.2/pycall/) -- HIGH confidence, feature comparison verified
- [ExponentialUtilities.jl (SciML docs)](https://docs.sciml.ai/ExponentialUtilities/stable/) -- HIGH confidence, Krylov expmv verified
- [OrdinaryDiffEq.jl split packages (SciML GitHub)](https://github.com/SciML/OrdinaryDiffEq.jl) -- MEDIUM confidence, split approach verified but exact version from training data
- [Aqua.jl (JuliaTesting, v0.8.14 Aug 2025)](https://github.com/JuliaTesting/Aqua.jl) -- HIGH confidence
- [Documenter.jl v1.16.1 (Nov 2025)](https://documenter.juliadocs.org/stable/) -- HIGH confidence
- [Julia 1.12 release blog (Oct 2025)](https://julialang.org/blog/2025/10/julia-1-12-highlights/) -- HIGH confidence
- [julia-actions GitHub org](https://github.com/julia-actions) -- HIGH confidence, CI actions verified
- [JuliaFormatter.jl](https://github.com/domluna/JuliaFormatter.jl) -- HIGH confidence
- [Arpack.jl threading limitation (Julia Discourse)](https://discourse.julialang.org/t/suggestions-needed-diagonalizing-large-hermitian-sparse-matrix/96580) -- MEDIUM confidence, community reports
- [FINUFFT.jl v3.x docs](https://ludvigak.github.io/FINUFFT.jl/latest/) -- HIGH confidence
- [QuantumOptics.jl architecture](https://docs.qojulia.org/) -- HIGH confidence, type system patterns
- [PrecompileTools.jl (JuliaLang)](https://github.com/JuliaLang/PrecompileTools.jl) -- HIGH confidence

---
*Stack research for: Julia quantum Lindbladian simulation (QuantumFurnace.jl)*
*Researched: 2026-02-13*
