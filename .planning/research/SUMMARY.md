# Project Research Summary

**Project:** QuantumFurnace.jl — Krylov-Based Lindbladian Spectral Gap Estimation
**Domain:** Matrix-free Krylov eigensolving for open quantum system spectral analysis (n=8-12 qubits)
**Researched:** 2026-02-20
**Confidence:** HIGH

## Executive Summary

QuantumFurnace.jl is a mature Julia package (26+ phases) for Lindbladian simulation of open quantum systems. The existing codebase provides dense Liouvillian construction and trajectory-based gap estimation, both limited to n<=6 qubits. The next milestone extends spectral gap estimation to n=8-12 qubits using matrix-free Krylov methods: the Lindbladian superoperator is never formed as an explicit matrix, only applied as a callable function L(rho). At n=12, the full dim^2 x dim^2 superoperator would be 4 petabytes; the Krylov approach stores only ~30-50 vectors of size dim^2 (7-12 GB total) and finds the spectral gap eigenvalue by repeated application of the Lindbladian action.

The critical structural insight enabling this without shift-invert is that all Lindbladian eigenvalues satisfy Re(lambda) <= 0, with the steady state at lambda=0 as the extremal eigenvalue under largest-real-part sorting. The spectral gap eigenvalue — typically an interior eigenvalue in other contexts — is exactly the second-most-extremal eigenvalue under `:LR` sorting. KrylovKit.jl's `eigsolve(..., :LR)` finds it directly. The existing precomputed jump data (NUFFT prefactors, OFT results, Bohr frequency dictionaries) is reused verbatim inside the matvec closure, making each matvec an O(dim^3 * n_jumps * n_freqs) operation rather than the O(dim^4) cost of assembling the full superoperator. At n=8 (dim=256), a full Krylov solve takes an estimated ~10 minutes where dense diagonalization is already impossible (64 GB matrix).

The primary risks are correctness bugs that produce plausible-looking but wrong eigenvalues: vectorization convention mismatches (reshape/transpose errors in the matvec), basis selection errors (using computational-basis jump operators instead of eigenbasis ones), and confusing the CPTP channel E with the Lindbladian generator L (which introduces catastrophic cancellation when computing (E(rho)-rho)/delta). All three are prevented by a mandatory round-trip test at n=4 before any eigsolve is attempted. Secondary risks are memory blowup from misconfigured krylovdim at n=12, and non-convergence for systems with clustered eigenvalues, both mitigated by pre-flight memory estimation and requesting 6-10 eigenvalues rather than just 2.

## Key Findings

### Recommended Stack

The existing stack requires exactly one new production dependency: **KrylovKit.jl v0.10.2** (released October 2025). This pure-Julia package accepts any callable `f(x) -> y` as the linear operator — no `LinearMap` wrapper, no explicit matrix, no Fortran dependency. It provides `eigsolve` for right eigenvectors with `:LR` targeting, `bieigsolve` (BiArnoldi) for simultaneous left+right eigenvectors needed for overlap diagnostics, and `exponentiate` for matrix exponential action. One new development dependency is recommended: **TimerOutputs.jl v0.5** for hierarchical profiling of multi-minute Krylov runs where `@benchmark` is impractical.

All existing dependencies remain unchanged. `Arpack.jl` stays for the dense n<=6 path. `LinearMaps.jl` stays for backward compatibility (the commented-out `linearmaps_liouv.jl` prototype was abandoned as "slow" because it recomputed OFT per iteration without NUFFT prefactor caching; the new implementation avoids this). `FINUFFT.jl` and `SparseArrays` continue to be used inside the matvec for NUFFT prefactors and Bohr-frequency indexing.

**Core technologies:**
- **KrylovKit.jl v0.10.2**: Matrix-free Arnoldi eigensolving — chosen over ArnoldiMethod.jl (no `bieigsolve`) and Arpack.jl (no callable interface without explicit matrix); API verified via official docs
- **TimerOutputs.jl v0.5**: Hierarchical profiling for multi-minute Krylov runs — complements existing BenchmarkTools for microbenchmarks
- **OpenBLAS (Julia default)**: Multi-threaded BLAS-3 for dim x dim matrix multiplications in matvec; benchmark with physical core count before switching to MKL
- **Existing precomputation stack**: FINUFFT/OFT/Bohr dictionaries reused verbatim inside the matvec closure, providing the per-frequency Lindbladian action at O(dim^3) cost per frequency per jump

### Expected Features

**Must have (table stakes — the milestone does not exist without these):**
- **Matrix-free `apply_lindbladian!` for all four domains (Bohr, Energy, Time, Trotter)** — the keystone; everything else depends on it; must include both forward L(rho) and adjoint L'(rho) for bieigsolve
- **KrylovKit `eigsolve` wrapper with `:LR` targeting** — finds steady state + gap mode without shift-invert; request howmany=6-10 for robustness against clustered eigenvalues
- **`KrylovGapResult` struct** — eigenvalues, spectral gap, convergence count, matvec count, residual norms, fixed point density matrix
- **Cross-validation against dense `eigen()` at n=4,6** — the correctness gate; no production run at n>=8 until round-trip test and eigsolve validation both pass
- **Memory estimation pre-flight check** — `krylovdim * 4^n * 16 * 1.5` bytes before any run; guard against OOM

**Should have (differentiators):**
- **Both L (Lindbladian) and E (CPTP channel) formulations** — cross-validates continuous vs discrete; gap from L equals -log(|mu_2(E)|)/delta up to O(delta^2)
- **`bieigsolve` for left+right eigenvectors at n<=6** — extends existing biorthogonal overlap diagnostics (DIAG-05) beyond the dense regime
- **Timing benchmarks with 4^n scaling extrapolation** — empirical timing at n=4,6,8 extrapolated to n=10,12 as go/no-go resource estimate
- **Adaptive krylovdim** — auto-increase from 30 to 50 to 100 if `info.converged < howmany`; prevents silent failure from wrong default

**Defer to subsequent milestones:**
- n=10, n=12 production runs (require cluster hardware and empirical parameter tuning from n=8 results)
- Sector-resolved gap computation (requires symmetry projection infrastructure not yet built)
- GPU acceleration (data transfer overhead dominates many small dim x dim matrix multiplications)
- Arnoldi-Lindblad time-evolution approach (higher complexity; only needed if `:LR` targeting fails, which is unlikely)

### Architecture Approach

The architecture adds a single new file `src/krylov.jl` (~350 lines) without modifying any existing source file. The design mirrors the established domain-dispatch pattern from `jump_workers.jl`: a `_lindbladian_jump_action!` function dispatched on domain type computes L_k(rho) directly on dim x dim density matrices rather than assembling the dim^2 x dim^2 superoperator. A `KrylovWorkspace` struct pre-allocates all scratch matrices once; the closure passed to KrylovKit captures this workspace and produces zero allocation on the hot path except for the unavoidable `copy(vec(d_rho))` stored per Krylov iteration. Precomputation calls `_precompute_data()` and `_precompute_coherent_total_B()` identically to `construct_lindbladian()`, ensuring exact basis and domain consistency.

**Major components:**
1. **`KrylovWorkspace{T}` struct** — 5 pre-allocated dim x dim scratch matrices (d_rho, A_w, LdagL, tmp1, tmp2) reused across all Krylov iterations; eliminates the 256 MB/call allocation problem at n=12
2. **`_lindbladian_jump_action!` (4 domain-dispatched methods)** — per-jump Lindbladian contribution L_k(rho) in operator form; calls shared `_accumulate_dissipator!` helper for the canonical `A rho A' - 0.5{A'A, rho}` formula
3. **`build_lindbladian_action()` closure factory** — returns the `f(v) -> w` callable consumed by KrylovKit; captures precomputed data, jump list, config, and workspace by closure
4. **`krylov_spectral_gap()` public API** — validates inputs, precomputes domain data, calls `eigsolve`, sorts results by |Re(lambda)|, returns `KrylovGapResult`
5. **`KrylovGapResult` struct (in structs.jl)** — distinct from `LindbladianResult` (no stored matrix), `ExactDiagnosticsResult` (dense only), and `SpectralGapResult` (trajectory-based)

### Critical Pitfalls

1. **Vectorization convention mismatch** — implementing the matvec with transposed reshape conventions produces wrong eigenvalues that appear numerically plausible. Prevention: mandatory round-trip test `norm(L_dense * vec(rho) - vec(L_matvec(rho))) < 1e-12` for 10 random density matrices at n=4 before any eigsolve call. Use only Julia's native `vec()` and `reshape(v, dim, dim)`; never use `permutedims` or `transpose` in vectorization steps.

2. **Implementing `(E(rho) - rho) / delta` instead of direct Lindblad formula** — catastrophic cancellation loses ~2 digits per order of delta, limiting achievable Krylov tolerance. Prevention: compute L(rho) = sum_k (A_k rho A_k' - 0.5{A_k'A_k, rho}) directly from the dissipator formula, exactly as `_vectorize_liouv_diss_and_add!` does in vector form.

3. **Basis mismatch between Krylov matvec and dense diagnostics** — using `jump.data` (computational basis) instead of `jump.in_eigenbasis`, or forgetting TrotterDomain uses a different eigenbasis. Prevention: call `_precompute_data(config, ham_or_trott)` identically to `construct_lindbladian()`; validate both BohrDomain and TrotterDomain at n=4; never access `jump.data` inside the Krylov path.

4. **Allocating inside the Krylov closure** — any scratch allocation inside the matvec function scales to gigabytes per Krylov solve. Prevention: all workspace pre-allocated in `KrylovWorkspace` before `eigsolve`; the only acceptable allocation inside the closure is the unavoidable `copy(vec(ws.d_rho))`.

5. **Not checking `info.converged`** — KrylovKit returns results silently even when zero eigenvalues converged; `info.converged < howmany` is not an error, just a field. Prevention: always assert `info.converged >= howmany` and fall back to increased `krylovdim` or different starting vector on failure.

## Implications for Roadmap

The feature dependency graph from FEATURES.md and the domain-structured build order from ARCHITECTURE.md both point to the same four-phase progression: matvec correctness first, eigensolver integration second, domain extension third, production scaling fourth. Each phase delivers a testable artifact that validates the previous phase before adding complexity.

### Phase 1: Core Matvec Infrastructure (EnergyDomain)

**Rationale:** EnergyDomain has the simplest A_w computation (direct `oft!` call, no Bohr-bucket iteration, no NUFFT). Building it first isolates the core dissipator logic from domain-specific complexity. The round-trip test at n=4 is the foundational correctness gate and must exist before any eigensolver is involved. The historical failure in `linearmaps_liouv.jl` was caused by recomputing OFT per iteration without caching; the existing NUFFT prefactor cache eliminates that bottleneck.

**Delivers:** `KrylovWorkspace{T}` struct, `_accumulate_dissipator!` and `_accumulate_dissipator_adjoint!` shared helpers, `_lindbladian_jump_action!` for EnergyDomain, `build_lindbladian_action()` closure factory, round-trip test passing at n=4 EnergyDomain (with and without coherent term)

**Addresses:** Matrix-free `apply_lindbladian!` (P0), preallocated workspace pattern (prerequisite for all other phases)

**Avoids:** Pitfalls 3 (vectorization convention), 4 (L vs E confusion), 7 (precision loss) — all addressed by the round-trip test before any eigsolve is attempted

### Phase 2: KrylovKit Integration and Gap Validation

**Rationale:** With the EnergyDomain matvec correct, integrating KrylovKit's `eigsolve` with `:LR` targeting is the next sequential dependency. This phase produces the first Krylov spectral gap estimate and validates it against `extract_leading_eigendata()` at n=4 and n=6. The convergence checking infrastructure and result struct are locked in here.

**Delivers:** KrylovKit.jl added to Project.toml with `[compat]` pin, `krylov_spectral_gap()` public API, `KrylovGapResult` struct in structs.jl, memory estimation function, cross-validation against dense `eigen()` at n=4,6 EnergyDomain with correct eigenvalue sorting and phase-independent eigenvector comparison

**Uses:** KrylovKit.jl v0.10.2 (new production dependency), TimerOutputs.jl (new dev dependency for Krylov profiling)

**Implements:** Architecture components 3-5 (closure factory, public API, result struct)

**Avoids:** Pitfalls 1 (wrong targeting — `:LR` justified by Lindbladian spectrum structure), 5 (basis mismatch — validated against dense), 6 (non-convergence — howmany=6-10 with adaptive krylovdim), 10 (wrong comparison — sort-by-real-part and residual-norm validation)

### Phase 3: Domain Extension and Full Cross-Validation

**Rationale:** With the eigensolver integration validated on EnergyDomain, extending to TrotterDomain, TimeDomain, and BohrDomain follows the established pattern. TrotterDomain is the physically relevant domain for quantum algorithm simulation — it uses a different eigenbasis and must be validated explicitly from the beginning, not as an afterthought. BohrDomain has a distinct two-operator dissipator loop structure requiring its own implementation.

**Delivers:** `_lindbladian_jump_action!` for TrotterDomain and TimeDomain (NUFFT prefactor path), `_lindbladian_jump_action!` for BohrDomain (Bohr-bucket iteration with generalized two-operator dissipator), coherent term `-i[B, rho]` integration, cross-validation at n=4,6 for all four domains (both KMS and GNS balance types), `bieigsolve` integration for left+right eigenvectors at n<=6 (P2 feature, enables overlap analysis)

**Avoids:** Pitfall 5 (basis mismatch — TrotterDomain test is required, not optional)

### Phase 4: Production Scaling and Performance

**Rationale:** After all domains are validated at n<=6, n=8 is the first genuinely Krylov-only regime (the dense Lindbladian would be 64 GB). This phase establishes the empirical performance baseline, tunes BLAS threading, and produces timing data for extrapolating to n=10,12. It also adds the discrete CPTP channel formulation as an independent cross-check.

**Delivers:** n=8 Krylov spectral gap results for all domains, timing benchmarks with 4^n scaling fit, resource estimation report for n=10,12, BLAS thread optimization with restore pattern (mirroring existing `test_threading.jl`), discrete CPTP channel `E(rho)` formulation with `:LM` targeting as cross-check, eigenvalue condition number `kappa(lambda)` computation for error bound reporting

**Uses:** TimerOutputs.jl for hierarchical profiling; `BLAS.set_num_threads()` with `try-finally` restore

**Avoids:** Pitfalls 2 (memory blowup — pre-flight guard established), 8 (BLAS thread contention — thread count optimized), 11 (convergence criteria — kappa(lambda) computed and reported)

### Phase Ordering Rationale

- Phase 1 before 2: the linear map must be correct before any eigensolver is used. The round-trip test is a hard gate that cannot be deferred.
- Phase 2 before 3: validate the eigensolver on the simplest domain before adding domain complexity that could obscure bugs. A bug found at n=4 EnergyDomain costs seconds; the same bug found for TrotterDomain after integrating all four domains takes hours to isolate.
- Phase 3 before 4: all domains must work at small scale before expensive production runs. A TrotterDomain basis bug found at n=4 costs seconds; found at n=10 it costs hours of compute time.
- The `:LR` targeting decision (from STACK research) eliminates shift-invert entirely, which removes a full category of architectural complexity — no GMRES inner solve, no nested Krylov, no linear system infrastructure.
- The preallocated workspace pattern (from ARCHITECTURE research) ensures zero hot-path allocation, which is the dominant performance risk at n=12 where a single scratch allocation is 256 MB.

### Research Flags

Phases likely needing deeper research during planning:

- **Phase 4 (Production Scaling at n=10,12):** Performance characteristics are estimated from O(dim^3) scaling extrapolation, not measured. Actual wall-clock times per matvec, memory behavior under Julia GC pressure, and optimal BLAS thread counts on specific cluster hardware are all empirical questions. Phase 4 should be treated as a research+implementation phase — measure first at n=8, then decide on n=10,12 feasibility.
- **Phase 3 (BohrDomain two-operator dissipator):** The BohrDomain Lindbladian action is structurally distinct: it uses two different operators (alpha_A_nu1 and A_nu2) in a generalized dissipator `A1 rho A2' - 0.5 A2' A1 rho - 0.5 rho A2' A1`, matching the existing `_vectorize_liouv_diss_and_add!(L, jump_1, jump_2_dag, scalar, ws)` call. The Krylov-side implementation of this generalized dissipator requires careful analysis of the Bohr-bucket loop structure in `jump_workers.jl` lines 11-45 before coding.

Phases with well-documented patterns (skip additional research):

- **Phase 1 (Core Matvec):** The Lindbladian dissipator formula is textbook physics. The `mul!`-based implementation is fully specified. The round-trip test design is unambiguous. No additional research needed.
- **Phase 2 (KrylovKit Integration):** The `eigsolve` API is thoroughly documented with clear examples. The `:LR` targeting decision is unambiguous given the Lindbladian spectrum structure (all Re(lambda) <= 0 means the gap eigenvalue is extremal under `:LR`). KrylovKit version is pinned to v0.10.2.

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Stack | HIGH | KrylovKit v0.10.2 API verified via official docs (October 2025 generation). `:LR` targeting justified mathematically from Lindbladian spectrum structure, confirmed by Arnoldi-Lindblad paper (peer-reviewed). Alternative packages (ArnoldiMethod.jl, Arpack.jl) explicitly evaluated and rejected with documented reasons. |
| Features | HIGH | Core features derived from direct codebase analysis (26 source files) and verified KrylovKit API. Feature prioritization follows hard technical dependencies (matvec before eigsolve; validation before production). Performance estimates for n=10,12 are LOW confidence (extrapolated from O(dim^3) scaling, not measured). |
| Architecture | HIGH | Component boundaries chosen via direct analysis of all existing source files. Domain dispatch pattern mirrors `jump_workers.jl` exactly. Zero modifications to existing files (except `structs.jl` addition and `QuantumFurnace.jl` includes). Historical failure mode in `linearmaps_liouv.jl` fully explained and addressed. |
| Pitfalls | HIGH | Most pitfalls identified from codebase analysis (existing BLAS thread management in `test_threading.jl`, abandoned prototype in `linearmaps_liouv.jl`) and KrylovKit GitHub issues (#9 memory overhead, #23/#38 degeneracy failures). Low-confidence estimates are explicitly flagged (timing at large n). |

**Overall confidence:** HIGH

### Gaps to Address

- **Timing at n=10,12:** The O(dim^3) extrapolation suggests hours per Krylov solve at n=12, but the actual cost depends on number of energy labels per jump, BLAS efficiency on specific hardware, and GC behavior. Empirical timing at n=8 (Phase 4) must be completed before committing to n=10,12 runs. Treat Phase 4 as empirical calibration, not just implementation.

- **KrylovKit memory overhead at scale:** GitHub issue #9 (2018) reports ~20x overhead vs Arpack for small problems; the 1.5x safety factor in the budget estimates may be too optimistic or too pessimistic for the large-vector regime at n=12. Validate with `@allocated` wrapping the eigsolve call at n=8 before scaling further.

- **Degenerate eigenvalue behavior at n>=8:** Systems with SZ conservation have eigenvalue multiplets whose density increases with system size. The PITFALLS research identifies this as a moderate risk (Pitfall 6) but the actual multiplet structure at n=8-12 for the target Heisenberg chain is unknown. The existing `detect_multiplets` infrastructure from `diagnostics.jl` could be adapted to predict convergence difficulty before a long run.

- **GNS balance type coverage:** FEATURES.md's cross-validation matrix includes GNS but ARCHITECTURE.md focuses on KMS examples. GNS uses a different `transition` function but the same domain dispatch structure. Should be straightforward but requires explicit test coverage in Phase 3.

## Sources

### Primary (HIGH confidence)

- [KrylovKit.jl official documentation](https://jutho.github.io/KrylovKit.jl/stable/) — API reference, `eigsolve`/`bieigsolve` signatures, algorithm parameters; generated October 10, 2025
- [KrylovKit.jl GitHub releases](https://github.com/Jutho/KrylovKit.jl/releases) — v0.10.2, October 11, 2025; confirms version and API stability
- [KrylovKit.jl eigenvalue problems](https://jutho.github.io/KrylovKit.jl/stable/man/eig/) — `:LR`/`:SR`/`:LM` targeting semantics; explicit no-shift-invert documentation
- [KrylovKit.jl algorithms reference](https://jutho.github.io/KrylovKit.jl/stable/man/algorithms/) — Arnoldi defaults (krylovdim=30, maxiter=100, tol=1e-12), BiArnoldi for bieigsolve
- QuantumFurnace.jl codebase (all 26 source files) — direct analysis for integration points, existing patterns, historical failure modes
- [ArnoldiMethod.jl documentation v0.4.0](https://julialinearalgebra.github.io/ArnoldiMethod.jl/dev/) — evaluated and rejected (no `bieigsolve`, no `exponentiate`)

### Secondary (MEDIUM confidence)

- [Arnoldi-Lindblad time evolution (Huybrechts & Roscilde, Quantum 2022)](https://quantum-journal.org/papers/q-2022-02-10-649/) — validates `:LR` approach for Lindbladian spectrum; peer-reviewed
- [Tensor Network Framework for Lindbladian Spectra (arXiv:2509.07709)](https://arxiv.org/abs/2509.07709) — complex-time Krylov methods for Lindbladian eigenvalues
- [KrylovKit.jl GitHub Issues #9, #23, #38](https://github.com/Jutho/KrylovKit.jl/issues/9) — memory allocation overhead and degenerate eigenvalue behavior; 2018 report, may not reflect v0.10 behavior
- [Krylov.jl performance tips](https://jso.dev/Krylov.jl/dev/tips/) — BLAS threading recommendations: use physical cores
- [QuantumToolbox.jl steadystate](https://qutip.org/QuantumToolbox.jl/stable/users_guide/steadystate) — validates `:SR`/`:LR` approach for Lindbladian steady state; different codebase but same physics
- [Three approaches for representing Lindblad dynamics (arXiv:1510.08634)](https://arxiv.org/pdf/1510.08634) — vectorization conventions documentation

### Tertiary (LOW confidence, needs empirical validation)

- O(dim^3) matvec timing extrapolations to n=10,12 — formula-based, not measured; requires n=8 benchmarks to calibrate
- n=12 full Krylov memory budget (~27 GB) — 1.5x KrylovKit overhead factor may be inaccurate for large-vector regime
- Krylov dimension requirements at n>=8 — system-dependent; empirical tuning required during Phase 4

---
*Research completed: 2026-02-20*
*Ready for roadmap: yes*
