# QuantumFurnace.jl Memory

## Completed Work

### Phase 43: Bi-exponential fitting (completed 2026-03-04)
- `BiexpFitResult` struct and `fit_biexponential_decay` in `src/fitting.jl`
- `estimate_mixing_time` extended with `model=:biexp` keyword in `src/mixing.jl`
- Extrapolation via `Roots.Bisection` (no closed-form for multi-exponential)
- Result: <0.001% error on synthetic data (vs 0.13% single-exp, vs 26% on real data)
- 1273 tests pass (1246 existing + 27 new)

## Key Architecture

- Fitting: `src/fitting.jl` (LsqFit.jl v0.15.x, Levenberg-Marquardt)
- Mixing estimation: `src/mixing.jl` (post-processing wrapper)
- Module: `src/QuantumFurnace.jl` — already has `using Roots`, `using LsqFit`
- Structs: `FitResult`, `BiexpFitResult`, `MixingTimeEstimate` (all immutable)
- Tests: `test/test_fitting.jl`, `test/test_mixing.jl`
- Diagnostic scripts: `scripts/diagnose_floor_v3.jl`, `scripts/mixing_time_extrapolate_verify.jl`

## Findings from Diagnostics (2026-03-03)
- Floor scales linearly with delta (exponent ~1.0 with coherent unitary): floor ≈ 0.068 × delta
- Coherent unitary is essential — without it, floor saturates at ~0.0017
- Single-exp overestimates offset C by ~47% → 26% error in extrapolated mixing time
- Root cause: sensitivity of `t = -ln((ε-C)/A)/gap` when ε ≈ C

## TrotterDomain Floor Analysis (2026-03-04)

### Two-component floor model: `floor = k_energy × δ + floor_Trotter(constant)`
- EnergyDomain k (scales with δ): n=3→0.059, n=4→0.094, n=5→0.142
- TrotterDomain constant floor: n=3→~8e-6, n=4→~8e-5 (10× larger!), n=5→~3.5e-5
- **n=4 "anomaly" explained**: at δ=0.0001, Trotter constant dominates → effective k≈0.879
- Trotter error in OFT/R-matrix is NOT the cause (~1e-8 for all n)
- Files: `src/errors.jl`, `scripts/trotter_error_comparison.jl`, `scripts/floor_superoperator.jl`
- Required δ to reach ε=1e-4: n=3→0.001, n=4→0.0001, n=5→0.0005

### Bohr Frequency Collision Root Cause (2026-03-04)
- **Even-n periodic chains are bipartite** → sublattice operator U=diag(1,-1,1,-1,...) exists
- U anti-commutes with hopping (XX+YY), commutes with on-site Z disorder
- In single-magnon sectors (Sz=±(n/2-1)): ZZ diagonal is constant → U creates spectral symmetry
- Result: E_k(Sz=+m) + E_{d+1-k}(Sz=-m) = C (exact cross-sector pairing, spread ~5e-16)
- This gives exact Bohr frequency collisions (gap reversal palindrome)
- n=4: 12 collisions (all from Sz=±1, single-magnon), n=6: 30 (all from Sz=±2)
- Odd n: ring is frustrated (not bipartite) → no U → no collisions
- Multi-magnon sectors: ZZ diagonal varies → symmetry breaks (no collisions from Sz=0)
- Scripts: `scripts/bohr_collision_sectors.jl`, `scripts/zz_disorder_symmetry_break.jl`

### ZZ Bond Disorder Symmetry Breaking Experiment (2026-03-04)
- Adding ZZ bond disorder (Σ_q ε_q Z_q Z_{q+1}) on top of existing Z disorder breaks bipartite pairing
- **Collisions eliminated**: n=4: 12→0, n=6: 30→0. Cross-sector pair sum spread 5e-16 → 1.4e-03
- ZZ-only disorder (without Z) fails: P=∏X_i commutes with all bilinear terms → Sz=±m degenerate → nu_min=0
- Need BOTH Z (breaks P) + ZZ (breaks bipartite pairing)
- **Trotter mismatch resolved (2026-03-04)**: HamHam now supports multiple disordering terms natively

## Bi-exp Verification Results

### Z-only disordering (2026-03-04, legacy)
- n=3 (δ=0.001): 1.35% error, PASS
- n=4 (δ=0.0001): 0.78% error, PASS (fit reports converged=false, still works)
- n=5 (δ=0.0005): 0.29% error, PASS
- Scripts: `scripts/biexp_mixing_verify.jl` (EnergyDomain), `scripts/biexp_mixing_verify_trotter.jl` (TrotterDomain)

### Z+ZZ disordering (2026-03-04, symmetry-breaking)
- All n values give <2% error when extrap time is sufficient
- δ=0.0005 for all n; floor ~3-6e-5, well below target 1e-4
- n=3: 0.22%, n=4: 1.3-1.6%, n=5: 0.12% (EnergyDomain)
- **Key finding**: bi-exp fit needs data ~3× past crossing to separate decay from floor
- Previous n=4 "failure" was simply insufficient extrapolation time (T=60 vs crossing at t≈59)
- Collisions broken for all n (min Bohr gap > 0)
- No qualitative even/odd difference — n=4 anomaly eliminated
- extrap_time_map: n=3→100, n=4→150, n=5→150
- Script: `scripts/biexp_mixing_verify_zzdisordering.jl`

## Hamiltonian Construction Notes
- `find_ideal_heisenberg`: optimizes disorder for max min-gap, `disordering_terms` kwarg (default `[[Z]]`)
- `HamHam` struct: `disordering_terms::Union{Vector{Vector{Matrix{Complex{T}}}}, Nothing}` — multiple terms
- `disordering_coeffs::Union{Vector{Vector{T}}, Nothing}` — per-term per-site coefficients
- Constructors: (1) no disorder, (2) multi-term, (2b) single-term convenience, (3) from NamedTuple
- NamedTuple constructor handles both legacy `disordering_term` and new `disordering_terms` keys
- `_construct_disordering_terms` has overloads for single-term and multi-term
- `_trotterize2` iterates over all disordering terms (each gets its own Trotter layer)
- `pad_term` handles 2-site terms with periodic wrapping
- Trotter: even-n uses 2 bond groups (clean), odd-n needs 3rd group for wrapping bond

## Thesis
- [Thesis Structure and Status](thesis_structure.md) — MSc thesis structure, completed/empty sections, original contributions (Props 5,7,9,10)
- [Thesis Numerics Plan](thesis_numerics_plan.md) — Full plan for Ch5 numerical data: plots, cross-checks of all analytical results, priority order
- [GQSP angle-finding decision](gqsp_angle_finding_decision.md) — Final: use GQSP (Motlagh-Wiebe) + Berntson-Sünderhauf FFT for Q, not optimization
- [Trotter sign convention mismatch](thesis_sign_convention_trotter.md) — 1_preliminaries.tex mixes $\ee^{+\ii Ht}$ (QPE) with $\ee^{-\ii Ht}$ (Trotter subsection); 2_methods.tex uses $+$. Unify in a later pass
- [TODO: bring R_b reflection into 2_methods.tex GQSP](todo_methods_gqsp_reflection.md) — prelim uses W = R_b U_B explicitly; methods chapter still treats U_B alone as the walk. Also resolve sin/cos convention.
- [Thesis Review chapter — quantum vs classical Gibbs sampling](thesis_review_chapter.md) — beads epic `qf-yk4`; framing: where could quantum Gibbs sampling first beat classical? Regime map + /lit-review corpus.

## Conventions
- [Circuits folder stays .tex-only](feedback_circuits_folder.md) — never leave .aux/.log/.pdf/.DS_Store in `drafts/circuits/`
- [Thesis .tex are read-only; drafts are readable Markdown](feedback_thesis_tex_readonly.md) — never Edit/Write `.tex` in `supplementary-informations/`; drafts go to `drafts/<slug>.md` as Markdown with inline `$math$`, not LaTeX dumps. Canonical example: `drafts/quantum-circuits-basics.md`
- [Plot PNGs may fail Read with API 400](feedback_image_read_fails.md) — verify generated plots via stdout diagnostics or ask the user; do not retry the same PNG after a "Could not process image" error

## Plotting
- [Thesis colour palette](reference_thesis_colors.md) — named colours (pinegreen, bordeaux, dustyplum, deepplum, aubergine, slateblue, sage, ochre, terracotta, dustyteal, mustard); refer to by name in plotting requests
- [Thesis gradient palettes](reference_gradient_palettes.md) — cold (mint→navy), warm (cream→mulberry), diverging (teal↔purple) gradients for heatmaps and continuous-scale plots; prefer over `:inferno`/`:viridis`

## Memory Persistence
- Memory lives in `.claude-memory/` at repo root, symlinked to `~/.claude/projects/.../memory/`
- Setup via `SessionStart` hook in `.claude/settings.json` (runs once per session)
- Migrated from `.planning/memory/` + `PreToolUse` hook on 2026-03-05
