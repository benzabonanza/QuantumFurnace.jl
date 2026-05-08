---
name: DLL Metropolis-type filter (qf-wmg, completed 2026-05-01)
description: Full Ding–Li–Lin Metropolis filter (Eq. 3.19–3.20) added alongside the existing Gaussian filter; BohrDomain + TimeDomain end-to-end with 2D NUFFT for the coherent G; 257 new tests, all passing
type: project
---

# DLL Metropolis-type filter — beads epic `qf-wmg` (10 atomic commits)

Closed 2026-05-01. Full implementation of the Metropolis-type Gevrey
filter from Ding–Li–Lin 2024 (Eq. 3.19–3.20) alongside the existing
`DLLGaussianFilter` (Eq. 3.21–3.22).

## Why the Metropolis filter

`f̂(ν) ≈ min{1, exp(-βν/2)}` on the flat top — the smoothed Metropolis
acceptance. The qualitative win over the Gaussian filter is that
`|α(ν, ν')| = |f̂(ν) f̂(ν')*|` stays **O(1)** at low T (β large) on the
negative-ν half of the flat top, whereas the Gaussian's Kossakowski entries
shrink ∝ `1/β` around the centre frequency `-1/β`. At β=10, n=3 fixture,
ν = -0.45: Metropolis 0.95 vs Gaussian 0.06 — **15× larger** for Metropolis.

## Mathematical content

```
u(x)   = exp(-√(1+x²)/4)                                  (Eq. 3.19)
q(ν)   = u(βν) · w(ν/S)                                   (Eq. 3.19)
f̂(ν)  = q(ν) · exp(-βν/4) ≈ min{1, exp(-βν/2)} on |ν| ≤ S/2  (Eq. 3.20)
```

`w` is the Hörmander mollifier: `w ≡ 1` on `|x| ≤ 1/2`, smooth decay to
0 at `|x| = 1`, identically 0 outside. Default `S = 2` puts the flat top
at `[-1, 1]` so the Bohr frequencies of our test fixtures
(`max|ν_BH| ≤ 0.9`) are covered and the bump is invisible to the
Lindbladian on the Bohr grid.

## Atomic commits

| Task   | Commit  | Topic                                                      |
|--------|---------|------------------------------------------------------------|
| qf-wmg.1  | 55e5d55 | Hörmander bump utility `_hormander_bump` + 22 tests        |
| qf-wmg.2  | 2e7e049 | `DLLMetropolisFilter` struct + `freq_kernel` + `q_weight`  |
| qf-wmg.3  | 6727311 | `time_kernel` via QuadGK + `filter_time_cutoff` (bisection)|
| qf-wmg.4  | 310e80a | Export + `validate_config!` routing (β + S checks)         |
| qf-wmg.5  | 498a13d | BohrDomain Gibbs FP + KMS skew + Kossakowski qual tests    |
| qf-wmg.6  | bef4206 | Refactor: extract `_dll_coherent_from_g_tt` helper         |
| qf-wmg.7  | 2485a27 | Coherent G TimeDomain via 2D NUFFT (Metropolis path)       |
| qf-wmg.8  | 9891002 | TimeDomain L + Bohr↔Time cross-check (correctness gate)    |
| qf-wmg.9  | 502594a | Diagnostic: Kossakowski heatmap Metropolis vs Gaussian     |
| qf-wmg.10 | (this)  | Memory + close epic                                        |

## Architecture

- `src/filters.jl`: adds `DLLMetropolisFilter{T}(beta; S=2)` next to `DLLGaussianFilter`. Hörmander bump utilities. Numerical `time_kernel` via QuadGK + doubling-search `filter_time_cutoff`.
- `src/dll.jl`: extends `dll_coherent_kernel_bohr` dispatch to the Union of both DLL filters. Adds `dll_coherent_op_time(::DLLMetropolisFilter, ...)` that tabulates `ĝ(ν, ν')` on `[-S, S]²` (default `Nν = 256`) + 2D type-3 NUFFT to time grid + calls the new `_dll_coherent_from_g_tt` helper. The helper is a pure refactor that the Gaussian path also uses.
- `src/misc_tools.jl`: extends `validate_config!` with `DLLMetropolisFilter.beta == Config.beta` and `S > 0` checks.
- `src/QuantumFurnace.jl`: exports `DLLMetropolisFilter`.

## Test coverage

Total +257 tests (full DLL suite: 1648/1648 passing):
- `test_dll_filter.jl`: 164 (struct/freq_kernel/q_weight) + 43 (time_kernel/cutoff) + 5 (validate_config) + 22 (Hörmander bump) = 234
- `test_dll_kossakowski.jl`: +183 (Metropolis subtests h2-h5)
- `test_dll_coherent.jl`: +12 (Metropolis subtests i-l)
- `test_dll_kms_db.jl`: +24 (Metropolis subtests f-i)

Key invariants verified at β ∈ {1, 5, 10}, n=3 disordered Heisenberg:
- KMS-DBC skew-symmetry α(ν,ν') = α(-ν',-ν) e^{-β(ν+ν')/2} (Eq. 4.7)
- Bohr↔Time agreement on G alone: ‖G_b − G_t‖ ≤ 1e-4 (Nν = 256)
- Bohr↔Time agreement on full L: ‖L_b − L_t‖ ≤ 1e-3
- TimeDomain Metropolis: KMS-DB up to quadrature (relative_norm ≤ 1e-3)

## Diagnostic script

`scripts/plot_kossakowski_metropolis_vs_gaussian.jl` produces a 2×3 grid
of log₁₀|α(ν,ν')| heatmaps (rows = filters, columns = β ∈ {1,5,10})
saved to `drafts/plots/kossakowski_metropolis_vs_gaussian.png`.

## Out of scope (deferred / not implemented)

- TrotterDomain DLL Metropolis (mirrors current `DLLGaussianFilter` scope)
- Custom user-supplied weight functions (a `DLLCustomFilter` would need a separate filter type)
- Closed-form `g(t, t')` for Metropolis (the (u, v) substitution does not separate)

## Files touched

- `src/filters.jl` (Metropolis struct + Hörmander)
- `src/dll.jl` (coherent G TimeDomain Metropolis + helper extraction)
- `src/misc_tools.jl` (validate_config!)
- `src/QuantumFurnace.jl` (export)
- `test/test_dll_filter.jl`, `test_dll_kossakowski.jl`, `test_dll_coherent.jl`, `test_dll_kms_db.jl`
- `scripts/plot_kossakowski_metropolis_vs_gaussian.jl` (diagnostic)
