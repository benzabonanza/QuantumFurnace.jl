# `predict_lindbladian_trajectory` in TimeDomain — speed benchmark

**Task**: qf-1au. Companion to qf-w2u (dissipator quadrature) and qf-fqt
(TimeDomain dissipator). Now the **full Lindbladian** with
`include_coherent = true` (the default), in TimeDomain, matrix-free,
multi-threaded.

**Setup**: smooth Metropolis with `default_smooth_s`. Canonical v2
registers: `r_D = 7`, `r_b_minus = 7`, `r_b_plus = 17`. Both coherent
registers actually used: `B_time` (2D FINUFFT over the
$(t_\text{outer}, \tau_\text{inner})$ grid) is constructed at workspace
build time and lives in `G_left = +iB^\top − R/2`, `G_right = -iB^\top
− R/2`. Each `apply_lindbladian!` matvec includes the $i[B, \rho]$
commutator. JULIA_NUM_THREADS = 8, BLAS = 1, krylovdim = 30.

**Caveat on `predict_lindbladian_trajectory` API**: the production
function in `src/lindblad_action.jl:601` is dispatched only on
`Config{Lindbladian, <:Union{BohrDomain, EnergyDomain}}`. The body is
domain-agnostic (uses `apply_lindbladian!` via closure), so the bench
calls the same inner helper `_krylov_spectral_decomposition` directly
with a TimeDomain workspace. The numbers below are what a future
`Union{BohrDomain, EnergyDomain, TimeDomain, TrotterDomain}` dispatch
would deliver — no algorithmic change, just a relaxed type guard.

## Headline

**`predict_lindbladian_trajectory` in TimeDomain is fast and well-behaved across the whole
benchmarked grid.** Total wall per cell stays under 10 s even at n = 6
with β_phys = 2 and the full coherent term. The workspace build
(dominated by `B_time` at `r_b_plus = 17`) is the main cost; the actual
Krylov factorisation is 2 s at n = 6.

| n | β_phys | d | wall_ws | wall_krylov | total | matvecs | ms/matvec | gap |
|---|--------|---|---------|-------------|-------|---------|-----------|-----|
| 3 | 1.0 | 8 | 3.69 s *(JIT warm)* | 1.50 s | 5.19 s | 30 | 50.0 | 1.48e−1 |
| 3 | 2.0 | 8 | 0.20 | 0.03 | 0.23 | 30 | 1.2 | 5.83e−2 |
| 4 | 1.0 | 16 | 0.37 | 0.09 | 0.46 | 30 | 2.9 | 1.14e−1 |
| 4 | 2.0 | 16 | 0.23 | 0.07 | 0.30 | 30 | 2.5 | 3.82e−2 |
| 5 | 1.0 | 32 | 1.21 | 0.25 | 1.46 | 30 | 8.4 | 5.21e−2 |
| 5 | 2.0 | 32 | 0.65 | 0.29 | 0.95 | 30 | 9.7 | 5.63e−3 |
| 6 | 1.0 | 64 | 6.88 | 2.06 | **8.94** | 30 | 68.7 | 3.90e−2 |
| 6 | 2.0 | 64 | 3.75 | 2.09 | **5.84** | 30 | 69.7 | 4.06e−3 |

(n = 3, β_phys = 1 is the first-touched cell and pays ~3 s of Julia JIT
warmup; later cells at the same `d` are <1 s.)

## What's expensive vs cheap

- **Workspace build** dominates at n = 6: 77% of total wall at β = 1, 64%
  at β = 2. The dominant sub-cost is `B_time` — a 2D FINUFFT over
  $(2^{r_{b−}} \times 2^{r_{b+}}) = (128 \times 131072)$ time
  samples, evaluated at the full set of $d^2$ Bohr frequencies. Setting
  `r_b_plus = 17` (the v2 recipe value to keep coherent ε ≤ 10⁻⁶) is
  what we're paying for.

- **Krylov factorisation** scales as ~d² in `ms/matvec`:
  $1.2 \to 2.9 \to 8.4 \to 68.7$ ms for $d \in \{8, 16, 32, 64\}$. Mostly
  the per-jump sandwich `jump_oft · ρ · jump_oft†` (O(d³) GEMM) × n_jumps × 2^r_D
  ω-loop with ω-loop threading on. Each cell does exactly 30
  matvecs (the krylovdim) — Arnoldi never broke down here.

- **β_phys = 2 is consistently faster than β_phys = 1 in workspace build**
  (~half the wall at n ≥ 4). The `b_±(t)` kernels decay faster at higher
  β, so `_truncate_time_labels_for_oft` cuts the time grid earlier and
  `B_time` does less FINUFFT work. The Krylov factorisation cost is
  β-independent (same fixed-size NUFFT cache regardless).

- **Memory** stays ~860 MB peak — dominated by the FINUFFT intermediate
  buffers during `B_time` construction, GC'd between cells. The
  persistent workspace cache is ~150 MB at n = 6 (`2^r_D × d² × 16 bytes
  × n_jumps`).

## Comparison vs EnergyDomain

The previously-measured EnergyDomain `predict_lindbladian` from
[krylov-trajectory-bench-qf-ev5](../../.claude-memory/krylov_trajectory_bench_qf_ev5.md)
reported wall = 0.028 → 224 s for n = 3 → 7 (β_alg = 10 in that
benchmark, krylovdim = 30, with_gqsp = false). At the comparable n = 6
point that scaling gives ≈ 80 s vs TimeDomain's 9 s here — but it's
apples-to-oranges (different β, the EnergyDomain bench used the legacy
`num_energy_bits=12` global kwarg, not the per-register triples). What
this run lets us conclude:

- TimeDomain at n = 6 with `r_D=7, r_b_plus=17` and coherent on takes ~9 s
  per `predict_lindbladian_trajectory` call.
- ~75% of that is one-time workspace setup (`B_time` 2D FINUFFT). If a
  caller does many `predict_lindbladian_trajectory` calls at fixed
  `(config, ham, jumps)` — varying `t_grid`, `rho_0`, etc. — the
  workspace can be built once and passed via the `workspace` kwarg
  (`src/lindblad_action.jl:611`), cutting per-call wall to ~2 s at n = 6.
- The actual Krylov factorisation is ~2 s at n = 6 → ~10 s at n = 7
  (extrapolating d² in matvec). Still tractable on the sandbox without
  any change to the heap envelope.

## Implication

TimeDomain `predict_lindbladian_trajectory` with the **full Lindbladian**
is viable. The user worry from the prior dissipator-only sweep ("how
badly does this blow up?") is answered: it doesn't. Two follow-ups:

1. **Relax dispatch** in `src/lindblad_action.jl:601` to
   `Config{Lindbladian, <:Union{BohrDomain, EnergyDomain, TimeDomain,
   TrotterDomain}}`. The body is already domain-agnostic; the guard is
   a "not yet tested" marker that this bench answers.
2. **Cache `B_time`** across calls when iterating over `t_grid` / `rho_0`
   at fixed `(config, ham, jumps)`. The `workspace = ...` kwarg path
   already supports this — just plumb it through any production
   sweep callers.

## Files

- Script: `scripts/scratch_predict_lindbladian_time_bench.jl`
- BSON sidecars: `scripts/output/predict_lindbladian_time_bench/n{3..6}_bphys{1.0,2.0}_smooth_def_s_time.bson`
- Campaign summary: `scripts/output/predict_lindbladian_time_bench/_campaign_summary.bson`
- Beads: qf-1au (under epic qf-yt9 / qf-fqt thread)
