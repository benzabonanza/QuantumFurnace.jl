# DLL Lindbladian — code+test audit (qf-fzj.1)

Read-only audit of the DLL Lindbladian surface (filters + Bohr/Time
operators + Kossakowski) covering Phase 51 (qf-3i8), the qf-hur perf
refactor, and the qf-wmg Metropolis filter. Source paths are absolute.

## 1. Public API surface

### Exports from `src/QuantumFurnace.jl`

| Group | Lines | Names |
|---|---|---|
| Filters (DLL-1) | `src/QuantumFurnace.jl:84-86` | `AbstractFilter`, `GaussianFilter`, `DLLGaussianFilter`, `DLLMetropolisFilter`, `time_kernel`, `freq_kernel`, `filter_time_cutoff` |
| DLL dissipator (DLL-2) | `src/QuantumFurnace.jl:88-89` | `dll_lindblad_op_bohr`, `dll_lindblad_op_time` |
| DLL coherent (DLL-3) | `src/QuantumFurnace.jl:91-92` | `dll_coherent_op_bohr`, `dll_coherent_op_time`, `dll_coherent_kernel_bohr` |
| DLL Kossakowski (DLL-4) | `src/QuantumFurnace.jl:94-95` | `dll_kossakowski_bohr` |

### `Config.filter` plumbing

- Field: `Config.filter::Union{Nothing, AbstractFilter} = nothing` at
  `src/structs.jl:126`. Docstring at lines `74-78`.
- Resolver: `_resolve_filter(config)` at `src/misc_tools.jl:117-126`.
  `nothing` → `GaussianFilter(config.sigma)`.
- Consumers in `src/`:
  - `src/furnace_utensils.jl:111` — DLL Bohr precompute
  - `src/furnace_utensils.jl:138, 143, 165, 185` — DLL Time and CKG
    Time/Trotter precomputes (calls `_truncate_time_labels_for_oft` and
    `_prepare_oft_nufft_prefactors`).
- Validation in `validate_config!`:
  - `src/misc_tools.jl:188-192` — DLL Gaussian β-match
  - `src/misc_tools.jl:194-203` — DLL Metropolis β-match + S>0
  - `src/misc_tools.jl:206-222` — DLL construction requires explicit
    filter; rejects `EnergyDomain` / `TrotterDomain` (deferred)
  - `src/misc_tools.jl:264-271` — special `_collect_config_errors!`
    overload for `Config{<:Any, TimeDomain, DLL}` (no `w0` requirement)

### Public-facing entry points

The two public DLL workflows are:

1. `Lindbladian` + `BohrDomain`+`DLL` → `_precompute_data` at
   `src/furnace_utensils.jl:107-112`, then `_jump_contribution!` at
   `src/jump_workers.jl:109-127`.
2. `Lindbladian` + `TimeDomain`+`DLL` → `_precompute_data` at
   `src/furnace_utensils.jl:134-158` (FINUFFT amortisation), then
   `_jump_contribution!` at `src/jump_workers.jl:150-169`.
3. Coherent term wired at `src/coherent.jl:52-70` for both domains.

`Thermalize` / `Trajectory` paths have **no DLL specialisation**
(checked: no DLL methods in `_precompute_coherent_unitary` at
`src/coherent.jl:88-140`). Per beads notes, DLL is currently
Lindbladian-only — not a bug, but worth noting.

## 2. Internal helpers

### Filters (`src/filters.jl`)

| Helper | file:line | Signature | Callers | Note |
|---|---|---|---|---|
| `q_weight(::DLLGaussianFilter, ν)` | `filters.jl:132` | `(f, nu) -> T` | NONE outside test (`test_dll_filter.jl`) | Documented as "DB weight before KMS factor"; only used in tests for round-trip equality `freq_kernel == q_weight·e^{-βν/4}`. |
| `_time_oft_prefactor_gaussian(::GaussianFilter)` | `filters.jl:108` | `(f) -> T` | `filter_time_cutoff(::GaussianFilter)` at line 119; `test_dll_filter.jl:107` | Single-call helper. |
| `_time_oft_prefactor_dll(::DLLGaussianFilter)` | `filters.jl:171` | `(f) -> T` | `filter_time_cutoff(::DLLGaussianFilter)` at line 182 | **Single-call helper, not used in tests** — could be inlined into `filter_time_cutoff`. |
| `_hormander_eta(t)` | `filters.jl:212` | `(t::T) -> T` | `_hormander_phi` at lines 229-230 | Internal stack. |
| `_hormander_phi(t)` | `filters.jl:223` | `(t::T) -> T` | `_hormander_bump` at line 248 | Internal stack. |
| `_hormander_bump(x)` | `filters.jl:241/252` | `(x::Real) -> T` | `q_weight(::DLLMetropolisFilter)` at `filters.jl:329`; `test_dll_filter.jl` (sanity tests) | Used. |
| `q_weight(::DLLMetropolisFilter, ν)` | `filters.jl:328` | `(f, nu) -> T` | `freq_kernel(::DLLMetropolisFilter)` at line 346; tests | Used. |

### DLL operators (`src/dll.jl`)

| Helper | file:line | Signature | Callers | Note |
|---|---|---|---|---|
| `_dll_J_quadrature(beta, s; kwargs)` | `dll.jl:239` | `(β, s) -> T` | `_dll_J_table` at line 266 | Internal stack. |
| `_dll_g_closed_form(beta, t, tp, J_val)` | `dll.jl:247` | `(β, t, t', J) -> Complex` | `dll_coherent_op_time(::DLLGaussianFilter)` at line 327 | **Single-call helper.** |
| `_dll_J_table(beta, time_labels; kwargs)` | `dll.jl:257` | `(β, ts) -> (J_vals, sum_to_index)` | `dll_coherent_op_time(::DLLGaussianFilter)` at line 322 | **Single-call helper.** |
| `_dll_coherent_from_g_tt(jumps, ham, g_tt, time_labels, τ)` | `dll.jl:466` | `(...) -> Matrix{CT}` | both `dll_coherent_op_time` overloads (lines 333 + 444) | **Two callers (Gaussian + Metropolis)** — confirms qf-wmg.6 refactor goal. |
| `dll_coherent_op_time_legacy(...)` | `dll.jl:542` | reference path | `test_dll_coherent.jl:198` only | **Reference-only**, kept per `feedback_fast_path_default.md`. |
| `_resolve_filter(config)` | `misc_tools.jl:125` | `(Config) -> AbstractFilter` | `furnace_utensils.jl:111, 138, 165, 185` (4 call sites) | Used. |

**Single-call helpers worth flagging** (not actionable, just visibility):

- `_time_oft_prefactor_dll` — only used by its own filter's
  `filter_time_cutoff`; could be a `let` or inline math.
- `_dll_g_closed_form` and `_dll_J_table` — both 1-call, but each is
  meaningfully named and improves readability of
  `dll_coherent_op_time(::DLLGaussianFilter)`. Recommend keeping.
- `_dll_J_quadrature` is a 1-call internal of `_dll_J_table` but is
  the one that actually invokes QuadGK — keep.

## 3. Dedup / simplification candidates

### Filter API minimality

Required interface (per `filters.jl:23-26` docstring): `time_kernel`,
`freq_kernel`, `filter_time_cutoff`. Plus `eltype` for the NUFFT path.

| Method | `GaussianFilter` | `DLLGaussianFilter` | `DLLMetropolisFilter` |
|---|---|---|---|
| `Base.eltype` | line 65 | line 66 | line 315 |
| `time_kernel` | line 88 | line 157 | line 373 |
| `freq_kernel` | line 98 | line 143 | line 345 |
| `filter_time_cutoff` | line 118 | line 181 | line 403 |
| `q_weight` | — (CKG has no q) | line 132 | line 328 |
| `_time_oft_prefactor_*` | line 108 | line 171 | — (Metropolis cutoff uses doubling search) |

The `AbstractFilter` interface is the smallest one fitting all three.
`q_weight` is **not part of `AbstractFilter`** — it's a DLL-only
introspection helper, defined only on the two DLL types and
non-exported. No legacy bloat on `DLLGaussianFilter`.

### `_dll_coherent_from_g_tt` reuse

Confirmed:
- `dll_coherent_op_time(::DLLGaussianFilter)` calls it at `dll.jl:333`.
- `dll_coherent_op_time(::DLLMetropolisFilter)` calls it at `dll.jl:444`.

Both are filter-agnostic from "Step 3" onward, exactly as the qf-wmg.6
refactor intended. **No duplication.**

### Default fast path

`Config()` defaults: `filter = nothing`, which `_resolve_filter` maps
to `GaussianFilter(config.sigma)`. The CKG path then routes via the
NUFFT prefactor stack in `src/nufft.jl::_prepare_oft_nufft_prefactors`
and `_jump_contribution!` for `Time/TrotterDomain`.

For DLL, `Config(...; filter=DLLGaussianFilter(beta))` routes:

- BohrDomain → `dll_lindblad_op_bohr` (closed-form, n³)
- TimeDomain → `oft_nufft_at_zero` slice (precompute) + `*=` in
  `_jump_contribution!` (`jump_workers.jl:165-167`)
- TimeDomain coherent → `_dll_coherent_from_g_tt` via NUFFT
- Metropolis TimeDomain coherent → 2D NUFFT over `[-S, S]²` then
  `_dll_coherent_from_g_tt`

`dll_coherent_op_time_legacy` (the O(Nt²·n³) reference at
`dll.jl:542`) is **never reachable** from a `Config(...)` default. Only
`test_dll_coherent.jl:198` calls it directly. Per
`feedback_fast_path_default.md` this is correct: legacy stays as a
test-only reference.

`dll_lindblad_op_time` (the explicit Riemann sum at `dll.jl:73`) is
**also only reachable from tests** (`test_dll_dissipator.jl:267`).
Production uses the precomputed `oft_nufft_at_zero` slice.

Both `dll_lindblad_op_time` and `dll_coherent_op_time_legacy` are
public-facing reference paths; consider noting in their docstrings
that they are not on the production path. (Currently
`dll_lindblad_op_time`'s docstring describes it as the canonical
Eq. 3.4 third-form computation without flagging that production
amortises via FINUFFT.)

### Bohr ↔ Time cross-checks

| File:line | Test name | Filter | β | Tol | What it compares |
|---|---|---|---|---|---|
| `test_dll_coherent.jl:101` | (c) Bohr ↔ Time agreement on G | DLLGaussian | {1,5,10} | `≤ 1e-3` | `G_b vs G_t` (G alone) |
| `test_dll_coherent.jl:270` | (k) Metropolis Bohr ↔ Time on G | DLLMetropolis | {1,5,10} | `≤ 1e-5` | `G_b vs G_t` (G alone) |
| `test_dll_coherent.jl:190` | (h) closed-form NUFFT == legacy | DLLGaussian | {1,5,10} | `≤ 1e-8` | NUFFT path vs legacy 2D-grid path |
| `test_dll_dissipator.jl:117` | (c) Bohr ↔ Time consistency | DLLGaussian | {1,5,10} | `≤ 1e-4` | full Liouvillian, n=2 toy |
| `test_dll_dissipator.jl:231` | (h) NUFFT slice == explicit Riemann | DLLGaussian | {1,5,10} | `≤ 1e-12` | jump-only operator (FINUFFT precision) |
| `test_dll_kms_db.jl:261` | (i) DLL Metropolis Bohr ↔ Time L | DLLMetropolis | {1,5,10} | `≤ 1e-4` | full Liouvillian, n=3 |
| `test_dll_kms_db.jl:290` | (j) Bohr ↔ Time converges with t0 | DLLMetropolis | {5, 10} | `≤ 1e-7` | grid-convergence sweep |

**Consolidation candidate.** The four "Bohr ↔ Time agreement" tests
(coherent.c, coherent.k, dissipator.c, kms_db.i) test **the same
identity** at different tolerances on different fixtures (n=2 toy vs
n=3 Heisenberg) and different filters (Gaussian vs Metropolis). They
could be a single parameterised testset with fixture/filter/tol
matrix:

```
(filter, fixture, β, what, tol) ∈
  {(Gaussian,   n=2-toy,    {1,5,10}, "L",          1e-4),
   (Gaussian,   n=3-heis,   {1,5,10}, "G alone",    1e-3),
   (Metropolis, n=3-heis,   {1,5,10}, "G alone",    1e-5),
   (Metropolis, n=3-heis,   {1,5,10}, "L",          1e-4)}
```

Recommend collapsing into one parameterised testset
`@testset "Bohr ↔ Time agreement (DLL)"` either in
`test_dll_dissipator.jl` (for L) or split G/L. The grid-convergence
test (kms_db.j) should stay separate — it's a structural witness
test, not an agreement test.

### KMS skew-symmetry checks

Asserted in:

- `test_dll_kossakowski.jl:136-153` — Eq. 4.7, DLLGaussian, β∈{1,5,10}
- `test_dll_kossakowski.jl:158-176` — Eq. 4.8 centrosymmetry
  (equivalent), DLLGaussian
- `test_dll_kossakowski.jl:239-257` — Eq. 4.7, DLLMetropolis,
  β∈{1,5,10}

Plus an **indirect KMS-DBC ratio check** in `test_dll_filter.jl:474`:
`f̂(ν)/f̂(-ν) = e^{-βν/2}` on the Metropolis flat top.

The Gaussian (ll. 136-153) and Metropolis (ll. 239-257) variants are
**near-identical code** with the filter type swapped. Could be a
single helper `function _check_kms_skew_symmetry(filter, ν_grid; β,
atol)` reused across both. Centrosymmetry (157-176) is
mathematically equivalent to skew-symmetry — testing both is
redundant but cheap; flag for consolidation but not removal.

The filter-level KMS-DBC ratio check at `test_dll_filter.jl:474` is
a different witness (single-ν, no Kossakowski) and stays.

## 4. Test redundancies

### `test/test_dll_filter.jl` (28 testsets, 603 lines)

Useful sweeps (β ∈ {1, 5, 10} per memory):

- `(h.2-h.7) q_weight, freq_kernel, support, symmetry` — sweep over
  `(β, S) ∈ {(5,2), (1,3), (10,1.5)}` — covers the `(β=10, β=1)`
  bracket plus a custom S; **good**.
- `(h.8) Metropolis asymptote` — uses `β=40, S=2.0` (single-point)
  for clean asymptotic; **fine, structural**.
- `(i.4) filter_time_cutoff bounds` — `β ∈ {1, 2, 5, 10}` × `tol ∈
  {1e-6, 1e-9, 1e-12}`; **good**.

Pure duplicates / cheap overlaps:

- **Letter collision**: `(g) DLL NUFFT prefactor spot check` at line
  157 and `(g) Hörmander bump` at line 335 share the letter "g";
  similarly `(h) validate_config! DLL beta mismatch` at line 210 and
  `(h) DLLMetropolisFilter freq_kernel + q_weight` at line 383.
  Cosmetic only — not a duplicate, but the section labels in the
  outer `@testset` are visually confusing. Suggest renaming the
  Metropolis block to `(j..)` and `(k..)`, or restructuring as nested
  Gaussian/Metropolis groups.
- `(h.5) bump invisible on flat top` and `(h.2) q_weight` overlap on
  the flat-top region — both check `q_weight(ν) == bare-u(βν)` on
  `|ν| ≤ S/2`. (h.5) is structural, (h.2) is endpoint coverage —
  keep both.
- `(c) Fourier round-trip via QuadGK` (line 40) and `(i.2) Fourier
  round-trip (time → freq)` (line 549) check the same identity for
  Gaussian/DLLGaussian and Metropolis respectively; OK as
  per-filter coverage.

### `test/test_dll_dissipator.jl` (8 testsets, 277 lines)

- `(a-d)` four tests on the same n=2 toy fixture sweeping β∈{1,5,10}.
  Each tests a distinct property (Gibbs FP, Bohr FP, Bohr↔Time
  agreement, dual trace), so they are not duplicates — but they all
  rebuild `_make_dll_config(domain; β)` and `_build_dll_toy_system(β)`
  per β. Could share a fixture cache.
- `(e-g)` validate_config! rejects: 3 distinct rejection paths
  (filter=nothing, TrotterDomain, EnergyDomain). Each tests a
  different error message; **keep**.
- `(h)` NUFFT slice vs explicit Riemann sum: `≤ 1e-12` per-jump
  comparison on n=3 fixture, β∈{1,5,10}. Structural test (FINUFFT
  precision floor); **keep**.

The fixture used in (a-d) is the **2-qubit XXZ toy**, while (h) uses
the n=3 disordered Heisenberg. Per `feedback_n3_minimum_test_size.md`
"tests should use n=3 minimum, not n=2", the toy is below the
recommended minimum. **Flag for migration**: the toy was chosen
because it's the only fixture where the bare DLL dissipator (without
G) preserves σ_β at machine precision (lines 5-12 docstring). Once G
is wired in (which it is in `_precompute_coherent_B`), the n=3 fixture
should work at the same tolerance — the toy may now be obsolete.

### `test/test_dll_coherent.jl` (13 testsets, 300 lines)

- `(a-d)` G isolation tests on n=3 Heisenberg, β∈{1,5,10}. Each
  testset sweeps β, no inner duplication.
- `(e-g)` single-jump Lindbladian tests, all at β=5 except (g)
  which is the β-sweep variant — (e) and (f) at fixed β are mostly
  consistency checks; **keep**.
- `(h)` legacy reference cross-check, β-sweep — structural, keep.
- `(i-l)` Metropolis variants of (a-d). **Mostly** parallel to the
  Gaussian (a-d):
  - (i) ↔ (a) — Hermiticity Bohr; same tolerance (`≤ 1e-15`)
  - (j) ↔ (b) — Hermiticity Time; same tolerance (`≤ 1e-6`)
  - (k) ↔ (c) — Bohr ↔ Time; tighter for Metropolis (`≤ 1e-5` vs `1e-3`)
  - (l) ↔ (d) — norm bounded; **different intent**: Gaussian (d)
    asserts `≤ 1.0`, Metropolis (l) asserts `1e-4 ≤ ‖G‖ ≤ 10` (the
    "stays O(1)" qualitative claim)

  The Hermiticity tests (i)/(a), (j)/(b) **could be merged** into a
  single parameterised testset over filter type. Same for
  Bohr↔Time agreement (c)/(k) — see §3 above.

### `test/test_dll_kms_db.jl` (10 testsets, 335 lines)

- `(a, c)` Bohr DLL Gaussian KMS-DB exact; β-sweep. Two distinct
  facts (relative norm + spectral gap match); **keep**.
- `(b)` Time DLL Gaussian KMS-DB up to quadrature; β-sweep.
- `(d, e)` β=10 fixed: dissipator-only fails, CKG also passes.
  Cross-construction comparison; **keep**.
- `(f, g)` Metropolis variants of (a, c) and (e); **same shape as
  Gaussian (a,c,e) on the Metropolis filter**. The Metropolis (f)
  combines (a)+(c) — slightly more compact than Gaussian.
- `(h)` Time DLL Metropolis KMS-DB up to quadrature — Metropolis
  variant of (b); β-sweep. **Direct duplicate of (b) with filter
  swapped.**
- `(i)` full Bohr↔Time Liouvillian agreement, Metropolis only —
  why no Gaussian variant? `test_dll_dissipator.jl:117` (c)
  partially covers Gaussian L Bohr↔Time on the n=2 toy at
  `≤ 1e-4`; this (i) is on n=3 at `≤ 1e-4`. Separate fixture, so
  not a duplicate, but consider parameterising.
- `(j)` t0 grid convergence sweep — structural, keep.

**Pattern**: the Gaussian/Metropolis split runs through 3-4 testsets
in each of `coherent.jl` (a/i, b/j, c/k, d/l) and `kms_db.jl` (a/f,
b/h, c/{f-tail}). A single per-property parameterised loop over
`filter ∈ {DLLGaussian, DLLMetropolis}` would cut duplication ~40%.

### `test/test_dll_kossakowski.jl` (12 testsets, 287 lines)

- `(a-c)` shape/Hermitian/PSD/rank-1/outer-product. β-sweep.
- `(d)` HamHam overload — single fixture; structural, keep.
- `(e)` CKG vs DLL rank gap — single β; structural, keep.
- `(f, g)` KMS skew-symmetry + centrosymmetry — β-sweep. Are
  mathematically equivalent (g follows from f); both checked at
  `≤ 1e-12`. Keep one as primary, demote the other to a
  consistency check or merge.
- `(h)` β-scaling — β-sweep. Structural.
- `(h2-h5)` Metropolis variants of (a, c, f) plus a qualitative
  contrast. **Same Gaussian/Metropolis duplication pattern as
  coherent.jl/kms_db.jl** — merge into per-property
  parameterised testsets.

## 5. Smells / surprises

### Filter methods: paper-faithful

I checked each filter method against the docstring's paper reference;
all match.

- `freq_kernel(::DLLGaussianFilter)` at `filters.jl:143` is
  `e^{1/8} exp(-(βν+1)²/8)` — Eq. 3.22.
- `time_kernel(::DLLGaussianFilter)` at `filters.jl:157` is the
  closed-form inverse FT — Eq. 3.3 / 3.22.
- `freq_kernel(::DLLMetropolisFilter)` at `filters.jl:345` is
  `q · e^{-βν/4}` — Eq. 3.20.
- `q_weight(::DLLMetropolisFilter)` at `filters.jl:328` is
  `exp(-√(1+(βν)²)/4) · w(ν/S)` — Eq. 3.19.
- `dll_coherent_kernel_bohr` at `dll.jl:136` is `(1/2i)·tanh·f̂·conj(f̂)`
  — Eq. 3.5.
- The Hörmander bump at `filters.jl:241` matches Assumption 15 (the
  standard `η/(η+η(1-·))` construction).

No method contradicts the paper.

### Dead code

None of the helpers I checked are dead — every `_dll_*` is reached
from at least the `dll_coherent_op_time` Gaussian or Metropolis
overloads.

`dll_coherent_kernel_bohr` is **exported** but used only by the
non-exported `dll_coherent_op_time_legacy` (test-only reference).
Either:
- demote to private (`_dll_coherent_kernel_bohr`) and keep it
  alongside the legacy reference, or
- remove from exports if nothing in `drafts/` / scripts depends on
  the public name.

`dll_lindblad_op_time` is **exported** but in production the time
NUFFT slice path bypasses it. The function is used only by tests
(`test_dll_dissipator.jl:267`). Safe to keep exported (it is the
canonical Eq. 3.4 third-form for documentation), but the docstring
should mention "production amortises this via FINUFFT in
`_precompute_data`" so future readers don't think this is the hot
path.

### Asymmetric β / ν tolerances

- `test_dll_coherent.jl:110` "`≤ 1e-3`" Bohr↔Time G (Gaussian) is
  loose; `test_dll_kms_db.jl:269` "`≤ 1e-4`" full L Bohr↔Time
  (Metropolis) is tighter. Both at default grid (Nt=4096). Likely
  fine — G alone has a different error structure than L.
- `test_dll_dissipator.jl:109` `(b) TimeDomain L[σ_β] ≤ 1e-4` and
  `test_dll_dissipator.jl:139` `(d) L†[I] ≤ 1e-4` — same tolerance.
  (`(a)` and `(c)` at `1e-10` are Bohr-domain — these are exact.)

No surprising asymmetry across β values within a single test;
tolerances are uniform across `_BETAS = (1.0, 5.0, 10.0)`. Where
β=10 is the worst case, the test bound is set to accommodate it.

### `S` parameter consistency

- Defined: `DLLMetropolisFilter{T<:AbstractFloat}(beta, S)` at
  `filters.jl:305-308`; convenience constructor at lines 310-311
  takes `S::Real = T(2)`.
- Used at:
  - `filters.jl:329` (`q_weight`): `_hormander_bump(nu / f.S)`
  - `filters.jl:375` (`time_kernel` integration limits)
  - `filters.jl:404` (`filter_time_cutoff` oscillation period
    `4π/S`)
  - `dll.jl:381, 389-390` (Metropolis 2D-NUFFT integration grid)
- Validated: `validate_config!` at `misc_tools.jl:200` only checks
  `S > 0`.

**Surprise / smell.** The docstring at `filters.jl:293-298` claims:

> `S = 2` keyword default. This must satisfy `S/2 ≥ max|ν_BH|` for
> the chosen Hamiltonian, otherwise the bump bites the Lindbladian
> (only the central region is Metropolis; outside `[-S/2, S/2]` the
> asymptote `min{1, exp(-βν/2)}` does not hold). **`validate_config!`
> emits a warning when this is violated.**

This is **not implemented** — `validate_config!` at
`misc_tools.jl:194-203` only checks `S > 0`. The `S/2 ≥ max|ν_BH|`
check requires the Hamiltonian and is therefore not currently
performed. Either:
- update the docstring to say "the caller must ensure `S/2 ≥
  max|ν_BH|`; otherwise Eq. 3.20 asymptote silently breaks", or
- add the check in `_precompute_data` where the Hamiltonian
  becomes available (would require a warning, not an error).

This is a documentation drift, not a bug. The default `S=2` is
correct for the test fixtures (`max|ν_BH| ≤ 0.9`).

### TrotterDomain DLL

`validate_config!` at `src/misc_tools.jl:217-221` rejects
`TrotterDomain` DLL with the message
`"DLL construction in TrotterDomain is deferred — not yet supported."`

I searched all DLL paths in `src/coherent.jl`, `src/jump_workers.jl`,
`src/dll.jl`, `src/filters.jl`: **no Trotter-specific DLL method
exists**. The dispatch on `Config{<:Any, TimeDomain, DLL}` (e.g.
`coherent.jl:65`, `jump_workers.jl:154`) does **not** match
`TrotterDomain`, so a stray DLL+TrotterDomain config can only reach
the runtime path via `validate_config!` failing. **No leak.**

The single `_collect_config_errors!` overload for
`Config{<:Any, TimeDomain, DLL}` at `misc_tools.jl:264-271` is the
only DLL-aware domain check (correctly skipping `w0`); the generic
`TrotterDomain` overload at `misc_tools.jl:273` would still apply
if an internal call ever reached it, but `validate_config!` rejects
DLL+TrotterDomain first at line 220.

### Other surprises

- `dll_coherent_op_bohr(jumps, hamiltonian, filter, beta)` accepts
  `beta` even though `filter.beta` is the source of truth. Docstring
  notes this at `dll.jl:171-174`: *"The `beta` argument is accepted
  for symmetry with `dll_coherent_op_time` but the kernel reads `β`
  from the filter; mismatch is the caller's responsibility (validated
  upstream by `validate_config!`)."* `validate_config!` enforces
  `filter.beta == config.beta`. Worth a passing thought: the explicit
  `beta` arg is now redundant — could be removed for clarity, or
  wrapped to assert `beta == filter.beta` defensively. Touches
  `coherent.jl:59` (call site).
- `dll_coherent_op_time(::DLLMetropolisFilter)` builds a fresh FINUFFT
  plan twice per call (`dll.jl:428` for the ν-grid → time grid, and
  `_dll_coherent_from_g_tt` builds another at `dll.jl:504`). At Nt
  ≈ 4096 the plan-construction overhead is negligible (10s of µs)
  but worth noting: a single shared planner could amortise.
- `_dll_J_table` builds a `Dict{T, Int}` keyed by `Float64` sums
  (`dll.jl:268`). On uniformly spaced `t0` grids the sums `t_m + t_n`
  are exactly representable and the dict lookup is safe; on a
  non-uniform grid this would be a hash-collision risk. Currently
  fine because production uses the uniform grid from
  `_precompute_labels` / `_truncate_time_labels_for_oft`. Not a bug,
  but a hidden assumption.
- `dll_coherent_op_time` (Gaussian path at `dll.jl:307`) uses
  `filter::DLLGaussianFilter` (concrete). The Metropolis overload at
  `dll.jl:369` uses `filter::DLLMetropolisFilter{T}`. **The legacy
  reference** `dll_coherent_op_time_legacy` at `dll.jl:546` uses
  `filter::AbstractFilter` (broadest) — could in theory be called
  with a CKG `GaussianFilter` and silently produce nonsense (the
  KMS factor would be missing). Document that legacy is DLL-filter-only
  or restrict the type.
- `Config.filter` is annotated `Union{Nothing, AbstractFilter}` —
  this works but precludes specialisation through `Config`. Fine
  because `_resolve_filter` resolves before any hot path; just
  noting this for completeness.

---

# Summary

The DLL surface is paper-faithful and well-isolated: the `_dll_*`
helpers are all reached, `_dll_coherent_from_g_tt` is correctly
shared between Gaussian and Metropolis, and the fast NUFFT path is
the production default with the legacy `O(Nt²·n³)` reference
gated to test-only. Main cleanup opportunities are (1) merging the
Gaussian/Metropolis test pairs into parameterised testsets (saves
~40% in `coherent.jl`, `kms_db.jl`, `kossakowski.jl`), (2) migrating
the n=2 toy fixture in `test_dll_dissipator.jl` to n=3 per the
project's minimum-test-size rule, and (3) fixing the
`DLLMetropolisFilter` docstring to match the actual `validate_config!`
behaviour (no `S/2 ≥ max|ν_BH|` check is performed).
