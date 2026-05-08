# GQSP audit (qf-fzj.1)

Read-only audit of the GQSP integration in QuantumFurnace.jl. All references are
absolute file:line. No source modified.

## 1. Public API surface

`src/QuantumFurnace.jl` does NOT export any GQSP symbol — neither the helpers
`_gqsp_block_encoding_alpha` / `_gqsp_apply_polynomial` (correct: they are `_`-prefixed
private helpers) nor the convenience block-encoding norm. The whole feature is reachable
through `Config` flags only.

Public-facing knobs (both fields of the exported `Config` struct):

- `Config.with_gqsp::Bool` — `src/structs.jl:122` (default `false`).
- `Config.gqsp_degree::Int` — `src/structs.jl:123` (default `1`).

Docstring entries describing them: `src/structs.jl:67-72`. The (non-exported) helpers
themselves live in `src/coherent.jl` and are reached via `QuantumFurnace.<name>` in tests
and scratch scripts.

Note: `Config` is exported at `src/QuantumFurnace.jl:63`. Any cleanup that renames or
removes the GQSP fields is therefore a public-facing API change.

## 2. Internal helpers

Every GQSP-related private helper lives in `src/coherent.jl`. The simulator's "α-helper"
plus the polynomial evaluator are the only GQSP-specific functions; the rest of the
machinery (`B_time`, `B_trotter`, `_compute_b_minus`, `_compute_b_plus*`,
`_compute_truncated_func`) is shared with the matrix-exp path.

| Helper | File:line | Signature | Call sites in `src/` |
|---|---|---|---|
| `_gqsp_block_encoding_alpha` | `src/coherent.jl:257-268` | `(jump::JumpOp, b_minus, b_plus, t0_sim::Real, gamma_norm_factor::Real) -> Real` | 3 (TimeDomain branch `coherent.jl:109`, TrotterDomain branch `coherent.jl:123`, trajectory branch `trajectories.jl:103-106`) |
| `_gqsp_apply_polynomial` | `src/coherent.jl:288-348` | `(B::AbstractMatrix{<:Complex}, alpha::Real, delta::Real, d::Int) -> Matrix` | 3 (same three sites: `coherent.jl:110`, `coherent.jl:124`, `trajectories.jl:107`) |

Both helpers have **3 call sites in src/**, none of which are single-call. Test-only
direct uses live in `test/test_gqsp_polynomial.jl` (qualified as `QuantumFurnace._gqsp_*`)
and one scratch script `scripts/scratch_gqsp_trace_drift.jl`.

**No single-call helpers to inline.** The functions are used three times each across two
files — extracting them was the right call. (See "Dedup" below for whether the three
call sites should themselves be merged into one.)

Worth noting (not flagged for inlining, just for context): `_gqsp_apply_polynomial`
internally has a `d == 1` fast path (`coherent.jl:302-312`) and a `d ≥ 2` Clenshaw branch
(`coherent.jl:317-347`). Both paths are exercised by tests. Default `gqsp_degree=1` makes
the fast path the production case; the Clenshaw branch is used for `gqsp_degree ≥ 2`
which only the regression test `test/test_gqsp_thermalize.jl:108-120` exercises.

## 3. Dedup / simplification candidates

### 3.1 The three `with_gqsp` branches are near-clones

The TimeDomain and TrotterDomain branches inside `_precompute_coherent_unitary`
(`src/coherent.jl:103-114` and `src/coherent.jl:116-128`) do **structurally identical**
work after building `B`: compute `α_be`, dispatch on `with_gqsp`, write `U_terms[k]`. The
two `with_gqsp` arms differ only in the third argument to `_gqsp_block_encoding_alpha`:
`config.t0` vs `trotter.t0`. Diff:

```
coherent.jl:108-113  if config.with_gqsp
                         α_be = _gqsp_block_encoding_alpha(jump, b_minus, b_plus, config.t0,  gamma_norm_factor)
                         U_terms[k] = _gqsp_apply_polynomial(B, α_be, delta, config.gqsp_degree)
                     else
                         U_terms[k] = exp(-1im * delta * Hermitian(B))
                     end

coherent.jl:122-127  if config.with_gqsp
                         α_be = _gqsp_block_encoding_alpha(jump, b_minus, b_plus, trotter.t0, gamma_norm_factor)
                         U_terms[k] = _gqsp_apply_polynomial(B, α_be, delta, config.gqsp_degree)
                     else
                         U_terms[k] = exp(-1im * delta * Hermitian(B))
                     end
```

Plus the trajectories.jl counterpart (`src/trajectories.jl:101-110`):

```
trajectories.jl:101-110  if config.with_gqsp
                             t0_sim = config.domain isa TrotterDomain ? ham_or_trott.t0 : config.t0
                             α_be = _gqsp_block_encoding_alpha(jumps[a], precomputed_data.b_minus, precomputed_data.b_plus,
                                                               t0_sim, precomputed_data.gamma_norm_factor)
                             per_op_U_B[a] = _gqsp_apply_polynomial(B_a, α_be, delta_eff, config.gqsp_degree)
                         else
                             per_op_U_B[a] = exp(-1im * delta_eff * Hermitian(B_a))
                         end
```

The same 5-line `if config.with_gqsp ... else ...` pattern (build α, apply polynomial,
else `exp(Hermitian)`) appears 3× in 2 files. **Cleanup candidate:** extract to a helper
of the form

```julia
# pseudocode signature
_coherent_unitary_step(jump, B, b_minus, b_plus, t0_sim, gamma_norm_factor, delta, config) -> Matrix
```

that internally branches on `config.with_gqsp`. The trajectories.jl branch already
computes `t0_sim` explicitly (`trajectories.jl:102`); applying the same pattern to
`_precompute_coherent_unitary` is what the helper should encapsulate — and it removes
the special-casing of `config.t0` vs `trotter.t0` from the caller.

### 3.2 α-helper logic — single source already, but its callers re-pluck precomputed-data

`_gqsp_block_encoding_alpha` itself is in one place (`coherent.jl:257`). Good.
What is duplicated is the **destructure+call** at the three call sites: each one
unpacks `b_minus`, `b_plus`, `gamma_norm_factor`, picks a `t0_sim`, and calls. Folding
this into the helper proposed in §3.1 removes that triplication.

### 3.3 Polynomial construction (Jacobi-Anger / Chebyshev) — no rebuild, but coefficient evaluation could be hoisted

`_gqsp_apply_polynomial` recomputes the Bessel coefficients `J_k(δα)` and the phase factor
`cis(-π/2 · k)` **on every call** (`coherent.jl:299, 303, 330`). For
`_precompute_coherent_unitary`, this is fine — called once per jump per simulation. For
`_build_trajectory_workspace` (`trajectories.jl:97-111`), also once per jump.

**Not flagged for cleanup:** the per-call cost is `O(d)` Bessel evaluations and is
dominated by the matmuls. A coefficient cache would only matter if the polynomial were
re-evaluated inside the trajectory hot loop, which it is not (per-op `U_B` is precomputed
into `per_op_U_B[a]` once at workspace-build time).

### 3.4 `trajectories.jl` GQSP branch — does it duplicate `coherent.jl` logic?

Yes, structurally — see §3.1. Functionally, the trajectories.jl branch differs in *what*
`B` it operates on:

- `_precompute_coherent_unitary` (`coherent.jl:88`) calls `B_time`/`B_trotter` directly
  to build a single-jump `B_k` in the eigenbasis (lines 106, 120) and immediately
  multiplies by `gamma_norm_factor`.
- `_build_trajectory_workspace` (`trajectories.jl:99`) goes through
  `_precompute_coherent_B` (the polymorphic helper that already wraps `B_time` /
  `B_trotter` / `B_bohr`), which **also** multiplies by `gamma_norm_factor` internally.

Both arrive at the same Hermitian `B` for the GQSP polynomial; the trajectories.jl path
also calls `hermitianize!(B_a)` (`trajectories.jl:100`), but `_precompute_coherent_unitary`
does NOT (`coherent.jl:106-107`, `coherent.jl:120-121`). See §5 for the smell.

### 3.5 `_compute_truncated_func` ℓ¹-norm path

`_gqsp_block_encoding_alpha` does `sum(abs, values(b_minus))` /
`sum(abs, values(b_plus))` (`coherent.jl:264-265`) on every call. The truncated `Dict`s
are built once per simulation at `furnace_utensils.jl:172-175`. The ℓ¹ sums are
recomputed from scratch for every jump (3, 9, 12 times per simulation depending on
`num_qubits`). Tiny constant-factor optimisation: cache the two ℓ¹ norms next to
`b_minus`/`b_plus` in the precomputed-data NamedTuple. Not flagged as a structural
duplication (single call-site definition), but worth a one-liner in the cleanup phase
if the helper is rewritten.

## 4. Test redundancies

### `test/test_gqsp_config.jl`

Top-level testset `GQSP config fields and validation (qf-63j.1)`
(`test_gqsp_config.jl:15`). Helper `_make_cfg` parameterises `(construction, domain,
with_gqsp, gqsp_degree)`. Five inner testsets:

| line | name | property asserted |
|---|---|---|
| 30 | Defaults: `with_gqsp=false`, `gqsp_degree=1` | default values + non-throw |
| 37 | Accepts `with_gqsp=true` in supported regime | KMS+TimeDomain, KMS+TrotterDomain, d=2, d=100 |
| 47 | Rejects with GNS | GNS construction throws |
| 51 | Rejects outside Time/TrotterDomain | EnergyDomain, BohrDomain throw |
| 56 | Rejects `gqsp_degree` out of [1,100] | 0, -3, 101 throw |

No duplication. Each branch of `validate_config!` (`misc_tools.jl:170-183`) corresponds
to exactly one assertion. **Useful sweep**, not redundant: domain sweep at line 39-44
exercises each accepted combination once.

### `test/test_gqsp_polynomial.jl`

Four top-level testsets, ten inner testsets total.

`@testset "GQSP polynomial _gqsp_apply_polynomial (qf-63j.2)"` at `test_gqsp_polynomial.jl:78`:

| line | name | property asserted |
|---|---|---|
| 80 | d=1 closed form vs hand formula | `f_1 = J_0 I − 2i J_1 (B/α)` exactly, 3 δ values |
| 95 | d=2 closed form vs hand formula | `f_2 = (J_0+2J_2) I − 2i J_1 (B/α) − 4 J_2 (B/α)²`, 3 δ values |
| 114 | Reduces to identity at δ=0 | `J_0(0)=1`, `J_k(0)=0` for any d, α |
| 126 | Slope-(d+1) δ-scaling | finite-difference slope ≈ d+1 for d∈{1,2,3} |
| 147 | Joint (B,α,δ)-invariance | `f_d` depends only on `(B/α, δα)`; (cB, cα, δ/c) |
| 166 | Tail-bound scaling | doubling α inflates error by ~2^(d+1) for d∈{1,2,3} |

`@testset "GQSP block-encoding norm _gqsp_block_encoding_alpha (qf-63j.2)"`
(`test_gqsp_polynomial.jl:187`): one block — synthetic `b_minus`/`b_plus` Dict fixture,
checks formula + 4× scaling under jump-doubling.

`@testset "GQSP operator/circuit equivalence: f_d(B/α) ≈ [L_d(W)]_anc=0"`
(`test_gqsp_polynomial.jl:211`): builds qubitization walk via `_qubit_block_encoding` /
`_qubit_walk` (locally defined ports of `scripts/scratch_gqsp_random_h.jl`), runs through
`d∈{1,2,3} × δ∈{1e-1, 1e-2}` = **6 cases**.

`@testset "n=3 Heisenberg slope-2 anchor (qf-63j.2)"` (`test_gqsp_polynomial.jl:235`):
end-to-end with `B_time` from the actual n=3 fixture, slope-2 check.

**Redundancy assessment:**
- The `d=1` and `d=2` closed-form testsets at lines 80 and 95 are useful as anchor
  formulas (they pin the sign convention and the `2(-i)^k` coefficient explicitly).
  **Not duplicate.**
- "Slope-(d+1)" (line 126) and "Tail-bound scaling" (line 166) test the **same**
  Bessel-tail estimate `‖f_d − exp(-iδB)‖ ∼ (δα)^{d+1}` on a random Hermitian `B` with
  `‖B‖ = 1`, just rotating which axis (δ or α) is varied. These are **two views of one
  phenomenon**. Not pure duplication (the tail-bound version is a stronger
  multiplicative constraint), but the cleanup pass could merge them into one
  parametrised `@testset` with `for axis in (:delta, :alpha)`.
- "Joint (B,α,δ)-invariance" (line 147) is **structural**: it pins the rescaling group
  symmetry. Not redundant with the slope tests.
- "n=3 Heisenberg anchor" (line 235) and "Slope-(d+1)" (line 126) **both** test slope-2
  at d=1, on different `B` (random Hermitian vs `B_time` for n=3 Heisenberg). The
  anchor test is meaningful as a regression (uses the simulator's truncated-func dicts
  + `B_time`), so it's a **useful sweep** — but the redundant inner d=1 slope check at
  line 134-144 (`d in (1, 2, 3)`) already covers d=1. **Mark for review:** is the n=3
  anchor adding signal beyond what slope-(d+1) at d=1 already gives? Yes, marginally —
  it stresses `B_time` end-to-end with real grid + truncated `Dict`s, not a synthetic
  `B`. Keep, but it's the weakest of the four top-level testsets.

### `test/test_gqsp_thermalize.jl`

Single top-level testset `GQSP integration in run_thermalize / run_trajectory (qf-63j.3)`
(`test_gqsp_thermalize.jl:23`). Seven inner testsets:

| line | name | property asserted |
|---|---|---|
| 44 | TimeDomain smoke + GQSP | `run_thermalize` returns valid `ThermalizeResults`, finite trace dist, tr≈1, hermiticity |
| 55 | TrotterDomain smoke + GQSP | same as above on TrotterDomain |
| 66 | TimeDomain regression vs matrix-exp | trace_distances match, ‖final_dm‖ within 5e-3 |
| 77 | TrotterDomain regression vs matrix-exp | same on TrotterDomain |
| 87 | Thermalization quality (longer T_mix=5.0) | both paths thermalize to td<0.15 AND |td_gqsp − td_exp| < 1e-3 |
| 108 | Higher gqsp_degree improves agreement | `err_d2 < err_d1` strictly |
| 122 | Trajectory regression vs matrix-exp | `run_trajectory` n=1, ‖rho_mean diff‖<5e-3 |

**Redundancy assessment:**
- Lines 44 and 55 (smoke tests on TimeDomain vs TrotterDomain) are **useful sweep**
  across the two GQSP-supported domains. Not duplicate.
- Lines 66 and 77 (regression vs matrix-exp on TimeDomain vs TrotterDomain) are also a
  **useful sweep**. Not duplicate.
- The pair of smoke (line 44) + regression (line 66) on TimeDomain is **partial
  duplication**: both run `run_thermalize` with `with_gqsp=true` on the same n=3 Heisenberg
  fixture at `delta=1e-3, mixing_time=0.05`. The smoke test asserts {finite, tr≈1,
  hermitian}; the regression asserts {td-match vs exp baseline}. Cleanup possibility:
  fold smoke into regression (the regression already runs both paths and checks both
  return `ThermalizeResults` implicitly via `r_exp.trace_distances`). One run of
  `run_thermalize` would cover both. **Mark: weak duplicate at the workflow level.** Same
  pattern on TrotterDomain (lines 55 and 77).
- Line 87 ("Thermalization quality") at T_mix=5.0 is **distinct** — actually exercises
  thermalization, not just per-step agreement. Not duplicate.
- Line 108 ("Higher gqsp_degree") is the only test exercising `gqsp_degree=2`, i.e. the
  only test that hits the Clenshaw branch in `_gqsp_apply_polynomial`. Definitely not
  redundant.
- Line 122 ("Trajectory regression") is the only `run_trajectory` test for GQSP; not
  duplicate with `run_thermalize` tests.

**Summary of test redundancies:**
- Soft duplicate: smoke + regression at the same `(domain, params)` pair (lines 44/66
  on TimeDomain, 55/77 on TrotterDomain). Each pair could collapse into a single,
  slightly richer `@testset`.
- Soft duplicate: slope-(d+1) and tail-bound scaling (test_gqsp_polynomial.jl:126 and
  166) test the same Bessel tail in two parametrisations.
- All other testsets are useful sweeps or distinct properties.

## 5. Smells / surprises

### 5.1 `Hermitian` wrapper around `B` is inconsistent before GQSP

In `_precompute_coherent_unitary` (the Thermalize-path coherent unitaries), the GQSP
branches at `coherent.jl:108-110` and `coherent.jl:122-124` pass `B` **raw** to
`_gqsp_apply_polynomial`. The matrix-exp branches `coherent.jl:112` and `coherent.jl:126`
wrap `B` in `Hermitian(B)` to enforce hermiticity (silently symmetrises `B` then takes
the matrix exp). In the Trajectory path (`trajectories.jl:99-110`), `hermitianize!(B_a)`
is called explicitly **before** the `if config.with_gqsp` branch, so both arms see a
hermitised `B_a`.

**Asymmetry:** the Thermalize path's GQSP branch operates on a possibly-non-Hermitian
`B` (as built by `B_time`/`B_trotter`). Numerical noise on the order of
quadrature-truncation error can introduce small anti-Hermitian components. The
docstring of `_gqsp_apply_polynomial` (`coherent.jl:283-284`) explicitly says "do not
Hermitianize the output," which is correct for the **output**, but says nothing about
the **input**. `B_time`'s docstring (`coherent.jl:143`) says "Has to be on a symmetric
time domain, otherwise it can't be Hermitian", relying on the symmetric-grid contract.

Recommendation for the cleanup phase (qf-fzj.2): make Thermalize-path GQSP branches
hermitianize-then-call (matching the Trajectory path), or document explicitly why the
matrix-exp side must hermitianize but the GQSP side need not. Currently this is an
**asymmetric pre-call invariant** that future maintainers will trip on.

### 5.2 Magic number cap `gqsp_degree ≤ 100`

`misc_tools.jl:180-182` caps `gqsp_degree` at 100. There is no comment explaining the
choice. The Clenshaw recurrence is `O(d)` matmuls and is numerically stable for moderate
`d`; the Bessel-tail bound shrinks geometrically. 100 is plausible but unmotivated. Flag
for physics-checker (qf-fzj.6): is `d=100` actually meaningful given that the
matrix-exp baseline has zero polynomial-truncation error? The default `gqsp_degree=1`
produces `O((δα)²)` error which already matches the `O(δ²)` Lie-Trotter splitting error,
so for the simulator the only relevant degrees are 1 and 2.

### 5.3 `_gqsp_block_encoding_alpha` ℓ¹ recomputation per jump

See §3.5. The two `sum(abs, values(b_minus))` / `sum(abs, values(b_plus))` calls
(`coherent.jl:264-265`) are recomputed once per jump per simulation. Tiny but
suboptimal; flag only if cleanup rewrites the helper signature.

### 5.4 Two `_get_truncated_indices` overloads with identical bodies

`coherent.jl:385-393` defines two overloads for `Real` and `Complex` `fvals` whose
bodies are byte-identical (`return findall(abs.(fvals) .>= atol)`). The `Complex`
overload is the only one actually exercised (the truncated-func returns
`Vector{ComplexF64}`). The `Real` overload looks like a defensive leftover from a
previous refactor. **Not strictly GQSP-specific** (used by `_compute_truncated_func`
which feeds both the matrix-exp path and the GQSP path), but flagged because the GQSP
α-helper depends on the same machinery. Cleanup candidate during qf-fzj.2: drop the
`Real` overload, or unify under a single `AbstractVector{<:Number}` signature.

### 5.5 `gqsp_degree=0` is rejected; identity case is unreachable

The `@assert d ≥ 1` in `_gqsp_apply_polynomial` (`coherent.jl:294`) and the
`gqsp_degree < 1` rejection in `validate_config!` (`misc_tools.jl:177`) together mean
the `d=0` case (`f_0 = J_0(δα) I`, no operator structure) is unreachable. This is
fine — `d=0` would give a useless coherent step — but it's worth noting as a
deliberate design choice with no comment.

### 5.6 BS vs MW convention — no contradiction in the Julia code

Per MEMORY.md note `gqsp_bs_vs_mw_convention.md`, the thesis figures use the BS
(transpose) form with `L_0` placed last. The Julia simulator code in `src/coherent.jl`
contains **no `L_0` reference** at all — it computes only the post-selected
anc=|0⟩ block via the closed-form Chebyshev expansion. The BS-vs-MW distinction is
only relevant in `src/python/gqsp/circuit.py` (which does the angle-finding POC) and
the thesis circuits. **Cleared:** no Julia simulator code contradicts the convention
note. Do **not** "fix" anything here.

The docstring of `_gqsp_apply_polynomial` (`coherent.jl:286`) cites the POCs by path:
`scripts/scratch_gqsp_B_n3.jl, scripts/scratch_gqsp_random_h.jl`. The latter is the
1-ancilla qubitization referenced in `test/test_gqsp_polynomial.jl:33`. Cross-references
are correct.

### 5.7 Docstring claims `Trajectory`-path but signature is `Config{Thermalize}`

`_precompute_coherent_unitary` is documented (`coherent.jl:73-87`) as "Precompute
per-jump coherent unitaries for Kraus thermalization" with signature
`config::Config{Thermalize}`. It is only called from `furnace.jl:190-191`
(`run_thermalize`). The Trajectory path does NOT use this helper — it has its own
inlined logic in `trajectories.jl:91-114`. Two parallel implementations of the same
3-step recipe (build B, dispatch on gqsp, store in vector). See §3.1 for the
extract-helper recommendation.

### 5.8 `_precompute_coherent_unitary` does not handle DLL+GQSP

`validate_config!` accepts `with_gqsp=true` for **both** `KMS` and `DLL` constructions
(`misc_tools.jl:170-176`: only requires `with_coherent(construction)`, which is true for
both). However, `_precompute_coherent_unitary` (`coherent.jl:88-140`) only branches
on the **domain**, not the **construction**. For a `Config{Thermalize, TimeDomain, DLL}`
with `with_gqsp=true`, the function would enter the TimeDomain branch
(`coherent.jl:103-114`) and call `B_time` + `_gqsp_block_encoding_alpha`, both of which
need `b_minus`, `b_plus`, `gamma_norm_factor`. But `_precompute_data` for
`(TimeDomain, DLL)` returns only `(filter, time_labels, t0, oft_nufft_at_zero)`
(`furnace_utensils.jl:152-157`) — no `b_minus`/`b_plus`/`gamma_norm_factor`.

Result: the destructure `(; b_minus, b_plus, gamma_norm_factor) = precomputed_data` at
`coherent.jl:104` would error at runtime. The configuration is "validatable" but
unrunnable. This is a **gap** between `validate_config!` and the implementation.
Either:
  - tighten `validate_config!` to also reject `with_gqsp=true ∧ DLL` (the easier fix), or
  - extend `_precompute_coherent_unitary` to dispatch on `Config{..., ..., DLL}` and
    use the DLL coherent operator `G` + a DLL-specific block-encoding norm.

This is **out of scope for cleanup-only** but should be flagged for qf-fzj.6 (physics-
check) and qf-fzj.7 (integration). MEMORY entry on `feedback_dll_paper_first.md`
applies: any DLL+GQSP integration must consult the paper before coding.

### 5.9 `delta_eff` differs between the two paths

- `_precompute_coherent_unitary` (Thermalize path): `delta = delta_scale * config.delta`
  where `delta_scale = (1.0 / p_jump)` from `furnace.jl:191` when
  `rescale_by_inv_prob=true`, else `1.0`. This is then passed straight into
  `_gqsp_apply_polynomial` (`coherent.jl:110, 124`).
- `_build_trajectory_workspace` (Trajectory path): `delta_eff = delta / p_jump`
  (`trajectories.jl:96`), unconditionally. No `delta_scale` knob.

The Trajectory path is **always** rescaled by `1/p_jump`. The Thermalize path is
**conditionally** rescaled. This is a real semantic difference between the two paths,
not a bug, but it's another reason §3.1's helper-extraction needs to **not** flatten
the two paths blindly: the unitary's effective δ depends on the simulation type. The
extracted helper should take `delta_eff` as an argument; the call sites compute it.

---

## Cross-references to MEMORY.md

- `gqsp_bs_vs_mw_convention.md`: confirmed — Julia simulator is convention-agnostic
  (no `L_0`); aligned with the note. (§5.6)
- `feedback_fast_path_default.md`: GQSP path is opt-in via `with_gqsp=true`; the
  matrix-exp path is the default. The `_gqsp_apply_polynomial` `d=1` fast path
  (`coherent.jl:301-312`, no matmul) is itself the default GQSP branch and aligns with
  the fast-path-default principle. (no contradiction)
- `feedback_dll_paper_first.md`: §5.8 (DLL+GQSP gap) requires consulting the paper
  before any cross-coupling fix.
- `feedback_n3_minimum_test_size.md`: all GQSP integration tests use
  `num_qubits=3` (`test_gqsp_thermalize.jl:29`, `test_gqsp_polynomial.jl:240`). Aligned.
- `feedback_beta_test_values.md`: `test_gqsp_thermalize.jl` uses the global `BETA`
  (=10.0 from `test_helpers.jl:74`) for all tests. **Fully aligned with β=10 rule.**
  No β sweep — but for a polynomial-truncation regression vs matrix-exp, the β
  dependence is not load-bearing (the comparison is at fixed B). No flag.
- `feedback_targeted_tests_during_integration.md`: applies to any cleanup PR — run
  only `test/test_gqsp_*.jl` between commits; full suite at end of phase.

---

## Summary

- GQSP surface is small: 0 exported names, 2 `Config` fields, 2 private helpers
  (`_gqsp_block_encoding_alpha`, `_gqsp_apply_polynomial`), 3 call sites each across
  `coherent.jl` (×2) and `trajectories.jl` (×1).
- Main cleanup target: the **3× near-clone `if config.with_gqsp ... else exp(...)`
  block** at `coherent.jl:108-113`, `coherent.jl:122-127`, and `trajectories.jl:101-110`
  — extract one helper that takes `B`, `delta_eff`, and the precomputed data.
- Smells worth flagging for the physics/integration phases: (a) Thermalize-path GQSP
  branch does **not** hermitianize `B` before calling `_gqsp_apply_polynomial` while
  the matrix-exp branch wraps in `Hermitian` and the Trajectory path calls
  `hermitianize!` (`coherent.jl:106-110` vs `trajectories.jl:100`); (b) `validate_config!`
  accepts `with_gqsp=true ∧ DLL` but `_precompute_coherent_unitary` would runtime-error
  on it (`misc_tools.jl:170-183` vs `coherent.jl:103-138`); (c) magic cap
  `gqsp_degree ≤ 100` with no comment.
