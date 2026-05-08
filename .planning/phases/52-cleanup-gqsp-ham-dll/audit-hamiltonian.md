# Audit: Hamiltonian generation (qf-fzj.1)

Read-only audit of the Hamiltonian-generation surface introduced/refactored in epic
`qf-k1u.5` (2D Heisenberg builders + multi-term disordering). Scope:
`src/hamiltonian.jl`, related struct fields in `src/structs.jl`, the
loader path in `src/misc_tools.jl`, the Trotter consumer in
`src/trotter_domain.jl`, the cached BSON families in `hamiltonians/`, and
`test/test_hamiltonian.jl`.

---

## 1. Public API surface

Names exported from `src/QuantumFurnace.jl` (lines 68, 73-74) related to
`HamHam` / Hamiltonian construction:

| Export | Defined at | User-callable? |
|---|---|---|
| `HamHam` (type) | `src/hamiltonian.jl:26` | yes (struct, plus 4 outer constructors) |
| `find_ideal_heisenberg` | `src/hamiltonian.jl:279` | yes (1D builder) |
| `find_ideal_2d_heisenberg` | `src/hamiltonian.jl:330` | yes (2D builder) |
| `load_hamiltonian` | `src/misc_tools.jl:13` | yes (BSON loader, **legacy filename only**) |
| `create_bohr_dict` | (defined elsewhere; called by HamHam ctors at lines 88, 145, 193) | yes |
| `pad_term` | `src/misc_tools.jl:376` | yes (helper used by both base and disordering paths) |

`HamHam` outer constructors (all callable from outside the module since `HamHam` is exported):
- `HamHam(terms, coeffs, num_qubits, beta; ...)` — no disorder, `src/hamiltonian.jl:60`
- `HamHam(terms, coeffs, dis_terms, dis_coeffs, num_qubits, beta; ...)` — multi-term disorder, `src/hamiltonian.jl:109`
- `HamHam(terms, coeffs, dis_term, dis_coeffs, num_qubits, beta; ...)` — single-term convenience, `src/hamiltonian.jl:167`
- `HamHam(raw::NamedTuple, beta)` — from `find_ideal_*` output, `src/hamiltonian.jl:189`

Private (not exported, but reached via `QuantumFurnace.<name>` in tests):
`_pad_two_site_op`, `_construct_2d_heisenberg_base`, `_optimize_disordered_heisenberg`,
`_construct_disordering_terms` (2 overloads), `_construct_base_ham`,
`_rescaling_and_shift_factors`, `_unpack_disordering_fields`, `_gibbs_in_eigen`,
`_load_hamiltonian_bson`.

---

## 2. Internal helpers — single-call helpers flagged for inlining

Format: `name | file:line | signature | callers`. **★** marks single-call
helpers that are inlining candidates.

### `_pad_two_site_op` — `src/hamiltonian.jl:446`
```
_pad_two_site_op(term::Vector{Matrix{ComplexF64}}, num_qubits::Int, q1::Int, q2::Int)
```
Callers:
- `src/hamiltonian.jl:507, 510, 516, 519` — all four bond-direction branches
  in `_construct_2d_heisenberg_base` (one helper, one consumer; 4 call sites).
- `test/test_hamiltonian.jl:19, 25, 29, 33, 34, 39–42` — direct unit tests.

Has its own argument validation block (lines 447–455) and a small `q1<q2`
canonicalisation (lines 458–464). Worth keeping as a named helper because
the q1/q2 ordering logic is non-trivial and the 2D builder calls it from 4
sites.

### `_construct_2d_heisenberg_base` — `src/hamiltonian.jl:491` ★
```
_construct_2d_heisenberg_base(Lx, Ly, terms, coeffs; periodic_x=true, periodic_y=true)
    -> Hermitian{ComplexF64, Matrix{ComplexF64}}
```
Callers:
- `src/hamiltonian.jl:342` — sole production caller, inside `find_ideal_2d_heisenberg`.
- `test/test_hamiltonian.jl:53, 58, 66, 108, 117, 119, 124, 126` — direct tests.

**Single production call site.** Could be inlined into `find_ideal_2d_heisenberg`,
but the unit-test surface is large (8 references) and the helper is well-documented.
Inlining would force tests to construct via the public `find_ideal_2d_heisenberg`
(which runs the optimisation loop for `batch_size` random realisations) — ~200×
slower for the bond-counting / Hermiticity assertions. **Keep as helper.**

### `_optimize_disordered_heisenberg` — `src/hamiltonian.jl:353`
```
_optimize_disordered_heisenberg(base_hamiltonian, terms, coeffs, num_qubits,
    disordering_terms; batch_size, periodic, disorder_strength)
```
Callers:
- `src/hamiltonian.jl:287` — from `find_ideal_heisenberg`
- `src/hamiltonian.jl:345` — from `find_ideal_2d_heisenberg`

Two production callers — exactly what was promised by the qf-k1u.5 refactor.
**Keep.**

### `_construct_disordering_terms` (single-term overload) — `src/hamiltonian.jl:528` ★
```
_construct_disordering_terms(term::Vector{Matrix{ComplexF64}},
    coeffs::Vector{Float64}, num_qubits::Int64) -> Hermitian{...}
```
Callers: **none in src or test**. The grep for `_construct_disordering_terms`
returns only call sites that pass `Vector{Vector{...}}` (multi-term):
`src/hamiltonian.jl:128` and `src/hamiltonian.jl:373`. The single-term overload
exists but is **dead code**. See §3 for details.

### `_construct_disordering_terms` (multi-term overload) — `src/hamiltonian.jl:543`
```
_construct_disordering_terms(terms::Vector{Vector{Matrix{ComplexF64}}},
    coeffs::Vector{Vector{Float64}}, num_qubits::Int64) -> Hermitian{...}
```
Callers (production):
- `src/hamiltonian.jl:128` — multi-term `HamHam` constructor.
- `src/hamiltonian.jl:373` — `_optimize_disordered_heisenberg` inner loop.

### `_construct_base_ham` — `src/hamiltonian.jl:416`
Callers: `src/hamiltonian.jl:73, 127, 285` (no-disorder ctor, multi-term ctor,
`find_ideal_heisenberg`); also `test/test_hamiltonian.jl:50` for the 1D-vs-2D
equivalence test. Three production callers — keep.

### `_rescaling_and_shift_factors` — `src/hamiltonian.jl:559`
Callers: `src/hamiltonian.jl:75, 131, 376` (all three Hamiltonian-construction
paths). Keep.

### `_unpack_disordering_fields` — `src/hamiltonian.jl:224` ★
```
_unpack_disordering_fields(raw::NamedTuple, ::Type{T}) -> (terms, coeffs)
```
Callers: `src/hamiltonian.jl:197` — sole call in the NamedTuple constructor.
**Single call site, ~20 lines.** Could be inlined into the NamedTuple
constructor; the docstring at `:217` notes its role as the legacy-format
adapter. Inlining is reasonable but the named helper makes the back-compat
shim discoverable. **Borderline — keep unless cleaning the legacy path.**

### `_gibbs_in_eigen` — `src/hamiltonian.jl:49`
Callers: `src/hamiltonian.jl:89, 146, 193` (all three HamHam ctor paths). Keep.

### `pad_term` — `src/misc_tools.jl:376`
Callers: 14+ sites in `src/`, `test/`, scripts (full grep above). Used by
both base-term construction and disordering-term construction; extensively
used outside the Hamiltonian-generation module. Keep.

---

## 3. Dedup / simplification candidates

### 3a. 1D vs 2D path sharing through `_optimize_disordered_heisenberg`

The kernel refactor (qf-k1u.5) is clean — both paths share:
- `_construct_disordering_terms` (multi-term overload) for the disorder Hamiltonian
- `_rescaling_and_shift_factors` for the spectrum-rescale step
- the `eigen` + `nu_min = minimum(diff(eigvals))` body

What differs:
- `find_ideal_heisenberg:285` calls `_construct_base_ham(terms, coeffs, num_qubits; periodic=periodic)`
- `find_ideal_2d_heisenberg:342` calls `_construct_2d_heisenberg_base(Lx, Ly, terms, coeffs; periodic_x, periodic_y)`

That asymmetry is unavoidable — the bond geometry differs. **No further dedup possible.**
The `terms = [[X,X], [Y,Y], [Z,Z]]` literal is duplicated at lines 284 and 341 —
trivial; could be a `const` but doesn't matter.

### 3b. The 4 `HamHam` constructors

Layout:
- Ctor (1) `HamHam(terms, coeffs, num_qubits, beta; ...)` lines 60–107 (47 lines)
- Ctor (2) `HamHam(terms, coeffs, dis_terms, dis_coeffs, num_qubits, beta; ...)` lines 109–164 (55 lines)
- Ctor (2b) `HamHam(terms, coeffs, dis_term, dis_coeffs, num_qubits, beta; ...)` lines 167–175 — pure
  delegation to (2). 9 lines, **already a thin wrapper**. ✓
- Ctor (3) `HamHam(raw::NamedTuple, beta)` lines 189–215 (27 lines)

Boilerplate shared between (1) and (2):
- Identical mixed-precision check (lines 67–71 ≡ 117–121). 5 lines duplicated.
- Identical eigen + bohr_freqs + bohr_dict + gibbs setup (lines 79–105 ≡ 135–163,
  modulo the disordering-coeff handling). ~25 lines duplicated.

(1) is *almost* "(2) with `disordering_terms = Matrix{ComplexF64}[]`". Compare:
- (1) calls `_construct_base_ham(terms, coeffs, num_qubits; periodic=periodic)` (line 73)
- (2) calls `_construct_base_ham(terms, coeffs, num_qubits)` (line 127, **no `periodic` kwarg passed!** — see §5)
  then `_construct_disordering_terms(...)` then sums

**Simplification candidate:** rewrite (1) as a thin wrapper around (2) with empty
`disordering_terms`/`disordering_coeffs`. The multi-term overload at
`:543` already iterates `for (term, term_coeffs) in zip(terms, coeffs)`; an
empty zip yields a zero matrix, so the math works. Stored
`disordering_terms` would become `Vector{...}[]` instead of `nothing`,
which **would change the BSON wire format** if anyone serialised a (1)-built
HamHam — but neither cached BSON family was generated via (1) (F1/F2/F3 all
go through `find_ideal_heisenberg` / `find_ideal_2d_heisenberg`), and
`_load_hamiltonian_bson` reads field 6 as `nothing`-or-vector and adapts.
**Net: safe to consolidate; saves ~30 LOC.**

### 3c. `HamHam(raw::NamedTuple, beta)` and `_unpack_disordering_fields`

`_unpack_disordering_fields` (`:224`) handles three cases:
1. New format: `raw.disordering_terms === nothing` → `(nothing, nothing)`
2. New format: `raw.disordering_terms` is a vector → unwrap each entry
3. Legacy format: `raw.disordering_term` (singular, with optional nothing) → wrap into a 1-element vector

The cached BSON files written by `generate_hamiltonians.jl` always go through
the new format (the `find_ideal_*` returns from `_optimize_disordered_heisenberg`
include `disordering_terms = disordering_terms` at line 405). The legacy
single-key NamedTuple is reachable from:
- `src/misc_tools.jl:_load_hamiltonian_bson` lines 59–71 — synthesises a NamedTuple
  with `disordering_term =` (singular) when reading **legacy** `heis_disordered_periodic_n*.bson` BSONs.
- `test/test_helpers.jl:_load_test_hamiltonian` lines 52–64 — duplicate of the above.

So the legacy branch in `_unpack_disordering_fields` is load-bearing for
`heis_disordered_periodic_n3..5.bson` (dated 2025-10-28) only. **Cannot remove**
without rewriting `_load_hamiltonian_bson` + `_load_test_hamiltonian` to map to
the new key, or migrating those three legacy BSON files.

**Minimal alternative:** at the loader, write `disordering_terms = [dt]` instead
of `disordering_term = dt`, eliminating the legacy branch from `_unpack_disordering_fields`.
That would also **remove the duplication between `_load_hamiltonian_bson` and `_load_test_hamiltonian`**
(see §5).

### 3d. `_construct_disordering_terms` single-term vs multi-term

The single-term overload at `:528` is called nowhere (full repo grep). It
implements
```julia
disordering_hamiltonian += coeffs[q] * pad_term(term, num_qubits, q)
```
which is exactly what the multi-term overload at `:543` does in its inner
`for q in 1:num_qubits` after fixing one `(term, term_coeffs)` pair. **The
single-term overload is dead code.** Remove it; if a single-term call site
ever returns, wrap as `[term]`/`[coeffs]` and call the multi-term overload.

### 3e. `_load_hamiltonian_bson` and `_load_test_hamiltonian` are near-duplicates

`src/misc_tools.jl:27–74` (48 lines) ≡ `test/test_helpers.jl:21–67` (47 lines).
Verbatim duplication except for the `init = @__MODULE__` vs `init = QuantumFurnace`
on lines 41/34, and the docstring. **Test should call `QuantumFurnace._load_hamiltonian_bson`**
instead of maintaining its own copy — same legacy back-compat coverage,
half the maintenance.

---

## 4. Test redundancies in `test/test_hamiltonian.jl`

8 testsets, 188 lines total. None of them sweep n=3,4,5,6 for bipartite-collision
verification — that work lives elsewhere (in scripts under `scripts/`, see
MEMORY.md "Bohr Frequency Collision Root Cause"). The test file is focused
on the qf-k1u.5 builders.

| # | Testset | At | Asserts | Parameters | Redundant? |
|---|---|---|---|---|---|
| 1 | `_pad_two_site_op: adjacent and non-adjacent placements` | `:17` | size, hermiticity, exact-kron equivalence, q1↔q2 symmetry | n∈{2,4} | n=2 violates `feedback_n3_minimum_test_size`; n=4 covers all cases. **Drop n=2 placement; n=4 suffices** |
| 2 | `_pad_two_site_op: argument validation` | `:38` | 4 ArgumentError throws | n=4 | Keep — covers each guard once |
| 3 | `_construct_2d_heisenberg_base: 2x1 and 1x2 lattices match 1D periodic n=2` | `:46` | 2D 2×1 ≡ 2D 1×2 ≡ 1D n=2 | (2,1)+(1,2) | n=2 violates the n=3 minimum rule. The 1D-n=2 reference itself uses the periodic-double-counting convention and is a degenerate case. **Re-base on (3,1)/(1,3) ≡ 1D n=3** |
| 4 | `_construct_2d_heisenberg_base: dimension and Hermiticity for several lattices` | `:63` | dim, hermiticity, traceless | (2,2),(2,3),(3,3),(2,5) | (2,2) is borderline (mixed n=2 directions). **Drop (2,2); keep (2,3),(3,3),(2,5)** |
| 5 | `_construct_2d_heisenberg_base: bond counting via Frobenius norm` | `:76` | exact closed-form `‖H‖_F²` | (2,2),(2,3),(3,3),(2,5) | Same observation. (2,2) duplicates the wrap-double-counting check that (2,3) already exercises (Lx=2 path). **Drop (2,2)** |
| 6 | `_construct_2d_heisenberg_base: open boundary disables wrap` | `:115` | OBC norm < PBC norm; explicit 3×3 OBC bond count | (2,2),(3,3) | (2,2) is the wrap-double-counting case from #4/#5; (3,3) is generic OBC. **Drop (2,2) check (line 117–121); keep 3×3 closed-form** |
| 7 | `find_ideal_2d_heisenberg: returns a valid raw NamedTuple` | `:133` | nu_min>0, sizes, periodic flag, dis_coeff shape, spectrum bounds, hermiticity | Lx=Ly=2 | n=4 OK; **representative end-to-end check, keep** |
| 8 | `find_ideal_2d_heisenberg: HamHam wrap end-to-end` | `:153` | HamHam type, gibbs trace=1 | Lx=Ly=2 | Overlaps #7 on data construction (uses `find_ideal_2d_heisenberg`). The novel assertion is the `HamHam(raw, β)` wrap. **Could merge into #7** |
| 9 | `find_ideal_2d_heisenberg: argument validation` | `:162` | 2 ArgumentErrors for L=0 | — | Keep |
| 10 | `find_ideal_heisenberg: still works after refactor` | `:168` | nu_min>0, dim=8, dis_terms=2 | n=3, [[Z],[Z,Z]] | Keep — pins the multi-term 1D path |
| 11 | `find_ideal_heisenberg: disorder_strength scales` | `:178` | per-coeff bound | n=3 | Keep |

**Coverage of the four `HamHam` constructors:** none of the constructors (1),
(2), (2b) are tested directly — they are reached only transitively via
`find_ideal_2d_heisenberg(...)` → returned NamedTuple → constructor (3) at
line 155 (`HamHam(raw, 1.0)`). Constructors (1), (2), (2b) have **no direct
test coverage at all** in `test/test_hamiltonian.jl`. The legacy
`_load_test_hamiltonian` path in `test_helpers.jl` exercises ctor (3) on
the legacy NamedTuple shape via the n=3,4 fixtures, indirectly.

**1D vs 2D builder overlap:** none — the 1D tests only assert `find_ideal_heisenberg`
top-level behaviour, not bond placement (that's tested in
`_construct_2d_heisenberg_base: 2x1 and 1x2 lattices match 1D periodic n=2`).

**Net redundancy:** mostly the (2,2)/n=2 cases (4 instances at lines 19, 53–60,
64, 117). Removing them would tighten the file by ~20 lines.

---

## 5. Smells / surprises

### 5a. Dead code: single-term `_construct_disordering_terms` overload (`src/hamiltonian.jl:528`)
Defined but never called. See §3d.

### 5b. Constructor (2) loses `periodic` keyword on its way to `_construct_base_ham`
Line 127: `base_hamiltonian = _construct_base_ham(terms, coeffs, num_qubits)` —
the function signature accepts `periodic::Bool = true` (line 417), so the
default is taken. But constructor (2) takes `periodic::Bool = true` as its
own kwarg (line 112) and **silently drops it**. Constructor (1) at line 73
correctly forwards it (`periodic=periodic`). For `periodic=false`, ctor (2)
would build a periodic base Hamiltonian — wrong.

This is reachable from any user code that calls
`HamHam(terms, coeffs, dis_terms, dis_coeffs, n, β; periodic=false)`; no
production caller does this today (`generate_hamiltonians.jl` always uses
periodic=true and goes through `find_ideal_*` anyway), but the ctor's
signature suggests it works. **Bug.**

### 5c. `load_hamiltonian` is broken for the new BSON families
`load_hamiltonian("heis_xxx", n; ...)` constructs the path
`heis_xxx_disordered_periodic_n<n>.bson` (line 16), which **does not exist**.
The only matching files in `hamiltonians/` are
`heis_xxx_zzdisordered_periodic_n*.bson` and `heis_xxx_clean_periodic_n*.bson`.
And `heis_xxz_2d_*` doesn't fit the format at all.

Effect:
- Scripts that use `load_hamiltonian("heis", n; ...)` only work for the legacy
  `heis_disordered_periodic_n3..5.bson` files (still on disk). All scratch
  scripts (`scratch_dll_dissipator_nufft.jl`, `scratch_benchmark_dll_vs_ckg.jl`,
  `scratch_dll_coherent_*.jl`, `scratch_dll_floor_*.jl`) hit this path.
- `scripts/scratch_smoke_dll_coherent.jl:7` calls `load_hamiltonian("heis_xxx", 2; ...)`
  → tries to read `heis_xxx_disordered_periodic_n2.bson`, which does not exist.
  This script must be currently broken.

`load_hamiltonian` was designed pre-qf-k1u.5 and never updated. It also still
documents itself as loading "legacy HamHam structs" — which would not match
the new NamedTuple-format BSONs even if the filename pattern were fixed.

### 5d. `_load_hamiltonian_bson` only handles the legacy struct serialisation
Lines 27–74. The function unpacks fields[1..14] of a serialised
`HamHam{Float64}`. Both family generators (`generate_hamiltonians.jl:67`,
`scripts/...`) write **NamedTuples**, not HamHam structs (the helper
`_save_and_report` does `BSON.bson(path, Dict(:hamiltonian => hamiltonian))`
where `hamiltonian` is the raw NamedTuple). Only the three pre-qf-k1u.5
files (`heis_disordered_periodic_n3..5.bson`, dated 2025-10-28) are in the
serialised-struct format. **The loader function is now misnamed and inflexible.**

### 5e. Duplicate BSON loader between `src/misc_tools.jl` and `test/test_helpers.jl`
See §3e. ~50 lines copy-pasted.

### 5f. Inconsistent argument validation between 1D and 2D builders
- `find_ideal_heisenberg` (`:279`): no validation of `num_qubits` (any
  `num_qubits ≥ 1` proceeds; `num_qubits = 0` would silently produce
  `2^0 × 2^0 = 1 × 1` matrices and probably crash inside `eigen`).
- `find_ideal_2d_heisenberg` (`:330`): explicitly throws `ArgumentError`
  for `Lx < 1` or `Ly < 1` at line 336.

The 2D path also asserts `length(terms) == length(coeffs)` inside
`_construct_2d_heisenberg_base:495`; the 1D path asserts the same inside
`_construct_base_ham:419`. So this is consistent.

### 5g. `_optimize_disordered_heisenberg` discards `terms` argument in returned NamedTuple
Lines 401–413: returns `terms = terms` literally. This means the **2D path
returns the bare `terms = [[X,X], [Y,Y], [Z,Z]]` list, not the lattice
geometry** — a downstream caller using only `raw.terms` cannot reconstruct
the 2D bond placement. This is fine for the rest of the pipeline because
`raw.matrix` carries the full operator, but if a future Trotterizer ever
relied on `raw.terms` to know which bond pairs are present, it would
silently treat the 2D system as a 1D chain. Currently no Trotterizer or
HamHam constructor uses `raw.terms` for bond-geometry inference (they read
`hamiltonian.data` directly), so this is **latent** and might surprise a
new contributor.

### 5h. `_construct_2d_heisenberg_base` periodic-wrap double-counting
Lines 506–512: when `Lx == 2` with periodic_x, both `i=1 → i=2` and the
wrap `Lx → 1` add the same bond (sites 1↔2), giving coefficient 2× per bond.
The docstring (lines 487–489) calls this out and says it matches the 1D
n=2 convention. **Documented but counter-intuitive** — easy to misread when
auditing physics. The bond-counting test at line 84 explicitly handles this
case in `expected_sumsq`. Not a bug, but worth flagging as physics surprise.

### 5i. `find_ideal_2d_heisenberg` site-index ordering: `(i-1)*Ly + (j-1) + 1`
Line 500. So adjacent **j**-values map to consecutive qubits and adjacent
**i**-values are separated by `Ly`. This is consistent across
`_construct_2d_heisenberg_base` and the test at `verify_hamiltonians.jl:70`.
The docstring (lines 310–312) covers it. Note: this is the *opposite* of
many 2D-physics conventions (row-major in `Lx`). Easy to flip when
reading code from another reference. **Documented but flag-worthy.**

### 5j. `HamHam` field `periodic` is a single Bool for both 1D and 2D
For 2D systems, `periodic = periodic_x && periodic_y` (line 346). A 2D
system with mixed BCs (periodic_x=true, periodic_y=false) is possible at
construction but **can't be round-tripped through `HamHam` losslessly**:
`find_ideal_2d_heisenberg(2, 3, [...]; periodic_x=true, periodic_y=false)` →
`raw.periodic = false` → indistinguishable from fully open or fully closed.
Latent — `generate_hamiltonians.jl` always uses periodic in both directions.

### 5k. Mixed-precision check is downward-only
Constructors (1) and (2) at lines 67–71 and 117–121 reject `T ∈ {Float16, Float32}`
when terms are `ComplexF64`, but quietly accept `T = Float64` and any
`T <: AbstractFloat` other than F16/F32 (e.g. `BigFloat`, `Float128`). The
error message says "downward mismatch errors, upward promotion allowed".
The constructor body for (1) at line 76 then forces the rescaled
Hamiltonian to be `Hermitian{ComplexF64, Matrix{ComplexF64}}` regardless of
`T`, then casts at the final `HamHam{T}(...)` (line 91). For `T=BigFloat`,
this means the eigen decomposition runs at Float64 precision and gets
truncated to `Complex{BigFloat}` storage — likely not what the user
intended. **Latent quality smell.**

### 5l. Unused/Any-typed fields in `HamHam`
None visible. All 15 fields are used. `bohr_dict::Dict{T, Vector{CartesianIndex{2}}}`
uses Float keys, which has the standard "near-equal floats hash apart" trap;
`create_bohr_dict` is presumably bucketed (defined elsewhere).

### 5m. BSON forward-compat: struct field renames are unsafe
`HamHam` is a positional struct — `BSON.parse` reads fields in order
`fields[1..14]` (legacy 14-field shape, see `_load_hamiltonian_bson:36–39`).
The current struct has 15 fields. Adding a field at the **end** is safe for
load-if-not-present; renaming any of `disordering_terms`, `disordering_coeffs`,
or any field referenced by index in `_load_hamiltonian_bson` would silently
read wrong data. The current loader hard-codes the legacy index 6 →
`disordering_term` (singular), index 7 → `disordering_coeffs`. Reordering
or merging fields requires updating both loaders. **Documented constraint
from the task description, recorded here for completeness.**

### 5n. Contradicts MEMORY.md? No.
The audit findings are consistent with MEMORY.md notes:
- "`HamHam` struct: `disordering_terms::Union{Vector{Vector{Matrix{Complex{T}}}}, Nothing}` — multiple terms" ✓
- "Constructors: (1) no disorder, (2) multi-term, (2b) single-term convenience, (3) from NamedTuple" ✓
- "NamedTuple constructor handles both legacy `disordering_term` and new `disordering_terms` keys" ✓
- "Trotter: even-n uses 2 bond groups (clean), odd-n needs 3rd group for wrapping bond" ✓
  (verified at `src/trotter_domain.jl:71–74`, the `is_bdr_strange` branch)
- "Z-only disorder hides bipartite collisions for even n" — see `generate_hamiltonians.jl:9`
  comment which warns the legacy `heis_disordered_periodic_n*.bson` exhibits this.

---

## Summary

- **Dead code:** single-term `_construct_disordering_terms` overload (`src/hamiltonian.jl:528`).
- **Bug:** ctor (2) drops `periodic` kwarg before calling `_construct_base_ham` (`:127`).
- **Broken loader:** `load_hamiltonian` filename pattern is stuck on the legacy `heis_*disordered_*.bson` shape; new F1/F2/F3 families and `scratch_smoke_dll_coherent.jl` are unreachable through it.
