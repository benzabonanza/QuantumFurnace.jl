# Phase 45: Rename `b` to `s` and add eta-regularized smooth Metropolis variant

## Goal

After this phase, the codebase consistently calls the smooth-Metropolis lower-limit
shift `s` (matching the thesis notation in `2_methods.tex` eq. 524 and
`metro-funcs.md`), and supports the previously-forbidden case **`a = 0, s > 0`**:
the eta-regularized smooth-Metropolis `b_+` from the thesis main case.

## Research Summary

### Where `b` is used today (verified)

The `b::Union{T, Nothing}` field on `Config` (in `src/structs.jl:94`) is referenced from:

| File | Lines | Usage |
|------|-------|-------|
| `src/structs.jl` | 55, 67-71, 94 | docstring, field declaration |
| `src/bohr_domain.jl` | 36, 44, 62, 70, 82, 90, 96-101, 112, 125, 130-135 | `_pick_alpha_*`, `create_alpha`, `create_alpha_gns` (`u_min = sqrt(beta * sigma^2 * (1 + b) / 2)`) |
| `src/coherent.jl` | 223 | `_compute_b_plus_smooth(t, beta, sigma, a, b)` and inside body |
| `src/energy_domain.jl` | 13, 15, 16, 31, 34, 36, 56, 62, 66, 115, 121, 126 | branch logic and `u_min = sqrt(config.beta * config.sigma^2 * config.b / 2)` |
| `src/furnace_utensils.jl` | 150 | `_select_b_plus_calculator` passes `config.b` |
| `src/misc_tools.jl` | 81-82, 93-95, 151-158, 218-229, 248-259 | filenames, `_print_press`, `validate_config!` rule |
| `src/results.jl` | 40 | `_config_to_dict` writes `:b` key (loop in `_dict_to_config_kwargs:104` reads it) |
| `src/simulation_time.jl` | 179-187 | `_determine_filter_info` uses `b_val`, returns Dict `:b => b_val` |
| `test/test_helpers.jl` | 201 | `make_config` passes `b = 0.4` |
| `test/test_simulation_time.jl` | 22 | constructs SimulationTimeBudget with `:b => 0.4` |
| `scripts/*.jl` (7 files) | various | keyword `b = 0.4` |

### Conflicting `b`/`s` identifiers in `src/coherent.jl`

These are unrelated to the smooth-Metro shift parameter and must NOT be renamed:

- `_compute_b_minus`, `_compute_b_plus`, `_compute_b_plus_metro`, `_compute_b_plus_smooth` — function names referring to the coherent operator B's `b_minus`/`b_plus` Fourier kernels
- `b_minus`, `b_plus` (Dict variables in `B_time`/`B_trotter`) — same kernels
- `b_plus_summand`, `b_t`, `b_s`, `B_time`, `B_trotter` — local variables

The local loop variable `for s in keys(b_plus)` (`B_time` line 122; `B_trotter` line 168 destructures `(s, b_s)`) IS a time key, NOT the new `s` shift parameter, but renaming
to `tau` improves readability after the rename (since after the rename `s` would be
a Config field name and confusing).

### `_compute_b_plus_metro` t=0 limit

Current code hard-codes `1 / (2*sqrt(2)*pi^2)`. This is the L'Hopital limit of the
formula at `t=0` ONLY when `sigma*beta = 1` AND `s = 0`. For the new variant we
must use the more general analytic limit:

```
lim_{t->0} (e^{-σ²β²(2t² + it)(1+s)} + 𝟙(|t|≤η)·i(2t + i)) / (t(2t + i))
       = (2 - σ²β²(1+s)) / (2√2 π²)
```

Verification: setting `σβ = 1`, `s = 0` gives `(2 - 1)/(2√2 π²) = 1/(2√2 π²)` —
exactly the legacy hard-coded value. For all our tests (`σβ = 1`) and most scripts,
the legacy value still holds at `s = 0`. For `scripts/biexp_mixing_verify_trotter.jl`
which uses `σ = 0.8/β` (so `σβ = 0.8`), the new formula is more accurate.

### `validate_config!` rule needs to relax

Current rule (`misc_tools.jl:151-158`):

```julia
if config.with_linear_combination && config.a == 0.0
    if config.b != 0.0
        push!(errors, "For linear combinations with b != 0, a must also be non-zero.")
    end
    if config.domain isa Union{TimeDomain, TrotterDomain} && with_coherent(...) && (isnothing(config.eta) || config.eta <= 0.0)
        push!(errors, "...eta must be > 0.")
    end
end
```

The inner `b != 0` check explicitly forbids the new `a = 0, s > 0` case. After the
rename + new variant we keep the eta requirement (still need eta when `a = 0` for
ALL `s ≥ 0`, since the eta regularization is what cures the t=0 singularity in the
plain Chen kernel) but drop the `b != 0` rejection.

### Beads issue tracking

There is no existing beads issue for this phase. Will create a new epic
`Phase 45: rename b → s and add a=0,s>0 smooth Metropolis` with child tasks per
commit.

## Tasks

### Task 1: Rename Config field `b` → `s` and update docstrings (Task A.1)

**Files**: `src/structs.jl`

**Action**:
- Line 94: `b::Union{T, Nothing} = nothing` → `s::Union{T, Nothing} = nothing`
- Line 55 docstring: `- \`a\` and \`b\`: Parameters for the linear combination type.` → `- \`a\` and \`s\`: Parameters for the linear combination type.`
- Lines 67-71: replace the case table with the four-row table from Task B (kept here in this same commit since the table inherently mentions both rename and the new `(0, >0)` row):
  ```
  ## Currently possible linear combinations:
  (a, s) =
  - (0, 0)   - plain Metropolis (kinky, eta-regularized in time domain)
  - (0, >0)  - smooth Metropolis (eta-regularized, kink-smoothed by s; thesis-main case)
  - (>0, 0)  - a-regularized smooth Metro (alternative regularization, no kink-smoothing)
  - (>0, >0) - a-regularized Glauberish (smooth in both senses)
  ```

**Note**: Since the test suite imports `Config` and tests pass `b = 0.4`, after
this commit alone the tests will FAIL. To keep "tests pass after every commit"
true, this task must be combined with Task 2-7 below into ONE atomic commit.

**Decision**: Bundle Task 1 + 2 + 3 + 4 + 5 + 6 + 7 into a single commit
`refactor(config): rename b → s parameter to match thesis notation`. Subdividing
breaks the cross-file references and would either fail tests or require throwaway
shim code. Justification: a pure rename (no behavior change) is one logical change.

### Task 2: Update `bohr_domain.jl` (Task A.2 contents)

**Files**: `src/bohr_domain.jl`

**Action**:
- `_pick_alpha_kms` (line 76-88): rename local `b = config.b` → `s = config.s`, and update the closure call `create_alpha(nu_1, nu_2, beta, sigma, a, b)` → `(..., a, s)`
- 2-arg `_pick_alpha` for KMS (line 60-66): `config.a, config.b` → `config.a, config.s`
- 2-arg `_pick_alpha` for GNS (line 68-74): same
- `_pick_alpha_gns` (line 106-118): same as KMS variant
- `create_alpha` (line 90-104): rename param `b` → `s` everywhere; `u_min = sqrt(beta * sigma^2 * (1 + b) / 2)` → `... (1 + s) / 2)`
- `create_alpha_gns` (line 125-138): same
- `create_f` (line 44-47): `b::Real` → `s::Real`, pass `s` instead of `b`
- `_pick_f` (line 30-42): rename local `b = config.b` → `s = config.s`, update closure

### Task 3: Update `coherent.jl`

**Files**: `src/coherent.jl`

**Action**:
- `_compute_b_plus_smooth(t, beta, sigma, a, b)` (line 223-226) → `(t, beta, sigma, a, s)`. Body: `exp(-a*b/2)` → `exp(-a*s/2)`, `(1 + b)` → `(1 + s)`. Local variable `b_vals` (line 224) is fine — that's the function output, not the shift parameter, but rename for consistency: leave alone, it's not ambiguous.
- `B_time` (line 105-149): rename inner loop var `for s in keys(b_plus)` → `for tau in keys(b_plus)` and `t_s = s * beta` → `t_tau = tau * beta`. The `b_plus[s]` access site (line 135) becomes `b_plus[tau]`.
- `B_trotter` (line 152-198): rename `for (s, b_s) in b_plus` → `for (tau, b_tau) in b_plus`. Update `s * beta` → `tau * beta` (line 169) and `b_s .* ...` → `b_tau .* ...` (line 183).

**Decision (Task A.3)**: Function names `_compute_b_plus_*`, `b_minus`, `b_plus`,
`b_plus_summand`, `B_time`, `B_trotter`: NOT renamed. These refer to the Fourier
kernels of the coherent B operator and have nothing to do with the shift `s`.

### Task 4: Update `energy_domain.jl`

**Files**: `src/energy_domain.jl`

**Action**:
- Replace all `config.b` (12 occurrences) with `config.s`.
- Update branch comments: `# No time singularity but kinky Metro in energy` (line 56) — this is the `s == 0, a != 0` case (i.e. `b == 0` was misleading; rename the comment too: `# a-regularized, no smoothing (s=0)`). Likewise line 62 — `# a-regularized + smoothed (s>0)`.
- The arithmetic `u_min = sqrt(config.beta * config.sigma^2 * config.b / 2)` (lines 16, 36, 66, 126) becomes `u_min = sqrt(config.beta * config.sigma^2 * config.s / 2)`.

### Task 5: Update `furnace_utensils.jl`

**Files**: `src/furnace_utensils.jl`

**Action**: line 150 — `(_compute_b_plus_smooth, (config.beta, config.sigma, config.a, config.b))` → `(... config.a, config.s)`. The dispatch logic itself is updated in Task 8 (Task B); this task is just the rename.

### Task 6: Update `misc_tools.jl`

**Files**: `src/misc_tools.jl`

**Action**:
- `_generate_filename(config::Config{Lindbladian})` (line 76-87): `b_str = "b=$(config.b)"` → `s_str = "s=$(config.s)"`. Update `join([..., a_str, b_str], "_")` → `join([..., a_str, s_str], "_")`. **User-facing change**: BSON filenames will now contain `s=0.4` instead of `b=0.4`.
- `_generate_filename(config::Config{Thermalize})` (line 89-101): same.
- `validate_config!` (line 151-158): keep the eta-requirement; remove the `b != 0` rejection (we'll add `a = 0, s > 0` support in Task 8). For now, in the rename-only commit, replace with:
  ```julia
  if config.with_linear_combination && config.a == 0.0
      # Note: the (a=0, s != 0) case is supported in eta-regularized smooth Metro (Task 8).
      if config.domain isa Union{TimeDomain, TrotterDomain} && with_coherent(config.construction) && (isnothing(config.eta) || config.eta <= 0.0)
          push!(errors, "For linear combinations in the KMS DB case with a=0 in TIME or TROTTER domain, eta must be > 0.")
      end
  end
  ```
  The error wording stays factual; the inner `b != 0` rejection is dropped because nothing now consumes `s != 0` in the `a = 0` branch — `_select_b_plus_calculator` falls through to the `_compute_b_plus_metro` path (which still ignores `s`). After Task 8 it will start using `s`.
- `_print_press` (line 218-235 and 248-267): `("b", config.b)` → `("s", config.s)`.

### Task 7: Update `results.jl`

**Files**: `src/results.jl`

**Action**: line 40 — `d[:b] = config.b` → `d[:s] = config.s`. In `_dict_to_config_kwargs` (line 104), the loop tuple `(:a, :b, :num_energy_bits, :t0, :w0, :eta, :num_trotter_steps_per_t0)` becomes `(:a, :s, :num_energy_bits, :t0, :w0, :eta, :num_trotter_steps_per_t0)`.

**Backward compatibility note**: Existing BSON results saved with `:b` key will fail to load. Since these are dev-internal results and we own all of them, this is acceptable. The user explicitly stated "deliberate user-facing change". Flag in commit message that legacy `.bson` files saved with `:b` key cannot be loaded after this commit.

### Task 8: Add eta-regularized smooth-Metro `b_+` (Task B)

**Files**: `src/coherent.jl`, `src/furnace_utensils.jl`, `src/structs.jl` (already updated in commit 1)

**Action 8a** — `src/coherent.jl`: replace `_compute_b_plus_metro(t, beta, sigma, eta)` with the four-arg version:
```julia
function _compute_b_plus_metro(t::Real, beta::Real, sigma::Real, eta::Real, s::Real=0.0)
    if abs(t) < 1e-12
        # Analytical t=0 limit (L'Hopital); reduces to 1/(2√2 π²) at σβ=1, s=0.
        return complex((2 - sigma^2 * beta^2 * (1 + s)) / (2 * sqrt(2) * pi^2))
    elseif abs(t) <= eta
        numerator = exp(-sigma^2 * beta^2 * (2 * t^2 + 1im * t) * (1 + s)) + 1im * (2 * t + 1im)
    else
        numerator = exp(-sigma^2 * beta^2 * (2 * t^2 + 1im * t) * (1 + s))
    end
    denominator = t * (2 * t + 1im)
    return (1 / (2 * sqrt(2) * pi^2)) * numerator / denominator
end
```

**Action 8b** — `src/furnace_utensils.jl` `_select_b_plus_calculator` (line 143-156):

```julia
function _select_b_plus_calculator(config::Config{<:Any, <:Any, KMS})
    if !config.with_linear_combination
        return (_compute_b_plus, (config.beta, config.gaussian_parameters[1], config.gaussian_parameters[2]))
    elseif config.a != 0.0
        # a-regularized smooth (covers Metro a>0,s=0 and Glauberish a>0,s>0)
        return (_compute_b_plus_smooth, (config.beta, config.sigma, config.a, config.s))
    else
        # eta-regularized smooth Metro (covers Chen plain s=0 and thesis-main s>0)
        s_val = something(config.s, 0.0)
        return (_compute_b_plus_metro, (config.beta, config.sigma, config.eta, s_val))
    end
end
```

The default `s::Real = 0.0` argument in 8a means the test `_compute_b_plus_metro(t, beta, sigma, eta)` (4-arg) still works on its own. Together with 8b passing 5 args now, the `s = 0` default is only relevant for direct callers.

### Task 9: Update tests

**Files**: `test/test_helpers.jl`, `test/test_simulation_time.jl`

**Action 9a** — `test/test_helpers.jl` line 201: `b = 0.4,` → `s = 0.4,` (inside `make_config`).

**Action 9b** — `test/test_simulation_time.jl` line 22: `Dict{Symbol,Float64}(:a => 0.333, :b => 0.4)` → `... :s => 0.4`.

**Note**: tasks 1–9 happen IN A SINGLE COMMIT. Splitting across commits would
require leaving stale `b = 0.4` in test_helpers temporarily, breaking tests. The
user's "split into a few atomic commits along logical lines" suggestion was
exploratory; in practice the rename touches struct + every consumer simultaneously.

### Task 10: Update scripts

**Files**: 7 scripts under `scripts/`

**Action**: each `b = 0.4` keyword arg → `s = 0.4`. Comment in `plot_transition_weights.jl:28` updated too.

**Decision**: Scripts are not run by tests — they're called manually. Updating them
at the same time is necessary because they will fail the next time someone runs
them with the new code. Bundle into the same commit as Task 1–9.

### Task 11: Add tests for new dispatch path

**Files**: `test/test_simulation_time.jl` or a new section in an appropriate place

**Action**:
1. Test that `_compute_b_plus_metro(t, β, σ, η, 0.0)` reproduces the legacy formula for several `t`. Use the explicit reference inline:
   ```julia
   function ref_b_plus_metro_legacy(t, beta, sigma, eta)
       if abs(t) < 1e-12
           return complex(1 / (2 * sqrt(2) * pi^2))
       elseif abs(t) <= eta
           numerator = exp(-sigma^2 * beta^2 * (2*t^2 + 1im*t)) + 1im * (2*t + 1im)
       else
           numerator = exp(-sigma^2 * beta^2 * (2*t^2 + 1im*t))
       end
       denominator = t * (2*t + 1im)
       return (1 / (2 * sqrt(2) * pi^2)) * numerator / denominator
   end
   # Skip t=0 in the ≈ comparison: legacy hard-code only matches at σβ=1.
   for t in [-1.0, -0.5, -0.05, 0.05, 0.5, 1.0]
       new = QuantumFurnace._compute_b_plus_metro(t, BETA, SIGMA, 0.05, 0.0)
       legacy = ref_b_plus_metro_legacy(t, BETA, SIGMA, 0.05)
       @test new ≈ legacy atol=TOL_EXACT
   end
   # Verify the t=0 case at σβ=1 (i.e. our test default)
   @test QuantumFurnace._compute_b_plus_metro(0.0, BETA, SIGMA, 0.05, 0.0) ≈ complex(1 / (2*sqrt(2)*pi^2)) atol=TOL_EXACT
   ```
2. Test that `_select_b_plus_calculator` returns `_compute_b_plus_metro` (with 4 fixed args including s_val) for `a = 0, s = 0.4`:
   ```julia
   config_a0_spos = Config(; sim=Thermalize(), domain=TimeDomain(), construction=KMS(),
       num_qubits=NUM_QUBITS, with_linear_combination=true,
       beta=BETA, sigma=SIGMA, a=0.0, s=0.4, eta=0.05,
       num_energy_bits=NUM_ENERGY_BITS, w0=W0, t0=T0,
       num_trotter_steps_per_t0=NUM_TROTTER_STEPS_PER_T0,
       mixing_time=1.0, delta=TEST_DELTA,
   )
   validate_config!(config_a0_spos)  # should NOT throw
   bp_fn, bp_args = QuantumFurnace._select_b_plus_calculator(config_a0_spos)
   @test bp_fn === QuantumFurnace._compute_b_plus_metro
   @test bp_args == (BETA, SIGMA, 0.05, 0.4)
   ```
3. Smoke test: precompute with this config and confirm `b_plus` Dict is non-empty.
   ```julia
   pd = QuantumFurnace._precompute_data(config_a0_spos, TEST_HAM)
   @test pd.b_plus !== nothing
   @test !isempty(pd.b_plus)
   ```

**Place these tests** in a new file `test/test_smooth_metro_eta.jl` and include
it in `test/runtests.jl`. New test file gives a clean, named test set rather
than scattering across other files.

### Commit ordering

1. **Commit 1**: `refactor(config): rename b → s parameter to match thesis notation` — all 10 tasks above (1–10). Pure rename, no behavior change. Tests pass.
2. **Commit 2**: `feat(coherent): add eta-regularized smooth Metropolis b_+ for a=0, s>0` — Tasks 8a + 8b. Behavior change: t=0 limit now general; `a=0, s>0` now valid. Tests still pass with old (s=0) parameters.
3. **Commit 3**: `test(coherent): cover a=0, s>0 dispatch path and s=0 reduction` — Task 11. Adds new test file.

After commit 1: same 1273 tests pass.
After commit 2: same 1273 tests pass (no test exercises the new path yet).
After commit 3: 1273 + N new tests pass where N is number of @test calls in
`test_smooth_metro_eta.jl` (estimated ~10).

## Must-Haves

Each is observable-from-code or test-runnable:

- [ ] `grep -rn "config\.b\b" src/ test/` returns 0 hits
- [ ] `grep -rn "::Real, b::Real" src/` returns 0 hits (function param `b::Real` from old `create_alpha` etc. — only `_compute_b_plus_smooth` had that signature)
- [ ] `grep -rn "\bb\s*=\s*0\." scripts/ test/` returns 0 hits (excluding loop body comments)
- [ ] `grep -n "config\.s\b" src/` returns >= 12 hits
- [ ] `grep -rn "validate_config" src/misc_tools.jl` no longer rejects `b != 0` / `s != 0` when `a == 0`
- [ ] `julia --project -e 'using Pkg; Pkg.test()'` passes after each commit (1273 → 1273 → 1273+N)
- [ ] New test file `test/test_smooth_metro_eta.jl` exists and runs
- [ ] `_compute_b_plus_metro` accepts 5 positional args and reduces to legacy form when 5th arg is 0.0
- [ ] `_select_b_plus_calculator` returns `_compute_b_plus_metro` with 4 args (including `s_val`) for the `a=0, s>0` config

## Risk Flags

- **BSON filename change** (intentional): old result files saved as `..._b=0.4_...bson` will not match new naming. User has confirmed acceptable.
- **BSON `:b` key change**: `_dict_to_config_kwargs` reading `(:a, :s, ...)` from a Dict that was saved with `:b` key will silently drop the `b` value, leaving `config.s = nothing`. Loaded configs will have wrong `s`. **Mitigation**: this is a dev-only repo; we accept this. If anyone needs to load legacy results, do it before committing.
- **t=0 limit change in `_compute_b_plus_metro`**: For `σβ ≠ 1`, the `s = 0` reduction at `t = 0` differs from the old hard-coded value. New value is mathematically correct; old was a coincidence at `σβ=1`. Affects scripts using `σ = 0.8/β` (in `biexp_mixing_verify_trotter.jl`). **Mitigation**: physics check — the t=0 limit only contributes one Dict entry; the L'Hopital formula is provably the analytic limit; the `_get_truncated_indices` step might drop or keep this one point but the integrated effect on simulation outputs is at most one more truncated time-step contribution.
- **Loop var rename `s → tau` in `B_time`/`B_trotter`**: pure local-variable rename, no exported names involved. Verifiable by tests.
- **Default arg `s::Real = 0.0` in `_compute_b_plus_metro`**: ensures direct 4-arg callers (e.g. anywhere outside `_select_b_plus_calculator`) keep working. Need to check that no test or script calls `_compute_b_plus_metro` directly with 4 args — if they do, behavior unchanged. **Verified by `grep`**: no direct calls outside `_select_b_plus_calculator` and tests using `bp_fn(t, args...)` form.
- **`with_linear_combination = true, a = 0, s = 0`** path: now reaches `_compute_b_plus_metro(beta, sigma, eta, 0.0)`. Since `s = 0.0` reduces exactly to the original formula, behavior unchanged for this legacy case.

## Test Strategy

- **Existing 1273 tests pass after every commit** — each commit is verified by `julia --project -e 'using Pkg; Pkg.test()'`.
- **New tests** in `test/test_smooth_metro_eta.jl`:
  1. `s = 0` reduction matches legacy form (8 points sampled across the kernel, including `t = 0`)
  2. Dispatch returns `_compute_b_plus_metro` for `a = 0, s > 0` config
  3. `validate_config!` does NOT throw for `a = 0, s > 0` (when `eta > 0`)
  4. `validate_config!` DOES throw for `a = 0, s > 0, eta = 0` (or nothing) in TimeDomain
  5. `_precompute_data` returns non-empty `b_plus` for `a = 0, s > 0` config
  6. End-to-end smoke: a one-step `run_thermalize` doesn't error with `a = 0, s > 0`

- **Physics check** (will spawn `physics-checker` agent if user requests): the
  thesis main case `a = 0, s > 0` should give a smooth `b_+(t)` (less kink than
  `s = 0` Chen Metro). The eta-regularized denominator handles `t → 0`; the
  `(1 + s)` in the exponent makes the Gaussian decay faster, reducing the
  kink-induced quadrature error.
