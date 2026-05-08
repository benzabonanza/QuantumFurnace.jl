---
name: Multi-rank DLL completion notes (qf-7go epic)
description: Outcome of the parametric multi-channel DLL filter epic — symmetrised translates, KMS DBC verified, τ_mix vs k headlines.
type: project
---

## Epic qf-7go: closed 2026-05-03

All 7 sub-issues closed. 7 atomic commits dd0f732..035bad1 (DLL-MR.1..5),
plus 64f6c04 (sweep harness) and follow-up plot+draft commit.

## Construction

`dll_multichannel_translates(base; centers, weights)` builds a
`DLLMultiChannelFilter` whose channels are symmetrised translates:
`q_ℓ(ν) = √(w/2)·[q_base(ν - ν_ℓ) + q_base(ν + ν_ℓ)]` (or `√w · q_base`
when ν_ℓ = 0). Each q_ℓ is real-even by construction → KMS Eq. 3.2
holds per channel → KMS-DBC of α^multi follows.

**KMS-DBC of the coherent term G**:
By Theorem 10 (Eq. 2.33), `G = -i tanh ∘ log(Δ^{1/4}) · (½ Σ_j L_j† L_j)`.
The map `V → -i tanh ∘ log(Δ^{1/4}) · V` is *linear* in V, so for the flat
multi-channel set {L_a^(ℓ)}, `G^multi = Σ_{a,ℓ} G_{a,ℓ}` — exactly what
`dll_coherent_op_bohr(::DLLMultiChannelFilter, β)` returns. Numerical
witness: `verify_detailed_balance` gives `‖A_anti‖₂ / ‖D‖₂ ≤ 1e-10` for
k ∈ {1, 2, 4} at β ∈ {1, 5, 10} (test_dll_multichannel_simulator.jl
testset (g.0)).

## Library API

- `DLLMultiChannelFilter{T, F<:AbstractFilter}(channels, beta)` — wraps k
  channels; constructor enforces β-consistency.
- `ShiftedSymmetricFilter{T, F<:AbstractFilter}(base, shift, weight, beta)`
  — single shifted-symmetric channel; closed-form q_l, f̂_l, f_l, cutoff.
- `dll_multichannel_translates(base; centers, weights)` — factory; rejects
  empty centers, non-positive weights, |center| > base.S/2 for Metropolis.
- BohrDomain ops `dll_lindblad_op_bohr(::DLLMultiChannelFilter)` returns
  `Vector{Matrix}` (no cross terms in dissipator); `dll_coherent_op_bohr`
  returns single sum-over-channels Matrix.
- TimeDomain analogues via NUFFT per channel; Vector for L, single Matrix
  for G.
- Simulator dispatch via `Config{Lindbladian, *, DLL}` with the
  multi-channel filter — `_jump_contribution!` loops per channel.
- Matrix-free Krylov path: `_accumulate_R_total_dll!` flattens
  per-channel L's into the existing `dll_lindblads::Vector{Matrix}`
  workspace slot; `apply_lindbladian!` already iterates this vector.

## Refactor: `oft_nufft_at_zero` → `oft_nufft_at_zero_list`

Time-domain DLL precompute (`_precompute_data` in
`src/furnace_utensils.jl`) now returns a **`Vector{Matrix}`** of per-channel
NUFFT prefactors. Length-1 for single-channel DLL filters (preserves
hot-path performance), length-k for multi-channel. `_jump_contribution!`
(TimeDomain DLL) loops over the list and accumulates per-channel
dissipators. **Test suite update**: `test_dll_dissipator.jl` test (h)
updated to `pre.oft_nufft_at_zero_list[1]`.

## Headline empirical (qf-7go.6 sweep, n=3, target_ε=1e-3)

DLL Metropolis (S=2, centers uniform on [0, S/2]):
| β  | k=1  | k=2  | k=4  | k=8  | k=1→k=8 |
|----|------|------|------|------|---------|
| 1  | 11.2 | 4.56 | 1.80 | 0.82 | 13.7×   |
| 5  | 18.1 | 13.6 | 5.06 | 2.22 | 8.2×    |
| 10 | 23.3 | 22.2 | 10.3 | 4.52 | 5.2×    |
| 20 | 35.4 | 35.3 | 26.4 | 12.9 | 2.7×    |

DLL Gaussian (centers uniform on [0, 1.0]):
| β  | k=1  | k=2  | k=4  | k=8  | k=1→k=8 |
|----|------|------|------|------|---------|
| 1  | 6.79 | 2.65 | 1.08 | 0.50 | 13.7×   |
| 5  | 11.8 | 11.4 | 4.70 | 2.02 | 5.9×    |
| 10 | 25.3 | 25.3 | 9.69 | 4.36 | 5.8×    |
| 20 | **264**  | **264**  | 20.8 | 13.3 | **19.9×** |

**Key takeaway**: The β=20 DLL Gaussian collapse documented in
`ckg_vs_dll_first_findings.md` is **fully recovered by adding rank** —
not a structural limitation of the Gaussian-type weighting, but an
artefact of insisting on rank 1.

## Revises §8.7 prediction

§8.7 of `drafts/ckg-vs-dll-comparison-findings.md` predicted multi-rank
would give a null result on τ_mix because H2 (diagonal-of-α drives the
rate) implied that adding channels at the same Metropolis shape would
not help. **§9 contradicts this**: shifted-translate channels enrich the
diagonal at multiple ν values, and the diagonal *shape* matters, not
just the per-frequency value at one ν. H2 is consistent — but the
optimal diagonal is not unique to CKG, it's recoverable inside DLL.

## Cost caveat

§5 / Remark 23 cost picture inherits per-channel: `‖G^multi‖ ≤ Σ_ℓ
‖G^(ℓ)‖ = O(k · β^{1.4})`. So the τ_mix gain at k=8 vs k=1 nets ~ratio/8
in total cost. β=20 Gaussian: 19.9×/8 = 2.5× total-cost reduction —
real but smaller than the τ_mix table alone suggests. Verdict §6 stands:
CKG remains the better choice for low-T large-system Gibbs sampling
under the [Ding et al. 2024] Table 1 cost metric.

## Files

- `src/filters.jl`: +DLLMultiChannelFilter, ShiftedSymmetricFilter,
  dll_multichannel_translates.
- `src/dll.jl`: multi-channel `dll_lindblad_op_bohr` (Vector return),
  `dll_coherent_op_bohr` (single sum), and TimeDomain analogues.
- `src/furnace_utensils.jl`: `_precompute_data` for TimeDomain DLL now
  returns `oft_nufft_at_zero_list::Vector{Matrix}`.
- `src/jump_workers.jl`: `_jump_contribution!` (Bohr+Time) loops per
  channel via new `_accumulate_dll_bohr_dissipator!` helper.
- `src/krylov_workspace.jl`: `_accumulate_R_total_dll!` flattens
  per-jump per-channel L's into workspace.
- `src/misc_tools.jl`: validate_config! per-channel β + S checks.
- `test/test_dll_multichannel_filter.jl`, `test_dll_multichannel_bohr.jl`,
  `test_dll_multichannel_time.jl`, `test_dll_multichannel_simulator.jl`:
  +130+365+32+34 tests = ~561 new tests across the four new files.
- `scripts/scratch_dll_multirank_taumix_sweep.jl`: 32-cell sweep harness.
- `scripts/plot_dll_multirank_taumix.jl`: τ_mix(k) panel plot.
- `drafts/dll-multirank-taumix-findings.md`: §9 thesis draft.
- `drafts/figures/numerics/dll_multirank_taumix.{png,pdf}`: figure.
- `scripts/output/dll_multirank/taumix_vs_k.bson`: sweep BSON.

## Future work (not in epic)

- Spectrum-adapted centers (place ν_ℓ at observed Bohr clusters), per-
  channel weights w_ℓ, mixed base shapes (some Metropolis, some Gaussian).
- Analytical prediction of optimal centers via the qf-mto Dirichlet-form
  framework — does the multi-rank DLL "diagonal shape" admit a closed
  form?
- n ≥ 4 sweep — current sweep is n=3 only.
- TrotterDomain support (deferred per qf-3i8 scope).
