---
name: feedback-check-gibbs-when-simulation-off
description: When a simulation result looks surprising or "off" (unexpected scaling, anomalous gap, weird τ_mix), check the Gibbs occupations p_i = exp(-β·ΔE_i)/Z at the (n, β) cell before reaching for algorithmic explanations.
metadata:
  type: feedback
---

# Check the Gibbs state first when a simulation looks off (set 2026-05-17)

**Rule:** When a simulation result is surprising — unexpected scaling, anomalous gap, weird τ_mix, "is this really low temperature?" intuition — *before* reaching for algorithmic explanations (parity traps, Krylov truncation, register sizing), **directly compute the Gibbs-state occupation profile** `p_i = exp(-β_phys · (E_i^phys − E_0^phys)) / Z` at the relevant (n, β_phys) cell.

**Why:** the physics interpretation of a gap_phys number depends critically on *which* states the Gibbs state actually populates. "β_phys = 2.5" sounds cold, but the actual coldness depends on the model's energy gap ΔE_1 in physical units — could mean p_0 = 99 % (essentially pure GS) for one system and p_0 = 28 % (broad multiplet) for another. The interpretation of the slowest L-mode is completely different in those two regimes.

**How to apply:** load the relevant `HamHam`, work in algorithm-frame eigenvalues (they're equivalent for Gibbs: `β_phys · ΔE_phys = β_alg · ΔE_alg`), compute `p_i` and the effective rank `1/Σ p_i²`. Boilerplate that takes 30 seconds:

```julia
using QuantumFurnace, Printf
raw = QuantumFurnace._parse_hamiltonian_bson(path_to_fixture)
ham = HamHam(raw; beta_phys = β_phys)
ΔE = (ham.eigvals .- minimum(ham.eigvals)) .* ham.rescaling_factor   # phys units
W  = exp.(-β_phys .* ΔE)
p  = W ./ sum(W)
@printf("p_0=%.4f  p_1=%.4f  eff_rank=%.2f  ΔE_1_phys=%.4f\n",
        p[1], p[2], 1/sum(p .^ 2), ΔE[2])
```

This came up in 2026-05-17 while interpreting the qf-e4z.23 even/odd-n split: the original "thermal cluster of low-energy states" mechanism we had attached to the even-n gap collapse was wrong — at β_phys=2.5 the even-n Gibbs is essentially pure GS (p_0 = 99.7 % for n=6). The actual mechanism (GS↔first-excited coherence bottleneck for even n; intra-multiplet relaxation for odd n) only became visible after checking the Gibbs occupations. See [[heisenberg-1d-even-odd-mechanism]] for the full worked example.

**General lesson:** for KMS-DB Lindbladian work, the Gibbs occupation profile is the cheap-to-compute, ground-truth physics input that disambiguates "thermal cluster vs unique GS vs degenerate multiplet" — three regimes that give qualitatively different gap mechanisms. Compute it before writing the physics interpretation, not after.
