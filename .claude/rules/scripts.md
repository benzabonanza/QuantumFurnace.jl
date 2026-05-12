---
paths:
  - "scripts/**/*.jl"
---

# Script Rules

- Scripts in `scripts/` are standalone runnable: `julia --project scripts/<name>.jl`
- Scripts get their own commits separate from library code in `src/`
- Include comments explaining the algorithm strategy (not every line)
- Include a runnable main block that demonstrates correctness
- Print key results to stdout
- Flag physics decisions with `# PHYSICS CHECK:` comments
- Use existing deps from the project; don't add new ones without justification
- **Coherent term ON by default.** When calling `apply_lindbladian!`, `construct_lindbladian`, `predict_lindbladian_trajectory`, `krylov_spectral_gap`, etc., keep the default `include_coherent = true` — that is the physical KMS Lindbladian. Only pass `include_coherent = false` when you are explicitly isolating the dissipator for a register-sizing or quadrature-error diagnostic, and the user has asked for that. When you do: (1) state it in the script header docstring, (2) name the script and outputs accordingly (`*_dissipator_*`, "dissipator-only sweep"), and (3) call out the caveat in any markdown summary / draft you produce. A dissipator-only sweep does NOT size the coherent registers `r_b_minus, r_b_plus`; see the [Coherent Term: ON by Default](../rules/julia-code.md#coherent-term-on-by-default) rule for the full reasoning.
