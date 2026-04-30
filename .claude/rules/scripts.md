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
