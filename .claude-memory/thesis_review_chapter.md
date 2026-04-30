---
name: Thesis Review chapter — quantum vs classical Gibbs sampling
description: Plan for the thesis "Review" chapter framing where quantum Gibbs sampling could first beat classical samplers; tracks /lit-review infrastructure and corpus-building plan
type: project
---

# Thesis "Review" Chapter — Quantum vs Classical Gibbs Sampling

**Beads epic**: `qf-yk4` — full plan and paper list lives there. Run `bd show qf-yk4` for the canonical description; child tasks `qf-yk4.1`–`qf-yk4.3`.

**Why**: Thesis needs a chapter ("Review") that maps the state of the art on Gibbs sampling — quantum and classical — and frames the central question: **where could quantum Gibbs sampling first beat a classical sampler?** Goal is a regime map (temperature × commutativity × locality × particle statistics) showing best known classical/quantum runtimes per cell and whether quantum advantage is proven, conjectured, ruled out, or unknown.

**How to apply**:
- When asked about the Review chapter, the literature corpus, or the quantum-vs-classical landscape — the plan is in `qf-yk4`. Use `/lit-review <paper>` to add a paper to the corpus.
- Reviews land in `drafts/literature-review/<slug>.md` with mandatory regime YAML frontmatter (`temperature`, `commutativity`, `locality`, `particle-statistics`, `paradigm`, `quantum-or-classical`, `result-type`, `key-scaling`, `related`).
- Synthesis skill (`/lit-synthesize`) is **not yet built** — task `qf-yk4.2`. Build it once the corpus has ~5–10 reviews.
- Final composition pass (task `qf-yk4.3`) is `thesis-writer` pulling the synthesized regime overviews into thesis prose anchored on regime-map figure(s).

**Scope of "classical side"**: practitioners' actual current methods — Glauber/Metropolis MCMC, cluster algorithms (Swendsen-Wang, Wolff), tensor networks (DMRG/PEPS finite-T, METTS), QMC/sign-problem-free results, parallel tempering. These need to be sourced beyond `supplementary-informations/`.

**Quantum-classical comparison framing** (mandatory for every review): why is the quantum version harder/easier here, and what is the parameter gap (poly / super-poly / unknown)? Source of the difference (sign problem, noncommutativity, KMS condition, locality of Lindblad operators) must be identified.

**Infrastructure** (already in place):
- `.claude/agents/lit-reviewer.md` — agent (opus, effort: high)
- `.claude/skills/lit-review/SKILL.md` — `/lit-review <paper>` skill (forks lit-reviewer)
- `drafts/literature-review/` — output directory
