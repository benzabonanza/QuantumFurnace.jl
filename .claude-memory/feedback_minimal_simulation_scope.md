---
name: feedback-minimal-simulation-scope
description: "Default to the minimum simulation that answers the physics question — do not inherit knobs from existing sweep scripts (multi-kdim, multi-seed, Bohr cross-check) by reflex."
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 2e0ec557-e0b4-44c1-92de-7a08e13be5d7
---

When the user asks for a small physics diagnostic, write a **minimal-scope** simulation that answers only the specific question — do not inherit knobs from existing sweep harnesses by reflex (multiple krylovdim values, multiple seeds, BohrDomain dense cross-check, trajectory τ_mix predictions, β-sweeps when a single β is asked) when they are not strictly needed for the answer.

**Why:** stated explicitly 2026-05-24 on qf-biz ("You sometimes simulate more than necessary because of existing sweep scripts, like many krylovdims, or what not, but I want you to only simulate what's needed to figure out the answer to these physics questions"). The user also added "I dont want a heavy simulation at the moment" — they value cycle-time and clarity of the per-question scope over completeness.

**How to apply:**

- Before writing any new diagnostic script, read the relevant existing sweep / driver script and explicitly list the knobs (krylovdim, howmany, register sizes, β grid, seeds, BohrDomain dense, trajectory predictor, multi-domain cross-checks).
- For each knob, ask: does this question genuinely need it, or is it inheritance from the bigger sweep? Default-trim — keep only what the question demands. Document the trimmed knobs in the script header docstring as "Skipped deliberately (per scope): …".
- A canonical CKG smooth-Metro krylov_spectral_gap at the validated kdim=40, howmany=4 is enough for almost any single-cell gap question — multiple kdim is a register-sizing diagnostic, not a routine cross-check.
- BohrDomain dense at n=6 costs ~165 s/cell — never run it as a "safety" cross-check when the same configuration already has a clean Bohr-vs-Energy validation in a prior issue (qf-1jj has this for n=6 ORDERED at r_D=7).
- Trajectory `predict_lindbladian_trajectory` is for τ_mix questions only — irrelevant to "is the gap due to mechanism X" questions.

Connected memories: [[canonical-taumix-setup-qf-e4z-30]] (default kdim, single-pass settings), [[feedback-more-data-points-for-scaling-claims]] (only fit asymptotic exponents with enough cells — orthogonal but related sentiment about not over-claiming from a sweep), [[feedback-check-gibbs-when-simulation-off]] (cheap Gibbs check first before invoking algorithmic explanations).
