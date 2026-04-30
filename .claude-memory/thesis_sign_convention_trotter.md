---
name: Thesis sign convention mismatch (Trotter / Hamiltonian evolution)
description: The preliminaries and CKG-methods chapters use opposite signs for exp(iHt) vs exp(-iHt); likely needs unification across the whole thesis.
type: project
originSessionId: 15982775-9c8b-4b29-8a7a-595d3e8cafce
---
The thesis currently has a sign-convention mismatch for Hamiltonian evolution that the user may want to unify in a later pass.

**Observation:**
- `supplementary-informations/1_preliminaries.tex` (QPE subsection and onward, including the new Trotterization subsection `sec:prelim-trotter`): uses $U = \ee^{+\ii H t_0}$ in the QPE circuit and $S_p(t) \approx \ee^{-\ii H t}$ in the Trotter subsection. (Note: the prelims itself is inconsistent — QPE uses `+` so that the readout register labels $\nu$ directly via (eq:qpe-unified); the Trotter subsection uses the more common `-` for product formulas.)
- `supplementary-informations/2_methods.tex` (CKG `Hamiltonian Trotterization` section, line 889 onward): uses $\ee^{+\ii H t}$ everywhere — see (eq:childs-trotter-bound).

**Why:** The OFT construction in Chen et al. 2023 uses $U = \ee^{+\ii H t_0}$ so that the Fourier kernel's phase matches the Bohr-frequency labelling in the discretised OFT; this propagates through the CKG generator construction. The Trotter literature (including Childs et al. 2021) universally writes product formulas for $\ee^{-\ii H t}$. The mismatch is cosmetic — Childs's commutator bound is symmetric under $t \mapsto -t$ — but visually jarring inside the thesis.

**How to apply:** The Trotter-subsection draft (drafts/trotterization.md as of 2026-04-21) now adopts $\ee^{-\ii H t}$ per Trotter-literature convention and does not paper over the mismatch in the body. The old v2 draft had a standalone explanatory sentence resolving it; at the user's request this was removed so the user can make a deliberate thesis-wide choice later. Two resolution options when doing that pass:
  (a) Unify to $\ee^{-\ii H t}$ everywhere (change `2_methods.tex:889` and propagate into the Childs bound display; also adjust the QPE circuit convention and rewire (eq:qpe-unified) so the readout register labels $-\nu$ or flip the filter phase).
  (b) Unify to $\ee^{+\ii H t}$ everywhere (change the Trotter subsection to use $S_p(t) \approx \ee^{+\ii H t}$; Childs's bound still applies since it is $t \mapsto -t$ symmetric). Probably less invasive.

Flag this when the user is doing a thesis-wide consistency sweep.
