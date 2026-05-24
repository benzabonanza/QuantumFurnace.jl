---
name: feedback-canonical-rho0-always-plus
description: "POLICY (2026-05-24, user-decreed) — ρ_0 = |+⟩⟨+|^⊗N for ALL simulations. NEVER I/d. NO exceptions."
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 2e0ec557-e0b4-44c1-92de-7a08e13be5d7
---

For **every** Lindbladian/channel trajectory simulation in this project — `predict_lindbladian_trajectory`, `predict_channel_trajectory`, any τ_mix computation, any thesis numerics — the initial state is:

$$\rho_0 = |{+}\rangle\langle{+}|^{\otimes N},
\quad\text{where}\quad |{+}\rangle = \frac{|0\rangle + |1\rangle}{\sqrt{2}}.$$

**Never use `I/d`** (maximally mixed state). Never use other "default" ρ_0 from scratch scripts (e.g. `:thermal_perturbed`) without explicit user instruction.

**Why:** the user decreed this 2026-05-24 after noticing the qf-1jj 2D TFIM sweep had used `I/d` (predating the qf-e4z.30 canonical-ρ_0 decision for 1D Heisenberg) and produced a misleading τ_mix vs Gap plot. They are frustrated by inconsistent ρ_0 across sweeps. Decreed: |+⟩⟨+|^⊗N is the SOLE canonical ρ_0 going forward.

**How to apply:**

1. Default any new sweep / trajectory driver to `ρ_0 = build_plus_state(num_qubits)` where:
   ```julia
   function build_plus_state(n)
       psi = ones(ComplexF64, 2^n) ./ sqrt(2.0^n)   # |+⟩^⊗N = ∑_x |x⟩ / √(2^N)
       return psi * psi'                            # rank-1 d × d density matrix
   end
   ```
2. **Audit any existing script before reusing it as a template** — if it has `rho_0 = I/d` or `Matrix{ComplexF64}(I, d, d) ./ d`, fix it.
3. When inheriting a sweep / driver pattern, explicitly check the `rho_0 = ...` line. Do not let the project default of `:maximally_mixed` in `_make_init_state` silently apply.
4. Don't use `_make_init_state(:thermal_perturbed, ...)` without explicit user instruction either — that's a different policy decision (it was discussed in qf-1jj methodology notes but never adopted).

**What this does NOT mean**:
- It does NOT mean we're avoiding all symmetry-protected initial states. |+⟩⟨+|^⊗N is still Z₂-even (`X|+⟩=|+⟩`), so it remains symmetry-protected against the Z₂-odd half of slow modes in Z₂-symmetric Hamiltonians (e.g. 2D TFIM). But it IS the canonical thesis choice (qf-e4z.30) and gives consistent τ_mix values across sweeps, which is the user's priority.

Connected: [[canonical-taumix-setup-qf-e4z-30]] (where this was first decided for 1D Heisenberg multiseed), [[heisenberg-1d-multiseed-even-odd-qf-e4z-23]] (where I/d caused the parity trap), [[qf-biz-2d-tfim-matrix-element-refuted]] (the qf-biz analysis that used the qf-1jj I/d sidecars and led to the user catching this).
