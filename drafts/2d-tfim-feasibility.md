# 2D Quantum TFIM Pipeline Smoke Test — Feasibility Note (qf-833)

**Status**: complete — all 3 sandbox cells finished in 16.3 min wall (1-hour budget).
**Driver**: `scripts/scratch_2d_tfim_gap_sweep.jl`
**Sidecars**: `scripts/output/sweep_2d_tfim/tfim_n{n}_Lx{Lx}_Ly{Ly}_h1_betaphys0.5_{Domain}_clean.bson`
**Summary**: `scripts/output/sweep_2d_tfim/sweep_2d_tfim_summary.bson`
**Log**: `scripts/output/sweep_2d_tfim/sweep.log`

> **Retrospective correction (2026-05-14, physics-checker review for qf-1jj)**:
> the "disordered phase" classification of the $(h = 1, \beta_{\rm phys} = 0.5)$
> operating point used throughout this note is **incorrect** — strictly,
> $T = 2.0 J$ at $h = J$ sits at most $\sim 3\%$ BELOW $T_c(h{=}1) \approx 2.07 J$
> per Hesselmann–Wessel 2016 (arXiv:1602.02096) SSE-QMC, putting us narrowly
> INSIDE the ordered "dome", not above it. The numerics ($\lambda_{\rm gap}$,
> $\tau_{\rm mix}$, Energy↔Bohr cross-checks, RSS) are all correct as
> measured — but the flat $\lambda_{\rm gap,phys} \approx 0.285$ across n=6
> and n=8 is **pre-asymptotic finite-size behaviour, NOT the expected
> $\Omega(1)$ disordered-phase gap**. The sandbox cells $n \le 8$ have
> $\sqrt{n}$ range $0.83$ — far too narrow to resolve the asymptotic
> stretched-exponential $\exp(-\gamma \sqrt{n})$ decay that's predicted
> for any Lindbladian quasi-local Gibbs sampler in the ordered phase
> (Lubetzky–Sly 2013 classical; Gamarnik–Kiani–Zlokapa 2024 quantum lower
> bound $T_{\rm mix} \ge 2^{n^{1/2 - o(1)}}$). Cite this note as
> *"narrowly inside the ordered phase at the critical line"*, not as
> *"disordered phase"*. The qf-1jj follow-up runs at $(h = 1.5,
> \beta_{\rm phys} = 2.0)$ which is unambiguously inside the dome
> (73% margin below $T_c$); compare against that. Phase-diagram
> reference: [[T_c(h) for 2D TFIM]].

## Headlines

- **All three sandbox cells completed** (n=4, 6, 8) inside 16.3 min — comfortably under the 1 h budget. n=8 dominated at 768 s; n=4 and n=6 are under 10 s each on EnergyDomain.
- **Energy↔Bohr cross-check passes at machine precision** at both n=4 and n=6: rel = 4.80e-14 and 8.20e-14. r_D = 7 with the unified ω-range is ~7 orders of magnitude over-converged for the 1e-7 gate. **No symmetry contingency triggered** — the qf-8fr Krylov-x0 fix (1e-10 GUE kick, [[krylov_x0_symmetric_bug_qf_8fr]]) survives clean 2D TFIM at β_phys = 0.5 without modification.
- **Pipeline is correct**: gap_phys = gap_alg × rescaling_factor lands at $\approx 0.29$ for both bulk-like cells (n=6, n=8). The n=4 (2×2) cell is a PBC-double-bond outlier ($g_{\rm phys} = 0.45$) because the wrap-around in $L_x = 2$ produces effective coupling 2J — it's barely 2D.
- **Dense BohrDomain `construct_lindbladian + eigvals` is the right cross-check path at n ≤ 6.** Matrix-free Krylov in BohrDomain is $\mathcal O(\text{krylovdim} \cdot d^4 \cdot n)$ wall and takes 15+ min at n=6. Dense build is $\mathcal O(d^4 n) + \mathcal O(d^6)$ — 190 s at n=6. Wired into the driver as `bohr_gap_dense`.

## Model

Clean 2D quantum transverse-field Ising model on an `Lx × Ly` square lattice with periodic BCs:

$$
H \;=\; -J \sum_{\langle ij \rangle} Z_i Z_j \;-\; h \sum_i X_i, \qquad J = h = 1.0,\; \beta_{\rm phys} = 0.5 .
$$

At $\beta_{\rm phys} = 0.5$ with $h = J$ the operating point sits moderate/borderline-T relative to $T_c(h=J)$ — probably in the disordered phase. Goal here is **pipeline correctness on 2D**, not "see exponential gap closing"; that lives entirely on the cluster.

Construction:

- Bond term — `_construct_2d_heisenberg_base(Lx, Ly, [[Z,Z]], [-J]; periodic_x=true, periodic_y=true)`.
- Transverse field — `_construct_disordering_terms([[X]], [-h .* ones(n)], n)` (uniform per-site coefficients; the "disordering" builder is just the canonical single-site-Pauli summer).
- Disorder is **OFF by default** (`disorder_strength = 0.0`). Subtask B turns on `disorder_strength = 1e-3` with a `find_typical`-style W₂-spectral-median selector if and only if the clean Energy↔Bohr cross-check fails the 1e-7 gate.

## Sandbox Phase A — Cells Attempted

Sandbox budget: 1 h wall, in increasing-n order, stop-on-overrun (per-cell estimated extrapolation in `script.main` checks `elapsed + estimated_next` against the budget before starting each cell).

| Order | n | (Lx, Ly) | Class | d  | d² | OFT cache @ r_D=7 | Status |
|---|---|---|---|---|---|---|---|
| 1 | 4 | 2×2 | smoke (PBC double-bonds) | 16  | 256   | 512 KB | _measured below_ |
| 2 | 6 | 2×3 | quasi-1D ladder         | 64  | 4 096 | 8 MB   | _measured below_ |
| 3 | 8 | 2×4 | quasi-1D ladder         | 256 | 65 536| 128 MB | _measured below — only if budget remains_ |

## Pipeline

- CKG smooth-Metropolis, $s = 0.25$, $a = 0$, coherent term ON (canonical thesis convention 2026-05-12).
- EnergyDomain, $r_D = 7$, $w_{0,D} = 2(\lVert H \rVert_{\rm op} + 8\sigma) / 2^{r_D}$ (unified ω-range), $t_{0,D} = 2\pi / (2^{r_D} w_{0,D})$, $\sigma = 1/\beta_{\rm alg}$.
- `krylov_spectral_gap` matrix-free, `krylovdim = 40`, `tol = 1e-10`, default `include_coherent = true`.
- `predict_lindbladian_trajectory` for $\tau_{\rm mix}$ at $\varepsilon_{\rm target} = 10^{-3}$, $t_{\max} = 500$, $N_{\rm grid} = 81$.
- BohrDomain cross-check (`krylov_spectral_gap` only, no trajectory) at $n = 4$ and $n = 6$.

JULIA_NUM_THREADS = 8, OPENBLAS_NUM_THREADS = 1.

## Symmetry / Krylov-x0 contingency (subtask B) — _not triggered_

Clean 2D TFIM has a Z₂ spin-flip symmetry $P = \bigotimes_i X_i$ that commutes with both $H_{\rm bond}$ and $H_{\rm field}$. The CKG Lindbladian inherits a dual symmetry $\rho \mapsto P \rho P$, so $\mathcal L$ is block-diagonal in the Z₂-even / Z₂-odd sectors of vec($\rho$). The maximally mixed state $I/d$ and the Gibbs state $\sigma_\beta$ are both Z₂-even, so an unperturbed seed $\text{vec}(I/d)$ cannot reach any odd-sector mode through Arnoldi.

The qf-8fr `_krylov_default_x0` adds $\varepsilon \cdot \text{vec}(H_{\rm GUE,traceless})$ with $\varepsilon = 10^{-10}$ — generic Hermitian, so has $\mathcal O(\varepsilon)$ overlap with each sector, including Z₂-odd. At KrylovKit `tol = 1e-10` this sits **right at the floor**: it works if the gap eigenmode is Z₂-even (e.g. in the disordered phase), but is one of the most dangerous regimes if the gap mode is Z₂-odd (e.g. ordered-phase tunneling).

At $\beta_{\rm phys} = 0.5,\, h = J$ the operating point is in (or very near) the disordered phase, and the gap eigenmode is Z₂-even. The **Energy↔Bohr gate passed at 4.8e-14 (n=4) and 8.2e-14 (n=6) — no contingency triggered**. The qf-8fr 1e-10 GUE kick is sufficient here without any modification. The dense BohrDomain `eigvals` path is also unaffected by the symmetry blind-spot (it gets the full spectrum), and it agrees with EnergyDomain Krylov, which independently confirms the Krylov path is locking onto the right mode.

This does **not** rule out the Z₂-odd contingency biting at ordered-phase / critical-fan operating points (deeper $\beta_{\rm phys}$, $h < J$). The contingency machinery is wired into the driver and will fire automatically if it is needed.

Contingency menu wired in `scripts/scratch_2d_tfim_gap_sweep.jl::main`:

1. **`disorder_strength = 1e-3`** ([[Z],[Z,Z]] family, find-typical batch 256) — **wired in the script.** Not exercised in this run.
2. Bump `_krylov_default_x0` perturbation amplitude from $10^{-10}$ to $10^{-6}$ — not wired, would need a one-line patch to `src/krylov_eigsolve.jl`.
3. Multi-vector / block Krylov with one seed lying explicitly in the Z₂-odd sector (`vec(M_z · I/d)` with $M_z = \sum_i Z_i$) — significant refactor, deferred.

## Measured Feasibility

All wall numbers measured inside the Docker sandbox running on the laptop, single Julia process, `JULIA_NUM_THREADS = 8`, `OPENBLAS_NUM_THREADS = 1`. Run date: 2026-05-14.

### Spectral-gap sweep at $\beta_{\rm phys} = 0.5$, $h = J = 1.0$

| n | (Lx, Ly) | Domain | $\beta_{\rm alg}$ | $R$ (rescale) | $\lambda_{\rm gap,alg}$ | $\lambda_{\rm gap,phys}$ | $\tau_{\rm mix}(\varepsilon{=}10^{-3})$ | wall (s) | RSS (MB) | conv |
|---|---|---|---|---|---|---|---|---|---|---|
| 4 | 2×2 | EnergyDomain | 19.0  | 38.0 | 1.1816e-02 | 0.449 | 48.9  | 5.74   | 769.7  | 5 (KrylovKit) |
| 4 | 2×2 | BohrDomain   | 19.0  | 38.0 | 1.1816e-02 | 0.449 | —     | 1.50   | 797.3  | dense_eigvals |
| 6 | 2×3 | EnergyDomain | 25.2  | 50.4 | 5.6694e-03 | 0.286 | 101   | 9.13   | 844.5  | 4 (KrylovKit) |
| 6 | 2×3 | BohrDomain   | 25.2  | 50.4 | 5.6694e-03 | 0.286 | —     | 190.39 | 1277.5 | dense_eigvals |
| 8 | 2×4 | EnergyDomain | 37.8  | 75.6 | 3.7662e-03 | 0.285 | 186   | 768.15 | 1113.9 | 4 (KrylovKit) |

(`conv` is `info.converged` from `KrylovKit.eigsolve`; `dense_eigvals` indicates the path took `construct_lindbladian + eigvals(L_dense)` — LAPACK `geev`, no convergence iterations.)

EnergyDomain matvec count: 55 (n=4), 86 (n=6), 103 (n=8). All under `krylovdim = 40` × 3 restarts.

### Energy↔Bohr cross-check (1e-7 gate)

| n | (Lx, Ly) | $\lambda_{\rm Energy}$ | $\lambda_{\rm Bohr}$ | $\lvert \Delta \rvert / \lvert \lambda_B \rvert$ | result |
|---|---|---|---|---|---|
| 4 | 2×2 | 1.1815990074e-02 | 1.1815990074e-02 | **4.80e-14** | ✓ PASS |
| 6 | 2×3 | 5.6694266844e-03 | 5.6694266844e-03 | **8.20e-14** | ✓ PASS |

The 1e-7 gate is satisfied with 7 orders of magnitude headroom on both cells. `r_D = 7` is therefore massively over-converged for this fixture — consistent with [[Quadrature register recipe v2 qf-yt9]] which found r_D = 6 sufficient for $\varepsilon = 10^{-9}$ reference on canonical 1D fixtures, and the 2D TFIM rescaled spectrum lives in the same $[0, 0.45]$ band by construction.

### Cell wall breakdown

EnergyDomain wall per cell splits into Krylov-gap ($t_{\rm gap}$) plus `predict_lindbladian_trajectory` ($t_{\rm traj}$). Dense BohrDomain wall is the `construct_lindbladian + eigvals` build.

| n | $t_{\rm gap}^{\rm Energy}$ (s) | $t_{\rm traj}^{\rm Energy}$ (s) | $t_{\rm gap}^{\rm Bohr}$ (s) |
|---|---|---|---|
| 4 |   4.26 |   1.48 |   1.50 |
| 6 |   6.20 |   2.93 | 190.39 |
| 8 | 557.37 | 210.78 |  — (forbidden at n > 6) |

EnergyDomain wall scales sublinearly in $d^2$ between n=4 and n=6 (4.26 → 6.20 s on a $d^2$ jump of 16×) because the Arnoldi step count is similar (~50–86) and per-step BLAS GEMV throughput saturates only at larger $d$. The n=6 → n=8 jump (6.20 → 557 s) is the expected $d^4 \cdot n$ matvec growth dominating.

### gap_phys vs n — pipeline correctness

After undoing the rescaling, $\lambda_{\rm gap,phys} = \lambda_{\rm gap,alg} \cdot R$:

- n = 4 (2×2): 0.449 — outlier, PBC double-bonds make the lattice effectively a thick 1D ring.
- n = 6 (2×3): 0.286
- n = 8 (2×4): 0.285

The two bulk-like cells (n=6, n=8) collapse to $\lambda_{\rm gap,phys} \approx 0.285$. At $\beta_{\rm phys} = 0.5$ with $h = J = 1$ we expect the disordered phase, so an $\Omega(1)$ physical gap with no system-size decay is the predicted behaviour. The sweep is too short to claim a slope, but the qualitative match — flat in $n$ once past the PBC outlier — is the pipeline-correctness signal asked of subtask C.

### Memory usage

Peak RSS in the sweep was 1.28 GB (n=6 BohrDomain dense path; the $d^2 \times d^2 = 4096 \times 4096$ complex Lindbladian takes $\sim 256$ MB, on top of workspace + OFT cache + KrylovKit scratch). n=8 EnergyDomain peaked at 1.11 GB. Comfortable within the sandbox 1.8 GB envelope.

## Cluster Run Plan (subtask D)

The cluster run plan picks up at the cells the sandbox cannot reach. **Not executed in this issue.**

Common config:

- CKG smooth-Metropolis, $s = 0.25$, $a = 0$, coherent ON, $J = h = 1.0$.
- EnergyDomain, $r_D = 7$, $w_{0,D}$ from unified ω-range, `with_linear_combination = true`.
- `krylov_spectral_gap(krylovdim = 40, tol = 1e-10)`, default `include_coherent = true`.
- `predict_lindbladian_trajectory` for $\tau_{\rm mix}$ at $\varepsilon_{\rm target} = 10^{-3}$.

OFT cache memory is the binding constraint at our register sizes: at $r_D = 7$ the cache stores `N_ω = 2^r_D = 128` complex `d × d` matrices per jump = $128 \cdot d^2 \cdot 16$ bytes per jump, times $3n$ jumps. The "OFT cache" column below is the dominant term; full RSS (workspace, eigenvectors, trajectory states) is typically $\sim 1.5\text{–}3\times$ that.

### Tier 1 — first bulk-2D pair

| n  | (Lx, Ly) | Class | d   | d²       | OFT cache | Total RSS est | Threads | Expected wall |
|----|---------|---|-----|----------|---|---|---|---|
|  9 | 3×3     | **bulk 2D (first)**  | 512 | 262 144 | 0.5 GB | 1.5–3 GB | 8 | ~30 min – 1 h (also feasible on laptop overnight) |
| 12 | 3×4     | bulk 2D rect         | 4096| 1.68e7  | 32 GB | 100–150 GB | 64 | 2–6 h |

`JULIA_NUM_THREADS = max-available`, `OPENBLAS_NUM_THREADS = 1`.

Sidecars and summary BSON follow the same naming as the sandbox sweep — drop them under `scripts/output/sweep_2d_tfim/cluster/`.

### Tier 2 — physics sweeps (after Tier 1 lands)

After the bulk-2D pair $(n=9, n=12)$ is in hand at the baseline operating point, map the dome and the quantum-critical fan:

- **$\beta_{\rm phys}$ sweep $\in \{0.5, 1.0, 2.0\}$** on $n = 9$ and $n = 12$.
  Maps ordered-phase gap-closing rate vs $T$. $\beta_{\rm phys} = 1.0$ should already sit inside the dome at $h = J$; $\beta_{\rm phys} = 2.0$ deep inside.
- **$h$ sweep $\in \{0, 0.5, 1.0, 2.0, 3.0\}$** at $\beta_{\rm phys} = 1.0$ on $n = 9$.
  Crosses the finite-$T$ transition line; gap should be $\Omega(1)$ outside the dome, exponential inside. $h = 0$ recovers the 2D classical Ising.
- **Quantum-critical fan**: $h \approx h_c \approx 3.044$, $\beta_{\rm phys} \in \{2, 4, 8\}$.
  Power-law gap closing with dynamical exponent $z = 1$, $\nu \approx 0.6301$ predicted (3D classical Ising universality).

Each sweep adds $\mathcal O(\text{a few})$ to $\mathcal O(\text{a dozen})$ cells; per-cell cost dominated by $n = 12$ at the long-$\tau$ end.

### Infeasible at $r_D = 7$ (informational)

| n  | (Lx, Ly) | d   | d²       | OFT cache | Note |
|----|---------|-----|----------|---|---|
| 15 | 3×5     | 32768 | 1.07e9 | ~2 TB | infeasible at $r_D = 7$ |
| 16 | 4×4     | 65536 | 4.29e9 | ~8 TB | infeasible at $r_D = 7$ |

A "$r_D = 6$, $w_{0,D}$ retuned" path could shave the cache by 2× but doubles the dissipator quadrature error; not pursued without a separate quadrature sweep on 2D fixtures.

## Open follow-ups (deferred from this issue)

- **Promotion of `build_tfim_2d_raw` from scratch into `src/hamiltonian.jl`** as a `find_typical_2d_tfim` peer of `find_typical_2d_heisenberg`. Wait until the cluster run validates the approach and we know whether the optional disorder is ever needed in production.
- **DLL comparison on 2D** — separate issue (probably blocked on multi-rank DLL on a 2D fixture).
- **TrotterDomain + GQSP channel cross-check on 2D** — separate issue.
- **Ordered-phase operating point** (e.g. $h = 0.5$, $\beta_{\rm phys} = 2.0$ on $n = 9$) — only on the cluster, gates "is the symmetry contingency actually needed?"; if the cluster Energy↔Bohr gate also fails there, escalate to contingency #2 (GUE amplitude bump in `_krylov_default_x0`).

## Cross-references

- qf-8fr — 1D classical Ising / disordered Heisenberg gap-vs-n cross-check, the parent issue this builds on.
- qf-833 — this issue.
- `drafts/classical-ising-crosscheck.md` — qf-8fr findings note (algorithm-frame vs physical-frame gap, full Krylov-x0 rationale).
- `.claude-memory/krylov_x0_symmetric_bug_qf_8fr.md` — Krylov-x0 fix details + theorems backing the Ω(1) physical gap.
