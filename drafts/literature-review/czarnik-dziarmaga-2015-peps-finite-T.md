---
paper: "Czarnik and Dziarmaga (2015)"
title: "Projected entangled pair states at finite temperature: Imaginary time evolution with ancillas"
arxiv: "1209.0454"
year: 2015
venue: "PRB"
pdf: "supplementary-informations/classical-review/czarnik-dziarmaga-2015-peps-finite-T.pdf"

temperature: [any-finite-T, near-critical]
commutativity: [non-commuting, stoquastic]
locality: [geometric-local, 2D-square-lattice]
particle-statistics: [spin]
hamiltonian-models: [transverse-field-Ising, 2D-quantum-Ising]

paradigm: [classical-tensor-network]
quantum-or-classical: [classical]

result-type: [algorithm-heuristic, benchmark]
key-scaling: "Per-time-step cost $O(M^3 D^4)$ for $M\\ge D^4$ (else $O(M^2 D^8)$); bond dimension $D$ heuristic, environment dimension $M$ grows as criticality is approached; Trotter step $d\\beta\\ge 10^{-6}\\beta_c$ near $\\beta_c$"

related: ["molnar-schuch-verstraete-cirac-2015-pepo-gibbs", "white-2009-metts", "verstraete-garcia-ripoll-cirac-2004-mpdo", "kuwahara-alhambra-anshu-2021-thermal-area-law-1d", "hastings-2006-solving-gapped-locally", "feiguin-white-2005-finite-T-dmrg", "bravyi-gosset-2017-quantum-ferromagnets"]
---

# Projected Entangled Pair States at Finite Temperature: Imaginary Time Evolution with Ancillas — Czarnik, Cincio & Dziarmaga 2015

**One-sentence takeaway**: PEPS-with-ancilla purification evolved in imaginary time gives the workhorse classical method for thermal states of *2D quantum* lattice systems, demonstrated on the 2D quantum Ising model in transverse field with a full $h$–$\beta$ phase diagram at minimal bond dimension $D=2$.

## Setting

- **Hamiltonian**: 2D quantum Ising model on the infinite square lattice,
  $$ \mathcal{H} = -\sum_{\langle s,s'\rangle} Z_s Z_{s'} - h \sum_s X_s \equiv \mathcal{H}_{ZZ} + \mathcal{H}_X $$
  Pauli operators, ferromagnetic nearest-neighbour coupling, transverse field $h$. Classical critical point $\beta_c = -\tfrac{1}{2}\ln(\sqrt{2}-1)\approx 0.441$ (at $h=0$); quantum critical point $h_c\approx 3.04$ (at $\beta^{-1}=0$).
- **Ansatz**: translationally-invariant PEPS with ancilla, one tensor $A^{ia}_{trbl}(\beta)$ per site (spin index $i$, ancilla index $a$, four bond indices $t,r,b,l$ each of dimension $D$). Purification: $\rho(\beta)\propto \mathrm{Tr}_{\mathrm{anc}}|\psi(\beta)\rangle\langle\psi(\beta)|$ with $|\psi(\beta)\rangle = U(\beta)|\psi(0)\rangle = e^{-\tfrac{1}{2}\beta\mathcal{H}}|\psi(0)\rangle$, and $|\psi(0)\rangle = \prod_s \tfrac{1}{\sqrt{S}}\sum_i|i_s,a_s=i_s\rangle$ a product of maximally-entangled spin–ancilla pairs giving $\rho(0)\propto \mathbf{1}$.
- **Goal**: classical, polynomial-in-$D$ algorithm for thermal expectation values $\langle\mathcal{O}\rangle_{\beta,h}$ across the entire 2D phase diagram, including across the quantum/thermal critical lines.
- **Temperature regime**: any $\beta>0$ in principle; the practical bottleneck is the environment-tensor dimension $M$, which must grow as criticality is approached.

## Main Results

- **Algorithm**: imaginary-time evolution of the PEPS-with-ancilla via a second-order Suzuki–Trotter split $U(d\beta) = U_X(d\beta/2)U_{ZZ}(d\beta)U_X(d\beta/2) + \mathcal{O}(d\beta^3)$, with bond-dimension truncation after each $U_{ZZ}$ step performed by a corner-matrix-renormalization-group (CTMRG) variant that constructs an *effective* environment (tensors $C, V, H$) of dimension $M$.
- **Per-step cost** (verbatim, p. 2): "polynomial in both $D$ and $M$. It is dominated by the calculation of $V'$ in Fig. 3 that scales like $M^3 D^4$ when $M\ge D^4$ or $M^2 D^8$ otherwise."
- **Classical limit ($h=0$)**: $D=2$ is exact at every $\beta$ (Eq. 9 doubles the bond dimension exactly from $D=1\to 2$); benchmark against Onsager's exact solution shows agreement across a wide $\beta$ range; convergence to the exact $\langle ZZ\rangle$, $\langle Z\rangle$ at $\beta\approx\beta_c$ requires $M\in\{8,16,32,56\}$ (Fig. 5). At $D=2$, $M=32$ the correlator $C_R = \langle Z_s Z_{s'}\rangle$ tracks the exact $C_R\sim R^{-1/4}$ over $R\in[1,1000]$ down to $\varepsilon = (\beta_c-\beta)/\beta_c = 10^{-4}$ (Fig. 6).
- **Quantum case ($h>0$)**: With a small longitudinal symmetry-breaking field $h_Z\in[10^{-3},10^{-10}]$ to smooth out the non-analyticity at the critical point, the magnetization $\langle Z\rangle$ at $h = \tfrac{2}{3}h_c$ converges to a non-analytic critical curve as $h_Z\to 0$; results for $D\in\{2,4,6\}$ with required $M\le 16$ are quantitatively close, supporting that "even the minimal non-trivial $D=2$ ($D=S$ in general) … converge to a non-analytic critical curve" (p. 5). A full $h$–$\beta$ phase diagram is produced at $D=2$, $M=24$, $h_Z=10^{-10}$, $d\beta=10^{-2}\beta_c$ (Fig. 8).
- **Failure mode at criticality**: for finite $M$ the renormalization is bound to accumulate unrecoverable errors crossing $\beta_c$; the explicit symmetry-breaking field $h_Z$ rounds the non-analyticity into a smooth crossover that finite-$M$ PEPS can integrate through.

## Method

A single-tensor PEPS purification on the doubled (spin+ancilla) Hilbert space is initialized at $\beta=0$ (maximally mixed) and propagated to inverse temperature $\beta$ by repeated application of a second-order Suzuki–Trotter step: each two-body $U_{ZZ}(d\beta)$ doubles the bond dimension from $D$ to $2D$, then a CTMRG-style truncation reduces it back to $D$ via an isometry $W$ chosen to be least distortive in the $\langle\psi_B|\psi_B\rangle$ metric, using converged "environment tensors" $C, V, H$ of effective dimension $M$ that approximate the rest-of-the-lattice contraction. Practical fixes: recycle the previous step's environment with $<1\%$ added noise to prevent it getting trapped in subspaces $M_{\rm eff}<M$, and add a tiny longitudinal $\delta\mathcal{H}=-h_Z\sum_s Z_s$ to round the critical non-analyticity into a finite-$M$-tractable crossover.

## Quantum vs Classical

- **Classical baseline this *defines***: this is the canonical *classical* finite-$T$ method for *2D quantum* lattice systems in the Tier-1 corpus's "2D / higher-$d$ quantum" cross-cutting cell. The corpus already labels it as such ("Czarnik–Dziarmaga (PEPS imaginary-time)"). Any quantum Gibbs-sampling claim of advantage on a 2D quantum lattice model has to beat this method's *empirical* cost — which is where the contest actually lives, since the rigorous PEPO bound (Molnar et al.) is much looser.
- **What this beats**: the only prior 2D options were (i) the higher-order SVD tensor product / PEPO ansatz of Xie et al. and Orús (rigorous bond-dimension bound but conservative), (ii) brute-force exact diagonalization on small finite tori, and (iii) QMC where applicable — but the 2D *quantum* TFIM is non-stoquastic only mildly (it *is* stoquastic in the standard basis once one chooses signs), so QMC is competitive here, not beaten outright. PEPS wins on access to the *thermodynamic limit* directly via the translationally-invariant ansatz.
- **What it doesn't do**: provide *rigorous* polynomial-time guarantees. The $D$ vs $\beta$ vs correlation-length cost trade-off is heuristic; no theorem says $D = D(\beta,\varepsilon)$ suffices. That rigorous backing is supplied separately by Molnar-Schuch-Verstraete-Cirac (2015), which proves $D = \exp(O(\beta))$ suffices for *any* local Hamiltonian in *any* dimension via a different (PEPO) construction, and by Hastings (2006) / Kuwahara-Alhambra-Anshu (2021) for the 1D analogue.
- **Where quantum *might* still win**: (i) frustrated 2D quantum models where the Trotter error and the $M$ blow-up at criticality compound — e.g. 2D Heisenberg antiferromagnet with frustration, or non-stoquastic 2D quantum lattice models, where the heuristic environment tensor renormalization may not converge. (ii) 3D quantum lattices, where the analogous ansatz becomes vastly more expensive. (iii) Deep low-$T$ regimes where $D$ must scale as $\exp(O(\beta))$ to capture long-range entanglement.
- **Where quantum *won't* win**: 2D stoquastic / quantum-Ising-like models in the disordered or weakly-correlated phases — empirically, $D=2$, $M\le 24$ already gives the full phase diagram. A quantum Gibbs sampler on this regime cell offers no asymptotic improvement.

## Implications for Quantum Advantage

- **Regime cell**: this paper anchors the **2D quantum × any-T × geometric-local × non-frustrated** cell with a heuristic but practically very effective classical method. Given that QMC also handles 2D stoquastic models efficiently (no sign problem on the standard basis here), this cell is *doubly* covered classically.
- **Reading for the thesis Review chapter**: the punch-line non-stoquastic × low-T cell remains the strongest quantum-advantage candidate precisely because *both* QMC *and* PEPS imaginary-time evolution degrade there — QMC from the sign problem, PEPS from the environment-tensor blow-up at low $T$ in models with high entanglement (which is also where rigorous PEPO bounds force $D = \exp(O(\beta))$).
- **Empirical bar to beat**: any quantum 2D Gibbs-sampler claim must improve over $O(M^3 D^4)$ per Trotter step at the $D, M$ values needed for the target observable's relative error in the 2D quantum Ising model — a low bar in absolute terms ($D=2$, $M\sim 30$) but a high bar in *scaling guarantees*, since this paper makes none.

## Open Questions / Limitations

- **Rigorous bond-dimension bound for PEPS-with-ancilla**: the $D = D(\beta,\varepsilon)$ scaling is left empirical. Molnar-Schuch-Verstraete-Cirac (2015) supply a rigorous $\exp(O(\beta))$ PEPO bond bound but their construction is different. Whether the imaginary-time-evolution route admits a comparable theorem is open. [CHECK]
- **Convergence at criticality**: explicitly noted that for $M$ finite the algorithm "is bound to accumulate unrecoverable errors near the critical point." The symmetry-breaking field $h_Z$ is a workaround, not a fix. There is no quantitative bound on how $M$ scales with $h_Z$ or with $(\beta_c-\beta)/\beta_c$.
- **Beyond Ising**: only the 2D quantum transverse-field Ising model is benchmarked. The paper does not test e.g. quantum Heisenberg antiferromagnet on the square lattice, where frustration could break the simple environment-tensor picture.
- **Real-time evolution / finite-$T$ dynamics**: the conclusion notes the algorithm can be modified for real-time evolution and (with fermionic swap gates) for fermionic 2D systems, but those extensions are deferred (cited as "P. Czarnik et al., in preparation"). Status today is independently developed — this paper itself does not address them.
- **Comparison vs PEPO ansatz** (Czarnik-Dziarmaga 2015b, "Variational tensor network renormalization in imaginary time", arXiv:1503.01077): the same authors followed up with a variational PEPO. Side-by-side accuracy/cost comparison vs the present method is not given here.

## Connections

- **Direct 1D analogue**: Verstraete, García-Ripoll, Cirac (2004) introduced the matrix product density operator (MPDO) purification + imaginary-time TDVP framework that this paper lifts to 2D. The 1D version is rigorous (Hastings 2006; Kuwahara-Alhambra-Anshu 2021); the 2D version is heuristic.
- **Direct rigorous counterpart**: Molnar-Schuch-Verstraete-Cirac (2015) prove $D = \exp(O(\beta))$ PEPO bond dimension suffices for any local Hamiltonian in any dimension. This paper is the *algorithmic* companion — what one actually computes when the rigorous bound is too pessimistic.
- **Sibling tensor-network finite-$T$ method (1D)**: White (2009) METTS samples MPS-friendly typical thermal states by alternating imaginary-time evolution and projective measurement, achieving $10^3$–$10^{10}$ speedup over MPDO purification in 1D — the lower-entanglement-per-sample idea is orthogonal to and could in principle be married with the 2D PEPS purification here, but 2D METTS is not in this paper.
- **Earlier 2D PEPO work**: Verstraete-García-Ripoll-Cirac (Ref. [10] in the paper) proposed a 2D PEPO using higher-order SVD; the present paper's contribution is the *imaginary-time-evolution* route via a single PEPS-with-ancilla, which is the natural lift of the MPS-with-ancilla 1D method.
- **CTMRG / corner matrix renormalization**: the environment-tensor construction adapts Baxter–Nishino corner matrix renormalization to 2D quantum; this is the same tool used in independent PEPS ground-state codes (Orús-Vidal).
- **vs sign-free QMC**: Syljuåsen-Sandvik directed-loop SSE and Prokof'ev-Svistunov worm algorithms cover the same 2D-quantum-stoquastic cell with rigorous mixing only in narrow regimes but excellent empirical performance. PEPS-imaginary-time is the alternative when one wants direct access to the *thermodynamic limit* without finite-size extrapolation, or when frustration/sign-problem makes QMC fail.
- **vs Bravyi-Gosset (2017)**: Bravyi-Gosset gives a *rigorous* polynomial-time classical FPRAS for stoquastic transverse-field XY/Heisenberg ferromagnets on *any* graph, including the 2D TFIM — so the cell that this paper covers heuristically is also covered with rigorous guarantees by Bravyi-Gosset. PEPS-imaginary-time wins on practical cost; Bravyi-Gosset wins on theoretical guarantee.
- **Quantum-side counterparts**: Chen et al. 2025 noncommutative KMS Gibbs sampler and Ding et al. 2024 KMS samplers are the quantum-algorithm pieces this regime-cell entry exists to be compared against in the thesis Review chapter.
