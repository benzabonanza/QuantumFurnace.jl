---
paper: "White (2009)"
title: "Minimally entangled typical thermal states"
arxiv: "0902.4475"
year: 2009
venue: "PRL 102, 190601"
pdf: "n/a"

temperature: [all]
commutativity: [non-commuting]
locality: [1D, geometric-local]
particle-statistics: [spin]
hamiltonian-models: [xxz, heisenberg-1d, transverse-ising-1d]

paradigm: [classical-tensor-network, mps, monte-carlo-hybrid]
quantum-or-classical: [classical]

result-type: [algorithm, sample-complexity]
key-scaling: "Per-sample cost $O(N m_0^3 \\beta N_\\tau)$ vs ancilla DMRG $O(N m_0^6 \\beta)$; speedup factor $m_0^3/N_\\tau \\sim 10^3$–$10^{10}$ at $m_0 \\in [50,5000]$, $N_\\tau \\in [10,100]$"

related: ["czarnik-dziarmaga-2015-peps-finite-T", "kuwahara-alhambra-anshu-2021-thermal-area-law", "molnar-schuch-verstraete-cirac-2015-pepo-gibbs"]
---

# Minimally Entangled Typical Thermal States — White 2009

**One-sentence takeaway**: A Markov chain on classical product states (CPS), with each step an imaginary-time evolution $e^{-\beta H/2}$ on an MPS followed by a projective measurement collapsing back to a CPS, samples thermal expectation values $\langle A\rangle_\beta$ in 1D at $10^3$–$10^{10}$ times the speed of ancilla/purification DMRG by exploiting that *typical* thermal states are far less entangled than a purification of $\rho_\beta$.

## Setting

- **What is sampled**: thermal expectation values $\langle A\rangle_\beta = \mathrm{Tr}(\rho_\beta A)$ with $\rho_\beta = e^{-\beta H}/Z$, written in the Monte-Carlo form
  $$ \langle A\rangle_\beta = \frac{1}{Z}\sum_i P(i)\,\langle\phi(i)|A|\phi(i)\rangle,\qquad P(i)=\langle i|e^{-\beta H}|i\rangle,\quad |\phi(i)\rangle = \frac{e^{-\beta H/2}|i\rangle}{\sqrt{P(i)}}, $$
  where $\{|i\rangle\}$ is any orthonormal product basis (CPS).
- **Hamiltonian class**: 1D local Hamiltonians (benchmarked on a 100-site spin-$\tfrac12$ Heisenberg chain; method applies broadly to 1D MPS-amenable models — XXZ, transverse-field Ising, Hubbard).
- **Access model**: classical algorithm; only requires that $e^{-\beta H/2}|i\rangle$ admits a faithful MPS representation at modest bond dimension $m$.
- **Temperature regime**: any $\beta>0$; in practice efficient as long as the *typical* states $|\phi(i)\rangle$ obey an area law with manageable bond dimension $m_0$. Low-$T$ pushes $m_0$ towards the *ground-state* MPS bond dimension (not its square).

## Main Results

- **METTS sampler (algorithm)**: alternate (i) imaginary-time evolution $|\phi(i)\rangle = e^{-\beta H/2}|i\rangle/\sqrt{P(i)}$ via tDMRG; (ii) sequential single-site projective measurement of $|\phi(i)\rangle$ in (random) local bases, yielding a new CPS $|i'\rangle$ sampled with the correct conditional probability. The chain has stationary measure $P(i)/Z$, so empirical averages converge to $\langle A\rangle_\beta$.
- **Per-sample cost**: $O(N m_0^3 \beta N_\tau)$ for METTS vs $O(N m_0^6 \beta)$ for ancilla/purification DMRG (Feiguin-White 2005), where $m_0$ is the bond dimension needed to represent a single METTS, $N$ the chain length, $N_\tau$ the number of samples per estimate.
- **Empirical speedup**: factor $m_0^3/N_\tau \approx 10^3$–$10^{10}$ over ancilla DMRG, taking $m_0\in[50,5000]$ and $N_\tau\in[10,100]$ from typical 1D benchmarks. Energy converges to ancilla-DMRG accuracy after only $5$–$10$ thermal steps from a random initial CPS (Fig. 1).
- **Bond-dimension separation**: METTS bond dimension $m$ starts at $1$ for $\beta\to 0$ and saturates to $m_0\sim$ ground-state $m$ as $\beta\to\infty$; ancilla DMRG saturates to $m_0^2$ (purification doubles entanglement entropy). This is the structural source of the speedup.
- **Variance / basis freedom**: measurement axes can be all the same, all different, fixed or random; the paper recommends *random* local axes to break ergodicity traps and reduce variance. Variance is controlled by the standard MCMC empirical variance of $\langle\phi(i)|A|\phi(i)\rangle$ over the chain.

## Method

The thermal trace is rewritten as a sum over CPS indices weighted by $P(i)/Z$. Two facts make this a classical algorithm: (a) typical $|\phi(i)\rangle$ are MPS-friendly because imaginary-time evolution from a product state cannot generate more entanglement than the *typical* (not maximal) thermal entropy density, in line with 1D thermal area laws; (b) the projective measurement on $|\phi(i)\rangle$ exactly samples a new CPS $|i'\rangle$ from the conditional distribution making the chain detailed-balanced w.r.t. $P(i)/Z$. Random measurement bases reduce autocorrelation and avoid getting stuck in symmetry sectors. Compared to ancilla DMRG, one trades a deterministic but doubly-entangled purification for a Monte-Carlo over singly-entangled MPS — an exponential saving in $m$ at the cost of $N_\tau$ samples.

## Quantum vs Classical

- **Baseline (classical)**: ancilla / purification finite-$T$ DMRG, Verstraete-García-Ripoll-Cirac (2004) and Feiguin-White (2005), at cost $O(N m_0^6 \beta)$ per evaluation. METTS is a strictly better classical algorithm in the same regime cell.
- **Baseline (quantum)**: Chen-Kastoryano-Brandão-Gilyén (2023) and Ding-Chen-Lin (2024) KMS detailed-balance Lindbladians on 1D local Hamiltonians, with mixing time bounds (under MLSI/spectral-gap assumptions) of polynomial-in-$n,\beta,1/\varepsilon$ form. For 1D gapped local Hamiltonians, recent work (Kastoryano-Brandão, Bardet et al., Capel et al.) gives quantum mixing times $\tilde O(\mathrm{poly}(n))$ at fixed $T$.
- **Gap**: classical METTS achieves polynomial-in-$N$ cost on 1D thermal expectation problems wherever the typical-state bond dimension is polynomial — i.e. the entire 1D area-law regime (Hastings 2006; Kuwahara-Alhambra-Anshu 2021). The corresponding quantum Gibbs sampler has the same polynomial scaling (often with worse constants and dimension-dependent prefactors). **No quantum advantage in $n$** in this cell. The gap, if any, is in $\beta$ and $\varepsilon$ prefactors and is empirically *favourable to classical*.
- **Source of the difference**: METTS exploits a structural fact — typical pure thermal states obey a tighter area law than the full mixed-state purification. The quantum sampler does not exploit this structure; it works at the level of the full $\rho_\beta$. Thus the *physical* reason classical wins here is the same reason 1D thermal physics is tractable: bounded entanglement of typical states.
- **Caveat**: METTS is a heuristic — there is no theorem giving $m_0(\beta,\varepsilon)$ or autocorrelation-time bounds. Rigorous backing comes separately from the 1D thermal area law (Hastings; Kuwahara-Alhambra-Anshu). The quantum samplers, by contrast, often come with proofs (modulo MLSI/gap assumptions). So *rigorously* the comparison is heuristic-classical vs proven-quantum, but *practically* METTS wins by orders of magnitude on every 1D benchmark.

## Implications for Quantum Advantage

- **Regime cell**: **all-T × non-commuting × 1D-geometric-local × spin**. METTS, together with the rigorous PEPO/MPS bounds of Hastings (2006), Molnar et al. (2015), and Kuwahara-Alhambra-Anshu (2021), make this the strongest cell *against* quantum advantage in the corpus map.
- **What this changes**: closes the upper-bound side of the 1D quantum thermal-sampling cell with a fast, general-purpose, empirically robust classical algorithm. Any quantum 1D Gibbs-sampling claim must beat $O(N m_0^3 \beta N_\tau)$ at the relevant $m_0, N_\tau$ — and given $m_0$ is polynomial in $N$ for gapped 1D systems, this is a high bar.
- **Promising or not**: 1D quantum Gibbs sampling is *unpromising* as a venue for first quantum advantage. The remaining quantum-leverage candidates within 1D are (i) critical 1D systems where $m_0$ grows polynomially with $N$ (ancilla then grows as $m_0^2$, but METTS still scales better), (ii) very low $T$ where sampling autocorrelation $N_\tau$ explodes — but these are also where the assumed quantum mixing-time bounds degrade. The cell is doubly hard for quantum.

## Open Questions / Limitations

- **No autocorrelation bound**: empirical claim of $5$–$10$ thermal steps to convergence is not backed by a mixing-time theorem. Pathological 1D models (e.g. near criticality, or at very low $T$ where the chain tunnels between symmetry-broken sectors) could in principle have long autocorrelation.
- **Critical regime**: the paper does not address criticality directly; bond dimension $m_0$ scales polynomially with $N$ (as $\sim N^{c/6}$ via CFT entanglement-entropy bounds), making METTS only polynomially-not-exponentially efficient there. Still beats ancilla, but no longer "trivial".
- **Failure of typicality interpretation**: White himself notes METTS "are very likely not 'typical'" in the strict random-state sense; they are *minimally* entangled, which is what makes them MPS-friendly but means they don't reflect the entanglement of generic thermal pure states.
- **2D and higher**: the original paper is 1D-only. 2D METTS with PEPS (Chen-Stoudenmire 2024, arXiv:2310.08533) is a separate development; the speedup over 2D PEPS purification (Czarnik-Dziarmaga 2015) is much less dramatic than in 1D because PEPS contraction is itself heuristic.
- **Real-time dynamics**: the paper notes that real-time evolution from METTS suffers entanglement growth like any other initial state — METTS does not help with finite-$T$ dynamics, only thermal averages.

## Connections

- **Direct rigorous backing for the 1D cell**: Kuwahara-Alhambra-Anshu (2021) prove a 1D thermal area law that justifies the polynomial $m_0(\beta)$ that METTS empirically observes. Together they make the 1D classical case airtight.
- **Direct rigorous companion in higher dimensions**: Molnar-Schuch-Verstraete-Cirac (2015) prove $D=\exp(O(\beta))$ PEPO bond-dimension suffices for $\rho_\beta$ in any dimension; 2D METTS (PEPS variant) is the algorithmic counterpart. Worse rigorous bound, better empirical performance.
- **2D-quantum sibling**: Czarnik-Dziarmaga (2015) PEPS-with-ancilla is the 2D analogue of *ancilla DMRG*; 2D METTS with PEPS sampling exists (arXiv:2310.08533) as the analogue of the present paper.
- **Direct improvement ancestor**: Verstraete-García-Ripoll-Cirac (2004) MPDO purification + Feiguin-White (2005) ancilla DMRG. METTS replaces both with a Monte-Carlo over MPS at the cost of variance.
- **Hybrid descendants**: Hybrid purification + sampling (Phys. Rev. B 101, 195119) and METTS with auxiliary MPS bases (arXiv:1910.03329) — refinements that further reduce sample variance.
- **Quantum-side counterparts (the comparison this enables)**: Chen et al. (2023, 2025) and Ding et al. (2024, 2025) KMS-detailed-balance Lindbladians. These are the quantum 1D Gibbs samplers the thesis Review chapter must measure against this paper.
- **Bravyi-Gosset (2017)**: rigorous classical FPRAS for stoquastic ferromagnetic XY/Heisenberg in any dimension covers the 1D Heisenberg ferromagnet rigorously; METTS covers the broader 1D non-stoquastic / antiferromagnetic / Hubbard regime heuristically but at much greater speed in practice.
