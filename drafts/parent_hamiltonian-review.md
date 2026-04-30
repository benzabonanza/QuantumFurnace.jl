# Draft Review: Parent Hamiltonian Subsection

**Reviewed**: `drafts/parent_hamiltonian.md`
**Date**: 2026-04-09
**Overall assessment**: The draft is well-structured and covers the right material at the right level. However, it contains one **critical misattribution** of the area law result, a **wrong citation location** for the mixing time bound, several notation inconsistencies with the thesis, and a completely wrong author list for the SDP hierarchy paper. Needs one round of targeted fixes before integration.

---

## Critical Issues

### Issue 1: Area law misattributed to wrong authors

- **Location**: Line 44, final paragraph: "Gilyen, Rouze and Stilck Franca \cite{gilyen2024quantum} proved an area law..."
- **Problem**: The paper at arXiv:2212.10061 is titled "Area law for steady states of detailed-balance local Lindbladians" and is by **Raz Firanko, Moshe Goldstein, and Itai Arad** -- not by Gilyen, Rouze, and Stilck Franca. The draft's writing notes already flag this as [CHECK], and the check reveals the attribution is wrong.
- **Why it matters**: Attributing a result to the wrong authors in a thesis is a serious scholarly error.
- **Suggested fix**: Replace `Gilyen, Rouze and Stilck~Franca \cite{gilyen2024quantum}` with `Firanko, Goldstein and Arad \cite{firanko2023area}` (using arXiv:2212.10061). The bib key `gilyen2024quantum` must be replaced with the correct entry. Note: the CKG paper (page 19) cites [Has07] (Hastings 2007) for the area law connection, which is a *different* area law (for ground states of 1D gapped Hamiltonians). The Firanko-Goldstein-Arad paper specifically proves an area law for *steady states* of detailed-balance Lindbladians, which is what the draft discusses. Verify the exact statement: the bound may differ from the stated $\tilde{\mathcal{O}}(\log^3(d)/\lambda)$ -- check the actual theorem in the paper.

### Issue 2: The $\tilde{\mathcal{O}}(\log^3(d)/\lambda)$ mutual information bound needs verification

- **Location**: Line 44: "the mutual information across any cut satisfies $I(L:R)_\sigma = \tilde{\mathcal{O}}(\log^3(d) / \lambda)$"
- **Problem**: This specific functional form needs to be verified against the actual theorem in Firanko-Goldstein-Arad (arXiv:2212.10061). The draft's writing notes already flag this as [CHECK]. The CKG paper's Appendix C discusses the area law but references Hastings [Has07] for ground states, not the steady-state version. The exact bound may have a different form or different parameters.
- **Why it matters**: Stating a wrong bound with a citation that doesn't match is a factual error.
- **Suggested fix**: Read arXiv:2212.10061 to extract the precise theorem statement, then update both the bound and the citation.

### Issue 3: SDP hierarchy paper misattributed to "Baez et al."

- **Location**: Line 35: "The recent SDP hierarchy of Baez et al.\ \cite{baez2024sdp}"
- **Problem**: The paper at arXiv:2411.03680 ("A Hierarchy of Spectral Gap Certificates for Frustration-Free Spin Systems") is by **Kshiti Sneh Rai, Ilya Kull, Patrick Emonts, Jordi Tura, Norbert Schuch, and Flavio Baccari** -- not "Baez et al." The draft's writing notes already flag this as [CHECK] and note the title confusion.
- **Why it matters**: Wrong author attribution in a thesis.
- **Suggested fix**: Replace `Baez et al.\ \cite{baez2024sdp}` with `Rai et al.\ \cite{rai2024hierarchy}` and create the correct bib entry for arXiv:2411.03680.

### Issue 4: "Bergamaschi and Chen" citation needs updating

- **Location**: Line 31: "Bergamaschi and Chen \cite{bergamaschi2024quantum}"
- **Problem**: The thesis already cites [BCL24] for "Bergamaschi, Chi-Fang Chen, and Yunchao Liu" (arXiv:2404.14639, "Quantum computational advantage with constant-temperature Gibbs sampling"). But the draft seems to refer to a *different* paper about rapid mixing and spectral gap stability. Looking at the CKG paper's reference list, [BC25] = "Thiago Bergamaschi and Chi-Fang Chen, Quantum spin chains thermalize at all temperatures, arXiv:2510.08533, 2025." The draft should cite this latter paper (about gap stability via frustration-freeness), not the quantum advantage paper.
- **Why it matters**: Citing the wrong paper by the same authors.
- **Suggested fix**: Use the correct reference -- arXiv:2510.08533 (Bergamaschi and Chen 2025, "Quantum spin chains thermalize at all temperatures") for the spectral gap stability result. Note the existing thesis bib key [BCL24] is for a different paper by partially overlapping authors.

---

## Notation Issues

### Issue 5: Discriminant notation mismatch

- **Location**: Lines 9-11, the discriminant formula
- **Thesis convention**: Equation (3.42) defines the discriminant as $\mathcal{D}(\sigma, \mathcal{L})$ using $\sigma$ for a general stationary state. The surrounding text on thesis p.34 introduces it generically. When specializing to Gibbs sampling, the thesis uses $\rho_\beta$ for the Gibbs state (eq. 3.32).
- **Draft uses**: $\mathcal{D}(\rho_\beta, \mathcal{L})$ on line 10, which is consistent with the thesis when specialized. However, the draft then switches to calling the result $\mathcal{H}_\beta$, which is not defined in the thesis preliminaries. The thesis text on p.34 calls the Hermitian part $\mathcal{H}(\rho, \mathcal{L})$ (not $\mathcal{H}_\beta$).
- **Fix**: The symbol $\mathcal{H}_\beta$ is fine as a shorthand for the specialization to Gibbs states, but it should be explicitly defined: "We write $\mathcal{H}_\beta := \mathcal{D}(\rho_\beta, \mathcal{L})$ for the discriminant specialized to the CKG Lindbladian." The CKG paper uses the same symbol $\mathcal{H}_\beta$ (p.6), so this is natural, but the thesis has not yet introduced it. Add a sentence establishing the shorthand.

### Issue 6: Purified Gibbs state notation

- **Location**: Line 13: $|\sqrt{\rho_\beta}\rangle$
- **Thesis convention**: The thesis does not yet introduce this notation (the Parent Hamiltonian section was empty). However, the CKG paper uses $|\sqrt{\rho_\beta}\rangle$ (eq. 1.9 on p.5).
- **Draft uses**: $|\sqrt{\rho_\beta}\rangle$ with the sum using $\ee^{-\beta E_i/2}$ and $|\psi_i^*\rangle$.
- **Fix**: This notation is consistent with the CKG paper. However, the thesis uses the custom macro `\ee` for Euler's number -- verify the draft's LaTeX will compile correctly with the thesis macros. Also: the draft writes $\ee^{-\beta E_i/2}$ in the sum but $\ee^{-\beta H}$ in the normalization, which is correct since $\tr[\ee^{-\beta H}] = \sum_i \ee^{-\beta E_i}$. No error here, but the notation $|\psi_i^*\rangle$ (complex conjugate) should be briefly explained since it is used without comment.

### Issue 7: $\lambda_\text{gap}$ vs $\lambda$ for spectral gap

- **Location**: Lines 17-19, eq. for $t_\text{mix}$
- **Thesis convention**: The thesis defines the spectral gap as $\lambda(\mathcal{L}) := -\lambda_2(\mathcal{L})$ (eq. 3.52, p.38). The mixing time bound on p.39 (eq. 3.55) uses $\lambda(\mathcal{L})$ in the denominator.
- **Draft uses**: $\lambda_\text{gap}(\mathcal{H}_\beta)$ -- a different symbol that refers to the gap of the *parent Hamiltonian*, not the Lindbladian.
- **Fix**: This is actually correct and intentional (the point of the section is that the Lindbladian gap equals the parent Hamiltonian gap), but the relationship should be made explicit. Add a sentence: "Since the discriminant transformation preserves the spectrum, $\lambda_\text{gap}(\mathcal{H}_\beta) = \lambda(\mathcal{L})$." This connects the notation to eq. (3.52).

### Issue 8: Doubled Hilbert space notation

- **Location**: Line 12: $\mathbb{C}^{2^n}\otimes\mathbb{C}^{2^n}$
- **Thesis convention**: The thesis uses $\mathcal{H}$ for the Hilbert space and $\mathcal{B}(\mathcal{H})$ for the space of bounded operators. It does not use $\mathbb{C}^{2^n}$ except implicitly.
- **Fix**: Minor. Could write "$\mathcal{H}\otimes\mathcal{H}$" for consistency, but $\mathbb{C}^{2^n}\otimes\mathbb{C}^{2^n}$ is fine for a qubit system. No change strictly needed.

---

## Missing Citations

- **Line 44** (Verstraete, Garcia-Ripoll, Cirac): The `verstraete2004matrix` citation is new and needs to be added to the bibliography. This is a well-known paper (PRL 93, 207204, 2004) and the attribution is correct.

- **Line 33** (Lin 2025): The `lin2025dissipative` citation (arXiv:2505.21308) is new. This paper is available in `supplementary-informations/`. The reference to "Section 3.1" for $\mathcal{O}(\log n)$-depth circuit constructions should be verified against the actual paper.

- **Lines 31-35**: Several new citations (`michalakis2013stability`, `anshu2016simple`, `gosset2016local`, `baez2024sdp` [should be `rai2024hierarchy`], `minganti2022arnoldi`) all need bib entries. The draft's writing notes correctly flag these.

- **Missing from draft**: The thesis already cites [RFA25] (Rouze, Stilck Franca, Alhambra 2025) and [SMBB25] (Smid et al. 2025) -- the draft uses `rouze2025efficient` and `vsmid2025polynomial` which will need to map to these existing bib keys.

---

## Consistency Issues

### The mixing time bound attribution

- **Location**: Line 16-17: "as shown in \cite[Appendix~C.4]{chen2023efficient}"
- **Issue**: The bound $\frac{\ln 2}{\lambda_\text{gap}(\mathcal{H}_\beta)} \leq t_\text{mix}(\mathcal{L}) \leq \frac{\ln(2\|\rho_\beta^{-1/2}\|)}{\lambda_\text{gap}(\mathcal{H}_\beta)}$ is stated in the CKG paper's Section C.4 (p.39), so the citation `[Appendix~C.4]{chen2023efficient}` is correct. However, the lower bound $\ln(2)/\lambda_\text{gap}$ and upper bound $\ln(2\|\rho^{-1/2}\|)/\lambda_\text{gap}$ were originally proven in the CKBG paper (Proposition II.2 on p.68 and Proposition E.5 on p.68-69). The CKG paper merely restates this result. For proper attribution, consider citing both, or at least the primary source.

### Frustration-freeness sign convention

- **Location**: Line 24-26
- **Issue**: The CKG paper's footnote 7 (p.6) notes: "Strictly speaking, originating from a Lindbladian, here the parent Hamiltonian is negative semi-definite, and the purified Gibbs state is the top-eigenstate. Introducing a global negative sign will make it the ground state." The draft writes $\mathcal{H}_\beta^a \geq 0$ on line 26, implying positive semi-definiteness, which corresponds to $-\mathcal{H}_\beta$ (the sign-flipped version). The CKG Proposition I.1 states $\mathcal{H}_\beta^a |\sqrt{\rho_\beta}\rangle = 0$ without specifying a sign convention for the individual terms. The draft should be explicit about whether it works with $\mathcal{H}_\beta$ or $-\mathcal{H}_\beta$. As written, the draft says the discriminant is non-positive (Lindbladian spectrum lies on non-positive real line, line 11) but then says $\mathcal{H}_\beta^a \geq 0$ (line 26), which requires flipping the overall sign. This should be made explicit.
- **Suggested fix**: Add a sentence: "Working with $-\mathcal{H}_\beta$ (equivalently, negating the Lindbladian spectrum to make it non-negative), we obtain a positive semidefinite parent Hamiltonian..." or explicitly define $\mathcal{H}_\beta := -\mathcal{D}(\rho_\beta, \mathcal{L})$ to flip the sign.

### Relationship to thesis eq. (3.55)

- The thesis already derives a mixing time upper bound in eq. (3.55): $t_\text{mix}(\varepsilon) \leq \frac{1}{\lambda(\mathcal{L})} \ln\left(\frac{1}{\varepsilon\sqrt{\lambda_\text{min}(\rho_\beta)}}\right)$. The draft's eq. on line 18 gives $t_\text{mix} \leq \frac{\ln(2\|\rho_\beta^{-1/2}\|)}{\lambda_\text{gap}(\mathcal{H}_\beta)}$. These should be the same bound (note $\|\rho_\beta^{-1/2}\| = \lambda_\text{min}(\rho_\beta)^{-1/2}$). The draft should reference eq. (3.55) and note that this is the same bound rewritten in terms of the parent Hamiltonian gap, to avoid the reader thinking it is a new result.

---

## Exposition Improvements

### Suggestion 1: Sign convention needs upfront clarification

- **Location**: Lines 11-12 and 24-26
- **Current**: The draft says the spectrum "lies on the non-positive real line" (line 11) but then switches to $\mathcal{H}_\beta^a \geq 0$ (line 26) without explaining the sign flip.
- **Suggested**: After the discriminant equation, add: "Since $\mathcal{L}$ generates a contraction semigroup, the spectrum of $\mathcal{H}_\beta$ is non-positive. For the frustration-free interpretation, we work with $-\mathcal{H}_\beta \geq 0$, whose ground state is $|\sqrt{\rho_\beta}\rangle$." Then use $-\mathcal{H}_\beta$ consistently in the frustration-free discussion, or define $\tilde{\mathcal{H}}_\beta := -\mathcal{H}_\beta$.
- **Rationale**: The reader will be confused by the sign inconsistency otherwise.

### Suggestion 2: Connect to thesis eq. (3.55) explicitly

- **Location**: Lines 16-19
- **Current**: The mixing time bound is stated and cited to [Appendix C.4]{chen2023efficient} as if new.
- **Suggested**: Add: "This recovers the spectral bound derived in Section~3.2 (eq.~3.55), now expressed in terms of the parent Hamiltonian gap."
- **Rationale**: The reader should see that the Parent Hamiltonian subsection is giving a new perspective on an existing result, not introducing a separate bound.

### Suggestion 3: Arnoldi/Lanczos complexity claims

- **Location**: Lines 42, the paragraph on Krylov methods
- **Current**: "For a general non-Hermitian Lindbladian, Arnoldi maintains a full upper-Hessenberg matrix requiring $\mathcal{O}(mn)$ storage and $\mathcal{O}(mn^2)$ operations per step, where $m$ is the Krylov dimension and $n$ the Hilbert space size. For the Hermitian parent Hamiltonian..., Lanczos iteration: a three-term recurrence needing only $\mathcal{O}(m + n)$ storage and $\mathcal{O}(mn)$ operations per step."
- **Issue**: The variable $n$ is overloaded. In the thesis, $n$ is the number of qubits (line 12 uses $2^n$ for the Hilbert space dimension). Here, "Hilbert space size" likely means $d^2 = 4^n$ (the doubled space). The Arnoldi storage claim of $\mathcal{O}(mn)$ is standard but should clarify that $n = d^2$ is the dimension of the vectorized space. Also, "operations per step" for Arnoldi is $\mathcal{O}(mn)$ for the orthogonalization (Gram-Schmidt against $m$ vectors of length $n$), plus one matrix-vector product. The draft says $\mathcal{O}(mn^2)$ which conflates the matrix-vector multiply cost (for dense matrices, $\mathcal{O}(n^2)$, but for sparse Lindbladians this is much less). The comparison is somewhat misleading since in practice both methods are dominated by the matrix-vector product (sparse matvec), not the recurrence overhead.
- **Suggested**: Use $d$ for the Hilbert space dimension to avoid overloading $n$. Clarify that the savings are in the recurrence overhead, not the matvec cost.

### Suggestion 4: Strengthen the DMRG paragraph

- **Location**: Lines 44, final paragraph
- **Current**: Mentions MPO representation and DMRG but only gives one citation (Verstraete et al. 2004).
- **Suggested**: The connection is: (1) area law $\Rightarrow$ efficient MPO representation, (2) frustration-freeness + MPO $\Rightarrow$ DMRG converges to global ground state. Step (1) relies on the (corrected) Firanko-Goldstein-Arad result. Step (2) is classical DMRG lore. Making these two steps explicit would strengthen the argument.

### Suggestion 5: Opening sentence is slightly misleading

- **Location**: Line 6: "The KMS detailed balance condition established in the previous subsections has a structural consequence that extends well beyond the algorithmic setting"
- **Current**: This suggests the parent Hamiltonian is a consequence of KMS detailed balance specifically.
- **Issue**: The discriminant transformation works for any detailed-balanced Lindbladian (GNS or KMS). What is special about the CKG/KMS case is frustration-freeness. The opening should be more precise.
- **Suggested**: "The KMS detailed balance condition established in the previous subsections, combined with the per-jump structure of the CKG Lindbladian, yields a structural consequence..."

---

## Connections and Enrichments

- **Quantum simulated annealing**: The CKG paper's Appendix C discusses quasi-adiabatic preparation of the purified Gibbs state via the gapped path $\mathcal{H}_{s\beta}$ for $s \in [0,1]$. This is briefly mentioned in the CKG abstract. Since the thesis plans an Algorithm subsection just before this one, mentioning the connection to quantum simulated annealing (Szegedy speedup) would provide a nice forward pointer. One sentence would suffice.

- **Bergamaschi-Chen all-temperatures result**: The recent result (arXiv:2510.08533) showing that 1D Hamiltonians admit a system-size-independent gap at ALL finite temperatures is a significant strengthening of the Rouze-Stilck Franca-Alhambra result (which only covers high temperature). If the draft mentions the former, it should note this distinction.

- **Connection to QuantumFurnace numerics**: The draft mentions (line 42) that Krylov methods on the parent Hamiltonian will be used "in our analysis." This should reference the specific thesis chapter (Part IV: QuantumFurnace.jl) where this is done. A forward reference would help the reader.

---

## What Works Well

- The three-part enumeration (spectral gap stability, detectability lemma, Knabe-type certificates) is well-organized and gives the reader a clear picture of the analytical payoffs.
- The connection between Arnoldi and Lanczos for the parent Hamiltonian is a nice practical observation that ties the abstract theory to the numerical work in the thesis.
- The overall length (~1.5 pages) is appropriate.
- The logical flow from discriminant definition $\to$ frustration-freeness $\to$ analytical tools $\to$ numerical methods is clean.

---

## Summary Scorecard

| Category | Rating | Notes |
|----------|--------|-------|
| Logical correctness | MINOR | Sign convention inconsistency (non-positive spectrum vs. PSD local terms); otherwise sound |
| Notation consistency | MINOR | $\mathcal{H}_\beta$ not yet defined in thesis; $n$ overloaded in Krylov paragraph |
| Citations | MAJOR | Area law misattributed to wrong authors; SDP hierarchy paper wrong authors; Bergamaschi-Chen cites wrong paper |
| Thesis consistency | MINOR | Should connect to eq. (3.55); sign flip needs reconciliation with thesis conventions |
| Exposition quality | MINOR | Good overall; sign convention and Krylov complexity claims need tightening |
