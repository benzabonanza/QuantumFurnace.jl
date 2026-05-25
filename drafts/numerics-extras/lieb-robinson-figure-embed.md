Insertion target: appendix of the numerics chapter (likely the locality /
operator-spread subsection). Drop the figure file into the thesis project's
`plots/` directory as `jump_spread_LR_n9.pdf`, then paste the block below.

# LaTeX embedding for the Lieb–Robinson appendix figure

The figure file to embed:

- `drafts/figures/numerics/jump_spread_LR_n9.pdf`  → copy into the thesis
  project's `plots/jump_spread_LR_n9.pdf`.

The wide aspect ratio (3 panels in a row, $\sim$ 1300:440) prefers
`width=\textwidth` placement and a top/bottom float.

## Figure block

```latex
\begin{figure}[t]
    \centering
    \includegraphics[width=\textwidth]{plots/jump_spread_LR_n9.pdf}
    \caption[Operator spread of $\sigma^x_5(t)$ on the 1D Heisenberg chain
    versus the Lieb--Robinson lightcone.]{%
        Per-site Pauli weight $w_k(t) = \tfrac{1}{d}\bigl[\|A(t)\|_{\mathrm{HS}}^2
        - \tfrac{1}{2}\|\mathrm{Tr}_k[A(t)]\|_{\mathrm{HS}}^2\bigr]$ of the
        Heisenberg-evolved operator $A(t) = U^\dagger \sigma^x_5\, U$ on the
        disordered 1D Heisenberg chain ($n = 9$, seed $= 46$, $Z+ZZ$
        $\varepsilon = 0.1$, periodic), at the two TimeDomain-relevant times
        $t_{\mathrm{median}}$ and $t_{\max}$ pulled from the truncated CKG
        OFT grid at each $\beta \in \{0.25, 0.5, 1\}$. Vertical red lines
        and the grey envelope mark the Lieb--Robinson lightcone with the
        tight spinon velocity $v_{\mathrm{LR}} = 2\pi$
        (1D~XXX,~Bethe~ansatz). At $\beta = 0.25$ the lightcone covers
        fewer than half the sites and $A(t_{\max})$ stays localised; at
        $\beta \geq 0.5$ the cone exceeds the chain length and the
        operator is effectively global. Locality intuition for the
        sampler is therefore cell-dependent on the canonical $\beta$ grid.}
    \label{fig:jump-spread-LR-n9}
\end{figure}
```

## Concise 3-line caption alternative

If the long-form caption above is too wordy for the appendix flow, use this
shorter version (matches the rest of the chapter's caption density):

```latex
\caption[Operator spread of $\sigma^x_5(t)$ vs.\ Lieb--Robinson lightcone.]{%
    Per-site Pauli weight $w_k(t)$ of $A(t) = U^\dagger \sigma^x_5 U$ on the
    1D Heisenberg chain ($n = 9$), at $t_{\mathrm{median}}$ and $t_{\max}$
    from the truncated CKG OFT grid for each $\beta$; vertical lines mark
    the tight Lieb--Robinson cone $v_{\mathrm{LR}} = 2\pi$. The lightcone
    fits inside the chain only at $\beta = 0.25$; at $\beta \geq 0.5$ the
    operator is effectively global.}
\label{fig:jump-spread-LR-n9}
```

## Cross-reference template

In the surrounding text:

```latex
Figure~\ref{fig:jump-spread-LR-n9} shows the Heisenberg-picture spread of a
single-site jump operator $\sigma^x_5$ on the $n = 9$ chain; see
\cite{BravyiHastingsVerstraete2006,Hastings2010LesHouches} for the
Lieb--Robinson bound and \cite{FaddeevTakhtajan1981} for the spinon
velocity that saturates the cone for 1D XXX.
```

(Bibtex entries for these citations live at the bottom of
`drafts/numerics-extras/lieb-robinson-notes.md`.)
