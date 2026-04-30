---
name: circuits folder stays .tex-only
description: Never leave .aux/.log/.pdf/.DS_Store or other build artefacts in drafts/circuits/
type: feedback
originSessionId: a3d4f276-d388-4028-ad8e-cd7d37bfa72d
---
`drafts/circuits/` must contain only `.tex` source files. No `.aux`, `.log`, `.pdf`, `.out`, `.DS_Store`, or other build artefacts.

**Why:** User maintains the folder as a clean source directory; thesis build happens elsewhere. Artefacts clutter the tree and pollute git status.

**How to apply:** When compiling a circuit to verify it builds, either (a) `cd` to a tmp dir and `pdflatex` with `-output-directory`, or (b) run `pdflatex` in-place and clean up (`rm *.aux *.log *.pdf`) immediately afterwards. Never leave them behind.
