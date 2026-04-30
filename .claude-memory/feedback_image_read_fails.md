---
name: Image Read fails with API 400 on plot PNGs
description: Reading generated plot PNGs (e.g. drafts/plots/*.png) via the Read tool can fail with "Could not process image" — verify plots by other means
type: feedback
originSessionId: fb027b0b-f0dc-4966-86ae-594c02506dce
---
The Read tool failed with `API Error: 400 {"type":"error","error":{"type":"invalid_request_error","message":"Could not process image"}}` when trying to view `drafts/plots/transition_weights.png` (qf-fx2.1, transition weights plot).

**Why:** Some plot PNGs produced by Plots.jl in this project (DPI 200, fontfamily "Computer Modern", PDF+PNG dual save) trigger this error in the multimodal vision pipeline. Cause is not fully diagnosed — possibly metadata, embedded fonts, or color profile.

**How to apply:**
- Do **not** rely on `Read` to visually inspect generated plot PNGs in `drafts/plots/`. If it fails once with that 400, do not retry the same file.
- Verify plots by other means: (a) read the script's stdout diagnostics (printed γ values, sanity asserts), (b) check the PDF/PNG file size is sensible, (c) ask the user to inspect the figure, (d) re-render to a different format (SVG) or different DPI if visual confirmation is essential.
- When generating new figures, print enough numerical diagnostics from the script that correctness can be judged without seeing the image.
