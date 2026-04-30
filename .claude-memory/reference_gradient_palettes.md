---
name: Thesis gradient palettes (cold / warm / diverging)
description: Three gradient palettes (cold, warm, diverging) provided by the user for heatmaps and continuous-scale plots — use these as Plots.cgrad sources instead of `:inferno` / `:viridis` defaults
type: reference
---

User-supplied gradient palettes for thesis plots. Use these whenever a script
needs a `cgrad`/colormap (heatmaps, density plots, continuous overlays) rather
than the matplotlib defaults. Sources are in
`supplementary-informations/colors-{warm,cold,opposite}-gradient.png`.

## Cold (sequential, light → dark)
Light mint → deep navy. Good for heatmaps where the *background* should read
as a soft pastel rather than pure black.

| # | Hex       | Notes                       |
|---|-----------|-----------------------------|
| 1 | `#B7E6A5` | light mint (zero-anchor)    |
| 2 | `#7CCBA2` | green-teal                  |
| 3 | `#46AEA0` | teal                        |
| 4 | `#089099` | deep teal                   |
| 5 | `#00718B` | petrol blue                 |
| 6 | `#045275` | dark teal-blue              |
| 7 | `#003147` | deep navy (peak)            |

```julia
const COLD_GRAD = cgrad(["#B7E6A5", "#7CCBA2", "#46AEA0", "#089099",
                         "#00718B", "#045275", "#003147"])
```

## Warm (sequential, light → dark)
Cream → deep magenta-purple. Most "elegant"; cream base reads as paper, peaks
read as ink. Good for matrices with sparse high-value structure (e.g. diagonal-
dominant Kossakowski).

| # | Hex       | Notes                       |
|---|-----------|-----------------------------|
| 1 | `#FCE1A4` | cream / butter              |
| 2 | `#FABF7B` | warm sand                   |
| 3 | `#F08F6E` | coral                       |
| 4 | `#E05C5C` | rose                        |
| 5 | `#D12959` | crimson                     |
| 6 | `#AB1866` | magenta                     |
| 7 | `#6E005F` | deep mulberry (peak)        |

```julia
const WARM_GRAD = cgrad(["#FCE1A4", "#FABF7B", "#F08F6E", "#E05C5C",
                         "#D12959", "#AB1866", "#6E005F"])
```

## Diverging (opposite, dark teal ↔ deep purple)
Use only for *signed* quantities centred on zero (e.g. real part of an off-
diagonal coherence). The light cream `#FCDE9C` is the zero-anchor.

| # | Hex       | Side       |
|---|-----------|------------|
| 1 | `#045275` | negative   |
| 2 | `#089099` |            |
| 3 | `#7CCBA2` |            |
| 4 | `#FCDE9C` | zero       |
| 5 | `#F0746E` |            |
| 6 | `#DC3977` |            |
| 7 | `#7C1D6F` | positive   |

```julia
const DIVERGING_GRAD = cgrad(["#045275", "#089099", "#7CCBA2", "#FCDE9C",
                              "#F0746E", "#DC3977", "#7C1D6F"])
```

## When to pick which
- **Sequential (warm or cold)** — non-negative quantities (Kossakowski α, |γ|², densities, populations).
- **Diverging** — quantities that go negative as well as positive (off-diagonal real parts, error w.r.t. a reference, perturbations).
- **Reversed** (`cgrad(...; rev=true)`) — when you need a dark zero-anchor (e.g. low values fade into background) but want to avoid the harsh pure black of `:inferno`.
