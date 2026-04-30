---
name: Thesis colour palette (named)
description: Named colour palette for thesis plots and figures — refer to these colours by name (e.g. "use bordeaux", "make γ_M pinegreen")
type: reference
originSessionId: fb027b0b-f0dc-4966-86ae-594c02506dce
---
User-defined LaTeX/thesis colour palette. When the user names one of these in a request ("plot γ_M in bordeaux", "use sage for the Gaussian"), use the matching hex.

In Julia `Plots.jl`, hex strings work directly: `color = "#7A2E39"`. RGB tuple form (0–1 floats) is also given for readability.

| Name        | Hex       | RGB (0–1)            | Notes / palette role                  |
|-------------|-----------|----------------------|---------------------------------------|
| pinegreen   | `#2D5A3D` | (0.176, 0.353, 0.239)| deep green                            |
| bordeaux    | `#7A2E39` | (0.478, 0.180, 0.224)| deep wine red                         |
| dustyplum   | `#8E6F8C` | (0.557, 0.435, 0.549)| muted plum                            |
| deepplum    | `#735874` | (0.451, 0.345, 0.455)| darker muted plum                     |
| aubergine   | `#4F3B5C` | (0.310, 0.231, 0.361)| dark purple                           |
| slateblue   | `#5C7794` | (0.361, 0.467, 0.580)| muted blue                            |
| sage        | `#8B9F7E` | (0.545, 0.624, 0.494)| muted green                           |
| ochre       | `#B89143` | (0.722, 0.569, 0.263)| warm yellow-brown                     |
| terracotta  | `#B5654A` | (0.710, 0.396, 0.290)| warm orange-red                       |
| dustyteal   | `#5F8B8E` | (0.373, 0.545, 0.557)| muted teal                            |
| mustard     | `#C9A86A` | (0.788, 0.659, 0.416)| pale warm yellow                      |

LaTeX source (matches the thesis):
```latex
\definecolor{pinegreen}{HTML}{2D5A3D}
\definecolor{bordeaux}{HTML}{7A2E39}
\definecolor{dustyplum}{HTML}{8E6F8C}
\definecolor{deepplum}{HTML}{735874}
\definecolor{aubergine}{HTML}{4F3B5C}
\definecolor{slateblue}{HTML}{5C7794}
\definecolor{sage}{HTML}{8B9F7E}
\definecolor{ochre}{HTML}{B89143}
\definecolor{terracotta}{HTML}{B5654A}
\definecolor{dustyteal}{HTML}{5F8B8E}
\definecolor{mustard}{HTML}{C9A86A}
```

**Pairing notes (palette character):**
- Warm pairs: `bordeaux`/`terracotta`, `ochre`/`mustard`
- Cool/muted: `slateblue`, `dustyteal`, `sage`, `pinegreen`
- Purple family: `aubergine` (darkest) → `deepplum` → `dustyplum` (lightest)

**Use in Plots.jl examples:**
```julia
plot(x, y; color = "#7A2E39")            # bordeaux
plot!(x, y2; color = "#2D5A3D")           # pinegreen
# or with parse(Colorant, ...) for richer Color types:
using Colors
c_bordeaux = parse(Colorant, "#7A2E39")
```

**"Greyer variant" recipe.** Blend 50/50 with a neutral grey of the *same luminance* L = (R+G+B)/3:
```
greyer = 0.5 * (r, g, b) + 0.5 * (L, L, L)   # L = (r+g+b)/3
```
Examples: sage `#8B9F7E` → `#8D9686` ; slateblue `#5C7794` → `#6A7786`. Pulls the colour toward neutral while preserving visual weight. **Caveat:** the user found this too "gloomy" for primary thesis plots — prefer the brighter recipe below unless the curve is truly meant to fade away.

**"Brighter variant" recipe.** HSV value-boost: scale all three RGB channels so the maximum channel hits 0.80, preserving hue and saturation:
```
ratio = 0.80 / max(r, g, b)
brighter = (ratio * r, ratio * g, ratio * b)
```
Examples: sage `#8B9F7E` → `#B2CCA1` ; slateblue `#5C7794` → `#7FA4CC`. The result is more luminous and cheerful, while keeping the named colour's hue identity. Good for secondary/recessive curves that still need to read as colourful, not muddy.

For stronger brightening, raise the target V to 0.85 or 0.90 (caps at 1.0).
