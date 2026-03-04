# Technology Stack: v2.2 Hamiltonian Simulation Time Counter

**Project:** QuantumFurnace.jl v2.2
**Researched:** 2026-03-04
**Confidence:** HIGH (all inputs are existing codebase quantities; no external library research needed)

## Executive Finding: No New Dependencies Required

The Hamiltonian simulation time counting functions are pure arithmetic on quantities the codebase already computes. They compute sums of `|t_j|` over QPE time grids, weighted by truncation masks. This is elementary Julia: `sum(abs, collection)`. No numerical libraries, no optimization, no special algorithms.

## Recommended Stack

### Core Framework

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| Julia stdlib only | >= 1.9 (existing requirement) | `sum`, `abs`, `ceil`, `length`, `log2`, `size` | All counting operations are stdlib arithmetic |

### Existing Dependencies Used (No Changes to Project.toml)

| Dependency | Already In Project | Used By Counting For |
|------------|-------------------|---------------------|
| SpecialFunctions.jl | Yes (`erfc`) | Smooth Metropolis transition weight evaluation (needed if computing energy truncation count) |
| LinearAlgebra | Yes (stdlib) | `size(ham.data, 1)` to extract dimension |

### Not Needed

| Technology | Why Not Needed |
|------------|----------------|
| LsqFit.jl | No curve fitting in time counting |
| Roots.jl | No root finding |
| QuadGK.jl | QPE grids are discrete finite sums, not continuous integrals |
| BSON.jl | No serialization of SimulationTimeBudget (lightweight transient result) |
| KrylovKit.jl | No eigensolves |
| FINUFFT | NUFFT computes OFT values, not time costs (category error to use it for counting) |
| SparseArrays | No matrix construction |
| Any new package | `sum(abs, v)` does not need a dependency |

## What the Counting Functions Compute

### OFT Ham Sim Time

The OFT computes `A(omega) = sum_j prefactor(t_j) * e^{iHt_j} * A * e^{-iHt_j}`. Each term requires Hamiltonian simulation for time `|t_j|`. The total per-energy OFT time:

```
OFT_time_per_energy = sum(|t_j| for t_j in full_QPE_time_grid)
                    = t0 * sum(|j| for j in -N/2 : N/2-1)
                    = t0 * N^2 / 4    (closed form, N = 2^r)
```

This uses the **full** QPE grid (what the quantum computer runs), NOT the Gaussian-truncated grid (a classical simulation optimization).

The total OFT cost per step sums over truncated energy grid points:
```
OFT_time_per_step = n_truncated_energies * OFT_time_per_energy
```

### B Coherent Term Ham Sim Time

For KMS construction, the B term involves double sums over truncated b_minus and b_plus time dictionaries:

```
B_cost = sum(4 * |s * beta| for s in keys(b_plus))     # inner sum
       + sum(2 * |t / sigma| for t in keys(b_minus))    # outer sum
```

The number of non-truncated entries depends on beta, sigma, transition weight type, and the time grid.

### Per-Step and Total Cost

```
per_step_cost = B_cost + 2 * OFT_time_per_step
total_cost = ceil(mixing_time / delta) * per_step_cost
```

The factor of 2 for OFT: one weak-measurement U and one controlled-U-dagger in the quantum circuit.

## Existing Internal Functions to Reuse

| Function | File | Used For |
|----------|------|----------|
| `_create_energy_labels(r, w0)` | `energy_domain.jl` | Full QPE energy grid |
| `_truncate_energy_labels(labels, config)` | `energy_domain.jl` | Truncated energy count (needs transition weight) |
| `_compute_truncated_func(func, labels, ...)` | `coherent.jl` | Truncated b_plus/b_minus dicts for B cost |
| `_compute_b_minus` | `coherent.jl` | b_minus function values |
| `_compute_b_plus` / `_compute_b_plus_metro` / `_compute_b_plus_smooth` | `coherent.jl` | b_plus function values (variant depends on transition weight) |
| `_select_b_plus_calculator(config)` | `furnace_utensils.jl` | Selects correct b_plus variant |
| `pick_transition(config)` / `pick_transition(config, w)` | `energy_domain.jl` | Transition weight for energy truncation |
| `with_coherent(construction)` | `structs.jl` | Whether B term contributes |

**Note on reuse strategy:** Some of these functions take a `Config` struct. The time counter can either:
1. Construct a minimal Config for reuse (coupling to Config type hierarchy)
2. Reimplement the logic inline (3-way dispatch on transition weight type)

Both are viable. Option (2) is cleaner for the public API but duplicates ~20 lines of formula code. Option (1) avoids duplication but requires assembling a Config. The ARCHITECTURE.md recommends option (2) for the public API but either works.

## Alternatives Considered

| Category | Recommended | Alternative | Why Not |
|----------|-------------|-------------|---------|
| Transition weight dispatch | Symbol keyword (`:gaussian`, `:metropolis`, `:smooth_metro`) | Reuse `pick_transition(Config, w)` | Would require constructing a Config just for weight evaluation |
| B term cost computation | Actually compute truncated dicts (cheap) | Estimate count analytically | Analytical estimate would be inaccurate; actual computation takes microseconds |
| Result struct | New `SimulationTimeBudget` immutable struct | NamedTuple | Struct provides type safety, docstring, and matches FitResult/MixingTimeEstimate pattern |
| Float precision | Float64 only | Parametric `{T<:AbstractFloat}` | Counting is microsecond arithmetic; Float64 is already overkill for sums of ~4000 terms |
| Caching | None (compute on the fly) | Cache results per parameter set | Computation takes microseconds; caching adds complexity for zero benefit |

## Integration with Module

```julia
# In src/QuantumFurnace.jl, after existing includes:
include("simulation_time.jl")  # New file

# New exports:
export SimulationTimeBudget, compute_simulation_time
```

No `using` additions needed. The file uses only existing internal functions and Julia stdlib.

## Version Compatibility

No changes to Julia version requirements or compat bounds:

| Requirement | Current | Impact |
|-------------|---------|--------|
| Julia >= 1.9 | Already set | No change |
| All [compat] entries | Unchanged | Counting functions touch no external packages |

## Summary: Dependency Changes

```toml
# Project.toml: NO CHANGES
# No new [deps], no new [compat], no new [extras]
```

**Total new external dependencies: 0**
**Total new internal dependencies: 0** (all called functions already exist)

## Sources

- `src/QuantumFurnace.jl` (import statements, include order)
- `src/energy_domain.jl` (transition weight formulas, `_create_energy_labels`, `_truncate_energy_labels`)
- `src/coherent.jl` (`B_time`, `B_trotter`, `_compute_b_minus`, `_compute_b_plus*`, `_compute_truncated_func`)
- `src/furnace_utensils.jl` (`_precompute_labels`, `_select_b_plus_calculator`)
- `src/structs.jl` (`Config{S,D,C,T}` fields, `with_coherent` trait)
- `src/hamiltonian.jl` (`HamHam{T}` struct fields)
- `src/misc_tools.jl` (`validate_config!` -- Fourier relation enforcement)
- `Project.toml` (current dependency list)
