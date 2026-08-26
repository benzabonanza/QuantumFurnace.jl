# QuantumFurnace.jl

QuantumFurnace.jl constructs and simulates detailed-balance Lindbladians for
quantum Gibbs sampling. It provides dense reference implementations,
matrix-free Krylov methods, and full-density-matrix channel simulations.

The package is under active development and is not yet registered in Julia's
General registry. Install it directly from GitHub:

```julia
import Pkg
Pkg.add(url="https://github.com/benzabonanza/QuantumFurnace.jl")
```

Start with [Finding a thermal state](generated/tutorial_thermalize.md), then
consult the [API reference](api.md) for the current public interfaces.
