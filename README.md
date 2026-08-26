# QuantumFurnace.jl

QuantumFurnace.jl is a Julia package for constructing and simulating quantum
Gibbs samplers based on detailed-balance Lindbladians. It supports dense
reference calculations, matrix-free Krylov methods, and full-density-matrix
channel simulations across Bohr, energy, time, and Trotter domains.

The package is under active development and should currently be treated as
pre-alpha research software.

## Features

- Construct KMS, GNS, and Ding--Li--Lin (DLL) Gibbs-sampling Lindbladians.
- Compare exact Bohr-domain constructions with energy-, time-, and
  Trotter-domain approximations.
- Compute fixed points and spectral gaps with dense or matrix-free methods.
- Simulate retained weak-measurement channels and track convergence to the
  Gibbs state.
- Build reproducible Heisenberg and transverse-field Ising Hamiltonians.

## Installation

QuantumFurnace.jl is not yet registered in Julia's General registry. Install
the development version directly from GitHub:

```julia
import Pkg
Pkg.add(url="https://github.com/benzabonanza/QuantumFurnace.jl")
```

Then load it with:

```julia
using QuantumFurnace
```

## Quick start

The following example loads a packaged three-qubit Hamiltonian, constructs
single-site Pauli jump operators, and runs a short KMS thermalisation
trajectory. The packaged Hamiltonian is rescaled, so `beta_alg` denotes the
algorithm-side inverse temperature.

```julia
using QuantumFurnace
using LinearAlgebra

n = 3
beta_alg = 10.0
ham = load_hamiltonian("heis", n; beta=beta_alg)

local_jumps = ([X], [Y], [Z])
jump_norm = sqrt(length(local_jumps) * n)
jumps = JumpOp[]

for local_jump in local_jumps, site in 1:n
    op = Matrix(pad_term(local_jump, n, site)) / jump_norm
    op_eig = ham.eigvecs' * op * ham.eigvecs
    push!(jumps, JumpOp(op, op_eig, op == transpose(op), ishermitian(op)))
end

config = Config(;
    sim=Thermalize(),
    domain=EnergyDomain(),
    construction=KMS(),
    num_qubits=n,
    with_linear_combination=true,
    beta=beta_alg,
    sigma=1 / beta_alg,
    a=beta_alg / 30,
    s=0.4,
    num_energy_bits_D=7,
    w0_D=0.05,
    mixing_time=0.1,
    delta=0.01,
    jump_selection=:sweep,
)

result = run_thermalize(jumps, config, ham; save_every=5)
result.trace_distances
```

`KMS()` includes its coherent correction by construction. Use
`beta_phys(ham, beta_alg)` to convert the example's inverse temperature back to
the un-rescaled physical Hamiltonian convention.

## Documentation

Tutorials, background material, and the API reference are available at
[benzabonanza.github.io/QuantumFurnace.jl](https://benzabonanza.github.io/QuantumFurnace.jl/).

## References

- C.-F. Chen, M. J. Kastoryano, F. G. S. L. Brandao, and A. Gilyen,
  "Quantum thermal state preparation," arXiv:2303.18224 (2023).
- C.-F. Chen, M. J. Kastoryano, and A. Gilyen, "An efficient and exact
  noncommutative quantum Gibbs sampler," arXiv:2311.09207 (2023).

## Citing

A formal software citation will accompany the first release. Until then,
please cite the repository URL together with the exact commit used in your
work.

## License

QuantumFurnace.jl is available under the MIT License.
