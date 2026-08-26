```@meta
EditURL = "../literate/tutorial_thermalize.jl"
```

# Finding a thermal state

This tutorial runs a short full-density-matrix simulation of a quantum Gibbs
sampler. The target state is
```math
\rho_\beta = \frac{e^{-\beta H}}{\operatorname{tr}(e^{-\beta H})}.
```
QuantumFurnace stores the example Hamiltonian with a rescaled spectrum, so
`beta_alg` below is the algorithm-side inverse temperature. The corresponding
physical value for the un-rescaled Hamiltonian is available through
`beta_phys(ham, beta_alg)`.

````@example tutorial_thermalize
using QuantumFurnace
using LinearAlgebra
````

## Load the system

Small reproducible Heisenberg fixtures are shipped with the package. Loading
one constructs its eigendecomposition and caches the Gibbs state at the
requested algorithm-side inverse temperature.

````@example tutorial_thermalize
n = 3
beta_alg = 10.0
ham = load_hamiltonian("heis", n; beta=beta_alg);
nothing #hide
````

## Construct jump operators

We use normalised single-site Pauli operators. `JumpOp` stores each operator
both in the computational basis and in the Hamiltonian eigenbasis used by the
energy-domain construction.

````@example tutorial_thermalize
local_jumps = ([X], [Y], [Z])
jump_norm = sqrt(length(local_jumps) * n)
jumps = JumpOp[]

for local_jump in local_jumps, site in 1:n
    op = Matrix(pad_term(local_jump, n, site)) / jump_norm
    op_eig = ham.eigvecs' * op * ham.eigvecs
    push!(jumps, JumpOp(op, op_eig, op == transpose(op), ishermitian(op)))
end
````

## Configure the sampler

`KMS()` includes the coherent correction required by this construction. The
energy domain uses a finite frequency grid for the dissipative integral; the
exact coherent term is constructed from the Hamiltonian's Bohr frequencies.
A short trajectory keeps this tutorial quick while still showing motion
towards the cached Gibbs state.

````@example tutorial_thermalize
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
);
nothing #hide
````

## Evolve and inspect the result

The default initial state is maximally mixed. The result stores the final
density matrix together with the sampled times and trace distances to the
Gibbs state.

````@example tutorial_thermalize
result = run_thermalize(jumps, config, ham; save_every=5);

result.time_steps
result.trace_distances
tr(result.final_dm)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

