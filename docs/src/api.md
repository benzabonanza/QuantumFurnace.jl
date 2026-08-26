# API reference

## Core types

```@docs
QuantumFurnace.Config
QuantumFurnace.HamHam
QuantumFurnace.JumpOp
QuantumFurnace.LindbladResults
QuantumFurnace.ThermalizeResults
```

## Hamiltonians

```@docs
QuantumFurnace.load_hamiltonian
QuantumFurnace.build_heis_1d
QuantumFurnace.build_tfim_2d
QuantumFurnace.beta_alg
QuantumFurnace.beta_phys
```

## Dynamics

```@docs
QuantumFurnace.construct_lindbladian
QuantumFurnace.run_lindblad
QuantumFurnace.run_thermalize
QuantumFurnace.predict_lindbladian_trajectory
QuantumFurnace.krylov_spectral_gap
QuantumFurnace.eigenmode_mixing_time
```

## Filters and validation

```@docs
QuantumFurnace.GaussianFilter
QuantumFurnace.DLLGaussianFilter
QuantumFurnace.DLLMetropolisFilter
QuantumFurnace.validate_config!
```
