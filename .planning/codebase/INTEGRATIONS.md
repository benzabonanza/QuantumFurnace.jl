# External Integrations

**Analysis Date:** 2026-02-13

## APIs & External Services

**None detected** - This is a pure computational physics package with no external API integrations.

## Data Storage

**File Storage:**
- Local filesystem only
- Precomputed Hamiltonians stored as BSON files in `hamiltonians/` directory
- Runtime results saved to `results/` and `simulations/` directories

**Databases:**
- Not used - All state is either computed in-memory or serialized to disk as BSON

**Caching:**
- Not applicable - No external caching layer used

## Authentication & Identity

**Not applicable** - No user authentication or identity management required. Package is self-contained computational software.

## Monitoring & Observability

**Error Tracking:**
- Not integrated - Errors handled via Julia's native exception system

**Logs:**
- Approach: Standard output via `@printf()` and `@warn` macros
- Output goes to `stdout`/`stderr`
- Example in `src/energy_domain.jl`:
  ```julia
  @warn "Lower bound cutoff not found for energies, using default range."
  ```

**Progress Tracking:**
- ProgressMeter used for long-running simulations to display progress bars
- Imported in `src/QuantumFurnace.jl`: `using ProgressMeter`

## CI/CD & Deployment

**Hosting:**
- Package registered in Julia package ecosystem
- Documentation hosted on GitHub Pages: https://tembence.github.io/QuantumFurnace.jl/
- GitHub repository: https://github.com/tembence/QuantumFurnace.jl

**CI Pipeline:**
- Documentation deployment configured in `docs/make.jl`
- Deploydocs to GitHub Pages on `main` branch (devbranch)
- No external CI service detected (GitHub Actions config may exist in `.github/` but not analyzed)

**HPC Submission:**
- SLURM batch script present: `simulations/run_julia.sbatch`
- Enables distributed computing via ClusterManagers on HPC clusters

## Environment Configuration

**Required env vars:**
- None required
- All configuration passed directly via struct initialization

**Configuration Method:**
- Configuration structs in `src/structs.jl`:
  - `LiouvConfig` and `LiouvConfigGNS` - For Liouvillian construction
  - `ThermalizeConfig` and `ThermalizeConfigGNS` - For thermalization simulation
  - `TrotterDomain`, `TimeDomain`, `BohrDomain`, `EnergyDomain` - Domain representations

**Example Configuration:**
```julia
config = ThermalizeConfig(
    num_qubits = num_qubits,
    with_coherent = with_coherent,
    with_linear_combination = with_linear_combination,
    domain = domain,
    beta = beta,
    a = a,
    b = b,
    num_energy_bits = num_energy_bits,
    w0 = w0,
    t0 = t0,
    mixing_time = mixing_time_bound,
    delta = delta,
)
```

**Secrets location:**
- Not applicable - No secrets or credentials used

## Webhooks & Callbacks

**Incoming:**
- Not applicable

**Outgoing:**
- Documentation auto-deployment via `deploydocs()` in `docs/make.jl` to GitHub Pages
- Integrated with Documenter.jl: `repo = "github.com/tembence/QuantumFurnace.jl.git"`

## Data I/O

**Reading:**
- Precomputed Hamiltonian matrices from BSON files via `load_hamiltonian()` in `src/misc_tools.jl`
  ```julia
  function load_hamiltonian(type::String, num_qubits::Int)
      project_root = Pkg.project().path |> dirname
      data_dir = joinpath(project_root, "hamiltonians")
      output_filename = join([type, "disordered", "periodic", "n$num_qubits"], "_") * ".bson"
      ham_path = joinpath(data_dir, output_filename)
      bson_hamiltonian_data = BSON.load(ham_path)
      return bson_hamiltonian_data[:hamiltonian]
  end
  ```

**Writing:**
- Simulation results saved as BSON via generated filenames from `generate_filename()` functions
- Filenames encode configuration: `alg_DB_TimeDomain_n=4_beta=10.0_B_a=0.0_b=0.0_mix=10.0.bson`
- Example usage in `src/kossakowski.jl`:
  ```julia
  BSON.load(ham_path)
  ```

## Integration Points Summary

**Self-contained Design:**
- No external API dependencies
- No authentication/authorization requirements
- No database connectivity
- All dependencies are scientific computing libraries
- Pure Julia package suitable for offline research environments

**Cluster Integration:**
- Via Julia Distributed system and ClusterManagers
- Enables embarrassingly-parallel trajectory simulations across HPC nodes

---

*Integration audit: 2026-02-13*
