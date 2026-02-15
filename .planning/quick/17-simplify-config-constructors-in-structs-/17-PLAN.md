---
phase: 17-simplify-config-constructors
plan: 01
type: execute
wave: 1
depends_on: []
files_modified: [src/structs.jl]
autonomous: true

must_haves:
  truths:
    - "LiouvConfigGNS keyword construction works identically to before"
    - "ThermalizeConfigGNS keyword construction works identically to before"
    - "with_coherent=true is still rejected at construction time for both GNS configs"
    - "Defaults (with_coherent=false, beta=1.0, sigma=0.1, etc.) match current behavior"
    - "Package loads and existing tests pass"
  artifacts:
    - path: "src/structs.jl"
      provides: "All 4 config structs with consistent @kwdef pattern"
      contains: "@kwdef struct LiouvConfigGNS"
    - path: "src/structs.jl"
      provides: "All 4 config structs with consistent @kwdef pattern"
      contains: "@kwdef struct ThermalizeConfigGNS"
  key_links:
    - from: "src/structs.jl (LiouvConfigGNS)"
      to: "src/bohr_domain.jl, src/energy_domain.jl, src/misc_tools.jl"
      via: "type dispatch on LiouvConfigGNS"
      pattern: "config::LiouvConfigGNS"
    - from: "src/structs.jl (ThermalizeConfigGNS)"
      to: "src/bohr_domain.jl, src/energy_domain.jl, src/misc_tools.jl"
      via: "type dispatch on ThermalizeConfigGNS"
      pattern: "config::ThermalizeConfigGNS"
---

<objective>
Simplify LiouvConfigGNS and ThermalizeConfigGNS from 3-part constructor patterns (manual struct + inner constructor + standalone keyword constructor) to clean @kwdef pattern with inner constructor validation.

Purpose: Make all 4 config structs use a consistent @kwdef pattern, reducing ~114 lines of boilerplate to ~54 lines while preserving identical public API and with_coherent validation.

Output: Simplified src/structs.jl with consistent constructor patterns across all config structs.
</objective>

<execution_context>
@./.claude/get-shit-done/workflows/execute-plan.md
@./.claude/get-shit-done/templates/summary.md
</execution_context>

<context>
@src/structs.jl
</context>

<tasks>

<task type="auto">
  <name>Task 1: Convert LiouvConfigGNS to @kwdef with inner constructor</name>
  <files>src/structs.jl</files>
  <action>
Replace lines 114-166 (the `LiouvConfigGNS` struct definition, inner constructor, and standalone keyword constructor) with a single `@kwdef struct` block plus a one-line bridging outer constructor.

The new code should be:

```julia
@kwdef struct LiouvConfigGNS{D <: AbstractDomain, T <: AbstractFloat} <: AbstractLiouvConfig{D,T}
    num_qubits::Int64
    with_coherent::Bool = false
    with_linear_combination::Bool
    domain::D
    beta::T = 1.0
    sigma::T = 0.1
    gaussian_parameters::Union{Tuple{T, T}, Tuple{Nothing, Nothing}} = (nothing, nothing)
    a::Union{T, Nothing} = nothing
    b::Union{T, Nothing} = nothing
    num_energy_bits::Union{Int, Nothing} = nothing
    t0::Union{T, Nothing} = nothing
    w0::Union{T, Nothing} = nothing
    eta::Union{T, Nothing} = nothing
    num_trotter_steps_per_t0::Union{Int, Nothing} = nothing

    function LiouvConfigGNS{D,T}(
        num_qubits, with_coherent, with_linear_combination, domain,
        beta, sigma, gaussian_parameters, a, b,
        num_energy_bits, t0, w0, eta, num_trotter_steps_per_t0
    ) where {D, T}
        with_coherent && error("GNS configs must have with_coherent=false")
        new{D,T}(num_qubits, with_coherent, with_linear_combination, domain,
               beta, sigma, gaussian_parameters, a, b,
               num_energy_bits, t0, w0, eta, num_trotter_steps_per_t0)
    end
end

# Bridge for @kwdef: unparameterized positional -> parameterized inner constructor
function LiouvConfigGNS(
    num_qubits::Int64, with_coherent::Bool, with_linear_combination::Bool, domain::D,
    beta::T, sigma::T, gaussian_parameters, a, b,
    num_energy_bits, t0, w0, eta, num_trotter_steps_per_t0
) where {D <: AbstractDomain, T <: AbstractFloat}
    LiouvConfigGNS{D,T}(num_qubits, with_coherent, with_linear_combination, domain,
                      beta, sigma, gaussian_parameters, a, b,
                      num_energy_bits, t0, w0, eta, num_trotter_steps_per_t0)
end
```

Key points:
- Defaults match the old keyword constructor: with_coherent=false, beta=1.0, sigma=0.1, gaussian_parameters=(nothing,nothing), optional fields default to nothing.
- Inner constructor keeps the `with_coherent && error(...)` validation.
- The bridging outer constructor is needed because @kwdef generates `LiouvConfigGNS(args...)` (unparameterized) but the inner constructor only defines `LiouvConfigGNS{D,T}(args...)`. Without the bridge, keyword construction fails with MethodError.
- Keep the docstring above the struct unchanged.
- Delete the old comment "# Keyword constructor for LiouvConfigGNS (replaces @kwdef)" and the standalone function.
  </action>
  <verify>
Run: `julia --project=. -e 'using QuantumFurnace; c = LiouvConfigGNS(num_qubits=2, with_linear_combination=false, domain=BohrDomain(), beta=1.0, sigma=0.1); println(typeof(c)); try LiouvConfigGNS(num_qubits=2, with_coherent=true, with_linear_combination=false, domain=BohrDomain(), beta=1.0, sigma=0.1) catch e; println("Validation OK: ", e.msg) end'`

Expected: prints `LiouvConfigGNS{BohrDomain, Float64}` then `Validation OK: GNS configs must have with_coherent=false`
  </verify>
  <done>LiouvConfigGNS uses @kwdef, keyword construction works with same defaults, with_coherent=true is rejected</done>
</task>

<task type="auto">
  <name>Task 2: Convert ThermalizeConfigGNS to @kwdef with inner constructor</name>
  <files>src/structs.jl</files>
  <action>
Replace the ThermalizeConfigGNS struct definition, inner constructor, and standalone keyword constructor (old lines 211-271) with a single `@kwdef struct` block plus a one-line bridging outer constructor.

The new code should be:

```julia
@kwdef struct ThermalizeConfigGNS{D <: AbstractDomain, T <: AbstractFloat} <: AbstractThermalizeConfig{D,T}
    num_qubits::Int64
    with_coherent::Bool = false
    with_linear_combination::Bool
    domain::D
    beta::T = 1.0
    sigma::T = 0.1
    gaussian_parameters::Union{Tuple{T, T}, Tuple{Nothing, Nothing}} = (nothing, nothing)
    a::Union{T, Nothing} = nothing
    b::Union{T, Nothing} = nothing
    num_energy_bits::Union{Int, Nothing} = nothing
    t0::Union{T, Nothing} = nothing
    w0::Union{T, Nothing} = nothing
    eta::Union{T, Nothing} = nothing
    num_trotter_steps_per_t0::Union{Int, Nothing} = nothing

    mixing_time::T
    delta::T

    function ThermalizeConfigGNS{D,T}(
        num_qubits, with_coherent, with_linear_combination, domain,
        beta, sigma, gaussian_parameters, a, b,
        num_energy_bits, t0, w0, eta, num_trotter_steps_per_t0,
        mixing_time, delta
    ) where {D, T}
        with_coherent && error("GNS configs must have with_coherent=false")
        new{D,T}(num_qubits, with_coherent, with_linear_combination, domain,
               beta, sigma, gaussian_parameters, a, b,
               num_energy_bits, t0, w0, eta, num_trotter_steps_per_t0,
               mixing_time, delta)
    end
end

# Bridge for @kwdef: unparameterized positional -> parameterized inner constructor
function ThermalizeConfigGNS(
    num_qubits::Int64, with_coherent::Bool, with_linear_combination::Bool, domain::D,
    beta::T, sigma::T, gaussian_parameters, a, b,
    num_energy_bits, t0, w0, eta, num_trotter_steps_per_t0,
    mixing_time::T, delta::T
) where {D <: AbstractDomain, T <: AbstractFloat}
    ThermalizeConfigGNS{D,T}(num_qubits, with_coherent, with_linear_combination, domain,
                           beta, sigma, gaussian_parameters, a, b,
                           num_energy_bits, t0, w0, eta, num_trotter_steps_per_t0,
                           mixing_time, delta)
end
```

Key points:
- Same pattern as LiouvConfigGNS: @kwdef + inner constructor + bridging outer.
- Has two extra fields (mixing_time, delta) with NO defaults (required, same as current).
- Defaults match old keyword constructor: with_coherent=false, beta=1.0, sigma=0.1, etc.
- Keep the docstring above the struct unchanged.
- Delete the old comment "# Keyword constructor for ThermalizeConfigGNS (replaces @kwdef)" and the standalone function.
  </action>
  <verify>
Run: `julia --project=. -e 'using QuantumFurnace; c = ThermalizeConfigGNS(num_qubits=2, with_linear_combination=false, domain=BohrDomain(), beta=1.0, sigma=0.1, mixing_time=10.0, delta=0.1); println(typeof(c)); try ThermalizeConfigGNS(num_qubits=2, with_coherent=true, with_linear_combination=false, domain=BohrDomain(), beta=1.0, sigma=0.1, mixing_time=10.0, delta=0.1) catch e; println("Validation OK: ", e.msg) end'`

Expected: prints `ThermalizeConfigGNS{BohrDomain, Float64}` then `Validation OK: GNS configs must have with_coherent=false`

Then run full package load check: `julia --project=. -e 'using QuantumFurnace; println("Package loads OK")'`
  </verify>
  <done>ThermalizeConfigGNS uses @kwdef, keyword construction works with same defaults, with_coherent=true is rejected, package loads cleanly</done>
</task>

</tasks>

<verification>
1. Package loads: `julia --project=. -e 'using QuantumFurnace; println("OK")'`
2. LiouvConfigGNS keyword construction with defaults works
3. ThermalizeConfigGNS keyword construction with defaults works
4. Both GNS configs reject with_coherent=true
5. All 4 config structs now use @kwdef pattern (grep for `@kwdef struct.*Config`)
6. No standalone keyword constructor functions remain for GNS configs
</verification>

<success_criteria>
- src/structs.jl has 4 config structs, all using @kwdef
- GNS configs are ~30 lines each (down from ~55 each), net reduction of ~50 lines
- Keyword API is identical: same field names, same defaults, same required fields
- with_coherent validation still fires at construction time
- Package compiles and loads without errors
</success_criteria>

<output>
After completion, create `.planning/quick/17-simplify-config-constructors-in-structs-/17-SUMMARY.md`
</output>
