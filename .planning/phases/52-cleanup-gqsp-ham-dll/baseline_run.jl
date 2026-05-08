using Test
using QuantumFurnace
using LinearAlgebra
using Printf

include("../../../test/test_helpers.jl")

const TARGETS = [
    "test_hamiltonian.jl",
    "test_dll_filter.jl",
    "test_dll_dissipator.jl",
    "test_dll_coherent.jl",
    "test_dll_kossakowski.jl",
    "test_dll_kms_db.jl",
    "test_gqsp_config.jl",
    "test_gqsp_polynomial.jl",
    "test_gqsp_thermalize.jl",
]

results = Tuple{String, Float64, Bool}[]

for target in TARGETS
    path = joinpath(@__DIR__, "..", "..", "..", "test", target)
    @info "Running $target"
    t0 = time()
    ok = true
    try
        Base.include(@__MODULE__, path)
    catch err
        @warn "Test errored" target err
        ok = false
    end
    push!(results, (target, time() - t0, ok))
end

println("\n=== Baseline timing (cleanup-epic targets) ===")
total = 0.0
for (target, t, ok) in results
    @printf "%-32s  %8.2f s  %s\n" target t (ok ? "OK" : "FAIL")
    total += t
end
@printf "%-32s  %8.2f s\n" "TOTAL" total
