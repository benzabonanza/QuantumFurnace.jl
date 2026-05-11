using Aqua
using QuantumFurnace

@testset "Aqua.jl Package Quality (TINF-03)" begin
    Aqua.test_all(QuantumFurnace;
        ambiguities = false,          # Disable: multiple dispatch may create legitimate ambiguities
        piracies = false,             # Disable: kron! on AbstractMatrix may be flagged
        # BenchmarkTools and Profile are intentionally listed in [deps] so that
        # `scripts/` (run with `--project`) can `using` them for benchmarking /
        # profiling; `src/` never imports them. Aqua flags this as "stale";
        # ignore those two specifically.
        stale_deps = (ignore = [:BenchmarkTools, :Profile],),
    )
end
