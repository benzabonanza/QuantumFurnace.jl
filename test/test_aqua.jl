using Aqua
using QuantumFurnace

@testset "Aqua.jl package quality" begin
    Aqua.test_all(QuantumFurnace;
        # BenchmarkTools and Profile are intentionally listed in [deps] so that
        # `scripts/` (run with `--project`) can `using` them for benchmarking /
        # profiling; `src/` never imports them. Aqua flags this as "stale";
        # ignore those two specifically.
        stale_deps = (ignore = [:BenchmarkTools, :Profile],),
    )
end
