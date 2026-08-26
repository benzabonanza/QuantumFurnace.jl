using Aqua
using QuantumFurnace

@testset "Aqua.jl package quality" begin
    Aqua.test_all(QuantumFurnace;
        # BenchmarkTools and Profile are intentionally listed in [deps] so that
        # `scripts/` (run with `--project`) can `using` them for benchmarking /
        # profiling; `src/` never imports them. Aqua flags this as "stale";
        # ignore those two specifically.
        stale_deps = (ignore = [:BenchmarkTools, :Profile],),
        # A cold, serialized CI precompile can take longer than Aqua's
        # 10-second default after the wrapper package has loaded. Keep the
        # persistent-task check enabled, but allow the process time to finish.
        persistent_tasks = (tmax = 120,),
    )
end
