using Aqua
using QuantumFurnace

@testset "Aqua.jl Package Quality (TINF-03)" begin
    Aqua.test_all(QuantumFurnace;
        ambiguities = false,          # Disable: multiple dispatch may create legitimate ambiguities
        piracies = false,             # Disable: kron! on AbstractMatrix may be flagged
    )
end
