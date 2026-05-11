using Test
using QuantumFurnace
using LinearAlgebra

# Shared n=3 fixture: HamHam loaded from disk at β_alg = BETA = 10.0.
# We re-load via the public NamedTuple-keyword constructor to exercise the
# new `HamHam(raw; beta_phys=...)` path against the legacy positional
# `HamHam(raw, beta)` form.
const _BPC_SRC_ROOT = dirname(@__DIR__)
const _BPC_HAM_PATH = joinpath(_BPC_SRC_ROOT, "hamiltonians",
                               "heis_disordered_periodic_n3.bson")

# Reuse the test-helpers loader (which already lifts the legacy 14-field BSON
# schema into the NamedTuple needed by HamHam(raw; beta_phys=...)).
const _BPC_HAM_ALG10 = QuantumFurnace._load_hamiltonian_bson(_BPC_HAM_PATH, 10.0)

@testset "qf-6vr Task 1 — β_phys / β_alg helpers + HamHam keyword constructor" begin

    rescale = _BPC_HAM_ALG10.rescaling_factor
    @test rescale > 1.0

    @testset "(a) beta_alg / beta_phys helpers on HamHam" begin
        β_phys = 0.5
        β_alg_expected = β_phys * rescale
        @test beta_alg(_BPC_HAM_ALG10, β_phys) ≈ β_alg_expected
        @test beta_phys(_BPC_HAM_ALG10, β_alg_expected) ≈ β_phys
        # Round-trip identity
        for β in (0.1, 1.0, 5.0, 10.0)
            @test beta_alg(_BPC_HAM_ALG10, beta_phys(_BPC_HAM_ALG10, β)) ≈ β
            @test beta_phys(_BPC_HAM_ALG10, beta_alg(_BPC_HAM_ALG10, β)) ≈ β
        end
        # Type stability — helpers return ham's parametric T
        T = eltype(_BPC_HAM_ALG10.eigvals)
        @test typeof(beta_alg(_BPC_HAM_ALG10, 0.5)) === T
        @test typeof(beta_phys(_BPC_HAM_ALG10, 5.0)) === T
    end

    @testset "(b) HamHam(raw; beta_phys=…) byte-identity with HamHam(raw, β_alg)" begin
        # Build a NamedTuple from the loaded ham (re-using the eigvals / data
        # already present is the cleanest test fixture).
        raw_nt = (
            matrix             = Matrix{ComplexF64}(_BPC_HAM_ALG10.data),
            terms              = _BPC_HAM_ALG10.base_terms,
            base_coeffs        = _BPC_HAM_ALG10.base_coeffs,
            disordering_terms  = _BPC_HAM_ALG10.disordering_terms,
            disordering_coeffs = _BPC_HAM_ALG10.disordering_coeffs,
            eigvals            = _BPC_HAM_ALG10.eigvals,
            eigvecs            = _BPC_HAM_ALG10.eigvecs,
            nu_min             = _BPC_HAM_ALG10.nu_min,
            shift              = _BPC_HAM_ALG10.shift,
            rescaling_factor   = _BPC_HAM_ALG10.rescaling_factor,
            periodic           = _BPC_HAM_ALG10.periodic,
        )

        for β_phys in (0.25, 1.0, 3.0)
            β_alg = β_phys * raw_nt.rescaling_factor
            ham_alg  = HamHam(raw_nt, β_alg)
            ham_phys = HamHam(raw_nt; beta_phys=β_phys)
            # Gibbs states must be byte-identical (same β_alg in
            # `_gibbs_in_eigen`).
            @test ham_phys.gibbs ≈ ham_alg.gibbs atol=1e-15
            @test isapprox(Matrix(ham_phys.gibbs), Matrix(ham_alg.gibbs);
                           atol=1e-15)
            # All other fields are loaded from raw — identical
            @test ham_phys.rescaling_factor == ham_alg.rescaling_factor
            @test ham_phys.eigvals == ham_alg.eigvals
        end
    end

    @testset "(c) gibbs_state(ham, β_alg) consistency with stored ham.gibbs" begin
        # gibbs_state keeps legacy semantics: β positional argument is β_alg
        # (matches the rescaled eigvals stored on the HamHam). Test against
        # the eigenbasis variant since the ham stores its Gibbs state in the
        # eigenbasis (diagonal), not the computational basis.
        β_alg = 10.0  # matches the ham was loaded at
        ρ_eigen = gibbs_state_in_eigen(_BPC_HAM_ALG10, β_alg)
        @test isapprox(ρ_eigen, Matrix(_BPC_HAM_ALG10.gibbs); atol=1e-12)
    end
end
