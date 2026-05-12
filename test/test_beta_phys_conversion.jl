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

    @testset "(d) Config accessors: beta_alg, beta_phys" begin
        # Minimal valid Config{Lindbladian, BohrDomain, KMS, Float64}: BohrDomain
        # has no Fourier register requirements, and a single Gaussian transition
        # is enough to satisfy `validate_config!`.
        rescale = _BPC_HAM_ALG10.rescaling_factor
        β_phys_val = 1.0
        β_alg_val = β_phys_val * rescale

        cfg = Config(
            sim = Lindbladian(), domain = BohrDomain(), construction = KMS(),
            num_qubits = 3, with_linear_combination = true,
            beta = β_alg_val, beta_phys = β_phys_val,
            sigma = 1.0 / β_alg_val,
            a = 0.0, s = 0.25,
        )
        @test beta_alg(cfg) == β_alg_val
        @test beta_phys(cfg) == β_phys_val

        # Config without β_phys: returns nothing
        cfg_noφ = Config(
            sim = Lindbladian(), domain = BohrDomain(), construction = KMS(),
            num_qubits = 3, with_linear_combination = true,
            beta = β_alg_val, sigma = 1.0 / β_alg_val, a = 0.0, s = 0.25,
        )
        @test beta_alg(cfg_noφ) == β_alg_val
        @test beta_phys(cfg_noφ) === nothing
    end

    @testset "(eX) fit_scaling reader: prefers :beta_phys, tracks beta_kind" begin
        # Synthetic vector of NamedTuples matching the qf-6vr sidecar schema.
        rescale = 3.0
        rows_phys = NamedTuple[]
        for n in 3:8, βp in (1.0, 2.0, 3.0)
            βa = βp * rescale
            push!(rows_phys, (
                n = n, beta = βa, beta_alg = βa, beta_phys = βp,
                rescaling_factor = rescale,
                mixing_time = Float64(n)^2 * βp^1.5 * (1.0 + 0.0),
                mixing_time_source = :extrapolated,
            ))
        end
        fits_phys = fit_scaling(rows_phys)
        @test fits_phys[:M0].beta_kind === :phys
        # Auto-mode prefers :beta_phys when present
        @test all(b -> b ∈ (1.0, 2.0, 3.0), fits_phys[:M0].beta_values)

        # Legacy rows with only :beta — reader picks :alg.
        rows_legacy = NamedTuple[]
        for n in 3:8, β in (5.0, 10.0, 20.0)
            push!(rows_legacy, (
                n = n, beta = β,
                mixing_time = Float64(n)^2 * β^1.5,
                mixing_time_source = :extrapolated,
            ))
        end
        fits_legacy = fit_scaling(rows_legacy)
        @test fits_legacy[:M0].beta_kind === :alg
        @test all(b -> b ∈ (5.0, 10.0, 20.0), fits_legacy[:M0].beta_values)

        # Explicit beta_kind=:alg on β_phys-aware data: reads :beta_alg.
        fits_force_alg = fit_scaling(rows_phys; beta_kind = :alg)
        @test fits_force_alg[:M0].beta_kind === :alg
        @test all(b -> b ∈ (3.0, 6.0, 9.0), fits_force_alg[:M0].beta_values)
    end

    @testset "(f) Lindbladian built via β_phys fixes the physical Gibbs state" begin
        # The hams cached on disk hold the *rescaled* spectrum so the simulator
        # stays inside [0, 0.45]; the un-rescaled (physical) Hamiltonian is
        #   H_phys = ham.rescaling_factor · H_alg + ham.shift · I
        # so its spectrum is
        #   eigvals_phys = rescaling_factor · eigvals_alg + shift.
        # The Gibbs state lives in the eigenbasis (shared between H_phys and
        # H_alg — the global shift commutes and the rescale is positive), and
        # in eigenbasis coordinates it reads
        #   ρ_phys(β_phys) ∝ diag(exp(-β_phys · eigvals_phys))
        # which is mathematically identical to
        #   ρ_alg(β_alg = β_phys · rescaling_factor)
        # because the `exp(-β_phys · shift · I)` factor cancels in the
        # normalisation.  The test below verifies this numerically and then
        # asserts that a Lindbladian assembled with `HamHam(raw; beta_phys=…)`
        #   (i)  has ρ_phys as a fixed point (both `construct_lindbladian`
        #        and the matrix-free `predict_lindbladian_trajectory` path),
        #   (ii) drives an arbitrary initial state toward ρ_phys.

        ham_alg = _BPC_HAM_ALG10
        rescale = ham_alg.rescaling_factor
        shift   = ham_alg.shift
        β_phys  = 1.0
        β_alg   = β_phys * rescale

        # ρ_phys directly from the un-rescaled Hamiltonian (in eigenbasis).
        eigvals_phys = rescale .* ham_alg.eigvals .+ shift
        w_phys = exp.(-β_phys .* eigvals_phys)
        w_phys ./= sum(w_phys)
        ρ_phys_eigen = Diagonal(w_phys)

        # Sanity: ρ_phys(β_phys) ≡ ρ_alg(β_alg) as a density matrix
        # (the analytical statement at the top of the testset).
        ρ_alg_eigen = gibbs_state_in_eigen(ham_alg, β_alg)
        @test isapprox(Matrix(ρ_phys_eigen), Matrix(ρ_alg_eigen); atol=1e-12)

        # Rebuild a HamHam from the legacy fixture via the β_phys keyword
        # constructor.  Stored .gibbs must be the physical Gibbs state in the
        # eigenbasis (since ρ_phys == ρ_alg as a density matrix).
        raw_nt = (
            matrix             = Matrix{ComplexF64}(ham_alg.data),
            terms              = ham_alg.base_terms,
            base_coeffs        = ham_alg.base_coeffs,
            disordering_terms  = ham_alg.disordering_terms,
            disordering_coeffs = ham_alg.disordering_coeffs,
            eigvals            = ham_alg.eigvals,
            eigvecs            = ham_alg.eigvecs,
            nu_min             = ham_alg.nu_min,
            shift              = ham_alg.shift,
            rescaling_factor   = ham_alg.rescaling_factor,
            periodic           = ham_alg.periodic,
        )
        ham_phys = HamHam(raw_nt; beta_phys = β_phys)
        @test isapprox(Matrix(ham_phys.gibbs), Matrix(ρ_phys_eigen); atol=1e-12)

        # BohrDomain CKG smooth-Metro Lindbladian, n=3.  BohrDomain is the
        # analytic construction (no Riemann-sum quadrature), so the assertions
        # below are limited by floating-point only — the right home for the
        # MATHEMATICAL fixed-point claim about ρ_phys.  An EnergyDomain
        # quadrature-precision spot check follows.
        cfg_b = Config(
            sim = Lindbladian(), domain = BohrDomain(), construction = KMS(),
            num_qubits = 3, with_linear_combination = true,
            beta = β_alg, beta_phys = β_phys,
            sigma = 1.0 / β_alg, a = 0.0, s = 0.25,
        )
        @test validate_config!(cfg_b, ham_phys) === nothing

        num_qubits = 3
        jump_norm  = sqrt(3 * num_qubits)
        jumps = JumpOp[]
        for pauli in (X, Y, Z), site in 1:num_qubits
            op    = Matrix(pad_term([pauli], num_qubits, site)) ./ jump_norm
            op_eb = ham_phys.eigvecs' * op * ham_phys.eigvecs
            push!(jumps, JumpOp(op, op_eb, op == transpose(op), op == op'))
        end

        # (i.a) dense `construct_lindbladian` fixes ρ_phys (eigenbasis).
        L_super = construct_lindbladian(jumps, cfg_b, ham_phys)
        ρ_phys_mat = Matrix{ComplexF64}(ρ_phys_eigen)
        residual_dense = norm(L_super * vec(ρ_phys_mat)) / norm(vec(ρ_phys_mat))
        @test residual_dense < 1e-10

        # (i.b) matrix-free `apply_lindbladian!` agrees with the dense path —
        # the same fact via the Krylov code-path used by the trajectory
        # predictor and the gap solver.
        ws = Workspace(cfg_b, ham_phys, jumps)
        apply_lindbladian!(ws, ρ_phys_mat, cfg_b, ham_phys)
        residual_matvec = norm(ws.scratch.rho_out) / norm(vec(ρ_phys_mat))
        @test residual_matvec < 1e-10

        # (ii) `predict_lindbladian_trajectory` from a generic initial state
        # converges to ρ_phys.  Horizon t=200 is generous at n=3, β_alg≈14 —
        # τ_mix ≈ 30 here (gap ~ 0.1) so 200/30 ≈ 6 e-folds is enough for the
        # captured eigenmode to dominate; the captured `rho_inf` is exact on
        # the Krylov subspace regardless of horizon.
        d = size(ham_phys.data, 1)
        rho_0 = Matrix{ComplexF64}(I(d) ./ d)
        t_grid = collect(range(0.0, 200.0, length=21))
        res = predict_lindbladian_trajectory(cfg_b, ham_phys, jumps, rho_0, t_grid;
                                              krylovdim = 30, tol = 1e-10)
        @test res.all_converged
        # The trajectory's `sigma_beta` reference must match the physical
        # Gibbs state in eigenbasis.
        @test isapprox(Matrix(res.sigma_beta), Matrix(ρ_phys_eigen); atol=1e-12)
        # The captured ρ_∞ from the Krylov eigendecomposition must agree with
        # ρ_phys.  BohrDomain is analytic so the residual is FP-only, but the
        # KrylovKit Arnoldi `tol = 1e-10` setting controls the leading-mode
        # null-space accuracy — atol≈5e-9 is the realistic floor here.
        @test isapprox(Matrix(res.rho_inf), Matrix(ρ_phys_eigen); atol=5e-9)
        # Distance at horizon is gap-limited; just assert it is shrinking.
        @test res.distances[end] < res.distances[1]
        @test res.distances[end] < 1e-3

        # EnergyDomain quadrature-precision spot check — the dense Lindbladian
        # built via the simulator-realistic Riemann-sum construction (which
        # production / Krylov sweeps actually consume) still fixes ρ_phys, but
        # only to the EnergyDomain quadrature precision at r_D=12 (~ 1e-5 at
        # this fixture, per qf-7xt).
        cfg_e = Config(
            sim = Lindbladian(), domain = EnergyDomain(), construction = KMS(),
            num_qubits = 3, with_linear_combination = true,
            beta = β_alg, beta_phys = β_phys,
            sigma = 1.0 / β_alg, a = 0.0, s = 0.25,
            num_energy_bits = 12, w0 = 0.05,
            t0 = 2π / (2^12 * 0.05),
            num_trotter_steps_per_t0 = 10,
        )
        L_super_e = construct_lindbladian(jumps, cfg_e, ham_phys)
        residual_energy = norm(L_super_e * vec(ρ_phys_mat)) / norm(vec(ρ_phys_mat))
        @test residual_energy < 1e-4

        @info "(f) β_phys fixed-point" rescale β_phys β_alg residual_dense residual_matvec residual_energy dist_end=res.distances[end]
    end

    @testset "(e) validate_config!(cfg, ham) consistency" begin
        rescale = _BPC_HAM_ALG10.rescaling_factor
        β_phys_val = 1.0
        β_alg_val = β_phys_val * rescale

        # Consistent — must succeed
        cfg_ok = Config(
            sim = Lindbladian(), domain = BohrDomain(), construction = KMS(),
            num_qubits = 3, with_linear_combination = true,
            beta = β_alg_val, beta_phys = β_phys_val,
            sigma = 1.0 / β_alg_val, a = 0.0, s = 0.25,
        )
        @test validate_config!(cfg_ok, _BPC_HAM_ALG10) === nothing

        # Inconsistent: β_alg differs from β_phys * rescale by > tolerance
        cfg_bad = Config(
            sim = Lindbladian(), domain = BohrDomain(), construction = KMS(),
            num_qubits = 3, with_linear_combination = true,
            beta = β_alg_val * 1.1,            # 10% off
            beta_phys = β_phys_val,
            sigma = 1.0 / β_alg_val, a = 0.0, s = 0.25,
        )
        @test_throws ArgumentError validate_config!(cfg_bad, _BPC_HAM_ALG10)

        # β_phys = nothing — skip the check (legacy-compatible path)
        cfg_legacy = Config(
            sim = Lindbladian(), domain = BohrDomain(), construction = KMS(),
            num_qubits = 3, with_linear_combination = true,
            beta = β_alg_val,
            sigma = 1.0 / β_alg_val, a = 0.0, s = 0.25,
        )
        @test validate_config!(cfg_legacy, _BPC_HAM_ALG10) === nothing
    end
end
