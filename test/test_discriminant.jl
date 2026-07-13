using Random: MersenneTwister
using LinearAlgebra: mul!

@testset "Discriminant (Phase 1: primitives)" begin

    # -----------------------------------------------------------------------
    # DiscriminantBuffers
    # -----------------------------------------------------------------------
    @testset "DiscriminantBuffers" begin
        bufs = DiscriminantBuffers(N3_DIM)

        @test bufs isa DiscriminantBuffers{ComplexF64}
        @test size(bufs.work1) == (N3_DIM, N3_DIM)
        @test size(bufs.work2) == (N3_DIM, N3_DIM)
        @test size(bufs.work3) == (N3_DIM, N3_DIM)

        # Parametric construction
        bufs32 = DiscriminantBuffers{ComplexF32}(4)
        @test bufs32 isa DiscriminantBuffers{ComplexF32}
        @test eltype(bufs32.work1) === ComplexF32
    end

    # -----------------------------------------------------------------------
    # gibbs_fractional_powers: identities
    # -----------------------------------------------------------------------
    @testset "gibbs_fractional_powers identities (N3 fixture)" begin
        powers = gibbs_fractional_powers(N3_GIBBS)

        @test powers isa NamedTuple
        @test propertynames(powers) == (:sigma_quarter, :sigma_inv_quarter, :sigma_half)
        @test length(powers.sigma_quarter) == N3_DIM
        @test length(powers.sigma_inv_quarter) == N3_DIM
        @test length(powers.sigma_half) == N3_DIM

        diag_real = real.(diag(N3_GIBBS))

        # sigma_quarter^4 = sigma  (after eps_trunc floor; for N3 fixture
        # the floor is not active, so equality should be exact-up-to-fp).
        @test isapprox(powers.sigma_quarter .^ 4, diag_real; atol=TOL_EXACT)

        # sigma_quarter * sigma_inv_quarter = 1
        @test isapprox(powers.sigma_quarter .* powers.sigma_inv_quarter,
                       ones(Float64, N3_DIM); atol=TOL_EXACT)

        # sigma_half = sigma_quarter^2
        @test isapprox(powers.sigma_half, powers.sigma_quarter .^ 2; atol=TOL_EXACT)

        # sigma_half^2 = sigma
        @test isapprox(powers.sigma_half .^ 2, diag_real; atol=TOL_EXACT)

        # All powers strictly positive (Gibbs has no zero eigenvalues at finite beta)
        @test all(powers.sigma_quarter .> 0)
        @test all(powers.sigma_inv_quarter .> 0)
        @test all(powers.sigma_half .> 0)
    end

    # -----------------------------------------------------------------------
    # apply_discriminant!: closure-style action against vec-form reference
    # -----------------------------------------------------------------------
    @testset "apply_discriminant! matches vec form (N3 KMS Bohr)" begin
        # Dense 3-qubit KMS-DB Lindbladian.
        config = make_config(Lindbladian(), BohrDomain(); num_qubits=3, construction=KMS())
        L_sparse = construct_lindbladian(N3_JUMPS, config, N3_HAM)
        L_dense  = Matrix{ComplexF64}(L_sparse)

        powers = gibbs_fractional_powers(N3_GIBBS)
        sq, sq_inv = powers.sigma_quarter, powers.sigma_inv_quarter

        # Reference materialised discriminant (column-stacking convention,
        # mirrors src/diagnostics.jl:236-249):
        #     D = (sigma^{-1/4} kron sigma^{-1/4}) * L * (sigma^{1/4} kron sigma^{1/4})
        d_left  = kron(sq_inv, sq_inv)
        d_right = kron(sq, sq)
        D_ref   = d_left .* L_dense .* d_right'

        # Closure-style action wrapping the dense Lindbladian.
        function lindblad_action!(out_mat, in_mat)
            mul!(vec(out_mat), L_dense, vec(in_mat))
            return out_mat
        end

        bufs = DiscriminantBuffers(N3_DIM)
        out  = Matrix{ComplexF64}(undef, N3_DIM, N3_DIM)

        # Random Hermitian + non-Hermitian inputs
        rng = MersenneTwister(0xC0FFEE)
        for trial in 1:5
            X = randn(rng, ComplexF64, N3_DIM, N3_DIM)
            apply_discriminant!(out, X, lindblad_action!, sq, sq_inv, bufs)
            ref = reshape(D_ref * vec(X), N3_DIM, N3_DIM)
            @test isapprox(out, ref; atol=TOL_EXACT, rtol=TOL_EXACT)
        end

        # Hermitian input
        X_h = Matrix(random_density_matrix(3))
        apply_discriminant!(out, X_h, lindblad_action!, sq, sq_inv, bufs)
        ref_h = reshape(D_ref * vec(X_h), N3_DIM, N3_DIM)
        @test isapprox(out, ref_h; atol=TOL_EXACT, rtol=TOL_EXACT)
    end

    # -----------------------------------------------------------------------
    # apply_discriminant!: zero allocations in the body itself
    # -----------------------------------------------------------------------
    @testset "apply_discriminant! body is allocation-free" begin
        # No-op Lindbladian action (just copy in -> out).  Proves any
        # allocation seen with a real Lindbladian closure comes from the
        # user closure, not from apply_discriminant! itself.
        noop_action!(o, x) = (copyto!(o, x); o)

        powers = gibbs_fractional_powers(N3_GIBBS)
        sq, sq_inv = powers.sigma_quarter, powers.sigma_inv_quarter

        bufs = DiscriminantBuffers(N3_DIM)
        out  = Matrix{ComplexF64}(undef, N3_DIM, N3_DIM)
        X    = Matrix(random_density_matrix(3))

        # Warmup -- compile method specialisations.
        apply_discriminant!(out, X, noop_action!, sq, sq_inv, bufs)

        allocs = @allocated apply_discriminant!(out, X, noop_action!, sq, sq_inv, bufs)
        @test allocs == 0
        @info "apply_discriminant! body allocations (noop closure)" allocs
    end

    # -----------------------------------------------------------------------
    # materialize_discriminant: kernel matches the broadcast/kron formula
    # -----------------------------------------------------------------------
    @testset "materialize_discriminant matches broadcast/kron form" begin
        config = make_config(Lindbladian(), BohrDomain(); num_qubits=3, construction=KMS())
        L_dense = Matrix{ComplexF64}(construct_lindbladian(N3_JUMPS, config, N3_HAM))

        powers = gibbs_fractional_powers(N3_GIBBS)
        sq, sq_inv = powers.sigma_quarter, powers.sigma_inv_quarter

        # Reference via the explicit broadcast / kron form (the same
        # formula compute_anti_hermitian_defect used pre-refactor).
        d_left  = kron(sq_inv, sq_inv)
        d_right = kron(sq, sq)
        D_ref   = d_left .* L_dense .* d_right'

        D_new = materialize_discriminant(L_dense, N3_GIBBS)
        @test isapprox(D_new, D_ref; atol=TOL_EXACT, rtol=TOL_EXACT)

        # In-place form writes into the supplied buffer.
        D_buf = similar(L_dense)
        materialize_discriminant!(D_buf, L_dense, N3_GIBBS)
        @test D_buf === materialize_discriminant!(D_buf, L_dense, N3_GIBBS)
        @test isapprox(D_buf, D_ref; atol=TOL_EXACT, rtol=TOL_EXACT)
    end

    # -----------------------------------------------------------------------
    # Spectrum invariance: D and L are similarity-equivalent (same spectrum)
    # -----------------------------------------------------------------------
    @testset "materialize_discriminant preserves spectrum (N3 KMS Bohr)" begin
        config = make_config(Lindbladian(), BohrDomain(); num_qubits=3, construction=KMS())
        L_dense = Matrix{ComplexF64}(construct_lindbladian(N3_JUMPS, config, N3_HAM))

        D = materialize_discriminant(L_dense, N3_GIBBS)

        spec_L = sort(eigvals(L_dense); by = z -> (real(z), imag(z)))
        spec_D = sort(eigvals(D);       by = z -> (real(z), imag(z)))

        @test maximum(abs.(spec_L .- spec_D)) < 1e-10
        @info "spectrum invariance" max_eig_err=maximum(abs.(spec_L .- spec_D))
    end

    # -----------------------------------------------------------------------
    # hermitian_antihermitian_split: structural identities
    # -----------------------------------------------------------------------
    @testset "hermitian_antihermitian_split identities" begin
        rng = MersenneTwister(0xBADCAFE)
        d2 = N3_DIM^2
        D = randn(rng, ComplexF64, d2, d2)

        H, A = hermitian_antihermitian_split(D)

        # Reconstruction
        @test isapprox(H + A, D; atol=TOL_EXACT, rtol=TOL_EXACT)

        # Hermiticity / anti-Hermiticity
        @test isapprox(H, H'; atol=TOL_EXACT, rtol=TOL_EXACT)
        @test isapprox(A, -A'; atol=TOL_EXACT, rtol=TOL_EXACT)

        # In-place version produces the same output.
        H_buf = similar(D)
        A_buf = similar(D)
        hermitian_antihermitian_split!(H_buf, A_buf, D)
        @test isapprox(H_buf, H; atol=TOL_EXACT, rtol=TOL_EXACT)
        @test isapprox(A_buf, A; atol=TOL_EXACT, rtol=TOL_EXACT)
    end

    # -----------------------------------------------------------------------
    # Refactor regression: compute_anti_hermitian_defect still produces a
    # consistent DefectResult on the N3 KMS Bohr fixture.
    # -----------------------------------------------------------------------
    @testset "compute_anti_hermitian_defect regression after refactor" begin
        config = make_config(Lindbladian(), BohrDomain(); num_qubits=3, construction=KMS())
        L_dense = Matrix{ComplexF64}(construct_lindbladian(N3_JUMPS, config, N3_HAM))

        # Compute it the OLD way inline (pre-refactor algorithm) and compare
        # to the new delegated implementation.  Any change in numerical
        # values would indicate the refactor altered behaviour.
        gibbs_diag      = real.(diag(Matrix(N3_GIBBS)))
        gibbs_diag_safe = max.(gibbs_diag, 1e-12)
        rq      = gibbs_diag_safe .^ 0.25
        rq_inv  = gibbs_diag_safe .^ (-0.25)
        D_old   = kron(rq_inv, rq_inv) .* L_dense .* kron(rq, rq)'
        H_old   = (D_old + D_old') / 2
        A_old   = (D_old - D_old') / 2

        A_norm_old = opnorm(A_old)
        H_eigs_old = sort(abs.(eigvals(Hermitian(H_old))))
        H_gap_old  = H_eigs_old[2]
        ratio_old  = A_norm_old / max(H_gap_old, 1e-30)

        result = compute_anti_hermitian_defect(L_dense, N3_GIBBS)

        @test isapprox(result.A_norm,       A_norm_old; atol=TOL_EXACT, rtol=TOL_EXACT)
        @test isapprox(result.H_gap,        H_gap_old;  atol=TOL_EXACT, rtol=TOL_EXACT)
        @test isapprox(result.defect_ratio, ratio_old;  atol=TOL_EXACT, rtol=TOL_EXACT)
    end

    # -----------------------------------------------------------------------
    # discriminant_spectrum: H-part eigenvalues match L spectrum for KMS-DB
    # -----------------------------------------------------------------------
    @testset "discriminant_spectrum on KMS-DB Lindbladian (N3 Bohr)" begin
        config = make_config(Lindbladian(), BohrDomain(); num_qubits=3, construction=KMS())
        L_dense = Matrix{ComplexF64}(construct_lindbladian(N3_JUMPS, config, N3_HAM))

        spec = discriminant_spectrum(L_dense, N3_GIBBS; n_modes=10)

        @test spec isa DiscriminantSpectrum
        @test spec.n_modes == 10
        @test length(spec.H_eigenvalues) == 10

        # Sorted ascending by |λ|
        @test issorted(abs.(spec.H_eigenvalues))

        # First eigenvalue ≈ 0 (steady-state purification)
        @test abs(spec.H_eigenvalues[1]) < 1e-10

        # H_gap matches |H_eigenvalues[2]|
        @test spec.H_gap ≈ abs(spec.H_eigenvalues[2])

        # For KMS-DB: H_part = D, so eigvals(H_part) ⊆ eigvals(L) ∈ ℝ.
        # Compare the leading H eigenvalues to the leading real(L) eigenvalues.
        L_eigs = eigvals(L_dense)
        L_real = sort(real.(L_eigs); by = abs)[1:10]
        # Sort both by |λ| and compare ascending.
        sorted_H = sort(spec.H_eigenvalues; by = abs)
        sorted_L = sort(L_real; by = abs)
        @test maximum(abs.(sorted_H .- sorted_L)) < 1e-10
        @info "discriminant_spectrum vs L spectrum (KMS-DB)" max_err=maximum(abs.(sorted_H .- sorted_L))
    end

    # -----------------------------------------------------------------------
    # Phase 4 DBV physics diagnostics: with vs without coherent term
    # -----------------------------------------------------------------------
    @testset "verify_detailed_balance physics diagnostics (n=3,4,5 KMS Bohr)" begin
        for n in (3, 4, 5)
            @testset "n=$n" begin
                sys   = make_test_system(; num_qubits=n)
                ham   = sys.hamiltonian
                jumps = sys.jumps
                gibbs = sys.gibbs

                config = make_config(Lindbladian(), BohrDomain(); num_qubits=n, construction=KMS())
                L_with = Matrix{ComplexF64}(construct_lindbladian(jumps, config, ham))
                L_no   = Matrix{ComplexF64}(construct_lindbladian(jumps, config, ham; include_coherent=false))

                # === DBV-A: with coherent → KMS-DB (Hermitian discriminant) ===
                res_with = verify_detailed_balance(L_with, gibbs)

                @test res_with.relative_norm < 1e-10
                @test res_with.fixed_point_residual < 1e-10
                @test res_with.is_kms_db
                @info "DBV-A n=$n with coherent" rel=res_with.relative_norm fp=res_with.fixed_point_residual

                # === DBV-B: parent-Hamiltonian gap == L spectral gap ===
                # For KMS-DB L, eigvals(D) = eigvals(L) ⊆ ℝ with one near-zero
                # eigenvalue (steady state) and the rest negative.  The
                # Hermitian-part gap of D therefore equals the spectral gap
                # of L (Chen et al. 2025 Thm I.3 framing).
                @test isapprox(res_with.hermitian_part_gap,
                               res_with.spectral_gap_L; rtol=1e-9)
                @info "DBV-B n=$n H-part gap = L gap" H_gap=res_with.hermitian_part_gap L_gap=res_with.spectral_gap_L

                # Cross-check via discriminant_spectrum (separate code path).
                spec = discriminant_spectrum(L_with, gibbs; n_modes=4)
                @test isapprox(spec.H_gap, res_with.hermitian_part_gap; rtol=1e-9)

                # === DBV-C: without coherent → measurably nonzero ||A|| ===
                res_no = verify_detailed_balance(L_no, gibbs)

                # Fixture-robust: dropping the Lamb-shift makes the DB residual
                # orders of magnitude larger than the with-coherent residual (which
                # is ~machine precision). Absolute magnitude is fixture-dependent
                # (seed-46: 1.2e-5 at n=3 .. 7e-4 at n=4), so compare to res_with.
                @test res_no.relative_norm > 100 * res_with.relative_norm
                @test !res_no.is_kms_db
                # Without the Lamb-shift coherent term, σ is no longer the
                # exact steady state of L (the dissipator alone satisfies
                # `D(σ) = i[B, σ]`, where B is the Lamb-shift), so the
                # fixed-point residual is also measurably nonzero -- a
                # second independent signal of KMS-DB violation.
                @test res_no.fixed_point_residual > 1e-6
                @info "DBV-C n=$n without coherent" rel=res_no.relative_norm fp=res_no.fixed_point_residual
            end
        end
    end

    # -----------------------------------------------------------------------
    # eps_trunc floor: defensive behaviour at synthetic low temperature
    # -----------------------------------------------------------------------
    @testset "gibbs_fractional_powers eps_trunc floor" begin
        # Synthetic 4x4 diagonal Gibbs with one entry below eps_trunc=1e-12.
        diag_vals = [0.5, 0.49, 1e-20, 0.01]
        rho = Hermitian(diagm(ComplexF64.(diag_vals)))

        # Default eps_trunc = 1e-12: third entry is floored.
        powers = gibbs_fractional_powers(rho)

        @test powers.sigma_inv_quarter[3] ≈ (1e-12)^(-0.25)
        @test isfinite(powers.sigma_inv_quarter[3])

        # Other entries untouched.
        for i in (1, 2, 4)
            @test powers.sigma_quarter[i] ≈ diag_vals[i]^0.25
            @test powers.sigma_inv_quarter[i] ≈ diag_vals[i]^(-0.25)
        end

        # Lower eps_trunc lets the small entry through.
        powers_loose = gibbs_fractional_powers(rho; eps_trunc=1e-30)
        @test powers_loose.sigma_inv_quarter[3] ≈ (1e-20)^(-0.25)
        @test powers_loose.sigma_inv_quarter[3] > powers.sigma_inv_quarter[3]
    end
end
