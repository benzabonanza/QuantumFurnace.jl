@testset "DLL multi-channel filter (qf-7go.1)" begin
    # -----------------------------------------------------------------------
    # (a) Constructor: k = 1 stores the channel verbatim, β-mismatch throws,
    #     empty channel vector throws, channel without `beta` throws.
    # -----------------------------------------------------------------------
    @testset "(a) constructor sanity" begin
        β = 5.0
        ch = [DLLMetropolisFilter(β; S = 2.0)]
        multi = DLLMultiChannelFilter(ch, β)
        @test multi.beta === β
        @test length(multi.channels) == 1
        @test multi.channels[1] === ch[1]
        # element-type ⇒ Complex{T} regardless of channel mix.
        @test eltype(multi) === ComplexF64

        # Two-channel mix of homogeneous Metropolis filters at matching β.
        ch2 = [DLLMetropolisFilter(β; S = 2.0), DLLMetropolisFilter(β; S = 2.0)]
        multi2 = DLLMultiChannelFilter(ch2, β)
        @test length(multi2.channels) == 2

        # Empty channels → ArgumentError.
        @test_throws ArgumentError DLLMultiChannelFilter(
            DLLMetropolisFilter{Float64}[], β)

        # Mismatched β → ArgumentError.
        @test_throws ArgumentError DLLMultiChannelFilter(
            [DLLMetropolisFilter(β; S = 2.0), DLLMetropolisFilter(β + 1; S = 2.0)],
            β,
        )

        # Channel without `.beta` field → ArgumentError.
        @test_throws ArgumentError DLLMultiChannelFilter(
            [GaussianFilter(0.5)], β)
    end

    # -----------------------------------------------------------------------
    # (b) Multi-channel q_weight / freq_kernel are literal sums over channels.
    # -----------------------------------------------------------------------
    @testset "(b) q_weight + freq_kernel = sum over channels" begin
        β = 2.0
        ch = [DLLMetropolisFilter(β; S = 2.0), DLLGaussianFilter(β)]
        # NB: this constructor path uses the inner ctor with F = AbstractFilter.
        multi = DLLMultiChannelFilter{Float64, AbstractFilter}(
            AbstractFilter[ch[1], ch[2]], β)

        for nu in (-0.7, -0.2, 0.0, 0.4, 1.1)
            expected_q  = sum(QuantumFurnace.q_weight(c, nu) for c in ch)
            expected_fk = sum(freq_kernel(c, nu) for c in ch)
            @test isapprox(QuantumFurnace.q_weight(multi, nu), expected_q;
                           atol = 1e-15)
            @test isapprox(freq_kernel(multi, nu), expected_fk; atol = 1e-15)
        end
    end

    # -----------------------------------------------------------------------
    # (c) Multi-channel time_kernel = sum of per-channel time_kernel
    #     (complex-valued; verify on both real and imaginary parts).
    # -----------------------------------------------------------------------
    @testset "(c) time_kernel = sum over channels" begin
        β = 1.5
        ch = [DLLGaussianFilter(β), DLLGaussianFilter(β)]  # k = 2 same channel
        multi = DLLMultiChannelFilter(ch, β)

        for t in (-3.0, -0.7, 0.0, 0.4, 2.1)
            ref = ComplexF64(0)
            for c in ch
                ref += time_kernel(c, t)
            end
            @test isapprox(time_kernel(multi, t), ref; atol = 1e-13)
            # Two identical channels ⇒ exactly 2× single-channel value.
            @test isapprox(time_kernel(multi, t), 2 * time_kernel(ch[1], t);
                           atol = 1e-13)
        end
    end

    # -----------------------------------------------------------------------
    # (d) k = 1 reduces to the single-channel filter at machine precision
    #     for q_weight / freq_kernel / time_kernel.
    # -----------------------------------------------------------------------
    @testset "(d) k = 1 reduces to single-channel filter" begin
        β = 5.0
        for ch in (DLLGaussianFilter(β), DLLMetropolisFilter(β; S = 2.0))
            multi = DLLMultiChannelFilter([ch], β)
            for nu in (-0.5, 0.0, 0.3, 0.9)
                @test QuantumFurnace.q_weight(multi, nu) ==
                      QuantumFurnace.q_weight(ch, nu)
                @test freq_kernel(multi, nu) == freq_kernel(ch, nu)
            end
            for t in (-2.0, 0.0, 1.0)
                @test time_kernel(multi, t) == time_kernel(ch, t)
            end
        end
    end

    # -----------------------------------------------------------------------
    # (e) filter_time_cutoff uses tol/k per channel and returns the max.
    #     The returned cutoff must satisfy |time_kernel(multi, tc)| ≤ tol
    #     (sum of per-channel residuals).
    # -----------------------------------------------------------------------
    @testset "(e) filter_time_cutoff bounds |time_kernel(multi, tc)| ≤ tol" begin
        β = 2.0
        ch = [DLLGaussianFilter(β), DLLGaussianFilter(β)]
        multi = DLLMultiChannelFilter(ch, β)
        # `slack = 4` matches existing DLLGaussianFilter cutoff tests
        # (test_dll_filter.jl:103); absorbs ULP loss in log/sqrt for the
        # closed-form Gaussian cutoff at the boundary.
        slack = 4.0
        for tol in (1e-6, 1e-9, 1e-12)
            tc = filter_time_cutoff(multi, tol)
            # Triangle inequality: each channel ≤ tol/k, so the sum ≤ tol (× slack).
            @test abs(time_kernel(multi, tc)) <= slack * tol
            # Equals max over per-channel cutoffs at the per-tol budget.
            expected_tc = maximum(filter_time_cutoff(c, tol / length(ch)) for c in ch)
            @test tc ≈ expected_tc
            @test 0 < tc < 1e6
        end
    end

    # -----------------------------------------------------------------------
    # (f) Type stability: Float64 channels ⇒ Float64-typed cutoff and
    #     ComplexF64 time kernel.
    # -----------------------------------------------------------------------
    @testset "(f) type stability" begin
        β = 2.0
        ch = [DLLMetropolisFilter(β; S = 2.0), DLLMetropolisFilter(β; S = 2.0)]
        multi = DLLMultiChannelFilter(ch, β)
        @test typeof(time_kernel(multi, 1.0)) === ComplexF64
        @test typeof(QuantumFurnace.q_weight(multi, 0.3)) === Float64
        @test typeof(freq_kernel(multi, 0.3)) === Float64
        @test typeof(filter_time_cutoff(multi, 1e-9)) === Float64
    end

    # =====================================================================
    # ShiftedSymmetricFilter + dll_multichannel_translates (qf-7go.5)
    # =====================================================================
    @testset "(g) ShiftedSymmetricFilter shift=0 reduces to scaled base" begin
        β = 5.0
        for base in (DLLGaussianFilter(β), DLLMetropolisFilter(β; S = 2.0))
            for w in (1.0, 0.4, 2.5)
                s = ShiftedSymmetricFilter(base, 0.0, w)
                for nu in (-0.5, 0.0, 0.3, 0.9)
                    @test isapprox(QuantumFurnace.q_weight(s, nu),
                                   sqrt(w) * QuantumFurnace.q_weight(base, nu);
                                   atol = 1e-15)
                    @test isapprox(freq_kernel(s, nu),
                                   sqrt(w) * freq_kernel(base, nu);
                                   atol = 1e-15)
                end
                for t in (-2.0, 0.0, 1.5)
                    @test isapprox(time_kernel(s, t),
                                   sqrt(w) * time_kernel(base, t);
                                   atol = 1e-13)
                end
            end
        end
    end

    @testset "(h) ShiftedSymmetricFilter shift≠0: q is real-even" begin
        β = 5.0
        for base in (DLLGaussianFilter(β), DLLMetropolisFilter(β; S = 2.0))
            for sh in (0.2, 0.5, 0.9)
                # base.S/2 = 1; centers must be ≤ 1.
                base isa DLLMetropolisFilter && abs(sh) > base.S/2 && continue
                s = ShiftedSymmetricFilter(base, sh, 1.0)
                for nu in (0.1, 0.3, 0.7, 1.4)
                    @test isapprox(QuantumFurnace.q_weight(s, nu),
                                   QuantumFurnace.q_weight(s, -nu);
                                   atol = 1e-15)
                end
            end
        end
    end

    @testset "(i) ShiftedSymmetricFilter shift≠0: closed-form q_l" begin
        β = 5.0
        base = DLLMetropolisFilter(β; S = 2.0)
        sh = 0.5
        w = 1.0
        s = ShiftedSymmetricFilter(base, sh, w)
        for nu in (-1.2, -0.5, 0.0, 0.3, 0.9, 1.4)
            expected = sqrt(w/2) * (
                QuantumFurnace.q_weight(base, nu - sh) +
                QuantumFurnace.q_weight(base, nu + sh))
            @test isapprox(QuantumFurnace.q_weight(s, nu), expected;
                           atol = 1e-15)
        end
    end

    @testset "(j) ShiftedSymmetricFilter shift≠0: time_kernel formula" begin
        # f_l(t) = √(w/2) · f_base(t) · 2 cosh(β·shift/4 + i·shift·t)
        β = 2.0
        base = DLLGaussianFilter(β)
        sh = 0.3
        w = 1.0
        s = ShiftedSymmetricFilter(base, sh, w)
        for t in (-2.0, -0.5, 0.0, 0.5, 2.0)
            z = complex(β * sh / 4, sh * t)
            expected = sqrt(w/2) * time_kernel(base, t) * (2 * cosh(z))
            @test isapprox(time_kernel(s, t), expected; atol = 1e-13)
        end
    end

    @testset "(k) dll_multichannel_translates: k=1 center=0 ≈ base × √w" begin
        β = 5.0
        base = DLLMetropolisFilter(β; S = 2.0)
        # k=1, center=[0.0], weight=[1.0]: each method matches base byte-for-byte.
        multi1 = dll_multichannel_translates(base; centers = [0.0])
        @test length(multi1.channels) == 1
        ch = multi1.channels[1]
        @test ch.shift == 0.0 && ch.weight == 1.0 && ch.base === base
        # multi-channel diagnostic methods evaluate to base's values.
        for nu in (-0.5, 0.0, 0.3, 0.9)
            @test isapprox(QuantumFurnace.q_weight(multi1, nu),
                           QuantumFurnace.q_weight(base, nu); atol = 1e-15)
        end
    end

    @testset "(l) dll_multichannel_translates: KMS skew-symmetry, k=2" begin
        for β in (1.0, 5.0, 10.0)
            base = DLLMetropolisFilter(β; S = 2.0)
            multi = dll_multichannel_translates(base; centers = [0.0, 0.4])
            # Sum of per-channel rank-1 α^(ℓ); each is KMS-skew-symmetric ⇒ sum is too.
            ν_grid = collect(range(-0.5, 0.5; length = 11))
            α = sum(
                QuantumFurnace.dll_kossakowski_bohr(c, ν_grid)
                for c in multi.channels
            )
            assert_kms_skew_symmetric(α, ν_grid, β; atol = 1e-12)
        end
    end

    @testset "(m) dll_multichannel_translates rejects out-of-bump centers" begin
        β = 5.0
        base = DLLMetropolisFilter(β; S = 2.0)
        # centers must satisfy |center| ≤ S/2 = 1.0.
        @test_throws ArgumentError dll_multichannel_translates(base; centers = [1.5])
        @test_throws ArgumentError dll_multichannel_translates(base; centers = [0.0, 1.001])
        # In-range succeeds.
        @test dll_multichannel_translates(base; centers = [0.0, 1.0]) isa DLLMultiChannelFilter
        # Empty centers throws.
        @test_throws ArgumentError dll_multichannel_translates(base; centers = Float64[])
        # Mismatched weights length throws.
        @test_throws ArgumentError dll_multichannel_translates(base;
                                                                centers = [0.0, 0.5],
                                                                weights = [1.0])
        # Non-positive weights throw.
        @test_throws ArgumentError dll_multichannel_translates(base;
                                                                centers = [0.0],
                                                                weights = [-0.1])
        # Non-DLL base (no .beta) throws.
        @test_throws ArgumentError dll_multichannel_translates(GaussianFilter(0.5))
    end
end
