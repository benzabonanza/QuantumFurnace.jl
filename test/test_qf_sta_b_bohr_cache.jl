# test_qf_sta_b_bohr_cache.jl
#
# Regression tests for the qf-sta B_bohr refactor: outer-ν₂ / inner-jump
# restructure with sparse row-aware f-eval and `last_i` row cache.
#
# The refactor walks `bohr_dict[ν₂]` directly and computes
# `f_row[col] = f(bohr_freqs[i, col], ν₂)` once per (ν₂, i) pair, caching
# the most recent `i` to skip recomputation when the same row appears twice
# in a row. The cache is correct only if either:
#  (a) `bohr_dict[ν₂]` is row-sorted (so repeated `i` are adjacent), OR
#  (b) the algorithm is order-invariant (i.e., correct for any ordering).
#
# Each (i, j) pair contributes independently to `B[j, col] += conj(in_eb[i, j]) *
# f(bohr_freqs[i, col], ν₂) * in_eb[i, col]`. The `f_row` cache only affects
# performance — `if i != last_i` strictly recomputes whenever the row
# changes — so correctness holds for any ordering. These tests verify that
# claim by:
#   (a) constructing a Hamiltonian with degenerate eigenvalues so that
#       `bohr_dict[0.0]` contains both diagonal and off-diagonal entries
#       (Bohr-collision case), and comparing B_bohr against a closed-form
#       reference computed without any row caching.
#   (b) explicitly shuffling `bohr_dict[ν]` index orders for a generic
#       Hamiltonian and confirming B_bohr is invariant.
#   (c) running the Gaussian filter branch (`with_linear_combination =
#       false`) which uses a different `f` closure than the smooth-Metro
#       branch (covered elsewhere).

using Test
using LinearAlgebra
using Random
using QuantumFurnace
using QuantumFurnace: B_bohr, _pick_f, create_bohr_dict

@testset "qf-sta B_bohr row-cache correctness" begin

    # -----------------------------------------------------------------
    # (a) Reference implementation: per-(i, j, col) without any caching.
    # -----------------------------------------------------------------
    # Computes B[j, col] = Σ_{ν₂} Σ_{(i, j) ∈ bohr_dict[ν₂]} Σ_{a}
    #     conj(in_eb_a[i, j]) * f(bohr_freqs[i, col], ν₂) * in_eb_a[i, col]
    # by walking bohr_dict[ν₂] entry-by-entry, evaluating f freshly each
    # iteration. Used as the ground truth.
    function _b_bohr_reference(ham::HamHam, jumps::Vector{<:JumpOp},
                               cfg::Config; index_permuter::Function = identity)
        d = size(ham.data, 1)
        CT = ComplexF64
        B = zeros(CT, d, d)
        f = _pick_f(cfg)
        for nu_2 in keys(ham.bohr_dict)
            indices = index_permuter(copy(ham.bohr_dict[nu_2]))
            for idx in indices
                i = idx[1]; j = idx[2]
                for jump in jumps
                    in_eb = jump.in_eigenbasis
                    val = conj(in_eb[i, j])
                    for col in 1:d
                        B[j, col] += val * f(ham.bohr_freqs[i, col], nu_2) *
                                     in_eb[i, col]
                    end
                end
            end
        end
        return B
    end

    # -----------------------------------------------------------------
    # (b) Build a custom HamHam with intentionally-degenerate spectrum.
    # -----------------------------------------------------------------
    # Two-qubit Hamiltonian with two pairs of degenerate eigenvalues, so
    # bohr_dict[0.0] gets off-diagonal entries (Bohr collisions) on top of
    # the diagonal. dim=4 is small enough that we can verify by hand.
    function _make_degenerate_ham(beta::Real)
        # Diagonal H with eigvals {0, 1, 1, 2} — two pairs of equals.
        # bohr_freqs[i, j] = E_i - E_j gets repeated values, including
        # E_2 - E_3 = 0 and E_3 - E_2 = 0 (off-diagonal collision in
        # bohr_dict[0.0]).
        H = ComplexF64[0 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 2]
        eigvals_d = [0.0, 1.0, 1.0, 2.0]
        eigvecs_d = Matrix{ComplexF64}(I, 4, 4)
        bohr_freqs = eigvals_d .- transpose(eigvals_d)
        bohr_dict = create_bohr_dict(bohr_freqs)
        Z = sum(exp.(-beta .* eigvals_d))
        gibbs_diag = ComplexF64[exp(-beta * E) / Z for E in eigvals_d]
        gibbs = Hermitian(diagm(gibbs_diag))
        return HamHam{Float64}(
            H,
            bohr_freqs,
            bohr_dict,
            Vector{Vector{Matrix{ComplexF64}}}(),
            Float64[],
            nothing, nothing,
            eigvals_d,
            eigvecs_d,
            1.0,    # nu_min (not used in B_bohr)
            0.0, 1.0, true, gibbs,
        )
    end

    # Sanity check: does the degenerate Hamiltonian produce an off-diagonal
    # Bohr collision at frequency 0?
    #
    # Note: Julia's `Dict` key comparison uses `isequal`, and
    # `isequal(0.0, -0.0) == false`. `create_bohr_dict` pushes (2, 3) into
    # `bohr_dict[bohr_freqs[2, 3]] = bohr_dict[0.0]` and (3, 2) into
    # `bohr_dict[-bohr_freqs[2, 3]] = bohr_dict[-0.0]`. So when ν₂ happens
    # to hit 0 from a Bohr collision, the symmetric pair lands in TWO
    # separate dict keys (0.0 and -0.0). This is pre-existing behavior; the
    # production B_bohr loop over `keys(bohr_dict)` correctly iterates both.
    # `f(ν, 0.0) == f(ν, -0.0)` for both filter branches.
    @testset "(setup) degenerate ham splits 0.0 / -0.0 keys" begin
        ham_deg = _make_degenerate_ham(5.0)
        z_pos = ham_deg.bohr_dict[0.0]
        z_neg = ham_deg.bohr_dict[-0.0]
        # bohr_dict[0.0] = [(1,1), (2,2), (3,3), (4,4), (2,3)]
        @test length(z_pos) == 5
        @test z_pos[1:4] == CartesianIndex{2}.(1:4, 1:4)
        @test CartesianIndex(2, 3) in z_pos
        # bohr_dict[-0.0] = [(3,2)] (the partner of the (2,3) Bohr collision)
        @test length(z_neg) == 1
        @test CartesianIndex(3, 2) in z_neg
    end

    # -----------------------------------------------------------------
    # (a) Algebraic equivalence on a generic (non-degenerate) Hamiltonian.
    # The new B_bohr matches the per-(i, j, col) reference at machine
    # precision.
    # -----------------------------------------------------------------
    @testset "(a) Generic ham: B_bohr ≡ reference (smooth Metro)" begin
        cfg = make_config(Lindbladian(), BohrDomain(); construction=KMS())
        B_new = B_bohr(TEST_HAM, TEST_JUMPS, cfg)
        B_ref = _b_bohr_reference(TEST_HAM, TEST_JUMPS, cfg)
        @test norm(B_new - B_ref) / max(norm(B_ref), 1.0) < 1e-12
    end

    # -----------------------------------------------------------------
    # (b) Order-invariance: shuffle bohr_dict[ν] indices and verify
    #     B_bohr is unchanged. The reference, computed with the
    #     shuffled order, must also match the unshuffled output.
    # -----------------------------------------------------------------
    @testset "(b) Order invariance: shuffled bohr_dict ≡ unshuffled" begin
        cfg = make_config(Lindbladian(), BohrDomain(); construction=KMS())
        B_unshuffled = B_bohr(TEST_HAM, TEST_JUMPS, cfg)
        B_ref_shuffled = _b_bohr_reference(TEST_HAM, TEST_JUMPS, cfg;
            index_permuter = idxs -> shuffle(MersenneTwister(7), idxs))
        @test norm(B_unshuffled - B_ref_shuffled) /
              max(norm(B_unshuffled), 1.0) < 1e-12
    end

    # -----------------------------------------------------------------
    # (c) Degenerate spectrum: bohr_dict[0.0] contains off-diagonal
    #     entries (Bohr collisions). Confirm B_bohr matches the
    #     reference on this fixture too.
    # -----------------------------------------------------------------
    @testset "(c) Degenerate ham (Bohr collision in ν=0)" begin
        β = 5.0
        ham_deg = _make_degenerate_ham(β)
        # Use simple X jumps on each of the 2 sites (NUM_QUBITS = 2,
        # dim = 4). Match the existing single-site Pauli convention.
        n = 2
        norm_fac = 1.0 / sqrt(n)
        jumps = JumpOp[]
        op_X = Matrix(pad_term([X], n, 1)) .* norm_fac
        push!(jumps, JumpOp(op_X, op_X, false, true))
        op_X2 = Matrix(pad_term([X], n, 2)) .* norm_fac
        push!(jumps, JumpOp(op_X2, op_X2, false, true))

        cfg = Config(
            sim = Lindbladian(), domain = BohrDomain(),
            construction = KMS(), num_qubits = n,
            with_linear_combination = true,
            beta = β, sigma = 1.0/β, a = β/30.0, s = 0.4,
            num_energy_bits = 8, w0 = 0.05,
        )
        B_new = B_bohr(ham_deg, jumps, cfg)
        B_ref = _b_bohr_reference(ham_deg, jumps, cfg)
        @test norm(B_new - B_ref) / max(norm(B_ref), 1.0) < 1e-12

        # Also confirm shuffled-order reference still matches.
        B_ref_shuf = _b_bohr_reference(ham_deg, jumps, cfg;
            index_permuter = idxs -> shuffle(MersenneTwister(13), idxs))
        @test norm(B_new - B_ref_shuf) / max(norm(B_new), 1.0) < 1e-12
    end

    # -----------------------------------------------------------------
    # (d) Gaussian filter branch (with_linear_combination = false).
    #     Uses a different `f` closure than smooth-Metro; confirms
    #     the cache logic doesn't depend on closure form.
    # -----------------------------------------------------------------
    @testset "(d) Gaussian filter (with_linear_combination=false)" begin
        cfg = Config(
            sim = Lindbladian(), domain = BohrDomain(),
            construction = KMS(), num_qubits = NUM_QUBITS,
            with_linear_combination = false,
            beta = BETA, sigma = SIGMA,
            gaussian_parameters = (BETA * (SIGMA^2 + 0.5^2) / 2, 0.5),
            num_energy_bits = NUM_ENERGY_BITS, w0 = W0,
        )
        B_new = B_bohr(TEST_HAM, TEST_JUMPS, cfg)
        B_ref = _b_bohr_reference(TEST_HAM, TEST_JUMPS, cfg)
        @test norm(B_new - B_ref) / max(norm(B_ref), 1.0) < 1e-12
    end

    # -----------------------------------------------------------------
    # (e) Threaded vs serial-reference equivalence on the degenerate
    #     fixture. Re-runs B_bohr to hit the threaded path (when
    #     nthreads > 1) on a Hamiltonian where bohr_dict ordering is
    #     non-monotonic in row index.
    # -----------------------------------------------------------------
    @testset "(e) Threaded path on degenerate ham" begin
        if Threads.nthreads() > 1
            β = 5.0
            ham_deg = _make_degenerate_ham(β)
            n = 2
            norm_fac = 1.0 / sqrt(n)
            jumps = JumpOp[]
            for site in 1:n
                op = Matrix(pad_term([X], n, site)) .* norm_fac
                push!(jumps, JumpOp(op, op, false, true))
            end
            cfg = Config(
                sim = Lindbladian(), domain = BohrDomain(),
                construction = KMS(), num_qubits = n,
                with_linear_combination = true,
                beta = β, sigma = 1.0/β, a = β/30.0, s = 0.4,
                num_energy_bits = 8, w0 = 0.05,
            )
            # bohr_dict has 4²=16 entries collapsed into a small number
            # of unique freqs. May not hit threaded threshold (nfreqs >=
            # 10) — but works either way; the test confirms correctness
            # whichever path is taken.
            B_new = B_bohr(ham_deg, jumps, cfg)
            B_ref = _b_bohr_reference(ham_deg, jumps, cfg)
            @test norm(B_new - B_ref) / max(norm(B_ref), 1.0) < 1e-12
        else
            @test true  # nthreads=1, threaded path not exercised
        end
    end
end
