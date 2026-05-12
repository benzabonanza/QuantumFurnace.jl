"""
Tests for the grid-independent γ-norm replacement (qf-etx).

qf-etx.1: `pick_gamma_sup(config)` returns the closed-form continuum supremum
of γ (= 1.0 for every standard family). The supremum is verified two ways:
(1) γ evaluated at its closed-form maximiser equals 1.0;
(2) γ on a 2^16-point fine grid never exceeds 1.0.

Subsequent regression tests (register invariance, BohrDomain ↔ EnergyDomain
agreement, Krylov / simulator route invariance, GQSP α_be invariant) are
added as separate `@testset` blocks below as the qf-etx.{2..7} sub-issues
land.
"""

using Test
using QuantumFurnace
using LinearAlgebra

@testset "qf-etx.1: pick_gamma_sup closed-form is correct continuum sup" begin
    N_FINE = 2^16
    BETA = 10.0
    SIGMA = 0.1
    SIGMA_G = 0.5
    KMS_OMEGA_G = BETA * (SIGMA^2 + SIGMA_G^2) / 2  # = 1.3
    GNS_OMEGA_G = BETA * SIGMA_G^2 / 2              # = 1.25

    # Standard wide grid (halfwidth ±10/β = ±1.0). Peaks of γ for KMS lie at
    # ω = -βσ²/2 = -0.05 for kinky / a-reg / smooth, and at ω = -ω_γ = -1.3
    # for Gaussian; we use case-specific grids so the peak is well-resolved.
    function fine_grid_centered(center::Real; halfwidth::Real=1.0)
        # Shift the grid so the analytical peak lies exactly on a sample point.
        return collect(center .+ range(-halfwidth, halfwidth; length=N_FINE))
    end

    @testset "KMS Gaussian (γ ≤ 1)" begin
        cfg = Config(
            sim = Lindbladian(), domain = BohrDomain(), construction = KMS(),
            num_qubits = 3, with_linear_combination = false,
            beta = BETA, sigma = SIGMA,
            gaussian_parameters = (KMS_OMEGA_G, SIGMA_G),
            num_energy_bits = 8, w0 = 0.05,
        )
        γ = pick_transition(cfg)
        # Closed-form maximiser: ω = -ω_γ.
        @test isapprox(γ(-KMS_OMEGA_G), 1.0; atol=1e-15)
        ω = fine_grid_centered(-KMS_OMEGA_G; halfwidth=2.0)
        @test maximum(γ.(ω)) <= 1.0 + 1e-15
        @test pick_gamma_sup(cfg) == 1.0
    end

    @testset "KMS kinky Metropolis (s=0, a=0)" begin
        cfg = Config(
            sim = Lindbladian(), domain = BohrDomain(), construction = KMS(),
            num_qubits = 3, with_linear_combination = true,
            beta = BETA, sigma = SIGMA, a = 0.0, s = 0.0,
            num_energy_bits = 8, w0 = 0.05,
        )
        γ = pick_transition(cfg)
        # Closed-form maximiser: γ ≡ 1 on ω ≤ -βσ²/2.
        peak = -BETA * SIGMA^2 / 2
        @test isapprox(γ(peak), 1.0; atol=1e-15)
        @test isapprox(γ(peak - 1.0), 1.0; atol=1e-15)
        ω = fine_grid_centered(peak; halfwidth=2.0)
        @test maximum(γ.(ω)) <= 1.0 + 1e-15
        @test pick_gamma_sup(cfg) == 1.0
    end

    @testset "KMS smooth Metropolis (s=0.25, a=0) — locked thesis fixture" begin
        cfg = Config(
            sim = Lindbladian(), domain = BohrDomain(), construction = KMS(),
            num_qubits = 3, with_linear_combination = true,
            beta = BETA, sigma = SIGMA, a = 0.0, s = 0.25,
            num_energy_bits = 8, w0 = 0.05,
        )
        γ = pick_transition(cfg)
        # Smooth Metropolis sup is attained as ω → -∞; on any finite grid the
        # sample sup is strictly less than 1, but extending the grid pulls
        # it toward 1. This is exactly the discretisation residue the fix
        # eliminates from `gamma_norm_factor`. The closed-form factor `exp(4·sqrtA·sqrtB)`
        # in the smoothing term overflows for |ω| ≳ 70, so we evaluate at a
        # moderate ω = -10 where γ is already within 1e-4 of the analytical sup.
        ω_grid = collect(range(-10.0, 0.0; length=N_FINE))
        sup_grid = maximum(γ.(ω_grid))
        @test sup_grid <= 1.0 + 1e-12
        @test sup_grid >= 0.999
        @test pick_gamma_sup(cfg) == 1.0
    end

    @testset "GNS Gaussian / kinky / smooth all give pick_gamma_sup = 1.0" begin
        cases = [
            (false, nothing, nothing, (GNS_OMEGA_G, SIGMA_G), -GNS_OMEGA_G),  # Gaussian
            (true,  0.0,     0.0,     (nothing, nothing),   0.0),             # kinky
            (true,  0.0,     0.25,    (nothing, nothing),   nothing),         # smooth
        ]
        for (with_lc, a, s, gp, peak) in cases
            cfg = Config(
                sim = Lindbladian(), domain = BohrDomain(), construction = GNS(),
                num_qubits = 3, with_linear_combination = with_lc,
                beta = BETA, sigma = SIGMA,
                gaussian_parameters = gp,
                a = a, s = s,
                num_energy_bits = 8, w0 = 0.05,
            )
            @test pick_gamma_sup(cfg) == 1.0
            γ = pick_transition(cfg)
            if peak !== nothing
                @test isapprox(γ(peak), 1.0; atol=1e-15)
            else
                # smooth Metro GNS — sup approached as ω → -∞; evaluate at a
                # moderate ω where γ is within 1e-3 of the sup but the
                # `exp(4·sqrtA·sqrtB)` factor has not overflowed.
                @test γ(-10.0) >= 0.999
            end
            # Sample sup never exceeds 1 on a moderate grid.
            ω = collect(range(-10.0, 5.0; length=N_FINE))
            @test maximum(γ.(ω)) <= 1.0 + 1e-12
        end
    end
end

# ---------------------------------------------------------------------------
# qf-etx.2: `_precompute_data` populates `gamma_norm_factor = 1.0` for every
# CKG branch — i.e. construction is grid-independent post-fix.
# ---------------------------------------------------------------------------
@testset "qf-etx.2: _precompute_data has grid-independent gamma_norm_factor" begin
    n = 3
    src_root = dirname(@__DIR__)
    ham_path = joinpath(src_root, "hamiltonians", "heis_xxx_zzdisordered_periodic_n$(n).bson")
    ham = QuantumFurnace._load_hamiltonian_bson(ham_path, 10.0)

    register_pairs = [(8, 0.05), (10, 0.025)]
    cases = [
        (false, nothing, nothing, (1.3, 0.5), "Gaussian"),
        (true,  0.0,     0.0,     (nothing, nothing), "kinky Metro"),
        (true,  0.0,     0.25,    (nothing, nothing), "smooth Metro a=0 s=0.25"),
    ]

    for (r_D, w0_D) in register_pairs, (with_lc, a, s, gp, label) in cases
        cfg = Config(
            sim = Lindbladian(), domain = BohrDomain(), construction = KMS(),
            num_qubits = n, with_linear_combination = with_lc,
            beta = 10.0, sigma = 0.1,
            gaussian_parameters = gp,
            a = a, s = s,
            num_energy_bits = r_D, w0 = w0_D,
        )
        pd = QuantumFurnace._precompute_data(cfg, ham)
        @test pd.gamma_norm_factor ≈ 1.0 atol=1e-15
    end
end

# ---------------------------------------------------------------------------
# qf-etx.3: BohrDomain `construct_lindbladian` is byte-identical across
# `(r_D, w0_D)` register choices. Pre-fix, the two builds disagreed by the
# ratio of their respective `1.0 / maximum(transition.(...))` samples.
# ---------------------------------------------------------------------------
@testset "qf-etx.3: BohrDomain construct_lindbladian register invariance" begin
    n = 3
    src_root = dirname(@__DIR__)
    ham_path = joinpath(src_root, "hamiltonians", "heis_xxx_zzdisordered_periodic_n$(n).bson")
    ham = QuantumFurnace._load_hamiltonian_bson(ham_path, 10.0)

    jump_paulis = [[X], [Y], [Z]]
    jump_norm = sqrt(length(jump_paulis) * n)
    jumps = JumpOp[]
    for pauli in jump_paulis, site in 1:n
        op = Matrix(pad_term(pauli, n, site)) ./ jump_norm
        op_eb = ham.eigvecs' * op * ham.eigvecs
        push!(jumps, JumpOp(op, op_eb, op == transpose(op), op == op'))
    end

    cases = [
        (false, nothing, nothing, (1.3, 0.5), "Gaussian"),
        (true,  0.0,     0.0,     (nothing, nothing), "kinky Metro"),
        (true,  0.0,     0.25,    (nothing, nothing), "smooth Metro"),
    ]
    for (with_lc, a, s, gp, label) in cases
        cfg_a = Config(
            sim = Lindbladian(), domain = BohrDomain(), construction = KMS(),
            num_qubits = n, with_linear_combination = with_lc,
            beta = 10.0, sigma = 0.1,
            gaussian_parameters = gp, a = a, s = s,
            num_energy_bits = 8, w0 = 0.05,
        )
        cfg_b = Config(
            sim = Lindbladian(), domain = BohrDomain(), construction = KMS(),
            num_qubits = n, with_linear_combination = with_lc,
            beta = 10.0, sigma = 0.1,
            gaussian_parameters = gp, a = a, s = s,
            num_energy_bits = 10, w0 = 0.025,
        )
        L_a = construct_lindbladian(jumps, cfg_a, ham)
        L_b = construct_lindbladian(jumps, cfg_b, ham)
        @test isapprox(L_a, L_b; atol=1e-13, rtol=1e-13)
    end
end

# ---------------------------------------------------------------------------
# qf-etx.4: BohrDomain ↔ EnergyDomain agreement WITHOUT the script-side
# `L_test/gnf_test - L_ref/gnf_ref` workaround. Pre-fix raw `‖L_eng-L_bohr‖`
# was dominated by the gnf mismatch (~5% of ‖L‖); post-fix it is the actual
# quadrature error which is at machine precision for KMS smooth Metro at
# R_REF=10 (per qf-7xt: smooth Metro converges at any R_REF ≥ 8).
#
# Uses KMS construction (the thesis-canonical construction) including the
# coherent (Lamb-shift) term. BohrDomain and EnergyDomain coherent terms
# both call `B_bohr` (exact Bohr frequencies) — they should agree at
# machine precision regardless of register size.
# ---------------------------------------------------------------------------
@testset "qf-etx.4: BohrDomain ↔ EnergyDomain agreement (no /gnf workaround)" begin
    n = 3
    src_root = dirname(@__DIR__)
    ham_path = joinpath(src_root, "hamiltonians", "heis_xxx_zzdisordered_periodic_n$(n).bson")
    ham = QuantumFurnace._load_hamiltonian_bson(ham_path, 10.0)

    jump_paulis = [[X], [Y], [Z]]
    jump_norm = sqrt(length(jump_paulis) * n)
    jumps = JumpOp[]
    for pauli in jump_paulis, site in 1:n
        op = Matrix(pad_term(pauli, n, site)) ./ jump_norm
        op_eb = ham.eigvecs' * op * ham.eigvecs
        push!(jumps, JumpOp(op, op_eb, op == transpose(op), op == op'))
    end

    # ω-range from the unified principle (`scripts/scratch_dissipative_quadrature.jl`):
    # ω_max = ‖H‖ + 8σ; full range = 2·ω_max.
    # R_REF=12 per qf-7xt.3 (kinky converges with slope -2 in 1/N, so r=12
    # gives ~1e-7; smooth Metro and Gaussian are at machine precision).
    R_REF = 12
    omega_range = 2.0 * (opnorm(ham.data) + 8 * 0.1)
    w0_ref = omega_range / 2^R_REF

    cases = [
        (false, nothing, nothing, (1.3, 0.5), "KMS Gaussian", 1e-12),
        (true,  0.0,     0.0,     (nothing, nothing), "KMS kinky Metro", 1e-6),
        (true,  0.0,     0.25,    (nothing, nothing), "KMS smooth Metro", 1e-12),
    ]
    for (with_lc, a, s, gp, label, tol) in cases
        cfg_bohr = Config(
            sim = Lindbladian(), domain = BohrDomain(), construction = KMS(),
            num_qubits = n, with_linear_combination = with_lc,
            beta = 10.0, sigma = 0.1,
            gaussian_parameters = gp, a = a, s = s,
            num_energy_bits = R_REF, w0 = w0_ref,
        )
        cfg_eng = Config(
            sim = Lindbladian(), domain = EnergyDomain(), construction = KMS(),
            num_qubits = n, with_linear_combination = with_lc,
            beta = 10.0, sigma = 0.1,
            gaussian_parameters = gp, a = a, s = s,
            num_energy_bits = R_REF, w0 = w0_ref,
        )
        L_bohr = construct_lindbladian(jumps, cfg_bohr, ham)
        L_eng  = construct_lindbladian(jumps, cfg_eng,  ham)
        diff_op = opnorm(L_eng - L_bohr)
        @test diff_op <= tol
    end
end

# ---------------------------------------------------------------------------
# qf-etx.5: Krylov route — `apply_lindbladian!` matvec parity vs the dense
# `construct_lindbladian * vec(ρ)` (already a property of the codebase, but
# we re-verify it survives the gnf fix), and `krylov_spectral_gap`
# register-invariance for BohrDomain.
# ---------------------------------------------------------------------------
@testset "qf-etx.5: Krylov route invariance" begin
    n = 3
    src_root = dirname(@__DIR__)
    ham_path = joinpath(src_root, "hamiltonians", "heis_xxx_zzdisordered_periodic_n$(n).bson")
    ham = QuantumFurnace._load_hamiltonian_bson(ham_path, 10.0)

    jump_paulis = [[X], [Y], [Z]]
    jump_norm = sqrt(length(jump_paulis) * n)
    jumps = JumpOp[]
    for pauli in jump_paulis, site in 1:n
        op = Matrix(pad_term(pauli, n, site)) ./ jump_norm
        op_eb = ham.eigvecs' * op * ham.eigvecs
        push!(jumps, JumpOp(op, op_eb, op == transpose(op), op == op'))
    end

    @testset "apply_lindbladian! matvec parity vs construct_lindbladian (BohrDomain, KMS smooth Metro)" begin
        cfg = Config(
            sim = Lindbladian(), domain = BohrDomain(), construction = KMS(),
            num_qubits = n, with_linear_combination = true,
            beta = 10.0, sigma = 0.1, a = 0.0, s = 0.25,
            num_energy_bits = 10, w0 = 0.05,
        )
        L = construct_lindbladian(jumps, cfg, ham)
        dim = size(ham.data, 1)
        ws = Workspace(cfg, ham, jumps)
        # Random Hermitian, trace-normalised input ρ
        rng_seed = 4242
        ρ_init = (m = randn(QuantumFurnace.MersenneTwister(rng_seed), ComplexF64, dim, dim);
                   QuantumFurnace.hermitianize!(m); m ./= tr(m); m)
        out_dense = reshape(L * vec(ρ_init), dim, dim)
        apply_lindbladian!(ws, ρ_init, cfg, ham)
        out_krylov = copy(ws.scratch.rho_out)
        @test isapprox(out_krylov, out_dense; atol=1e-12, rtol=1e-12)
    end

    @testset "krylov_spectral_gap register invariance (BohrDomain, KMS smooth Metro)" begin
        cfg_a = Config(
            sim = Lindbladian(), domain = BohrDomain(), construction = KMS(),
            num_qubits = n, with_linear_combination = true,
            beta = 10.0, sigma = 0.1, a = 0.0, s = 0.25,
            num_energy_bits = 8, w0 = 0.05,
        )
        cfg_b = Config(
            sim = Lindbladian(), domain = BohrDomain(), construction = KMS(),
            num_qubits = n, with_linear_combination = true,
            beta = 10.0, sigma = 0.1, a = 0.0, s = 0.25,
            num_energy_bits = 10, w0 = 0.025,
        )
        res_a = krylov_spectral_gap(cfg_a, ham, jumps; krylovdim=20, tol=1e-12)
        res_b = krylov_spectral_gap(cfg_b, ham, jumps; krylovdim=20, tol=1e-12)
        # BohrDomain has no register-grid dependence post-fix; spectral gaps
        # must agree to machine precision.
        @test isapprox(res_a.spectral_gap, res_b.spectral_gap; atol=1e-10, rtol=1e-10)
    end
end

# ---------------------------------------------------------------------------
# qf-etx.6: Simulator routes (run_thermalize, predict_channel_trajectory)
# register invariance for BohrDomain. Pre-fix, two BohrDomain configs at
# different `(r_D, w0_D)` produced different `gamma_norm_factor` values,
# which fed into `jump_weight_scaling` in the trajectory simulator and
# polluted otherwise grid-independent results. Post-fix the two register
# choices must produce byte-identical final ρ trajectories.
# ---------------------------------------------------------------------------
@testset "qf-etx.6: Simulator routes register invariance" begin
    n = 3
    src_root = dirname(@__DIR__)
    ham_path = joinpath(src_root, "hamiltonians", "heis_xxx_zzdisordered_periodic_n$(n).bson")
    ham = QuantumFurnace._load_hamiltonian_bson(ham_path, 10.0)

    jump_paulis = [[X], [Y], [Z]]
    jump_norm = sqrt(length(jump_paulis) * n)
    jumps = JumpOp[]
    for pauli in jump_paulis, site in 1:n
        op = Matrix(pad_term(pauli, n, site)) ./ jump_norm
        op_eb = ham.eigvecs' * op * ham.eigvecs
        push!(jumps, JumpOp(op, op_eb, op == transpose(op), op == op'))
    end

    @testset "run_thermalize register invariance (BohrDomain, KMS smooth Metro)" begin
        delta = 0.01
        mixing_time = 0.05
        cfg_a = Config(
            sim = Thermalize(), domain = BohrDomain(), construction = KMS(),
            num_qubits = n, with_linear_combination = true,
            beta = 10.0, sigma = 0.1, a = 0.0, s = 0.25,
            num_energy_bits = 8, w0 = 0.05,
            mixing_time = mixing_time, delta = delta,
        )
        cfg_b = Config(
            sim = Thermalize(), domain = BohrDomain(), construction = KMS(),
            num_qubits = n, with_linear_combination = true,
            beta = 10.0, sigma = 0.1, a = 0.0, s = 0.25,
            num_energy_bits = 10, w0 = 0.025,
            mixing_time = mixing_time, delta = delta,
        )
        res_a = run_thermalize(jumps, cfg_a, ham)
        res_b = run_thermalize(jumps, cfg_b, ham)
        @test isapprox(res_a.final_dm, res_b.final_dm; atol=1e-12, rtol=1e-12)
        @test isapprox(res_a.trace_distances, res_b.trace_distances; atol=1e-10, rtol=1e-10)
    end

    @testset "predict_channel_trajectory register invariance (BohrDomain, KMS smooth Metro)" begin
        delta = 0.01
        cfg_a = Config(
            sim = Thermalize(), domain = BohrDomain(), construction = KMS(),
            num_qubits = n, with_linear_combination = true,
            beta = 10.0, sigma = 0.1, a = 0.0, s = 0.25,
            num_energy_bits = 8, w0 = 0.05,
            mixing_time = 0.5, delta = delta,
        )
        cfg_b = Config(
            sim = Thermalize(), domain = BohrDomain(), construction = KMS(),
            num_qubits = n, with_linear_combination = true,
            beta = 10.0, sigma = 0.1, a = 0.0, s = 0.25,
            num_energy_bits = 10, w0 = 0.025,
            mixing_time = 0.5, delta = delta,
        )
        rho0 = Matrix{ComplexF64}(I(2^n) / 2^n)
        k_grid = collect(0:5:50)
        res_a = predict_channel_trajectory(cfg_a, ham, jumps, rho0, k_grid; krylovdim=20, tol=1e-12)
        res_b = predict_channel_trajectory(cfg_b, ham, jumps, rho0, k_grid; krylovdim=20, tol=1e-12)
        # Distances along the trajectory must match (trace distance to ρ_inf).
        @test isapprox(res_a.distances, res_b.distances; atol=1e-10, rtol=1e-10)
        @test isapprox(res_a.spectral_gap, res_b.spectral_gap; atol=1e-10, rtol=1e-10)
    end
end

# ---------------------------------------------------------------------------
# qf-etx.7: GQSP block-encoding invariant `‖B_a / α_a‖_op ≤ 1` after the fix.
# The bound was originally proved assuming continuum γ ≤ 1; post-fix
# `gamma_norm_factor = 1.0 / pick_gamma_sup(config) = 1.0` realises that
# continuum bound exactly. `‖B_a‖_op` and `α_a` both inherit the same scalar
# γ_nf, so the ratio is invariant under the fix; the meaningful regression
# is that the invariant continues to hold across filter families and
# `gqsp_degree` choices.
# ---------------------------------------------------------------------------
@testset "qf-etx.7: GQSP α_be block-encoding invariant" begin
    n = 3
    src_root = dirname(@__DIR__)
    ham_path = joinpath(src_root, "hamiltonians", "heis_xxx_zzdisordered_periodic_n$(n).bson")
    ham = QuantumFurnace._load_hamiltonian_bson(ham_path, 10.0)

    jump_paulis = [[X], [Y], [Z]]
    jump_norm = sqrt(length(jump_paulis) * n)
    jumps = JumpOp[]
    for pauli in jump_paulis, site in 1:n
        op = Matrix(pad_term(pauli, n, site)) ./ jump_norm
        op_eb = ham.eigvecs' * op * ham.eigvecs
        push!(jumps, JumpOp(op, op_eb, op == transpose(op), op == op'))
    end

    R = 10
    w0 = 0.05
    t0 = 2π / (2^R * w0)

    # Three KMS filter families on TimeDomain (the GQSP path).
    cases = [
        (false, nothing, nothing, (1.3, 0.5), "KMS Gaussian"),
        (true,  0.0,     0.0,     (nothing, nothing), "KMS kinky Metro"),
        (true,  0.0,     0.25,    (nothing, nothing), "KMS smooth Metro"),
    ]
    for (with_lc, a, s, gp, label) in cases
        for d_gqsp in [1, 2, 4]
            cfg = Config(
                sim = Thermalize(), domain = TimeDomain(), construction = KMS(),
                num_qubits = n, with_linear_combination = with_lc,
                beta = 10.0, sigma = 0.1, a = a, s = s,
                gaussian_parameters = gp,
                num_energy_bits = R, w0 = w0, t0 = t0,
                eta = 1e-3,
                mixing_time = 0.01, delta = 0.01,
                with_gqsp = true, gqsp_degree = d_gqsp,
            )
            pd = QuantumFurnace._precompute_data(cfg, ham)
            t0_outer = register_t0_b_minus(cfg)
            t0_inner = register_t0_b_plus(cfg)
            for jump in jumps
                B = QuantumFurnace.B_time([jump], ham, pd.b_minus, pd.b_plus,
                    t0_outer, t0_inner, cfg.beta, cfg.sigma)
                rmul!(B, pd.gamma_norm_factor)
                α_be = QuantumFurnace._gqsp_block_encoding_alpha(jump,
                    pd.b_minus, pd.b_plus, t0_outer, t0_inner, pd.gamma_norm_factor)
                @test opnorm(B) / α_be <= 1.0 + 1e-10
            end
        end
    end
end

# ---------------------------------------------------------------------------
# qf-96o: default_smooth_s β-scaling
# ---------------------------------------------------------------------------
@testset "qf-96o: default_smooth_s preserves absolute smoothing width σ·√s = 0.05" begin
    @testset "calibration point (β=10, σ=1/β) returns s = 0.25" begin
        @test QuantumFurnace.default_smooth_s(10.0, 0.1) == 0.25
    end

    @testset "constant σ·√s along σ = 1/β across β-sweep" begin
        for β in (2.0, 5.0, 10.0, 20.0, 50.0)
            σ = 1.0 / β
            s = QuantumFurnace.default_smooth_s(β, σ)
            @test σ * sqrt(s) ≈ 0.05 atol = 1e-12
            @test s ≈ (β / 20.0)^2 atol = 1e-12     # equivalent form on σ = 1/β
        end
    end

    @testset "constant σ·√s off σ = 1/β (σ = c/β, c ∈ {0.25, 2, 3})" begin
        β = 10.0
        for c in (0.25, 0.5, 2.0, 3.0)
            σ = c / β
            s = QuantumFurnace.default_smooth_s(β, σ)
            @test σ * sqrt(s) ≈ 0.05 atol = 1e-12
        end
    end

    @testset "monotone β-scaling at fixed c = σβ" begin
        # s increases with β when σ = c/β (σ shrinks → need more s).
        @test QuantumFurnace.default_smooth_s(20.0, 0.05) >
              QuantumFurnace.default_smooth_s(10.0, 0.10) >
              QuantumFurnace.default_smooth_s( 5.0, 0.20)
    end
end

# ---------------------------------------------------------------------------
# qf-nq5: validate_config! enforces the (a, s) Metropolis taxonomy.
# Kinky Metropolis is exactly (s = 0, a = 0); smooth Metropolis is
# (s > 0, any a ≥ 0). The (s = 0, a > 0) combination is rejected.
# ---------------------------------------------------------------------------
@testset "qf-nq5: validate_config! (a, s) Metropolis taxonomy" begin
    function _taxonomy_cfg(; a, s, construction=KMS(), domain=BohrDomain())
        Config(
            sim = Lindbladian(), domain = domain, construction = construction,
            num_qubits = 3, with_linear_combination = true,
            beta = 10.0, sigma = 0.1, a = a, s = s,
            num_energy_bits = 8, w0 = 0.05,
        )
    end
    # Allowed: kinky (s = a = 0)
    validate_config!(_taxonomy_cfg(a = 0.0, s = 0.0))
    # Allowed: smooth with a = 0, s > 0 (thesis-default branch)
    validate_config!(_taxonomy_cfg(a = 0.0, s = 0.25))
    # Allowed: smooth with a > 0, s > 0
    validate_config!(_taxonomy_cfg(a = 0.333, s = 0.4))
    # Rejected: a-regularised but unsmoothed (s = 0, a > 0)
    @test_throws ArgumentError validate_config!(_taxonomy_cfg(a = 1.0, s = 0.0))
    # Same rejection on GNS construction
    @test_throws ArgumentError validate_config!(_taxonomy_cfg(a = 1.0, s = 0.0, construction = GNS()))
end
