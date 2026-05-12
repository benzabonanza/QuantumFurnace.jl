#!/usr/bin/env julia
# qf-bm1 / Phase 1+2 diagnostic — Non-Hermitian jumps end-to-end KMS-DB check.
#
# Strategy:
#   * Build a (σ⁺, σ⁻) jump pair on a single site of the n=3 disordered
#     Heisenberg fixture (paired to satisfy the Kossakowski α-skew-symmetry
#     requirement; cf. memory: feedback_non_hermitian_jumps.md).
#   * Construct the CKG-KMS Lindbladian with smooth-Metropolis filter in
#     each of {Bohr, Energy, Time, Trotter} domains.
#   * Run `verify_detailed_balance` and report `relative_norm`, the
#     anti-Hermitian discriminant magnitude.
#   * Cross-compare against a Hermitian baseline (single X jump on the
#     same site) to disambiguate quadrature / discretisation noise from a
#     genuine KMS-DB violation.
#
# Hypothesis under test:
#   * If `B_bohr` already contains the structure required by the
#     non-Hermitian extension (per-jump A†(ω₁) f(ω₁,ω₂) A(ω₂) summed over
#     each jump in the paired set), then BohrDomain should hit ≤ 1e-10.
#   * Energy/Time inherit Bohr correctness via discretisation; the gap is
#     pure quadrature error.
#   * Trotter inherits Time + Trotter-error.
#
# # PHYSICS CHECK: paired (σ⁺, σ⁻) is the minimal physical fixture.
# Single non-Hermitian jumps cannot satisfy KMS-DB because the
# α(ω₁,ω₂) = α(-ω₂,-ω₁) e^{-β(ω₁+ω₂)/2} skew-symmetry needs both
# (A, A†) Bohr decompositions in the jump set.

using QuantumFurnace
using LinearAlgebra
using Printf

# Resolve fixture path the same way `test_helpers.jl` does, so this
# script runs from any working directory.
const SOURCE_ROOT = dirname(@__DIR__)

#* PAULIs ------------------------------------------------------------------
# σ⁺ = (X + iY)/2 = [0 1; 0 0]; σ⁻ = (X - iY)/2 = [0 0; 1 0]; σ⁺' = σ⁻.
const SIGMA_PLUS  = ComplexF64[0 1; 0 0]
const SIGMA_MINUS = ComplexF64[0 0; 1 0]

#* Build n=3 system ---------------------------------------------------------
function _load_n3_ham(beta)
    ham_path = joinpath(SOURCE_ROOT, "hamiltonians", "heis_xxx_zzdisordered_periodic_n3.bson")
    return QuantumFurnace._load_hamiltonian_bson(ham_path, Float64(beta))
end

# `JumpOp[]` followed by `push!` would normally produce `Vector{JumpOp}`,
# but Julia's type-narrowing can sometimes concretise to
# `Vector{JumpOp{Matrix{ComplexF64}}}` (which is a strict subtype, not
# accepted by the `Vector{JumpOp}` method signature). We construct the
# vector with explicit `Vector{JumpOp}(undef, n)` then assign in-place,
# which guarantees the abstract eltype across all callers.
function _wrap_jumps(ops_with_data::Vector)
    jumps = Vector{JumpOp}(undef, length(ops_with_data))
    for (i, j) in pairs(ops_with_data)
        jumps[i] = j
    end
    return jumps
end

"""
    build_nh_pair_jumps(ham; site=1, normalize=true) -> Vector{JumpOp}

Single-site (σ⁺, σ⁻) jump pair on the requested site. With `normalize=true`
the pair is rescaled by 1/√2 so the L¹ jump-norm sum is comparable to a
single Hermitian Pauli. `hermitian = false` is set explicitly. `orthogonal`
is false here; the operators are not symmetric under transpose.
"""
function build_nh_pair_jumps(ham::HamHam; site::Int=1, normalize::Bool=true,
                              basis::AbstractMatrix = ham.eigvecs)
    n = Int(log2(size(ham.data, 1)))
    @assert 1 <= site <= n
    norm_fac = normalize ? 1.0 / sqrt(2) : 1.0
    ops = JumpOp[]
    for op2x2 in (SIGMA_PLUS, SIGMA_MINUS)
        op = Matrix(pad_term([op2x2], n, site)) .* norm_fac
        op_eb = basis' * op * basis
        push!(ops, JumpOp(op, op_eb, false, false))
    end
    return _wrap_jumps(ops)
end

"""
    build_hermitian_baseline(ham; site=1, normalize=true) -> Vector{JumpOp}

Hermitian X jump on the same site, with the same overall norm (1/√2 if
`normalize=true`). Used as the rel-norm baseline.
"""
function build_hermitian_baseline(ham::HamHam; site::Int=1, normalize::Bool=true,
                                   basis::AbstractMatrix = ham.eigvecs)
    n = Int(log2(size(ham.data, 1)))
    norm_fac = normalize ? 1.0 / sqrt(2) : 1.0
    op = Matrix(pad_term([X], n, site)) .* norm_fac
    op_eb = basis' * op * basis
    return _wrap_jumps([JumpOp(op, op_eb, false, true)])
end

"""
    build_unpaired_nh(ham; site=1, normalize=true) -> Vector{JumpOp}

Single non-Hermitian jump (σ⁺ alone) — Q1 violator. Should produce a
genuinely-non-zero KMS-DB residual that nothing in the production code
can fix.
"""
function build_unpaired_nh(ham::HamHam; site::Int=1, normalize::Bool=true,
                            basis::AbstractMatrix = ham.eigvecs)
    n = Int(log2(size(ham.data, 1)))
    norm_fac = normalize ? 1.0 / sqrt(2) : 1.0
    op = Matrix(pad_term([SIGMA_PLUS], n, site)) .* norm_fac
    op_eb = basis' * op * basis
    return _wrap_jumps([JumpOp(op, op_eb, false, false)])
end

#* Config factory -----------------------------------------------------------
# Wider grid than the 4-qubit test default — the smooth-Metropolis filter
# at β=1 is ~5× wider than at β=10, and a 12-bit Eb=12 register at w0=0.05
# (range [-1, +1]) cuts off the tail. Enlarging num_energy_bits=14 and
# w0=0.1 gives [-2.05, +2.05] which fully encloses the integrand at β=1.
const NUM_ENERGY_BITS = 14
const W0 = 0.1
const T0 = 2pi / (2^NUM_ENERGY_BITS * W0)
const NUM_TROTTER_STEPS_PER_T0 = 10

function nh_config(domain; beta::Real, num_qubits::Int=3)
    Config(;
        sim = Lindbladian(),
        domain = domain,
        construction = KMS(),
        num_qubits = num_qubits,
        with_linear_combination = true,
        beta = Float64(beta),
        sigma = 1.0 / Float64(beta),
        a = beta / 30.0,
        s = 0.4,
        num_energy_bits = NUM_ENERGY_BITS,
        w0 = W0,
        t0 = T0,
        num_trotter_steps_per_t0 = NUM_TROTTER_STEPS_PER_T0,
    )
end

#* Per-domain construction --------------------------------------------------
"""
    trotter_gibbs_target(ham, trotter) -> Hermitian{ComplexF64}

KMS-DB target for a TrotterDomain Lindbladian. Gibbs of the **original H**
(not H_T) expressed in the Trotter eigenbasis. The TrotterDomain
construction is designed to drive toward `gibbs(H)`; Trotter approximation
only affects the unitary evolution, not the fixed point. Empirical M-sweep
(scripts/scratch_trotter_M_sweep.jl) confirms `verify_detailed_balance(L, gibbs_H_in_T)
≈ 1e-10` independent of M.

`ham.gibbs` is diagonal in the H eigenbasis. We transform to computational
basis via `U_H * gibbs * U_H'`, then to Trotter eigenbasis via `U_T' · ... · U_T`.
"""
function trotter_gibbs_target(ham::HamHam, trotter::TrottTrott)
    gibbs_comp = ham.eigvecs * ham.gibbs * ham.eigvecs'
    return Hermitian(trotter.eigvecs' * gibbs_comp * trotter.eigvecs)
end

function build_lindbladian_for_domain(jumps_builder::Function, ham, beta, domain;
                                       allow_unpaired_nonhermitian::Bool=false)
    cfg = nh_config(domain; beta=beta, num_qubits=Int(log2(size(ham.data, 1))))
    if domain isa TrotterDomain
        # Use make_trotter_for_config to honour qf-d0w shared-δt₀ scheme
        # (per memory: qf_d0w_shared_delta_t0.md). The Trotter basis is
        # close to but not identical to the Hamiltonian eigenbasis, so
        # rebuild jumps in `trotter.eigvecs` AND the Gibbs target must be
        # the Trotter-Hamiltonian Gibbs (not ham.gibbs).
        trotter = make_trotter_for_config(ham, cfg)
        jumps = jumps_builder(ham; basis=trotter.eigvecs)
        L = Matrix{ComplexF64}(construct_lindbladian(jumps, cfg, ham; trotter=trotter,
            allow_unpaired_nonhermitian=allow_unpaired_nonhermitian))
        gibbs_T = trotter_gibbs_target(ham, trotter)
        return L, gibbs_T
    else
        jumps = jumps_builder(ham)
        L = Matrix{ComplexF64}(construct_lindbladian(jumps, cfg, ham;
            allow_unpaired_nonhermitian=allow_unpaired_nonhermitian))
        return L, Hermitian(ham.gibbs)
    end
end

#* Main diagnostic ----------------------------------------------------------
function run_diagnostic(beta::Real)
    println("=" ^ 78)
    @printf("β = %g — n=3 single-site jumps on site 1\n", beta)
    println("=" ^ 78)

    ham = _load_n3_ham(beta)
    gibbs = ham.gibbs

    # `(label, builder, allow_unpaired)` — the Q1 fixture uses an unpaired
    # non-Hermitian jump on purpose to exhibit the validation-bypass + the
    # KMS-DB violation it produces.
    fixtures = [
        ("Hermitian baseline (X)",       build_hermitian_baseline, false),
        ("Paired non-Hermitian (σ⁺,σ⁻)", build_nh_pair_jumps,     false),
        ("Unpaired σ⁺ (Q1 violator)",    build_unpaired_nh,        true),
    ]
    domains = [BohrDomain(), EnergyDomain(), TimeDomain(), TrotterDomain()]

    for (label, builder, allow) in fixtures
        println()
        println(label)
        for dom in domains
            try
                L, gibbs_target = build_lindbladian_for_domain(builder, ham, beta, dom;
                    allow_unpaired_nonhermitian=allow)
                res = verify_detailed_balance(L, gibbs_target; atol=1e-10)
                @printf("  %-12s  rel_norm = %.3e   fp_resid = %.3e   D_norm = %.3e\n",
                        string(typeof(dom).name.name),
                        res.relative_norm, res.fixed_point_residual, res.discriminant_norm)
            catch e
                @printf("  %-12s  ERROR: %s\n", string(typeof(dom).name.name),
                        first(sprint(showerror, e), 200))
            end
        end
    end
    println()
end

"""
    krylov_vs_dense_paired_nh(beta) -> nothing

Empirically verify the Krylov matvec route matches the dense
construct_lindbladian for paired (σ⁺, σ⁻) at the requested β. This pins
down whether downstream Krylov methods (krylov_spectral_gap,
lindblad_action_integrate, predict_lindbladian_trajectory) inherit
correctness for non-Hermitian jumps from the audit-clean dispatch in
`src/krylov_matvec.jl`.
"""
function krylov_vs_dense_paired_nh(beta::Real)
    println("=" ^ 78)
    @printf("Krylov ↔ dense matvec equivalence for paired NH at β = %g\n", beta)
    println("=" ^ 78)

    ham = _load_n3_ham(beta)
    dim = size(ham.data, 1)

    using_random_density_matrix = QuantumFurnace.random_density_matrix
    rho = Matrix(using_random_density_matrix(Int(log2(dim))))

    for dom in (EnergyDomain(), TimeDomain())
        cfg = nh_config(dom; beta=beta, num_qubits=Int(log2(dim)))
        jumps = build_nh_pair_jumps(ham)

        L_dense = Matrix{ComplexF64}(construct_lindbladian(jumps, cfg, ham))
        result_dense = reshape(L_dense * vec(rho), dim, dim)

        ws = QuantumFurnace.Workspace(cfg, ham, jumps)
        result_krylov = copy(apply_lindbladian!(ws, rho, cfg, ham))

        # Adjoint
        ws_adj = QuantumFurnace.Workspace(cfg, ham, jumps)
        result_dense_adj = reshape(L_dense' * vec(rho), dim, dim)
        result_krylov_adj = copy(apply_adjoint_lindbladian!(ws_adj, rho, cfg, ham))

        @printf("  %-12s  fwd ‖dense - krylov‖ = %.3e   adj ‖dense - krylov‖ = %.3e\n",
                string(typeof(dom).name.name),
                norm(result_dense - result_krylov),
                norm(result_dense_adj - result_krylov_adj))
    end
    println()
end

if abspath(PROGRAM_FILE) == @__FILE__
    for β in (1.0, 5.0, 10.0)
        run_diagnostic(β)
    end
    for β in (5.0, 10.0)
        krylov_vs_dense_paired_nh(β)
    end
end
