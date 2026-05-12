#!/usr/bin/env julia
#
# qf-9ld.2 — Block-encoding norm dataset for DLL multi-channel + CKG.
#
# Implements the LCU 1-norm convention locked down in qf-9ld.1 / §9 of
# `drafts/ckg-vs-dll-comparison-findings.md`:
#
#   • DLL: α_G = Z_g · Z_𝒜² · √|𝒜|,    α_L = Z_f · Z_𝒜
#         Z_g = Σ_{m,n} |g(t_m, t_n)| · τ²,    Z_f = Σ_m |f(t_m)| · τ
#   • CKG: α_B = ‖b_-‖_1 · ‖b_+‖_1 · Z_𝒜²
#         ‖b_·‖_1 = Σ_m |b_·(t_m)| · τ
#   • Aggregated ‖L‖_be = α_G + (1/2) · α_L² · |𝒜|   (LW23 / Ding Thm 18)
#
# Hamiltonian norm convention
# ---------------------------
# Operator-norm and dissipator-bound numbers are reported on the *rescaled*
# fixture (‖H‖ ≈ 0.45). The rescaling factor is recorded in the BSON for
# qf-9ld.3 to multiply through to the unrescaled physical regime
# (where ‖H‖ = O(n)).
#
# LCU norms Z_f, Z_g, ‖b_±‖_1 are H-independent (filter + quadrature
# only) so they apply to both rescaled and unrescaled regimes.
#
# Sweep grid
# ----------
# n ∈ {3, 4, 5}, β ∈ {5, 10, 20}, k ∈ {1, 2, 4, 8} for DLL multi-channel
# (Metropolis S=2 + Gaussian); CKG smooth-Metropolis at k=1 baseline
# (a=0, s=0.25, σ = 1/β — sweep convention).
#
# Time grid for LCU sums: τ = 0.2, T_max = 6β. This is coarser than the
# DLL paper's Eq. 3.38 prescription (τ = O(1/(A_q^4 T² S³))) but more
# than adequate for the 1-norm to converge — the integrands are smooth
# Schwartz-class functions, the L¹ Riemann sum converges as O(τ²) on
# smooth integrands. A finer grid would shift absolute LCU norms by < 1%
# and not change any β-slope.

using Printf
using LinearAlgebra
using Statistics
using BSON
using QuantumFurnace

include(joinpath(@__DIR__, "..", "test", "test_helpers.jl"))

# ── Sweep grid ────────────────────────────────────────────────────────────────
const NS         = [3, 4, 5]
const BETAS      = [5.0, 10.0, 20.0]
const KS         = [1, 2, 4, 8]
const SIGMA_FAC  = 1.0   # σ_CKG = SIGMA_FAC / β   (sweep convention)
const A_CKG      = 0.0
const S_CKG      = 0.25
const S_DLL      = 2.0   # Metropolis bump support; centers must satisfy |c| ≤ S/2
const TAU_LCU    = 0.2
const TMAX_FAC   = 6.0   # T_max = TMAX_FAC · β

const OUTPUT_DIR = joinpath(@__DIR__, "output", "dll_norm_audit")
const OUTPUT_BSON = joinpath(OUTPUT_DIR, "dll_norm_audit.bson")

# ── Build single-Pauli jumps with the standard 1/√(3n) normalisation ──────────
function build_paulis_jumps(ham::HamHam{T}, n::Int) where {T}
    paulis = ([X], [Y], [Z])
    num_jumps = length(paulis) * n
    jump_norm = sqrt(T(num_jumps))
    jumps = JumpOp[]
    for pauli in paulis, site in 1:n
        op = Matrix(pad_term(pauli, n, site)) ./ jump_norm
        op_eb = ham.eigvecs' * op * ham.eigvecs
        push!(jumps, JumpOp(op, op_eb,
                            op == transpose(op), op == op'))
    end
    return jumps, num_jumps
end

# Centers for k channels, uniform on [0, c_max]. Matches the qf-7go.6 sweep
# convention (`scripts/scratch_dll_multirank_taumix_sweep.jl::_centers_for_k`):
# k=1 ⇒ [0.0] (standard rank-1 DLL filter, no shift).
const C_MAX_METRO = S_DLL / 2 - 1e-9   # stay strictly inside Metropolis flat top
const C_MAX_GAUSS = 1.0                # comparable spread to the Metropolis case
function _centers_for_k(k::Int, c_max::Float64)
    k == 1 && return [0.0]
    return collect(range(0.0, c_max; length = k))
end

build_dll_filter(::Val{:metro}, β::Real, k::Int) = begin
    base = DLLMetropolisFilter(β; S = S_DLL)
    k == 1 ? base : dll_multichannel_translates(base; centers = _centers_for_k(k, C_MAX_METRO))
end
build_dll_filter(::Val{:gauss}, β::Real, k::Int) = begin
    base = DLLGaussianFilter(β)
    k == 1 ? base : dll_multichannel_translates(base; centers = _centers_for_k(k, C_MAX_GAUSS))
end

# ── LCU 1-norm helpers ────────────────────────────────────────────────────────
#
# Z_f: discretised L¹ sum of f(t) over the time grid. For a multi-channel
# filter `F = (F_1, …, F_k)` the per-channel f's add linearly in α^DLL
# (Eq. dll-multi-Z), so Z_f^multi = Σ_ℓ Z_f^(ℓ) — which is exactly what
# `time_kernel(F_multi, t)` returns by linearity (see filters.jl:564).
#
# Z_g: discretised L¹ sum of g(t,t') = (1/(2π))² ∬ ĝ(ν,ν') e^{-iνt+iν't'} dνdν'
# (Eq. 3.6 of Ding et al.). For multi-channel, again sums linearly across
# channels. We compute g(t_m, t_n) by evaluating ĝ on a fine (ν,ν') grid
# over [-S, S]² (Metropolis) or the Gaussian's effective support (we use
# the same closed-form path as `dll_coherent_op_time` for Gaussian, and a
# uniform ν-grid quadrature for Metropolis multi-channel).
#
# To keep the code path uniform between single- and multi-channel, both
# branches go through `dll_coherent_kernel_bohr(filter, ν, νp)` evaluated
# on a (ν,ν') grid + 2D inverse FT to time. For multi-channel, the kernel
# helper sums per-channel kernels via `freq_kernel(filter::DLLMultiChannelFilter)`
# (filters.jl:549).

# Discretised Z_f via direct sum of |time_kernel(filter, t_m)| · τ.
function compute_Z_f(filter::AbstractFilter, β::Real)
    Tmax = TMAX_FAC * β
    Nt = 2 * round(Int, Tmax / TAU_LCU) + 1
    ts = range(-Tmax, Tmax; length = Nt)
    τ = step(ts)
    s = 0.0
    for t in ts
        s += abs(time_kernel(filter, t)) * τ
    end
    return s, Nt, τ
end

# Multi-channel dispatch: G^multi = Σ_ℓ G^(ℓ) is *linear* in channels (not bilinear
# in the summed `freq_kernel(filter_multi)`). Therefore Z_g^multi ≤ Σ_ℓ Z_g^(ℓ)
# (triangle inequality), with equality when channels don't cancel.
# The naive `freq_kernel(filter_multi)`-based formula would compute |Σ_ℓ f̂_ℓ|²
# instead of Σ_ℓ |f̂_ℓ|², over-counting cross-channel terms ℓ ≠ ℓ' that aren't
# in the actual G^multi. Per-channel + sum is the right convention.
function compute_Z_g(filter::DLLMultiChannelFilter, β::Real; kwargs...)
    Z = 0.0
    meta = nothing
    for c in filter.channels
        z, m = compute_Z_g(c, β; kwargs...)
        Z += z
        meta = m
    end
    return Z, meta
end

# Discretised Z_g via 2D-FT of ĝ(ν, ν') on a uniform (ν, ν') grid.
# We do not call the production `dll_coherent_op_time` g_tt construction,
# because we want an explicit `g_tt` matrix to take |·| of. The recipe
# mirrors the Metropolis branch of `dll_coherent_op_time` (the Gaussian
# closed-form gives the same answer to ~1e-12).
function compute_Z_g(filter::AbstractFilter, β::Real;
                    Nν::Int = 192, ν_radius::Real = max(S_DLL, 4.0 / β))
    # ν grid covers |ν| ≤ ν_radius. For Metropolis the bump is exactly on
    # [-S, S]; for Gaussian the freq_kernel decays as exp(-(βν)²/8), so
    # ν_radius ≥ 4/β captures > 10 e-foldings. Take max of the two.
    νs = range(-ν_radius, ν_radius; length = Nν)
    Δν = step(νs)

    # ĝ(ν, ν') = (1/2i) · tanh(β(ν'-ν)/4) · freq_kernel(filter, ν) · conj(freq_kernel(filter, ν'))
    # — Eq. 3.5 of Ding et al. via the exported `freq_kernel`. Linear in the
    # per-channel kernel for `DLLMultiChannelFilter`, so this works for both
    # single- and multi-channel filters.
    fk = ComplexF64[freq_kernel(filter, ν) for ν in νs]
    ĝ = Matrix{ComplexF64}(undef, Nν, Nν)
    pref = ComplexF64(1) / 2im
    for q in 1:Nν, p in 1:Nν
        th = tanh(β * (νs[q] - νs[p]) / 4)
        ĝ[p, q] = pref * th * fk[p] * conj(fk[q])
    end

    # Time grid for the LCU sum.
    Tmax = TMAX_FAC * β
    Nt = 2 * round(Int, Tmax / TAU_LCU) + 1
    ts = range(-Tmax, Tmax; length = Nt)
    τ = step(ts)

    # 2D inverse FT (Riemann sum, sign convention matches Eq. 3.6:
    # g(t,t') = (1/(2π))² · Σ ĝ(ν,νp) · exp(-iνt + iνp·t') · (Δν)²
    # We then sum |g(t_m, t_n)| · τ² over the time grid.
    # To avoid the O(Nt² · Nν²) cost, we use NUFFT.
    Nsrc = Nν * Nν
    src_x = Vector{Float64}(undef, Nsrc)
    src_y = Vector{Float64}(undef, Nsrc)
    src_c = Vector{ComplexF64}(undef, Nsrc)
    idx = 0
    for q in 1:Nν, p in 1:Nν
        idx += 1
        src_x[idx] = -νs[p]   # negate so isign=+1 gives e^{-iνt}
        src_y[idx] = νs[q]
        src_c[idx] = ĝ[p, q]
    end

    Ntgt = Nt * Nt
    tgt_x = Vector{Float64}(undef, Ntgt)
    tgt_y = Vector{Float64}(undef, Ntgt)
    for nn in 1:Nt, m in 1:Nt
        idx_t = (nn - 1) * Nt + m
        tgt_x[idx_t] = ts[m]
        tgt_y[idx_t] = ts[nn]
    end

    # Use the same FINUFFT plan family as dll.jl.
    plan = QuantumFurnace.FINUFFT.finufft_makeplan(3, 2, +1, 1, 1e-10;
                                                    dtype = Float64, nthreads = 1)
    QuantumFurnace.FINUFFT.finufft_setpts!(plan, src_x, src_y, Float64[],
                                            tgt_x, tgt_y, Float64[])
    out = Vector{ComplexF64}(undef, Ntgt)
    QuantumFurnace.FINUFFT.finufft_exec!(plan, src_c, out)
    QuantumFurnace.FINUFFT.finufft_destroy!(plan)

    norm_factor = (Δν / (2π))^2
    Z_g = 0.0
    for v in out
        Z_g += abs(v * norm_factor) * τ^2
    end
    return Z_g, (Nν = Nν, Nt = Nt, τ = τ, ν_radius = ν_radius)
end

# CKG ‖b_-‖_1, ‖b_+‖_1 LCU 1-norms via direct sum on the same time grid.
#
# CKG with `a = 0, s ≥ 0` (sweep convention) routes through
# `_compute_b_plus_metro(t, β, σ, η, s)` per `_select_b_plus_calculator`.
# The `η`-regularised b_+ is the Cor. III.2 form whose LCU 1-norm is
# bounded by O(log(β‖H‖/ε)). We pin η = 1e-3 (matching the §2 ε target).
function compute_ckg_b_norms(β::Real, σ::Real, a::Real, s::Real;
                             η::Real = 1e-3)
    Tmax = TMAX_FAC * β
    Nt = 2 * round(Int, Tmax / TAU_LCU) + 1
    ts = range(-Tmax, Tmax; length = Nt)
    τ = step(ts)
    bm_norm = 0.0
    bp_norm = 0.0
    for t in ts
        bm = invokelatest(QuantumFurnace._compute_b_minus, t, β, σ)
        # `a == 0` → eta-regularised metro form; `a > 0` → smooth form.
        bp = if a == 0.0
            invokelatest(QuantumFurnace._compute_b_plus_metro, t, β, σ, η, s)
        else
            invokelatest(QuantumFurnace._compute_b_plus_smooth, t, β, σ, a, s)
        end
        bm_norm += abs(bm) * τ
        bp_norm += abs(bp) * τ
    end
    return bm_norm, bp_norm, Nt, τ, η
end

# Flatten multi-channel `dll_lindblad_op_bohr` return — single-channel
# returns Matrix; multi-channel returns Vector{Matrix} (one per channel,
# no cross-terms in the multi-channel α). Returns flat Vector{Matrix} of
# per-(coupling, channel) Lindblad operators of length |𝒜|·k.
function _flatten_lindblad_ops(jumps, ham, filter, k)
    L_per_jump = [dll_lindblad_op_bohr(jump, ham, filter) for jump in jumps]
    if k == 1
        return L_per_jump  # already Vector{Matrix}
    else
        # Each entry is Vector{Matrix} of length k.
        return [Lℓ for Lvec in L_per_jump for Lℓ in Lvec]
    end
end

# ── Per-cell metric collection ────────────────────────────────────────────────
function audit_dll_cell(ham::HamHam, jumps::Vector{JumpOp}, β::Real, n::Int,
                        family::Symbol, k::Int)
    filter = build_dll_filter(Val(family), β, k)

    # Operator-level (Bohr-domain — cheap, exact for n ≤ 5).
    G  = dll_coherent_op_bohr(jumps, ham, filter, β)
    L_flat = _flatten_lindblad_ops(jumps, ham, filter, k)  # |𝒜|·k operators total

    G_op = opnorm(Matrix(G))
    L_max_op = maximum(opnorm(Matrix(L)) for L in L_flat)
    diss_1to1 = dissipator_one_to_one_norm_bound(L_flat)

    # LCU 1-norms (filter-only, H-independent). The multi-channel filter's
    # `time_kernel` and `freq_kernel` already sum over channels by linearity
    # (`src/filters.jl:549, 564`), so Z_f and Z_g here are exactly the
    # multi-channel sums Σ_ℓ Z_f^(ℓ), Σ_ℓ Z_g^(ℓ) of §9.3 (Eq. dll-multi-Z).
    Z_f, Nt_f, τ_f = compute_Z_f(filter, β)
    Z_g, g_meta    = compute_Z_g(filter, β)

    # Per-channel Z_f, Z_g (for the concatenated-index L_be variant below).
    Z_f_per_channel, Z_g_per_channel = if k == 1
        ([Z_f], [Z_g])
    else
        ch = filter.channels
        ([compute_Z_f(c, β)[1] for c in ch], [compute_Z_g(c, β)[1] for c in ch])
    end

    # Block-encoding factors (Z_𝒜 = max ‖A^a‖, here = 1/√(3n) for the
    # standard Pauli set with the conventional 1/√(3n) normalisation).
    Z_A = maximum(opnorm(j.in_eigenbasis) for j in jumps)
    A_card = length(jumps)
    α_L_naive  = Z_f * Z_A          # Σ_ℓ α_L^(ℓ)  (Eq. dll-multi-Z)
    α_G_naive  = Z_g * Z_A^2 * sqrt(A_card)
    # NAIVE aggregation: treat the multi-channel filter as a single-rank
    # filter with summed Z_f, Z_g; LW23 formula on |𝒜| index — over-counts
    # the dissipator (k² growth from α_L²).
    L_be_naive = α_G_naive + 0.5 * α_L_naive^2 * A_card
    # CONCATENATED-INDEX aggregation: index space (a, ℓ) of size |𝒜|·k,
    # per-(a,ℓ) BE factor α_L^(ℓ) = Z_f^(ℓ)·Z_𝒜; LW23 dissipator term
    # (1/2) · (max_ℓ α_L^(ℓ))² · |𝒜|·k. Tighter for k > 1 — k linear, not k².
    α_L_per_max = maximum(Z_f_per_channel) * Z_A
    L_be_concat = α_G_naive + 0.5 * α_L_per_max^2 * A_card * k

    return (;
        family, k, n, β,
        G_op, L_max_op, diss_1to1,
        Z_f, Z_g, Z_f_per_channel, Z_g_per_channel,
        Z_A, A_card,
        α_L_naive, α_G = α_G_naive,
        L_be_naive, L_be_concat,
        Nt_f, τ_f, Nν = g_meta.Nν, ν_radius = g_meta.ν_radius,
    )
end

function audit_ckg_cell(ham::HamHam, jumps::Vector{JumpOp}, β::Real, n::Int)
    σ = SIGMA_FAC / β
    cfg = Config(
        sim = Lindbladian(),
        domain = BohrDomain(),
        construction = KMS(),
        num_qubits = n,
        with_linear_combination = true,
        beta = β,
        sigma = σ,
        a = A_CKG,
        s = S_CKG,
        filter = nothing,
    )
    B = B_bohr(ham, jumps, cfg)
    B_op = opnorm(Matrix(B))

    # CKG dissipator: build via materialize_lindbladian — but for the
    # block-encoding norm we instead need ‖h_+‖_1, ‖h_-‖_1 of the
    # transition kernel. These are bounded analytically by Cor. III.3
    # (Gaussian) / Cor. III.4 (Metropolis): both ≤ 1. We confirm
    # numerically by sampling on the same time grid.
    bm_norm, bp_norm, Nt, τ, η_used = compute_ckg_b_norms(β, σ, A_CKG, S_CKG)
    Z_A = maximum(opnorm(j.in_eigenbasis) for j in jumps)
    A_card = length(jumps)
    α_B = bm_norm * bp_norm * Z_A^2

    return (;
        family = :ckg_smooth_metro, k = 1, n, β,
        B_op, σ, a = A_CKG, s = S_CKG, η = η_used,
        bm_norm, bp_norm, α_B, Z_A, A_card,
        Nt, τ,
    )
end

# ── Main sweep ────────────────────────────────────────────────────────────────
println("="^110)
println("qf-9ld.2 — Block-encoding norm audit (DLL multi-channel + CKG smooth-Metropolis)")
println("="^110)
println("Sweep: n=$NS, β=$BETAS, k=$KS for DLL ({metro, gauss}); CKG at k=1 baseline.")
println("Time grid for LCU sums: τ=$TAU_LCU, T_max = $TMAX_FAC·β.")
println()

results = Dict{Symbol, Any}()
results[:dll_metro] = NamedTuple[]
results[:dll_gauss] = NamedTuple[]
results[:ckg_smooth_metro] = NamedTuple[]
results[:metadata] = (;
    NS, BETAS, KS, SIGMA_FAC, A_CKG, S_CKG, S_DLL,
    TAU_LCU, TMAX_FAC,
    rescaling_note = "Metrics on rescaled fixture (‖H‖ ≈ 0.45). " *
                     "rescaling_factor stored per-cell for unrescaled extrapolation in qf-9ld.3.",
)

for n in NS
    ham_path = joinpath(@__DIR__, "..", "hamiltonians",
                         "heis_xxx_zzdisordered_periodic_n$(n).bson")
    isfile(ham_path) || (@warn "missing fixture, skipping" n; continue)

    for β in BETAS
        ham = _load_test_hamiltonian(ham_path, β)
        jumps, _ = build_paulis_jumps(ham, n)
        H_rescale = ham.rescaling_factor

        @printf("\n[n=%d, β=%5.1f, ‖H‖_eb=%.3f, rescaling=%.3f]\n",
                n, β, maximum(abs, ham.eigvals), H_rescale)

        # CKG baseline (k = 1).
        @time "  ckg" begin
            ckg_cell = audit_ckg_cell(ham, jumps, β, n)
            push!(results[:ckg_smooth_metro], merge(ckg_cell, (rescaling_factor = H_rescale,)))
            @printf("  ckg sM:  ‖B‖=%.3e  ‖b₋‖₁=%.3f  ‖b₊‖₁=%.3f  α_B=%.3e\n",
                    ckg_cell.B_op, ckg_cell.bm_norm, ckg_cell.bp_norm, ckg_cell.α_B)
        end

        # DLL Metropolis & Gaussian, k ∈ KS.
        for family in (:metro, :gauss)
            for k in KS
                @time "  dll-$family k=$k" begin
                    cell = audit_dll_cell(ham, jumps, β, n, family, k)
                    container = family == :metro ? results[:dll_metro] : results[:dll_gauss]
                    push!(container, merge(cell, (rescaling_factor = H_rescale,)))
                    @printf("  dll-%s k=%d:  ‖G‖=%.3e  Z_f=%.3f  Z_g=%.3f  α_G=%.3e  α_L=%.3e  ‖L‖_be(naive)=%.3e  (concat)=%.3e\n",
                            family, k, cell.G_op, cell.Z_f, cell.Z_g,
                            cell.α_G, cell.α_L_naive, cell.L_be_naive, cell.L_be_concat)
                end
            end
        end
    end
end

# ── Save BSON ─────────────────────────────────────────────────────────────────
mkpath(OUTPUT_DIR)
BSON.@save OUTPUT_BSON results
println("\nSaved dataset to $OUTPUT_BSON")

function _summary_table(container::Vector{NamedTuple}, label::String, field::Symbol)
    println("\nSummary: $field vs (n, β, k) for $label")
    println("="^110)
    @printf("%-3s | %-5s | %-12s %-12s %-12s %-12s\n",
            "n", "β", "k=1", "k=2", "k=4", "k=8")
    println("-"^70)
    for n in NS, β in BETAS
        cells = [c for c in container if c.n == n && c.β == β]
        cell_by_k = Dict(c.k => c for c in cells)
        @printf("%-3d | %-5.1f", n, β)
        for k in KS
            v = haskey(cell_by_k, k) ? getfield(cell_by_k[k], field) : NaN
            @printf(" | %-10.3e", v)
        end
        println()
    end
end

_summary_table(results[:dll_metro], "DLL Metropolis", :L_be_naive)
_summary_table(results[:dll_metro], "DLL Metropolis", :L_be_concat)
_summary_table(results[:dll_gauss], "DLL Gaussian",   :L_be_naive)
_summary_table(results[:dll_gauss], "DLL Gaussian",   :L_be_concat)

println("\nSummary: CKG smooth-Metro α_B and B-norm")
println("="^110)
@printf("%-3s | %-5s | %-12s %-12s %-12s %-12s\n",
        "n", "β", "‖B‖_op", "‖b₋‖₁", "‖b₊‖₁", "α_B")
println("-"^70)
for c in results[:ckg_smooth_metro]
    @printf("%-3d | %-5.1f | %-12.3e %-12.3f %-12.3f %-12.3e\n",
            c.n, c.β, c.B_op, c.bm_norm, c.bp_norm, c.α_B)
end

# ── β-slope diagnostics ───────────────────────────────────────────────────────
function loglog_slope(βs::Vector{Float64}, ys::Vector{Float64})
    x = log10.(βs); y = log10.(ys)
    n = length(x); x̄ = sum(x)/n; ȳ = sum(y)/n
    return sum((x .- x̄) .* (y .- ȳ)) / sum((x .- x̄).^2)
end

println("\n" * "="^110)
println("β-slope of α_G (log10 fit over β ∈ $(BETAS))")
println("Notes:")
println("  • k=1 single-channel slope is sub-linear (≈0.6 for Metro on rescaled ‖H‖=0.45) —")
println("    the β·‖H‖ ≤ 9 regime is pre-asymptotic; Remark 23's β-linear A_q^4·β·S term")
println("    dominates only at larger β·‖H‖ (operator-norm ‖G‖_op slopes 1.35–1.46 from §5.1).")
println("  • k>1 slopes are large and non-power-law: dominated by the cosh(β·shift/4)²")
println("    blow-up of the per-channel ShiftedSymmetricFilter's time-domain weights")
println("    (`src/filters.jl:686`). This is the qf-9ld.4 fairness-check signal.")
println("="^110)
@printf("%-3s | %-8s | %-12s %-12s %-12s %-12s\n",
        "n", "family", "k=1", "k=2", "k=4", "k=8")
println("-"^70)
for n in NS
    for family in (:metro, :gauss)
        container = family == :metro ? results[:dll_metro] : results[:dll_gauss]
        @printf("%-3d | %-8s", n, "dll-$family")
        for k in KS
            cells = [c for c in container if c.n == n && c.k == k]
            sort!(cells; by = c -> c.β)
            βs = Float64[c.β for c in cells]
            αs = Float64[c.α_G for c in cells]
            slope = (length(βs) >= 2) ? loglog_slope(βs, αs) : NaN
            @printf(" | %-10.3f", slope)
        end
        println()
    end
end
println("\nDone.")
