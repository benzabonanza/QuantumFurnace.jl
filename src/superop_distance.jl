# ---------------------------------------------------------------------------
# qf-72g: Matrix-free trace-norm distance between two superoperator propagators
# ---------------------------------------------------------------------------
#
# A general tool for measuring how close two (matrix-free) superoperator actions
# stay along a trajectory, in the operationally-meaningful TRACE norm. Either arm
# may be a Lindbladian action `L(ρ)` (propagated as `e^{tL}`) or a CPTP channel
# action `Φ_δ(ρ)` (propagated as `Φ_δ^k`). Any pairing works: L vs L (e.g.
# EnergyDomain quadrature-error checks at different r_D / r_bp), L vs Φ_δ (the
# qf-72g ideal-vs-implemented comparison), or Φ_δ vs Φ_δ.
#
# WHY TRACE NORM (qf-72g decision): the trace distance ½‖ρ−σ‖₁ = sup over POVMs of
# total-variation distance bounds the error in EVERY observable / mixing-time
# readout. It is matrix-free (one d×d SVD per snapshot), needs no adjoint map, and
# lower-bounds the induced 1→1 norm hence the diamond distance. Unlike the
# Hilbert–Schmidt norm it has no built-in √d = 2^{n/2} blind spot on the very axis
# (n) we make claims about.
#
# WHY EXACT PROPAGATION: the quantity of interest E = Φ_δ − e^{δL} is SMALL (it is
# an implementation error). Reconstructing each trajectory by a truncated Krylov
# spectral expansion (as `predict_*_trajectory` does) would inject a
# krylovdim-truncation error that is independent of — and at early times larger
# than — E, corrupting the measurement. So each arm is propagated as exactly as
# possible: the Lindbladian arm via `lindblad_action_integrate` (adaptive
# `KrylovKit.exponentiate`, tight tol) and the channel arm by exact iteration of
# its matvec.
#
# BASIS RECONCILIATION (correctness crux): the EnergyDomain Lindbladian matvec
# evolves ρ in the Hamiltonian eigenbasis (`ham.eigvecs`); the TrotterDomain
# channel matvec evolves ρ in the Trotter eigenbasis (`trotter.eigvecs`).
# Subtracting raw outputs is meaningless. Each `PropagatorArm` carries the unitary
# `V` (its working-basis vectors in the common/computational basis): ρ₀ is rotated
# IN via `V' ρ₀ V`, the arm evolves, and every snapshot is rotated OUT via
# `V ρ V'` to the common basis before the trace distance. Because the trace norm
# is unitarily invariant this is well defined; and because each arm is rotated by
# its OWN `V` (never the other's) the genuine O(δ²) fixed-point offset between the
# two maps is preserved, not gauged away.

"""
    PropagatorArm(apply!; kind, delta=NaN, basis=nothing, label="")

One arm of a [`propagator_trace_distance`](@ref) / [`propagator_fixed_point_distance`](@ref)
comparison: a matrix-free superoperator action plus the metadata needed to turn it
into a propagator and bring its output into a common basis.

# Arguments
- `apply!`: in-place matvec `apply!(out, x) -> out` acting on `d×d` matrices in the
  arm's WORKING basis. For a Lindbladian arm this computes `L(x)`; for a channel
  arm one full `Φ_δ(x)` step. (Reuse `apply_lindbladian!` / `apply_delta_channel!`
  via a closure, or supply [`lindbladian_arm`](@ref) / [`channel_arm`](@ref).)

# Keywords
- `kind::Symbol`: `:lindbladian` (propagate as `e^{tL}` via Krylov exponentiate)
  or `:channel` (propagate as `Φ_δ^k` by exact iteration).
- `delta::Real`: the channel step `δ` (required, `>0`, for `:channel`; ignored for
  `:lindbladian`). Used to map a physical time `t` to the integer step `k = t/δ`.
- `basis`: the unitary `V` whose columns are the working-basis vectors expressed in
  the common basis, so `ρ_common = V ρ_work V'`. `nothing` means the working basis
  IS the common basis (identity rotation) — appropriate e.g. when both arms are
  EnergyDomain Lindbladians compared in the shared `ham.eigvecs` basis.
- `label::AbstractString`: a human-readable tag carried into the result.
"""
struct PropagatorArm{F}
    apply!::F
    kind::Symbol
    delta::Float64
    basis::Union{Nothing, Matrix{ComplexF64}}
    label::String
end

function PropagatorArm(apply!::F; kind::Symbol, delta::Real = NaN,
                       basis::Union{Nothing, AbstractMatrix} = nothing,
                       label::AbstractString = "") where {F}
    kind in (:lindbladian, :channel) ||
        throw(ArgumentError("PropagatorArm kind must be :lindbladian or :channel (got :$kind)"))
    if kind === :channel
        (isfinite(delta) && delta > 0) ||
            throw(ArgumentError("a :channel PropagatorArm requires delta > 0 (got $delta)"))
    end
    b = basis === nothing ? nothing : Matrix{ComplexF64}(basis)
    return PropagatorArm{F}(apply!, kind, Float64(delta), b, String(label))
end

# Rotate a common-basis operator into the arm's working basis: V' X V.
_rotate_in(arm::PropagatorArm, X::AbstractMatrix) =
    arm.basis === nothing ? Matrix{ComplexF64}(X) : arm.basis' * X * arm.basis

# Rotate a working-basis operator back to the common basis: V X V'.
_rotate_out(arm::PropagatorArm, X::AbstractMatrix) =
    arm.basis === nothing ? X : arm.basis * X * arm.basis'

"""
    _iterate_channel_states(apply!, rho_0_work, k_grid; hermitize, renormalize)
        -> (states, ks_sorted, matvecs)

Exact iteration of a channel matvec: returns `Φ_δ^k(rho_0_work)` for every `k` in
`k_grid` (a `Dict`-indexed snapshot during a single 0→max(k_grid) sweep). No Krylov
truncation. `states[i]` corresponds to `ks_sorted[i]` (sorted unique `k`). Optional
defensive Hermitisation + trace renormalisation each step removes the channel's
≪1e-12 non-CPTP assembly drift (the genuine O(δ²) fixed-point offset is Hermitian
and trace-1, so it is NOT removed).
"""
function _iterate_channel_states(apply!::F, rho_0_work::Matrix{ComplexF64},
                                 k_grid::AbstractVector{<:Integer};
                                 hermitize::Bool = true,
                                 renormalize::Bool = true) where {F}
    d = size(rho_0_work, 1)
    ks = sort(unique(Int.(k_grid)))
    ks[1] >= 0 || throw(ArgumentError("k_grid must be non-negative (got min $(ks[1]))"))
    pos = Dict(k => i for (i, k) in enumerate(ks))
    states = Vector{Matrix{ComplexF64}}(undef, length(ks))
    rho = copy(rho_0_work)
    out = Matrix{ComplexF64}(undef, d, d)
    haskey(pos, 0) && (states[pos[0]] = copy(rho))
    matvecs = 0
    @inbounds for k in 1:ks[end]
        apply!(out, rho)
        matvecs += 1
        copyto!(rho, out)
        if hermitize
            for j in 1:d, i in 1:d
                rho[i, j] = (rho[i, j] + conj(rho[j, i])) / 2
            end
        end
        if renormalize
            tr_now = real(tr(rho))
            tr_now != 0 && (rho ./= tr_now)
        end
        haskey(pos, k) && (states[pos[k]] = copy(rho))
    end
    return states, ks, matvecs
end

# Propagate one arm and return its trajectory in the COMMON basis, aligned to the
# original `t_grid` order, plus (matvecs, converged).
#
# TIME CONVENTION (must match across arms): both arms interpret `t_grid` as
# ABSOLUTE physical time measured from `t = 0`, where the snapshot at `t = 0` is
# `rho_0` itself. The channel arm enforces this directly via `k = round(t/δ)`
# (so `t = 0 ⟺ k = 0 ⟺ rho_0`). The Lindbladian arm delegates to
# `lindblad_action_integrate`, which sets `states[1] = rho_0` and then integrates
# the RELATIVE steps `t_grid[i+1] - t_grid[i]`; that coincides with absolute time
# ONLY when `t_grid[1] == 0`. `propagator_trace_distance` therefore requires
# `t_grid[1] == 0` so the two arms are sampled at the same physical times.
function _propagate_arm(arm::PropagatorArm, rho_0::AbstractMatrix,
                        t_grid::AbstractVector{<:Real};
                        krylovdim::Integer, tol::Real,
                        hermitize::Bool, renormalize::Bool,
                        commensurate_atol::Real)
    rho_w = _rotate_in(arm, rho_0)
    if arm.kind === :lindbladian
        res = lindblad_action_integrate(
            arm.apply!, rho_w, zero(rho_w), t_grid;
            krylovdim = Int(krylovdim), tol = tol, save_states = true,
        )
        states_w = res.states
        return [_rotate_out(arm, s) for s in states_w], res.total_matvecs, res.all_converged
    else  # :channel
        k_of_t = Vector{Int}(undef, length(t_grid))
        @inbounds for (i, t) in enumerate(t_grid)
            kk = float(t) / arm.delta
            kr = round(Int, kk)
            abs(kk - kr) <= commensurate_atol * max(1.0, abs(kk)) || throw(ArgumentError(
                "channel arm '$(arm.label)': t=$t is not an integer multiple of δ=$(arm.delta) " *
                "(t/δ = $kk). Build t_grid as k_grid .* δ."))
            k_of_t[i] = kr
        end
        states_sorted, ks, mv = _iterate_channel_states(
            arm.apply!, rho_w, k_of_t; hermitize = hermitize, renormalize = renormalize)
        posmap = Dict(k => i for (i, k) in enumerate(ks))
        states_w = [states_sorted[posmap[k]] for k in k_of_t]
        return [_rotate_out(arm, s) for s in states_w], mv, true
    end
end

"""
    propagator_trace_distance(arm_A, arm_B, rho_0, t_grid; kwargs...) -> NamedTuple

Trace-distance between the trajectories of two superoperator propagators starting
from a common initial state `rho_0` (given in the COMMON / computational basis),
evaluated at the physical times in `t_grid`:

    T_k = ½ ‖ Φ_A(t_k)(rho_0) − Φ_B(t_k)(rho_0) ‖₁ ,

where each arm is propagated EXACTLY in its own working basis and rotated back to
the common basis before subtraction (see [`PropagatorArm`](@ref)). `t_grid` must be
sorted ascending and START at `t = 0` (`t_grid[1] == 0`): both arms read it as
absolute time from `t = 0` with the `t = 0` snapshot equal to `rho_0`, and the
Lindbladian arm's integrator measures time relative to `t_grid[1]`, so a nonzero
first node would silently desynchronise the two arms (the channel arm at `k = t/δ`
vs the Lindbladian arm at `t − t_grid[1]`). For a `:channel` arm every `t` in
`t_grid` must additionally be an integer multiple of that arm's `δ`.

This is the general engine behind the qf-72g ideal-Lindbladian-vs-implemented-
channel error-norm plot, and a reusable tool for quadrature-error checks (e.g.
two EnergyDomain Lindbladians at different `r_D`) at sizes where dense `d²×d²`
superoperators are infeasible.

# Keywords
- `krylovdim::Integer = 30`, `tol::Real = 1e-12`: forwarded to the Lindbladian
  arm's `lindblad_action_integrate` (Krylov `exponentiate`). The tight default tol
  keeps the ideal arm's reconstruction error well below the implementation error
  being measured.
- `hermitize::Bool = true`, `renormalize::Bool = true`: defensive cleanup of the
  CHANNEL arm's per-step output (removes ≪1e-12 non-CPTP drift; does not touch the
  genuine O(δ²) offset). No effect on the Lindbladian arm, which always
  re-Hermitises + renormalises.
- `save_states::Bool = false`: also return both common-basis trajectories.
- `commensurate_atol::Real = 1e-9`: tolerance for the channel `t = k·δ` check (also
  used for the `t_grid[1] == 0` check).

# Returns
NamedTuple with `t`, `trace_distances` (the `T_k`), `max_distance`, `argmax_t`,
`argmax_index`, `per_step_distance` (the `T_k` at the smallest `t_k > 0`, i.e. T₁
when `t_grid` is `k·δ`), `matvecs` (`(A, B)`), `converged` (`(A, B)`), `labels`
(`(A, B)`), and `states` (`(A, B)` common-basis trajectories) when `save_states`.
`trace_distances[1]` is `0` to machine precision (since `t_grid[1] == 0` ⇒ both
arms snapshot `rho_0`) — a built-in basis-correctness sanity check (a nonzero value
flags a rotation bug).
"""
function propagator_trace_distance(
    arm_A::PropagatorArm, arm_B::PropagatorArm,
    rho_0::AbstractMatrix, t_grid::AbstractVector{<:Real};
    krylovdim::Integer = 30, tol::Real = 1e-12,
    hermitize::Bool = true, renormalize::Bool = true,
    save_states::Bool = false, commensurate_atol::Real = 1e-9,
)
    d = size(rho_0, 1)
    size(rho_0, 2) == d || throw(ArgumentError("rho_0 must be square"))
    isempty(t_grid) && throw(ArgumentError("t_grid must be non-empty"))
    issorted(t_grid) || throw(ArgumentError("t_grid must be sorted ascending"))
    # Both arms read t_grid as absolute time from t=0; the Lindbladian arm's
    # integrator measures time relative to t_grid[1] (states[1]=rho_0), so a
    # nonzero first node desynchronises the arms. Require t_grid[1]==0 (the
    # documented T_0=0 invariant). See the TIME CONVENTION note on _propagate_arm.
    abs(float(t_grid[1])) <= commensurate_atol || throw(ArgumentError(
        "t_grid must start at 0 (got t_grid[1] = $(t_grid[1])); both arms measure " *
        "absolute time from t=0 with the t=0 snapshot equal to rho_0. Prepend 0.0."))

    statesA, mvA, cvA = _propagate_arm(arm_A, rho_0, t_grid;
        krylovdim = krylovdim, tol = tol, hermitize = hermitize,
        renormalize = renormalize, commensurate_atol = commensurate_atol)
    statesB, mvB, cvB = _propagate_arm(arm_B, rho_0, t_grid;
        krylovdim = krylovdim, tol = tol, hermitize = hermitize,
        renormalize = renormalize, commensurate_atol = commensurate_atol)

    T = [trace_distance_nh(statesA[i], statesB[i]) for i in eachindex(t_grid)]
    imax = argmax(T)

    # per-step T₁: trace distance at the smallest strictly-positive time.
    i_first_pos = findfirst(>(0), t_grid)
    per_step = i_first_pos === nothing ? T[1] : T[i_first_pos]

    return (
        t                 = collect(float.(t_grid)),
        trace_distances   = T,
        max_distance      = T[imax],
        argmax_t          = float(t_grid[imax]),
        argmax_index      = imax,
        per_step_distance = per_step,
        matvecs           = (mvA, mvB),
        converged         = (cvA, cvB),
        labels            = (arm_A.label, arm_B.label),
        states            = save_states ? (statesA, statesB) : nothing,
    )
end

# Generator matvec of an arm in its WORKING basis:
#   :lindbladian -> G = L  (apply! itself)
#   :channel     -> G = (Φ_δ − I)/δ  (finite-difference generator)
# The channel generator form shares Φ_δ's eigenVECTORS exactly (Φ_δ R = μ R ⇒
# G R = ((μ−1)/δ) R) but spreads the eigenvalues from the |μ|≈1 cluster to the
# O(1)-separated λ^eff = (μ−1)/δ, so the Re≈0 steady mode and the slow gap mode
# are well-conditioned at ANY δ (a direct |μ|≈1 read on Φ_δ degenerates as δ→0).
# Shared by `_arm_fixed_point` (steady = Re≈0 mode) and
# `slow_subspace_generator_distance` (the per-arm generator whose mismatch
# M = G_test − G_ref is projected onto the reference slow subspace).
function _arm_generator(arm::PropagatorArm)
    arm.kind === :channel || return arm.apply!
    return let f = arm.apply!, δ = arm.delta
        (out::AbstractMatrix, x::AbstractMatrix) -> (f(out, x); @. out = (out - x) / δ; out)
    end
end

# Leading eigenmode (steady state) of an arm, hermitised + trace-normalised, in
# the COMMON basis. Reuses `_krylov_spectral_decomposition`; the steady mode is
# the robust, seed-insensitive mode (unlike the slow gap mode).
#
# For a :channel arm we extract the fixed point from the GENERATOR form
# G = (Φ_δ − I)/δ (sorted :lindbladian) via `_arm_generator`, NOT from Φ_δ
# directly — see that helper for why this stays well-conditioned as δ → 0. This is
# essential for the qf-72g floor-controllability study, which sweeps δ down to ~1e-5.
function _arm_fixed_point(arm::PropagatorArm, rho_seed::AbstractMatrix;
                          krylovdim::Integer, tol::Real)
    seed_w = _rotate_in(arm, rho_seed)
    dim = size(seed_w, 1)
    decomp = _krylov_spectral_decomposition(
        _arm_generator(arm), Matrix{ComplexF64}(seed_w), dim;
        krylovdim = Int(krylovdim), tol = tol, sort_mode = :lindbladian)
    return _rotate_out(arm, decomp.rho_inf), decomp.matvec_count
end

"""
    propagator_fixed_point_distance(arm_A, arm_B; kwargs...) -> NamedTuple

Trace distance `½‖σ_A − σ_B‖₁` between the fixed points (steady states) of two
propagators — the `t → ∞` limit of [`propagator_trace_distance`](@ref). Each fixed
point is the leading (steady) eigenmode of the arm, brought to the common basis.
For a `:channel` arm the fixed point is extracted from the generator form
`(Φ_δ − I)/δ` rather than `Φ_δ`, so it stays well-conditioned as `δ → 0` (where all
`μ_i → 1` and a direct |μ|≈1 extraction degenerates) — essential for the qf-72g
floor-controllability study, which sweeps `δ`. For the qf-72g comparison this is the
channel's asymptotic Gibbs-state offset (the `floor_distance` the channel sweeps
report), computable at every `n` even when iterating to mixing is infeasible.

Note the generator form trades the |μ|≈1 degeneracy for an O(1/δ) amplification of
the matvec / threaded-reduction noise in the extracted `σ_A` (the conditioning of
the `(Φ_δ − I)/δ` Re≈0 eigenvector scales as `1/(δ·gap)`): empirically the
deviation from the dense channel fixed point is ~`7e-6` at `δ = 1e-3` and grows
~linearly to ~`1e-3` at `δ = 1e-4`. This is the inherent small-δ conditioning of
the quantity, not a method error (a basis/extraction bug would be O(1)); a δ-sweep
that goes well below `1e-4` should treat `σ_A` as accurate to ~`O(δ⁻¹·10⁻⁹)`.

# Keywords
- `rho_seed = nothing`: Arnoldi seed in the common basis (defaults to `I/d`, which
  is basis-invariant; the steady mode is captured regardless of the seed).
- `krylovdim::Integer = 40`, `tol::Real = 1e-10`: forwarded to the eigen-decomposition.

# Returns
NamedTuple with `distance`, `sigma_A`, `sigma_B` (common-basis fixed points),
`matvecs` (`(A, B)`), `labels` (`(A, B)`).
"""
function propagator_fixed_point_distance(
    arm_A::PropagatorArm, arm_B::PropagatorArm;
    rho_seed::Union{Nothing, AbstractMatrix} = nothing,
    d::Union{Nothing, Integer} = nothing,
    krylovdim::Integer = 40, tol::Real = 1e-10,
)
    dim = if rho_seed !== nothing
        size(rho_seed, 1)
    elseif arm_A.basis !== nothing
        size(arm_A.basis, 1)
    elseif d !== nothing
        Int(d)
    else
        throw(ArgumentError("provide rho_seed, or a basis on arm_A, or d to fix the dimension"))
    end
    seed = rho_seed === nothing ? Matrix{ComplexF64}(I(dim) / dim) : Matrix{ComplexF64}(rho_seed)

    σA, mvA = _arm_fixed_point(arm_A, seed; krylovdim = krylovdim, tol = tol)
    σB, mvB = _arm_fixed_point(arm_B, seed; krylovdim = krylovdim, tol = tol)
    return (
        distance = trace_distance_nh(σA, σB),
        sigma_A  = σA,
        sigma_B  = σB,
        matvecs  = (mvA, mvB),
        labels   = (arm_A.label, arm_B.label),
    )
end

# ---------------------------------------------------------------------------
# qf-e4z.48: ROBUST fixed-point extraction + fixed-point-vs-Gibbs distance
# ---------------------------------------------------------------------------
#
# The complement of `slow_subspace_generator_distance` (which measures the slow-RATE
# drift / gap distortion): here we measure how far the channel's STATIONARY STATE has
# drifted from the target Gibbs state, ½‖σ_Φ − ρ_β‖₁. The two are physically distinct
# axes — the gap distortion is blind to coherent (Lamb-shift) quadrature error that
# nonetheless drifts the fixed point (qf-e4z.48). Extracting the channel fixed point is
# made robust here by three lessons from qf-e4z.48:
#   • The STEADY mode (unlike the slow gap mode, qf-e4z.37) is the robust dominant mode,
#     but the generator-form Krylov `_arm_fixed_point` is SEED-sensitive: I/d is
#     near-orthogonal to the sharply-peaked steady eigenvector at large β_alg. Seeding
#     with ρ_β (the channel's fixed point ≈ ρ_β) fixes it (residual 7e-10 vs 2.6e-7).
#   • The EXACT dense superoperator is valid for the fixed point even though the
#     ℝ-linear hermitizing channel matvec injects d(d−1)/2 spurious μ=0 eigenvalues:
#     those live in the antisymmetric kernel; the physical μ≈1 fixed point is the
#     Hermitian eigenvector and is untouched. It is the ground truth at small d.
#   • A reference-free stationarity residual ½‖Φ(σ)−σ‖₁ certifies any extracted σ.

# Stationarity residual of a candidate fixed point `σ_common` (COMMON basis) under one
# application of the arm's action — a reference-free certificate (0 at the fixed point):
#   :channel     → ½‖Φ_δ(σ) − σ‖₁   (one channel step)
#   :lindbladian → ½‖𝓛(σ)‖₁          (generator action; trace-norm of 𝓛(σ))
function _arm_stationarity_residual(arm::PropagatorArm, sigma_common::AbstractMatrix)
    x = Matrix{ComplexF64}(_rotate_in(arm, sigma_common))
    out = Matrix{ComplexF64}(undef, size(x, 1), size(x, 2))
    arm.apply!(out, x)
    Y = Matrix{ComplexF64}(_rotate_out(arm, out))
    return arm.kind === :channel ?
        trace_distance_nh(Y, Matrix{ComplexF64}(sigma_common)) :
        trace_norm_nh(Y) / 2
end

# Exact dense fixed point of an arm via `build_dense_superoperator` (small `d` only):
# the eigenvector whose eigenvalue is closest to the arm's fixed-point target (μ=1 for a
# :channel, λ=0 for a :lindbladian generator), reshaped, Hermitised, trace-normalised,
# and rotated to the COMMON basis. For a :channel the ℝ-linear hermitizing matvec injects
# d(d−1)/2 spurious μ=0 eigenvalues in the antisymmetric kernel (qf-e4z.45), but the
# physical fixed point is the Hermitian μ≈1 eigenvector and is untouched — the returned
# anti-Hermitian fraction (~1e-13 when clean) is a built-in check. Returns
# `(σ_common, eigval, antiherm_frac, steady_gap)`, where `steady_gap` is the distance of
# the SECOND-closest eigenvalue to the fixed-point target — the spectral separation of the
# steady state from the rest (×δ for a channel). `steady_gap ≈ 0` flags a DEGENERATE steady
# space (symmetry block / decoherence-free subspace), where the single returned σ is an
# arbitrary member and the fixed point is not unique.
function _dense_arm_fixed_point(arm::PropagatorArm, d::Integer)
    S = build_dense_superoperator(arm.apply!, Int(d))
    F = eigen(S)
    target = arm.kind === :channel ? one(ComplexF64) : zero(ComplexF64)
    order = sortperm(abs.(F.values .- target))
    i1 = order[1]
    eval1 = F.values[i1]
    steady_gap = length(order) >= 2 ? abs(F.values[order[2]] - target) : Inf
    R = reshape(Vector{ComplexF64}(F.vectors[:, i1]), Int(d), Int(d))
    nR = norm(R)
    antiherm = nR == 0 ? 0.0 : norm((R .- R') ./ 2) / nR
    R = (R .+ R') ./ 2
    trR = real(tr(R)); trR != 0 && (R ./= trR)
    return Matrix{ComplexF64}(_rotate_out(arm, R)), eval1, antiherm, steady_gap
end

# Physical-validity certificate on an extracted fixed point (qf-dee #7): the steady
# state of a CPTP channel — and any Gibbs fixed point — is a density matrix, hence
# PSD. A minimum eigenvalue below `-psd_tol` is NOT float noise: it flags a
# non-physical fixed point (a parameter choice where the assembled map is not
# completely positive, or a corrupted extraction). Returns the min eigenvalue of the
# Hermitian part and warns when it is too negative.
function _fixed_point_min_eigval(sigma::AbstractMatrix, label::AbstractString;
                                 psd_tol::Real = 1e-10)
    λmin = minimum(real, eigvals(Hermitian(Matrix{ComplexF64}((sigma .+ sigma') ./ 2))))
    λmin < -psd_tol && @warn "arm_fixed_point: extracted fixed point ($label) is not " *
        "PSD — min eigenvalue $(round(λmin, sigdigits = 3)) < -psd_tol = " *
        "$(-float(psd_tol)). A steady state must be a valid density matrix; suspect a " *
        "non-CP parameter choice or a corrupted extraction, not float noise." maxlog = 1
    return λmin
end

"""
    arm_fixed_point(arm; seed=nothing, method=:auto, krylovdim=120, tol=1e-10,
                    residual_gate=1e-7, dense_max_dim=4096, psd_tol=1e-10) -> NamedTuple

Robustly extract the fixed point (steady state) of a [`PropagatorArm`](@ref), returned in
the COMMON basis, with a reference-free stationarity certificate. This is the robust
companion to [`fixed_point_gibbs_distance`](@ref) and a hardened alternative to the bare
generator-form extraction in [`propagator_fixed_point_distance`](@ref) (qf-e4z.48).

# Methods (`method`)
- `:krylov` — generator-form Arnoldi (`_arm_fixed_point` on `𝓛` / `(Φ_δ−I)/δ`). Fast and
  scalable, but the steady-mode overlap depends on `seed`: pass `ρ_β` for a `:channel`
  arm at large `β` (I/d is a poor seed for the sharply-peaked steady state). For a
  `:lindbladian` arm do NOT seed with its exact fixed point `ρ_β` — `𝓛(ρ_β)≈0` makes
  Arnoldi break down; use `I/d` (the default when `seed=nothing`).
- `:dense` — exact `build_dense_superoperator` eigendecomposition (ground truth; valid for
  the fixed point despite the channel's spurious kernel eigenvalues). Cost `O(d²)` matvecs
  + an `O(d⁶)` eigensolve, so only for small `d` (`d ≤ √dense_max_dim`).
- `:auto` (default) — try `:krylov`; accept it iff the stationarity residual `≤
  residual_gate`; otherwise escalate to `:dense` when `d ≤ dense_max_dim` (else return the
  ungated Krylov result with a warning). A Krylov breakdown also triggers the dense fallback.

# Uniqueness assumption
The extractors return the SINGLE mode closest to the fixed-point target. This is the unique
steady state for a primitive (ergodic) channel / Lindbladian — true for the disordered-
Heisenberg CKG sampler (full-rank gapped Gibbs fixed point). If the generator has a
DEGENERATE steady space (an exact symmetry block / decoherence-free subspace — e.g. an
ordered-phase doublet), the returned `sigma` is one arbitrary member; the `:dense` path
reports `steady_gap` (the second-closest eigenvalue's distance to the target) and warns when
it is below `degeneracy_tol`, so the caller can detect this.

# Returns
NamedTuple `(sigma, residual, eigval, antiherm_frac, steady_gap, min_eigval, method, matvecs,
converged)`: `sigma` the COMMON-basis fixed point; `residual` the [`_arm_stationarity_residual`](@ref)
certificate (`≈3e-9` for a clean dense channel fixed point); `eigval` the dense eigenvalue
(`μ≈1` channel / `λ≈0` Lindbladian; `NaN` for Krylov); `antiherm_frac` the dense
eigenvector's anti-Hermitian fraction (`NaN` for Krylov); `steady_gap` the steady-state
spectral separation (`NaN` for Krylov; small ⇒ degenerate, see above); `min_eigval` the
minimum eigenvalue of the Hermitised `sigma` (a PSD certificate — `≥ -psd_tol` for a
valid density-matrix fixed point; a `@warn` fires below `-psd_tol`, qf-dee #7); `method`
the path actually taken (`:krylov`, `:dense`, or `:krylov_ungated`).
"""
function arm_fixed_point(arm::PropagatorArm;
        seed::Union{Nothing, AbstractMatrix} = nothing,
        method::Symbol = :auto, krylovdim::Integer = 120, tol::Real = 1e-10,
        residual_gate::Real = 1e-7, dense_max_dim::Integer = 4096,
        degeneracy_tol::Real = 1e-9, psd_tol::Real = 1e-10)
    method in (:auto, :krylov, :dense) ||
        throw(ArgumentError("method must be :auto, :krylov or :dense (got :$method)"))
    d = arm.basis !== nothing ? size(arm.basis, 1) :
        (seed !== nothing ? size(seed, 1) :
         throw(ArgumentError("provide a seed, or an arm with a basis, to fix the dimension")))
    seed_m = seed === nothing ? Matrix{ComplexF64}(I(d) / d) : Matrix{ComplexF64}(seed)

    dense_result() = let (σ, μ, ah, sgap) = _dense_arm_fixed_point(arm, d)
        sgap < degeneracy_tol && @warn "arm_fixed_point: near-degenerate steady space " *
            "($(arm.label)): 2nd-closest eigenvalue is $(round(sgap, sigdigits = 3)) from the " *
            "fixed-point target (< degeneracy_tol = $degeneracy_tol). The returned σ is one " *
            "arbitrary member of a degenerate steady space — the fixed point is not unique." maxlog = 1
        (sigma = σ, residual = _arm_stationarity_residual(arm, σ), eigval = μ,
         antiherm_frac = ah, steady_gap = sgap,
         min_eigval = _fixed_point_min_eigval(σ, arm.label; psd_tol = psd_tol),
         method = :dense, matvecs = d * d, converged = true)
    end
    method === :dense && return dense_result()

    σk = nothing; mvk = 0; resk = NaN; broke = false
    try
        σk, mvk = _arm_fixed_point(arm, seed_m; krylovdim = Int(krylovdim), tol = tol)
        resk = _arm_stationarity_residual(arm, σk)
    catch err
        method === :krylov && rethrow(err)
        broke = true
        @warn "arm_fixed_point: Krylov broke down ($(arm.label)); escalating to dense" err maxlog = 1
    end

    # min_eigval of the (Krylov) σ — NaN when Krylov broke down (σk undefined).
    mineig_k() = broke ? NaN : _fixed_point_min_eigval(σk, arm.label; psd_tol = psd_tol)

    method === :krylov && return (sigma = σk, residual = resk, eigval = NaN,
        antiherm_frac = NaN, steady_gap = NaN, min_eigval = mineig_k(),
        method = :krylov, matvecs = mvk, converged = !broke)

    # :auto — accept a certified Krylov result, else escalate to dense when feasible.
    (!broke && isfinite(resk) && resk <= residual_gate) && return (sigma = σk,
        residual = resk, eigval = NaN, antiherm_frac = NaN, steady_gap = NaN,
        min_eigval = mineig_k(), method = :krylov, matvecs = mvk, converged = true)
    d <= dense_max_dim && return dense_result()
    @warn "arm_fixed_point: Krylov residual $(resk) > gate $(residual_gate) but d=$d > " *
          "dense_max_dim=$(dense_max_dim) — returning UNGATED Krylov fixed point (approximate)." maxlog = 1
    return (sigma = σk, residual = resk, eigval = NaN, antiherm_frac = NaN,
            steady_gap = NaN, min_eigval = mineig_k(), method = :krylov_ungated,
            matvecs = mvk, converged = !broke)
end

"""
    fixed_point_gibbs_distance(arm, gibbs; seed=:auto, kwargs...) -> NamedTuple

Robust trace distance `½‖σ_arm − ρ_β‖₁` between an arm's fixed point and an EXPLICIT
Gibbs state `gibbs` (`ρ_β`), both in the COMMON / computational basis — the qf-e4z.48
fixed-point-quality measure, designed to be read off alongside
[`slow_subspace_generator_distance`](@ref)'s gap distortion in the same sweep (so both
the dynamical and the thermodynamic error axes plot on one figure).

Build `ρ_β` in the computational basis with `gibbs_state(ham, β_alg)` (NOT `ham.gibbs`,
which is stored in the eigenbasis). The fixed point is extracted by [`arm_fixed_point`](@ref)
(robust, certified). `seed=:auto` picks the qf-e4z.48-validated seed by arm kind: `ρ_β`
for a `:channel` arm (its fixed point ≈ `ρ_β`, so overlap is high), `I/d` for a
`:lindbladian` arm (seeding with its exact fixed point `ρ_β` would break Arnoldi down).
Pass an explicit `seed` matrix or `nothing` (= `I/d`) to override; remaining `kwargs`
(`method`, `krylovdim`, `tol`, `residual_gate`, `dense_max_dim`) forward to
`arm_fixed_point`.

# Returns
NamedTuple `(distance, sigma, residual, eigval, antiherm_frac, steady_gap, min_eigval,
method, matvecs, converged, label)` — `distance` is the headline `½‖σ_arm − ρ_β‖₁`; the
rest is the `arm_fixed_point` certificate (check `residual` ≲ `residual_gate` to trust
`distance`, and `min_eigval ≥ -psd_tol` that the fixed point is a valid density matrix).
"""
function fixed_point_gibbs_distance(arm::PropagatorArm, gibbs::AbstractMatrix;
        seed::Union{Symbol, Nothing, AbstractMatrix} = :auto, kwargs...)
    size(gibbs, 1) == size(gibbs, 2) || throw(ArgumentError("gibbs must be square"))
    g = Matrix{ComplexF64}(gibbs)
    chosen = seed === :auto ? (arm.kind === :channel ? g : nothing) :
             (seed === nothing ? nothing : Matrix{ComplexF64}(seed))
    fp = arm_fixed_point(arm; seed = chosen, kwargs...)
    return (distance = trace_distance_nh(fp.sigma, g), sigma = fp.sigma,
            residual = fp.residual, eigval = fp.eigval, antiherm_frac = fp.antiherm_frac,
            steady_gap = fp.steady_gap, min_eigval = fp.min_eigval, method = fp.method,
            matvecs = fp.matvecs, converged = fp.converged, label = arm.label)
end

# Phase α (|α| = 1) that makes αR maximally Hermitian. For a non-degenerate
# real-eigenvalue mode R = e^{iθ}H (H Hermitian) of a KMS-normal generator, R† =
# e^{-2iθ}R, so α² = ⟨R,R†⟩_HS / ‖R‖² = e^{-2iθ} and α = e^{-iθ} gives αR = H. The
# four-way branch resolves the √ sign and rejects the anti-Hermitian (i·α) branch.
# A degenerate-cluster mode is a complex eigen-combination that NO single phase
# Hermitizes — the caller detects this via the residual ‖aherm(αR)‖/‖R‖.
function _hermitizing_phase(R::AbstractMatrix)
    z = dot(vec(R), vec(R'))                 # ⟨R, R†⟩_HS
    n2 = real(dot(vec(R), vec(R)))           # ‖R‖²_HS
    (abs(z) < 1e-300 || n2 == 0) && return one(ComplexF64)
    α = sqrt(z / n2)
    α /= abs(α)                              # enforce |α| = 1
    best = α
    bestval = -1.0
    for c in (α, -α, im * α, -im * α)
        v = norm((c .* R .+ (c .* R)') ./ 2)
        if v > bestval
            bestval = v
            best = c
        end
    end
    return best
end

"""
    slow_subspace_generator_distance(ref_arm, test_arm, rho_seed; kwargs...) -> NamedTuple

One scalar error per `(n, β)` measuring how far a TEST propagator's generator is
from a REFERENCE propagator's, restricted to the reference's slow (long-time)
subspace — the qf-e4z.45 channel-vs-Lindbladian operator-level distance. With
`(ref = lindbladian_arm, test = channel_arm)` this is "how far is the deployed CKG
channel `Φ_δ` from the ideal KMS-DB Lindbladian `𝓛`, on the modes that survive to
`τ_mix`", computed WITHOUT ever eigensolving the channel.

# The measure
Let `M = G_test − G_ref` be the generator mismatch (units of rate), where the
per-arm generator `G` is `𝓛` itself for a `:lindbladian` arm and the
finite-difference `(Φ_δ − I)/δ` for a `:channel` arm (`_arm_generator`). Let
`{R_k, L_k}` be the reference generator's `K` slowest biorthonormal right/left
eigenpairs (`⟨L_j, R_k⟩_HS = δ_jk`), ordered by `|Re λ|` ascending. Project `M`
onto that subspace and take the spectral norm of the small `K×K` matrix:

    M_jk = ⟨L_j | M | R_k⟩_HS = tr(L_j' (G_test − G_ref) R_k),
    ε_slow = ‖M_jk‖₂  (`eps_slow`, units of rate),  ε_slow / |λ₂| (`eps_slow_rel_gap`).

`K = 1` (`num_slow_modes = 1`, default) reduces to the relative gap shift
`|⟨L₂|M|R₂⟩| / |λ₂|` — first-order biorthogonal perturbation theory for the gap
eigenvalue. **`K = 1` is the robust, validated deliverable.** With
`include_stationary = true` a `k = 0` stationary rung (the steady mode, `λ ≈ 0`) is
prepended, subsuming the qf-72g fixed-point bias.

# Hermiticity requirement — why `K > 1` is only conditionally reliable (qf-e4z.45)
A `:channel` test arm's matvec (`apply_delta_channel!`) is Hermiticity-PRESERVING
and defensively hermitizes its output, hence is **ℝ-linear**: it acts faithfully
only on **Hermitian** inputs (on a non-Hermitian `X` it returns `Φ_δ(½(X+X†))`).
The reference's slow modes are Hermitian *up to a phase* (a KMS-normal `𝓛` has real
slow spectrum, so each non-degenerate eigenmode is `e^{iθ}H`), so this routine
**phase-fixes** every retained `R_k → αR_k` (Hermitian) and `L_k → αL_k` (preserving
`⟨L_j,R_k⟩ = δ_jk`) before applying the generators. The non-degenerate gap mode is
thereby made exactly Hermitian, so `K = 1` is clean. A *degenerate cluster*, however,
returns arbitrary complex eigen-combinations no single phase can Hermitize; if
`num_slow_modes` retains a PARTIAL cluster the modes stay non-Hermitian, the channel
hermitization corrupts the mismatch (an `O(1/δ)` artifact), and the projection onto a
partial degenerate eigenspace is itself ill-defined. The returned `max_antiherm_frac`
(the worst post-phase-fix anti-Hermitian fraction) flags this: it is `~1e-13` when the
retained modes are Hermitian (trust the result) and `O(1)` when a cluster is cut (a
`@warn` fires; do not trust `eps_slow`). Use `K > 1` only when `max_antiherm_frac`
stays small; otherwise the gap shift (`K = 1`) is the reliable scalar.

NOTE: build a dense `:channel` superoperator via `build_dense_superoperator(arm.apply!,
d)` ONLY for cross-checks on Hermitian inputs — feeding the non-Hermitian `E_ij`
basis through the hermitizing matvec makes `E_ij`/`E_ji` map to identical columns,
injecting `≈ d(d−1)/2` spurious zero eigenvalues (an artifact, not channel physics;
the true `Φ_δ` has all `|μ| ≈ 1`). Use the Hermitian operator basis to materialize it.

# Why PT on the reference modes, not an independent channel gap (qf-e4z.37 / qf-9lp)
The true channel-vs-`𝓛` gap difference is `O(δ)` and tiny (dense `gap_Φ/gap_𝓛 ∈
[1.0003, 1.0007]`), but extracting it by an INDEPENDENT matrix-free channel
eigensolve is fragile: `Φ_δ`'s slow eigenVECTOR is `~1e-5`-orthogonal to every
physical/random seed, so Krylov picks the wrong mode (11–47 % inflation,
krylovdim-independent). The generator form fixes the eigenvalue clustering but NOT
the seed-overlap (same eigenvectors). PT sidesteps both: it applies `M` to the
REFERENCE's slow modes, which ARE robust (`𝓛` is KMS-normal — its `|+⟩` Pass-1 gap
equals the true operator gap to `≤1.3e-6`, qf-e4z.43), with no channel eigensolve,
no `O(δ)` cancellation of two `O(gap)` eigenvalues, and no wrong-mode risk. PT IS
"the generator-form channel solve seeded by the ideal slow mode, to first order".

# Basis contract (correctness crux)
`R_k`, `L_k` are returned by the reference arm's Krylov decomposition in the
reference WORKING basis, rotated ONCE to the COMMON basis via `_rotate_out` (HS
biorthonormality is unitarily invariant, so it survives the rotation). When applying
`G_test` / `G_ref` to a common-basis `R_k`, each arm rotates IN with its OWN basis,
applies its generator, and rotates OUT — so a reference in `ham.eigvecs`
(EnergyDomain) and a test in `trotter.eigvecs` (TrotterDomain) are reconciled
correctly. Identical arms (or the same operator in two different working bases) give
`ε_slow = 0` to machine precision — the rotation-bug canary.

# Arguments
- `ref_arm::PropagatorArm`: source of the slow eigenbasis (typically the ideal
  `lindbladian_arm`; must be KMS-normal for the modes to be robust).
- `test_arm::PropagatorArm`: the propagator whose generator is compared (typically
  `channel_arm`; may live in a DIFFERENT working basis).
- `rho_seed::AbstractMatrix`: Arnoldi seed for the reference decomposition, in the
  COMMON basis. It must OVERLAP the reference's slow modes — pass the canonical
  `|+⟩⟨+|^⊗N` (validated to capture the gap for the disordered-Heisenberg CKG `𝓛`,
  qf-e4z.43). `I/d` is symmetry-protected and silently misses the (traceless) gap
  mode (qf-8fr), so it is NOT a valid seed here.

# Keywords
- `num_slow_modes::Integer = 1`: `K`, the number of slowest NON-stationary modes
  retained. `K = 1` is the deliverable; `K > 1` only if `max_antiherm_frac` stays
  small (see the Hermiticity section). `krylovdim` must exceed `K + 1`.
- `include_stationary::Bool = false`: prepend the `k = 0` steady mode (gives a
  `(K+1)×(K+1)` block whose lowest rung is the fixed-point bias).
- `krylovdim::Integer = 60`: Arnoldi subspace size for the reference decomposition
  (must be `> K + 1`; the slow modes converge well below this, qf-5vx).
- `tol::Real = 1e-10`: forwarded to the reference decomposition.
- `antiherm_tol::Real = 1e-6`: threshold on `max_antiherm_frac` above which a
  `@warn` fires (retained modes non-Hermitian ⇒ a degenerate cluster is cut ⇒
  `eps_slow` is corrupted; trust only `K = 1`).
- `pt_spacing_frac::Real = 0.5`: threshold on `eps_slow / |λ₂−λ₃|` above which a
  `@warn` fires, FOR THE `K = 1` DELIVERABLE ONLY (qf-dee #4). First-order PT for the
  gap eigenvalue needs the perturbation `≪` the gap mode's spacing to its nearest
  non-stationary neighbour `|λ₂−λ₃|`, not merely `≪ |λ₂|`; this is the spacing
  PRECONDITION (complementing the `max_antiherm_frac` cut-cluster guard, which catches
  the symptom). The canonical n=3 deliverable sits at `≈ 0.11` (dense-validated, so the
  default 0.5 is silent there); `→ 1` flags genuine near-degeneracy. For `K > 1` block
  validity is governed by `max_antiherm_frac` and no spacing warn fires.

# Returns
NamedTuple with `eps_slow` (`‖M_jk‖₂`, rate units), `eps_slow_rel_gap`
(`eps_slow / |λ₂|`, the dimensionless headline), `eps_slow_rel_block`
(`eps_slow / ‖ΠLΠ‖₂`), `lambda_neighbor_spacing` (`|λ₂−λ₃|`, the gap mode's spacing to
its nearest non-stationary neighbour; `NaN` if `<3` modes resolved), `eps_slow_rel_neighbor`
(`eps_slow / |λ₂−λ₃|`, the first-order-PT validity ratio — `≪ 1` to trust the gap shift),
`M` (the `K×K` mismatch block), `M_diagonal` (per-mode
first-order rate shifts), `eigenvalues_ref` (retained reference `λ`), `gap_ref`
(`|Re λ₂|`), `block_norm_ref` (`max|λ|` over retained modes), `num_slow_modes`,
`include_stationary`, `max_antiherm_frac` (worst post-phase-fix anti-Hermitian
fraction of the retained modes — `~1e-13` ⇒ trust, `O(1)` ⇒ cut cluster, untrusted),
`ref_gen_residual` (`max|⟨L_j|G_ref|R_k⟩ − λ_k δ_jk|`, a near-machine-precision
convergence canary for the reference modes), `matvecs` (`(ref_decomp, test_gen,
ref_gen)` counts), `converged` (reference Arnoldi did not break down), and `labels`
(`(ref, test)`).
"""
function slow_subspace_generator_distance(
    ref_arm::PropagatorArm, test_arm::PropagatorArm,
    rho_seed::AbstractMatrix;
    num_slow_modes::Integer = 1,
    include_stationary::Bool = false,
    krylovdim::Integer = 60,
    tol::Real = 1e-10,
    antiherm_tol::Real = 1e-6,
    pt_spacing_frac::Real = 0.5,
)
    d = size(rho_seed, 1)
    size(rho_seed, 2) == d || throw(ArgumentError("rho_seed must be square"))
    K = Int(num_slow_modes)
    K >= 1 || throw(ArgumentError("num_slow_modes must be ≥ 1 (got $K)"))
    Int(krylovdim) > K + 1 || throw(ArgumentError(
        "krylovdim ($krylovdim) must exceed num_slow_modes + 1 ($(K + 1)) to resolve " *
        "the steady + $K slow modes"))

    # 1. Reference slow eigenpairs: a single forward Arnoldi on the reference
    #    GENERATOR, seeded by rho_seed rotated into the reference working basis.
    #    `_krylov_spectral_decomposition` returns biorthonormal R_modes / L_modes
    #    sorted with the steady state at index 1, the gap mode at index 2, … .
    gen_ref = _arm_generator(ref_arm)
    seed_ref_work = Matrix{ComplexF64}(_rotate_in(ref_arm, rho_seed))
    decomp = _krylov_spectral_decomposition(
        gen_ref, seed_ref_work, d;
        krylovdim = Int(krylovdim), tol = tol, sort_mode = :lindbladian)

    m = length(decomp.R_modes)
    m >= K + 1 || throw(ArgumentError(
        "reference Krylov decomposition returned only $m modes; need ≥ $(K + 1) " *
        "(steady + K=$K slow). Increase krylovdim (got $krylovdim)."))

    # Decomposition idx 1 = steady (design k=0), idx 2 = gap (design k=1). Retain
    # the K slowest NON-stationary modes (idx 2:K+1), optionally prepend the steady.
    idxs = include_stationary ? collect(1:(K + 1)) : collect(2:(K + 1))
    Ksel = length(idxs)

    # Rotate retained modes ONCE to the COMMON basis (HS biorthonormality is
    # unitarily invariant ⇒ ⟨L_j,R_k⟩ = δ_jk survives the rotation).
    R_common = [Matrix{ComplexF64}(_rotate_out(ref_arm, decomp.R_modes[i])) for i in idxs]
    L_common = [Matrix{ComplexF64}(_rotate_out(ref_arm, decomp.L_modes[i])) for i in idxs]
    Λ = decomp.eigenvalues[idxs]

    # Phase-fix each retained mode to be Hermitian: a :channel test arm is ℝ-linear
    # (Hermiticity-preserving + defensive hermitize), so it acts faithfully only on
    # HERMITIAN inputs. A KMS-normal reference has real slow spectrum ⇒ each
    # non-degenerate eigenmode is e^{iθ}H; rotate by α (|α|=1) so αR_k is Hermitian
    # and counter-rotate L_k (⟨L_j,R_k⟩ = δ_jk is phase-invariant: both ×α). The gap
    # mode becomes exactly Hermitian (K=1 clean). A cut degenerate cluster cannot be
    # Hermitized by any single phase — `max_antiherm_frac` then stays O(1) and warns.
    antiherm_fracs = Vector{Float64}(undef, Ksel)
    for kk in 1:Ksel
        α = _hermitizing_phase(R_common[kk])
        R_common[kk] = α .* R_common[kk]
        L_common[kk] = α .* L_common[kk]
        antiherm_fracs[kk] = norm((R_common[kk] .- R_common[kk]') ./ 2) / norm(R_common[kk])
    end
    max_antiherm_frac = maximum(antiherm_fracs)
    # The corruption is specific to an ℝ-linear (Hermiticity-hermitizing) TEST arm —
    # a :lindbladian test arm is ℂ-linear, so non-Hermitian modes are harmless there.
    if test_arm.kind === :channel && max_antiherm_frac > antiherm_tol
        @warn "slow_subspace_generator_distance: retained reference mode(s) remain " *
              "non-Hermitian after phase-fixing (max anti-Hermitian fraction " *
              "$(round(max_antiherm_frac, sigdigits = 3))). The test arm's matvec is " *
              "Hermiticity-preserving (ℝ-linear), so the projected mismatch is " *
              "corrupted on non-Hermitian modes — this usually means num_slow_modes = " *
              "$K cuts a degenerate cluster of the reference. Trust only the gap shift " *
              "(num_slow_modes = 1)." maxlog = 1
    end

    # 2. M_jk = ⟨L_j | (G_test − G_ref) | R_k⟩_HS. Apply each generator to the
    #    common-basis R_k via its OWN arm (rotate in → generate → rotate out), then
    #    contract with L_j. ⟨L_j|G_ref|R_k⟩ = λ_k δ_jk exactly (Arnoldi relation),
    #    captured separately as a convergence diagnostic.
    gen_test = _arm_generator(test_arm)
    M = zeros(ComplexF64, Ksel, Ksel)
    Gref_proj = zeros(ComplexF64, Ksel, Ksel)
    out_t = Matrix{ComplexF64}(undef, d, d)
    out_r = Matrix{ComplexF64}(undef, d, d)
    mv_test = 0
    mv_ref_gen = 0
    for (kk, Rk) in enumerate(R_common)
        gen_test(out_t, Matrix{ComplexF64}(_rotate_in(test_arm, Rk)))
        mv_test += 1
        GtestRk = _rotate_out(test_arm, out_t)
        gen_ref(out_r, Matrix{ComplexF64}(_rotate_in(ref_arm, Rk)))
        mv_ref_gen += 1
        GrefRk = _rotate_out(ref_arm, out_r)
        MRk = GtestRk .- GrefRk
        for (jj, Lj) in enumerate(L_common)
            M[jj, kk] = dot(vec(Lj), vec(MRk))
            Gref_proj[jj, kk] = dot(vec(Lj), vec(GrefRk))
        end
    end

    # 3. Scalars. The gap mode is the slowest non-stationary retained mode: it sits
    #    at retained-index 1 (steady excluded) or 2 (steady prepended).
    eps_slow = opnorm(M, 2)
    gap_idx = include_stationary ? 2 : 1
    gap_ref = abs(real(Λ[gap_idx]))                 # gap RATE = |Re λ₂| (krylov_spectral_gap convention)
    block_norm = maximum(abs.(Λ))                   # ‖ΠLΠ‖₂ = max|λ| (diag in eigenbasis; |λ| not |Re λ|, valid for a non-normal reference too)
    ref_gen_residual = maximum(abs.(Gref_proj .- Diagonal(Λ)))
    # ⟨L_j|G_ref|R_k⟩ = λ_k δ_jk holds to the Arnoldi-relation residual when the
    # reference modes are converged. A large residual means the reference modes are
    # NOT well-resolved (e.g. a non-KMS-normal / channel arm used as REFERENCE, whose
    # slow eigenvector is seed-orthogonal, qf-9lp) ⇒ the whole measure is unreliable.
    if ref_gen_residual > sqrt(tol)
        @warn "slow_subspace_generator_distance: reference generator is not diagonal " *
              "on its own slow modes (max|⟨L_j|G_ref|R_k⟩ − λ_k δ_jk| = " *
              "$(round(ref_gen_residual, sigdigits = 3)) > √tol). The reference arm's " *
              "slow modes are not well-resolved — use a KMS-normal reference (the ideal " *
              "Lindbladian), not a channel, and a seed overlapping its slow modes." maxlog = 1
    end

    # qf-dee #4: first-order PT for the gap eigenvalue is valid only when the
    # perturbation is ≪ the gap mode's spacing to its NEAREST non-stationary NEIGHBOR
    # |λ₂−λ₃| — not merely ≪ |λ₂|. (`max_antiherm_frac` catches a CUT degenerate
    # cluster, i.e. the symptom; this checks the spacing PRECONDITION.) The reference
    # decomposition is sorted steady=1, gap=2, neighbor=3, so the neighbor spacing is
    # |λ₂−λ₃| from the FULL spectrum (independent of which modes were retained).
    λ_neighbor_spacing = length(decomp.eigenvalues) >= 3 ?
        abs(decomp.eigenvalues[2] - decomp.eigenvalues[3]) : NaN
    # Warn only for the K=1 gap-shift deliverable: there λ₃ IS the gap mode's first
    # excluded neighbour and ε_slow IS the gap shift, so the check is exactly "is the
    # first-order gap shift ≪ the gap-to-neighbour spacing". For K>1 the retained block
    # may already include λ₃ and block validity is governed by max_antiherm_frac instead.
    if K == 1 && isfinite(λ_neighbor_spacing) && λ_neighbor_spacing > 0 &&
       eps_slow > pt_spacing_frac * λ_neighbor_spacing
        @warn "slow_subspace_generator_distance: perturbation ε_slow = " *
              "$(round(eps_slow, sigdigits = 3)) is not ≪ the gap mode's spacing to its " *
              "nearest neighbor |λ₂−λ₃| = $(round(λ_neighbor_spacing, sigdigits = 3)) " *
              "(ratio $(round(eps_slow / λ_neighbor_spacing, sigdigits = 3)) > " *
              "pt_spacing_frac = $(pt_spacing_frac)). First-order PT for the gap " *
              "eigenvalue needs the perturbation ≪ the neighbor spacing, not just ≪ " *
              "|λ₂|; the gap shift may be unreliable (gap and neighbor can mix)." maxlog = 1
    end

    return (
        eps_slow            = eps_slow,
        eps_slow_rel_gap    = gap_ref > 0 ? eps_slow / gap_ref : NaN,
        eps_slow_rel_block  = block_norm > 0 ? eps_slow / block_norm : NaN,
        lambda_neighbor_spacing = λ_neighbor_spacing,
        eps_slow_rel_neighbor   = (isfinite(λ_neighbor_spacing) && λ_neighbor_spacing > 0) ?
                                  eps_slow / λ_neighbor_spacing : NaN,
        M                   = M,
        M_diagonal          = diag(M),
        eigenvalues_ref     = Λ,
        gap_ref             = gap_ref,
        block_norm_ref      = block_norm,
        num_slow_modes      = K,
        include_stationary  = include_stationary,
        max_antiherm_frac   = max_antiherm_frac,
        ref_gen_residual    = ref_gen_residual,
        matvecs             = (ref_decomp = decomp.matvec_count, test_gen = mv_test, ref_gen = mv_ref_gen),
        converged           = decomp.converged,
        labels              = (ref_arm.label, test_arm.label),
    )
end

# ---------------------------------------------------------------------------
# Convenience arm constructors (reuse Workspace + the canonical matvecs)
# ---------------------------------------------------------------------------

"""
    lindbladian_arm(config::Config{Lindbladian}, hamiltonian, jumps;
                    basis=hamiltonian.eigvecs, workspace=nothing,
                    label="ideal L") -> PropagatorArm

A `:lindbladian` [`PropagatorArm`](@ref) wrapping `apply_lindbladian!` for a
`Config{Lindbladian}` (EnergyDomain / BohrDomain / …). `jumps` must be expressed in
the matvec's working basis (e.g. `_jumps_in_basis(ham, n, ham.eigvecs)` for
EnergyDomain). `basis` is the unitary mapping that working basis to the common
basis — for EnergyDomain it is `hamiltonian.eigvecs`. Builds (or reuses) a
`Workspace`; the coherent term is ON (canonical KMS Lindbladian).
"""
function lindbladian_arm(config::Config{Lindbladian}, hamiltonian::HamHam,
                         jumps::Vector{JumpOp};
                         basis::Union{Nothing, AbstractMatrix} = hamiltonian.eigvecs,
                         workspace::Union{Nothing, Workspace} = nothing,
                         label::AbstractString = "ideal L")
    ws = workspace === nothing ? Workspace(config, hamiltonian, jumps) : workspace
    apply! = let ws = ws, config = config, ham = hamiltonian
        (out::AbstractMatrix, x::AbstractMatrix) -> begin
            apply_lindbladian!(ws, Matrix{ComplexF64}(x), config, ham)
            copyto!(out, ws.scratch.rho_out)
            return out
        end
    end
    return PropagatorArm(apply!; kind = :lindbladian, basis = basis, label = label)
end

"""
    channel_arm(config::Config{Thermalize}, hamiltonian, jumps, trotter;
                basis=trotter.eigvecs, workspace=nothing,
                label="implemented Φ_δ") -> PropagatorArm

A `:channel` [`PropagatorArm`](@ref) wrapping the faithful `apply_delta_channel!`
for a `Config{Thermalize}` (the deployed Strang/GQSP/OFT/weak-measurement channel,
bit-identical to `run_thermalize :sweep`). For `TrotterDomain` pass the `trotter`
object and build `jumps = _jumps_in_basis(ham, n, trotter.eigvecs)`; `basis` is then
`trotter.eigvecs`. `delta` is read from `config.delta`.
"""
function channel_arm(config::Config{Thermalize}, hamiltonian::HamHam,
                     jumps::Vector{JumpOp},
                     trotter::Union{Nothing, AbstractTrotter} = nothing;
                     basis::Union{Nothing, AbstractMatrix} =
                         trotter === nothing ? nothing : trotter.eigvecs,
                     workspace::Union{Nothing, Workspace} = nothing,
                     label::AbstractString = "implemented Φ_δ")
    config.jump_selection === :sweep || throw(ArgumentError(
        "channel_arm requires config.jump_selection = :sweep (got :$(config.jump_selection))"))
    config.delta === nothing && throw(ArgumentError("channel_arm requires config.delta to be set"))
    ws = workspace === nothing ? Workspace(config, hamiltonian, jumps; trotter = trotter) : workspace
    apply! = let ws = ws, config = config, ham = hamiltonian
        (out::AbstractMatrix, x::AbstractMatrix) -> begin
            apply_delta_channel!(ws, Matrix{ComplexF64}(x), config, ham)
            copyto!(out, ws.scratch.rho_next)
            return out
        end
    end
    return PropagatorArm(apply!; kind = :channel, delta = float(config.delta),
                         basis = basis, label = label)
end

# ---------------------------------------------------------------------------
# qf-e4z.50: anti-Hermitian norm of the quantum-discriminant of a generator —
# the channel-vs-Lindbladian KMS-detailed-balance-violation error axis.
# ---------------------------------------------------------------------------
#
# The KMS quantum discriminant of a generator G w.r.t. the Gibbs state σ is
#     D(G)(X) = σ^{-1/4} G(σ^{1/4} X σ^{1/4}) σ^{-1/4}                (discriminant.jl).
# KMS detailed balance ⇔ D(G) is HS-self-adjoint ⇔ its anti-Hermitian part
#     A = (D − D†_HS)/2
# vanishes. For the ideal KMS-DB Lindbladian 𝓛, A = 0 exactly (only quadrature
# residual survives). For the IMPLEMENTED channel Φ_δ, the effective generator
# G_eff = (Φ_δ − I)/δ breaks detailed balance through Trotter/GQSP/δ-splitting and
# ‖A‖_op > 0 measures that violation (units of rate, like the spectral gap). This
# is the 3rd error axis alongside `slow_subspace_generator_distance` (gap
# distortion) and `fixed_point_gibbs_distance` (Gibbs drift).
#
# Matrix-free recipe (scales to n=9/10 like the other two axes):
#  • A is anti-Hermitian (A†_HS = −A), so ‖A‖_op = its largest singular value, got
#    from `hs_operator_norm_krylov` (GKL) with forward `A` and adjoint `−A` — both
#    built from the FORWARD generator G and its HS-ADJOINT G† via `apply_discriminant!`
#    (no d²×d² superoperator is ever formed).
#  • Work in the basis where σ is DIAGONAL (the Hamiltonian eigenbasis), so the σ^{±1/4}
#    conjugations are length-d diagonal scalings. Conditioning κ_σ = (σ_max/σ_min)^{1/2}
#    = e^{β_alg·ΔE/2}; benign at the operating point (≈5 at n=3 / β_phys=0.5).

"""
    discriminant_antiherm_norm(gen!, gen_adj!, sigma_quarter, sigma_inv_quarter, d;
                               krylovdim=30, tol=1e-12, max_retries=3,
                               compute_discriminant_norm=false) -> NamedTuple

Matrix-free operator-2 norm of the anti-Hermitian part `A = (D − D†_HS)/2` of the
KMS quantum discriminant `D(G)` of a generator `G`, w.r.t. a Gibbs state `σ` that is
DIAGONAL in the working basis (passed as the length-`d` fractional-power vectors
`sigma_quarter = σ^{1/4}`, `sigma_inv_quarter = σ^{-1/4}` from
[`gibbs_fractional_powers`](@ref)).

`gen!(out, X)` writes the forward generator action `G(X)` into `out`; `gen_adj!(out, X)`
writes the HS-adjoint `G†(X)`. Both act on `d×d` matrices in the working basis. The
discriminant and its HS-adjoint are
`D(X)  = σ^{-1/4} G(σ^{1/4} X σ^{1/4}) σ^{-1/4}` and
`D†(X) = σ^{1/4} G†(σ^{-1/4} X σ^{-1/4}) σ^{1/4}`, assembled via
[`apply_discriminant!`](@ref); `A = (D − D†)/2` is anti-Hermitian so its largest
singular value (the returned `antiherm_norm`) is `‖A‖_op = ‖D − D†‖_op / 2`,
computed by [`hs_operator_norm_krylov`](@ref).

For a KMS-DB `G` (e.g. the ideal Lindbladian) `antiherm_norm ≈ 0` — the GKL
self-consistency floor handling returns the noise-floor lower bound there.

# Returns
NamedTuple with `antiherm_norm` (`‖A‖_op`), `discriminant_norm` (`‖D‖_op`, only when
`compute_discriminant_norm`, else `NaN`), `relative` (`antiherm_norm/discriminant_norm`
or `NaN`), and `conditioning` (`κ_σ = (σ_max/σ_min)^{1/2}`).
"""
function discriminant_antiherm_norm(
    gen!::FG, gen_adj!::FA,
    sigma_quarter::AbstractVector{<:Real},
    sigma_inv_quarter::AbstractVector{<:Real},
    d::Integer;
    krylovdim::Integer = 30, tol::Real = 1e-12, max_retries::Integer = 3,
    compute_discriminant_norm::Bool = false,
) where {FG, FA}
    sq    = collect(Float64, sigma_quarter)
    sqinv = collect(Float64, sigma_inv_quarter)
    length(sq) == d || throw(ArgumentError("sigma_quarter length $(length(sq)) ≠ d=$d"))

    bufs_D  = DiscriminantBuffers(d)
    bufs_Dd = DiscriminantBuffers(d)
    bufD  = Matrix{ComplexF64}(undef, d, d)
    bufDd = Matrix{ComplexF64}(undef, d, d)

    # D(X) = σ^{-1/4} G(σ^{1/4} X σ^{1/4}) σ^{-1/4}
    D!(out, X) = apply_discriminant!(out, Matrix{ComplexF64}(X), gen!, sq, sqinv, bufs_D)
    # D†(X) = σ^{1/4} G†(σ^{-1/4} X σ^{-1/4}) σ^{1/4}  (swap the two scale vectors)
    Dadj!(out, X) = apply_discriminant!(out, Matrix{ComplexF64}(X), gen_adj!, sqinv, sq, bufs_Dd)

    # A(X) = (D(X) − D†(X))/2 ; A is anti-Hermitian ⇒ A†_HS = −A.
    A!(out, X) = begin
        D!(bufD, X); Dadj!(bufDd, X)
        @inbounds @. out = (bufD - bufDd) / 2
        out
    end
    negA!(out, X) = (A!(out, X); @inbounds @. out = -out; out)

    antiherm = hs_operator_norm_krylov(A!, negA!, d;
                                       tol = tol, krylovdim = Int(krylovdim),
                                       max_retries = Int(max_retries))

    discr_norm = NaN
    if compute_discriminant_norm
        discr_norm = hs_operator_norm_krylov(D!, Dadj!, d;
                                             tol = tol, krylovdim = Int(krylovdim),
                                             max_retries = Int(max_retries))
    end

    # κ_σ = (σ_max/σ_min)^{1/2} = (max σ^{1/4} / min σ^{1/4})^2
    κσ = (maximum(sq) / minimum(sq))^2
    return (
        antiherm_norm     = antiherm,
        discriminant_norm = discr_norm,
        relative          = compute_discriminant_norm ? antiherm / discr_norm : NaN,
        conditioning      = κσ,
    )
end

# σ^{±1/4} diagonal (Hamiltonian eigenbasis) from a Gibbs spectrum exp(-β_alg E)/Z.
function _gibbs_quarter_powers(hamiltonian::HamHam, beta_alg::Real)
    g = exp.(-float(beta_alg) .* hamiltonian.eigvals)
    g ./= sum(g)
    powers = gibbs_fractional_powers(Hermitian(Matrix(Diagonal(ComplexF64.(g)))))
    return powers.sigma_quarter, powers.sigma_inv_quarter
end

"""
    channel_discriminant_antiherm_norm(config::Config{Thermalize}, hamiltonian, jumps,
                                       trotter=nothing; krylovdim=30, tol=1e-12,
                                       max_retries=3, workspace=nothing,
                                       compute_discriminant_norm=false) -> NamedTuple

`‖A[D(G_eff)]‖_op` for the implemented CKG channel `Φ_δ`, where the effective
generator is `G_eff = (Φ_δ − I)/δ` — the channel-side KMS-detailed-balance-violation
error (qf-e4z.50), matrix-free to n=9/10. The genuine ℂ-linear channel and its
HS-adjoint are used ([`apply_delta_channel!`](@ref) / [`apply_adjoint_delta_channel!`](@ref),
both `hermitize=false`); the discriminant is taken in the Hamiltonian eigenbasis
(σ diagonal), with the TrotterDomain channel rotated in/out by
`W = trotter.eigvecs' · ham.eigvecs`. Hermitian (Pauli) jumps only.

`jumps` must be expressed in the channel's working basis (`_jumps_in_basis(ham, n,
trotter.eigvecs)` for TrotterDomain). Pass a prebuilt `workspace` to skip the
`Workspace` build. See [`discriminant_antiherm_norm`](@ref) for the returned fields.
"""
function channel_discriminant_antiherm_norm(
    config::Config{Thermalize}, hamiltonian::HamHam, jumps::Vector{JumpOp},
    trotter::Union{Nothing, AbstractTrotter} = nothing;
    krylovdim::Integer = 30, tol::Real = 1e-12, max_retries::Integer = 3,
    workspace::Union{Nothing, Workspace} = nothing,
    compute_discriminant_norm::Bool = false,
)
    config.delta === nothing && throw(ArgumentError(
        "channel_discriminant_antiherm_norm requires config.delta to be set"))
    config.jump_selection === :sweep || throw(ArgumentError(
        "channel_discriminant_antiherm_norm requires config.jump_selection = :sweep " *
        "(got :$(config.jump_selection))"))
    ws = workspace === nothing ?
        Workspace(config, hamiltonian, jumps; trotter = trotter) : workspace
    d = size(hamiltonian.data, 1)
    δ = float(config.delta)
    sq, sqinv = _gibbs_quarter_powers(hamiltonian, config.beta)

    # Rotation ham eigenbasis → channel working (trotter) basis: X_trott = W X_ham W'.
    Wrot = trotter === nothing ?
        Matrix{ComplexF64}(I, d, d) :
        Matrix{ComplexF64}(trotter.eigvecs)' * Matrix{ComplexF64}(hamiltonian.eigvecs)
    Xtr = Matrix{ComplexF64}(undef, d, d)
    rot = Matrix{ComplexF64}(undef, d, d)

    chan_gen!(out, Xham) = begin
        Xh = Matrix{ComplexF64}(Xham)
        mul!(rot, Wrot, Xh); mul!(Xtr, rot, Wrot')              # ham → trotter
        apply_delta_channel!(ws, Xtr, config, hamiltonian; hermitize = false)
        mul!(rot, Wrot', ws.scratch.rho_next); mul!(out, rot, Wrot)  # trotter → ham
        @inbounds @. out = (out - Xh) / δ                       # G_eff = (Φ−I)/δ
        out
    end
    chan_gen_adj!(out, Xham) = begin
        Xh = Matrix{ComplexF64}(Xham)
        mul!(rot, Wrot, Xh); mul!(Xtr, rot, Wrot')
        apply_adjoint_delta_channel!(ws, Xtr, config, hamiltonian; hermitize = false)
        mul!(rot, Wrot', ws.scratch.rho_next); mul!(out, rot, Wrot)
        @inbounds @. out = (out - Xh) / δ                       # G_eff† = (Φ†−I)/δ
        out
    end

    return discriminant_antiherm_norm(chan_gen!, chan_gen_adj!, sq, sqinv, d;
        krylovdim = krylovdim, tol = tol, max_retries = max_retries,
        compute_discriminant_norm = compute_discriminant_norm)
end

"""
    lindbladian_discriminant_antiherm_norm(config::Config{Lindbladian}, hamiltonian,
                                           jumps; krylovdim=30, tol=1e-12,
                                           max_retries=3, workspace=nothing,
                                           compute_discriminant_norm=false) -> NamedTuple

`‖A[D(𝓛)]‖_op` for the ideal KMS-DB Lindbladian — the BASELINE / noise-floor of the
[`channel_discriminant_antiherm_norm`](@ref) error axis. For a KMS-DB `𝓛` this is the
residual detailed-balance violation from EnergyDomain quadrature only (≈ machine
precision at r_D=8), so it bounds the floor against which the channel value is read.
`jumps` in the Lindbladian working basis (`_jumps_in_basis(ham, n, ham.eigvecs)` for
EnergyDomain). See [`discriminant_antiherm_norm`](@ref) for the returned fields.
"""
function lindbladian_discriminant_antiherm_norm(
    config::Config{Lindbladian}, hamiltonian::HamHam, jumps::Vector{JumpOp};
    krylovdim::Integer = 30, tol::Real = 1e-12, max_retries::Integer = 3,
    workspace::Union{Nothing, Workspace} = nothing,
    compute_discriminant_norm::Bool = false,
)
    ws = workspace === nothing ? Workspace(config, hamiltonian, jumps) : workspace
    d = size(hamiltonian.data, 1)
    sq, sqinv = _gibbs_quarter_powers(hamiltonian, config.beta)

    lind_gen!(out, X) = (apply_lindbladian!(ws, Matrix{ComplexF64}(X), config, hamiltonian);
                         copyto!(out, ws.scratch.rho_out); out)
    lind_gen_adj!(out, X) = (apply_adjoint_lindbladian!(ws, Matrix{ComplexF64}(X), config, hamiltonian);
                             copyto!(out, ws.scratch.rho_out); out)

    return discriminant_antiherm_norm(lind_gen!, lind_gen_adj!, sq, sqinv, d;
        krylovdim = krylovdim, tol = tol, max_retries = max_retries,
        compute_discriminant_norm = compute_discriminant_norm)
end
