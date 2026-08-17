# Computational basis: site 1 is the leftmost Kronecker factor and spins are
# `+1` for up and `-1` for down.
const PX = ComplexF64[0 1; 1 0]
const PY = ComplexF64[0 -im; im 0]
const PZ = ComplexF64[1 0; 0 -1]
const I2 = ComplexF64[1 0; 0 1]

"Kronecker-pad a single-site operator O at site q on n sites (site 1 leftmost)."
function one_site(O::Matrix{ComplexF64}, q::Int, n::Int)
    ops = [k == q ? O : I2 for k in 1:n]
    return foldl(kron, ops)
end

"Kronecker-pad a two-site operator O⊗O at sites (a,b) on n sites."
function two_site(O::Matrix{ComplexF64}, a::Int, b::Int, n::Int)
    @assert a != b
    ops = [(k == a || k == b) ? O : I2 for k in 1:n]
    return foldl(kron, ops)
end

# Model data are stored in physical units; spins are integer signs.

# Operator kind codes
const ID      = Int8(0)
const EXCH    = Int8(1)   # Heisenberg exchange (bond); diagonal H1 or off-diag H2
const FIELD   = Int8(2)   # Heisenberg longitudinal field (site), diagonal
const ZZDIS   = Int8(3)   # Heisenberg ZZ disorder (bond), diagonal
const ISING   = Int8(4)   # TFIM Ising bond, diagonal
const FIELDC  = Int8(5)   # TFIM field constant (site), diagonal
const FIELDX  = Int8(6)   # TFIM field flip (site), off-diagonal
const ZDIS    = Int8(7)   # TFIM Z disorder (site), diagonal
const ZZDIS_T = Int8(8)   # TFIM ZZ disorder (bond), diagonal

struct Op
    kind::Int8
    idx::Int32       # bond index (into model.bonds) or site index
    offdiag::Bool
end
const IDOP = Op(ID, Int32(0), false)

is_id(op::Op)      = op.kind == ID
is_offdiag(op::Op) = (op.kind == EXCH && op.offdiag) || op.kind == FIELDX
is_diag(op::Op)    = !is_id(op) && !is_offdiag(op)


# Convert a package Hamiltonian back to the physical frame.
function fixture_H_phys(raw)
    d = size(raw.matrix, 1)
    return (raw.matrix .- raw.shift .* Matrix{ComplexF64}(I, d, d)) .* raw.rescaling_factor
end

"""
    sse_reconstruction_error(model, raw) -> Float64

Return the maximum entrywise mismatch between an SSE model and a package
Hamiltonian after undoing the algorithm-frame shift and rescaling.
"""
function sse_reconstruction_error end

# Exact small-system reference in the physical frame.
struct ExactRef
    energy::Float64
    mz2::Float64
    zz::Dict{Tuple{Int,Int},Float64}     # <Z_i Z_j>
end


"""
    sse_exact_reference(model, beta_phys; pairs) -> ExactRef

Compute exact thermal energy, magnetisation, and requested correlations.

# Arguments
- `beta_phys`: Physical inverse temperature.
- `pairs`: Site pairs for which to return `Z_i Z_j` correlations.

# Returns
An `ExactRef` containing dense Gibbs-state observables.
"""
function sse_exact_reference end

function _exact_reference(H::Matrix{ComplexF64}, n::Int, beta_phys::Float64; pairs::Vector{Tuple{Int,Int}})
    validate_observable_inputs(n, beta_phys, pairs)
    Hh = Hermitian(H)
    vals, vecs = eigen(Hh)
    w = exp.(-beta_phys .* (vals .- minimum(vals)))
    w ./= sum(w)
    rho = vecs * Diagonal(ComplexF64.(w)) * vecs'
    energy = real(tr(rho * H))
    Mz = zeros(ComplexF64, size(H))
    for q in 1:n
        Mz .+= one_site(PZ, q, n)
    end
    mz2 = real(tr(rho * (Mz * Mz))) / n^2
    zzd = Dict{Tuple{Int,Int},Float64}()
    for (i, j) in pairs
        zzd[(i, j)] = real(tr(rho * two_site(PZ, i, j, n)))
    end
    return ExactRef(energy, mz2, zzd)
end

# Statistical estimators.
"Block-jackknife mean and error of the ratio mean(X)/mean(Y). For Y≡1, plain jackknife of mean(X)."
function jackknife_ratio(X::Vector{Float64}, Y::Vector{Float64}; nblocks::Int=40)
    N = length(X)
    length(Y) == N || throw(DimensionMismatch("X and Y must have equal lengths."))
    N >= 2 || throw(ArgumentError("jackknife_ratio requires at least two samples."))
    2 <= nblocks <= N || throw(ArgumentError("nblocks must lie between 2 and the sample count."))
    nb = nblocks
    edges = round.(Int, range(0, N; length=nb + 1))
    sumX = sum(X); sumY = sum(Y)
    iszero(sumY) && throw(ArgumentError("ratio estimator has zero total denominator."))
    Rb = Float64[]
    for b in 1:nb
        lo, hi = edges[b] + 1, edges[b + 1]
        sx = sum(@view X[lo:hi]); sy = sum(@view Y[lo:hi])
        denominator = sumY - sy
        iszero(denominator) && throw(ArgumentError(
            "ratio estimator has a zero leave-one-block-out denominator."))
        push!(Rb, (sumX - sx) / denominator)
    end
    Rjk = mean(Rb)
    err = sqrt((nb - 1) / nb * sum((Rb .- Rjk) .^ 2))
    Rfull = sumX / sumY
    return Rfull, err
end

"""
    tau_int(series; c=6.0, max_pairs=20_000_000, label="series")

Estimate the integrated autocorrelation time using a finite Sokal window.
Direct autocovariances are capped at `max_pairs` multiply-add pairs, preventing
the unresolved-window case from degrading to quadratic cost. A warning is
emitted when no self-consistent window is found. Returns `(tau, err, window)`.
"""
function tau_int(
    series::Vector{Float64};
    c::Float64=6.0,
    max_pairs::Int=20_000_000,
    label::AbstractString="series",
)
    N = length(series)
    N >= 2 || throw(ArgumentError("tau_int requires at least two samples."))
    isfinite(c) && c > 0 || throw(ArgumentError("c must be finite and > 0."))
    max_pairs > 0 || throw(ArgumentError("max_pairs must be positive."))
    all(isfinite, series) || throw(ArgumentError("autocorrelation samples must be finite."))

    m = mean(series)
    centered = series .- m
    v = dot(centered, centered) / N
    if v ≤ 0
        return 0.5, 0.0, 0
    end

    max_lag = min(N - 1, max(1, fld(max_pairs, N)))
    tau = 0.5
    W = max_lag
    resolved = false
    for t in 1:max_lag
        rho_t = dot(@view(centered[1:(N - t)]), @view(centered[(t + 1):N])) /
            ((N - t) * v)
        tau += rho_t
        if t >= c * max(tau, 0.5)
            W = t
            resolved = true
            break
        end
    end
    tau = max(tau, 0.5)
    if !resolved
        @warn(
            "No self-consistent autocorrelation window was resolved; the reported tau is truncated.",
            label=label,
            N=N,
            max_lag=max_lag,
            tau=tau,
        )
    end
    err = tau * sqrt(2 * (2 * W + 1) / N)     # Madras–Sokal
    return tau, err, W
end

"Choose blocks much longer than the pilot autocorrelation time."
function blocking_nblocks(N::Int, tau::Float64; label::AbstractString="SSE observables")
    N >= 2 || throw(ArgumentError("blocking requires at least two samples."))
    isfinite(tau) && tau >= 0.5 || throw(ArgumentError(
        "pilot autocorrelation time must be finite and >= 0.5."))

    # Sandvik, arXiv:1101.3281, Eqs. (46)--(48): bins must be much longer
    # than tau. A factor of ten is conservative while retaining enough bins.
    target_block_length = max(1, ceil(Int, 10 * tau))
    available_blocks = fld(N, target_block_length)
    if available_blocks < 20
        @warn(
            "SSE measurement run is too short for 20 autocorrelation-sized jackknife blocks.",
            label=label,
            N=N,
            tau=tau,
            target_block_length=target_block_length,
            available_blocks=available_blocks,
        )
    end
    return clamp(available_blocks, 2, min(N, 1_000))
end

function validate_observable_inputs(
    n::Int,
    beta_phys::Float64,
    pairs::Vector{Tuple{Int,Int}},
)
    n >= 1 || throw(ArgumentError("model size must be positive."))
    isfinite(beta_phys) && beta_phys > 0 || throw(ArgumentError(
        "beta_phys must be finite and > 0."))
    for (i, j) in pairs
        1 <= i <= n || throw(ArgumentError("pair index $i lies outside 1:$n."))
        1 <= j <= n || throw(ArgumentError("pair index $j lies outside 1:$n."))
        i != j || throw(ArgumentError("pair sites must be distinct (got ($i, $j))."))
    end
    return nothing
end

function validate_run_inputs(
    n::Int,
    beta_phys::Float64,
    pairs::Vector{Tuple{Int,Int}},
    nsweeps::Int,
    nwarm::Int,
)
    validate_observable_inputs(n, beta_phys, pairs)
    nsweeps >= 2 || throw(ArgumentError("nsweeps must be at least 2."))
    nwarm >= 0 || throw(ArgumentError("nwarm must be nonnegative."))
    return nothing
end


# SSE state and diagonal updates.
mutable struct SSEState
    n::Int
    spins::Vector{Int}        # |alpha(0)>, ±1
    opstring::Vector{Op}
    M::Int
    n_op::Int
end

"""
    SSEResult

Thermal observables, jackknife errors, autocorrelation estimates, and update
diagnostics returned by [`run_sse`](@ref).
"""
struct SSEResult
    energy::Float64
    energy_err::Float64
    mz2::Float64
    mz2_err::Float64
    zz::Dict{Tuple{Int,Int},Tuple{Float64,Float64}}
    avg_sign::Float64
    avg_sign_err::Float64
    tau_E::Tuple{Float64,Float64}
    tau_mz::Tuple{Float64,Float64}          # τ_int of |m_z| (within-sector fluctuations)
    tau_mz_signed::Tuple{Float64,Float64}   # τ_int of signed m_z — the order-parameter /
                                            # sector-tunneling slow mode (ordered phase)
    mz_mean::Float64                        # ⟨signed m_z⟩ — freezing/ergodicity diagnostic:
                                            # ≈0 if both sectors visited, ≠0 if trapped in one
    n_signflips::Int                        # # sign changes of m_z over the measured sweeps
                                            # (tunnelling-event count; →0 when frozen)
    loop_accept::Float64
    mean_r::Float64
end

"Count sign changes of a (signed) series — the number of sector-tunnelling events."
function count_signflips(arr::Vector{Float64})
    c = 0
    previous_sign = 0
    @inbounds for value in arr
        current_sign = value > 0 ? 1 : value < 0 ? -1 : 0
        current_sign == 0 && continue
        previous_sign != 0 && current_sign != previous_sign && (c += 1)
        previous_sign = current_sign
    end
    return c
end
