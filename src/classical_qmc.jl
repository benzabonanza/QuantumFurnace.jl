"""
    ClassicalQMC

Sign-free Stochastic Series Expansion (SSE) quantum Monte Carlo thermal sampler for the
two QuantumFurnace thesis fixtures (1D disordered AFM Heisenberg, 2D disordered TFIM). This
is the *classical* Gibbs-sampling baseline that competes with the KMS Lindbladian quantum
sampler (beads qf-h23; design note `drafts/classical-gibbs-sampling/method.md`).

Algorithm (Sandvik, arXiv:1909.10591): diagonal update Eq. 19; deterministic operator-loop
[Sandvik 1999] for the Heisenberg model; Swendsen–Wang cluster [Sandvik 2003] for the TFIM.

Strategy:
- Sample in the **physical frame** (bare `J=1`, physical `h`, `β_phys`); `ρ_β` is
  frame-independent, so observables cross-check against the fixture's exact dense `ρ_β`.
- A model is reconstructed in physical units by replaying the `build_heis_1d` /
  `build_tfim_2d` RNG recipe; [`sse_reconstruction_error`](@ref) checks the dense `H_phys`
  matches the fixture's un-rescaled Hamiltonian to ~1e-10 (exercised in the test suite — it
  certifies the *exact reference* targets the right model; the SSE weights are certified
  separately by the MC-vs-exact agreement).
- Energy estimator `⟨H⟩ = C_shift − ⟨n⟩/β` with `C_shift` the sum of the SSE operator
  constants (worked out per model below).

This module is self-contained (it rebuilds the fixture Hamiltonians from scratch rather than
importing the dense builders), wrapped as a sub-module so its many generic internal names
(`Op`, `EXCH`, `FIELD`, …) stay isolated from the main `QuantumFurnace` namespace.

Public API: [`build_sse_heis_model`](@ref), [`build_sse_tfim_model`](@ref), [`run_sse`](@ref),
[`sse_exact_reference`](@ref), [`sse_reconstruction_error`](@ref), [`SSEResult`](@ref).

Honest scope: the frustrated odd-`n` *periodic* Heisenberg sign sector is NOT solved — the
deterministic operator-loop is non-ergodic there (true `⟨s⟩<1` not captured); use OBC
(always bipartite ⇒ sign-free) or a worm/directed-loop update (beads qf-h23.3.8). Even-`n`
PBC, all OBC, and the TFIM are fully sign-free and validated.
"""
module ClassicalQMC

using LinearAlgebra
using Random
using Statistics

export SSEResult, run_sse,
    build_sse_heis_model, build_sse_tfim_model,
    sse_exact_reference, sse_reconstruction_error

# ----------------------------------------------------------------------------------------------
# Pauli matrices and dense operator builders (computational basis; site 1 = leftmost kron factor,
# matching pad_term / _pad_two_site_op in src/hamiltonian.jl). Spin convention: +1 = up (Z=+1),
# -1 = down (Z=-1).
# ----------------------------------------------------------------------------------------------
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

# ----------------------------------------------------------------------------------------------
# Model representation (physical units). Operator constants follow the worked-out algebra in the
# plan. Spins are ±1 integers.
# ----------------------------------------------------------------------------------------------

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

struct HeisModel
    n::Int
    bonds::Vector{Tuple{Int,Int}}     # exchange bonds (with duplicates for wrap / n=2)
    J::Float64
    field::Vector{Float64}            # hz per site (length n)
    zz::Vector{Float64}               # Kzz per bond (aligned with bonds)
    frustrated::Vector{Bool}          # per bond: intra-sublattice (Marshall sign -1)
    C_shift::Float64
    H_phys::Matrix{ComplexF64}
    periodic::Bool
end

struct TfimModel
    n::Int
    Lx::Int
    Ly::Int
    bonds::Vector{Tuple{Int,Int}}     # Ising bonds (with duplicates for wrap)
    J::Float64
    h::Float64
    zdis::Vector{Float64}             # eps_z per site
    zz::Vector{Float64}               # eps_zz per bond (aligned with bonds)
    C_shift::Float64
    H_phys::Matrix{ComplexF64}
    periodic::Bool
end

"2-colour the bond graph by BFS; a bond is frustrated iff its endpoints share a colour."
function frustrated_bonds(n::Int, bonds::Vector{Tuple{Int,Int}})
    color = fill(0, n)
    adj = [Tuple{Int,Int}[] for _ in 1:n]      # (neighbour, bondsign placeholder)
    for (i, j) in bonds
        push!(adj[i], (j, 0))
        push!(adj[j], (i, 0))
    end
    for s0 in 1:n
        color[s0] != 0 && continue
        color[s0] = 1
        queue = [s0]
        while !isempty(queue)
            u = popfirst!(queue)
            for (v, _) in adj[u]
                if color[v] == 0
                    color[v] = 3 - color[u]
                    push!(queue, v)
                end
            end
        end
    end
    return [color[i] == color[j] for (i, j) in bonds]
end

"""
    build_sse_heis_model(n; seed, periodic, disorder_strength, J=1.0) -> HeisModel

Reconstruct the 1D Heisenberg fixture in physical units by replaying the `build_heis_1d`
RNG recipe (MersenneTwister(seed); [Z] per-site draw THEN [Z,Z] per-bond draw, each `rand!`
then `.*= disorder_strength`). Builds the dense `H_phys` for the reconstruction cross-check.
"""
function build_sse_heis_model(n::Int; seed::Int, periodic::Bool, disorder_strength::Float64, J::Float64=1.0)
    # Exchange bonds: pad_term places [O,O] at (q,q+1), wrap (n,1) when periodic.
    bonds = Tuple{Int,Int}[]
    for q in 1:(periodic ? n : n - 1)
        push!(bonds, q < n ? (q, q + 1) : (n, 1))
    end

    # Disorder draw (order: field hz first, then bond Kzz) — matches build_heis_1d exactly.
    rng = MersenneTwister(seed)
    hz = zeros(Float64, n); rand!(rng, hz); hz .*= disorder_strength
    Kzz_site = zeros(Float64, n); rand!(rng, Kzz_site); Kzz_site .*= disorder_strength
    # zz coeff for bond enumerated at position q is Kzz_site[q]
    zz = Float64[Kzz_site[q < n ? q : n] for q in 1:(periodic ? n : n - 1)]

    # Dense H_phys = J sum (XX+YY+ZZ) + sum hz Z + sum Kzz ZZ
    d = 2^n
    H = zeros(ComplexF64, d, d)
    for (i, j) in bonds
        H .+= J .* (two_site(PX, i, j, n) .+ two_site(PY, i, j, n) .+ two_site(PZ, i, j, n))
    end
    for q in 1:n
        H .+= hz[q] .* one_site(PZ, q, n)
    end
    for (b, (i, j)) in enumerate(bonds)
        H .+= zz[b] .* two_site(PZ, i, j, n)
    end

    frust = frustrated_bonds(n, bonds)
    C_shift = J * length(bonds) + sum(abs, hz) + sum(abs, zz)
    return HeisModel(n, bonds, J, hz, zz, frust, C_shift, H, periodic)
end

"""
    build_sse_tfim_model(Lx, Ly; seed, h, disorder_strength, J=1.0, periodic=true) -> TfimModel

Reconstruct the 2D TFIM fixture in physical units (build_tfim_2d recipe). Bonds: right + up
neighbour with periodic wrap (doubled on L=2); field -h X per site; disorder [Z] then [Z,Z].
"""
function build_sse_tfim_model(Lx::Int, Ly::Int; seed::Int, h::Float64, disorder_strength::Float64,
                          J::Float64=1.0, periodic::Bool=true)
    n = Lx * Ly
    site_index(i, j) = (i - 1) * Ly + (j - 1) + 1

    rng = MersenneTwister(seed)
    ez = zeros(Float64, n); rand!(rng, ez); ez .*= disorder_strength
    ezz_site = zeros(Float64, n); rand!(rng, ezz_site); ezz_site .*= disorder_strength

    # Build bonds and aligned per-bond ezz in the SAME enumeration (right + up neighbour).
    bonds = Tuple{Int,Int}[]
    zz = Float64[]
    for i in 1:Lx, j in 1:Ly
        c = ezz_site[site_index(i, j)]
        if i < Lx
            push!(bonds, (site_index(i, j), site_index(i + 1, j))); push!(zz, c)
        elseif periodic && Lx > 1
            push!(bonds, (site_index(Lx, j), site_index(1, j))); push!(zz, c)
        end
        if j < Ly
            push!(bonds, (site_index(i, j), site_index(i, j + 1))); push!(zz, c)
        elseif periodic && Ly > 1
            push!(bonds, (site_index(i, Ly), site_index(i, 1))); push!(zz, c)
        end
    end

    d = 2^n
    H = zeros(ComplexF64, d, d)
    for (i, j) in bonds
        H .+= (-J) .* two_site(PZ, i, j, n)
    end
    for q in 1:n
        H .+= (-h) .* one_site(PX, q, n)
    end
    for q in 1:n
        H .+= ez[q] .* one_site(PZ, q, n)
    end
    for (b, (i, j)) in enumerate(bonds)
        H .+= zz[b] .* two_site(PZ, i, j, n)
    end

    C_shift = J * length(bonds) + h * n + sum(abs, ez) + sum(abs, zz)
    return TfimModel(n, Lx, Ly, bonds, J, h, ez, zz, C_shift, H, periodic)
end

# ----------------------------------------------------------------------------------------------
# Reconstruction cross-check vs fixture: H_phys ?= (raw.matrix - raw.shift*I)*raw.rescaling_factor
# ----------------------------------------------------------------------------------------------
function fixture_H_phys(raw)
    d = size(raw.matrix, 1)
    return (raw.matrix .- raw.shift .* Matrix{ComplexF64}(I, d, d)) .* raw.rescaling_factor
end

"""
    sse_reconstruction_error(model, raw) -> Float64

Maximum absolute entrywise difference between the SSE model's reconstructed physical
Hamiltonian `model.H_phys` and the fixture's un-rescaled Hamiltonian
`(raw.matrix - raw.shift*I) * raw.rescaling_factor`, where `raw` is the NamedTuple returned
by `build_heis_1d` / `build_tfim_2d`. A value `≤ 1e-10` certifies that the
exact reference targets exactly the fixture model. Not called by `run_sse` (the model is
self-contained and carries no fixture handle); it is exercised in `test/test_classical_qmc.jl`.
"""
sse_reconstruction_error(model::HeisModel, raw) = maximum(abs.(model.H_phys .- fixture_H_phys(raw)))
sse_reconstruction_error(model::TfimModel, raw) = maximum(abs.(model.H_phys .- fixture_H_phys(raw)))

# ----------------------------------------------------------------------------------------------
# Exact thermal reference from dense H_phys (physical frame).
# ----------------------------------------------------------------------------------------------
struct ExactRef
    energy::Float64
    mz2::Float64
    zz::Dict{Tuple{Int,Int},Float64}     # <Z_i Z_j>
end

"""
    sse_exact_reference(model, beta_phys; pairs) -> ExactRef

Exact dense thermal reference (`⟨H⟩`, `⟨m_z²⟩`, `⟨Z_iZ_j⟩` for the requested `pairs`) from
`ρ_β = e^{-β_phys H_phys}/Z` of the model's reconstructed physical Hamiltonian. The
ground-truth the QMC estimates are cross-checked against at small `n`.
"""
sse_exact_reference(model::HeisModel, beta_phys::Float64; pairs::Vector{Tuple{Int,Int}}) =
    _exact_reference(model.H_phys, model.n, beta_phys; pairs=pairs)
sse_exact_reference(model::TfimModel, beta_phys::Float64; pairs::Vector{Tuple{Int,Int}}) =
    _exact_reference(model.H_phys, model.n, beta_phys; pairs=pairs)

function _exact_reference(H::Matrix{ComplexF64}, n::Int, beta_phys::Float64; pairs::Vector{Tuple{Int,Int}})
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

# ----------------------------------------------------------------------------------------------
# Statistics: block jackknife (with optional sign reweighting) + Sokal integrated autocorr time.
# ----------------------------------------------------------------------------------------------
"Block-jackknife mean and error of the ratio mean(X)/mean(Y). For Y≡1, plain jackknife of mean(X)."
function jackknife_ratio(X::Vector{Float64}, Y::Vector{Float64}; nblocks::Int=40)
    N = length(X)
    nb = min(nblocks, N)
    edges = round.(Int, range(0, N; length=nb + 1))
    sumX = sum(X); sumY = sum(Y)
    Rb = Float64[]
    for b in 1:nb
        lo, hi = edges[b] + 1, edges[b + 1]
        sx = sum(@view X[lo:hi]); sy = sum(@view Y[lo:hi])
        push!(Rb, (sumX - sx) / (sumY - sy))
    end
    Rjk = mean(Rb)
    err = sqrt((nb - 1) / nb * sum((Rb .- Rjk) .^ 2))
    Rfull = sumX / sumY
    return Rfull, err
end

"Integrated autocorrelation time via Sokal automatic windowing. Returns (tau, err, window)."
function tau_int(series::Vector{Float64}; c::Float64=6.0)
    N = length(series)
    m = mean(series)
    v = sum((series .- m) .^ 2) / N
    if v ≤ 0
        return 0.5, 0.0, 0
    end
    tau = 0.5
    W = N - 1
    for t in 1:(N - 1)
        rho_t = 0.0
        @inbounds for k in 1:(N - t)
            rho_t += (series[k] - m) * (series[k + t] - m)
        end
        rho_t /= (N - t) * v
        tau += rho_t
        if t ≥ c * tau
            W = t
            break
        end
    end
    err = tau * sqrt(2 * (2 * W + 1) / N)     # Madras–Sokal
    return tau, err, W
end

# ----------------------------------------------------------------------------------------------
# SSE configuration + generic diagonal update.
# ----------------------------------------------------------------------------------------------
mutable struct SSEState
    n::Int
    spins::Vector{Int}        # |alpha(0)>, ±1
    opstring::Vector{Op}
    M::Int
    n_op::Int
end

"Diagonal matrix element <alpha|H_op|alpha> on the current local spins."
@inline function diag_element(kind::Int8, idx::Int32, spins::Vector{Int}, m::HeisModel)
    if kind == EXCH
        (i, j) = m.bonds[idx]
        return spins[i] != spins[j] ? 2.0 * m.J : 0.0
    elseif kind == FIELD
        hz = m.field[idx]
        return abs(hz) - hz * spins[idx]
    elseif kind == ZZDIS
        (i, j) = m.bonds[idx]
        K = m.zz[idx]
        return abs(K) - K * spins[i] * spins[j]
    end
    return 0.0
end

"Off-diagonal propagation: apply the operator's spin flip to `state`."
@inline function apply_offdiag!(state::Vector{Int}, op::Op, m::HeisModel)
    if op.kind == EXCH
        (i, j) = m.bonds[op.idx]
        state[i] = -state[i]; state[j] = -state[j]
    end
    return nothing
end

"Build the diagonal-operator menu (nonzero-coupling instances). Returns Vector{(kind,idx)}."
function heis_menu(m::HeisModel)
    menu = Tuple{Int8,Int32}[]
    for b in eachindex(m.bonds)
        push!(menu, (EXCH, Int32(b)))
    end
    for i in 1:m.n
        m.field[i] != 0 && push!(menu, (FIELD, Int32(i)))
    end
    for b in eachindex(m.bonds)
        m.zz[b] != 0 && push!(menu, (ZZDIS, Int32(b)))
    end
    return menu
end

"One diagonal update sweep (Sandvik Eq.19). Mutates st.opstring and st.n_op; propagates spins."
function diagonal_update!(st::SSEState, m::HeisModel, menu::Vector{Tuple{Int8,Int32}},
                          beta::Float64, rng::AbstractRNG)
    Nd = length(menu)
    state = st.spins                      # working copy = |alpha(0)>; mutated then restored
    work = copy(state)
    @inbounds for p in 1:st.M
        op = st.opstring[p]
        if is_id(op)
            (kind, idx) = menu[rand(rng, 1:Nd)]
            elem = diag_element(kind, idx, work, m)
            if elem > 0
                if rand(rng) < beta * Nd * elem / (st.M - st.n_op)
                    st.opstring[p] = Op(kind, idx, false)
                    st.n_op += 1
                end
            end
        elseif is_diag(op)
            elem = diag_element(op.kind, op.idx, work, m)
            if rand(rng) < (st.M - st.n_op + 1) / (beta * Nd * elem)
                st.opstring[p] = IDOP
                st.n_op -= 1
            end
        else
            apply_offdiag!(work, op, m)
        end
    end
    return nothing
end

# ----------------------------------------------------------------------------------------------
# Heisenberg deterministic operator-loop (Sandvik 1991) with sequential heat-bath disorder
# correction (full-recompute of the disorder weight per loop — exact, cheap at tiny n).
# Leg layout per exchange vertex v (legs 4(v-1)+1..+4): 1=bot-i, 2=bot-j, 3=top-i, 4=top-j.
# Vertex pairing for the SU(2) deterministic loop: 1<->2, 3<->4 (flip a bottom/top pair toggles
# the operator diag<->offdiag, keeping the antiparallel-only vertex valid).
# ----------------------------------------------------------------------------------------------
@inline function vertex_partner(leg::Int)
    l = ((leg - 1) % 4) + 1
    base = leg - l
    pl = l == 1 ? 2 : l == 2 ? 1 : l == 3 ? 4 : 3
    return base + pl
end

"Disorder log-weight: propagate alpha0 and sum log(element) over FIELD + ZZDIS operators."
function disorder_logweight(spins::Vector{Int}, opstring::Vector{Op}, m::HeisModel)
    work = copy(spins)
    lw = 0.0
    @inbounds for p in eachindex(opstring)
        op = opstring[p]
        if op.kind == FIELD
            hz = m.field[op.idx]
            e = abs(hz) - hz * work[op.idx]
            lw += (e > 0 ? log(e) : -Inf)
        elseif op.kind == ZZDIS
            (i, j) = m.bonds[op.idx]
            K = m.zz[op.idx]
            e = abs(K) - K * work[i] * work[j]
            lw += (e > 0 ? log(e) : -Inf)
        elseif is_offdiag(op)
            apply_offdiag!(work, op, m)
        end
    end
    return lw
end

"One operator-loop update sweep. Returns (n_loops, n_flipped, sum_r) for diagnostics."
function loop_update!(st::SSEState, m::HeisModel, rng::AbstractRNG)
    n = st.n
    opstring = st.opstring
    has_disorder = any(!=(0), m.field) || any(!=(0), m.zz)

    # Collect exchange vertices and their leg spins (propagate to record bottom/top spins).
    exch_pos = Int[]
    for p in 1:st.M
        opstring[p].kind == EXCH && push!(exch_pos, p)
    end
    V = length(exch_pos)

    # Free-site bookkeeping handled after loops; if no exchange ops, just flip free spins.
    if V == 0
        free_site_flips!(st, m, has_disorder, rng)
        return (0, 0, 0.0)
    end

    legspin = Vector{Int}(undef, 4V)
    link = fill(-1, 4V)
    first_leg = fill(-1, n)
    last_leg = fill(-1, n)

    work = copy(st.spins)
    vv = 0
    @inbounds for p in 1:st.M
        op = opstring[p]
        if op.kind == EXCH
            vv += 1
            base = 4 * (vv - 1)
            (i, j) = m.bonds[op.idx]
            legspin[base + 1] = work[i]; legspin[base + 2] = work[j]
            if op.offdiag
                work[i] = -work[i]; work[j] = -work[j]
            end
            legspin[base + 3] = work[i]; legspin[base + 4] = work[j]
        elseif is_offdiag(op)
            apply_offdiag!(work, op, m)
        end
    end

    # Build worldline linked list (legs of consecutive exchange ops on each site, periodic in tau).
    vv = 0
    @inbounds for p in 1:st.M
        opstring[p].kind != EXCH && continue
        vv += 1
        base = 4 * (vv - 1)
        (i, j) = m.bonds[opstring[p].idx]
        for (sbot, stop, site) in ((base + 1, base + 3, i), (base + 2, base + 4, j))
            if last_leg[site] == -1
                first_leg[site] = sbot
            else
                link[last_leg[site]] = sbot
                link[sbot] = last_leg[site]
            end
            last_leg[site] = stop
        end
    end
    for s in 1:n
        if last_leg[s] != -1
            link[last_leg[s]] = first_leg[s]
            link[first_leg[s]] = last_leg[s]
        end
    end

    # Find loops (alternating vertex-partner and worldline-link moves).
    visited = falses(4V)
    loops = Vector{Vector{Int}}()
    for g in 1:4V
        visited[g] && continue
        legs = Int[]
        leg = g
        while true
            visited[leg] = true; push!(legs, leg)
            pp = vertex_partner(leg)
            visited[pp] = true; push!(legs, pp)
            leg = link[pp]
            leg == g && break
        end
        push!(loops, legs)
    end

    # Reconstruct (alpha0, exchange op types) from current legspins.
    function reconstruct!()
        @inbounds for s in 1:n
            first_leg[s] != -1 && (st.spins[s] = legspin[first_leg[s]])
        end
        vc = 0
        @inbounds for p in 1:st.M
            opstring[p].kind != EXCH && continue
            vc += 1
            base = 4 * (vc - 1)
            offd = legspin[base + 1] != legspin[base + 3]
            opstring[p] = Op(EXCH, opstring[p].idx, offd)
        end
    end

    n_flipped = 0
    sum_r = 0.0
    W0 = has_disorder ? disorder_logweight(st.spins, opstring, m) : 0.0
    for legs in loops
        if !has_disorder
            if rand(rng) < 0.5
                @inbounds for g in legs
                    legspin[g] = -legspin[g]
                end
                n_flipped += 1
            end
            sum_r += 1.0
            continue
        end
        # tentative flip
        @inbounds for g in legs
            legspin[g] = -legspin[g]
        end
        reconstruct!()
        W1 = disorder_logweight(st.spins, opstring, m)
        r = exp(W1 - W0)
        sum_r += r
        if rand(rng) < r / (1 + r)
            W0 = W1
            n_flipped += 1
        else
            @inbounds for g in legs
                legspin[g] = -legspin[g]
            end
            reconstruct!()
        end
    end
    reconstruct!()

    free_site_flips!(st, m, has_disorder, rng)
    return (length(loops), n_flipped, sum_r)
end

"Flip free-site boundary spins (sites with no exchange operators) by heat-bath on the field/zz weight."
function free_site_flips!(st::SSEState, m::HeisModel, has_disorder::Bool, rng::AbstractRNG)
    has_exch = falses(st.n)
    for p in 1:st.M
        op = st.opstring[p]
        if op.kind == EXCH
            (i, j) = m.bonds[op.idx]
            has_exch[i] = true; has_exch[j] = true
        end
    end
    for s in 1:st.n
        has_exch[s] && continue
        if !has_disorder
            rand(rng) < 0.5 && (st.spins[s] = -st.spins[s])
        else
            W0 = disorder_logweight(st.spins, st.opstring, m)
            st.spins[s] = -st.spins[s]
            W1 = disorder_logweight(st.spins, st.opstring, m)
            r = exp(W1 - W0)
            if !(rand(rng) < r / (1 + r))
                st.spins[s] = -st.spins[s]   # revert
            end
        end
    end
    return nothing
end

# ----------------------------------------------------------------------------------------------
# Driver + measurement (Heisenberg).
# ----------------------------------------------------------------------------------------------
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
    @inbounds for i in 2:length(arr)
        (arr[i] * arr[i-1] < 0) && (c += 1)
    end
    return c
end

"Configuration sign s = (-1)^(# off-diagonal exchange ops on frustrated bonds)."
function config_sign(st::SSEState, m::HeisModel)
    cnt = 0
    for p in 1:st.M
        op = st.opstring[p]
        if op.kind == EXCH && op.offdiag && m.frustrated[op.idx]
            cnt += 1
        end
    end
    return isodd(cnt) ? -1.0 : 1.0
end

"""
    run_sse(model, beta_phys; pairs, nsweeps=200_000, nwarm=40_000, seed=12345) -> SSEResult

Run the SSE sampler on a Heisenberg or TFIM `model` (from [`build_sse_heis_model`] /
[`build_sse_tfim_model`]) at physical inverse temperature `beta_phys`, returning an immutable
[`SSEResult`] with thermal observables (`⟨H⟩`, `⟨m_z²⟩`, `⟨Z_iZ_j⟩` for `pairs`), the
average sign, integrated autocorrelation times for the energy, `|m_z|`, and the signed `m_z`
(the order-parameter / sector-tunneling slow mode in an ordered phase), and update
diagnostics. `nwarm` equilibration sweeps (with `M` adaptation) precede `nsweeps` measured
sweeps; `seed` seeds the QMC RNG (independent of the fixture disorder seed). β is `β_phys`.

Heisenberg uses the deterministic operator-loop; TFIM uses the Swendsen–Wang cluster by
default (`update=:cluster`), or a local single-segment Metropolis baseline (`update=:local`,
TFIM only) that exhibits critical slowing down in the ordered phase — the classical-vs-quantum
mixing-time contrast. For the frustrated odd-`n` periodic Heisenberg the loop is non-ergodic in
the sign sector (`avg_sign` stays 1; true `⟨s⟩<1` not captured) — see the module docstring.
"""
function run_sse(m::HeisModel, beta_phys::Float64; pairs::Vector{Tuple{Int,Int}},
                        nsweeps::Int=200_000, nwarm::Int=40_000, seed::Int=12345)
    rng = MersenneTwister(seed)
    menu = heis_menu(m)
    st = SSEState(m.n, rand(rng, (-1, 1), m.n), Op[IDOP for _ in 1:20], 20, 0)

    # Warmup with M adaptation.
    nmax = 1
    for sweep in 1:nwarm
        diagonal_update!(st, m, menu, beta_phys, rng)
        loop_update!(st, m, rng)
        nmax = max(nmax, st.n_op)
        if sweep % 200 == 0
            Mnew = max(st.M, ceil(Int, 1.3 * nmax))
            if Mnew > st.M
                append!(st.opstring, [IDOP for _ in 1:(Mnew - st.M)])
                st.M = Mnew
            end
        end
    end

    # Measurement.
    n_arr = Float64[]; mz2_arr = Float64[]; absmz_arr = Float64[]; mz_arr = Float64[]
    sign_arr = Float64[]
    zz_arr = Dict(pr => Float64[] for pr in pairs)
    loops_total = 0; flips_total = 0; r_total = 0.0
    for _ in 1:nsweeps
        diagonal_update!(st, m, menu, beta_phys, rng)
        (nl, nf, sr) = loop_update!(st, m, rng)
        loops_total += nl; flips_total += nf; r_total += sr
        s = config_sign(st, m)
        push!(sign_arr, s)
        push!(n_arr, st.n_op)
        mz = sum(st.spins) / m.n
        push!(mz2_arr, mz^2)
        push!(absmz_arr, abs(mz))
        push!(mz_arr, mz)
        for pr in pairs
            push!(zz_arr[pr], st.spins[pr[1]] * st.spins[pr[2]])
        end
    end

    # Sign-reweighted estimators: <O> = <O*s>/<s>.
    avg_sign, avg_sign_err = jackknife_ratio(sign_arr, ones(length(sign_arr)))
    nbar, nbar_err = jackknife_ratio(n_arr .* sign_arr, sign_arr)
    energy = m.C_shift - nbar / beta_phys
    energy_err = nbar_err / beta_phys
    mz2, mz2_err = jackknife_ratio(mz2_arr .* sign_arr, sign_arr)
    zzres = Dict{Tuple{Int,Int},Tuple{Float64,Float64}}()
    for pr in pairs
        zzres[pr] = jackknife_ratio(zz_arr[pr] .* sign_arr, sign_arr)
    end

    tauE = tau_int(n_arr); taumz = tau_int(absmz_arr); taumzs = tau_int(mz_arr)
    loop_accept = loops_total > 0 ? flips_total / loops_total : 0.0
    mean_r = loops_total > 0 ? r_total / loops_total : 1.0
    return SSEResult(energy, energy_err, mz2, mz2_err, zzres, avg_sign, avg_sign_err,
                     (tauE[1], tauE[2]), (taumz[1], taumz[2]), (taumzs[1], taumzs[2]),
                     sum(mz_arr) / length(mz_arr), count_signflips(mz_arr), loop_accept, mean_r)
end

# ----------------------------------------------------------------------------------------------
# TFIM: SSE with Swendsen–Wang cluster update (Sandvik 2003). Operators: ISING bond (diag,
# element 2J on parallel), FIELDC (diag site, element h), FIELDX (off-diag site, element h,
# flips the spin), ZDIS / ZZDIS_T diagonal disorder ride-alongs. Stoquastic ⇒ sign-free (s≡1).
# Cluster: segments (cut by field ops) unioned across Ising bonds; flip each cluster (heat-bath
# on the tiny Z/ZZ disorder), toggling FIELDC↔FIELDX at segment boundaries (weight-neutral).
# ----------------------------------------------------------------------------------------------

@inline function tfim_diag_element(kind::Int8, idx::Int32, spins::Vector{Int}, m::TfimModel)
    if kind == ISING
        (i, j) = m.bonds[idx]
        return spins[i] == spins[j] ? 2.0 * m.J : 0.0
    elseif kind == FIELDC
        return m.h
    elseif kind == ZDIS
        ez = m.zdis[idx]
        return abs(ez) - ez * spins[idx]
    elseif kind == ZZDIS_T
        (i, j) = m.bonds[idx]
        K = m.zz[idx]
        return abs(K) - K * spins[i] * spins[j]
    end
    return 0.0
end

@inline function apply_offdiag_tfim!(state::Vector{Int}, op::Op)
    op.kind == FIELDX && (state[op.idx] = -state[op.idx])
    return nothing
end

function tfim_menu(m::TfimModel)
    menu = Tuple{Int8,Int32}[]
    for b in eachindex(m.bonds)
        push!(menu, (ISING, Int32(b)))
    end
    for i in 1:m.n
        push!(menu, (FIELDC, Int32(i)))
    end
    for i in 1:m.n
        m.zdis[i] != 0 && push!(menu, (ZDIS, Int32(i)))
    end
    for b in eachindex(m.bonds)
        m.zz[b] != 0 && push!(menu, (ZZDIS_T, Int32(b)))
    end
    return menu
end

function tfim_diagonal_update!(st::SSEState, m::TfimModel, menu::Vector{Tuple{Int8,Int32}},
                               beta::Float64, rng::AbstractRNG)
    Nd = length(menu)
    work = copy(st.spins)
    @inbounds for p in 1:st.M
        op = st.opstring[p]
        if is_id(op)
            (kind, idx) = menu[rand(rng, 1:Nd)]
            elem = tfim_diag_element(kind, idx, work, m)
            if elem > 0 && rand(rng) < beta * Nd * elem / (st.M - st.n_op)
                st.opstring[p] = Op(kind, idx, false)
                st.n_op += 1
            end
        elseif is_diag(op)
            elem = tfim_diag_element(op.kind, op.idx, work, m)
            if rand(rng) < (st.M - st.n_op + 1) / (beta * Nd * elem)
                st.opstring[p] = IDOP
                st.n_op -= 1
            end
        else
            apply_offdiag_tfim!(work, op)   # FIELDX
        end
    end
    return nothing
end

function tfim_disorder_logweight(spins::Vector{Int}, opstring::Vector{Op}, m::TfimModel)
    work = copy(spins)
    lw = 0.0
    @inbounds for p in eachindex(opstring)
        op = opstring[p]
        if op.kind == ZDIS
            ez = m.zdis[op.idx]
            e = abs(ez) - ez * work[op.idx]
            lw += (e > 0 ? log(e) : -Inf)
        elseif op.kind == ZZDIS_T
            (i, j) = m.bonds[op.idx]
            K = m.zz[op.idx]
            e = abs(K) - K * work[i] * work[j]
            lw += (e > 0 ? log(e) : -Inf)
        elseif op.kind == FIELDX
            work[op.idx] = -work[op.idx]
        end
    end
    return lw
end

"""
    tfim_offdiag_logweight(spins, opstring, m) -> Float64

Log-weight of the terms that change under a worldline-segment flip: the Ising-bond operators
(`log(2J)` when the bonded spins are parallel, `-Inf` when a flip would force them antiparallel
— the constraint that *locks* segments together) plus the diagonal disorder. Field operators
(`FIELDC`/`FIELDX`, both element `h`) are weight-neutral under a flip and omitted, but `FIELDX`
still propagates the running spin. Used by the local single-segment Metropolis update.
"""
function tfim_offdiag_logweight(spins::Vector{Int}, opstring::Vector{Op}, m::TfimModel)
    work = copy(spins)
    lw = 0.0
    @inbounds for p in eachindex(opstring)
        op = opstring[p]
        if op.kind == ISING
            (i, j) = m.bonds[op.idx]
            lw += (work[i] == work[j] ? log(2 * m.J) : -Inf)
        elseif op.kind == ZDIS
            ez = m.zdis[op.idx]
            e = abs(ez) - ez * work[op.idx]
            lw += (e > 0 ? log(e) : -Inf)
        elseif op.kind == ZZDIS_T
            (i, j) = m.bonds[op.idx]
            K = m.zz[op.idx]
            e = abs(K) - K * work[i] * work[j]
            lw += (e > 0 ? log(e) : -Inf)
        elseif op.kind == FIELDX
            work[op.idx] = -work[op.idx]
        end
    end
    return lw
end

uf_find(parent::Vector{Int}, x::Int) = begin
    while parent[x] != x
        parent[x] = parent[parent[x]]
        x = parent[x]
    end
    x
end
function uf_union!(parent::Vector{Int}, a::Int, b::Int)
    ra = uf_find(parent, a); rb = uf_find(parent, b)
    ra != rb && (parent[ra] = rb)
    return nothing
end

"One Swendsen–Wang cluster update sweep. Returns (n_clusters, n_flipped)."
function tfim_cluster_update!(st::SSEState, m::TfimModel, rng::AbstractRNG)
    n = st.n; M = st.M; ops = st.opstring
    has_disorder = any(!=(0), m.zdis) || any(!=(0), m.zz)

    # Field-op positions per site (in increasing p), defining segment cuts.
    field_pos = [Int[] for _ in 1:n]
    @inbounds for p in 1:M
        op = ops[p]
        (op.kind == FIELDC || op.kind == FIELDX) && push!(field_pos[op.idx], p)
    end

    # Global segment ids. Site s has nsegs_site segments; wrap segment index = nsegs_site.
    seg_start = zeros(Int, n); nsegs_site = zeros(Int, n); nseg = 0
    for s in 1:n
        F = length(field_pos[s])
        ns = F == 0 ? 1 : F
        seg_start[s] = nseg; nsegs_site[s] = ns; nseg += ns
    end

    seg_index = function (s::Int, p::Int)
        F = length(field_pos[s])
        F == 0 && return 1
        cnt = 0
        @inbounds for q in field_pos[s]
            q < p ? (cnt += 1) : break
        end
        return cnt == 0 ? F : cnt
    end

    # Segment spins from alpha0 propagated through the field ops.
    seg_spin = Vector{Int}(undef, nseg)
    for s in 1:n
        F = length(field_pos[s])
        wrap = seg_start[s] + (F == 0 ? 1 : F)
        seg_spin[wrap] = st.spins[s]
        if F > 0
            cur = st.spins[s]
            @inbounds for k in 1:(F - 1)
                opk = ops[field_pos[s][k]]
                cur = opk.kind == FIELDX ? -cur : cur
                seg_spin[seg_start[s] + k] = cur
            end
        end
    end

    # Union segments connected by Ising bonds (both sites must flip together).
    parent = collect(1:nseg)
    @inbounds for p in 1:M
        op = ops[p]
        if op.kind == ISING
            (i, j) = m.bonds[op.idx]
            uf_union!(parent, seg_start[i] + seg_index(i, p), seg_start[j] + seg_index(j, p))
        end
    end

    # Cluster membership.
    clusters = Dict{Int,Vector{Int}}()
    for sgid in 1:nseg
        r = uf_find(parent, sgid)
        push!(get!(clusters, r, Int[]), sgid)
    end

    function reconstruct_tfim!()
        @inbounds for s in 1:n
            wrap = seg_start[s] + nsegs_site[s]
            st.spins[s] = seg_spin[wrap]
        end
        @inbounds for s in 1:n
            F = length(field_pos[s])
            F == 0 && continue
            for k in 1:F
                below = k == 1 ? seg_start[s] + F : seg_start[s] + (k - 1)
                above = seg_start[s] + k
                isx = seg_spin[below] != seg_spin[above]
                ops[field_pos[s][k]] = Op(isx ? FIELDX : FIELDC, Int32(s), isx)
            end
        end
    end

    n_clusters = length(clusters)
    n_flipped = 0
    W0 = has_disorder ? tfim_disorder_logweight(st.spins, ops, m) : 0.0
    for (_, segs) in clusters
        if !has_disorder
            if rand(rng) < 0.5
                @inbounds for g in segs
                    seg_spin[g] = -seg_spin[g]
                end
                n_flipped += 1
            end
            continue
        end
        @inbounds for g in segs
            seg_spin[g] = -seg_spin[g]
        end
        reconstruct_tfim!()
        W1 = tfim_disorder_logweight(st.spins, ops, m)
        r = exp(W1 - W0)
        if rand(rng) < r / (1 + r)
            W0 = W1; n_flipped += 1
        else
            @inbounds for g in segs
                seg_spin[g] = -seg_spin[g]
            end
            reconstruct_tfim!()
        end
    end
    reconstruct_tfim!()
    return (n_clusters, n_flipped)
end

"""
    tfim_local_update!(st, m, rng) -> (n_segments, n_flipped)

LOCAL off-diagonal update: the no-cluster baseline that exhibits **critical slowing down** in
the ordered phase. Uses the identical worldline-segment decomposition as the cluster update but
flips each segment *alone* (not grouped into Swendsen–Wang clusters), accepting via Metropolis
on the full off-diagonal weight [`tfim_offdiag_logweight`]. An Ising bond that a lone flip would
turn antiparallel vetoes the move (`-Inf`), so in the ordered phase segments are locked together
and single-segment flips are mostly rejected — the magnetization barrier the cluster jumps in
one move. Correct (Metropolis detailed balance) but slow; for the classical-vs-quantum
mixing-time contrast only (`run_sse(...; update=:local)`).
"""
function tfim_local_update!(st::SSEState, m::TfimModel, rng::AbstractRNG)
    n = st.n; M = st.M; ops = st.opstring

    # Segment machinery — identical to tfim_cluster_update!.
    field_pos = [Int[] for _ in 1:n]
    @inbounds for p in 1:M
        op = ops[p]
        (op.kind == FIELDC || op.kind == FIELDX) && push!(field_pos[op.idx], p)
    end
    seg_start = zeros(Int, n); nsegs_site = zeros(Int, n); nseg = 0
    for s in 1:n
        F = length(field_pos[s]); ns = F == 0 ? 1 : F
        seg_start[s] = nseg; nsegs_site[s] = ns; nseg += ns
    end
    seg_spin = Vector{Int}(undef, nseg)
    for s in 1:n
        F = length(field_pos[s])
        wrap = seg_start[s] + (F == 0 ? 1 : F)
        seg_spin[wrap] = st.spins[s]
        if F > 0
            cur = st.spins[s]
            @inbounds for k in 1:(F - 1)
                opk = ops[field_pos[s][k]]
                cur = opk.kind == FIELDX ? -cur : cur
                seg_spin[seg_start[s] + k] = cur
            end
        end
    end
    function reconstruct_tfim!()
        @inbounds for s in 1:n
            wrap = seg_start[s] + nsegs_site[s]
            st.spins[s] = seg_spin[wrap]
        end
        @inbounds for s in 1:n
            F = length(field_pos[s]); F == 0 && continue
            for k in 1:F
                below = k == 1 ? seg_start[s] + F : seg_start[s] + (k - 1)
                above = seg_start[s] + k
                isx = seg_spin[below] != seg_spin[above]
                ops[field_pos[s][k]] = Op(isx ? FIELDX : FIELDC, Int32(s), isx)
            end
        end
    end

    # Single-segment Metropolis (the only difference from the cluster update).
    n_flipped = 0
    W0 = tfim_offdiag_logweight(st.spins, ops, m)
    @inbounds for g in randperm(rng, nseg)
        seg_spin[g] = -seg_spin[g]
        reconstruct_tfim!()
        W1 = tfim_offdiag_logweight(st.spins, ops, m)
        if rand(rng) < exp(W1 - W0)
            W0 = W1; n_flipped += 1
        else
            seg_spin[g] = -seg_spin[g]
            reconstruct_tfim!()
        end
    end
    reconstruct_tfim!()
    return (nseg, n_flipped)
end

function run_sse(m::TfimModel, beta_phys::Float64; pairs::Vector{Tuple{Int,Int}},
                  nsweeps::Int=200_000, nwarm::Int=40_000, seed::Int=12345,
                  update::Symbol=:cluster)
    offdiag! = update === :local ? tfim_local_update! :
               update === :cluster ? tfim_cluster_update! :
               throw(ArgumentError("update must be :cluster or :local, got :$update"))
    rng = MersenneTwister(seed)
    menu = tfim_menu(m)
    st = SSEState(m.n, rand(rng, (-1, 1), m.n), Op[IDOP for _ in 1:20], 20, 0)
    nmax = 1
    for sweep in 1:nwarm
        tfim_diagonal_update!(st, m, menu, beta_phys, rng)
        offdiag!(st, m, rng)
        nmax = max(nmax, st.n_op)
        if sweep % 200 == 0
            Mnew = max(st.M, ceil(Int, 1.3 * nmax))
            Mnew > st.M && (append!(st.opstring, [IDOP for _ in 1:(Mnew - st.M)]); st.M = Mnew)
        end
    end
    n_arr = Float64[]; mz2_arr = Float64[]; absmz_arr = Float64[]; mz_arr = Float64[]
    zz_arr = Dict(pr => Float64[] for pr in pairs)
    cl_total = 0; fl_total = 0
    for _ in 1:nsweeps
        tfim_diagonal_update!(st, m, menu, beta_phys, rng)
        (nc, nf) = offdiag!(st, m, rng)
        cl_total += nc; fl_total += nf
        push!(n_arr, st.n_op)
        mz = sum(st.spins) / m.n
        push!(mz2_arr, mz^2); push!(absmz_arr, abs(mz)); push!(mz_arr, mz)
        for pr in pairs
            push!(zz_arr[pr], st.spins[pr[1]] * st.spins[pr[2]])
        end
    end
    ones_ = ones(length(n_arr))
    nbar, nbar_err = jackknife_ratio(n_arr, ones_)
    energy = m.C_shift - nbar / beta_phys
    energy_err = nbar_err / beta_phys
    mz2, mz2_err = jackknife_ratio(mz2_arr, ones_)
    zzres = Dict{Tuple{Int,Int},Tuple{Float64,Float64}}()
    for pr in pairs
        zzres[pr] = jackknife_ratio(zz_arr[pr], ones_)
    end
    tauE = tau_int(n_arr); taumz = tau_int(absmz_arr); taumzs = tau_int(mz_arr)
    accept = cl_total > 0 ? fl_total / cl_total : 0.0
    return SSEResult(energy, energy_err, mz2, mz2_err, zzres, 1.0, 0.0,
                     (tauE[1], tauE[2]), (taumz[1], taumz[2]), (taumzs[1], taumzs[2]),
                     sum(mz_arr) / length(mz_arr), count_signflips(mz_arr), accept, 1.0)
end

end # module ClassicalQMC
