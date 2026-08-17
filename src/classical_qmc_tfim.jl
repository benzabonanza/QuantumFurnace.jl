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
    periodic_x::Bool
    periodic_y::Bool
end

"""
    build_sse_tfim_model(Lx, Ly; seed, h, disorder_strength, J=1.0,
        periodic_x=true, periodic_y=true) -> TfimModel

Reconstruct the 2D TFIM SSE model in physical units.

# Keywords
- `seed`: Seed for the fixture-compatible disorder draws.
- `h`: Transverse-field strength.
- `disorder_strength`: Scale of the longitudinal and bond disorder.
- `J`: Ising coupling.
- `periodic_x`, `periodic_y`: Whether to include boundary-crossing bonds in
  each lattice direction.

# Returns
A `TfimModel` containing the SSE data and dense physical Hamiltonian.
"""
function build_sse_tfim_model(Lx::Int, Ly::Int; seed::Int, h::Float64, disorder_strength::Float64,
                          J::Float64=1.0, periodic_x::Bool=true, periodic_y::Bool=true)
    Lx >= 1 || throw(ArgumentError("Lx must be at least 1."))
    Ly >= 1 || throw(ArgumentError("Ly must be at least 1."))
    isfinite(J) && J > 0 || throw(ArgumentError(
        "J must be finite and > 0 for the ferromagnetic SSE decomposition."))
    isfinite(h) && h > 0 || throw(ArgumentError(
        "h must be finite and > 0 for the transverse-field SSE decomposition."))
    isfinite(disorder_strength) && disorder_strength >= 0 || throw(ArgumentError(
        "disorder_strength must be finite and >= 0."))

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
        elseif periodic_x && Lx > 1
            push!(bonds, (site_index(Lx, j), site_index(1, j))); push!(zz, c)
        end
        if j < Ly
            push!(bonds, (site_index(i, j), site_index(i, j + 1))); push!(zz, c)
        elseif periodic_y && Ly > 1
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
    return TfimModel(
        n, Lx, Ly, bonds, J, h, ez, zz, C_shift, H, periodic_x, periodic_y)
end

sse_reconstruction_error(model::TfimModel, raw) = maximum(abs.(model.H_phys .- fixture_H_phys(raw)))

sse_exact_reference(model::TfimModel, beta_phys::Float64; pairs::Vector{Tuple{Int,Int}}) =
    _exact_reference(model.H_phys, model.n, beta_phys; pairs=pairs)

function validate_model(m::TfimModel)
    m.n == m.Lx * m.Ly || throw(DimensionMismatch("TFIM n must equal Lx * Ly."))
    m.Lx >= 1 && m.Ly >= 1 || throw(ArgumentError("TFIM dimensions must be positive."))
    isfinite(m.J) && m.J > 0 || throw(ArgumentError(
        "TFIM SSE requires finite ferromagnetic J > 0."))
    isfinite(m.h) && m.h > 0 || throw(ArgumentError(
        "TFIM SSE requires finite transverse field h > 0."))
    length(m.zdis) == m.n || throw(DimensionMismatch("TFIM site-disorder length must equal n."))
    length(m.zz) == length(m.bonds) || throw(DimensionMismatch(
        "TFIM bond-disorder length must match bonds."))
    all(isfinite, m.zdis) && all(isfinite, m.zz) || throw(ArgumentError(
        "TFIM SSE couplings must be finite."))
    isfinite(m.C_shift) || throw(ArgumentError("TFIM SSE energy shift must be finite."))
    for (i, j) in m.bonds
        1 <= i <= m.n && 1 <= j <= m.n && i != j || throw(ArgumentError(
            "TFIM bonds must join distinct sites in 1:$(m.n)."))
    end
    return nothing
end

# TFIM cluster update. Field operators cut worldline segments; Ising vertices
# join segments into clusters, and diagonal disorder supplies heat-bath weights.

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

Return the log-weight affected by a TFIM worldline-segment flip.

Ising bonds contribute `log(2J)` for parallel spins and `-Inf` otherwise;
field operators are weight-neutral but still propagate spins.
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

# Build the deterministic TFIM worldline segmentation shared by the cluster
# and local updates. The returned arrays retain the original global segment-id
# convention: each site's final id is its segment crossing the imaginary-time
# boundary.
function _tfim_segment_layout(st::SSEState)
    n = st.n
    ops = st.opstring

    field_positions = [Int[] for _ in 1:n]
    @inbounds for p in 1:st.M
        op = ops[p]
        (op.kind == FIELDC || op.kind == FIELDX) &&
            push!(field_positions[op.idx], p)
    end

    segment_starts = zeros(Int, n)
    segment_counts = zeros(Int, n)
    n_segments = 0
    for site in 1:n
        n_field = length(field_positions[site])
        count = n_field == 0 ? 1 : n_field
        segment_starts[site] = n_segments
        segment_counts[site] = count
        n_segments += count
    end

    segment_spins = Vector{Int}(undef, n_segments)
    for site in 1:n
        n_field = length(field_positions[site])
        wrap = segment_starts[site] + (n_field == 0 ? 1 : n_field)
        segment_spins[wrap] = st.spins[site]
        if n_field > 0
            current = st.spins[site]
            @inbounds for k in 1:(n_field - 1)
                field_op = ops[field_positions[site][k]]
                current = field_op.kind == FIELDX ? -current : current
                segment_spins[segment_starts[site] + k] = current
            end
        end
    end

    return field_positions, segment_starts, segment_counts, segment_spins
end

@inline function _tfim_segment_index(field_positions::Vector{Int}, p::Int)
    n_field = length(field_positions)
    n_field == 0 && return 1
    count = 0
    @inbounds for q in field_positions
        q < p ? (count += 1) : break
    end
    return count == 0 ? n_field : count
end

function _reconstruct_tfim_worldlines!(
    st::SSEState,
    field_positions::Vector{Vector{Int}},
    segment_starts::Vector{Int},
    segment_counts::Vector{Int},
    segment_spins::Vector{Int},
)
    @inbounds for site in 1:st.n
        wrap = segment_starts[site] + segment_counts[site]
        st.spins[site] = segment_spins[wrap]
    end
    @inbounds for site in 1:st.n
        n_field = length(field_positions[site])
        n_field == 0 && continue
        for k in 1:n_field
            below = k == 1 ?
                segment_starts[site] + n_field : segment_starts[site] + (k - 1)
            above = segment_starts[site] + k
            is_x = segment_spins[below] != segment_spins[above]
            st.opstring[field_positions[site][k]] =
                Op(is_x ? FIELDX : FIELDC, Int32(site), is_x)
        end
    end
    return nothing
end

"One Swendsen–Wang cluster update sweep. Returns (n_clusters, n_flipped)."
function tfim_cluster_update!(st::SSEState, m::TfimModel, rng::AbstractRNG)
    n = st.n; M = st.M; ops = st.opstring
    has_disorder = any(!=(0), m.zdis) || any(!=(0), m.zz)

    field_pos, seg_start, nsegs_site, seg_spin = _tfim_segment_layout(st)
    nseg = length(seg_spin)

    # Union segments connected by Ising bonds (both sites must flip together).
    parent = collect(1:nseg)
    @inbounds for p in 1:M
        op = ops[p]
        if op.kind == ISING
            (i, j) = m.bonds[op.idx]
            uf_union!(
                parent,
                seg_start[i] + _tfim_segment_index(field_pos[i], p),
                seg_start[j] + _tfim_segment_index(field_pos[j], p),
            )
        end
    end

    # Cluster membership.
    clusters = Dict{Int,Vector{Int}}()
    for sgid in 1:nseg
        r = uf_find(parent, sgid)
        push!(get!(clusters, r, Int[]), sgid)
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
        _reconstruct_tfim_worldlines!(
            st, field_pos, seg_start, nsegs_site, seg_spin)
        W1 = tfim_disorder_logweight(st.spins, ops, m)
        r = exp(W1 - W0)
        if rand(rng) < r / (1 + r)
            W0 = W1; n_flipped += 1
        else
            @inbounds for g in segs
                seg_spin[g] = -seg_spin[g]
            end
            _reconstruct_tfim_worldlines!(
                st, field_pos, seg_start, nsegs_site, seg_spin)
        end
    end
    _reconstruct_tfim_worldlines!(
        st, field_pos, seg_start, nsegs_site, seg_spin)
    return (n_clusters, n_flipped)
end

"""
    tfim_local_update!(st, m, rng) -> (n_segments, n_flipped)

Apply one local single-segment Metropolis update to the TFIM worldlines.

This correct but slow baseline exposes critical slowing down; the cluster
update is the production default.
"""
function tfim_local_update!(st::SSEState, m::TfimModel, rng::AbstractRNG)
    ops = st.opstring
    field_pos, seg_start, nsegs_site, seg_spin = _tfim_segment_layout(st)
    nseg = length(seg_spin)

    # Single-segment Metropolis (the only difference from the cluster update).
    n_flipped = 0
    W0 = tfim_offdiag_logweight(st.spins, ops, m)
    @inbounds for g in randperm(rng, nseg)
        seg_spin[g] = -seg_spin[g]
        _reconstruct_tfim_worldlines!(
            st, field_pos, seg_start, nsegs_site, seg_spin)
        W1 = tfim_offdiag_logweight(st.spins, ops, m)
        if rand(rng) < exp(W1 - W0)
            W0 = W1; n_flipped += 1
        else
            seg_spin[g] = -seg_spin[g]
            _reconstruct_tfim_worldlines!(
                st, field_pos, seg_start, nsegs_site, seg_spin)
        end
    end
    _reconstruct_tfim_worldlines!(
        st, field_pos, seg_start, nsegs_site, seg_spin)
    return (nseg, n_flipped)
end

function run_sse(m::TfimModel, beta_phys::Float64; pairs::Vector{Tuple{Int,Int}},
                  nsweeps::Int=200_000, nwarm::Int=40_000, seed::Int=12345,
                  update::Symbol=:cluster)
    offdiag! = update === :local ? tfim_local_update! :
               update === :cluster ? tfim_cluster_update! :
               throw(ArgumentError("update must be :cluster or :local, got :$update"))
    validate_model(m)
    validate_run_inputs(m.n, beta_phys, pairs, nsweeps, nwarm)

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
    tauE = tau_int(n_arr; label="expansion order")
    taumz = tau_int(absmz_arr; label="absolute magnetisation")
    taumzs = tau_int(mz_arr; label="signed magnetisation")
    tau_mz2 = tau_int(mz2_arr; label="squared magnetisation")
    tau_zz = [tau_int(zz_arr[pr]; label="correlation $pr") for pr in pairs]
    pilot_tau = maximum((tauE[1], taumz[1], taumzs[1], tau_mz2[1],
                         (tau[1] for tau in tau_zz)...))
    nblocks = blocking_nblocks(nsweeps, pilot_tau; label="TFIM observables")

    ones_ = ones(length(n_arr))
    nbar, nbar_err = jackknife_ratio(n_arr, ones_; nblocks=nblocks)
    energy = m.C_shift - nbar / beta_phys
    energy_err = nbar_err / beta_phys
    mz2, mz2_err = jackknife_ratio(mz2_arr, ones_; nblocks=nblocks)
    zzres = Dict{Tuple{Int,Int},Tuple{Float64,Float64}}()
    for pr in pairs
        zzres[pr] = jackknife_ratio(zz_arr[pr], ones_; nblocks=nblocks)
    end
    accept = cl_total > 0 ? fl_total / cl_total : 0.0
    return SSEResult(energy, energy_err, mz2, mz2_err, zzres, 1.0, 0.0,
                     (tauE[1], tauE[2]), (taumz[1], taumz[2]), (taumzs[1], taumzs[2]),
                     sum(mz_arr) / length(mz_arr), count_signflips(mz_arr), accept, 1.0)
end
