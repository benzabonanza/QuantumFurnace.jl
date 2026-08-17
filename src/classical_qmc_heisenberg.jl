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

Reconstruct the 1D Heisenberg SSE model in physical units.

# Keywords
- `seed`: Seed for the fixture-compatible disorder draws.
- `periodic`: Whether to include the wraparound bond.
- `disorder_strength`: Scale of the longitudinal and bond disorder.
- `J`: Exchange coupling.

# Returns
An `HeisModel` containing the SSE data and dense physical Hamiltonian.
"""
function build_sse_heis_model(n::Int; seed::Int, periodic::Bool, disorder_strength::Float64, J::Float64=1.0)
    n >= 2 || throw(ArgumentError("n must be at least 2."))
    isfinite(J) && J > 0 || throw(ArgumentError(
        "J must be finite and > 0 for the antiferromagnetic SSE decomposition."))
    isfinite(disorder_strength) && disorder_strength >= 0 || throw(ArgumentError(
        "disorder_strength must be finite and >= 0."))

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

    # Math: $H = J sum_(ij)(X_i X_j + Y_i Y_j + Z_i Z_j)
    # + sum_i h_i Z_i + sum_(ij) K_(ij) Z_i Z_j$.
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

sse_reconstruction_error(model::HeisModel, raw) = maximum(abs.(model.H_phys .- fixture_H_phys(raw)))

sse_exact_reference(model::HeisModel, beta_phys::Float64; pairs::Vector{Tuple{Int,Int}}) =
    _exact_reference(model.H_phys, model.n, beta_phys; pairs=pairs)

function validate_model(m::HeisModel)
    m.n >= 2 || throw(ArgumentError("Heisenberg SSE models require at least two sites."))
    isfinite(m.J) && m.J > 0 || throw(ArgumentError(
        "Heisenberg SSE requires finite antiferromagnetic J > 0."))
    length(m.field) == m.n || throw(DimensionMismatch("Heisenberg field length must equal n."))
    length(m.zz) == length(m.bonds) || throw(DimensionMismatch(
        "Heisenberg bond-disorder length must match bonds."))
    length(m.frustrated) == length(m.bonds) || throw(DimensionMismatch(
        "Heisenberg frustration flags must match bonds."))
    all(isfinite, m.field) && all(isfinite, m.zz) || throw(ArgumentError(
        "Heisenberg SSE couplings must be finite."))
    isfinite(m.C_shift) || throw(ArgumentError("Heisenberg SSE energy shift must be finite."))
    for (i, j) in m.bonds
        1 <= i <= m.n && 1 <= j <= m.n && i != j || throw(ArgumentError(
            "Heisenberg bonds must join distinct sites in 1:$(m.n)."))
    end
    any(m.frustrated) && throw(ArgumentError(
        "Heisenberg SSE is restricted to bipartite exchange graphs; odd periodic chains are unsupported."))
    return nothing
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

# Heisenberg operator-loop update. Each exchange vertex stores bottom-site
# legs 1--2 and top-site legs 3--4; deterministic pairing is 1<->2, 3<->4.
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

Run the SSE sampler at physical inverse temperature `beta_phys`.

# Keywords
- `pairs`: Site pairs for reported `Z_i Z_j` correlations.
- `nsweeps`: Number of measured sweeps.
- `nwarm`: Number of equilibration sweeps.
- `seed`: Sampling seed, independent of the model's disorder seed.

# Returns
An `SSEResult` with observables and sampling diagnostics.

The Heisenberg method is not valid for the frustrated odd-size periodic sign
sector described in the module documentation.
"""
function run_sse(m::HeisModel, beta_phys::Float64; pairs::Vector{Tuple{Int,Int}},
                        nsweeps::Int=200_000, nwarm::Int=40_000, seed::Int=12345)
    validate_model(m)
    validate_run_inputs(m.n, beta_phys, pairs, nsweeps, nwarm)

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

    tauE = tau_int(n_arr; label="expansion order")
    taumz = tau_int(absmz_arr; label="absolute magnetisation")
    taumzs = tau_int(mz_arr; label="signed magnetisation")
    tau_mz2 = tau_int(mz2_arr; label="squared magnetisation")
    tau_zz = [tau_int(zz_arr[pr]; label="correlation $pr") for pr in pairs]
    pilot_tau = maximum((tauE[1], taumz[1], taumzs[1], tau_mz2[1],
                         (tau[1] for tau in tau_zz)...))
    nblocks = blocking_nblocks(nsweeps, pilot_tau; label="Heisenberg observables")

    # Sign-reweighted estimators: <O> = <O*s>/<s>.
    avg_sign, avg_sign_err = jackknife_ratio(
        sign_arr, ones(length(sign_arr)); nblocks=nblocks)
    nbar, nbar_err = jackknife_ratio(n_arr .* sign_arr, sign_arr; nblocks=nblocks)
    energy = m.C_shift - nbar / beta_phys
    energy_err = nbar_err / beta_phys
    mz2, mz2_err = jackknife_ratio(mz2_arr .* sign_arr, sign_arr; nblocks=nblocks)
    zzres = Dict{Tuple{Int,Int},Tuple{Float64,Float64}}()
    for pr in pairs
        zzres[pr] = jackknife_ratio(
            zz_arr[pr] .* sign_arr, sign_arr; nblocks=nblocks)
    end

    loop_accept = loops_total > 0 ? flips_total / loops_total : 0.0
    mean_r = loops_total > 0 ? r_total / loops_total : 1.0
    return SSEResult(energy, energy_err, mz2, mz2_err, zzres, avg_sign, avg_sign_err,
                     (tauE[1], tauE[2]), (taumz[1], taumz[2]), (taumzs[1], taumzs[2]),
                     sum(mz_arr) / length(mz_arr), count_signflips(mz_arr), loop_accept, mean_r)
end
