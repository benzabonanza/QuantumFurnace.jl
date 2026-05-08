# scripts/scratch_omega_threading_threshold.jl
# qf-in3.4: numerically sweep OMEGA_THREAD_THRESHOLD for the Lindbladian
# matvec (apply_lindbladian!) and the Channel rho_jump accumulator
# (_accumulate_rho_jump!), to choose an evidence-backed crossover.
#
# Approach: time the threaded vs serial helpers directly on a Workspace at
# varying effective work-list sizes. The threaded path has fixed startup
# overhead (task spawn, BLAS toggle, per-task KrylovScratch allocation);
# serial has none. The crossover is where work-per-label × n_labels exceeds
# threading overhead.
#
# Run:
#   JULIA_NUM_THREADS=4 julia --project scripts/scratch_omega_threading_threshold.jl

using QuantumFurnace
using LinearAlgebra
using Random
using Printf

include(joinpath(@__DIR__, "..", "test", "test_helpers.jl"))

# Lightweight @elapsed-based timing — minimum of `samples` runs.
function _bench(f, samples::Int=200)
    best = Inf
    for _ in 1:samples
        t = @elapsed f()
        best = min(best, t)
    end
    return best * 1e3  # ms
end

# ---------------------------------------------------------------------------
# Helper: force-call serial vs threaded apply_lindbladian! variants on a
# truncated energy-labels view. We do NOT mutate the Workspace; we run the
# private helpers explicitly with a chosen subset of energy_labels.
# ---------------------------------------------------------------------------

function _serial_apply!(ws, rho, config, ham, energy_labels_view, prefactor, inv_4sigma2)
    sc = ws.scratch
    G_left = ws.G_left
    G_right = ws.G_right
    fill!(sc.rho_out, 0)
    mul!(sc.rho_out, G_left, rho)
    mul!(sc.rho_out, rho, G_right, 1.0, 1.0)

    bohr_freqs = ham.bohr_freqs
    @inbounds for (k, eigenbasis) in enumerate(ws.jump_eigenbases)
        is_herm = ws.jump_hermitian[k]
        if is_herm
            for w_raw in energy_labels_view
                w_raw > 1e-12 && continue
                w = abs(w_raw)
                QuantumFurnace.oft!(sc.jump_oft, eigenbasis, bohr_freqs, w, inv_4sigma2)
                scalar_w = prefactor * QuantumFurnace.pick_transition(config, w)
                QuantumFurnace._accumulate_sandwich_scratch!(sc.rho_out, sc.jump_oft, rho, scalar_w,
                                                             sc.sandwich_tmp, sc.sandwich_out)
                if w > 1e-12
                    scalar_neg = prefactor * QuantumFurnace.pick_transition(config, -w)
                    QuantumFurnace._accumulate_sandwich_adj_scratch!(sc.rho_out, sc.jump_oft, rho, scalar_neg,
                                                                    sc.sandwich_tmp, sc.sandwich_out)
                end
            end
        else
            for w in energy_labels_view
                QuantumFurnace.oft!(sc.jump_oft, eigenbasis, bohr_freqs, w, inv_4sigma2)
                scalar_w = prefactor * QuantumFurnace.pick_transition(config, w)
                QuantumFurnace._accumulate_sandwich_scratch!(sc.rho_out, sc.jump_oft, rho, scalar_w,
                                                             sc.sandwich_tmp, sc.sandwich_out)
            end
        end
    end
    return sc.rho_out
end

function _threaded_apply!(ws, rho, config, ham, energy_labels_view, prefactor, inv_4sigma2)
    sc = ws.scratch
    G_left = ws.G_left
    G_right = ws.G_right
    fill!(sc.rho_out, 0)
    mul!(sc.rho_out, G_left, rho)
    mul!(sc.rho_out, rho, G_right, 1.0, 1.0)

    QuantumFurnace._apply_lindbladian_threaded_energy!(
        sc, rho, ws.jump_eigenbases, ws.jump_hermitian,
        ham.bohr_freqs, Vector{Float64}(energy_labels_view),
        config, prefactor, inv_4sigma2; adjoint=false)
    return sc.rho_out
end

# ---------------------------------------------------------------------------
# Build a high-resolution Workspace and sweep the labels-window size.
# ---------------------------------------------------------------------------

function run_threshold_sweep(num_qubits::Int)
    Threads.nthreads() == 1 && error("Re-run with JULIA_NUM_THREADS>=2")

    sys = make_test_system(; num_qubits=num_qubits)
    ham = sys.hamiltonian
    jumps = sys.jumps

    config = make_config(Lindbladian(), EnergyDomain(); construction=KMS(),
                         num_qubits=num_qubits)
    ws = Workspace(config, ham, jumps)
    rho = Matrix(random_density_matrix(num_qubits))
    prefactor = ws.oft_domain_prefactor * ws.gamma_norm_factor
    inv_4sigma2 = 1.0 / (4 * config.sigma^2)

    base = ws.energy_labels
    full = vcat(base, base, base, base, base, base, base, base)

    candidates = [5, 10, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200, 300, length(full)]
    println("\n# OMEGA_THREAD_THRESHOLD sweep (Energy) — n=$num_qubits, dim=$(2^num_qubits), nthreads=$(Threads.nthreads())")
    println("# n_labels   serial_ms   threaded_ms   ratio (s/t)   recommend")
    println("# -------- ----------- ------------- -------------- ------------")

    for N in candidates
        N = min(N, length(full))
        labels = full[1:N]
        _serial_apply!(ws, rho, config, ham, labels, prefactor, inv_4sigma2)
        _threaded_apply!(ws, rho, config, ham, labels, prefactor, inv_4sigma2)
        ts = _bench(() -> _serial_apply!(ws, rho, config, ham, labels, prefactor, inv_4sigma2))
        tt = _bench(() -> _threaded_apply!(ws, rho, config, ham, labels, prefactor, inv_4sigma2))
        ratio = ts / tt
        rec = ratio >= 1.0 ? "thread" : "serial"
        @printf("  %5d   %9.4f   %11.4f   %10.3fx   %s\n", N, ts, tt, ratio, rec)
    end
end

# TimeDomain sweep: NUFFT-lookup hot loop is ~2× lighter than the OFT
# `exp()` per matrix entry of EnergyDomain, so the crossover should shift up.
function _serial_apply_time!(ws, rho, config, ham, energy_labels_view, prefactor)
    sc = ws.scratch
    nufft = ws.oft_nufft_prefactors
    nufft_data = nufft.data
    nufft_idx = nufft.energy_to_index
    G_left = ws.G_left
    G_right = ws.G_right
    fill!(sc.rho_out, 0)
    mul!(sc.rho_out, G_left, rho)
    mul!(sc.rho_out, rho, G_right, 1.0, 1.0)
    @inbounds for (k, eigenbasis) in enumerate(ws.jump_eigenbases)
        is_herm = ws.jump_hermitian[k]
        for w_raw in energy_labels_view
            if is_herm && w_raw > 1e-12
                continue
            end
            w = is_herm ? abs(w_raw) : w_raw
            li_idx = is_herm ? nufft_idx[w] : findfirst(==(w_raw), ws.energy_labels)
            li_idx === nothing && continue
            mat = @view nufft_data[:, :, li_idx]
            @. sc.jump_oft = eigenbasis * mat
            scalar_w = prefactor * QuantumFurnace.pick_transition(config, w)
            QuantumFurnace._accumulate_sandwich_scratch!(sc.rho_out, sc.jump_oft, rho, scalar_w,
                                                         sc.sandwich_tmp, sc.sandwich_out)
            if is_herm && w > 1e-12
                scalar_neg = prefactor * QuantumFurnace.pick_transition(config, -w)
                QuantumFurnace._accumulate_sandwich_adj_scratch!(sc.rho_out, sc.jump_oft, rho, scalar_neg,
                                                                sc.sandwich_tmp, sc.sandwich_out)
            end
        end
    end
    return sc.rho_out
end

function _threaded_apply_time!(ws, rho, config, ham, energy_labels_view, prefactor)
    sc = ws.scratch
    nufft = ws.oft_nufft_prefactors
    G_left = ws.G_left
    G_right = ws.G_right
    fill!(sc.rho_out, 0)
    mul!(sc.rho_out, G_left, rho)
    mul!(sc.rho_out, rho, G_right, 1.0, 1.0)
    QuantumFurnace._apply_lindbladian_threaded_timetrot!(
        sc, rho, ws.jump_eigenbases, ws.jump_hermitian,
        nufft.data, nufft.energy_to_index,
        Vector{Float64}(energy_labels_view), config, prefactor; adjoint=false)
    return sc.rho_out
end

function run_threshold_sweep_time(num_qubits::Int)
    sys = make_test_system(; num_qubits=num_qubits)
    ham = sys.hamiltonian
    jumps = sys.jumps

    config = make_config(Lindbladian(), TimeDomain(); construction=KMS(),
                         num_qubits=num_qubits)
    ws = Workspace(config, ham, jumps)
    rho = Matrix(random_density_matrix(num_qubits))
    prefactor = ws.oft_domain_prefactor * ws.gamma_norm_factor
    base = ws.energy_labels
    full = vcat(base, base, base, base, base, base, base, base)
    candidates = [5, 10, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200, 300, length(full)]
    println("\n# OMEGA_THREAD_THRESHOLD sweep (Time)   — n=$num_qubits, dim=$(2^num_qubits), nthreads=$(Threads.nthreads())")
    println("# n_labels   serial_ms   threaded_ms   ratio (s/t)   recommend")
    println("# -------- ----------- ------------- -------------- ------------")
    for N in candidates
        N = min(N, length(full))
        labels = full[1:N]
        _serial_apply_time!(ws, rho, config, ham, labels, prefactor)
        _threaded_apply_time!(ws, rho, config, ham, labels, prefactor)
        ts = _bench(() -> _serial_apply_time!(ws, rho, config, ham, labels, prefactor))
        tt = _bench(() -> _threaded_apply_time!(ws, rho, config, ham, labels, prefactor))
        ratio = ts / tt
        rec = ratio >= 1.0 ? "thread" : "serial"
        @printf("  %5d   %9.4f   %11.4f   %10.3fx   %s\n", N, ts, tt, ratio, rec)
    end
end

for n in (3, 4, 5)
    run_threshold_sweep(n)
end
for n in (3, 4, 5)
    run_threshold_sweep_time(n)
end
