#=
Phase B prototype: NUFFT-based DLL TimeDomain Lindblad operator at ω = 0.

The DLL paper Eq. 3.4 (third form) on the simulator's truncated time grid is

    L_a[i, j] = A_eb[i, j] · τ · Σ_m f(t_m) · cis((λ_i − λ_j) · t_m)

This is exactly a 1D type-3 NUFFT of `f(t_m)` evaluated at non-uniform Bohr-
frequency targets `{λ_i − λ_j}`, which is precisely what
`_prepare_oft_nufft_prefactors` computes when called with `energy_labels = [0]`
(then `cis(-omega · t_m) = 1`). FINUFFT precision: ε = 1e-12.

Verification: NUFFT path agrees with explicit `dll_lindblad_op_time` to ε ≤ 1e-12.

Run: julia --project scripts/scratch_dll_dissipator_nufft.jl
=#

using QuantumFurnace
using LinearAlgebra
using Printf

# ----- Setup -------------------------------------------------------------
n = 5
beta = 5.0
ham = load_hamiltonian("heis", n; beta=beta)

X = ComplexF64[0 1; 1 0]
Z = ComplexF64[1 0; 0 -1]
Id = ComplexF64[1 0; 0 1]
function _pad1(op)
    out = ComplexF64.(op)
    for _ in 2:n
        out = kron(out, Id)
    end
    return out
end
A1 = _pad1(X)
A2 = _pad1(Z)
B1 = ham.eigvecs' * A1 * ham.eigvecs
B2 = ham.eigvecs' * A2 * ham.eigvecs
jumps = JumpOp[JumpOp(A1, B1, true, true), JumpOp(A2, B2, true, true)]

filter = DLLGaussianFilter(beta)

# Build the simulator's time grid (matches src/furnace_utensils.jl path).
N_BITS = 12
W0 = 0.05
T0 = 2π / (2^N_BITS * W0)
N = 2^N_BITS
raw_time_labels = collect((-N÷2):(N÷2 - 1)) .* T0
time_labels = QuantumFurnace._truncate_time_labels_for_oft(raw_time_labels, 2/beta; filter=filter)
@printf "Truncated grid: Nt = %d  (raw N = %d)\n" length(time_labels) N

# ----- Reference: explicit Riemann sum via dll_lindblad_op_time ----------
println("\nBuilding L_a via explicit triple loop `dll_lindblad_op_time`...")
# Warm-up
_ = QuantumFurnace.dll_lindblad_op_time(jumps[1], ham, time_labels, filter, T0)
L_explicit = Matrix{ComplexF64}[]
t_explicit = @elapsed begin
    for jump in jumps
        push!(L_explicit, QuantumFurnace.dll_lindblad_op_time(jump, ham, time_labels, filter, T0))
    end
end
@printf "  explicit Riemann: %.3f s for %d jumps\n" t_explicit length(jumps)

# ----- NUFFT path: single-slice at ω = 0 ---------------------------------
println("\nBuilding L_a via NUFFT slice at ω = 0...")
T = eltype(ham.eigvals)
# Warm-up to avoid measuring JIT compilation
_ = QuantumFurnace._prepare_oft_nufft_prefactors(
    ham.bohr_freqs, time_labels, T[zero(T)], filter; eps=1e-12,
)
t_setup = @elapsed begin
    nufft = QuantumFurnace._prepare_oft_nufft_prefactors(
        ham.bohr_freqs, time_labels, T[zero(T)], filter; eps=1e-12,
    )
end
pf_zero = nufft.data[:, :, 1]
@printf "  NUFFT setup (one ω-slice): %.4f s\n" t_setup

# Also warmup explicit
_ = QuantumFurnace.dll_lindblad_op_time(jumps[1], ham, time_labels, filter, T0)

dim = size(ham.data, 1)
L_buf = Matrix{ComplexF64}(undef, dim, dim)
# Warm-up the elementwise multiply
@. L_buf = jumps[1].in_eigenbasis * pf_zero * T0

L_nufft = Matrix{ComplexF64}[]
t_nufft = @elapsed begin
    for jump in jumps
        L_a = Matrix{ComplexF64}(undef, dim, dim)
        @. L_a = jump.in_eigenbasis * pf_zero * T0
        push!(L_nufft, L_a)
    end
end
@printf "  NUFFT mul (per jump): %.6f s for %d jumps\n" t_nufft length(jumps)

# ----- Cross-check ----------------------------------------------------------
println("\nCross-check NUFFT vs explicit:")
for k in 1:length(jumps)
    err = opnorm(L_explicit[k] - L_nufft[k])
    rel = err / opnorm(L_explicit[k])
    @printf "  jump %d:  ‖L_exp − L_nufft‖_op = %.3e   relative = %.3e\n" k err rel
end

# ----- Speedup ------------------------------------------------------------
total_explicit = t_explicit
total_nufft = t_setup + t_nufft
@printf "\nTotal: explicit %.3f s   vs   NUFFT %.3f s   (%.1fx speedup)\n" total_explicit total_nufft total_explicit / total_nufft

println("\nDONE")
