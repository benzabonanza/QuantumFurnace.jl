"""
scratch_qf_bm1_verify.jl

Probes edge cases for `validate_jump_pairing` (qf-bm1 Q1).

Cases:
1. Empty jump set -> no error.
2. All Hermitian -> no error.
3. Multiple paired sets ([sigma+_1, sigma-_1, sigma+_2, sigma-_2]) -> no error.
4. Self-paired (A = A_dag, e.g., a Hermitian matrix marked hermitian=false) -> no error.
5. Three identical (A, A, A_dag) -> currently expected: no error (A and the second A
   should both find the third A_dag as their partner).
6. Mismatched pair on different sites (sigma+_1 + sigma-_2) -> error (NOT a valid pair).
7. Float64 noise: A and 1.0e-13*A_dag (within atol=1e-12) -> no error;
                 1.0e-10*A_dag (outside atol) -> error.
8. Single non-Hermitian -> error.
9. allow_unpaired_nonhermitian=true bypass -> no error even when invalid.
"""

using QuantumFurnace
using LinearAlgebra
using Test

const _S_PLUS  = ComplexF64[0 1; 0 0]
const _S_MINUS = ComplexF64[0 0; 1 0]
const _X = ComplexF64[0 1; 1 0]

# Helper: build a JumpOp with arbitrary data for testing the validation
# (we don't care about basis here — validation only inspects `.data` and `.hermitian`).
function _mkjump(data::AbstractMatrix; hermitian::Bool=false)
    eb = copy(data)  # eigenbasis copy is irrelevant for validate_jump_pairing
    return JumpOp(Matrix{ComplexF64}(data), Matrix{ComplexF64}(eb), false, hermitian)
end

"""
Build a 2-site Pauli operator for site 1 or 2, in a 4-dim Hilbert space.
op2 should be 2x2; pad on remaining site with I_2.
"""
function _twosite(op2::AbstractMatrix; site::Int)
    I2 = Matrix{ComplexF64}(I, 2, 2)
    if site == 1
        return kron(op2, I2)
    else
        return kron(I2, op2)
    end
end

passed = String[]
failed = String[]

function _record(name::AbstractString, ok::Bool, info::String="")
    if ok
        push!(passed, name * (isempty(info) ? "" : " ($(info))"))
        println("  PASS  $name")
        isempty(info) || println("        $info")
    else
        push!(failed, name * (isempty(info) ? "" : " ($(info))"))
        println("  FAIL  $name")
        isempty(info) || println("        $info")
    end
end

println("====================================================================")
println("qf-bm1 validate_jump_pairing edge-case probe")
println("====================================================================")

# -----------------------------------------------------------------
# Case 1: Empty jump set
# -----------------------------------------------------------------
println("\n[1] Empty jump set")
try
    res = validate_jump_pairing(JumpOp[])
    _record("empty jump set", res === nothing, "returned $(res)")
catch e
    _record("empty jump set", false, "threw: $(e)")
end

# -----------------------------------------------------------------
# Case 2: All Hermitian (with hermitian=true)
# -----------------------------------------------------------------
println("\n[2] All Hermitian, hermitian=true")
try
    jumps = JumpOp[_mkjump(_X; hermitian=true)]
    res = validate_jump_pairing(jumps)
    _record("single Hermitian (hermitian=true)", res === nothing)
catch e
    _record("single Hermitian (hermitian=true)", false, "threw: $(e)")
end

# -----------------------------------------------------------------
# Case 3: Multiple paired sets — sigma+_1, sigma-_1, sigma+_2, sigma-_2
# -----------------------------------------------------------------
println("\n[3] Two paired sets on different sites")
try
    jumps = JumpOp[
        _mkjump(_twosite(_S_PLUS;  site=1)),
        _mkjump(_twosite(_S_MINUS; site=1)),
        _mkjump(_twosite(_S_PLUS;  site=2)),
        _mkjump(_twosite(_S_MINUS; site=2)),
    ]
    res = validate_jump_pairing(jumps)
    _record("paired (s+_1, s-_1, s+_2, s-_2)", res === nothing)
catch e
    _record("paired (s+_1, s-_1, s+_2, s-_2)", false, "threw: $(e)")
end

# -----------------------------------------------------------------
# Case 4: Self-paired — Hermitian matrix marked hermitian=false (e.g. by mistake)
# A single such jump satisfies A = A_dag, but partner search excludes self (j == k).
# Expected: no other jump matches -> ERROR.
# -----------------------------------------------------------------
println("\n[4] Self-paired single Hermitian marked hermitian=false")
println("    Expected: ERROR (validation excludes self-match j==k)")
try
    jumps = JumpOp[_mkjump(_X; hermitian=false)]
    res = validate_jump_pairing(jumps)
    _record("self-paired single A=A_dag, hermitian=false",
        false, "should have errored, returned $(res)")
catch e
    if isa(e, ArgumentError)
        _record("self-paired single A=A_dag, hermitian=false",
            true, "correctly rejected (single jump cannot self-pair, j==k excluded)")
    else
        _record("self-paired single A=A_dag, hermitian=false",
            false, "wrong exception type: $(typeof(e))")
    end
end

# -----------------------------------------------------------------
# Case 4b: Two copies of A=A_dag, both hermitian=false
# Expected: each finds the OTHER -> no error
# -----------------------------------------------------------------
println("\n[4b] Two copies of A=A_dag, both hermitian=false")
println("     Expected: each finds the other -> no error")
try
    jumps = JumpOp[_mkjump(_X; hermitian=false), _mkjump(_X; hermitian=false)]
    res = validate_jump_pairing(jumps)
    _record("two copies of A=A_dag, hermitian=false", res === nothing)
catch e
    _record("two copies of A=A_dag, hermitian=false", false, "threw: $(e)")
end

# -----------------------------------------------------------------
# Case 5: Three identical (A, A, A_dag) where A is non-Hermitian
# Expected: BOTH As at index 1 and 2 should find A_dag at index 3. The implementation
# breaks on first match — but j != k, so first match for k=1 is j=3, first match for
# k=2 is also j=3. Both find the same A_dag. Should return nothing.
# -----------------------------------------------------------------
println("\n[5] Three jumps: (A, A, A_dag) with A non-Hermitian")
println("    Expected: both As find the same A_dag at j=3 -> no error")
try
    A = _twosite(_S_PLUS; site=1)
    Adag = A'
    jumps = JumpOp[_mkjump(A), _mkjump(A), _mkjump(Adag)]
    res = validate_jump_pairing(jumps)
    _record("(A, A, A_dag)", res === nothing,
        "both unpaired As should match A_dag at j=3")
catch e
    _record("(A, A, A_dag)", false, "threw: $(e)")
end

# -----------------------------------------------------------------
# Case 6: Mismatched pair on different sites: sigma+_1 + sigma-_2 (NOT a valid pair).
# (sigma+_1)_dag = sigma-_1, NOT sigma-_2. So both jumps are unpaired.
# -----------------------------------------------------------------
println("\n[6] Mismatched pair (s+_1, s-_2) — different sites")
println("    Expected: ERROR (s+_1 has no s-_1 partner)")
try
    jumps = JumpOp[
        _mkjump(_twosite(_S_PLUS;  site=1)),
        _mkjump(_twosite(_S_MINUS; site=2)),
    ]
    res = validate_jump_pairing(jumps)
    _record("(s+_1, s-_2)",
        false, "should have errored, returned $(res)")
catch e
    if isa(e, ArgumentError)
        _record("(s+_1, s-_2)", true, "correctly rejected as unpaired")
    else
        _record("(s+_1, s-_2)", false, "wrong exception: $(typeof(e))")
    end
end

# -----------------------------------------------------------------
# Case 7a: Float64 noise within atol — A and (1+1e-13)*A_dag
# atol = 1e-12 default. Diff is O(1e-13) -> within atol -> should pass.
# -----------------------------------------------------------------
println("\n[7a] Float64 noise within atol=1e-12 — A vs (1 + 1e-13)*A_dag")
try
    A = _twosite(_S_PLUS; site=1)
    Adag_noisy = (1.0 + 1e-13) * A'
    jumps = JumpOp[_mkjump(A), _mkjump(Adag_noisy)]
    res = validate_jump_pairing(jumps)  # default atol=1e-12
    _record("A vs (1+1e-13)*A_dag, atol=1e-12", res === nothing)
catch e
    _record("A vs (1+1e-13)*A_dag, atol=1e-12", false, "threw: $(e)")
end

# -----------------------------------------------------------------
# Case 7b: Float64 noise outside atol — A and (1+1e-10)*A_dag
# Diff is O(1e-10) > atol=1e-12 -> should ERROR.
# -----------------------------------------------------------------
println("\n[7b] Float64 noise outside atol=1e-12 — A vs (1 + 1e-10)*A_dag")
println("     Expected: ERROR (diff exceeds default atol)")
try
    A = _twosite(_S_PLUS; site=1)
    Adag_noisy = (1.0 + 1e-10) * A'
    jumps = JumpOp[_mkjump(A), _mkjump(Adag_noisy)]
    res = validate_jump_pairing(jumps; atol=1e-12)
    _record("A vs (1+1e-10)*A_dag, atol=1e-12",
        false, "should have errored")
catch e
    if isa(e, ArgumentError)
        _record("A vs (1+1e-10)*A_dag, atol=1e-12", true, "correctly rejected")
    else
        _record("A vs (1+1e-10)*A_dag, atol=1e-12", false, "wrong exception: $(typeof(e))")
    end
end

# -----------------------------------------------------------------
# Case 7c: Same as 7b but explicit atol=1e-9 -> should PASS
# -----------------------------------------------------------------
println("\n[7c] Float64 noise (1e-10) but atol=1e-9 — should PASS")
try
    A = _twosite(_S_PLUS; site=1)
    Adag_noisy = (1.0 + 1e-10) * A'
    jumps = JumpOp[_mkjump(A), _mkjump(Adag_noisy)]
    res = validate_jump_pairing(jumps; atol=1e-9)
    _record("A vs (1+1e-10)*A_dag, atol=1e-9", res === nothing)
catch e
    _record("A vs (1+1e-10)*A_dag, atol=1e-9", false, "threw: $(e)")
end

# -----------------------------------------------------------------
# Case 8: Single non-Hermitian jump — should ERROR
# -----------------------------------------------------------------
println("\n[8] Single non-Hermitian sigma+ — should ERROR")
try
    jumps = JumpOp[_mkjump(_twosite(_S_PLUS; site=1))]
    res = validate_jump_pairing(jumps)
    _record("single non-Hermitian", false, "should have errored")
catch e
    if isa(e, ArgumentError)
        _record("single non-Hermitian", true, "correctly rejected")
    else
        _record("single non-Hermitian", false, "wrong exception: $(typeof(e))")
    end
end

# -----------------------------------------------------------------
# Case 9: allow_unpaired_nonhermitian=true bypass
# -----------------------------------------------------------------
println("\n[9] allow_unpaired_nonhermitian=true bypass")
try
    jumps = JumpOp[_mkjump(_twosite(_S_PLUS; site=1))]
    res = validate_jump_pairing(jumps; allow_unpaired_nonhermitian=true)
    _record("bypass with allow_unpaired_nonhermitian=true", res === nothing)
catch e
    _record("bypass with allow_unpaired_nonhermitian=true", false, "threw: $(e)")
end

# -----------------------------------------------------------------
# Case 10: A jump and its (mistyped) Hermitian partner — A non-Hermitian,
# A_dag in the set but marked hermitian=true. Validation looks at jumps[k].hermitian
# only. The non-Hermitian A still searches for ANY jump whose data ≈ A_dag.
# It should find the marked-Hermitian A_dag and return nothing.
# -----------------------------------------------------------------
println("\n[10] A non-Herm, A_dag marked hermitian=true (mixed)")
println("     Expected: pass (validation only checks data, not the partner's flag)")
try
    A = _twosite(_S_PLUS; site=1)
    jumps = JumpOp[
        _mkjump(A; hermitian=false),
        _mkjump(A'; hermitian=true),
    ]
    res = validate_jump_pairing(jumps)
    _record("mixed flags (s+ NH, s- marked Herm)", res === nothing)
catch e
    _record("mixed flags (s+ NH, s- marked Herm)", false, "threw: $(e)")
end

# -----------------------------------------------------------------
# Case 11: Differently-sized non-Hermitian operators
# A is 4x4 and "partner" is 2x2 with same data — sizes mismatch. Should error.
# -----------------------------------------------------------------
println("\n[11] Size mismatch — A 4x4, fake partner 2x2")
println("     Expected: ERROR (different shapes)")
try
    A = _twosite(_S_PLUS; site=1)  # 4x4
    fake = _S_MINUS  # 2x2
    jumps = JumpOp[_mkjump(A), _mkjump(fake)]
    res = validate_jump_pairing(jumps)
    _record("size mismatch 4x4 vs 2x2", false, "should have errored")
catch e
    if isa(e, ArgumentError)
        _record("size mismatch 4x4 vs 2x2", true, "correctly rejected")
    else
        _record("size mismatch 4x4 vs 2x2", false, "wrong exception: $(typeof(e))")
    end
end

# -----------------------------------------------------------------
# Case 12: All-zero non-Hermitian jump (edge of all-zero matrix is its own dag)
# A=zeros, A_dag=zeros -> A and A_dag indistinguishable, partner search wants
# "any other jump whose data ≈ A_dag". With ONE jump, no other -> error.
# -----------------------------------------------------------------
println("\n[12] Single all-zero non-Hermitian jump")
println("     Expected: ERROR (no other jump to pair with)")
try
    jumps = JumpOp[_mkjump(zeros(ComplexF64, 4, 4); hermitian=false)]
    res = validate_jump_pairing(jumps)
    _record("single zeros 4x4 non-Herm", false, "should have errored (no other)")
catch e
    if isa(e, ArgumentError)
        _record("single zeros 4x4 non-Herm", true, "correctly rejected")
    else
        _record("single zeros 4x4 non-Herm", false, "wrong exception: $(typeof(e))")
    end
end

# -----------------------------------------------------------------
# Case 12b: Two all-zero non-Hermitian jumps -> they "pair" with each other.
# Edge case: a degenerate but technically "paired" set.
# -----------------------------------------------------------------
println("\n[12b] Two all-zero non-Hermitian jumps")
println("      Expected: pass (each pairs with the other)")
try
    jumps = JumpOp[
        _mkjump(zeros(ComplexF64, 4, 4); hermitian=false),
        _mkjump(zeros(ComplexF64, 4, 4); hermitian=false),
    ]
    res = validate_jump_pairing(jumps)
    _record("two zeros 4x4 non-Herm", res === nothing,
        "degenerate but algorithmically paired")
catch e
    _record("two zeros 4x4 non-Herm", false, "threw: $(e)")
end

# -----------------------------------------------------------------
# Summary
# -----------------------------------------------------------------
println("\n====================================================================")
println("SUMMARY: $(length(passed)) PASS, $(length(failed)) FAIL")
println("====================================================================")
for s in passed
    println("  PASS: $s")
end
for s in failed
    println("  FAIL: $s")
end
exit(length(failed) == 0 ? 0 : 1)
