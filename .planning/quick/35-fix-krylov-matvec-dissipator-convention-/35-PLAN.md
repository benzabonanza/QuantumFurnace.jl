---
phase: quick-35
plan: 01
type: execute
wave: 1
depends_on: []
files_modified:
  - src/krylov_matvec.jl
  - test/test_krylov_matvec.jl
autonomous: true
must_haves:
  truths:
    - "Krylov matvec dissipator matches dense kron-based convention for complex jump operators"
    - "All existing round-trip tests still pass (real/Hermitian jumps unaffected)"
    - "New test with complex non-Hermitian jump operator passes round-trip vs dense"
  artifacts:
    - path: "src/krylov_matvec.jl"
      provides: "Fixed dissipator accumulation functions"
    - path: "test/test_krylov_matvec.jl"
      provides: "Complex jump operator round-trip test"
  key_links:
    - from: "src/krylov_matvec.jl"
      to: "src/qi_tools.jl"
      via: "dissipator convention agreement"
      pattern: "conj.*rho.*transpose"
---

<objective>
Fix all 6 Krylov matvec dissipator helper functions to use the dense kron-based convention
`conj(J) * rho * J^T` instead of the textbook `J * rho * J'` for the sandwich term.

Purpose: The dense vectorized code (`_vectorize_liouv_diss_and_add!`) uses kron(J, conj(J))
which un-vectorizes to `conj(J) rho J^T`. The Krylov matvec currently uses `L rho L'`.
These disagree for general complex (non-Hermitian) jump operators. The existing tests use
Pauli X/Y/Z which are all Hermitian, masking the bug.

Output: Fixed `krylov_matvec.jl` + new test with complex non-Hermitian jump operator.
</objective>

<execution_context>
@./.claude/get-shit-done/workflows/execute-plan.md
@./.claude/get-shit-done/templates/summary.md
</execution_context>

<context>
@src/krylov_matvec.jl
@src/qi_tools.jl (lines 64-103: _vectorize_liouv_diss_and_add! -- the target convention)
@test/test_krylov_matvec.jl
@test/test_helpers.jl
@src/krylov_workspace.jl (KrylovWorkspace struct: scratch buffers tmp1, tmp2, LdagL)
</context>

<tasks>

<task type="auto">
  <name>Task 1: Fix all 6 dissipator helper functions in krylov_matvec.jl</name>
  <files>src/krylov_matvec.jl</files>
  <action>
Fix the sandwich term in all 6 dissipator accumulation functions to match the dense kron convention.
The anticommutator terms are correct in all functions and must NOT be changed.

The key identity: the dense code uses `kron(J, conj(J))` which un-vectorizes to `conj(J) * rho * J^T`.
For the HS adjoint, the adjoint of `X rho Y` is `X^H rho Y^H`, so the adjoint sandwich of
`conj(J) rho J^T` is `conj(J)^H rho (J^T)^H = J^T rho conj(J)`.

For the `adj_L` variants (operator is `M = L^H`), apply the same convention with M in place of J:
`conj(M) rho M^T = conj(L^H) rho (L^H)^T = L^T rho conj(L)`.

Strategy for each function: compute anticommutator terms FIRST (using LdagL), then reuse LdagL
as scratch for the sandwich result, avoiding extra allocations.

**1. `_accumulate_dissipator!` (forward, operator L):**
Change sandwich from `L * rho * L'` to `conj(L) * rho * L^T`.
New BLAS sequence (reordered: anticommutator first, then sandwich):
```
# L'L -> ws.LdagL  (unchanged)
BLAS.gemm!('C', 'N', CT, L_op, L_op, ZT, ws.LdagL)
# Term 2: -0.5 * scalar * L'L * rho
BLAS.gemm!('N', 'N', CT, ws.LdagL, rho, ZT, ws.tmp1)
BLAS.axpy!(T(-0.5 * scalar), ws.tmp1, out)
# Term 3: -0.5 * scalar * rho * L'L
BLAS.gemm!('N', 'N', CT, rho, ws.LdagL, ZT, ws.tmp1)
BLAS.axpy!(T(-0.5 * scalar), ws.tmp1, out)
# Term 1: scalar * conj(L) * rho * L^T  (LdagL now free as scratch)
@. ws.tmp2 = conj(L_op)                                    # tmp2 = conj(L)
BLAS.gemm!('N', 'N', CT, ws.tmp2, rho, ZT, ws.tmp1)        # tmp1 = conj(L) * rho
BLAS.gemm!('N', 'T', CT, ws.tmp1, L_op, ZT, ws.LdagL)     # LdagL = conj(L) * rho * L^T
BLAS.axpy!(T(scalar), ws.LdagL, out)
```

**2. `_accumulate_dissipator_adj_L!` (forward, operator L' i.e. L^H):**
Change sandwich from `L' * rho * L` to `conj(L^H) * rho * (L^H)^T = L^T * rho * conj(L)`.
Anticommutator: `(L^H)^H * L^H = L * L'` (unchanged -- current code already computes L*L').
```
# L * L' -> ws.LdagL  (unchanged: (L')^H * L' = L * L')
BLAS.gemm!('N', 'C', CT, L_op, L_op, ZT, ws.LdagL)
# Term 2: -0.5 * scalar * L L' * rho
BLAS.gemm!('N', 'N', CT, ws.LdagL, rho, ZT, ws.tmp1)
BLAS.axpy!(T(-0.5 * scalar), ws.tmp1, out)
# Term 3: -0.5 * scalar * rho * L L'
BLAS.gemm!('N', 'N', CT, rho, ws.LdagL, ZT, ws.tmp1)
BLAS.axpy!(T(-0.5 * scalar), ws.tmp1, out)
# Term 1: scalar * L^T * rho * conj(L)  (LdagL now free)
@. ws.tmp2 = conj(L_op)                                    # tmp2 = conj(L)
BLAS.gemm!('T', 'N', CT, L_op, rho, ZT, ws.tmp1)           # tmp1 = L^T * rho
BLAS.gemm!('N', 'N', CT, ws.tmp1, ws.tmp2, ZT, ws.LdagL)  # LdagL = L^T * rho * conj(L)
BLAS.axpy!(T(scalar), ws.LdagL, out)
```

**3. `_accumulate_adjoint_dissipator!` (HS adjoint, operator L):**
Change sandwich from `L' * rho * L` to HS adjoint of `conj(L) rho L^T` = `J^T * rho * conj(J)`.
Wait -- more carefully: forward sandwich is `conj(L) rho L^T`. HS adjoint of `X rho Y` is `X^H rho Y^H`.
X = conj(L), X^H = conj(L)^H = L^T.  Y = L^T, Y^H = (L^T)^H = conj(L).
So adjoint sandwich: `L^T * rho * conj(L)`.
Anticommutator: L'L (unchanged, same as forward).
```
# L'L -> ws.LdagL  (unchanged)
BLAS.gemm!('C', 'N', CT, L_op, L_op, ZT, ws.LdagL)
# Term 2: -0.5 * scalar * L'L * rho  (unchanged)
BLAS.gemm!('N', 'N', CT, ws.LdagL, rho, ZT, ws.tmp1)
BLAS.axpy!(T(-0.5 * scalar), ws.tmp1, out)
# Term 3: -0.5 * scalar * rho * L'L  (unchanged)
BLAS.gemm!('N', 'N', CT, rho, ws.LdagL, ZT, ws.tmp1)
BLAS.axpy!(T(-0.5 * scalar), ws.tmp1, out)
# Term 1: scalar * L^T * rho * conj(L)  (LdagL now free)
@. ws.tmp2 = conj(L_op)                                    # tmp2 = conj(L)
BLAS.gemm!('T', 'N', CT, L_op, rho, ZT, ws.tmp1)           # tmp1 = L^T * rho
BLAS.gemm!('N', 'N', CT, ws.tmp1, ws.tmp2, ZT, ws.LdagL)  # LdagL = L^T * rho * conj(L)
BLAS.axpy!(T(scalar), ws.LdagL, out)
```

**4. `_accumulate_adjoint_dissipator_adj_L!` (HS adjoint, operator L' i.e. L^H):**
Forward for L^H gives sandwich: `L^T * rho * conj(L)` (from function 2 above).
HS adjoint of that: X = L^T, X^H = conj(L). Y = conj(L), Y^H = L^T.
So adjoint sandwich: `conj(L) * rho * L^T`.
Anticommutator: L * L' (unchanged, same as function 2).
```
# L * L' -> ws.LdagL  (unchanged)
BLAS.gemm!('N', 'C', CT, L_op, L_op, ZT, ws.LdagL)
# Term 2: -0.5 * scalar * L L' * rho
BLAS.gemm!('N', 'N', CT, ws.LdagL, rho, ZT, ws.tmp1)
BLAS.axpy!(T(-0.5 * scalar), ws.tmp1, out)
# Term 3: -0.5 * scalar * rho * L L'
BLAS.gemm!('N', 'N', CT, rho, ws.LdagL, ZT, ws.tmp1)
BLAS.axpy!(T(-0.5 * scalar), ws.tmp1, out)
# Term 1: scalar * conj(L) * rho * L^T  (LdagL now free)
@. ws.tmp2 = conj(L_op)                                    # tmp2 = conj(L)
BLAS.gemm!('N', 'N', CT, ws.tmp2, rho, ZT, ws.tmp1)        # tmp1 = conj(L) * rho
BLAS.gemm!('N', 'T', CT, ws.tmp1, L_op, ZT, ws.LdagL)     # LdagL = conj(L) * rho * L^T
BLAS.axpy!(T(scalar), ws.LdagL, out)
```

**5. `_accumulate_dissipator_2op!` (forward, two-operator BohrDomain):**
Dense convention: `kron(A, B_dag^T) vec(rho) = vec(B_dag^T * rho * A^T)`.
Anticommutator: `B_dag * A` with `kron(B_dag*A, I)` and `kron(I, (B_dag*A)^T)`.
So: sandwich `B_dag^T * rho * A^T`, anticomm uses `(B_dag*A)` with 'N' and `(B_dag*A)^T`.

Currently: sandwich uses `B_dag' * rho * A'` (with 'C' flags). Anticomm uses `(B_dag*A)'`.
Fix: change 'C' to 'T' in sandwich and anticommutator transpose terms.
```
# B_dag * A -> ws.LdagL  (unchanged)
BLAS.gemm!('N', 'N', CT, B_dag, A, ZT, ws.LdagL)
# Term 2: -0.5 * scalar * (B_dag*A)^T * rho  (change 'C' to 'T')
BLAS.gemm!('T', 'N', CT, ws.LdagL, rho, ZT, ws.tmp1)
BLAS.axpy!(T(-0.5 * scalar), ws.tmp1, out)
# Term 3: -0.5 * scalar * rho * (B_dag*A)^T  (change 'C' to 'T')
BLAS.gemm!('N', 'T', CT, rho, ws.LdagL, ZT, ws.tmp1)
BLAS.axpy!(T(-0.5 * scalar), ws.tmp1, out)
# Term 1: scalar * B_dag^T * rho * A^T  (change 'C' to 'T')
BLAS.gemm!('T', 'N', CT, B_dag, rho, ZT, ws.tmp1)
BLAS.gemm!('N', 'T', CT, ws.tmp1, A, ZT, ws.tmp2)
BLAS.axpy!(T(scalar), ws.tmp2, out)
```

**6. `_accumulate_adjoint_dissipator_2op!` (HS adjoint, two-operator BohrDomain):**
HS adjoint of `X rho Y` is `X^H rho Y^H`.
Forward: X = B_dag^T, Y = A^T. X^H = conj(B_dag), Y^H = conj(A).
So adjoint sandwich: `conj(B_dag) * rho * conj(A)`.
Anticomm forward: `(B_dag*A)^T * rho` and `rho * (B_dag*A)^T`.
Adjoint: `conj(B_dag*A) * rho` and `rho * conj(B_dag*A)`.

Currently the function uses all 'N' flags (designed for real operators where conj = identity).
For complex: need to conjugate B_dag*A and the sandwich operators.
```
# B_dag * A -> ws.LdagL
BLAS.gemm!('N', 'N', CT, B_dag, A, ZT, ws.LdagL)
# conj(B_dag*A) -- conjugate LdagL in-place for anticommutator
# But we need the conjugated version for terms 2 and 3.
# Do anticommutator terms first, then sandwich.
# Actually, store conj(LdagL) into tmp2, use that for anticomm:
@. ws.tmp2 = conj(ws.LdagL)                                # tmp2 = conj(B_dag*A)
# Term 2: -0.5 * scalar * conj(B_dag*A) * rho
BLAS.gemm!('N', 'N', CT, ws.tmp2, rho, ZT, ws.tmp1)
BLAS.axpy!(T(-0.5 * scalar), ws.tmp1, out)
# Term 3: -0.5 * scalar * rho * conj(B_dag*A)
BLAS.gemm!('N', 'N', CT, rho, ws.tmp2, ZT, ws.tmp1)
BLAS.axpy!(T(-0.5 * scalar), ws.tmp1, out)
# Term 1: scalar * conj(B_dag) * rho * conj(A)
@. ws.tmp2 = conj(B_dag)                                   # tmp2 = conj(B_dag)
BLAS.gemm!('N', 'N', CT, ws.tmp2, rho, ZT, ws.tmp1)        # tmp1 = conj(B_dag) * rho
@. ws.LdagL = conj(A)                                      # LdagL = conj(A), free to reuse
BLAS.gemm!('N', 'N', CT, ws.tmp1, ws.LdagL, ZT, ws.tmp2)  # tmp2 = conj(B_dag) * rho * conj(A)
BLAS.axpy!(T(scalar), ws.tmp2, out)
```

Update all docstrings to document the corrected convention. The key formula change:
- Forward 1-op: `D(L)(rho) = conj(L) rho L^T - 0.5{L'L, rho}` (matches kron(L, conj(L)))
- Adjoint 1-op: `D*(L)(rho) = L^T rho conj(L) - 0.5{L'L, rho}`
- Forward adj_L: `D(L')(rho) = L^T rho conj(L) - 0.5{LL', rho}` (operator is L^H)
- Adjoint adj_L: `D*(L')(rho) = conj(L) rho L^T - 0.5{LL', rho}`
- Forward 2-op: `D(A,B_dag)(rho) = B_dag^T rho A^T - 0.5{(B_dag A)^T, rho}` (note: anticomm uses transpose)
- Adjoint 2-op: `D*(A,B_dag)(rho) = conj(B_dag) rho conj(A) - 0.5{conj(B_dag A), rho}`

IMPORTANT: The `@.` broadcast conjugation is zero-allocation when writing to a pre-allocated
same-typed matrix. Verify this in the allocation test.
  </action>
  <verify>
Run existing test suite:
```
cd /Users/bence/code/QuantumFurnace.jl && julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -40
```
All 19 existing krylov matvec tests must still pass (Pauli X/Y/Z are Hermitian/real, so
conj(J) rho J^T = J rho J' for them). Zero-allocation tests must still pass.
  </verify>
  <done>
All 6 dissipator functions use the kron-based convention. Existing tests pass unchanged.
  </done>
</task>

<task type="auto">
  <name>Task 2: Add complex jump operator round-trip test</name>
  <files>test/test_krylov_matvec.jl</files>
  <action>
Add a new testset "Round-trip: matvec vs dense with complex non-Hermitian jump (EnergyDomain)"
after the existing testset 19 (BohrDomain duality).

The test must use a jump operator that has both real and imaginary parts AND is NOT Hermitian,
so that `conj(J) rho J^T != J rho J'`. A good choice is a random complex matrix (not
constructed from Paulis).

Test structure:
1. Create a single complex non-Hermitian jump operator for the 4-qubit test system:
   ```julia
   # Create a random complex non-Hermitian jump operator
   rng = Random.MersenneTwister(42)  # deterministic
   raw_jump = (randn(rng, ComplexF64, DIM, DIM)) ./ sqrt(DIM)
   jump_in_eigen = TEST_HAM.eigvecs' * raw_jump * TEST_HAM.eigvecs
   # Not orthogonal, not Hermitian
   complex_jump = JumpOp(raw_jump, jump_in_eigen, false, false)
   complex_jumps = [complex_jump]
   ```

2. Test round-trip for EnergyDomain (forward):
   ```julia
   config = make_liouv_config(EnergyDomain(); with_coherent=true)
   L_dense = construct_lindbladian(complex_jumps, config, TEST_HAM)
   ws = KrylovWorkspace(config, TEST_HAM, complex_jumps)
   for _ in 1:10
       rho = Matrix(random_density_matrix(NUM_QUBITS))
       v_dense = L_dense * vec(rho)
       L_rho = apply_lindbladian!(ws, rho, config, TEST_HAM)
       @test norm(v_dense - vec(L_rho)) < 1e-12
   end
   ```

3. Test round-trip for EnergyDomain (adjoint):
   ```julia
   for _ in 1:10
       rho = Matrix(random_density_matrix(NUM_QUBITS))
       v_adj_dense = L_dense' * vec(rho)
       L_adj_rho = apply_adjoint_lindbladian!(ws, rho, config, TEST_HAM)
       @test norm(v_adj_dense - vec(L_adj_rho)) < 1e-12
   end
   ```

4. Test adjoint duality with complex jump:
   ```julia
   for _ in 1:5
       X = Matrix(random_density_matrix(NUM_QUBITS))
       Y = Matrix(random_density_matrix(NUM_QUBITS))
       L_Y = copy(apply_lindbladian!(ws, Y, config, TEST_HAM))
       lhs = tr(X' * L_Y)
       Lstar_X = copy(apply_adjoint_lindbladian!(ws, X, config, TEST_HAM))
       rhs = tr(Lstar_X' * Y)
       @test abs(lhs - rhs) < 1e-11
   end
   ```

Add `using Random` at the top of the test file (after the existing `using` statements) if not
already present.

Also add a similar test for TimeDomain with the same complex jump to ensure the Time/Trotter
code path is also covered:
   ```julia
   @testset "Round-trip: complex jump (TimeDomain)" begin
       config_td = make_liouv_config(TimeDomain(); with_coherent=true)
       L_dense_td = construct_lindbladian(complex_jumps, config_td, TEST_HAM)
       ws_td = KrylovWorkspace(config_td, TEST_HAM, complex_jumps)
       for _ in 1:10
           rho = Matrix(random_density_matrix(NUM_QUBITS))
           @test norm(L_dense_td * vec(rho) - vec(apply_lindbladian!(ws_td, rho, config_td, TEST_HAM))) < 1e-12
           @test norm(L_dense_td' * vec(rho) - vec(apply_adjoint_lindbladian!(ws_td, rho, config_td, TEST_HAM))) < 1e-12
       end
   end
   ```

NOTE: BohrDomain uses a fundamentally different code path (2-op dissipator with eigenbasis-derived
operators that are real-valued by construction). A BohrDomain test with complex jumps would require
a different test setup. The 1-op EnergyDomain and TimeDomain tests are sufficient to validate
the convention fix.
  </action>
  <verify>
Run full test suite:
```
cd /Users/bence/code/QuantumFurnace.jl && julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -50
```
New complex jump tests must pass. All existing tests must still pass.
  </verify>
  <done>
New tests with a complex non-Hermitian jump operator pass for EnergyDomain and TimeDomain,
confirming the Krylov matvec matches the dense kron-based convention for general complex operators.
All existing tests pass unchanged.
  </done>
</task>

</tasks>

<verification>
1. All 19+ existing Krylov matvec tests pass (real/Hermitian jumps)
2. New complex non-Hermitian jump tests pass (EnergyDomain forward + adjoint + duality, TimeDomain forward + adjoint)
3. Zero-allocation tests still pass for EnergyDomain and TimeDomain
4. Full project test suite passes: `julia --project -e 'using Pkg; Pkg.test()'`
</verification>

<success_criteria>
- Krylov matvec produces identical results to dense kron-based Lindbladian for complex non-Hermitian jump operators
- Adjoint duality tr(X' L(Y)) = tr(L*(X)' Y) holds for complex jumps
- No regressions in existing tests
- Zero-allocation property preserved
</success_criteria>

<output>
After completion, create `.planning/quick/35-fix-krylov-matvec-dissipator-convention-/35-SUMMARY.md`
</output>
