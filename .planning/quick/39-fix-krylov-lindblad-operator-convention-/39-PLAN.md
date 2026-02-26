---
phase: quick-39
plan: 01
type: execute
wave: 1
depends_on: []
files_modified:
  - src/krylov_matvec.jl
  - src/krylov_workspace.jl
  - src/qi_tools.jl
autonomous: true

must_haves:
  truths:
    - "All Lindblad sandwich terms use L*rho*L' convention (not conj(L)*rho*L^T)"
    - "Dense Lindbladian and Krylov matvec remain mutually consistent"
    - "Krylov Lindbladian convention matches Thermalize jump_workers convention"
    - "All existing tests pass (round-trip, adjoint duality, cross-validation)"
  artifacts:
    - path: "src/krylov_matvec.jl"
      provides: "Updated sandwich helpers and apply_lindbladian! functions"
      contains: "L_op.*rho"
    - path: "src/krylov_workspace.jl"
      provides: "Updated R_total accumulation with L'L convention using stored operator"
      contains: "_accumulate_R_total!"
    - path: "src/qi_tools.jl"
      provides: "Updated dense Lindbladian vectorization matching L*rho*L'"
      contains: "_vectorize_liouv_diss_and_add!"
  key_links:
    - from: "src/krylov_matvec.jl"
      to: "src/qi_tools.jl"
      via: "Same operator convention for sandwich terms"
      pattern: "L_op.*rho|jump.*rho"
    - from: "src/krylov_matvec.jl"
      to: "src/jump_workers.jl"
      via: "Matching L*rho*L' convention with Thermalize path"
      pattern: "jump_oft.*rho.*jump_oft'"
---

<objective>
Fix Lindblad jump operator convention across the entire codebase to consistently use
L*rho*L' (where L is the stored operator `in_eigenbasis` / `jump_oft`), matching the
Thermalize path in jump_workers.jl.

Currently, the dense Lindbladian (qi_tools.jl) and Krylov matvec (krylov_matvec.jl)
use `conj(L)*rho*L^T`, while the Thermalize path uses `L*rho*L'`. These differ for
complex non-Hermitian operators. The fix changes the Krylov and dense code to match
the Thermalize convention.

Purpose: Ensure spectral gap estimates from the Krylov Lindbladian path match the
physics of the Thermalize simulation path. Without this fix, the Krylov-computed
fixed point differs from the Thermalize fixed point for complex jump operators.

Output: Updated src/krylov_matvec.jl, src/krylov_workspace.jl, src/qi_tools.jl
</objective>

<execution_context>
@./.claude/get-shit-done/workflows/execute-plan.md
@./.claude/get-shit-done/templates/summary.md
</execution_context>

<context>
@.planning/STATE.md
@src/krylov_matvec.jl
@src/krylov_workspace.jl
@src/krylov_eigsolve.jl
@src/qi_tools.jl
@src/jump_workers.jl
@test/test_krylov_matvec.jl
@test/test_krylov_crossvalidation.jl
</context>

<tasks>

<task type="auto">
  <name>Task 1: Fix dense Lindbladian vectorization convention in qi_tools.jl</name>
  <files>src/qi_tools.jl</files>
  <action>
In `_vectorize_liouv_diss_and_add!` (single-operator version, ~line 64), change from
`conj(J)*rho*J^T` convention to `J*rho*J'` convention:

Current code:
```julia
@. jump_conj = conj(jump)
_kron!(L_target, jump, jump_conj, scalar)        # kron(J, conj(J)) => conj(J)*rho*J^T
mul!(jump_dag_jump, jump', jump)                  # J'J
_kron!(L_target, jump_dag_jump, Id, -0.5 * scalar)
_kron!(L_target, Id, transpose(jump_dag_jump), -0.5 * scalar)
```

Change to:
```julia
@. jump_conj = conj(jump)
_kron!(L_target, jump_conj, jump, scalar)        # kron(conj(J), J) => J*rho*conj(J)^T = J*rho*J'
mul!(jump_dag_jump, jump', jump)                  # J'J (same -- Hermitian, convention-independent)
_kron!(L_target, Id, jump_dag_jump, -0.5 * scalar)   # kron(I, J'J) => J'J*rho*I = J'J*rho
_kron!(L_target, transpose(jump_dag_jump), Id, -0.5 * scalar)  # kron((J'J)^T, I) => I*rho*(J'J) = rho*J'J
```

Explanation of the new kron terms:
- `kron(conj(J), J) * vec(rho) = vec(J * rho * conj(J)^T) = vec(J * rho * J')` since `conj(J)^T = J'`. This is the desired `L*rho*L'`.
- `kron(I, J'J) * vec(rho) = vec(J'J * rho)`. This is `-0.5 * J'J * rho`.
- `kron((J'J)^T, I) * vec(rho) = vec(rho * ((J'J)^T)^T) = vec(rho * J'J)`. This is `-0.5 * rho * J'J`.
- Together: `J*rho*J' - 0.5{J'J, rho}` -- the standard Lindblad dissipator.

For the two-operator version (`_vectorize_liouv_diss_and_add!` with jump_1, jump_2, ~line 86):

Current code:
```julia
_kron!(L_target, jump_1, transpose(jump_2), scalar)       # kron(J1, J2^T) => J2^T*rho*J1^T
mul!(jump2_jump1, jump_2, jump_1)                          # J2*J1
_kron!(L_target, jump2_jump1, Id, -0.5 * scalar)
_kron!(L_target, Id, transpose(jump2_jump1), -0.5 * scalar)
```

Change to:
```julia
_kron!(L_target, transpose(jump_2), jump_1, scalar)       # kron(J2^T, J1) => J1*rho*(J2^T)^T = J1*rho*J2
mul!(jump2_jump1, jump_2, jump_1)                          # J2*J1
_kron!(L_target, Id, jump2_jump1, -0.5 * scalar)           # J2*J1*rho
_kron!(L_target, transpose(jump2_jump1), Id, -0.5 * scalar)  # rho*J2*J1
```

Note: The two-operator docstring says `L = J1 * X * J2 - 0.5 * (J2 * J1 * X + X * J2 * J1)`.
With the new kron ordering, `kron(J2^T, J1)*vec(rho) = vec(J1*rho*(J2^T)^T) = vec(J1*rho*J2)`.
This gives the dissipator `J1*rho*J2 - 0.5{J2*J1, rho}` which matches the docstring.

Also update the docstring on the single-operator version to clarify the convention is `J*rho*J' - 0.5{J'J, rho}`.
  </action>
  <verify>
Run the Krylov matvec round-trip tests to confirm the dense Lindbladian and Krylov matvec
(which will be updated in Task 2) are still mutually consistent:
```bash
cd /Users/bence/code/QuantumFurnace.jl && julia --project -e 'using Pkg; Pkg.test(test_args=["krylov_matvec"])'
```
Note: tests may fail until Task 2 is complete (Krylov side must also be updated).
  </verify>
  <done>
Dense Lindbladian vectorization uses J*rho*J' - 0.5{J'J, rho} convention via
kron(conj(J), J) for sandwich and kron(I, J'J)/kron((J'J)^T, I) for anticommutator.
  </done>
</task>

<task type="auto">
  <name>Task 2: Fix Krylov sandwich helpers and R_total accumulation</name>
  <files>src/krylov_matvec.jl, src/krylov_workspace.jl</files>
  <action>
**Part A: Fix sandwich helpers in krylov_matvec.jl**

All four single-operator sandwich functions need swapping. The pattern: replace
`conj(L)*rho*L^T` with `L*rho*L'`, and the adjoint versions accordingly.

1. `_accumulate_sandwich!` (~line 38): Forward sandwich for positive frequency.

Current (conj(L)*rho*L^T):
```julia
@. ws.tmp2 = conj(L_op)
BLAS.gemm!('N', 'N', CT, ws.tmp2, rho, ZT, ws.tmp1)       # conj(L)*rho
BLAS.gemm!('N', 'T', CT, ws.tmp1, L_op, ZT, ws.LdagL)     # conj(L)*rho*L^T
```

Change to (L*rho*L'):
```julia
BLAS.gemm!('N', 'N', CT, L_op, rho, ZT, ws.tmp1)           # L*rho
BLAS.gemm!('N', 'C', CT, ws.tmp1, L_op, ZT, ws.LdagL)      # L*rho*L'
```
No need for ws.tmp2 anymore. Update docstring to say `scalar * L * rho * L'`.

2. `_accumulate_sandwich_adj_L!` (~line 63): Forward sandwich for negative-frequency
   Hermitian partner (M = L').

Current (L^T*rho*conj(L)):
```julia
@. ws.tmp2 = conj(L_op)
BLAS.gemm!('T', 'N', CT, L_op, rho, ZT, ws.tmp1)           # L^T*rho
BLAS.gemm!('N', 'N', CT, ws.tmp1, ws.tmp2, ZT, ws.LdagL)   # L^T*rho*conj(L)
```

Change to (L'*rho*L):
```julia
BLAS.gemm!('C', 'N', CT, L_op, rho, ZT, ws.tmp1)           # L'*rho
BLAS.gemm!('N', 'N', CT, ws.tmp1, L_op, ZT, ws.LdagL)      # L'*rho*L
```
No need for ws.tmp2. Update docstring to say `scalar * L' * rho * L`.

3. `_accumulate_adjoint_sandwich!` (~line 89): HS adjoint of forward sandwich.

The HS adjoint of `L*rho*L'` is `L'*rho*L`. This is the same computation as the
new `_accumulate_sandwich_adj_L!`.

Current (L^T*rho*conj(L)):
```julia
@. ws.tmp2 = conj(L_op)
BLAS.gemm!('T', 'N', CT, L_op, rho, ZT, ws.tmp1)
BLAS.gemm!('N', 'N', CT, ws.tmp1, ws.tmp2, ZT, ws.LdagL)
```

Change to (L'*rho*L):
```julia
BLAS.gemm!('C', 'N', CT, L_op, rho, ZT, ws.tmp1)           # L'*rho
BLAS.gemm!('N', 'N', CT, ws.tmp1, L_op, ZT, ws.LdagL)      # L'*rho*L
```
Update docstring to say `scalar * L' * rho * L` (HS adjoint of `L*rho*L'`).

4. `_accumulate_adjoint_sandwich_adj_L!` (~line 115): HS adjoint of negative-freq sandwich.

The HS adjoint of `L'*rho*L` is `L*rho*L'`. This is the same computation as the
new `_accumulate_sandwich!`.

Current (conj(L)*rho*L^T):
```julia
@. ws.tmp2 = conj(L_op)
BLAS.gemm!('N', 'N', CT, ws.tmp2, rho, ZT, ws.tmp1)
BLAS.gemm!('N', 'T', CT, ws.tmp1, L_op, ZT, ws.LdagL)
```

Change to (L*rho*L'):
```julia
BLAS.gemm!('N', 'N', CT, L_op, rho, ZT, ws.tmp1)           # L*rho
BLAS.gemm!('N', 'C', CT, ws.tmp1, L_op, ZT, ws.LdagL)      # L*rho*L'
```
Update docstring to say `scalar * L * rho * L'` (HS adjoint of `L'*rho*L`).

5. `_accumulate_sandwich_2op!` (~line 279): BohrDomain forward 2-operator sandwich.

Current (B_dag^T*rho*A^T):
```julia
BLAS.gemm!('T', 'N', CT, B_dag, rho, ZT, ws.tmp1)      # B_dag^T*rho
BLAS.gemm!('N', 'T', CT, ws.tmp1, A, ZT, ws.tmp2)       # B_dag^T*rho*A^T
```

The Bohr domain dissipator should be `A*rho*B_dag` (matching the docstring convention
`J1*X*J2` from qi_tools.jl). Actually, looking at the BohrDomain `_vectorize_liouv_diss_and_add!`
call: `_vectorize_liouv_diss_and_add!(L_target, alpha_A_nu1, A_nu_2_dag, gamma_norm_factor, ws)`.
With the NEW qi_tools convention, this computes `alpha_A*rho*A_nu2_dag` (J1*rho*J2).

So the 2op sandwich should change to `A*rho*B_dag` (where A=alpha_A=jump_1, B_dag=A_nu2_dag=jump_2).

Change to:
```julia
BLAS.gemm!('N', 'N', CT, A, rho, ZT, ws.tmp1)           # A*rho
BLAS.gemm!('N', 'N', CT, ws.tmp1, B_dag, ZT, ws.tmp2)   # A*rho*B_dag
```
Update docstring.

6. `_accumulate_adjoint_sandwich_2op!` (~line 304): HS adjoint of 2op sandwich.

The HS adjoint of `A*rho*B_dag` is `B_dag'*rho*A'`.

Change to:
```julia
BLAS.gemm!('C', 'N', CT, B_dag, rho, ZT, ws.tmp1)       # B_dag'*rho
BLAS.gemm!('N', 'C', CT, ws.tmp1, A, ZT, ws.tmp2)       # B_dag'*rho*A'
```
No need for ws.tmp2 = conj(...) scratch. Update docstring.

**Part B: Fix R_total accumulation in krylov_workspace.jl**

The R_total comment says "physics convention: R = sum rate^2 * L' * L". Since L is the
stored operator (with the new convention), L'L is `jump_oft' * jump_oft`, which is what
the code already does:

```julia
mul!(LdagL, jump_oft', jump_oft)   # L'L
```

This is `J'J` where J = stored operator. With the new convention where L = J (the stored
operator), this correctly computes `L'L`. So the R_total computation does NOT need changing.

The negative-frequency partner has L_neg = L', so L_neg'*L_neg = L*L'. The code does:
```julia
mul!(LdagL, jump_oft, jump_oft')   # L*L' for negative freq
```
This is also already correct.

Similarly for the BohrDomain R_total (`mul!(R, A_nu2_dag, jump_oft, ...)`), this needs
to match the new convention. Currently it computes `A_nu2_dag * alpha_A` which matches
the anticommutator `J2*J1` in the 2-op dissipator. Since we changed the 2-op kron terms
to use `kron(I, J2*J1)` and `kron((J2*J1)^T, I)`, this R_total is already correct.

**Verify R_total consistency**: The R_total should equal `sum rate^2 * J'J` for the
single-operator case. The existing code does `mul!(LdagL, jump_oft', jump_oft)` which
computes `J'J`. This is correct and does not change.

**Part C: Update comments and docstrings in krylov_matvec.jl**

Update the inline comments in `apply_lindbladian!` and `apply_adjoint_lindbladian!`
(all domain variants) to say `L * rho * L'` instead of `conj(L) * rho * L^T`.
Specifically:
- Line 169: "Sandwich-only loop: sum_i scalar_i * L_i * rho * L_i'"
- Line 237: "Adjoint sandwich-only loop: sum_i scalar_i * L_i' * rho * L_i"
- Line 481: Same pattern for TimeDomain/TrotterDomain
- Line 551: Same for adjoint TimeDomain/TrotterDomain
  </action>
  <verify>
Run the full Krylov test suite (matvec round-trips, cross-validation, eigsolve):
```bash
cd /Users/bence/code/QuantumFurnace.jl && julia --project -e 'using Pkg; Pkg.test(test_args=["krylov_matvec", "krylov_crossvalidation", "krylov_eigsolve"])'
```
All tests must pass. The round-trip tests compare Krylov matvec against the dense
Lindbladian (updated in Task 1), so both sides now use `L*rho*L'` and should agree.
  </verify>
  <done>
All 4 single-operator sandwich helpers, both 2-operator sandwich helpers, and all
comments/docstrings updated to L*rho*L' convention. R_total accumulation verified
unchanged (already correct). All Krylov tests pass.
  </done>
</task>

<task type="auto">
  <name>Task 3: Run full test suite and verify convention consistency</name>
  <files></files>
  <action>
Run the complete test suite to ensure no regressions:
```bash
cd /Users/bence/code/QuantumFurnace.jl && julia --project -e 'using Pkg; Pkg.test()'
```

If any tests fail, diagnose whether the failure is:
1. Convention mismatch (some code path still uses old convention) -- fix it
2. Numerical tolerance issue from changed convention -- adjust threshold
3. Pre-existing failure unrelated to this change -- document and skip

Key tests to watch:
- test_krylov_matvec.jl: Round-trip forward/adjoint for all domains (Energy, Time,
  Trotter, Bohr), complex non-Hermitian jump tests (testsets 20-23), adjoint duality,
  zero allocations
- test_krylov_crossvalidation.jl: Krylov vs dense gap matching, L-vs-E convergence
- test_krylov_eigsolve.jl: Eigsolve results
- Any thermalization tests that cross-reference Lindbladian results
  </action>
  <verify>
Full test suite passes:
```bash
cd /Users/bence/code/QuantumFurnace.jl && julia --project -e 'using Pkg; Pkg.test()'
```
  </verify>
  <done>
All tests pass. Dense Lindbladian, Krylov matvec, and Thermalize path all
consistently use L*rho*L' convention for Lindblad jump operators.
  </done>
</task>

</tasks>

<verification>
1. Round-trip tests pass: `apply_lindbladian!(ws, rho, config, ham)` matches
   `construct_lindbladian(jumps, config, ham) * vec(rho)` for all domains
2. Adjoint duality holds: `tr(X' * L(Y)) == tr(L*(X)' * Y)` within 1e-11
3. Cross-validation passes: Krylov spectral gap matches dense eigen() within 1e-8
4. L-vs-E convergence: Channel-derived gap converges to Lindbladian gap with order >= 0.9
5. Complex non-Hermitian jump tests pass (testsets 20-23) -- these are the most
   sensitive to convention differences
6. Zero-allocation property maintained for Krylov matvec hot path
</verification>

<success_criteria>
- All Krylov matvec round-trip tests pass
- All cross-validation tests pass
- All existing tests pass (full test suite)
- Convention is consistently L*rho*L' across dense Lindbladian, Krylov matvec,
  and Thermalize paths
</success_criteria>

<output>
After completion, create `.planning/quick/39-fix-krylov-lindblad-operator-convention-/39-SUMMARY.md`
</output>
