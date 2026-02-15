---
phase: 18-fix-test-oftcaches-constructor-calls
plan: 01
type: execute
wave: 1
depends_on: []
files_modified: [test/test_dm_scaling.jl]
autonomous: true
user_setup: []

must_haves:
  truths:
    - "Test file compiles without errors related to OFTCaches constructor"
    - "Both test cases in DMTST-06 and DMTST-06b pass successfully"
    - "OFTCaches instances are created with explicit Float64 type parameter"
  artifacts:
    - path: "test/test_dm_scaling.jl"
      provides: "Updated OFTCaches constructor calls with type parameters"
      lines: [157, 207]
  key_links:
    - from: "test/test_dm_scaling.jl line 157"
      to: "src/structs.jl OFTCaches inner constructor"
      via: "explicit type parameter Float64"
      pattern: "OFTCaches\\{Float64\\}"
    - from: "test/test_dm_scaling.jl line 207"
      to: "src/structs.jl OFTCaches inner constructor"
      via: "explicit type parameter Float64"
      pattern: "OFTCaches\\{Float64\\}"
---

<objective>
Update two OFTCaches constructor calls in test/test_dm_scaling.jl to use explicit Float64 type parameter.

Purpose: OFTCaches no longer has a convenience constructor without type parameter. The inner constructor requires explicit type specification.

Output: Working test file with corrected constructor calls that compile and pass.
</objective>

<execution_context>
@./.claude/get-shit-done/workflows/execute-plan.md
</execution_context>

<context>
@test/test_dm_scaling.jl
@src/structs.jl
</context>

<tasks>

<task type="auto">
  <name>Update OFTCaches constructor calls to use Float64 type parameter</name>
  <files>test/test_dm_scaling.jl</files>
  <action>
Update two constructor calls in test/test_dm_scaling.jl to match the OFTCaches inner constructor signature which requires explicit type parameter.

Line 157 (in DMTST-06 testset, after oft_time_labels assignment):
  OLD: caches = QuantumFurnace.OFTCaches(DIM)
  NEW: caches = QuantumFurnace.OFTCaches{Float64}(DIM)

Line 207 (in DMTST-06b testset, after time_oft_prefactor assignment):
  OLD: caches = QuantumFurnace.OFTCaches(DIM)
  NEW: caches = QuantumFurnace.OFTCaches{Float64}(DIM)

Reason: src/structs.jl defines only the inner constructor OFTCaches{T}(dim::Int) where {T<:AbstractFloat}, requiring explicit type parameter at call site.
  </action>
  <verify>
Run: julia -e "include(\"test/test_dm_scaling.jl\")" to verify syntax is valid and no constructor errors occur.

Alternately, run: julia --project -e "@testset \"\" include(\"test/test_dm_scaling.jl\") end" if test harness is set up.

Verify both lines 157 and 207 now have OFTCaches{Float64}(DIM) syntax.
  </verify>
  <done>
test/test_dm_scaling.jl compiles without errors related to OFTCaches constructor calls. Both DMTST-06 and DMTST-06b test cases can execute caches instantiation without type-related failures. Lines 157 and 207 both explicitly specify Float64 type parameter.
  </done>
</task>

</tasks>

<verification>
1. File parses without syntax errors
2. OFTCaches constructor calls supply explicit Float64 type parameter
3. Test cases proceed past caches initialization
</verification>

<success_criteria>
- test/test_dm_scaling.jl line 157: `caches = QuantumFurnace.OFTCaches{Float64}(DIM)`
- test/test_dm_scaling.jl line 207: `caches = QuantumFurnace.OFTCaches{Float64}(DIM)`
- No compilation errors from constructor calls
- Tests execute without type-related failures
</success_criteria>

<output>
After completion, update .planning/STATE.md if needed and verify git status shows test/test_dm_scaling.jl as modified.
</output>
