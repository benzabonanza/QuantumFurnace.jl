---
status: testing
phase: 35-workspace-and-channel-consolidation
source: 35-01-SUMMARY.md, 35-02-SUMMARY.md
started: 2026-02-27T10:00:00Z
updated: 2026-02-27T10:05:00Z
---

## Current Test
<!-- OVERWRITE each test - shows where we are -->

number: 2
name: Workspace{Krylov} constructs from Config{Lindbladian}
expected: |
  `Workspace(config)` where config is `Config{Lindbladian,...}` returns a `Workspace{Krylov,...}` with KrylovScratch.
awaiting: user response

## Tests

### 1. Full test suite passes
expected: All tests pass (1199 expected) confirming the workspace consolidation introduced no regressions
result: pass

### 2. Workspace{Krylov} constructs from Config{Lindbladian}
expected: `Workspace(config)` where config is a `Config{Lindbladian,...}` returns a `Workspace{Krylov,...}` with KrylovScratch. The old `KrylovWorkspace` alias also works.
result: [pending]

### 3. Workspace{Thermalize} constructs from Config{Thermalize}
expected: `Workspace(config)` where config is a `Config{Thermalize,...}` returns a `Workspace{Krylov,...}` with ThermalizeScratch (DM path uses Krylov workspace with ThermalizeScratch).
result: [pending]

### 4. Workspace{Trajectory} constructs via _build_trajectory_workspace
expected: `QuantumFurnace._build_trajectory_workspace(config, ham, jumps)` returns a `Workspace{Trajectory,...}` with TrajectoryScratch and flat per-operator Kraus vectors (ws.Rs, ws.K0s, ws.U_residuals).
result: [pending]

### 5. Near-zero allocation matvec hot path
expected: Krylov matvec allocation tests pass with allocations <= 512 bytes (MATVEC_ALLOC_BUDGET). Run the krylov_matvec tests specifically.
result: [pending]

### 6. Thread-safe workspace copies
expected: `_copy_workspace_for_thread` creates independent workspace copies that share immutable data but have separate scratch buffers. Workspace independence tests pass.
result: [pending]

### 7. Old type names removed from source
expected: `KrausScratch` struct, `LindbladianWorkspace` struct, `TrajectoryWorkspace` struct, `TrajectoryFramework` struct, and `PerOperatorKraus` struct no longer exist in src/. Only `KrylovWorkspace` remains as a const alias.
result: [pending]

### 8. SC type parameter unnecessary
expected: Workspace should not need a 5th type parameter SC for scratch -- it should be inferable from the other parameters (S, D, C, T)
result: issue
reported: "I dont like that we have now a scratch parameter for the structs, SC, is this really necessary, or can we just infer everything already from previous parameters? I think we can."
severity: major

### 9. Krylov singleton should not exist
expected: No new simulator singleton type -- KrylovSpectrum already exists and should be used instead of introducing a separate Krylov type
result: issue
reported: "I never wanted a new simulator type, but now we have Krylov as well, please get rid of it and consolidate everything with KrylovSpectrum. KrylovSpectrum has to be enough."
severity: major

### 10. Id should be stored in workspace
expected: Identity matrix Id should be a field in the workspace, not passed as argument to _jump_contributions()
result: issue
reported: "You are passing Id as argument in _jump_contributions(), Cant you put it into the workspace thats why its there to use it for these things?"
severity: minor

### 11. _TransitionWrap should be removed
expected: No hacky wrapper types -- use the existing pick_transition logic based on configs and their parametrization
result: issue
reported: "Get rid of the _TransitionWrap, I dont want these hacky solutions in the code; just keep the old logic of how we used pick_transition based on configs and their parametrization."
severity: major

### 12. KrylovWorkspace backward-compat alias should be removed
expected: No backward-compatibility aliases -- this is internal code not used by others yet
result: issue
reported: "Get rid of this line and your tries to make things backward compatible, its only my code right now, not used by others yet. const KrylovWorkspace = Workspace"
severity: minor

## Summary

total: 12
passed: 1
issues: 5
pending: 6
skipped: 0

## Gaps

- truth: "Workspace type parameters should be minimal and inferable"
  status: failed
  reason: "User reported: SC 5th type parameter is unnecessary, should be inferable from S,D,C,T"
  severity: major
  test: 8
  root_cause: ""
  artifacts: []
  missing: []
  debug_session: ""

- truth: "No new simulator singleton types beyond existing ones"
  status: failed
  reason: "User reported: Krylov singleton should not exist, KrylovSpectrum is sufficient"
  severity: major
  test: 9
  root_cause: ""
  artifacts: []
  missing: []
  debug_session: ""

- truth: "Identity matrix should be stored in workspace for reuse"
  status: failed
  reason: "User reported: Id passed as argument to _jump_contributions() instead of stored in workspace"
  severity: minor
  test: 10
  root_cause: ""
  artifacts: []
  missing: []
  debug_session: ""

- truth: "No hacky wrapper types -- use existing dispatch patterns"
  status: failed
  reason: "User reported: _TransitionWrap is hacky, use existing pick_transition logic"
  severity: major
  test: 11
  root_cause: ""
  artifacts: []
  missing: []
  debug_session: ""

- truth: "No backward-compatibility aliases for internal code"
  status: failed
  reason: "User reported: const KrylovWorkspace = Workspace alias is unnecessary"
  severity: minor
  test: 12
  root_cause: ""
  artifacts: []
  missing: []
  debug_session: ""
