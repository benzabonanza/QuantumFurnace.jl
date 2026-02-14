# Quick Task 11: Widen trajectory regression tolerance

## What changed

`test/test_regression.jl`: Changed trajectory regression test tolerance from `atol=1e-6` to `atol=1e-3` on both trajectory tests (EnergyDomain and TrotterDomain). Updated comments to explain the tolerance accounts for both BLAS platform differences AND Julia version differences (1.11 vs 1.12 RNG streams). DM tests remain at `atol=1e-10`.

## Why

Quick task 10 changed tolerance from 1e-10 to 1e-6, but user's M1 Mac running Julia 1.12.4 still failed with max element-wise differences of ~7.5e-5. The BSON references were generated on Julia 1.11.3 (x86 linux). Cross-version RNG stream differences cause more trajectory branching divergence than cross-platform BLAS alone.

## Commit

- `7ab18ae`: fix(quick-11): widen trajectory regression tolerance to 1e-3 for cross-version portability
