"""
    ClassicalQMC

Self-contained Stochastic Series Expansion sampler for the package's 1D
Heisenberg and 2D TFIM models.

Sampling uses physical-frame parameters and reports energy, magnetisation,
correlations, signs, and autocorrelation diagnostics. The implementation uses
diagonal updates with an operator loop for Heisenberg and cluster or local
updates for TFIM.

Odd-size periodic Heisenberg systems are outside the validated sign-free
regime; use open boundaries or a directed-loop implementation there.
"""
module ClassicalQMC

using LinearAlgebra
using Random
using Statistics

export SSEResult, run_sse,
    build_sse_heis_model, build_sse_tfim_model,
    sse_exact_reference, sse_reconstruction_error


include("classical_qmc_common.jl")
include("classical_qmc_heisenberg.jl")
include("classical_qmc_tfim.jl")

end # module ClassicalQMC
