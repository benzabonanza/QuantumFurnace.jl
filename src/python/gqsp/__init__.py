"""GQSP (Generalized Quantum Signal Processing) tools.

Implements the analytical pipeline (Jacobi-Anger Laurent polynomial,
Berntson-Sünderhauf complementary polynomial, Motlagh-Wiebe angle
extraction) and a Qiskit POC circuit realising e^{-i delta H} from a
block encoding of H.

Conventions follow `1_preliminaries.tex` (cos: cos theta = lambda/alpha,
target e^{-i delta alpha cos theta}).

This is a Python port of the Julia scripts/scratch_gqsp_*.jl pipeline.
"""

from .jacobi_anger import jacobi_anger_coeffs
from .berntson_sunderhauf import complementary_polynomial_bs
from .motlagh_wiebe import gqsp_extract_angles
from .circuit import (
    build_block_encoding_one_anc,
    build_walk,
    gqsp_circuit,
)

__all__ = [
    "jacobi_anger_coeffs",
    "complementary_polynomial_bs",
    "gqsp_extract_angles",
    "build_block_encoding_one_anc",
    "build_walk",
    "gqsp_circuit",
]
