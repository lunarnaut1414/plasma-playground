"""plasmaplay — shared utilities for the plasma-playground experiments.

Small, dependency-light helpers reused across experiments: physical constants,
analytic magnetic-field models, particle pushers, and plotting conveniences.
"""

from . import (
    constants, fields, pushers, plotting, integrators, diagnostics, solvers, pic,
    fvm, propulsion,
)

__all__ = [
    "constants", "fields", "pushers", "plotting", "integrators", "diagnostics",
    "solvers", "pic", "fvm", "propulsion",
]
__version__ = "0.1.0"
