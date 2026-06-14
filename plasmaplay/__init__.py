"""plasmaplay — shared utilities for the plasma-playground experiments.

Small, dependency-light helpers reused across experiments: physical constants,
analytic magnetic-field models, particle pushers, and plotting conveniences.
"""

from . import (
    constants, fields, pushers, plotting, integrators, diagnostics, solvers, pic,
    fvm, propulsion, guiding_center, dispersion, tokamak, transport, animate,
    equilibrium_metrics, operating_limits, tearing, cylinder_mhd, reduced_mhd,
)

__all__ = [
    "constants", "fields", "pushers", "plotting", "integrators", "diagnostics",
    "solvers", "pic", "fvm", "propulsion", "guiding_center", "dispersion", "tokamak",
    "transport", "animate", "equilibrium_metrics", "operating_limits", "tearing",
    "cylinder_mhd", "reduced_mhd",
]
__version__ = "0.1.0"
