"""V13 — field-line tracing & rotational transform.

In the screw-pinch field (fields.screw_pinch) field lines lie on cylinders
r = const and wind rigidly, so the rotational transform has the closed form
ι = twist · L / (2π Bz). We trace a line, build a Poincaré section at planes
spaced one period L apart, and check (a) the measured ι matches theory and
(b) the surface "closes" — every puncture sits at the same radius.
"""

import numpy as np

from plasmaplay import fields
from plasmaplay.diagnostics import (
    poincare_section,
    rotational_transform,
    trace_field_line,
)

BZ = 1.0
TWIST = 0.2
L = 2 * np.pi                       # axial period
IOTA_THEORY = TWIST * L / (2 * np.pi * BZ)   # = 0.2


def test_v13_rotational_transform_matches_theory():
    B = fields.screw_pinch(Bz=BZ, twist=TWIST)
    crossings = poincare_section(B, x0=[0.3, 0.0, 0.0], period=L,
                                 n_crossings=30, ds=0.01)
    iota = rotational_transform(crossings)
    assert np.isclose(iota, IOTA_THEORY, rtol=1e-2)


def test_v13_flux_surface_closes():
    # Punctures of a single field line should all land on one radius (a clean
    # nested flux surface), and the start radius is preserved.
    B = fields.screw_pinch(Bz=BZ, twist=TWIST)
    r0 = 0.3
    crossings = poincare_section(B, x0=[r0, 0.0, 0.0], period=L,
                                 n_crossings=30, ds=0.01)
    radii = np.hypot(crossings[:, 0], crossings[:, 1])
    assert np.allclose(radii, r0, rtol=1e-3)
    assert radii.std() / radii.mean() < 1e-3


def test_v13_iota_independent_of_radius_for_linear_profile():
    # Linear B_theta = twist*r is shearless: ι is the same on every surface.
    B = fields.screw_pinch(Bz=BZ, twist=TWIST)
    iotas = []
    for r0 in (0.2, 0.4, 0.6):
        c = poincare_section(B, x0=[r0, 0.0, 0.0], period=L,
                             n_crossings=20, ds=0.01)
        iotas.append(rotational_transform(c))
    assert np.allclose(iotas, IOTA_THEORY, rtol=1e-2)


def test_trace_field_line_stays_on_cylinder():
    # Arc-length tracing should keep r constant along the whole line.
    B = fields.screw_pinch(Bz=BZ, twist=TWIST)
    pts = trace_field_line(B, x0=[0.3, 0.0, 0.0], ds=0.02, n_steps=2000)
    r = np.hypot(pts[:, 0], pts[:, 1])
    assert np.allclose(r, 0.3, atol=1e-3)
