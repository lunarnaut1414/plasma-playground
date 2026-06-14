"""Unit tests for plasmaplay.fields.

The uniform/zero factories are simple but underpin every pusher test, so we pin
their contract (shape, position-independence). The magnetic mirror gets the most
attention: its on-axis profile, its symmetry, the inward-pointing radial field,
and — most important physically — that it satisfies div(B) = 0, which is the
whole reason its functional form looks the way it does.
"""

import numpy as np

from plasmaplay import fields
from plasmaplay.constants import MU_0

SAMPLE_POINTS = [
    np.array([0.0, 0.0, 0.0]),
    np.array([0.1, -0.2, 0.3]),
    np.array([-1.0, 2.0, -0.5]),
]


def _numerical_divergence(field, p, h=1e-5):
    """Central-difference div(B) at point p."""
    div = 0.0
    for i in range(3):
        pp = p.astype(float).copy(); pp[i] += h
        pm = p.astype(float).copy(); pm[i] -= h
        div += (field(pp)[i] - field(pm)[i]) / (2 * h)
    return div


# --- uniform / zero -------------------------------------------------------

def test_uniform_B_is_constant_and_well_shaped():
    B = fields.uniform_B(2.5)
    for p in SAMPLE_POINTS:
        out = np.asarray(B(p))
        assert out.shape == (3,)
        np.testing.assert_allclose(out, [0.0, 0.0, 2.5])


def test_uniform_E_is_constant_and_well_shaped():
    E = fields.uniform_E(Ex=1.0, Ey=-2.0, Ez=3.0)
    for p in SAMPLE_POINTS:
        np.testing.assert_allclose(np.asarray(E(p)), [1.0, -2.0, 3.0])


def test_zero_field():
    Z = fields.zero_field()
    for p in SAMPLE_POINTS:
        np.testing.assert_allclose(np.asarray(Z(p)), [0.0, 0.0, 0.0])


# --- magnetic mirror ------------------------------------------------------

def test_mirror_on_axis_center_value():
    B0 = 1.3
    mirror = fields.magnetic_mirror(B0=B0, mirror_ratio=4.0, length=2.0)
    np.testing.assert_allclose(mirror([0.0, 0.0, 0.0]), [0.0, 0.0, B0])


def test_mirror_throat_value_equals_ratio_times_B0():
    B0, ratio, length = 1.0, 4.0, 2.0
    mirror = fields.magnetic_mirror(B0=B0, mirror_ratio=ratio, length=length)
    # at z = length/2 (the throat) the on-axis field should be B0 * mirror_ratio
    Bz_throat = mirror([0.0, 0.0, length / 2])[2]
    assert np.isclose(Bz_throat, B0 * ratio, rtol=1e-12)


def test_mirror_field_is_minimum_at_center():
    mirror = fields.magnetic_mirror(B0=1.0, mirror_ratio=3.0, length=2.0)
    Bz0 = mirror([0.0, 0.0, 0.0])[2]
    for z in (0.1, 0.4, 0.8, 1.0):
        assert mirror([0.0, 0.0, z])[2] > Bz0


def test_mirror_axial_symmetry():
    mirror = fields.magnetic_mirror(B0=1.0, mirror_ratio=3.0, length=2.0)
    for z in (0.2, 0.5, 0.9):
        assert np.isclose(mirror([0, 0, z])[2], mirror([0, 0, -z])[2])


def test_mirror_radial_field_points_inward_where_field_rises():
    # for z > 0 the field is increasing, so B_r must point toward the axis
    # (restoring) — i.e. opposite in sign to the radial position.
    mirror = fields.magnetic_mirror(B0=1.0, mirror_ratio=3.0, length=2.0)
    Bx = mirror([0.05, 0.0, 0.3])[0]
    assert Bx < 0          # x > 0 -> B_x < 0
    By = mirror([0.0, 0.05, 0.3])[1]
    assert By < 0          # y > 0 -> B_y < 0


def test_mirror_is_divergence_free():
    # div(B) = 0 is the physical constraint the radial term was chosen to satisfy.
    mirror = fields.magnetic_mirror(B0=1.0, mirror_ratio=4.0, length=2.0)
    for p in ([0.02, 0.01, 0.1], [0.05, -0.03, 0.4], [-0.04, 0.02, -0.6]):
        div = _numerical_divergence(mirror, np.array(p))
        assert abs(div) < 1e-6


# --- V11: Biot-Savart of a circular loop vs the analytic on-axis field ----

def test_v11_biot_savart_loop_on_axis():
    # On the axis of a loop of radius a carrying current I:
    #   B_z(z) = mu0 I a^2 / (2 (a^2 + z^2)^{3/2}),  B_x = B_y = 0.
    a, I = 0.5, 1000.0
    loop = fields.circular_loop(a, I, n_segments=400)
    for z in (0.0, 0.1, 0.3, 0.7, 1.5):
        B = loop([0.0, 0.0, z])
        Bz_analytic = MU_0 * I * a**2 / (2.0 * (a**2 + z**2) ** 1.5)
        assert np.isclose(B[2], Bz_analytic, rtol=1e-3)   # < 0.1% on axis
        assert abs(B[0]) < 1e-6 * abs(Bz_analytic) + 1e-12
        assert abs(B[1]) < 1e-6 * abs(Bz_analytic) + 1e-12


def test_v11_loop_center_value_and_sign():
    # At the center, B_z = mu0 I / (2a), and CCW current gives +z field.
    a, I = 0.5, 1000.0
    loop = fields.circular_loop(a, I, n_segments=400)
    Bz_center = loop([0.0, 0.0, 0.0])[2]
    assert np.isclose(Bz_center, MU_0 * I / (2.0 * a), rtol=1e-3)
    assert Bz_center > 0


# --- helical stellarator: current-free rotational transform (Track E) -------
def _curl(B, p, h=1e-5):
    """Numerical curl of a field B at point p (central differences)."""
    p = np.asarray(p, float)
    d = np.eye(3) * h
    j = [(np.asarray(B(p + d[k])) - np.asarray(B(p - d[k]))) / (2 * h) for k in range(3)]
    # curl = (dBz/dy - dBy/dz, dBx/dz - dBz/dx, dBy/dx - dBx/dy)
    return np.array([j[1][2] - j[2][1], j[2][0] - j[0][2], j[0][1] - j[1][0]])


def test_stellarator_is_current_free():
    """The vacuum stellarator field is curl-free (no current in the plasma region),
    and the loop integral of B around the axis vanishes (no net plasma current) — the
    defining contrast with a current-driven tokamak / screw pinch."""
    B = fields.helical_stellarator(eps=0.4, l=2, h=1.0)
    for p in ([0.4, 0.0, 0.0], [0.2, 0.3, 0.5], [-0.3, 0.25, 1.0]):
        assert np.max(np.abs(_curl(B, p))) < 1e-4          # curl B ~ 0 (current-free)
    th = np.linspace(0, 2 * np.pi, 240, endpoint=False)
    r = 0.4
    circ = sum(np.dot(B([r * np.cos(t), r * np.sin(t), 0.0])[:2],
                      np.array([-r * np.sin(t), r * np.cos(t)]) * (2 * np.pi / 240))
               for t in th)
    assert abs(circ) < 1e-9                                  # no net enclosed current


def test_stellarator_has_transform_growing_with_shaping():
    """A current-free stellarator still twists its field lines (iota != 0), and the
    transform GROWS with the helical shaping amplitude (the 2nd-order geometric origin,
    iota ~ eps^2) — twist from geometry, not from current."""
    from plasmaplay.diagnostics import poincare_section, rotational_transform
    L = 2 * np.pi

    def iota(eps):
        pts = poincare_section(fields.helical_stellarator(eps=eps), x0=[0.5, 0.0, 0.0],
                               period=L, n_crossings=30, ds=0.012)
        return abs(rotational_transform(pts))
    assert iota(0.5) > 5e-3                                  # nonzero transform
    assert iota(0.6) > iota(0.3)                             # grows with shaping
