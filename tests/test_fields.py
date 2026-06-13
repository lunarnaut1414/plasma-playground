"""Unit tests for plasmaplay.fields.

The uniform/zero factories are simple but underpin every pusher test, so we pin
their contract (shape, position-independence). The magnetic mirror gets the most
attention: its on-axis profile, its symmetry, the inward-pointing radial field,
and — most important physically — that it satisfies div(B) = 0, which is the
whole reason its functional form looks the way it does.
"""

import numpy as np

from plasmaplay import fields

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
