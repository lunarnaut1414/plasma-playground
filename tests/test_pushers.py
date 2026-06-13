"""V1 & V2 — Boris pusher validation.

V1: a charged particle in a uniform B field gyrates at fixed energy on a circle
    of radius r_L. The Boris algorithm should conserve energy to ~machine
    precision (its defining property) and trace the analytic Larmor radius.
V2: in crossed E and B the guiding center drifts at exactly v = E/B, independent
    of charge and mass.
"""

import numpy as np

from plasmaplay import constants as k
from plasmaplay import fields, pushers


def _run_uniform_B(n_orbits=20, steps_per_orbit=60):
    """Proton in a 1 T field along z, launched perpendicular."""
    B0 = 1.0
    v_perp = 1.0e5
    omega_c = k.gyrofrequency(k.e, B0, k.m_p)
    T_c = 2 * np.pi / omega_c
    dt = T_c / steps_per_orbit
    t, pos, vel = pushers.boris_push(
        position=[0, 0, 0], velocity=[v_perp, 0, 0],
        charge=k.e, mass=k.m_p,
        E_func=fields.zero_field(), B_func=fields.uniform_B(B0),
        dt=dt, n_steps=n_orbits * steps_per_orbit,
    )
    return B0, v_perp, t, pos, vel


# --- V1 -------------------------------------------------------------------

def test_v1_energy_conservation():
    _, _, _, _, vel = _run_uniform_B()
    speed2 = np.sum(vel**2, axis=1)
    rel_drift = (speed2.max() - speed2.min()) / speed2.mean()
    assert rel_drift < 1e-10        # Boris is exactly energy-conserving


def test_v1_gyroradius():
    B0, v_perp, _, pos, _ = _run_uniform_B()
    r_L = k.gyroradius(v_perp, k.e, B0, k.m_p)
    # The orbit is a circle; its center is the trajectory centroid (integer
    # number of full orbits), and every point sits one r_L from it.
    center = pos.mean(axis=0)
    radii = np.linalg.norm(pos - center, axis=1)
    assert np.isclose(radii.mean(), r_L, rtol=1e-2)
    assert radii.std() / radii.mean() < 1e-2     # genuinely circular


# --- V2 -------------------------------------------------------------------

def test_v2_ExB_drift():
    B0, Ey = 1.0, 1.0e4
    omega_c = k.gyrofrequency(k.e, B0, k.m_p)
    T_c = 2 * np.pi / omega_c
    steps_per_orbit, n_orbits = 60, 20
    dt = T_c / steps_per_orbit
    _, _, vel = pushers.boris_push(
        position=[0, 0, 0], velocity=[1.0e5, 0, 0],
        charge=k.e, mass=k.m_p,
        E_func=fields.uniform_E(Ey=Ey), B_func=fields.uniform_B(B0),
        dt=dt, n_steps=n_orbits * steps_per_orbit,
    )
    # E = Ey ŷ, B = Bz ẑ  ->  v_drift = (E×B)/B² = (Ey/Bz) x̂.
    # Averaging v_x over an integer number of gyro-periods removes the
    # gyration, leaving the drift.
    v_drift_measured = vel[:, 0].mean()
    v_drift_theory = Ey / B0
    assert np.isclose(v_drift_measured, v_drift_theory, rtol=1e-2)
    # Drift is purely in x: no net y-velocity.
    assert abs(vel[:, 1].mean()) < 1e-2 * v_drift_theory
