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


# --- unit-level behavior of boris_push ------------------------------------

def _push(position, velocity, charge, mass, E, B, dt, n):
    return pushers.boris_push(position, velocity, charge, mass, E, B, dt, n)


def test_output_shapes_and_initial_conditions():
    n = 50
    t, pos, vel = _push([1, 2, 3], [4, 5, 6], k.e, k.m_p,
                        fields.zero_field(), fields.uniform_B(1.0), 1e-9, n)
    assert t.shape == (n + 1,)
    assert pos.shape == (n + 1, 3)
    assert vel.shape == (n + 1, 3)
    # the initial state is stored verbatim at index 0
    np.testing.assert_allclose(pos[0], [1, 2, 3])
    np.testing.assert_allclose(vel[0], [4, 5, 6])
    np.testing.assert_allclose(t[:3], [0, 1e-9, 2e-9])


def test_zero_fields_is_straight_line():
    # no E, no B  ->  constant velocity, exactly.
    v0 = np.array([1.0e4, -2.0e4, 3.0e4])
    dt, n = 1e-7, 100
    _, pos, vel = _push([0, 0, 0], v0, k.e, k.m_p,
                        fields.zero_field(), fields.zero_field(), dt, n)
    for i in range(n + 1):
        np.testing.assert_allclose(vel[i], v0, rtol=1e-12)
        np.testing.assert_allclose(pos[i], v0 * dt * i, rtol=1e-12, atol=1e-9)


def test_pure_E_field_uniform_acceleration():
    # B = 0, uniform E  ->  velocity grows linearly: v(t) = v0 + (q/m) E t.
    Ex = 500.0
    dt, n = 1e-9, 200
    t, _, vel = _push([0, 0, 0], [0, 0, 0], k.e, k.m_p,
                      fields.uniform_E(Ex=Ex), fields.zero_field(), dt, n)
    a = k.e * Ex / k.m_p
    np.testing.assert_allclose(vel[:, 0], a * t, rtol=1e-12)
    np.testing.assert_allclose(vel[:, 1:], 0.0, atol=1e-20)


def test_charge_sign_reverses_gyration():
    # positive and negative charges gyrate in opposite senses; check the first
    # step's transverse velocity flips sign.
    args = dict(E=fields.zero_field(), B=fields.uniform_B(1.0), dt=1e-10, n=1)
    _, _, vel_pos = _push([0, 0, 0], [1e5, 0, 0], k.e, k.m_p, **args)
    _, _, vel_neg = _push([0, 0, 0], [1e5, 0, 0], -k.e, k.m_p, **args)
    assert vel_pos[1, 1] != 0
    assert np.sign(vel_pos[1, 1]) == -np.sign(vel_neg[1, 1])


def test_returns_to_start_after_one_period():
    # a full gyro-period should bring the particle back to its start (Boris has
    # only a tiny frequency error, vanishing as steps/orbit increases).
    B0, v_perp = 1.0, 1.0e5
    omega_c = k.gyrofrequency(k.e, B0, k.m_p)
    steps = 1000
    dt = (2 * np.pi / omega_c) / steps
    _, pos, _ = _push([0, 0, 0], [v_perp, 0, 0], k.e, k.m_p,
                      fields.zero_field(), fields.uniform_B(B0), dt, steps)
    r_L = k.gyroradius(v_perp, k.e, B0, k.m_p)
    assert np.linalg.norm(pos[-1] - pos[0]) < 1e-3 * r_L


def test_list_and_array_inputs_agree():
    args = (k.e, k.m_p, fields.zero_field(), fields.uniform_B(1.0), 1e-9, 30)
    _, p_list, v_list = _push([0, 0, 0], [1e5, 0, 0], *args)
    _, p_arr, v_arr = _push(np.zeros(3), np.array([1e5, 0.0, 0.0]), *args)
    np.testing.assert_allclose(p_list, p_arr)
    np.testing.assert_allclose(v_list, v_arr)
