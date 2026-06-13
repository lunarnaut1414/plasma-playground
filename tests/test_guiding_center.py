"""Validation of the guiding-center integrator (experiment 02).

Three checks, matching the analytic drift catalogue:
  - in crossed E and B the guiding center moves at exactly v_E = E/B,
  - in a grad-B field it drifts at the analytic grad-B velocity, and
  - the guiding-center trajectory matches the gyro-averaged full Boris orbit
    (the whole point of the approximation) in the adiabatic regime.
"""

import numpy as np

from plasmaplay import constants as k
from plasmaplay import fields, guiding_center as gc, pushers


def test_ExB_drift_exact():
    B0, Ey = 1.0, 1.0e3
    # E×B drift = (E×B)/B² = (Ey/Bz) x̂, independent of charge/mass/energy.
    t, pos, _ = gc.gc_push(
        position=[0, 0, 0], v_par=0.0, mu=0.0,
        charge=k.e, mass=k.m_p,
        E_func=fields.uniform_E(Ey=Ey), B_func=fields.uniform_B(B0),
        dt=1e-7, n_steps=200,
    )
    v_meas = (pos[-1, 0] - pos[0, 0]) / t[-1]
    assert np.isclose(v_meas, Ey / B0, rtol=1e-4)
    assert abs(pos[-1, 1]) < 1e-6                  # no y-motion


def test_gradB_drift_matches_analytic():
    B0, L = 1.0, 0.5
    v_perp = 1.0e5
    Bmag0 = B0                                     # |B| at x=0
    mu = gc.magnetic_moment(v_perp, k.m_p, Bmag0)
    t, pos, _ = gc.gc_push(
        position=[0, 0, 0], v_par=0.0, mu=mu,
        charge=k.e, mass=k.m_p,
        E_func=fields.zero_field(), B_func=fields.gradient_B(B0, L),
        dt=1e-7, n_steps=300,
    )
    # analytic grad-B drift speed in +y: v = m v⊥² / (2 q B0 L)
    v_analytic = k.m_p * v_perp**2 / (2 * k.e * B0 * L)
    v_meas = (pos[-1, 1] - pos[0, 1]) / t[-1]
    assert np.isclose(v_meas, v_analytic, rtol=2e-3)
    assert abs(pos[-1, 0]) < 1e-3 * abs(pos[-1, 1])   # drift is in y


def test_guiding_center_matches_gyroaveraged_orbit():
    # Run the full Boris orbit and the guiding center in the same grad-B field;
    # their drift (net displacement) should agree in the adiabatic regime.
    B0, L = 1.0, 2.0
    v_perp = 5.0e4
    charge, mass = k.e, k.m_p
    omega_c = k.gyrofrequency(charge, B0, mass)
    T_c = 2 * np.pi / omega_c
    dt = T_c / 80
    n_orbits = 40
    n = n_orbits * 80

    _, pos_full, _ = pushers.boris_push(
        position=[0, 0, 0], velocity=[v_perp, 0, 0],
        charge=charge, mass=mass,
        E_func=fields.zero_field(), B_func=fields.gradient_B(B0, L),
        dt=dt, n_steps=n,
    )
    mu = gc.magnetic_moment(v_perp, mass, B0)
    _, pos_gc, _ = gc.gc_push(
        position=[0, 0, 0], v_par=0.0, mu=mu,
        charge=charge, mass=mass,
        E_func=fields.zero_field(), B_func=fields.gradient_B(B0, L),
        dt=dt, n_steps=n,
    )
    # The gyroradius dwarfs the per-orbit drift, so extract each guiding center
    # by averaging the y-position over one gyro-period (80 steps) at the start
    # and end of the run; the drift is the difference.
    per = 80
    drift_full = pos_full[-per:, 1].mean() - pos_full[:per, 1].mean()
    drift_gc = pos_gc[-per:, 1].mean() - pos_gc[:per, 1].mean()
    assert np.isclose(drift_full, drift_gc, rtol=0.05)
    assert np.sign(drift_full) == np.sign(drift_gc)
