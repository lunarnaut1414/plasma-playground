"""V8 — 1-D ideal MHD solver validation, plus converter/flux unit tests.

Three physics checks of increasing specificity:
  - Sod hydro shock tube (B=0) vs the well-known exact star-region values,
  - a circularly-polarized Alfvén wave (an exact MHD solution) propagating at
    v_A = Bx/√ρ with its amplitude preserved,
  - the Brio-Wu MHD shock tube: positivity, conservation, and the signature
    By sign reversal across the structure.
"""

import numpy as np

from plasmaplay import fvm

GAMMA = 2.0


# --- unit-level: converters & flux ----------------------------------------

def test_prim_cons_roundtrip():
    W = np.array([[1.0, 0.2, -0.1, 0.3, 0.5, 0.4, -0.2],
                  [0.5, -0.3, 0.0, 0.1, 2.0, -0.6, 0.3]])
    U = fvm.prim_to_cons(W, Bx=0.7, gamma=GAMMA)
    W2 = fvm.cons_to_prim(U, Bx=0.7, gamma=GAMMA)
    np.testing.assert_allclose(W, W2, atol=1e-14)


def test_uniform_state_is_steady():
    n = 100
    W0 = np.tile([1.0, 0.3, 0.0, 0.0, 1.0, 0.5, 0.0], (n, 1))
    Wf = fvm.solve_mhd_1d(W0, dx=1.0 / n, Bx=0.7, gamma=GAMMA,
                          t_end=0.1, bc="periodic")
    np.testing.assert_allclose(Wf, W0, atol=1e-12)


def test_fast_speed_reduces_to_sound_speed_without_field():
    # with B = 0 the fast magnetosonic speed is just the sound speed sqrt(gamma p/rho)
    W = np.array([[2.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0]])
    cf = fvm.fast_speed(W, Bx=0.0, gamma=GAMMA)
    assert np.isclose(cf[0], np.sqrt(GAMMA * 4.0 / 2.0))


# --- V8a: Sod hydro shock tube vs known exact solution --------------------

def test_v8_sod_shock_tube_star_region():
    n, dx = 800, 1.0 / 800
    x = (np.arange(n) + 0.5) * dx
    W0 = np.zeros((n, 7))
    left = x < 0.5
    W0[left] = [1.0, 0, 0, 0, 1.0, 0, 0]
    W0[~left] = [0.125, 0, 0, 0, 0.1, 0, 0]
    Wf = fvm.solve_mhd_1d(W0, dx, Bx=0.0, gamma=1.4, t_end=0.2)
    star = Wf[(x > 0.58) & (x < 0.68)]            # the contact/star plateau
    # exact Sod star values: p* = 0.30313, u* = 0.92745
    assert np.isclose(np.median(star[:, 4]), 0.30313, rtol=1e-2)
    assert np.isclose(np.median(star[:, 1]), 0.92745, rtol=1e-2)


# --- V8b: circularly-polarized Alfven wave propagates at v_A --------------

def test_v8_alfven_wave_speed_and_amplitude():
    n, L = 256, 1.0
    dx = L / n
    x = (np.arange(n) + 0.5) * dx
    A, k, rho, Bx, p = 1e-3, 2 * np.pi / L, 1.0, 1.0, 0.1
    By0, Bz0 = A * np.cos(k * x), A * np.sin(k * x)
    W0 = np.zeros((n, 7))
    W0[:, 0] = rho
    W0[:, 4] = p
    W0[:, 5], W0[:, 6] = By0, Bz0
    W0[:, 2], W0[:, 3] = -By0 / np.sqrt(rho), -Bz0 / np.sqrt(rho)  # +x Alfven wave

    v_A = Bx / np.sqrt(rho)
    t = 0.25
    Wf = fvm.solve_mhd_1d(W0, dx, Bx=Bx, gamma=GAMMA, t_end=t, bc="periodic")

    shift = int(round(v_A * t / dx))              # cells the wave should travel
    # the actual best-fit shift recovers the speed
    errs = [np.sum((Wf[:, 5] - np.roll(By0, s)) ** 2) for s in range(n)]
    assert int(np.argmin(errs)) == shift          # speed = v_A
    amp_ratio = np.max(np.abs(Wf[:, 5])) / A
    assert 0.9 < amp_ratio <= 1.02                # amplitude preserved (low diffusion)


# --- V8c: Brio-Wu MHD shock tube ------------------------------------------

def test_v8_brio_wu_structure():
    n, dx = 800, 1.0 / 800
    x = (np.arange(n) + 0.5) * dx
    W0 = np.zeros((n, 7))
    left = x < 0.5
    W0[left] = [1.0, 0, 0, 0, 1.0, 1.0, 0]
    W0[~left] = [0.125, 0, 0, 0, 0.1, -1.0, 0]
    mass0 = W0[:, 0].sum() * dx
    Wf = fvm.solve_mhd_1d(W0, dx, Bx=0.75, gamma=GAMMA, t_end=0.1)

    assert (Wf[:, 0] > 0).all() and (Wf[:, 4] > 0).all()    # positivity
    assert Wf[:, 5].min() < -0.1 and Wf[:, 5].max() > 0.1   # By reverses sign
    # waves don't reach the ends in t=0.1, so far states are untouched ...
    assert np.isclose(Wf[0, 0], 1.0, rtol=1e-3)
    assert np.isclose(Wf[-1, 0], 0.125, rtol=1e-3)
    # ... and total mass is conserved
    assert np.isclose(Wf[:, 0].sum() * dx, mass0, rtol=1e-3)
