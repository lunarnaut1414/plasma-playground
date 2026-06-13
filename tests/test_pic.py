"""Unit tests for the PIC weighting kernels (cloud-in-cell) and loaders.

These pin the properties the PIC cycle relies on: charge conservation on deposit,
partition of unity (a constant field gathers back to the constant), and linear
exactness of CIC interpolation. The physics-level validation of the full loop
(ω_pe, Landau damping) lives with experiment 03 (V5–V7).
"""

import numpy as np

from plasmaplay import pic

L = 2 * np.pi
NG = 64


def test_deposit_conserves_total_charge():
    rng = np.random.default_rng(0)
    x = rng.uniform(0, L, 1000)
    q = 0.5
    rho = pic.cic_deposit(x, q, NG, L)
    # integral of density = total charge
    assert np.isclose(rho.sum() * (L / NG), 1000 * q)


def test_deposit_on_grid_point_is_local():
    # a particle sitting exactly on grid node j deposits entirely to node j
    dx = L / NG
    x = np.array([5 * dx])
    rho = pic.cic_deposit(x, 1.0, NG, L)
    assert np.isclose(rho[5], 1.0 / dx)
    assert np.isclose(rho.sum(), 1.0 / dx)   # nothing leaked elsewhere


def test_interpolate_constant_field_partition_of_unity():
    field = np.full(NG, 3.7)
    rng = np.random.default_rng(1)
    x = rng.uniform(0, L, 500)
    vals = pic.cic_interpolate(field, x, L)
    np.testing.assert_allclose(vals, 3.7, atol=1e-12)


def test_interpolate_is_linearly_exact():
    # CIC is linear, so a field sampling a linear ramp interpolates exactly
    # (away from the periodic wrap). Build E = a*x on the grid.
    dx = L / NG
    xs_grid = np.arange(NG) * dx
    a = 2.0
    field = a * xs_grid
    x = np.linspace(0.1, L - 0.1, 50)        # avoid the wrap cell
    vals = pic.cic_interpolate(field, x, L)
    np.testing.assert_allclose(vals, a * x, atol=1e-10)


def test_deposit_interpolate_are_adjoint_weights():
    # gather then a matching scatter use identical weights -> momentum conserving.
    # Check the weight identity: sum_i W(i) = 1 for any particle.
    dx = L / NG
    x = np.array([3.27 * dx, 0.0, (NG - 0.001) * dx])
    for xi in x:
        rho = pic.cic_deposit(np.array([xi]), 1.0, NG, L) * dx  # back to charge
        assert np.isclose(rho.sum(), 1.0)       # weights sum to one


def test_loaders_shapes_and_statistics():
    rng = np.random.default_rng(2)
    x, v = pic.load_maxwellian(10000, L, v_thermal=1.5, drift=0.3, rng=rng)
    assert x.shape == v.shape == (10000,)
    assert x.min() >= 0 and x.max() < L
    assert np.isclose(v.mean(), 0.3, atol=0.05)
    assert np.isclose(v.std(), 1.5, rtol=0.05)


def test_two_stream_has_two_beams():
    rng = np.random.default_rng(3)
    x, v = pic.load_two_stream(10000, L, v_beam=2.0, v_thermal=0.1, rng=rng)
    assert np.isclose(v.mean(), 0.0, atol=0.05)
    # roughly half above zero, half below
    assert abs((v > 0).mean() - 0.5) < 0.02


def test_perturbation_creates_density_mode():
    # seeding mode-1 should put power in the mode-1 Fourier component of density
    rng = np.random.default_rng(4)
    x, _ = pic.load_maxwellian(20000, L, v_thermal=1.0, rng=rng)
    xp = pic.perturb_positions(x, L, amplitude=0.1, mode=1)
    rho = pic.cic_deposit(xp, 1.0, NG, L)
    spectrum = np.abs(np.fft.fft(rho - rho.mean()))
    assert spectrum[1] == spectrum[1:NG // 2].max()    # mode 1 dominates
