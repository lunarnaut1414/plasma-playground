"""V9 — kinetic Langmuir dispersion; V15 — ω–k spectrum recovery.

V9 solves the kinetic dispersion relation and checks the real frequency and the
Landau damping rate against the standard literature benchmark values. V15 checks
the ω–k FFT diagnostic recovers a known plane wave's (k, ω).
"""

import numpy as np

from plasmaplay.diagnostics import omega_k_spectrum
from plasmaplay.dispersion import (
    bohm_gross,
    langmuir_dispersion,
    plasma_dispersion_function,
)


# --- V9: kinetic Langmuir dispersion --------------------------------------

# (k λ_D, ω_r/ω_pe, γ/ω_pe) — classic electron-plasma-wave Landau benchmarks
LANDAU_BENCHMARKS = [
    (0.5, 1.4156, -0.1533),
    (0.4, 1.2850, -0.0661),
    (0.3, 1.1598, -0.0126),
]


def test_v9_langmuir_frequency_and_damping():
    for kld, wr_ref, g_ref in LANDAU_BENCHMARKS:
        wr, g = langmuir_dispersion(kld)
        assert np.isclose(wr, wr_ref, rtol=1e-3)
        assert np.isclose(g, g_ref, rtol=2e-2)
        assert g < 0                              # waves damp


def test_v9_damping_weakens_at_long_wavelength():
    # Landau damping is exponentially weak as k λ_D -> 0
    _, g_short = langmuir_dispersion(0.5)
    _, g_long = langmuir_dispersion(0.2)
    assert abs(g_long) < abs(g_short) < 1.0
    assert abs(g_long) < 1e-3


def test_v9_real_part_near_bohm_gross_for_small_k():
    # at small k λ_D the kinetic frequency approaches the Bohm-Gross fluid value
    wr, _ = langmuir_dispersion(0.15)
    assert np.isclose(wr, bohm_gross(0.15), rtol=2e-2)


def test_plasma_dispersion_function_known_value():
    # Z(0) = i√π
    assert np.isclose(plasma_dispersion_function(0.0), 1j * np.sqrt(np.pi))


def test_plasma_dispersion_function_vs_plasmapy():
    import pytest
    pytest.importorskip("plasmapy")
    from plasmapy.dispersion import plasma_dispersion_func
    for z in (0.5, 1.0 + 0.3j, -0.7j):
        assert np.isclose(plasma_dispersion_function(z), plasma_dispersion_func(z),
                          rtol=1e-8)


# --- V15: omega-k spectrum recovers a known plane wave ---------------------

def test_v15_omega_k_recovers_plane_wave():
    n_x, n_t = 128, 256
    L, T = 2 * np.pi, 10.0
    dx, dt = L / n_x, T / n_t
    x = np.arange(n_x) * dx
    t = np.arange(n_t) * dt
    k0 = 2 * np.pi * 5 / L          # mode 5
    omega0 = 2 * np.pi * 8 / T      # 8 periods over T
    field = np.cos(k0 * x[None, :] - omega0 * t[:, None])

    k, omega, power = omega_k_spectrum(field, dx, dt)
    iw, ik = np.unravel_index(np.argmax(power), power.shape)
    assert np.isclose(abs(k[ik]), k0, atol=1.5 * (k[1] - k[0]))
    assert np.isclose(abs(omega[iw]), omega0, atol=1.5 * (omega[1] - omega[0]))
