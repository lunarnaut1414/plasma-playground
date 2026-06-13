"""Physics validation of the full PIC loop (experiment 03): V5.

V5 — a perturbed *cold* plasma oscillates at the plasma frequency ω_pe,
independent of wavelength. We run in normalized units (ε0 = 1) where the
macroparticle charge/mass are chosen so ω_pe = 1 exactly, then measure the
oscillation frequency of the seeded Fourier mode from its time-spectrum peak.

(V6 Landau damping and V7 two-stream growth live with the experiment, since they
need more care with noise.)
"""

import numpy as np

from plasmaplay import pic
from plasmaplay.diagnostics import dominant_frequency


def _omega_pe_one_setup(n_particles, L):
    """Choose q, m so that ω_pe = 1 in normalized (ε0 = 1) units."""
    q = -L / n_particles      # -> mean |charge density| = 1
    m = L / n_particles       # -> q/m = -1  => ω_pe^2 = 1
    return q, m


def _mode_history(sim, mode, dt, n_steps):
    """Record the dominant real component of the mode-`mode` E amplitude."""
    comp = np.empty(n_steps, dtype=complex)
    for n in range(n_steps):
        E = sim.step(dt)
        comp[n] = np.fft.fft(E)[mode]
    # the spatial phase puts the signal in real or imag depending on seeding;
    # pick whichever carries the oscillation.
    return comp.real if comp.real.std() > comp.imag.std() else comp.imag


def test_v5_cold_plasma_oscillates_at_omega_pe():
    rng = np.random.default_rng(0)
    L, N, NG = 2 * np.pi, 20000, 64
    x, v = pic.load_maxwellian(N, L, v_thermal=0.0, rng=rng)   # cold
    x = pic.perturb_positions(x, L, amplitude=0.05, mode=1)
    q, m = _omega_pe_one_setup(N, L)
    sim = pic.ElectrostaticPIC1D(L, NG, x, v, q, m, eps0=1.0)

    dt, n_steps = 0.05, 3000
    signal = _mode_history(sim, 1, dt, n_steps)
    omega = dominant_frequency(signal, dt)
    assert np.isclose(omega, 1.0, rtol=0.05)        # ω_pe = 1


def test_v5_frequency_independent_of_mode():
    # cold oscillation frequency is ω_pe for any wavelength (no thermal disp.)
    rng = np.random.default_rng(1)
    L, N, NG = 2 * np.pi, 20000, 64
    q, m = _omega_pe_one_setup(N, L)
    dt, n_steps = 0.05, 3000
    for mode in (1, 2, 3):
        x, v = pic.load_maxwellian(N, L, v_thermal=0.0, rng=rng)
        x = pic.perturb_positions(x, L, amplitude=0.05, mode=mode)
        sim = pic.ElectrostaticPIC1D(L, NG, x, v, q, m, eps0=1.0)
        signal = _mode_history(sim, mode, dt, n_steps)
        omega = dominant_frequency(signal, dt)
        assert np.isclose(omega, 1.0, rtol=0.05)
