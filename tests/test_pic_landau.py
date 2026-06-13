"""Physics validation of the PIC loop, continued: V6 (Landau) & V7 (two-stream).

These collective instabilities are the headline results of experiment 03. PIC
measurements of their rates carry real particle noise, so the tolerances are
deliberately generous (~25-30%); the point is that the simulation reproduces the
right *physics* and the right ballpark rate, not 3-decimal agreement.

Smaller particle counts than the experiment keep the suite fast.
"""

import numpy as np

from plasmaplay import pic

L = 2 * np.pi
NG = 128


def _mode1_amplitude_history(sim, dt, n_steps):
    amp = np.empty(n_steps)
    for n in range(n_steps):
        E = sim.step(dt)
        amp[n] = np.abs(np.fft.fft(E)[1])
    return amp


def _omega_pe_one(N):
    return -L / N, L / N


def _fit_rate(t, amp, t0, t1):
    mask = (t >= t0) & (t <= t1)
    return np.polyfit(t[mask], np.log(amp[mask]), 1)[0]


def test_v6_landau_damping_rate():
    # k λ_D = 0.5 benchmark: kinetic theory gives γ ≈ 0.1533 ω_pe.
    rng = np.random.default_rng(0)
    N = 100000
    x, v = pic.load_maxwellian(N, L, v_thermal=0.5, rng=rng)
    x = pic.perturb_positions(x, L, amplitude=0.01, mode=1)
    q, m = _omega_pe_one(N)
    sim = pic.ElectrostaticPIC1D(L, NG, x, v, q, m, eps0=1.0)

    dt, n_steps = 0.05, 600
    t = np.arange(n_steps) * dt
    amp = _mode1_amplitude_history(sim, dt, n_steps)

    gamma = -_fit_rate(t, amp, 2.0, 11.0)
    assert np.isclose(gamma, 0.1533, rtol=0.35)
    assert gamma > 0          # it genuinely damps


def _two_stream_growth_rate(k, v0, omega_pe=1.0):
    a = k * v0
    wb2 = 0.5 * omega_pe**2
    u = np.roots([1.0, -(2 * a**2 + 2 * wb2), a**4 - 2 * wb2 * a**2])
    omegas = np.concatenate([np.sqrt(u + 0j), -np.sqrt(u + 0j)])
    return float(np.max(omegas.imag))


def test_v7_two_stream_growth_rate():
    v0 = 0.6                  # k v0 = 0.6 < ω_pe -> mode 1 unstable
    gamma_theory = _two_stream_growth_rate(1.0, v0)
    assert gamma_theory > 0.3     # near the max growth rate ω_pe/(2√2)≈0.354

    rng = np.random.default_rng(0)
    N = 100000
    x, v = pic.load_two_stream(N, L, v_beam=v0, v_thermal=0.03, rng=rng)
    x = pic.perturb_positions(x, L, amplitude=0.01, mode=1)
    q, m = _omega_pe_one(N)
    sim = pic.ElectrostaticPIC1D(L, NG, x, v, q, m, eps0=1.0)

    dt, n_steps = 0.05, 500
    t = np.arange(n_steps) * dt
    amp = _mode1_amplitude_history(sim, dt, n_steps)

    gamma = _fit_rate(t, amp, 4.0, 13.0)
    assert np.isclose(gamma, gamma_theory, rtol=0.30)


def test_v7_growth_rate_formula_marginal_stability():
    # the symmetric cold two-stream is stable once k v0 >= ω_pe
    assert _two_stream_growth_rate(1.0, 1.0) < 1e-9     # marginal
    assert _two_stream_growth_rate(1.0, 1.5) < 1e-9     # stable
    assert _two_stream_growth_rate(1.0, 0.5) > 0.0      # unstable
