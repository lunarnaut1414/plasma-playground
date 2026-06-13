"""Experiment 03 — 1-D electrostatic Particle-in-Cell.

The first *self-consistent, many-particle* experiment: the particles make the
field that pushes them. Here that buys us collective behavior with no analogue in
single-particle land.

Cases (choose with --case):
  cold        a perturbed cold plasma oscillates at the plasma frequency ω_pe (F1)
  landau      a warm Langmuir wave damps collisionlessly — Landau damping   (F2)
  twostream   two counter-streaming beams go unstable — two-stream           (F2)

We run in normalized units (ε0 = 1) with the macroparticle charge/mass chosen so
ω_pe = 1, which makes every frequency/rate directly comparable to theory.

Run:
    python run.py --case cold [--save]
    python run.py --case landau
    python run.py --case twostream
    python run.py --case all --save
"""

from __future__ import annotations

import argparse

import matplotlib.pyplot as plt
import numpy as np

from plasmaplay import pic, plotting
from plasmaplay.diagnostics import dominant_frequency

L = 2 * np.pi          # domain length -> fundamental mode k = 1
NG = 64                # grid cells


def omega_pe_one(n_particles):
    """q, m giving ω_pe = 1 in ε0 = 1 units (see module docstring)."""
    return -L / n_particles, L / n_particles


def _energies(sim):
    ke = 0.5 * sim.m * np.sum(sim.v ** 2)
    fe = sim.field_energy()
    return ke, fe


# --- F1: cold plasma oscillation -----------------------------------------

def run_cold(save=False):
    print("\n--- cold plasma oscillation (F1) ---")
    rng = np.random.default_rng(0)
    N = 20000
    # A *truly* cold plasma (v_th = 0) has zero Debye length, so the PIC grid is
    # always under-resolved (dx > π λ_D) and a numerical grid-heating instability
    # grows. A small thermal spread sets λ_D ≈ 0.05 > dx/π and stabilizes it,
    # with a negligible Bohm-Gross shift to the frequency. (Resolve λ_D!)
    v_th = 0.05
    x, v = pic.load_maxwellian(N, L, v_thermal=v_th, rng=rng)
    x = pic.perturb_positions(x, L, amplitude=0.05, mode=1)
    q, m = omega_pe_one(N)
    sim = pic.ElectrostaticPIC1D(L, NG, x, v, q, m, eps0=1.0)

    dt, n_steps = 0.05, 2000
    t = np.arange(n_steps) * dt
    mode1 = np.empty(n_steps)
    ke = np.empty(n_steps)
    fe = np.empty(n_steps)
    for n in range(n_steps):
        E = sim.step(dt)
        mode1[n] = np.fft.fft(E)[1].imag      # mode-1 amplitude (sine phase)
        ke[n], fe[n] = _energies(sim)

    omega = dominant_frequency(mode1, dt)
    print(f"  measured ω = {omega:.4f}   (ω_pe = 1.000)")
    total = ke + fe
    print(f"  energy drift = {np.ptp(total) / total.mean():.2e}")

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.5))
    a1.plot(t, mode1, lw=0.8)
    a1.set_xlabel("t  [1/ω_pe]"); a1.set_ylabel("mode-1 E amplitude")
    a1.set_title(f"Cold oscillation at ω_pe  (measured ω={omega:.3f})")
    a2.plot(t, fe, label="field"); a2.plot(t, ke - ke[0], label="kinetic − KE₀")
    a2.plot(t, fe + ke - ke[0], "k--", lw=0.8, label="total")
    a2.set_xlabel("t  [1/ω_pe]"); a2.set_ylabel("energy")
    a2.set_title("Energy exchange (field ↔ kinetic)"); a2.legend(fontsize=8)
    fig.tight_layout()
    _finish(fig, "cold_oscillation.png", save)


def _loglinear_rate(t, amp, t0, t1):
    """Fit log(amp) ~ rate*t over the window [t0, t1]; return the slope."""
    mask = (t >= t0) & (t <= t1)
    slope = np.polyfit(t[mask], np.log(amp[mask]), 1)[0]
    return slope


# --- F2: Landau damping --------------------------------------------------

# Landau damping of a Langmuir wave at k λ_D = 0.5 (the standard benchmark):
# the kinetic dispersion gives ω_r ≈ 1.4156 ω_pe and γ ≈ 0.1533 ω_pe.
LANDAU_GAMMA = 0.1533


def run_landau(save=False):
    print("\n--- Landau damping (F2) ---")
    rng = np.random.default_rng(0)
    N, NG_l = 200000, 128
    v_th = 0.5                       # k λ_D = 0.5 at k = 1
    x, v = pic.load_maxwellian(N, L, v_thermal=v_th, rng=rng)
    x = pic.perturb_positions(x, L, amplitude=0.01, mode=1)
    q, m = omega_pe_one(N)
    sim = pic.ElectrostaticPIC1D(L, NG_l, x, v, q, m, eps0=1.0)

    dt, n_steps = 0.05, 700
    t = np.arange(n_steps) * dt
    amp = np.empty(n_steps)
    for n in range(n_steps):
        E = sim.step(dt)
        amp[n] = np.abs(np.fft.fft(E)[1])

    gamma = -_loglinear_rate(t, amp, 2.0, 12.0)
    print(f"  measured γ = {gamma:.4f}   (theory γ = {LANDAU_GAMMA:.4f})")

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.semilogy(t, amp, lw=0.8, label="mode-1 |E| (PIC)")
    ax.semilogy(t, amp[40] * np.exp(-LANDAU_GAMMA * (t - t[40])), "r--",
                label=f"theory e^(−γt), γ={LANDAU_GAMMA}")
    ax.set_xlabel("t  [1/ω_pe]"); ax.set_ylabel("mode-1 |E|")
    ax.set_title("Landau damping (collisionless)"); ax.legend()
    ax.set_ylim(amp.min() * 0.5, amp.max() * 2)
    fig.tight_layout()
    _finish(fig, "landau_damping.png", save)


# --- F2: two-stream instability ------------------------------------------

def two_stream_growth_rate(k, v0, omega_pe=1.0):
    """Max growth rate of the symmetric cold two-stream dispersion.

    1 = ω_b²/(ω−kv0)² + ω_b²/(ω+kv0)²  with ω_b² = ω_pe²/2.  Reduces to a
    quadratic in u = ω²; a negative real root means ω is imaginary -> growth.
    """
    a = k * v0
    wb2 = 0.5 * omega_pe ** 2
    u = np.roots([1.0, -(2 * a**2 + 2 * wb2), a**4 - 2 * wb2 * a**2])
    omegas = np.concatenate([np.sqrt(u + 0j), -np.sqrt(u + 0j)])
    return float(np.max(omegas.imag))


def run_twostream(save=False):
    print("\n--- two-stream instability (F2) ---")
    rng = np.random.default_rng(0)
    N, NG_t = 200000, 128
    v0, v_th = 0.6, 0.03            # k v0 = 0.6 < ω_pe -> mode 1 unstable
    x, v = pic.load_two_stream(N, L, v_beam=v0, v_thermal=v_th, rng=rng)
    x = pic.perturb_positions(x, L, amplitude=0.01, mode=1)
    q, m = omega_pe_one(N)
    sim = pic.ElectrostaticPIC1D(L, NG_t, x, v, q, m, eps0=1.0)

    dt, n_steps = 0.05, 600
    t = np.arange(n_steps) * dt
    amp = np.empty(n_steps)
    snaps = {}
    for n in range(n_steps):
        E = sim.step(dt)
        amp[n] = np.abs(np.fft.fft(E)[1])
        if n in (0, n_steps - 1):
            snaps[n] = (sim.x.copy(), sim.v.copy())

    gamma = _loglinear_rate(t, amp, 4.0, 14.0)
    gamma_th = two_stream_growth_rate(1.0, v0)
    print(f"  measured γ = {gamma:.4f}   (theory γ = {gamma_th:.4f})")

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 5))
    a1.semilogy(t, amp, lw=0.9, label="mode-1 |E| (PIC)")
    a1.semilogy(t, amp[80] * np.exp(gamma_th * (t - t[80])), "r--",
                label=f"theory e^(+γt), γ={gamma_th:.3f}")
    a1.set_xlabel("t  [1/ω_pe]"); a1.set_ylabel("mode-1 |E|")
    a1.set_title("Two-stream growth"); a1.legend()
    a1.set_ylim(amp[amp > 0].min() * 0.5, amp.max() * 2)

    xf, vf = snaps[n_steps - 1]
    idx = rng.choice(len(xf), 8000, replace=False)
    a2.scatter(xf[idx], vf[idx], s=1, alpha=0.3)
    a2.set_xlabel("x"); a2.set_ylabel("v")
    a2.set_title("Phase space at saturation (the 'hole')")
    fig.tight_layout()
    _finish(fig, "two_stream.png", save)


def _finish(fig, name, save):
    if save:
        out = plotting.ensure_outputs_dir(__file__)
        path = out / name
        fig.savefig(path, dpi=150)
        print(f"  saved -> {path}")
    plt.show()


def main(case="cold", save=False):
    print("=" * 64)
    print("Experiment 03 — 1-D electrostatic PIC")
    print("=" * 64)
    if case in ("cold", "all"):
        run_cold(save)
    if case in ("landau", "all"):
        run_landau(save)
    if case in ("twostream", "all"):
        run_twostream(save)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--case", default="cold",
                   choices=["cold", "landau", "twostream", "all"])
    p.add_argument("--save", action="store_true")
    args = p.parse_args()
    main(case=args.case, save=args.save)
