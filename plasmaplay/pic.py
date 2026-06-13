"""Electrostatic particle-in-cell (PIC) building blocks, 1-D.

The PIC cycle each step:

    deposit charge to grid  ->  solve for E (Poisson)  ->  gather E to particles
    ->  push particles (leapfrog)

This module provides the particle<->grid weighting (cloud-in-cell), distribution
loaders, and a small `ElectrostaticPIC1D` stepper that wires those together with
the spectral field solver from `solvers.py`. It is deliberately plain NumPy and
periodic; `numba` can accelerate the deposit/gather later if needed.
"""

from __future__ import annotations

import numpy as np
from scipy.special import erfinv

from .constants import EPSILON_0
from .solvers import solve_efield_1d


# --- particle <-> grid weighting (cloud-in-cell) --------------------------

def cic_deposit(positions, charges, n_grid, L):
    """Deposit particle charge onto a periodic grid (cloud-in-cell).

    Each particle splits linearly between its two nearest grid points. Returns
    charge *density* (charge / length) of shape (n_grid,). `charges` may be a
    scalar (same for all) or a per-particle array.
    """
    dx = L / n_grid
    xp = np.asarray(positions, dtype=float) / dx          # position in cells
    i = np.floor(xp).astype(int)
    frac = xp - i
    i0 = i % n_grid
    i1 = (i + 1) % n_grid
    q = np.asarray(charges, dtype=float)
    # bincount is the fast vectorized scatter (much quicker than np.add.at)
    w0 = q * (1.0 - frac)
    w1 = q * frac
    rho = (np.bincount(i0, weights=np.broadcast_to(w0, i0.shape), minlength=n_grid)
           + np.bincount(i1, weights=np.broadcast_to(w1, i1.shape), minlength=n_grid))
    return rho / dx


def cic_interpolate(grid_field, positions, L):
    """Gather a grid field to particle positions (cloud-in-cell, periodic).

    The exact linear adjoint of `cic_deposit`: same weights, so momentum-
    conserving when paired with it. Returns one value per particle.
    """
    n_grid = len(grid_field)
    dx = L / n_grid
    xp = np.asarray(positions, dtype=float) / dx
    i = np.floor(xp).astype(int)
    frac = xp - i
    i0 = i % n_grid
    i1 = (i + 1) % n_grid
    return grid_field[i0] * (1.0 - frac) + grid_field[i1] * frac


# --- distribution loaders -------------------------------------------------

def load_maxwellian(n_particles, L, v_thermal, drift=0.0, rng=None, quiet=True):
    """Load particles uniform in x with a Maxwellian (normal) velocity spread.

    quiet=True uses a low-noise "quiet start": positions evenly spaced and
    velocities placed on the inverse Maxwellian CDF (exact marginal), then the
    velocities are shuffled to decorrelate them from position. This dramatically
    cuts the particle noise that would otherwise swamp e.g. Landau damping.
    """
    rng = np.random.default_rng() if rng is None else rng
    if quiet:
        x = (np.arange(n_particles) + 0.5) / n_particles * L
        quant = (np.arange(n_particles) + 0.5) / n_particles
        v = drift + v_thermal * np.sqrt(2.0) * erfinv(2.0 * quant - 1.0)
        rng.shuffle(v)
    else:
        x = rng.uniform(0.0, L, n_particles)
        v = rng.normal(drift, v_thermal, n_particles)
    return x, v


def perturb_positions(x, L, amplitude, mode=1):
    """Displace positions to seed a density perturbation δn/n ≈ amplitude·cos(kx)."""
    k = 2.0 * np.pi * mode / L
    return (x - (amplitude / k) * np.sin(k * x)) % L


def load_two_stream(n_particles, L, v_beam, v_thermal=0.0, rng=None, quiet=True):
    """Two counter-streaming beams (±v_beam), n_particles total."""
    half = n_particles // 2
    x1, v1 = load_maxwellian(half, L, v_thermal, drift=+v_beam, rng=rng, quiet=quiet)
    x2, v2 = load_maxwellian(half, L, v_thermal, drift=-v_beam, rng=rng, quiet=quiet)
    return np.concatenate([x1, x2]), np.concatenate([v1, v2])


# --- the PIC stepper ------------------------------------------------------

class ElectrostaticPIC1D:
    """Minimal periodic 1-D electrostatic PIC.

    Parameters
    ----------
    L, n_grid : domain length and number of grid cells.
    positions, velocities : initial particle state (length N).
    charge, mass : per-macroparticle charge and mass (q/m sets the dynamics).
    eps0 : permittivity (use 1.0 for normalized runs where ω_pe is set by choice
           of charge/mass — see experiment 03).
    """

    def __init__(self, L, n_grid, positions, velocities, charge, mass,
                 eps0=EPSILON_0):
        self.L = L
        self.n_grid = n_grid
        self.dx = L / n_grid
        self.x = np.mod(np.asarray(positions, dtype=float), L)
        self.v = np.asarray(velocities, dtype=float).copy()
        self.q = charge
        self.m = mass
        self.eps0 = eps0

    def efield(self):
        """Current electric field on the grid."""
        rho = cic_deposit(self.x, self.q, self.n_grid, self.L)
        return solve_efield_1d(rho, self.dx, self.eps0)

    def step(self, dt):
        """Advance one leapfrog step; return the grid E used this step."""
        E_grid = self.efield()
        E_p = cic_interpolate(E_grid, self.x, self.L)
        self.v += (self.q / self.m) * E_p * dt
        self.x = np.mod(self.x + self.v * dt, self.L)
        return E_grid

    def field_energy(self):
        """Electrostatic field energy ½ ε0 ∫E² dx."""
        E = self.efield()
        return 0.5 * self.eps0 * np.sum(E ** 2) * self.dx


def mode_amplitude(grid_field, mode):
    """Complex Fourier amplitude of a given integer mode of a grid field."""
    return np.fft.fft(grid_field)[mode]
