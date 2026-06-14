"""Validation suite for plasmaplay.animate (the gif/visualization foundation).

We don't assert on gif *bytes*; we assert on the **frame data** (which is the
physics) and that the gif pipeline actually writes a non-empty file.

Falsifiable checks:
  * `make_frames` samples a field on a time grid with the right count/shape.
  * Frames of an analytic 1-D diffusing Gaussian conserve mass and decay in peak
    as 1/√(1+2Dt/σ₀²) — a real property, so this catches a broken sampler/grid.
  * `torus_surface` points satisfy the implicit torus equation.
  * `animate_profiles` / `animate_cross_section` / `animate_torus_3d` each write a
    non-empty .gif with the requested frame pipeline.
"""

import numpy as np

from plasmaplay import animate as anim


# --- analytic reference: 1-D diffusion of a Gaussian -----------------------
def _diffusing_gaussian(x, t, D=0.5, sigma0=1.0, mass=1.0):
    """c(x,t) = mass/√(2π σ²) · exp(-x²/2σ²),  σ²(t) = σ₀² + 2Dt.
    Exact solution of ∂c/∂t = D ∂²c/∂x² for a Gaussian initial condition."""
    sig2 = sigma0**2 + 2.0 * D * t
    return mass / np.sqrt(2.0 * np.pi * sig2) * np.exp(-x**2 / (2.0 * sig2))


# --- make_frames -----------------------------------------------------------
def test_make_frames_shape_and_count():
    x = np.linspace(-8, 8, 129)
    times = np.linspace(0, 4, 25)
    frames = anim.make_frames(lambda t: _diffusing_gaussian(x, t), times)
    assert frames.shape == (25, 129)


def test_diffusion_frames_conserve_mass():
    """Total mass ∫c dx is invariant under diffusion -> frames must agree across time."""
    x = np.linspace(-30, 30, 2001)            # wide grid so the tails are captured
    times = np.linspace(0, 6, 13)
    frames = anim.make_frames(lambda t: _diffusing_gaussian(x, t, D=0.7), times)
    mass = np.trapezoid(frames, x, axis=1)
    assert np.allclose(mass, mass[0], rtol=1e-3)


def test_diffusion_peak_decay_law():
    """Peak amplitude must fall as 1/√(1+2Dt/σ₀²)."""
    x = np.linspace(-30, 30, 4001)
    D, sigma0 = 0.7, 1.0
    times = np.array([0.0, 1.0, 3.0, 6.0])
    frames = anim.make_frames(lambda t: _diffusing_gaussian(x, t, D=D, sigma0=sigma0), times)
    peak = frames.max(axis=1)
    expected = peak[0] / np.sqrt(1.0 + 2.0 * D * times / sigma0**2)
    assert np.allclose(peak, expected, rtol=2e-3)


# --- torus geometry --------------------------------------------------------
def test_torus_surface_on_torus():
    R0, a = 3.0, 1.0
    X, Y, Z = anim.torus_surface(R0, a, n_u=50, n_v=30)
    residual = (np.sqrt(X**2 + Y**2) - R0) ** 2 + Z**2 - a**2
    assert np.max(np.abs(residual)) < 1e-9


# --- gif pipeline (writes to a tmp file, asserts non-empty) ----------------
def test_animate_profiles_writes_gif(tmp_path):
    x = np.linspace(-6, 6, 81)
    times = np.linspace(0, 3, 12)
    frames = anim.make_frames(lambda t: _diffusing_gaussian(x, t), times)
    out = anim.animate_profiles(x, frames, times, path=tmp_path / "p.gif",
                                ylabel="c", title="diffusion", fps=10, dpi=60)
    assert out.exists() and out.stat().st_size > 0


def test_animate_cross_section_writes_gif(tmp_path):
    rho = np.linspace(0, 1, 41)
    times = np.linspace(0, 3, 10)
    # a peaked profile that decays in time
    frames = anim.make_frames(lambda t: np.exp(-(rho / 0.4) ** 2) * np.exp(-0.3 * t), times)
    out = anim.animate_cross_section(rho, frames, times, path=tmp_path / "x.gif",
                                     clabel="T", title="disk", fps=10, dpi=60)
    assert out.exists() and out.stat().st_size > 0


def test_animate_torus_3d_writes_gif(tmp_path):
    edge = np.linspace(1.0, 5.0, 8)
    out = anim.animate_torus_3d(edge, path=tmp_path / "torus.gif",
                                title="torus", fps=8, dpi=60)
    assert out.exists() and out.stat().st_size > 0
