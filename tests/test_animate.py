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

import matplotlib
import numpy as np
from PIL import Image

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from plasmaplay import animate as anim  # noqa: E402


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


def test_animate_profiles_shade_between(tmp_path):
    """Two series with a shaded gap band (the Ti-Te equipartition lag) animate cleanly."""
    x = np.linspace(0, 1, 20)
    te = np.array([(10 + k) * (1 - x**2) for k in range(6)])
    ti = np.array([(18 + k) * (1 - x**2) for k in range(6)])
    frames = np.stack([te, ti], axis=1)                       # (6, 2, 20)
    out = anim.animate_profiles(x, frames, np.arange(6.0), path=tmp_path / "two.gif",
                                labels=["Te", "Ti"], colors=["#22d3ee", "#ff9f45"],
                                shade_between=(0, 1), shade_label="Ti-Te", fps=6, dpi=60)
    assert out.exists() and Image.open(out).n_frames == 6


def test_animate_cross_section_writes_gif(tmp_path):
    rho = np.linspace(0, 1, 41)
    times = np.linspace(0, 3, 10)
    # a peaked profile that decays in time
    frames = anim.make_frames(lambda t: np.exp(-(rho / 0.4) ** 2) * np.exp(-0.3 * t), times)
    out = anim.animate_cross_section(rho, frames, times, path=tmp_path / "x.gif",
                                     clabel="T", title="disk", fps=10, dpi=60)
    assert out.exists() and out.stat().st_size > 0


def test_animate_phase_track_writes_gif(tmp_path):
    t = np.linspace(0, 10, 20)
    n = 1e20 * (1 + 0.1 * t)
    T = 2 + 2 * t
    out = anim.animate_phase_track(n, T, t, path=tmp_path / "track.gif",
                                   color=np.linspace(0, 0.1, 20), xlabel="n",
                                   ylabel="T", clabel="f_He", band=(10.0, 20.0),
                                   fps=8, dpi=60)
    assert out.exists() and out.stat().st_size > 0


def test_animate_torus_3d_writes_gif(tmp_path):
    edge = np.linspace(1.0, 5.0, 8)
    out = anim.animate_torus_3d(edge, path=tmp_path / "torus.gif",
                                title="torus", fps=8, dpi=60)
    assert out.exists() and out.stat().st_size > 0


def test_nested_torus_surface_radii(tmp_path):
    """Each nested flux surface is a torus of tube radius rho*a (it satisfies the
    torus equation at that radius), and the animation writes a gif."""
    R0, a = 3.0, 1.0
    for rho in (0.3, 0.7, 1.0):
        X, Y, Z = anim.torus_surface(R0, rho * a, n_u=40, n_v=24)
        resid = (np.sqrt(X**2 + Y**2) - R0) ** 2 + Z**2 - (rho * a) ** 2
        assert np.max(np.abs(resid)) < 1e-9
    rho_levels = np.array([0.2, 0.6, 1.0])
    T_rt = np.array([[10.0, 5.0, 1.0], [20.0, 8.0, 1.0], [4.0, 4.0, 1.0]])  # incl. a crash
    out = anim.animate_torus_nested(rho_levels, T_rt, path=tmp_path / "nested.gif",
                                    fps=6, dpi=60, n_u=30, n_v=18)
    assert out.exists() and out.stat().st_size > 0


def test_apply_house_style_dark_and_light():
    """The shared style sets a dark background + light text when dark, no-op otherwise."""
    fig, ax = plt.subplots()
    txt = anim.apply_house_style(fig, ax, dark=True)
    assert txt == anim.HOUSE_FG
    assert fig.get_facecolor()[:3] != (1.0, 1.0, 1.0)        # not white
    plt.close(fig)
    fig2, ax2 = plt.subplots()
    assert anim.apply_house_style(fig2, [ax2], dark=False) == "black"
    plt.close(fig2)


def test_animate_poloidal_field_writes_gif(tmp_path):
    """A masked 2-D field on an (R,Z) grid animates with a colorbar + boundary outline."""
    R = np.linspace(2.0, 4.0, 16); Z = np.linspace(-1.5, 1.5, 18)
    RR, ZZ = np.meshgrid(R, Z, indexing="ij")
    mask = ((RR - 3.0) ** 2 + ZZ**2) < 1.0
    frames = np.array([(5 + 3 * k) * np.exp(-((RR - 3) ** 2 + ZZ**2)) for k in range(5)])
    out = anim.animate_poloidal_field(R, Z, frames, np.arange(5.0), path=tmp_path / "p.gif",
                                      mask=mask, clabel="T", fps=5, dpi=60)
    assert out.exists() and Image.open(out).n_frames == 5


def test_animate_operating_space_writes_gif(tmp_path):
    """Several tracks sweep an (n,T) operating diagram with a band + Greenwald vline;
    the dark-themed movie writes a multi-frame gif."""
    t = np.linspace(0, 1, 8)
    tracks = [
        {"x": 5e19 + 1e19 * t, "y": 5 + 0 * t, "label": "L", "color": "#22d3ee"},
        {"x": 7e19 + 1e20 * t, "y": 22 - 20 * t, "label": "disrupt", "color": "#c0c5cf"},
    ]
    out = anim.animate_operating_space(tracks, t, path=tmp_path / "ops.gif",
                                       vlines=[(2.2e20, "n_G")], band=(10.0, 25.0),
                                       fps=6, dpi=60)
    assert out.exists() and Image.open(out).n_frames == 8


def test_torus_field_lines_on_surface():
    """Circular-tube field lines lie on a torus of radius a*rscale and there are n_lines."""
    R0, a, rscale = 3.0, 1.0, 1.1
    lines = anim.torus_field_lines(R0, a, iota=0.5, n_lines=3, n_tor=2.0, rscale=rscale)
    assert len(lines) == 3
    for X, Y, Z in lines:
        resid = (np.sqrt(X**2 + Y**2) - R0) ** 2 + Z**2 - (a * rscale) ** 2
        assert np.max(np.abs(resid)) < 1e-9


def test_poloidal_bp_quiver_is_circulating():
    """B_p arrows are perpendicular to the radius (purely circulating, no radial part)."""
    X, Y, U, V = anim.poloidal_bp_quiver(rings=(0.5,), n_ang=12)
    radial = X * U + Y * V                      # dot of position with arrow
    assert np.allclose(radial, 0.0, atol=1e-12)


def test_animate_stellarator_3d_writes_gif(tmp_path):
    """The two-panel stellarator burn (twisty torus + l=2 elliptical bullseye) animates a
    full profile and writes a multi-frame gif."""
    rho = np.linspace(0.0, 1.0, 12)
    T_rt = np.array([(15.0 + 2.0 * k) * (1.0 - rho**2) + 0.5 for k in range(5)])
    out = anim.animate_stellarator_3d(rho, T_rt, times=np.arange(5.0),
                                      path=tmp_path / "stell.gif", fps=5, dpi=60,
                                      n_u=40, n_v=18, n_s=16)
    assert out.exists() and Image.open(out).n_frames == 5


def test_animate_discharge_3d_writes_gif(tmp_path):
    """The two-panel 3-D discharge (glowing torus + bullseye T(rho)) animates a full
    profile, accepts a per-frame crash count, and writes a multi-frame gif."""
    rho = np.linspace(0.0, 1.0, 12)
    T_rt = np.array([(20.0 - 5.0 * k % 7) * (1.0 - rho**2) + 0.5 for k in range(6)])
    crashes = np.array([0, 1, 0, 2, 0, 1])
    out = anim.animate_discharge_3d(rho, T_rt, times=np.arange(6.0), crashes=crashes,
                                    path=tmp_path / "disch.gif", fps=6, dpi=60,
                                    n_u=30, n_v=18)
    assert out.exists() and out.stat().st_size > 0
    assert Image.open(out).n_frames == 6
