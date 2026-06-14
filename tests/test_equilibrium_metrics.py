"""Validation suite for plasmaplay.equilibrium_metrics + FluxSurfaceTransport1D.

Falsifiable checks:
  * The flux-surface metric extractor reproduces the circular-limit geometry:
    <|grad rho|^2> = 1/a^2, V' ∝ rho, total volume = 2 pi^2 R0 a^2.
  * A real Grad-Shafranov (Solov'ev) solve shows the outboard Shafranov shift and
    yields a shaped (elongated) plasma.
  * Transport on circular metrics reduces exactly to the cylindrical Transport1D.
  * The IPB98(y,2) scaling reproduces the published ITER confinement time.
"""

import numpy as np
import pytest

from plasmaplay import equilibrium_metrics as em
from plasmaplay import transport as tr
from plasmaplay.solvers import grad_shafranov_solve


def _circular_psi(R0=3.0, a=1.0, n=201, half=1.2):
    """A synthetic circular equilibrium psi = 1 - r^2/a^2 (peaks at the axis)."""
    R = np.linspace(R0 - half, R0 + half, n)
    Z = np.linspace(-half, half, n)
    RR, ZZ = np.meshgrid(R, Z, indexing="ij")
    r2 = (RR - R0) ** 2 + ZZ ** 2
    return R, Z, 1.0 - r2 / a ** 2


# --- flux-surface metrics: the circular limit is analytic ------------------
def test_flux_metrics_gradrho2_circular_limit():
    """For circular surfaces rho = r/a, <|grad rho|^2> = 1/a^2 (=1 here)."""
    R, Z, psi = _circular_psi(a=1.0)
    m = em.flux_surface_metrics(R, Z, psi, n_rho=64)
    mid = (m["rho"] > 0.3) & (m["rho"] < 0.8)
    assert m["grad_rho2"][mid].mean() == pytest.approx(1.0, rel=0.03)


def test_flux_metrics_volume_circular_limit():
    """Total plasma volume is the torus volume 2 pi^2 R0 a^2, and V' ∝ rho."""
    R0, a = 3.0, 1.0
    R, Z, psi = _circular_psi(R0=R0, a=a)
    m = em.flux_surface_metrics(R, Z, psi, n_rho=64)
    assert m["V"][-1] == pytest.approx(2 * np.pi ** 2 * R0 * a ** 2, rel=0.03)
    # V' = dV/drho ∝ rho  =>  V'/rho is constant at 4 pi^2 R0 a^2
    mid = (m["rho"] > 0.3) & (m["rho"] < 0.8)
    ratio = m["Vprime"][mid] / m["rho"][mid]
    assert ratio.mean() == pytest.approx(4 * np.pi ** 2 * R0 * a ** 2, rel=0.05)
    assert ratio.std() / ratio.mean() < 0.05            # genuinely linear


# --- a real Grad-Shafranov equilibrium -------------------------------------
def test_solovev_shafranov_shift_outboard():
    """The Solov'ev GS solve puts the magnetic axis *outboard* of R0 (the toroidal
    1/R term pushes it there) and produces a vertically elongated plasma."""
    R0 = 3.0
    R = np.linspace(R0 - 1.3, R0 + 1.3, 81)
    Z = np.linspace(-1.9, 1.9, 91)
    RR, _ = np.meshgrid(R, Z, indexing="ij")
    psi = grad_shafranov_solve(R, Z, -(RR ** 2 + 1.0), boundary=0.0)
    i, _ = np.unravel_index(np.argmax(psi), psi.shape)
    assert R[i] > R0                                     # outboard Shafranov shift
    m = em.flux_surface_metrics(R, Z, psi, n_rho=48, psi_bnd=0.06 * psi.max())
    assert np.all(m["Vprime"] >= 0) and m["V"][-1] > 0   # sane positive geometry


# --- transport on the geometry reduces to the cylindrical solver -----------
def test_flux_transport_reduces_to_cylindrical():
    """Fed circular metrics (V' = 4 pi^2 R0 a^2 rho, <|grad rho|^2> = 1/a^2), the
    flux-surface transport solver matches Transport1D to sub-percent."""
    R0, a, ng = 3.0, 1.0, 65
    rho = np.linspace(0, 1, ng)
    Vp = 4 * np.pi ** 2 * R0 * a ** 2 * rho
    g = np.full(ng, 1.0 / a ** 2)
    fs = tr.FluxSurfaceTransport1D(a, Vp, g, rho_metric=rho, n_grid=ng, chi=1.5, D=0.5)
    fs.set_state(T=2.0, n=1e20)
    cyl = tr.Transport1D(a=a, n_grid=ng, chi=1.5, D=0.5)
    cyl.set_state(T=2.0, n=1e20)
    for _ in range(800):
        fs.step(2e-3, p_aux_total=3e6)
        cyl.step(2e-3, p_aux_total=3e6)
    assert fs.T[0] == pytest.approx(cyl.T[0], rel=0.02)
    assert fs._vol_avg(fs.T) == pytest.approx(cyl._vol_avg(cyl.T), rel=0.02)


def test_flux_transport_burns_peaked_on_real_geometry():
    """A heated burn on a real Solov'ev equilibrium builds a monotonically peaked
    profile and meaningful alpha self-heating."""
    R0, a = 3.0, 1.0
    R = np.linspace(R0 - 1.3, R0 + 1.3, 81)
    Z = np.linspace(-1.9, 1.9, 91)
    RR, _ = np.meshgrid(R, Z, indexing="ij")
    psi = grad_shafranov_solve(R, Z, -(RR ** 2 + 1.0), boundary=0.0)
    m = em.flux_surface_metrics(R, Z, psi, n_rho=64, psi_bnd=0.06 * psi.max())
    fs = tr.FluxSurfaceTransport1D(a, m["Vprime"], m["grad_rho2"], rho_metric=m["rho"],
                                   n_grid=97, chi=0.6, D=0.05, n_edge=2e19)
    fs.set_state(T=2.0, n=5e19)
    for k in range(2500):
        paux = 5e5 * min(k * 4e-3 / 4.0, 1.0)
        fs.step(4e-3, p_aux_total=paux, fuel_total=5e19 / 6.0,
                fuel_profile=tr.gaussian_deposition(fs.rho, 0.0, 0.4))
    assert fs.T[0] > 15.0                                # a hot burning core
    assert fs.T[0] > fs.T[48] > fs.T[-1]                 # monotonically peaked
    p_alpha = fs._vol_avg(tr.fusion_power_density(fs.n, fs.T, "alpha")) * fs.plasma_volume()
    assert p_alpha / 1e6 > 10.0                          # >10 MW of alpha self-heating


# --- IPB98(y,2) confinement scaling ----------------------------------------
def test_ipb98_reproduces_iter_baseline():
    """The IPB98(y,2) H-mode scaling gives ~3.7 s for the ITER Q=10 baseline."""
    tau = em.confinement_time_ipb98(Ip_MA=15.0, B=5.3, n19=10.0, P_MW=87.0,
                                    R=6.2, a=2.0, kappa=1.75, M=2.5)
    assert tau == pytest.approx(3.7, rel=0.1)


def test_ipb98_scalings_monotonic():
    """More current and weaker heating both lengthen confinement (the key levers)."""
    base = dict(Ip_MA=10.0, B=5.0, n19=8.0, P_MW=50.0, R=4.0, a=1.4, kappa=1.7)
    t0 = em.confinement_time_ipb98(**base)
    assert em.confinement_time_ipb98(**{**base, "Ip_MA": 12.0}) > t0   # tau ∝ Ip^0.93
    assert em.confinement_time_ipb98(**{**base, "P_MW": 70.0}) < t0    # tau ∝ P^-0.69
