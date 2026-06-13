"""T0 — the 3-D tokamak equilibrium field (`plasmaplay.tokamak`).

The bridge from the experiment-04 poloidal flux ψ(R, Z) to a 3-D Cartesian B(x).
Every check here is one of the three falsifiable T0 validations from
`docs/3D_TOKAMAK_GUIDE.md`:

  * ∇·B ≈ 0 numerically (a magnetic field must be divergence-free),
  * |B| ∝ 1/R for the vacuum toroidal field (the defining tokamak falloff),
  * on the magnetic axis the field is purely toroidal (∇ψ = 0 there),

plus that the field composes with the existing pushers/tracers unchanged.
"""

import numpy as np

from plasmaplay.diagnostics import trace_field_line
from plasmaplay.solvers import grad_shafranov_solve
from plasmaplay.tokamak import (
    divergence,
    equilibrium_field,
    safety_factor,
    solovev_F,
    to_cartesian,
    to_cylindrical,
    toroidal_poincare,
    vacuum_F,
)

R0, B0 = 1.0, 2.0


def _equilibrium_grid(n=121):
    """A Solov'ev-style fixed-boundary equilibrium (experiment 04's setup)."""
    R = np.linspace(R0 - 0.55, R0 + 0.55, n)
    Z = np.linspace(-0.7, 0.7, n)
    RR, _ = np.meshgrid(R, Z, indexing="ij")
    source = -(1.0 * RR**2 + 1.0)            # c_p R² + c_0, both positive
    psi = grad_shafranov_solve(R, Z, source, boundary=0.0)
    return R, Z, psi


# --- the field construction matches the analytic B_R, B_Z, B_φ -----------

def test_field_matches_analytic_flux_derivatives():
    # ψ = R² Z  ->  B_R = -(1/R)∂ψ/∂Z = -R,  B_Z = (1/R)∂ψ/∂R = 2Z.
    # Both are linear, so grid samples and bilinear interpolation are exact.
    R = np.linspace(0.6, 1.6, 81)
    Z = np.linspace(-0.6, 0.6, 81)
    RR, ZZ = np.meshgrid(R, Z, indexing="ij")
    psi = RR**2 * ZZ
    field = equilibrium_field(R, Z, psi, vacuum_F(R0, B0))

    for (r, phi, z) in [(1.0, 0.0, 0.2), (1.2, 1.1, -0.3), (0.9, -2.0, 0.1)]:
        pos = to_cartesian(r, phi, z)
        B = field(pos)
        BR_exp, BZ_exp, Bphi_exp = -r, 2.0 * z, R0 * B0 / r
        Bx = BR_exp * np.cos(phi) - Bphi_exp * np.sin(phi)
        By = BR_exp * np.sin(phi) + Bphi_exp * np.cos(phi)
        np.testing.assert_allclose(B, [Bx, By, BZ_exp], atol=1e-9)


# --- |B| ∝ 1/R for the vacuum toroidal field -----------------------------

def test_vacuum_field_is_one_over_R():
    # Flat ψ -> B_R = B_Z = 0 exactly; B = (R0 B0 / R) φ̂, so |B| R = R0 B0.
    R = np.linspace(0.6, 1.6, 41)
    Z = np.linspace(-0.5, 0.5, 41)
    psi = np.zeros((R.size, Z.size))
    field = equilibrium_field(R, Z, psi, vacuum_F(R0, B0))

    for r in (0.7, 0.95, 1.25, 1.5):
        for phi in (0.0, 1.3, -2.5):
            B = field(to_cartesian(r, phi, 0.0))
            np.testing.assert_allclose(np.linalg.norm(B) * r, R0 * B0, rtol=1e-12)
            # and it is purely toroidal (no poloidal component)
            bhat_phi = np.array([-np.sin(phi), np.cos(phi), 0.0])
            np.testing.assert_allclose(B @ bhat_phi, np.linalg.norm(B), rtol=1e-12)


# --- ∇·B ≈ 0 in a real equilibrium ---------------------------------------

def test_divergence_free_in_equilibrium():
    R, Z, psi = _equilibrium_grid(n=161)
    field = equilibrium_field(R, Z, psi, vacuum_F(R0, B0))
    h = 2.0 * (R[1] - R[0])                  # span a grid cell, not intra-cell

    rng = np.random.default_rng(0)
    divs, scales = [], []
    for _ in range(40):
        r = rng.uniform(R0 - 0.35, R0 + 0.35)
        z = rng.uniform(-0.45, 0.45)
        phi = rng.uniform(-np.pi, np.pi)
        pos = to_cartesian(r, phi, z)
        divs.append(abs(divergence(field, pos, h=h)))
        # characteristic |∂B/∂x| at this point, to judge "small" relatively
        Bp = field(pos + np.array([h, 0, 0]))
        Bm = field(pos - np.array([h, 0, 0]))
        scales.append(np.linalg.norm(Bp - Bm) / (2 * h))
    rel = np.array(divs) / np.maximum(np.array(scales), 1e-12)
    assert np.max(rel) < 0.05                # div ≪ the field's own gradient scale


# --- on the magnetic axis the field is purely toroidal -------------------

def test_axis_is_purely_toroidal():
    R, Z, psi = _equilibrium_grid(n=161)
    i, j = np.unravel_index(np.argmax(psi), psi.shape)   # axis = ψ extremum
    R_axis, Z_axis = R[i], Z[j]
    field = equilibrium_field(R, Z, psi, vacuum_F(R0, B0))

    for phi in (0.0, 0.8, -1.7):
        B = field(to_cartesian(R_axis, phi, Z_axis))
        bhat_phi = np.array([-np.sin(phi), np.cos(phi), 0.0])
        B_tor = B @ bhat_phi
        B_pol = np.linalg.norm(B - B_tor * bhat_phi)
        assert B_pol / abs(B_tor) < 0.02     # ∇ψ ≈ 0 at the axis -> no poloidal B


# --- F-profile helpers and coordinate transforms -------------------------

def test_solovev_F_reduces_to_vacuum():
    psi = np.linspace(-1.0, 1.0, 11)
    np.testing.assert_allclose(solovev_F(R0 * B0, 0.0)(psi), vacuum_F(R0, B0)(psi))


def test_solovev_F_profile():
    # F(ψ) = sqrt(F0² + 2 FF'(ψ - ψ_b))
    F0, FFp, psi_b = 2.0, 0.5, 0.0
    psi = np.array([0.0, 1.0, 2.0])
    expected = np.sqrt(F0**2 + 2 * FFp * (psi - psi_b))
    np.testing.assert_allclose(solovev_F(F0, FFp, psi_b)(psi), expected)


def test_coordinate_roundtrip():
    rng = np.random.default_rng(1)
    for _ in range(20):
        R = rng.uniform(0.5, 2.0)
        phi = rng.uniform(-np.pi, np.pi)
        Zc = rng.uniform(-1.0, 1.0)
        Rb, phib, Zb = to_cylindrical(to_cartesian(R, phi, Zc))
        np.testing.assert_allclose([Rb, phib, Zb], [R, phi, Zc], atol=1e-12)


# --- the field drops straight into the existing tracer (T1 preview) -------

def test_field_line_stays_bounded_on_a_flux_surface():
    R, Z, psi = _equilibrium_grid(n=161)
    field = equilibrium_field(R, Z, psi, solovev_F(R0 * B0, 0.0))
    i, j = np.unravel_index(np.argmax(psi), psi.shape)
    R_axis, Z_axis = R[i], Z[j]
    # seed just off the axis; the line should wind on a nested surface, staying
    # inside the plasma region (it must not run off to the grid edge).
    x0 = to_cartesian(R_axis + 0.1, 0.0, Z_axis)
    pts = trace_field_line(field, x0, ds=0.01, n_steps=4000)
    Rline = np.hypot(pts[:, 0], pts[:, 1])
    assert Rline.min() > R0 - 0.55 and Rline.max() < R0 + 0.55
    assert np.all(np.abs(pts[:, 2]) < 0.7)


# --- T1: q-profile & toroidal Poincaré -----------------------------------

def _circular_equilibrium(Rc=10.0, B0c=1.0, a=1.0, Bth0=0.05, n=241):
    """Large-aspect-ratio circular equilibrium with uniform current.

    ψ = C·r², r² = (R−Rc)² + Z², C = Rc·Bθ0/(2a) gives B_θ = Bθ0·(r/a), so the
    safety factor is uniform: q = a·B0 /(Rc·Bθ0). Flux surfaces are exact circles.
    """
    C = Rc * Bth0 / (2.0 * a)
    R = np.linspace(Rc - 1.2, Rc + 1.2, n)
    Z = np.linspace(-1.2, 1.2, n)
    RR, ZZ = np.meshgrid(R, Z, indexing="ij")
    psi = C * ((RR - Rc) ** 2 + ZZ ** 2)
    field = equilibrium_field(R, Z, psi, vacuum_F(Rc, B0c))
    q_analytic = a * B0c / (Rc * Bth0)
    return field, (Rc, 0.0), q_analytic


def test_safety_factor_matches_large_aspect_analytic():
    # ds=0.05 is as accurate as ds=0.02 here (field lines are exact circles, so
    # RK4 arc-length tracing is converged) but ~2.5× faster — see the T1 memo.
    field, axis, q_exact = _circular_equilibrium()      # q_exact = 2.0
    for r in (0.4, 0.7, 1.0):
        q = safety_factor(field, (axis[0] + r, 0.0), axis, n_poloidal=6, ds=0.05)
        assert abs(abs(q) - q_exact) / q_exact < 0.03   # within 3% of analytic q


def test_safety_factor_converges_with_trace_length():
    field, axis, _ = _circular_equilibrium()
    q4 = safety_factor(field, (axis[0] + 0.7, 0.0), axis, n_poloidal=4, ds=0.05)
    q10 = safety_factor(field, (axis[0] + 0.7, 0.0), axis, n_poloidal=10, ds=0.05)
    assert abs(q4 - q10) / abs(q10) < 0.02              # q is a converged number


def _sheared_circular_equilibrium(Rc=10.0, B0c=1.0, a=1.0, q0=1.0, n=261):
    """Large-aspect circular equilibrium with q(r) ≈ q₀(1+(r/a)²) (sheared).

    ψ(r) = (B0c a²/2q₀) ln(1+(r/a)²). Because q varies with radius it is
    irrational on essentially every surface, so a traced field line fills its
    circle *densely* — the honest test of "closed nested surface". (The uniform
    q = 2 equilibrium above is rational: a line closes after 2 toroidal turns and
    punctures φ=0 at just 2 points, which would pass the std/mean check trivially.
    See docs/T1_QPROFILE_POINCARE.md.)
    """
    R = np.linspace(Rc - 1.2, Rc + 1.2, n)
    Z = np.linspace(-1.2, 1.2, n)
    RR, ZZ = np.meshgrid(R, Z, indexing="ij")
    r2 = (RR - Rc) ** 2 + ZZ ** 2
    psi = (B0c * a**2 / (2.0 * q0)) * np.log1p(r2 / a**2)
    return equilibrium_field(R, Z, psi, vacuum_F(Rc, B0c)), (Rc, 0.0)


def test_poincare_is_a_closed_nested_surface():
    field, axis = _sheared_circular_equilibrium()
    # one surface: a line on it densely fills a single circle — every one of many
    # distinct punctures sits at the same minor radius.
    pc = toroidal_poincare(field, (axis[0] + 0.7, 0.0), n_punctures=60, ds=0.06)
    minor = np.hypot(pc[:, 0] - axis[0], pc[:, 1])
    poloidal = np.arctan2(pc[:, 1], pc[:, 0] - axis[0])
    assert minor.std() / minor.mean() < 2e-3
    assert np.ptp(poloidal) > 5.0          # punctures spread around the circle (≈2π)
    # two surfaces are nested: the inner launch stays inside the outer
    inner = toroidal_poincare(field, (axis[0] + 0.4, 0.0), n_punctures=30, ds=0.06)
    r_inner = np.hypot(inner[:, 0] - axis[0], inner[:, 1]).mean()
    assert r_inner < minor.mean()


def test_safety_factor_shear_increases_outward():
    # a real (shaped, tight-aspect) Solov'ev equilibrium has magnetic shear:
    # |q| grows monotonically from the axis toward the edge.
    R, Z, psi = _equilibrium_grid(n=161)
    i, j = np.unravel_index(np.argmax(psi), psi.shape)
    axis = (R[i], Z[j])
    field = equilibrium_field(R, Z, psi, vacuum_F(R0, B0))
    qs = [abs(safety_factor(field, (axis[0] + r, axis[1]), axis,
                            n_poloidal=8, ds=0.05))
          for r in (0.1, 0.2, 0.3, 0.4)]
    assert all(qs[k] < qs[k + 1] for k in range(len(qs) - 1))
    assert qs[0] > 0.5                                  # finite q on/near the axis
