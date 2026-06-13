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

from plasmaplay.constants import e, m_p
from plasmaplay.diagnostics import trace_field_line
from plasmaplay.guiding_center import gc_push, magnetic_moment
from plasmaplay.pushers import boris_push
from plasmaplay.solvers import grad_shafranov_solve
from plasmaplay.tokamak import (
    divergence,
    equilibrium_field,
    helical_perturbation,
    safety_factor,
    solovev_F,
    superpose,
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


# --- T2: banana orbits, trapping boundary, μ invariance ------------------

# A large-aspect circular equilibrium: |B| = Rc B0 / R is an exact magnetic
# mirror around each flux surface, so the trapped/passing boundary is the clean
# analytic λ_c = √(2ε/(1+ε)). ZERO = no electric field.
_TRAP_RC, _TRAP_B0 = 10.0, 1.0
def _ZERO(x):
    return np.zeros(3)


def _trap_field(n=201):
    Bth0, aa = 0.05, 1.0
    C = _TRAP_RC * Bth0 / (2.0 * aa)
    R = np.linspace(_TRAP_RC - 1.2, _TRAP_RC + 1.2, n)
    Z = np.linspace(-1.2, 1.2, n)
    RR, ZZ = np.meshgrid(R, Z, indexing="ij")
    psi = C * ((RR - _TRAP_RC) ** 2 + ZZ ** 2)
    return equilibrium_field(R, Z, psi, vacuum_F(_TRAP_RC, _TRAP_B0))


def _gc_orbit(field, r, lam, energy_eV=1000.0, T=2e-3, nst=2000):
    """Guiding-center orbit launched at the outboard midplane (Rc+r, 0)."""
    v = np.sqrt(2.0 * energy_eV * e / m_p)
    x0 = to_cartesian(_TRAP_RC + r, 0.0, 0.0)
    B0 = np.linalg.norm(field(x0))
    vpar0 = lam * v
    vperp = np.sqrt(max(v * v - vpar0 * vpar0, 0.0))
    mu = magnetic_moment(vperp, m_p, B0)
    _, pos, vpar = gc_push(x0, vpar0, mu, e, m_p, _ZERO, field, T / nst, nst)
    return pos, vpar


def _is_trapped(field, r, lam):
    _, vpar = _gc_orbit(field, r, lam)
    return np.sign(vpar).min() != np.sign(vpar).max()   # v∥ reverses -> bounced


def test_trapping_boundary_matches_sqrt_epsilon():
    # The mirror prediction λ_c = √(2ε/(1+ε)) must separate trapped from passing:
    # a pitch 15% below it traps (bounces), 15% above it passes. Checked at three
    # minor radii — which also shows the boundary (= trapped fraction) grows with ε.
    field = _trap_field()
    for r in (0.3, 0.5, 0.7):
        eps = r / _TRAP_RC
        lam_c = np.sqrt(2 * eps / (1 + eps))
        assert _is_trapped(field, r, 0.85 * lam_c)          # below boundary -> trapped
        assert not _is_trapped(field, r, 1.15 * lam_c)      # above boundary -> passing


def test_trapped_bounces_passing_circulates():
    field = _trap_field()
    r = 0.5
    # deeply trapped: v∥ reverses and the orbit is a banana — poloidal angle stays
    # bounded (it never gets around to the inboard side).
    pos_t, vpar_t = _gc_orbit(field, r, 0.10)
    theta_t = np.arctan2(pos_t[:, 2], np.hypot(pos_t[:, 0], pos_t[:, 1]) - _TRAP_RC)
    assert np.sign(vpar_t).min() != np.sign(vpar_t).max()
    assert np.ptp(theta_t) < np.pi                          # confined poloidally
    # passing: v∥ keeps its sign and the particle circulates a full poloidal turn.
    pos_p, vpar_p = _gc_orbit(field, r, 0.6)
    theta_p = np.unwrap(np.arctan2(pos_p[:, 2],
                                   np.hypot(pos_p[:, 0], pos_p[:, 1]) - _TRAP_RC))
    assert np.sign(vpar_p).min() == np.sign(vpar_p).max()
    assert abs(theta_p[-1] - theta_p[0]) > 2 * np.pi        # completes a circuit


def test_mu_is_adiabatic_invariant():
    # μ = m v⊥²/(2|B|) is conserved over many gyro-orbits along a full Boris orbit.
    field = _trap_field()
    v = np.sqrt(2.0 * 1000.0 * e / m_p)
    x0 = to_cartesian(_TRAP_RC + 0.5, 0.0, 0.0)
    Bv = field(x0); B0 = np.linalg.norm(Bv); bhat = Bv / B0
    vpar0 = 0.12 * v; vperp = np.sqrt(v * v - vpar0 * vpar0)
    e1 = np.cross(bhat, [0, 0, 1.0]); e1 /= np.linalg.norm(e1)
    vel = vpar0 * bhat + vperp * e1
    wc = e * B0 / m_p
    _, pos, vels = boris_push(x0, vel, e, m_p, _ZERO, field, 2 * np.pi / wc / 30, 12000)
    s = slice(0, len(pos), 50)
    Bm = np.array([np.linalg.norm(field(p)) for p in pos[s]])
    bh = np.array([field(p) / np.linalg.norm(field(p)) for p in pos[s]])
    vp = vels[s]
    vparr = np.einsum("ij,ij->i", vp, bh)
    mu = m_p * (np.sum(vp**2, axis=1) - vparr**2) / (2 * Bm)
    assert (mu.max() - mu.min()) / mu.mean() < 0.01         # μ steady to <1%


# --- T3: break axisymmetry — magnetic islands & stochasticity ------------

# A sheared circular equilibrium (Rc=5) tuned so q = 2 sits at minor radius
# r = 0.5: an m/n = 2/1 helical perturbation is then resonant there and tears
# the surface into a two-island chain. Aspect ratio 5 (not 10) keeps the field-
# line traces fast — see docs/T3_MAGNETIC_ISLANDS.md.
_T3_RC, _T3_Q0 = 5.0, 1.6
_T3_AXIS = (_T3_RC, 0.0)
def _T3_ENV(r):  # Gaussian envelope localising the perturbation at the resonance
    return np.exp(-((r - 0.5) ** 2) / (2 * 0.3**2))


def _t3_base_field(n=261):
    R = np.linspace(_T3_RC - 1.15, _T3_RC + 1.15, n)
    Z = np.linspace(-1.15, 1.15, n)
    RR, ZZ = np.meshgrid(R, Z, indexing="ij")
    psi = (1.0 * 1.0 / (2.0 * _T3_Q0)) * np.log1p(((RR - _T3_RC) ** 2 + ZZ ** 2) / 1.0)
    return equilibrium_field(R, Z, psi, vacuum_F(_T3_RC, 1.0))


def _island_width(base, delta, start, npunc=70, ds=0.07):
    field = superpose(base, helical_perturbation(delta, 2, 1, _T3_AXIS, envelope=_T3_ENV))
    pc = toroidal_poincare(field, start, n_punctures=npunc, ds=ds)
    minor = np.hypot(pc[:, 0] - _T3_RC, pc[:, 1])
    return minor.max() - minor.min(), pc


def test_resonance_is_at_q_equals_m_over_n():
    base = _t3_base_field()
    q = abs(safety_factor(base, (_T3_RC + 0.5, 0.0), _T3_AXIS, n_poloidal=6, ds=0.06))
    assert abs(q - 2.0) < 0.05                     # the q = m/n = 2 surface is at r ≈ 0.5


def test_perturbation_stays_divergence_free():
    # the helical perturbation is built from a flux, so ∇·B = 0 must survive it.
    base = _t3_base_field()
    field = superpose(base, helical_perturbation(2e-3, 2, 1, _T3_AXIS, envelope=_T3_ENV))
    h = 2.0 * (1.15 * 2 / 260)
    for r, phi, z in [(0.5, 0.0, 0.0), (0.45, 1.0, 0.1), (0.55, -2.0, -0.05)]:
        BR = np.linalg.norm(field(to_cartesian(_T3_RC + r, phi, z)))
        d = abs(divergence(field, to_cartesian(_T3_RC + r, phi, z), h=h))
        assert d / BR < 1e-3                       # ∇·B ≪ |B|


def test_island_opens_centered_on_the_resonant_surface():
    # The 2/1 perturbation tears the q=2 surface into an island chain: a line
    # started at its X-point fills a wide band *centred on r = 0.5* (where q=2),
    # not somewhere else. An unperturbed surface there would be a thin circle.
    base = _t3_base_field()
    W_res, pc = _island_width(base, 1e-3, (_T3_RC, 0.5))
    minor = np.hypot(pc[:, 0] - _T3_RC, pc[:, 1])
    assert W_res > 0.1                                          # a real, wide island
    assert abs(minor.mean() - 0.5) < 0.06                      # straddles q=2 surface
    # an irrational surface well inside stays a thin KAM curve (no island there)
    W_off, _ = _island_width(base, 1e-3, (_T3_RC + 0.85, 0.0))
    assert W_off < W_res                                        # resonance is the wide one


def test_island_width_scales_as_sqrt_delta():
    # W ∝ √(perturbation amplitude): quadrupling δ doubles the island width.
    base = _t3_base_field()
    W1, _ = _island_width(base, 5e-4, (_T3_RC, 0.5))
    W4, _ = _island_width(base, 2e-3, (_T3_RC, 0.5))
    assert 1.8 < W4 / W1 < 2.2                                  # √4 = 2
