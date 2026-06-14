"""The Kadomtsev sawtooth cycle (rung B3) — periodic m=1 core reconnection.

B1 found *when* the m=1/n=1 internal kink goes unstable (q(0) < 1); B3 is what
happens *repeatedly* once it does — the **sawtooth oscillation** that paces the core
of every ohmically-heated tokamak:

  1. ohmic current (lower resistivity in the hot core) peaks J(r) and drives q(0)
     below 1 on the resistive timescale tau_R;
  2. the m=1 internal kink reconnects the core — Kadomtsev's "complete reconnection"
     — **flattening** the temperature (and q -> 1) inside the mixing radius r_mix and
     crashing the central temperature;
  3. q(0) relaxes back above 1, the core re-heats and the current re-peaks, and the
     whole thing repeats. The period scales with tau_R.

This module models that on a 1-D cylinder, reusing the B1 q-profile / kink trigger
(`cylinder_mhd`). The re-peaking is *genuine* resistive induction-equation diffusion
of the poloidal field (the current peaks where the resistivity is low); the crash is
the Kadomtsev flatten. The clean, exact invariant we check is **thermal-energy
conservation** across the crash (the reconnection redistributes heat, it does not
destroy it); the helical flux psi* is computed and flattened in the core (the
reconnection signature). The sawtooth period scales with tau_R = a^2 / eta.

Normalisation: minor radius a = 1, B_z = R = 1 (so q = r / B_theta). Resistivity eta
sets the resistive time tau_R = a^2 / eta.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import cumulative_trapezoid

from .cylinder_mhd import screw_pinch_q


def helical_flux(r, q):
    """The (1,1) helical flux psi*(r) = int_0^r B_theta (1 - q) dr'  (B_z = R = 1).

    Rises from 0 on axis, peaks at the q = 1 surface (where the integrand vanishes),
    then falls — the non-monotonic profile whose reconnection is the sawtooth crash.
    """
    r = np.asarray(r, dtype=float)
    q = np.asarray(q, dtype=float)
    b_th = np.where(r > 1e-12, r / np.maximum(q, 1e-12), 0.0)
    return cumulative_trapezoid(b_th * (1.0 - q), r, initial=0.0)


def mixing_radius(r, q):
    """Kadomtsev mixing radius r_mix: the outer radius where psi* returns to its axis
    value (0). Reconnection mixes the core out to here. Returns None if q(0) >= 1
    (no q = 1 surface -> no crash)."""
    q = np.asarray(q, dtype=float)
    if q[0] >= 1.0:
        return None
    psi = helical_flux(r, q)
    i1 = int(np.argmax(psi))                       # the q = 1 surface (psi* maximum)
    after = np.where(psi[i1:] <= 0.0)[0]
    return float(r[-1]) if after.size == 0 else float(r[i1 + after[0]])


def kadomtsev_flatten(r, field, r_mix):
    """Flatten `field` inside r_mix to its area-weighted mean, conserving the
    cross-section integral int field * r dr over [0, r_mix] **exactly** (so flattening
    a temperature conserves its thermal energy). Outside r_mix the field is untouched.
    """
    r = np.asarray(r, dtype=float)
    out = np.asarray(field, dtype=float).copy()
    inside = r <= r_mix
    if inside.sum() < 2:
        return out
    w = r[inside]
    mean = np.trapezoid(out[inside] * w, r[inside]) / np.trapezoid(w, r[inside])
    out[inside] = mean
    return out


def q_from_temperature(rho, T, q_edge=3.0, alpha=1.5):
    """Model safety-factor profile q(rho) from a Spitzer-ohmic current J ~ T^alpha.

    A hot core conducts better (Spitzer sigma ~ T^{3/2}), so the ohmic current peaks
    where the plasma is hot; a peaked current lowers q on axis. With B_theta(rho) =
    (1/rho) int_0^rho J rho' drho' and q ~ rho/B_theta, normalised so q(1) = q_edge
    (fixed total current),

        q(rho) = q_edge * G(1) * rho^2 / G(rho),   G(rho) = int_0^rho (T/T_edge)^alpha rho' drho'.

    The point of the coupling (Track C): as the burn peaks the core temperature, q(0)
    falls — and when it crosses 1 the m=1 internal kink fires a sawtooth crash.
    """
    rho = np.asarray(rho, dtype=float)
    T = np.asarray(T, dtype=float)
    j = np.maximum(T / max(T[-1], 1e-9), 1e-6) ** alpha    # ohmic current shape
    G = cumulative_trapezoid(j * rho, rho, initial=0.0)
    G = np.maximum(G, 1e-30)
    q = q_edge * G[-1] * rho ** 2 / G
    q[0] = q[1]                                            # finite on-axis value
    return q


def resistive_relaxation_time(T_core_keV, tau_ref=2.0, T_ref_keV=20.0):
    """Sawtooth recovery time = the resistive current-redistribution time near the axis.

    After a crash the central current profile is reset (q(0) ~ 1); the next crash can
    only fire once the *current* has resistively diffused back into a peaked, kink-
    unstable profile. That time is the local resistive diffusion time tau_R =
    mu0 a_1^2 / eta. Spitzer resistivity eta ~ T^{-3/2}, so

        tau_R(T) = tau_ref * (T_core / T_ref)^{3/2}.

    A HOTTER core conducts better -> current diffuses *slower* -> a LONGER sawtooth
    period and a bigger build-up between crashes ("monster sawteeth"). This is why a
    burning reactor shows a few large crashes, not many small ones — the opposite of a
    naive "crash whenever q(0)<1 every step" model. `tau_ref` is the period at the
    reference core temperature `T_ref_keV` (a device-scale calibration, not first
    principles). Returns seconds.
    """
    return tau_ref * (max(float(T_core_keV), 1e-3) / T_ref_keV) ** 1.5


def external_q_profile(rho, q_axis=1.4, q_edge=2.6):
    """A stellarator's safety factor: set by the EXTERNAL coils, not plasma current.

    In a stellarator the rotational transform comes from the 3-D coil geometry (see
    `fields.helical_stellarator`) with ~zero net plasma current, so q(rho) is fixed by
    the device — it does NOT respond to the burning core the way the current-driven
    tokamak q does (`q_from_temperature`). Modelled as a smooth monotone profile with
    q_axis on axis and q_edge at the edge. Kept **above 1 everywhere** (q_axis > 1), so
    there is no q = 1 surface — hence no m=1 internal kink, no sawteeth, no disruptions:
    the inherently steady-state operation that distinguishes a stellarator from a tokamak.
    """
    rho = np.asarray(rho, dtype=float)
    return q_axis + (q_edge - q_axis) * rho ** 2


def crash_profiles(rho, n, T, r_mix):
    """Apply a Kadomtsev crash to a (n, T) burn state inside r_mix, conserving BOTH
    the particle content (int n rho drho) and the thermal energy (int 3 n T rho drho)
    exactly. Returns (n_new, T_new): density and energy density are each flattened to
    their area-weighted means inside r_mix, then T = energy/(3 n)."""
    n_new = kadomtsev_flatten(rho, n, r_mix)               # conserves particles
    w_new = kadomtsev_flatten(rho, 3.0 * n * T, r_mix)     # conserves thermal energy
    T_new = np.asarray(T, dtype=float).copy()
    inside = np.asarray(rho) <= r_mix
    T_new[inside] = w_new[inside] / (3.0 * np.maximum(n_new[inside], 1e9))
    return n_new, T_new


def sawtooth_event(rho, n, T, *, q_edge=3.0, alpha=1.5, q_trigger=0.93,
                   min_radius=0.05):
    """Fire a sawtooth crash IF the Spitzer q-profile of this state has q(0) below the
    trigger (the m=1 internal kink reaches finite amplitude slightly below q(0)=1, the
    observed crash threshold ~0.9-0.95).

    Returns (n, T, crashed): the post-crash profiles (particle- and energy-conserving)
    and whether a crash occurred. The reusable MHD-event trigger for the coupled
    discharge — it leaves the state untouched when the core is kink-stable. The
    trigger margin is what makes crashes finite-amplitude and well-separated rather
    than firing every step at the marginal point."""
    q = q_from_temperature(rho, T, q_edge, alpha)
    if q[0] >= q_trigger:
        return n, T, False
    r_mix = mixing_radius(rho, q)
    if r_mix is None or r_mix < min_radius:
        return n, T, False
    n2, T2 = crash_profiles(rho, n, T, r_mix)
    return n2, T2, True


class SawtoothCycle:
    """A 1-D ohmic sawtooth: resistive current peaking + Kadomtsev crashes.

    The poloidal field B_theta(r) evolves by the resistive induction equation with a
    core-peaked conductivity (low resistivity in the hot centre), so the ohmic steady
    current peaks on axis and drives q(0) below 1; when it crosses the trigger the core
    reconnects (Kadomtsev), q is reset toward 1 and the temperature is flattened inside
    r_mix conserving thermal energy. The temperature is heated centrally and re-peaks
    between crashes — the sawtooth waveform. tau_R = a^2 / eta.
    """

    def __init__(self, *, eta=0.5, nr=97, q0_init=1.2, nu=1.0, chi=0.02,
                 heat=0.6, eta_peaking=8.0, q_trigger=0.99, q_reset=1.05):
        self.r = np.linspace(0.0, 1.0, nr)
        self.dr = self.r[1] - self.r[0]
        self.eta = eta
        self.chi = chi
        self.heat = heat
        self.q_trigger = q_trigger
        self.q_reset = q_reset
        # resistivity LOW in the hot core, rising outward (Spitzer-like) -> the ohmic
        # steady current peaks on axis, driving q(0) down toward q_ss(0) < 1.
        self.D_r = eta * (1.0 + eta_peaking * self.r ** 2)   # tau_R ~ a^2 / eta
        q_init = screw_pinch_q(self.r, q0_init, nu)
        self.bth = np.where(self.r > 1e-12, self.r / q_init, 0.0)
        self.B_edge = self.bth[-1]                          # fixed total current
        self.T = 1.0 - 0.8 * self.r ** 2                    # peaked starting profile
        self.t = 0.0
        self.crash_times = []

    # --- geometry from the poloidal field -----------------------------------
    def q(self):
        q = np.where(self.r > 1e-9, self.r / np.maximum(self.bth, 1e-12), 0.0)
        q[0] = self.q0()                                    # finite on-axis value
        return q

    def q0(self):
        """Central safety factor (B_theta ~ r/q0 near the axis -> q0 = r/B_theta)."""
        return float(self.r[1] / max(self.bth[1], 1e-12))

    def _current(self, bth):
        """Axial current density J_z = (1/r) d(r B_theta)/dr."""
        J = np.gradient(self.r * bth, self.r) / np.maximum(self.r, 1e-9)
        J[0] = J[1]
        return J

    def step(self, dt):
        # genuine resistive induction (implicit): dB_theta/dt = d/dr( D(r) J_z ),
        # J peaks where D is low (the core) -> q(0) falls to its ohmic steady value <1.
        self.bth = self._diffuse_bth(self.bth, dt)
        # temperature: weak transport + central heating
        self.T = self._diffuse_T(self.T, self.chi, dt)
        self.T = self.T + dt * self.heat * np.exp(-(self.r / 0.35) ** 2)
        self.t += dt
        if self.q0() < self.q_trigger:                      # Kadomtsev crash
            self._crash()
        return self

    def _diffuse_bth(self, u, dt):
        """Backward-Euler step of dB_theta/dt = d/dr[ D(r) (1/r) d/dr(r B_theta) ].

        Finite-volume on the staggered faces (E = D J, J = (1/r) d(r B_theta)/dr),
        tridiagonal in B_theta. Axis B_theta=0, edge held at B_edge (fixed current).
        """
        n = u.size
        r, dr = self.r, self.dr
        rf = 0.5 * (r[:-1] + r[1:])
        Df = 0.5 * (self.D_r[:-1] + self.D_r[1:])
        a = Df / (rf * dr)                                   # face coefficient (len n-1)
        lo = np.zeros(n); di = np.ones(n); up = np.zeros(n); rhs = u.copy()
        for i in range(1, n - 1):
            lo[i] = -dt / dr * a[i - 1] * r[i - 1]
            up[i] = -dt / dr * a[i] * r[i + 1]
            di[i] = 1.0 + dt / dr * r[i] * (a[i] + a[i - 1])
        di[0] = 1.0; up[0] = 0.0; rhs[0] = 0.0              # axis
        di[-1] = 1.0; lo[-1] = 0.0; rhs[-1] = self.B_edge   # fixed total current
        return _thomas(lo, di, up, rhs)

    def _diffuse_T(self, f, coeff, dt):
        """Implicit cylindrical diffusion of the temperature (edge held at 0.05)."""
        n = f.size
        r, dr = self.r, self.dr
        rf = 0.5 * (r[:-1] + r[1:])
        rc = r.copy(); rc[0] = dr / 8.0
        lo = np.zeros(n); di = np.ones(n); up = np.zeros(n); rhs = f.copy()
        for i in range(1, n - 1):
            w = dt * coeff / (rc[i] * dr ** 2)
            lo[i] = -w * rf[i - 1]; up[i] = -w * rf[i]
            di[i] = 1.0 + w * (rf[i - 1] + rf[i])
        w0 = dt * coeff / (rc[0] * dr ** 2)
        di[0] = 1.0 + w0 * rf[0]; up[0] = -w0 * rf[0]
        di[-1] = 1.0; lo[-1] = 0.0; rhs[-1] = 0.05
        return _thomas(lo, di, up, rhs)

    def _crash(self):
        r_mix = mixing_radius(self.r, self.q())
        if r_mix is None or r_mix < 2 * self.dr:
            return
        # flatten T inside r_mix (conserves thermal energy); reconnect q -> q_reset
        self.T = kadomtsev_flatten(self.r, self.T, r_mix)
        inside = self.r <= r_mix
        self.bth = self.bth.copy()
        self.bth[inside] = self.r[inside] / self.q_reset    # q = q_reset (~1) in core
        self.crash_times.append(self.t)

    def thermal_energy(self):
        """Cross-section thermal energy int T * r dr (the crash-conserved invariant)."""
        return float(np.trapezoid(self.T * self.r, self.r))

    def run(self, t_end, dt, *, record_every=1):
        """Integrate, returning dict: t, T0 (core T), q0, and the crash times."""
        ts, T0, q0s = [], [], []
        for k in range(int(round(t_end / dt))):
            self.step(dt)
            if k % record_every == 0:
                ts.append(self.t); T0.append(self.T[0]); q0s.append(self.q0())
        return {"t": np.array(ts), "T0": np.array(T0), "q0": np.array(q0s),
                "crashes": np.array(self.crash_times)}


def _thomas(lower, diag, upper, rhs):
    """Thomas tridiagonal solve, O(N)."""
    n = diag.size
    c = upper.copy(); d = rhs.astype(float).copy(); b = diag.copy()
    c[0] /= b[0]; d[0] /= b[0]
    for i in range(1, n):
        m = b[i] - lower[i] * c[i - 1]
        c[i] = upper[i] / m if i < n - 1 else 0.0
        d[i] = (d[i] - lower[i] * d[i - 1]) / m
    x = d
    for i in range(n - 2, -1, -1):
        x[i] -= c[i] * x[i + 1]
    return x
