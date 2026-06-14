"""Burning-plasma transport — the discharge arc, not the fluid motion.

"Ignition -> steady state -> fuel injection" is an *energy- and particle-balance*
story that plays out on the confinement timescale (~seconds), NOT an MHD/CFD story
(microseconds). The right model is a transport one: evolve the radial temperature
and density profiles in time, with sources (auxiliary heating, fusion alpha
self-heating, fuelling) and sinks (transport, radiation). This is the toy-model
cousin of integrated codes like ASTRA / RAPTOR / TRANSP.

Two rungs live here, sharing one set of physics terms:

  F0  burn_0d        — 0-D (volume-averaged) power balance. Two coupled ODEs for
                       the plasma energy W and density n. This is the Lawson /
                       POPCON picture and the *validation anchor* for F2: watch it
                       ignite (alpha heating overtakes losses -> thermal runaway)
                       and settle to a burning steady state.

  F2  Transport1D    — 1-D radial transport. Evolve T(rho, t) and n(rho, t) with a
                       diffusion equation, deposition-profile sources, and the same
                       reaction physics. Backward-Euler (implicit) diffusion so the
                       stiff conduction term doesn't dictate a tiny timestep;
                       sources are explicit.

Physics conventions in this module:
  * Temperature T is carried in **keV** in the public API (the natural fusion
    unit); densities n in m^-3; power densities in W/m^3; times in seconds.
  * Single temperature, T_e = T_i = T, and quasineutral n_e = n_i = n with a
    50:50 D-T mix (n_D = n_T = n/2). Z_eff is a scalar knob for radiation.
  * Geometry is a large-aspect-ratio circular column: the flux-surface volume
    element is V'(rho) proportional to rho, so the transport divergence is the
    cylindrical (1/rho) d/drho ( rho * flux ).

The reaction rate is the Bosch-Hale parametrization (H.-S. Bosch & G.M. Hale,
Nucl. Fusion 32, 611 (1992)) — the standard, ~0.3% accurate D-T <sigma v> over
0.2-100 keV. It is the nonlinearity that makes ignition happen.
"""

from __future__ import annotations

import numpy as np

from .constants import e as E  # elementary charge (J per eV)

# ---------------------------------------------------------------------------
# Energetics of the D-T reaction
# ---------------------------------------------------------------------------
E_ALPHA_MEV = 3.5      # alpha particle — stays in the plasma and heats it
E_NEUTRON_MEV = 14.1   # neutron — escapes (this is the wall/blanket power)
E_FUSION_MEV = E_ALPHA_MEV + E_NEUTRON_MEV  # 17.6 MeV total per D-T reaction
_MEV_J = 1.0e6 * E

# Bremsstrahlung coefficient: P_brem [W/m^3] = C_BREM * Z_eff * n_e^2 * sqrt(T_keV)
# (standard relativistic-corrected-free form; T in keV, n in m^-3).
C_BREM = 5.35e-37


# ---------------------------------------------------------------------------
# Bosch-Hale D-T reactivity
# ---------------------------------------------------------------------------
# Coefficients for the D(t,n)4He channel, Bosch & Hale (1992), Table VII.
_BG = 34.3827            # Gamow energy^(1/2)  [keV^(1/2)]
_MRC2 = 1.124656e6       # reduced-mass energy m_r c^2  [keV]
_C = (
    1.17302e-9, 1.51361e-2, 7.51886e-2, 4.60643e-3,
    1.35000e-2, -1.06750e-4, 1.36600e-5,
)


def reactivity_dt(T_keV):
    """D-T fusion reactivity <sigma v> in m^3/s as a function of T (keV).

    Bosch-Hale parametrization, valid ~0.2-100 keV. Vectorized over T. Below the
    fit range the rate is negligible; we clamp tiny/zero T to a floor so the burn
    integrators don't divide by zero. The output is in **m^3/s** (the original
    formula gives cm^3/s; we convert).
    """
    T = np.asarray(T_keV, dtype=float)
    T = np.maximum(T, 1e-3)  # floor: below ~0.2 keV the rate is effectively zero anyway
    c1, c2, c3, c4, c5, c6, c7 = _C
    theta = T / (1.0 - (T * (c2 + T * (c4 + T * c6)))
                 / (1.0 + T * (c3 + T * (c5 + T * c7))))
    xi = (_BG ** 2 / (4.0 * theta)) ** (1.0 / 3.0)
    sigv_cm3 = c1 * theta * np.sqrt(xi / (_MRC2 * T ** 3)) * np.exp(-3.0 * xi)
    return sigv_cm3 * 1.0e-6  # cm^3/s -> m^3/s


# ---------------------------------------------------------------------------
# Power densities (W/m^3) — the source/sink terms shared by F0 and F2
# ---------------------------------------------------------------------------
def fusion_power_density(n, T_keV, which="alpha"):
    """Fusion power density (W/m^3) for a 50:50 D-T plasma of total ion density n.

    n_D = n_T = n/2, so the reaction-rate density is (n/2)^2 <sigma v>. `which`
    selects which slice of the 17.6 MeV is counted:
      "alpha"   -> 3.5 MeV  (the part that self-heats the plasma)
      "neutron" -> 14.1 MeV (the part that leaves — fusion power to the blanket)
      "total"   -> 17.6 MeV
    """
    rate = 0.25 * np.asarray(n, dtype=float) ** 2 * reactivity_dt(T_keV)
    energy = {"alpha": E_ALPHA_MEV, "neutron": E_NEUTRON_MEV, "total": E_FUSION_MEV}[which]
    return rate * energy * _MEV_J


def bremsstrahlung_density(n, T_keV, z_eff=1.0):
    """Bremsstrahlung radiation loss density (W/m^3): C_BREM Z_eff n^2 sqrt(T)."""
    n = np.asarray(n, dtype=float)
    return C_BREM * z_eff * n ** 2 * np.sqrt(np.maximum(T_keV, 0.0))


def _w_of(n, T_keV):
    """Thermal energy density W = (3/2)(n_e T_e + n_i T_i) = 3 n T, in J/m^3."""
    return 3.0 * np.asarray(n, dtype=float) * np.asarray(T_keV, dtype=float) * 1e3 * E


def _T_of(n, w):
    """Invert W = 3 n T for T in keV (guards small n)."""
    n = np.asarray(n, dtype=float)
    return w / (3.0 * np.maximum(n, 1e10) * 1e3 * E)


# ---------------------------------------------------------------------------
# F0 — 0-D burn dynamics (the Lawson / POPCON anchor)
# ---------------------------------------------------------------------------
def burn_0d(n0, T0, *, tau_E, p_aux, tau_p=None, fuel_rate=0.0, z_eff=1.0,
            t_end=8.0, dt=1e-3):
    """Volume-averaged burning-plasma evolution — two coupled ODEs.

        dW/dt = P_aux + P_alpha(n,T) - P_brem(n,T) - W/tau_E
        dn/dt = S_fuel - n/tau_p

    W = 3 n T is the thermal energy density; the energy confinement time tau_E
    rolls all transport losses into one number (the POPCON convention). Ignition
    is the moment the alpha self-heating P_alpha overtakes the losses so the
    plasma stays hot with P_aux removed.

    Parameters
    ----------
    n0 : initial ion density [m^-3].
    T0 : initial temperature [keV].
    tau_E : energy confinement time [s] (scalar, or callable t -> tau_E).
    p_aux : auxiliary heating power density [W/m^3] (scalar, or callable t -> p_aux).
        Make it a callable to ramp heating on then off and watch the plasma stay lit.
    tau_p : particle confinement time [s]; default = 3*tau_E (particles are
        confined a few times longer than energy, a common rule of thumb).
    fuel_rate : fuelling source S_fuel [m^-3 s^-1] (scalar, or callable t).
    z_eff : effective charge for bremsstrahlung.
    t_end, dt : integration window and step [s]. Explicit RK4.

    Returns
    -------
    dict with arrays: t, T (keV), n, W, p_alpha, p_brem, p_aux, Q (= P_fusion/P_aux,
    inf when unheated), and triple (n T tau_E, the Lawson product in keV s / m^3).
    """
    tau_p = 3.0 * tau_E if tau_p is None else tau_p

    def as_f(x):
        return x if callable(x) else (lambda t, _x=x: _x)

    tauE_f, paux_f, fuel_f = as_f(tau_E), as_f(p_aux), as_f(fuel_rate)

    def deriv(t, w, n):
        T = _T_of(n, w)
        p_a = fusion_power_density(n, T, "alpha")
        p_b = bremsstrahlung_density(n, T, z_eff)
        dw = paux_f(t) + p_a - p_b - w / tauE_f(t)
        dn = fuel_f(t) - n / tau_p
        return dw, dn

    nstep = int(round(t_end / dt))
    t = np.empty(nstep + 1)
    W = np.empty(nstep + 1)
    N = np.empty(nstep + 1)
    W[0], N[0], t[0] = _w_of(n0, T0), float(n0), 0.0

    for k in range(nstep):
        tk, wk, nk = t[k], W[k], N[k]
        k1w, k1n = deriv(tk, wk, nk)
        k2w, k2n = deriv(tk + dt / 2, wk + dt / 2 * k1w, nk + dt / 2 * k1n)
        k3w, k3n = deriv(tk + dt / 2, wk + dt / 2 * k2w, nk + dt / 2 * k2n)
        k4w, k4n = deriv(tk + dt, wk + dt * k3w, nk + dt * k3n)
        W[k + 1] = wk + dt / 6 * (k1w + 2 * k2w + 2 * k3w + k4w)
        N[k + 1] = max(nk + dt / 6 * (k1n + 2 * k2n + 2 * k3n + k4n), 1e10)
        t[k + 1] = tk + dt

    T = _T_of(N, W)
    p_alpha = fusion_power_density(N, T, "alpha")
    p_fus = fusion_power_density(N, T, "total")
    p_brem = bremsstrahlung_density(N, T, z_eff)
    p_aux_arr = np.array([paux_f(tk) for tk in t])
    tauE_arr = np.array([tauE_f(tk) for tk in t])
    with np.errstate(divide="ignore", invalid="ignore"):
        Q = np.where(p_aux_arr > 0, p_fus / p_aux_arr, np.inf)
    return {
        "t": t, "T": T, "n": N, "W": W,
        "p_alpha": p_alpha, "p_fusion": p_fus, "p_brem": p_brem, "p_aux": p_aux_arr,
        "Q": Q, "triple": N * T * tauE_arr,
    }


# ---------------------------------------------------------------------------
# F2 — 1-D radial transport
# ---------------------------------------------------------------------------
def gaussian_deposition(rho, center, width):
    """A normalized radial deposition profile: Gaussian in rho, unit volume integral.

    Used for heating and fuelling source shapes. Normalized so that the
    volume integral  integral( profile * 2 rho drho ) = 1 (cylindrical V' ~ rho),
    i.e. multiply by a total power/particle rate to get a source density.
    """
    rho = np.asarray(rho, dtype=float)
    g = np.exp(-0.5 * ((rho - center) / width) ** 2)
    norm = np.trapezoid(g * 2.0 * rho, rho)
    return g / norm


class Transport1D:
    """1-D radial transport of a burning D-T plasma (the F2 rung).

    Evolves T(rho, t) [keV] and n(rho, t) [m^-3] on a uniform grid rho in [0, 1]
    (rho = r/a, a = minor radius). Per unit volume:

        (3/2) d(nT)/dt = (1/rho) d/drho( rho * n chi dT/drho )
                          + p_aux + p_alpha - p_brem
        dn/dt          = (1/rho) d/drho( rho * D dn/drho ) + S_fuel

    chi (heat) and D (particle) are prescribed diffusivities [m^2/s] (the F2
    "anomalous transport is an input" assumption; F3+ would compute them). The
    diffusion operators are advanced implicitly (backward Euler, tridiagonal
    solve) so conduction doesn't force a tiny timestep; the reaction/heating
    sources are explicit. Boundary: zero-gradient at the axis (rho=0), fixed
    edge values T_edge / n_edge at rho=1 (a crude pedestal).

    The temperature update is linearized per step (energy form): we step the
    energy density W = 3 n T implicitly in T at fixed n, then add sources, then
    recover T. Good enough for the burn arc; not a stiff Newton solve.
    """

    def __init__(self, a, n_grid=129, *, chi=1.0, D=0.4, z_eff=1.0,
                 T_edge=0.1, n_edge=2e19):
        self.a = float(a)
        self.rho = np.linspace(0.0, 1.0, n_grid)
        self.drho = self.rho[1] - self.rho[0]
        self.chi = chi
        self.D = D
        self.z_eff = z_eff
        self.T_edge = T_edge
        self.n_edge = n_edge
        self.t = 0.0
        # state — set by set_state()
        self.T = np.full(n_grid, T_edge)
        self.n = np.full(n_grid, n_edge)

    def set_state(self, T, n):
        """Set initial profiles (arrays over rho, or scalars)."""
        self.T = np.broadcast_to(np.asarray(T, float), self.rho.shape).astype(float).copy()
        self.n = np.broadcast_to(np.asarray(n, float), self.rho.shape).astype(float).copy()
        self.T[-1] = self.T_edge
        self.n[-1] = self.n_edge
        return self

    # --- implicit diffusion of one scalar field with a coefficient field -----
    def _diffuse(self, f, coeff, dt, f_edge):
        """Backward-Euler step of  df/dt = (1/rho) d/drho( rho * coeff * df/drho ).

        Cylindrical (rho-weighted) finite volume on the uniform grid. Returns the
        updated field with zero-gradient axis BC and Dirichlet edge value f_edge.
        `coeff` may be scalar or an array over rho (we use face-averaged values).
        """
        N = f.size
        rho, dr = self.rho, self.drho
        a2 = self.a ** 2  # rho is normalized; physical Laplacian carries 1/a^2
        coeff = np.broadcast_to(np.asarray(coeff, float), f.shape)
        # face radii and face diffusivities (i+1/2)
        rho_face = 0.5 * (rho[:-1] + rho[1:])           # length N-1
        c_face = 0.5 * (coeff[:-1] + coeff[1:])
        # cell "volume" weight ~ rho (axis cell uses rho ~ dr/4 effective)
        rho_cell = rho.copy()
        rho_cell[0] = dr / 8.0  # small positive weight at the axis singularity

        lower = np.zeros(N)
        diag = np.ones(N)
        upper = np.zeros(N)
        rhs = f.copy()

        for i in range(1, N - 1):
            w = dt / (rho_cell[i] * dr ** 2 * a2)
            aw = w * rho_face[i - 1] * c_face[i - 1]   # to i-1
            ae = w * rho_face[i] * c_face[i]           # to i+1
            lower[i] = -aw
            upper[i] = -ae
            diag[i] = 1.0 + aw + ae
        # axis: zero-gradient -> mirror, only the outward face contributes
        w0 = dt / (rho_cell[0] * dr ** 2 * a2)
        ae0 = w0 * rho_face[0] * c_face[0]
        diag[0] = 1.0 + ae0
        upper[0] = -ae0
        # edge: Dirichlet
        diag[-1] = 1.0
        lower[-1] = 0.0
        rhs[-1] = f_edge

        return _thomas(lower, diag, upper, rhs)

    def step(self, dt, *, p_aux_total=0.0, aux_profile=None,
             fuel_total=0.0, fuel_profile=None):
        """Advance the plasma by dt seconds.

        p_aux_total : auxiliary heating power per unit length-cubed? No — it is the
            volume-averaged auxiliary power density scale [W/m^3] multiplied into
            `aux_profile` (a normalized gaussian_deposition). Pass aux_profile to
            shape where the power goes; default is broad central heating.
        fuel_total : fuelling source scale [m^-3 s^-1] times `fuel_profile`.
            A deeply-peaked fuel_profile is a pellet; an edge-peaked one is a gas puff.
        """
        rho = self.rho
        if aux_profile is None:
            aux_profile = gaussian_deposition(rho, 0.0, 0.35)
        if fuel_profile is None:
            fuel_profile = gaussian_deposition(rho, 0.0, 0.3)

        n, T = self.n, self.T

        # --- particle transport (implicit) + fuelling (explicit) ---
        n_new = self._diffuse(n, self.D, dt, self.n_edge)
        n_new = n_new + dt * fuel_total * fuel_profile
        n_new[-1] = self.n_edge
        n_new = np.maximum(n_new, 1e16)

        # --- energy: implicit conduction of T, then explicit power sources ---
        # diffuse temperature (heat diffusivity chi acts on T directly here)
        T_diff = self._diffuse(T, self.chi, dt, self.T_edge)
        # power balance updates the energy density at fixed (new) density
        p_alpha = fusion_power_density(n, T, "alpha")
        p_brem = bremsstrahlung_density(n, T, self.z_eff)
        p_aux = p_aux_total * aux_profile
        dW = dt * (p_aux + p_alpha - p_brem)          # J/m^3 added this step
        # convert the conduction-updated T into energy, add source energy, recover T
        W = 3.0 * n_new * T_diff * 1e3 * E + dW
        T_new = W / (3.0 * np.maximum(n_new, 1e16) * 1e3 * E)
        T_new = np.maximum(T_new, self.T_edge)
        T_new[-1] = self.T_edge

        self.n, self.T, self.t = n_new, T_new, self.t + dt
        return self

    # --- volume-integrated diagnostics ---------------------------------------
    def _vol_avg(self, field):
        """Volume average over the circular cross-section: <f> = int f rho drho / int rho drho."""
        return np.trapezoid(field * self.rho, self.rho) / np.trapezoid(self.rho, self.rho)

    def diagnostics(self):
        """Return volume-integrated scalars: <T>, <n>, P_alpha, P_fusion, P_brem, triple, tau_E*."""
        p_alpha = self._vol_avg(fusion_power_density(self.n, self.T, "alpha"))
        p_fus = self._vol_avg(fusion_power_density(self.n, self.T, "total"))
        p_brem = self._vol_avg(bremsstrahlung_density(self.n, self.T, self.z_eff))
        W = self._vol_avg(3.0 * self.n * self.T * 1e3 * E)
        return {
            "t": self.t,
            "T_avg": self._vol_avg(self.T), "T0": self.T[0],
            "n_avg": self._vol_avg(self.n), "n0": self.n[0],
            "P_alpha": p_alpha, "P_fusion": p_fus, "P_brem": p_brem, "W": W,
        }


def _thomas(lower, diag, upper, rhs):
    """Thomas algorithm — solve a tridiagonal system in O(N). Arrays are 1-D of
    equal length; lower[0] and upper[-1] are ignored."""
    n = diag.size
    c = upper.copy()
    d = rhs.copy().astype(float)
    b = diag.copy().astype(float)
    c[0] /= b[0]
    d[0] /= b[0]
    for i in range(1, n):
        m = b[i] - lower[i] * c[i - 1]
        c[i] = upper[i] / m if i < n - 1 else 0.0
        d[i] = (d[i] - lower[i] * d[i - 1]) / m
    x = d
    for i in range(n - 2, -1, -1):
        x[i] -= c[i] * x[i + 1]
    return x
