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
from .constants import m_e as M_E, m_p as M_P  # electron / proton mass [kg]

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
def reaction_rate_dt(n_dt, T_keV):
    """D-T reactions per m^3 per second for a 50:50 mix of fuel-ion density `n_dt`.

    n_D = n_T = n_dt/2, so the reaction-rate density is (n_dt/2)^2 <sigma v>. This
    is the single quantity that drives *both* the alpha heating and the helium-ash
    production / fuel burnup (one He nucleus + 2 fuel ions consumed per reaction).
    """
    n = np.asarray(n_dt, dtype=float)
    return 0.25 * n ** 2 * reactivity_dt(T_keV)


def fusion_power_density(n, T_keV, which="alpha"):
    """Fusion power density (W/m^3) for a 50:50 D-T plasma of fuel-ion density n.

    n_D = n_T = n/2, so the reaction-rate density is (n/2)^2 <sigma v>. `which`
    selects which slice of the 17.6 MeV is counted:
      "alpha"   -> 3.5 MeV  (the part that self-heats the plasma)
      "neutron" -> 14.1 MeV (the part that leaves — fusion power to the blanket)
      "total"   -> 17.6 MeV
    """
    energy = {"alpha": E_ALPHA_MEV, "neutron": E_NEUTRON_MEV, "total": E_FUSION_MEV}[which]
    return reaction_rate_dt(n, T_keV) * energy * _MEV_J


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
# F1 — 0-D with helium ash, fuel dilution/burnup, and a soft beta-limit
# ---------------------------------------------------------------------------
MU_0 = 4.0e-7 * np.pi   # vacuum permeability [H/m]


def beta_thermal(w_density, B):
    """Thermal beta = 2 mu0 p / B^2, the ratio of plasma to magnetic pressure.

    For an ideal gas W = (3/2) p, so the thermal pressure is p = (2/3) W. `B` is the
    on-axis toroidal field [T]. Returned as a fraction (multiply by 100 for percent).
    """
    p = (2.0 / 3.0) * np.asarray(w_density, dtype=float)
    return 2.0 * MU_0 * p / B ** 2


def troyon_limit(beta_N, Ip_MA, a, B):
    """Troyon beta limit as a *fraction*: beta_max[%] = beta_N * Ip[MA] / (a[m] B[T]).

    beta_N is the normalized-beta coefficient in percent-meter-tesla-per-megaamp
    (the Troyon value is ~2.8; machines push ~3-4). Divided by 100 to return a
    fraction so it compares directly with `beta_thermal`.
    """
    return beta_N * Ip_MA / (a * B) / 100.0


def _w_of_ash(n_dt, n_he, T_keV):
    """Energy density W = (3/2)(n_e + n_i) T with ash present.

    Quasineutrality with doubly-charged helium: n_e = n_dt + 2 n_he; the total ion
    density is n_i = n_dt + n_he. So W = (3/2)(2 n_dt + 3 n_he) T (T in keV).
    """
    return 1.5 * (2.0 * n_dt + 3.0 * n_he) * T_keV * 1e3 * E


def _T_of_ash(n_dt, n_he, w):
    """Invert `_w_of_ash` for T in keV (guards a vanishing heat-capacity)."""
    denom = 1.5 * (2.0 * np.asarray(n_dt, float) + 3.0 * np.asarray(n_he, float)) * 1e3 * E
    return w / np.maximum(denom, 1e-12)


def z_eff_with_ash(n_dt, n_he):
    """Z_eff = sum(n_i Z_i^2) / n_e for a D-T + He(2+) plasma.

    Z_eff = (n_dt * 1^2 + n_he * 2^2) / (n_dt + 2 n_he). Ash raises Z_eff, which
    raises the bremsstrahlung loss — a second penalty for letting ash accumulate.
    """
    n_dt = np.asarray(n_dt, float)
    n_he = np.asarray(n_he, float)
    n_e = n_dt + 2.0 * n_he
    return (n_dt + 4.0 * n_he) / np.maximum(n_e, 1e10)


def burn_0d_ash(n_dt0, T0, *, tau_E, p_aux, B, tau_p=None, tau_he=None,
                fuel_rate=0.0, beta_limit=None, beta_stiffness=80.0,
                t_end=30.0, dt=1e-3):
    """0-D burn with He ash, fuel dilution + burnup, and a soft beta-limit (F1).

    Three coupled ODEs for the fuel-ion density n_DT, the helium-ash density n_He,
    and the energy density W = (3/2)(n_e + n_i) T:

        dn_DT/dt = S_fuel - n_DT/tau_p - 2 R_fus     (2 fuel ions consumed / reaction)
        dn_He/dt =          R_fus      - n_He/tau_he (ash born / reaction, pumped on tau_he)
        dW/dt    = P_aux + P_alpha - P_brem - W/tau_E_eff

    with R_fus = reaction_rate_dt(n_DT, T) and quasineutrality n_e = n_DT + 2 n_He.
    Bremsstrahlung uses the ash-raised Z_eff. The **soft beta-limit** degrades the
    energy confinement as the pressure approaches the limit,

        tau_E_eff = tau_E / (1 + beta_stiffness * max(0, beta/beta_limit - 1)^2),

    so the operating point is pinned near beta_limit instead of running away to the
    ~80 keV point of the limit-free model — landing it in the real ~10-25 keV band.
    The penalty is **one-sided**: below the limit confinement is untouched (ignition
    proceeds normally), and it bites only once beta exceeds beta_limit, acting as a
    soft stability wall the pressure cannot climb past.

    Why this rung matters: real machines must *pump* the helium ash and *refuel* the
    burned D-T continuously, and they operate beta-limited. Turn `tau_he` short and
    `beta_limit=None` to recover something close to the F0 picture.

    Parameters mirror `burn_0d`, plus: `B` on-axis field [T] (needed for beta),
    `tau_he` ash particle-confinement time [s] (default 5*tau_E), `beta_limit`
    fraction (None disables it), `beta_stiffness` the soft-limit exponent.

    Returns a dict of arrays: t, T (keV), n_DT, n_He, n_e, f_He (= n_He/n_e), z_eff,
    beta, tau_E_eff, p_alpha, p_fusion, p_brem, p_aux, Q (= P_fus/P_aux), triple.
    """
    tau_p = 3.0 * tau_E if tau_p is None else tau_p
    tau_he = 5.0 * tau_E if tau_he is None else tau_he

    def as_f(x):
        return x if callable(x) else (lambda t, _x=x: _x)

    tauE_f, paux_f, fuel_f = as_f(tau_E), as_f(p_aux), as_f(fuel_rate)

    def tauE_eff(t, w):
        tE = tauE_f(t)
        if beta_limit is None:
            return tE
        excess = max(0.0, beta_thermal(w, B) / beta_limit - 1.0)
        return tE / (1.0 + beta_stiffness * excess ** 2)

    def deriv(t, n_dt, n_he, w):
        T = _T_of_ash(n_dt, n_he, w)
        R = reaction_rate_dt(n_dt, T)
        n_e = n_dt + 2.0 * n_he
        z = float(z_eff_with_ash(n_dt, n_he))
        p_a = R * E_ALPHA_MEV * _MEV_J
        p_b = bremsstrahlung_density(n_e, T, z)
        dw = paux_f(t) + p_a - p_b - w / tauE_eff(t, w)
        d_ndt = fuel_f(t) - n_dt / tau_p - 2.0 * R
        d_nhe = R - n_he / tau_he
        return d_ndt, d_nhe, dw

    nstep = int(round(t_end / dt))
    t = np.empty(nstep + 1)
    NDT = np.empty(nstep + 1)
    NHE = np.empty(nstep + 1)
    W = np.empty(nstep + 1)
    NDT[0], NHE[0], W[0], t[0] = float(n_dt0), 0.0, _w_of_ash(n_dt0, 0.0, T0), 0.0

    for k in range(nstep):
        tk, a0, b0, w0 = t[k], NDT[k], NHE[k], W[k]
        k1 = deriv(tk, a0, b0, w0)
        k2 = deriv(tk + dt / 2, a0 + dt / 2 * k1[0], b0 + dt / 2 * k1[1], w0 + dt / 2 * k1[2])
        k3 = deriv(tk + dt / 2, a0 + dt / 2 * k2[0], b0 + dt / 2 * k2[1], w0 + dt / 2 * k2[2])
        k4 = deriv(tk + dt, a0 + dt * k3[0], b0 + dt * k3[1], w0 + dt * k3[2])
        NDT[k + 1] = max(a0 + dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]), 1e10)
        NHE[k + 1] = max(b0 + dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]), 0.0)
        W[k + 1] = max(w0 + dt / 6 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]), 1e-6)
        t[k + 1] = tk + dt

    T = _T_of_ash(NDT, NHE, W)
    n_e = NDT + 2.0 * NHE
    R = reaction_rate_dt(NDT, T)
    p_alpha = R * E_ALPHA_MEV * _MEV_J
    p_fus = R * E_FUSION_MEV * _MEV_J
    z = z_eff_with_ash(NDT, NHE)
    p_brem = bremsstrahlung_density(n_e, T, z)
    p_aux_arr = np.array([paux_f(tk) for tk in t])
    beta = beta_thermal(W, B)
    tauE_arr = np.array([tauE_eff(tk, wk) for tk, wk in zip(t, W)])
    with np.errstate(divide="ignore", invalid="ignore"):
        Q = np.where(p_aux_arr > 0, p_fus / p_aux_arr, np.inf)
    return {
        "t": t, "T": T, "n_DT": NDT, "n_He": NHE, "n_e": n_e,
        "f_He": NHE / np.maximum(n_e, 1e10), "z_eff": z, "beta": beta,
        "tau_E_eff": tauE_arr, "p_alpha": p_alpha, "p_fusion": p_fus,
        "p_brem": p_brem, "p_aux": p_aux_arr, "Q": Q, "triple": n_e * T * tauE_arr,
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


# ---------------------------------------------------------------------------
# F2.5 — two temperatures (Te, Ti): Spitzer equipartition + a heating mix
# ---------------------------------------------------------------------------
# Single-temperature transport hides three facts: fusion is an *ion* reaction
# (the burn rate and alpha power follow T_i), bremsstrahlung is an *electron*
# loss (it follows T_e), and the two species exchange energy only collisionally,
# on the Spitzer equipartition time. When the heating is lopsided — neutral beams
# and RF deposit on ions, fusion alphas slow down mostly on electrons — the two
# temperatures separate, and a beam-heated plasma runs T_i > T_e. These helpers
# add that physics; the coupling term Q_Delta is the validation anchor.
def coulomb_logarithm(n_e, T_e_keV):
    """Electron-ion Coulomb logarithm lnLambda (NRL Plasma Formulary, T_e>10 eV).

    lnLambda = 24 - ln( sqrt(n_e[cm^-3]) / T_e[eV] ), with n_e in m^-3 and T_e in
    keV on input. It is weakly varying (~15-20 for fusion plasmas) and sets the
    overall scale of every Coulomb collision rate.
    """
    n_cm3 = np.asarray(n_e, dtype=float) * 1e-6
    T_eV = np.asarray(T_e_keV, dtype=float) * 1e3
    return 24.0 - np.log(np.sqrt(n_cm3) / T_eV)


def collision_frequency_ei(n_e, T_e_keV, z_eff=1.0, coulomb_log=None):
    """Spitzer electron collision frequency nu_e [s^-1] (NRL Plasma Formulary).

        nu_e = 2.91e-6 * Z * n_e[cm^-3] * lnLambda * T_e[eV]^(-3/2)

    This is the electron-ion momentum collision rate (~ 1/tau_e), the clock that
    times collisional energy exchange between the species. Vectorized over inputs.
    """
    n_cm3 = np.asarray(n_e, dtype=float) * 1e-6
    T_eV = np.asarray(T_e_keV, dtype=float) * 1e3
    lnL = coulomb_logarithm(n_e, T_e_keV) if coulomb_log is None else coulomb_log
    return 2.91e-6 * z_eff * n_cm3 * lnL * T_eV ** (-1.5)


def equipartition_power(n_e, T_e_keV, T_i_keV, *, z_eff=1.0, mu_i=2.5,
                        coulomb_log=None):
    """Collisional electron->ion energy-exchange power density [W/m^3] (Braginskii).

        Q_Delta = 3 (m_e/m_i) n_e nu_ei k_B (T_e - T_i)

    Positive when the electrons are hotter (they heat the ions); it flips sign and
    vanishes at T_e = T_i. m_i = mu_i * m_p (mu_i ~ 2.5 amu for a 50:50 D-T mix).
    This is the only term that drags T_e and T_i toward each other.
    """
    n_e = np.asarray(n_e, dtype=float)
    nu = collision_frequency_ei(n_e, T_e_keV, z_eff, coulomb_log)
    m_ratio = M_E / (mu_i * M_P)
    dT = np.asarray(T_e_keV, dtype=float) - np.asarray(T_i_keV, dtype=float)
    return 3.0 * m_ratio * n_e * nu * dT * 1e3 * E


def equipartition_time(n_e, T_e_keV, *, z_eff=1.0, mu_i=2.5, coulomb_log=None):
    """e-folding time [s] of the temperature difference (T_e - T_i) at n_e = n_i.

    Substituting Q_Delta into the two energy equations (each species has heat
    capacity 3/2 n) gives d(T_e - T_i)/dt = -(T_e - T_i)/tau_eq with

        tau_eq = 1 / ( 4 (m_e/m_i) nu_ei ).

    It scales as T_e^{3/2}/n_e: hot, thin plasmas equilibrate slowly, which is
    exactly why they can sustain T_i != T_e long enough to matter.
    """
    nu = collision_frequency_ei(n_e, T_e_keV, z_eff, coulomb_log)
    m_ratio = M_E / (mu_i * M_P)
    return 1.0 / (4.0 * m_ratio * nu)


def two_temperature_relax_0d(n_e, Te0, Ti0, *, z_eff=1.0, mu_i=2.5,
                             t_end=1.0, dt=1e-4):
    """0-D collisional relaxation of (T_e, T_i) toward a common temperature.

    No heating, no transport — only the equipartition exchange Q_Delta. At fixed
    n_e = n_i the total energy (3/2) n (T_e + T_i) is conserved and T_e, T_i relax
    to their mean on ~tau_eq. This is the clean validation anchor for the coupling:
    the measured difference-decay rate must match `equipartition_time`. Forward
    Euler (the exchange is gentle for dt << tau_eq). Returns dict: t, T_e, T_i [keV].
    """
    nstep = int(round(t_end / dt))
    t = np.empty(nstep + 1)
    Te = np.empty(nstep + 1)
    Ti = np.empty(nstep + 1)
    Te[0], Ti[0], t[0] = float(Te0), float(Ti0), 0.0
    cap = 1.5 * float(n_e) * 1e3 * E      # heat capacity per species [J/m^3/keV]
    for k in range(nstep):
        q = equipartition_power(n_e, Te[k], Ti[k], z_eff=z_eff, mu_i=mu_i)
        Te[k + 1] = Te[k] - dt * q / cap
        Ti[k + 1] = Ti[k] + dt * q / cap
        t[k + 1] = t[k] + dt
    return {"t": t, "T_e": Te, "T_i": Ti}


class TwoTempTransport1D(Transport1D):
    """1-D transport with separate electron and ion temperatures (the F2.5 rung).

    Evolves T_e(rho,t), T_i(rho,t) and n(rho,t). Each temperature channel diffuses
    with its own heat diffusivity, the two exchange energy collisionally, and they
    receive different heating:

        (3/2) d(n T_e)/dt = div(chi_e grad T_e) + p_aux_e + f_ae p_alpha - p_brem - Q
        (3/2) d(n T_i)/dt = div(chi_i grad T_i) + p_aux_i + (1-f_ae) p_alpha       + Q

    with Q = equipartition_power(n, T_e, T_i). The split that single-T transport
    cannot represent: p_alpha and the burn use the *ion* temperature (fusion is an
    ion reaction); p_brem uses the *electron* temperature; fusion alphas slow down
    predominantly on electrons (f_ae ~ 0.8-0.9), while neutral-beam / RF-ion heating
    lands on the ions — so a beam-heated burning plasma settles at T_i > T_e until
    equipartition closes the gap. Mirrors `Transport1D`'s implicit-diffusion /
    explicit-source split and reuses its tridiagonal solver.
    """

    def __init__(self, a, n_grid=129, *, chi_e=1.0, chi_i=0.5, D=0.4, z_eff=1.0,
                 mu_i=2.5, Te_edge=0.1, Ti_edge=0.1, n_edge=2e19):
        super().__init__(a, n_grid, chi=chi_e, D=D, z_eff=z_eff,
                         T_edge=Te_edge, n_edge=n_edge)
        self.chi_e = chi_e
        self.chi_i = chi_i
        self.mu_i = mu_i
        self.Te_edge = Te_edge
        self.Ti_edge = Ti_edge
        self.Te = np.full(n_grid, Te_edge)
        self.Ti = np.full(n_grid, Ti_edge)

    def set_state(self, Te, Ti, n):
        """Set initial electron/ion temperature and density profiles (arrays or scalars)."""
        self.Te = np.broadcast_to(np.asarray(Te, float), self.rho.shape).astype(float).copy()
        self.Ti = np.broadcast_to(np.asarray(Ti, float), self.rho.shape).astype(float).copy()
        self.n = np.broadcast_to(np.asarray(n, float), self.rho.shape).astype(float).copy()
        self.Te[-1] = self.Te_edge
        self.Ti[-1] = self.Ti_edge
        self.n[-1] = self.n_edge
        self.T = self.Te                    # keep the parent single-T attribute sane
        return self

    def step(self, dt, *, p_aux_e_total=0.0, p_aux_i_total=0.0, aux_e_profile=None,
             aux_i_profile=None, frac_alpha_e=0.85, fuel_total=0.0, fuel_profile=None):
        """Advance T_e, T_i, n by dt seconds.

        p_aux_e_total / p_aux_i_total : electron- and ion-channel heating power-density
            scales [W/m^3] multiplied into their (normalized) deposition profiles —
            e.g. RF/ohmic to electrons, neutral-beam to ions.
        frac_alpha_e : fraction of the 3.5 MeV alpha power deposited on electrons
            (the rest on ions); ~0.85 for a hot D-T plasma.
        fuel_total / fuel_profile : as in `Transport1D.step`.
        """
        rho = self.rho
        if aux_e_profile is None:
            aux_e_profile = gaussian_deposition(rho, 0.0, 0.35)
        if aux_i_profile is None:
            aux_i_profile = gaussian_deposition(rho, 0.0, 0.35)
        if fuel_profile is None:
            fuel_profile = gaussian_deposition(rho, 0.0, 0.3)

        n, Te, Ti = self.n, self.Te, self.Ti

        # particle transport (implicit) + fuelling (explicit)
        n_new = self._diffuse(n, self.D, dt, self.n_edge)
        n_new = n_new + dt * fuel_total * fuel_profile
        n_new[-1] = self.n_edge
        n_new = np.maximum(n_new, 1e16)

        # conduction (implicit) of each temperature channel with its own chi
        Te_diff = self._diffuse(Te, self.chi_e, dt, self.Te_edge)
        Ti_diff = self._diffuse(Ti, self.chi_i, dt, self.Ti_edge)

        # explicit sources/sinks: fusion follows T_i, brem follows T_e, Q couples them
        p_alpha = fusion_power_density(n, Ti, "alpha")
        p_brem = bremsstrahlung_density(n, Te, self.z_eff)
        q_ei = equipartition_power(n, Te, Ti, z_eff=self.z_eff, mu_i=self.mu_i)
        p_e = p_aux_e_total * aux_e_profile + frac_alpha_e * p_alpha - p_brem - q_ei
        p_i = p_aux_i_total * aux_i_profile + (1.0 - frac_alpha_e) * p_alpha + q_ei

        cap = 1.5 * n_new * 1e3 * E          # heat capacity per species [J/m^3/keV]
        Te_new = Te_diff + dt * p_e / cap
        Ti_new = Ti_diff + dt * p_i / cap
        Te_new = np.maximum(Te_new, self.Te_edge)
        Ti_new = np.maximum(Ti_new, self.Ti_edge)
        Te_new[-1] = self.Te_edge
        Ti_new[-1] = self.Ti_edge

        self.n, self.Te, self.Ti, self.T = n_new, Te_new, Ti_new, Te_new
        self.t += dt
        return self

    def diagnostics(self):
        """Volume-integrated two-temperature scalars: <Te>, <Ti>, on-axis values,
        <n>, and the (ion-temperature) fusion / (electron) brem powers."""
        return {
            "t": self.t,
            "Te_avg": self._vol_avg(self.Te), "Te0": self.Te[0],
            "Ti_avg": self._vol_avg(self.Ti), "Ti0": self.Ti[0],
            "n_avg": self._vol_avg(self.n), "n0": self.n[0],
            "P_alpha": self._vol_avg(fusion_power_density(self.n, self.Ti, "alpha")),
            "P_fusion": self._vol_avg(fusion_power_density(self.n, self.Ti, "total")),
            "P_brem": self._vol_avg(bremsstrahlung_density(self.n, self.Te, self.z_eff)),
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
