"""Experiment 09 — Burning plasma: ignition -> steady state -> fuel injection.

The whole-discharge arc of a fusion plasma, as a *transport* problem (the right
model for the seconds-long confinement timescale — not MHD/CFD, which lives at
microseconds). We evolve the radial temperature and density profiles in time with
auxiliary heating, fusion alpha self-heating, bremsstrahlung, transport losses,
and fuelling — the toy cousin of integrated codes like ASTRA / RAPTOR.

Two modes:

  zerod   (F0) 0-D Lawson / POPCON burn: two coupled ODEs. The clearest view of
          *ignition* — alpha heating overtakes losses, the temperature runs away
          and settles at a thermally-stable burning point. The validation anchor.

  burn    (F2, default) 1-D radial transport. Three scripted phases:
            1. IGNITION       — ramp auxiliary heating until alpha-heating takes over
            2. STEADY BURN    — heating off; the plasma stays lit on its own (Q -> inf)
            3. FUEL INJECTION — a deep pellet pulse; density and fusion power respond
          Produces time traces, profile snapshots, and a poloidal cross-section
          montage (the "watch it burn" picture).

Run:
    python run.py                 # 1-D burn arc (three phases)
    python run.py --mode zerod    # 0-D ignition / Lawson demo
    python run.py --save          # write figures to ./outputs/
"""

from __future__ import annotations

import argparse

import matplotlib.pyplot as plt
import numpy as np

from plasmaplay import plotting, transport as tr

# Geometry of the toy device (large-aspect-ratio circular torus).
R0 = 3.0        # major radius [m]
A_MINOR = 1.0   # minor radius [m]
PLASMA_VOLUME = 2.0 * np.pi**2 * R0 * A_MINOR**2   # 2 pi^2 R0 a^2  [m^3]

# Confinement / transport inputs (the F2 "transport is prescribed" assumption).
# chi is tuned so the effective energy confinement tau_E ~ a^2/(5.78 chi) ~ 1 s,
# enough to cross the Lawson threshold at n ~ 1e20 and stay lit on alpha power.
CHI = 0.10      # heat diffusivity [m^2/s]
D_PART = 0.04   # particle diffusivity [m^2/s]
N_TARGET = 1.0e20   # operating density [m^-3]
TAU_P = 6.0     # particle confinement time used to set the holding fuel rate [s]

# Scripted-discharge timeline [s] and actuator levels (shared by the run and the
# cross-section re-simulation so the scenario is defined exactly once).
T_IGN_END = 4.0       # auxiliary heating ramps over [0, T_IGN_END] then switches off
T_STEADY_END = 14.0   # self-sustained burn over [T_IGN_END, T_STEADY_END]
T_PELLET = 14.0       # a deep pellet is injected at T_PELLET (for 0.2 s)
T_BURN_END = 22.0
DT = 2e-3
P_AUX_PEAK = 6.0e5    # peak auxiliary heating power density [W/m^3]
PELLET_RATE = 3.0e20  # pellet particle source [m^-3 s^-1] during the 0.2 s pulse


# ---------------------------------------------------------------------------
# F0 — 0-D ignition / Lawson
# ---------------------------------------------------------------------------
def run_zerod(save=False):
    print("\n--- 0-D burn dynamics (F0): ignition & the Lawson point ---")
    tau_E = 2.0
    fuel = N_TARGET / TAU_P            # baseline fuelling holds the density

    def paux(t):                       # ignition kick, then OFF
        return 3.0e5 if t < 4.0 else 0.0

    r = tr.burn_0d(N_TARGET, 2.0, tau_E=tau_E, p_aux=paux, fuel_rate=fuel,
                   tau_p=TAU_P, t_end=30.0, dt=1e-3)

    p_loss = r["W"] / tau_E + r["p_brem"]
    cross = np.where(r["p_alpha"] > p_loss)[0]
    t_ign = r["t"][cross[0]] if cross.size else np.nan
    print(f"  ignition (alpha overtakes losses) at t = {t_ign:.2f} s")
    print(f"  triple product there: {r['triple'][cross[0]]:.2e} keV s / m^3  (Lawson ~3e21)")
    print(f"  steady state: T = {r['T'][-1]:.1f} keV, n = {r['n'][-1]:.2e} m^-3")
    print(f"  power balance check  P_alpha/(P_loss+P_brem) = "
          f"{r['p_alpha'][-1] / p_loss[-1]:.4f}  (should be ~1)")
    print(f"  self-sustained after heating off: {r['p_aux'][-1] == 0 and r['T'][-1] > 10}")

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    ax[0].plot(r["t"], r["T"], color="crimson")
    ax[0].axvspan(0, 4, color="gold", alpha=0.2, label="aux heating on")
    ax[0].axvline(t_ign, ls="--", color="k", lw=0.8)
    ax[0].annotate("ignition", (t_ign, r["T"][cross[0]]), textcoords="offset points",
                   xytext=(8, -4))
    ax[0].set(xlabel="t [s]", ylabel="T [keV]", title="Temperature: ignition then burn")
    ax[0].legend()

    ax[1].plot(r["t"], r["p_alpha"] * PLASMA_VOLUME / 1e6, label="alpha (self) heat")
    ax[1].plot(r["t"], p_loss * PLASMA_VOLUME / 1e6, label="losses (transport+brem)")
    ax[1].plot(r["t"], r["p_aux"] * PLASMA_VOLUME / 1e6, label="aux heating")
    ax[1].set(xlabel="t [s]", ylabel="power [MW]", title="Power balance")
    ax[1].legend()
    fig.tight_layout()
    _finish(fig, save, "burn_0d_ignition.png")


# ---------------------------------------------------------------------------
# F2 — 1-D burn arc
# ---------------------------------------------------------------------------
def _make_sim(n_grid):
    """Build the simulator and the fixed deposition shapes for the scenario."""
    sim = tr.Transport1D(A_MINOR, n_grid=n_grid, chi=CHI, D=D_PART,
                         T_edge=0.1, n_edge=2e19)
    sim.set_state(T=2.0, n=N_TARGET)
    shapes = {
        "aux": tr.gaussian_deposition(sim.rho, 0.0, 0.35),     # central heating
        "hold": tr.gaussian_deposition(sim.rho, 0.0, 0.40),    # broad gas-puff fuelling
        "pellet": tr.gaussian_deposition(sim.rho, 0.35, 0.12), # deep, localized pellet
    }
    return sim, shapes


def _advance(sim, shapes):
    """Apply one scripted scenario step (the three phases live here, once)."""
    t = sim.t
    hold_rate = N_TARGET / TAU_P                       # holds the density
    p_aux = P_AUX_PEAK * (0.3 + 0.7 * t / T_IGN_END) if t < T_IGN_END else 0.0
    fuel_total, fuel_profile = hold_rate, shapes["hold"]
    if T_PELLET <= t < T_PELLET + 0.2:                 # 0.2 s pellet pulse
        fuel_total, fuel_profile = hold_rate + PELLET_RATE, shapes["pellet"]
    sim.step(DT, p_aux_total=p_aux, aux_profile=shapes["aux"],
             fuel_total=fuel_total, fuel_profile=fuel_profile)
    return p_aux


def run_burn(save=False, n_grid=129):
    print("\n--- 1-D burning plasma (F2): ignition -> steady -> fuel injection ---")
    sim, shapes = _make_sim(n_grid)

    snaps = {}
    hist = {k: [] for k in ("t", "T0", "Tavg", "n0", "navg", "Pfus", "Palpha", "Paux")}

    for _ in range(int(round(T_BURN_END / DT))):
        t = sim.t
        p_aux = _advance(sim, shapes)

        d = sim.diagnostics()
        hist["t"].append(t)
        hist["T0"].append(d["T0"]); hist["Tavg"].append(d["T_avg"])
        hist["n0"].append(d["n0"]); hist["navg"].append(d["n_avg"])
        hist["Pfus"].append(d["P_fusion"] * PLASMA_VOLUME / 1e6)   # MW
        hist["Palpha"].append(d["P_alpha"] * PLASMA_VOLUME / 1e6)
        hist["Paux"].append(p_aux * PLASMA_VOLUME / 1e6)

        # snapshot profiles at representative moments (ordered by time)
        for label, tt in (("heating", T_IGN_END * 0.6), ("ignited", T_IGN_END),
                          ("steady", T_STEADY_END - 0.1), ("post-pellet", T_PELLET + 1.5)):
            if label not in snaps and sim.t >= tt:
                snaps[label] = (sim.rho.copy(), sim.T.copy(), sim.n.copy(), tt)

    hist = {k: np.asarray(v) for k, v in hist.items()}

    i_steady = np.argmin(np.abs(hist["t"] - (T_STEADY_END - 0.2)))
    print(f"  ignition phase peak T0 = {hist['T0'][:int(T_IGN_END / DT)].max():.1f} keV")
    print(f"  steady burn: T0 = {hist['T0'][i_steady]:.1f} keV, "
          f"P_fusion = {hist['Pfus'][i_steady]:.1f} MW (aux off -> Q = inf)")
    dn = (hist['navg'][np.argmin(np.abs(hist['t'] - (T_PELLET + 0.3)))]
          - hist['navg'][np.argmin(np.abs(hist['t'] - (T_PELLET - 0.1)))])
    print(f"  pellet at t={T_PELLET}s raised <n> by {dn:.2e} m^-3; "
          f"P_fusion settles to {hist['Pfus'][-1]:.1f} MW")

    _plot_burn_traces(hist, T_IGN_END, T_STEADY_END, T_PELLET, save)
    _plot_burn_profiles(snaps, save)
    _plot_cross_sections(snaps, save, n_grid)


def _plot_burn_traces(h, t1, t2, tp, save):
    fig, ax = plt.subplots(3, 1, figsize=(9, 9), sharex=True)
    for a in ax:
        a.axvspan(0, t1, color="gold", alpha=0.18)
        a.axvspan(t1, t2, color="seagreen", alpha=0.10)
        a.axvline(tp, color="purple", ls="--", lw=1.0)
    ax[0].plot(h["t"], h["T0"], label="T core")
    ax[0].plot(h["t"], h["Tavg"], label="<T>", ls="--")
    ax[0].set(ylabel="T [keV]", title="Burning-plasma discharge: ignition | steady | pellet")
    ax[0].legend(loc="upper right")
    ax[1].plot(h["t"], h["n0"], label="n core", color="navy")
    ax[1].plot(h["t"], h["navg"], label="<n>", ls="--", color="navy")
    ax[1].set(ylabel="n [m$^{-3}$]"); ax[1].legend(loc="upper right")
    ax[2].plot(h["t"], h["Pfus"], label="P fusion", color="crimson")
    ax[2].plot(h["t"], h["Palpha"], label="P alpha (self-heat)", color="darkorange")
    ax[2].plot(h["t"], h["Paux"], label="P aux", color="gray")
    ax[2].set(xlabel="t [s]", ylabel="power [MW]"); ax[2].legend(loc="upper right")
    # phase labels
    ax[0].text(t1/2, ax[0].get_ylim()[1]*0.9, "IGNITION", ha="center", fontsize=9)
    ax[0].text((t1+t2)/2, ax[0].get_ylim()[1]*0.9, "STEADY BURN", ha="center", fontsize=9)
    ax[0].text((t2+h["t"][-1])/2, ax[0].get_ylim()[1]*0.9, "+ FUEL", ha="center", fontsize=9)
    fig.tight_layout()
    _finish(fig, save, "burn_1d_traces.png")


def _plot_burn_profiles(snaps, save):
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    order = ["heating", "ignited", "steady", "post-pellet"]
    for label in order:
        if label not in snaps:
            continue
        rho, T, n, tt = snaps[label]
        ax[0].plot(rho, T, label=f"{label} (t={tt:.1f}s)")
        ax[1].plot(rho, n, label=f"{label} (t={tt:.1f}s)")
    ax[0].set(xlabel=r"$\rho = r/a$", ylabel="T [keV]", title="Temperature profiles")
    ax[1].set(xlabel=r"$\rho = r/a$", ylabel="n [m$^{-3}$]", title="Density profiles")
    for a in ax:
        a.legend()
    fig.tight_layout()
    _finish(fig, save, "burn_1d_profiles.png")


def _plot_cross_sections(snaps, save, n_grid):
    """Map each captured T(rho) snapshot onto a poloidal disk — the 'watch it burn'
    view. The model is 1-D (flux functions), so each surface rho = const is an
    isotherm; we revolve the profile to fill the circular cross-section."""
    order = [k for k in ("heating", "ignited", "steady", "post-pellet") if k in snaps]
    theta = np.linspace(0, 2 * np.pi, 160)
    vmax = max(snaps[k][1].max() for k in order)

    fig, axes = plt.subplots(1, len(order), figsize=(4 * len(order), 4.4),
                             subplot_kw={"aspect": "equal"})
    if len(order) == 1:
        axes = [axes]
    for ax, label in zip(axes, order):
        rho, T, _n, tt = snaps[label]
        RR, TT = np.meshgrid(rho, theta)
        X, Y = RR * np.cos(TT), RR * np.sin(TT)
        field = np.broadcast_to(T, RR.shape)
        pc = ax.contourf(X, Y, field, levels=40, cmap="inferno", vmin=0, vmax=vmax)
        ax.set_title(f"{label}\nt = {tt:.1f} s, T0 = {T[0]:.0f} keV", fontsize=10)
        ax.set_xticks([]); ax.set_yticks([])
    fig.colorbar(pc, ax=axes, label="T [keV]", shrink=0.8)
    fig.suptitle("Poloidal cross-section: the plasma igniting and burning", y=1.02)
    _finish(fig, save, "burn_1d_cross_sections.png")


def _finish(fig, save, name):
    if save:
        out = plotting.ensure_outputs_dir(__file__) / name
        fig.savefig(out, dpi=130, bbox_inches="tight")
        print(f"  saved {out}")
    else:
        plt.show()
    plt.close(fig)


def main(mode="burn", save=False, n_grid=129):
    print("=" * 64)
    print("Experiment 09 — burning plasma (ignition -> steady -> fuelling)")
    print("=" * 64)
    if mode == "zerod":
        run_zerod(save=save)
    else:
        run_burn(save=save, n_grid=n_grid)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--mode", choices=["burn", "zerod"], default="burn")
    p.add_argument("--save", action="store_true", help="write figures to ./outputs/")
    p.add_argument("--n-grid", type=int, default=129)
    args = p.parse_args()
    main(mode=args.mode, save=args.save, n_grid=args.n_grid)
