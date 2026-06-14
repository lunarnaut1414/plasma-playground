"""Animation helpers — turn a time-series of fields into a `.gif` for the showcase.

The experiments in this playground produce time-evolving radial profiles (T(ρ,t),
n(ρ,t), ψ(r,θ,t), …). These helpers render those into animations:

  * `animate_profiles`       — line plot(s) over r, animated in time.
  * `animate_cross_section`  — a radial profile revolved into a poloidal disk
                               heatmap, animated in time (the "watch it burn" view).
  * `animate_torus_3d`       — a rotating 3-D torus surface colored by a field.

Design: the **frame-data construction is separated from the gif writing** so the
physics can be unit-tested without rendering bytes. `make_frames` samples a field
function on a time grid into a plain array; the `animate_*` functions take arrays
and write a gif via matplotlib's `PillowWriter` (pure Python — no external binary).

Gifs go to (gitignored) `outputs/`. Keep them modest: ~60–120 frames, dpi ≈ 90.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg", force=False)  # headless-safe; respects an already-set backend
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter  # noqa: E402


def make_frames(field_fn, times):
    """Sample `field_fn(t)` on `times` into a stacked array of shape (n_times, ...).

    Pure and testable: no plotting. `field_fn` returns an array (any shape) for a
    scalar time `t`; the results are stacked along a new leading time axis.
    """
    times = np.asarray(times, dtype=float)
    return np.stack([np.asarray(field_fn(t), dtype=float) for t in times], axis=0)


def torus_surface(R0, a, n_u=80, n_v=40):
    """(X, Y, Z) of a torus surface: tube radius `a`, major radius `R0`.

    u = toroidal angle (around the machine), v = poloidal angle (around the tube).
    Returned arrays have shape (n_u, n_v). Every point satisfies the implicit torus
    equation (√(X²+Y²) − R0)² + Z² = a², which is what the unit test checks.
    """
    u = np.linspace(0.0, 2.0 * np.pi, n_u)
    v = np.linspace(0.0, 2.0 * np.pi, n_v)
    U, V = np.meshgrid(u, v, indexing="ij")
    X = (R0 + a * np.cos(V)) * np.cos(U)
    Y = (R0 + a * np.cos(V)) * np.sin(U)
    Z = a * np.sin(V)
    return X, Y, Z


def _prepare(path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def animate_profiles(x, frames, times=None, *, path, labels=None, xlabel="r/a",
                     ylabel="", title="", fps=20, dpi=90, ylim=None):
    """Animate one or several radial line profiles over time, save to `path` (gif).

    `frames` is (n_t, n_x) for a single curve, or (n_t, n_series, n_x) for several
    curves drawn together (e.g. Te and Ti). Returns the saved Path.
    """
    frames = np.asarray(frames, dtype=float)
    x = np.asarray(x, dtype=float)
    if frames.ndim == 2:
        frames = frames[:, None, :]            # -> (n_t, 1, n_x)
    n_t, n_series, _ = frames.shape
    labels = labels or [None] * n_series
    if ylim is None:
        ylim = (min(0.0, frames.min()), frames.max() * 1.08 + 1e-30)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    lines = [ax.plot([], [], label=labels[s])[0] for s in range(n_series)]
    ax.set(xlim=(x.min(), x.max()), ylim=ylim, xlabel=xlabel, ylabel=ylabel)
    if any(lbl is not None for lbl in labels):
        ax.legend(loc="upper right")
    ttl = ax.set_title(title)

    def draw(i):
        for s, ln in enumerate(lines):
            ln.set_data(x, frames[i, s])
        if times is not None:
            ttl.set_text(f"{title}   t = {times[i]:.2f} s")
        return (*lines, ttl)

    anim = FuncAnimation(fig, draw, frames=n_t, blit=False)
    out = _prepare(path)
    anim.save(str(out), writer=PillowWriter(fps=fps), dpi=dpi)
    plt.close(fig)
    return out


def animate_phase_track(x, y, times=None, *, path, color=None, xlabel="", ylabel="",
                        title="", clabel="", cmap="viridis", xlim=None, ylim=None,
                        fps=20, dpi=90, logx=False):
    """Animate a 2-D trajectory (x[i], y[i]) as a growing, fading tail with a head.

    The reusable "phase-space movie" used for burn (n,T) ignition tracks: a point
    sweeps the plane over time, leaving a tail; if `color` (one scalar per frame) is
    given the head is colored by it with a colorbar (e.g. ash fraction). Returns the
    saved Path. Pure rendering — the *physics* lives in the (x,y,color) arrays.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n_t = x.size
    xlim = xlim or (x.min(), x.max() * 1.05 + 1e-30)
    ylim = ylim or (min(0.0, y.min()), y.max() * 1.08 + 1e-30)
    norm = None
    if color is not None:
        color = np.asarray(color, dtype=float)
        norm = plt.Normalize(vmin=float(color.min()), vmax=float(color.max()) + 1e-30)

    fig, ax = plt.subplots(figsize=(6.5, 5))
    ax.set(xlim=xlim, ylim=ylim, xlabel=xlabel, ylabel=ylabel)
    if logx:
        ax.set_xscale("log")
    (tail,) = ax.plot([], [], color="0.6", lw=1.2, alpha=0.8)
    head = ax.scatter([x[0]], [y[0]], s=80, zorder=3,
                      c=([color[0]] if color is not None else "crimson"),
                      cmap=cmap, norm=norm, edgecolors="k", linewidths=0.6)
    if color is not None:
        fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), ax=ax, label=clabel,
                     shrink=0.85)
    ttl = ax.set_title(title)

    def draw(i):
        tail.set_data(x[: i + 1], y[: i + 1])
        head.set_offsets([[x[i], y[i]]])
        if color is not None:
            head.set_array(np.array([color[i]]))
        if times is not None:
            ttl.set_text(f"{title}   t = {times[i]:.1f} s")
        return tail, head, ttl

    anim = FuncAnimation(fig, draw, frames=n_t, blit=False)
    out = _prepare(path)
    anim.save(str(out), writer=PillowWriter(fps=fps), dpi=dpi)
    plt.close(fig)
    return out


def animate_cross_section(rho, frames, times=None, *, path, title="",
                          cmap="inferno", clabel="", fps=20, dpi=90,
                          vmin=0.0, vmax=None, n_theta=160):
    """Revolve a radial profile into a poloidal disk and animate it over time.

    `frames` is (n_t, n_rho). Each frame is broadcast over poloidal angle (the
    model is 1-D in ρ, so ρ = const is an isotherm) and drawn as a filled disk.
    """
    frames = np.asarray(frames, dtype=float)
    rho = np.asarray(rho, dtype=float)
    vmax = frames.max() if vmax is None else vmax
    theta = np.linspace(0.0, 2.0 * np.pi, n_theta)
    RR, TT = np.meshgrid(rho, theta)
    X, Y = RR * np.cos(TT), RR * np.sin(TT)
    levels = np.linspace(vmin, vmax + 1e-30, 40)

    fig, ax = plt.subplots(figsize=(5.2, 5.0))
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    sm = plt.cm.ScalarMappable(cmap=cmap,
                               norm=plt.Normalize(vmin=vmin, vmax=vmax))
    fig.colorbar(sm, ax=ax, label=clabel, shrink=0.85)

    def draw(i):
        ax.clear()
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        field = np.broadcast_to(frames[i], RR.shape)
        ax.contourf(X, Y, field, levels=levels, cmap=cmap, extend="both")
        t = f"   t = {times[i]:.2f} s" if times is not None else ""
        ax.set_title(f"{title}{t}")

    anim = FuncAnimation(fig, draw, frames=len(frames), blit=False)
    out = _prepare(path)
    anim.save(str(out), writer=PillowWriter(fps=fps), dpi=dpi)
    plt.close(fig)
    return out


def animate_operating_space(tracks, times=None, *, path, xlabel="", ylabel="",
                            title="", vlines=(), band=None, xlim=None, ylim=None,
                            fps=20, dpi=90, logx=False):
    """Animate several trajectories sweeping together through a 2-D operating diagram.

    `tracks` is a list of dicts {x, y, label, color} that share one time grid (so
    frame i shows each track up to its i-th point with a head marker). `vlines` is a
    list of (x, label) verticals (e.g. the Greenwald density limit); `band` is an
    optional (y_lo, y_hi) shaded horizontal strip (e.g. the burning temperature
    band). The reusable "operation-modes" movie: L-mode / H-mode / disruption tracks
    on one (n, T) plane. Returns the saved Path.
    """
    tracks = list(tracks)
    n_t = len(np.asarray(tracks[0]["x"]))
    all_x = np.concatenate([np.asarray(t["x"], float) for t in tracks])
    all_y = np.concatenate([np.asarray(t["y"], float) for t in tracks])
    xlim = xlim or (all_x.min(), all_x.max() * 1.05 + 1e-30)
    ylim = ylim or (min(0.0, all_y.min()), all_y.max() * 1.1 + 1e-30)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.set(xlim=xlim, ylim=ylim, xlabel=xlabel, ylabel=ylabel)
    if logx:
        ax.set_xscale("log")
    if band is not None:
        ax.axhspan(band[0], band[1], color="0.85", alpha=0.6, zorder=0,
                   label="burning band")
    for xv, lbl in vlines:
        ax.axvline(xv, ls="--", color="crimson", lw=1.1)
        ax.text(xv, ylim[1] * 0.96, lbl, color="crimson", ha="right", va="top",
                rotation=90, fontsize=8)
    lines, heads = [], []
    for tk in tracks:
        (ln,) = ax.plot([], [], color=tk.get("color"), lw=1.8, label=tk.get("label"))
        hd = ax.scatter([], [], s=60, color=tk.get("color"), edgecolors="k",
                        linewidths=0.6, zorder=3)
        lines.append(ln); heads.append(hd)
    ax.legend(loc="upper left", fontsize=8)
    ttl = ax.set_title(title)

    def draw(i):
        for tk, ln, hd in zip(tracks, lines, heads):
            x = np.asarray(tk["x"], float); y = np.asarray(tk["y"], float)
            ln.set_data(x[: i + 1], y[: i + 1])
            hd.set_offsets([[x[i], y[i]]])
        if times is not None:
            ttl.set_text(f"{title}   t = {times[i]:.1f} s")
        return (*lines, *heads, ttl)

    anim = FuncAnimation(fig, draw, frames=n_t, blit=False)
    out = _prepare(path)
    anim.save(str(out), writer=PillowWriter(fps=fps), dpi=dpi)
    plt.close(fig)
    return out


def animate_poloidal_field(R, Z, frames, times=None, *, path, mask=None, title="",
                           cmap="inferno", clabel="", fps=20, dpi=90,
                           vmin=0.0, vmax=None, levels=40):
    """Animate a 2-D field on a real (R, Z) poloidal cross-section over time.

    Unlike `animate_cross_section` (which revolves a 1-D profile into a *circular*
    disk), this draws the field on the actual equilibrium grid — so a shaped /
    D-shaped / outboard-shifted plasma renders with its true geometry. `frames` is
    (n_t, nR, nZ); `mask` (nR, nZ bool) blanks everything outside the plasma. Used
    to render T(rho, t) mapped onto the flux surfaces of a Grad-Shafranov solve.
    """
    R = np.asarray(R, dtype=float)
    Z = np.asarray(Z, dtype=float)
    frames = np.asarray(frames, dtype=float)
    RR, ZZ = np.meshgrid(R, Z, indexing="ij")
    vmax = np.nanmax(frames) if vmax is None else vmax
    lv = np.linspace(vmin, vmax + 1e-30, levels)
    if mask is not None:
        frames = np.where(mask[None], frames, np.nan)

    fig, ax = plt.subplots(figsize=(4.6, 5.6))
    ax.set_aspect("equal")
    ax.set_xlabel("R [m]"); ax.set_ylabel("Z [m]")
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    fig.colorbar(sm, ax=ax, label=clabel, shrink=0.85)

    def draw(i):
        ax.clear()
        ax.set_aspect("equal")
        ax.set_xlabel("R [m]"); ax.set_ylabel("Z [m]")
        ax.contourf(RR, ZZ, frames[i], levels=lv, cmap=cmap, extend="both")
        t = f"   t = {times[i]:.1f} s" if times is not None else ""
        ax.set_title(f"{title}{t}")

    anim = FuncAnimation(fig, draw, frames=len(frames), blit=False)
    out = _prepare(path)
    anim.save(str(out), writer=PillowWriter(fps=fps), dpi=dpi)
    plt.close(fig)
    return out


def animate_torus_nested(rho_levels, T_rt, times=None, *, path, R0=3.0, a=1.0,
                         n_u=70, n_v=36, cmap="inferno", title="", fps=16, dpi=90,
                         vmin=0.0, vmax=None, rotate=True):
    """Animate NESTED toroidal flux surfaces each colored by its own temperature.

    `rho_levels` is a 1-D array of normalized radii (0 < rho <= 1); `T_rt` is
    (n_t, n_levels) the temperature on each surface over time. Each flux surface is
    drawn as a torus of tube radius rho*a, semi-transparent, outer surfaces first so
    the hot core shows through — the honest 3-D successor to `animate_torus_3d`'s
    single-color stand-in: you see the radial T structure, and a sawtooth crash
    visibly flattens (cools) the core surfaces. The camera azimuth sweeps if `rotate`.
    """
    rho_levels = np.asarray(rho_levels, dtype=float)
    T_rt = np.asarray(T_rt, dtype=float)
    n_t = T_rt.shape[0]
    vmax = float(T_rt.max()) if vmax is None else vmax
    norm = plt.Normalize(vmin=vmin, vmax=max(vmax, vmin + 1e-9))
    cmap_obj = matplotlib.colormaps[cmap]
    # precompute each surface geometry (outer -> inner) and a per-surface alpha
    order = np.argsort(rho_levels)[::-1]
    surfaces = [(torus_surface(R0, rho_levels[j] * a, n_u, n_v), j) for j in order]
    alphas = {j: 0.25 + 0.6 * (1.0 - rho_levels[j]) for j in order}   # core more opaque

    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(111, projection="3d")
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    fig.colorbar(sm, ax=ax, label="T [keV]", shrink=0.7)

    def draw(i):
        ax.clear(); ax.set_axis_off()
        for (X, Y, Z), j in surfaces:
            ax.plot_surface(X, Y, Z, color=cmap_obj(norm(T_rt[i, j])), rstride=2,
                            cstride=2, linewidth=0, antialiased=False,
                            alpha=alphas[j], shade=True)
        rng = R0 + a
        ax.set_xlim(-rng, rng); ax.set_ylim(-rng, rng); ax.set_zlim(-rng, rng)
        ax.set_box_aspect((1, 1, 1))
        if rotate:
            ax.view_init(elev=32, azim=(360.0 * i / max(n_t, 1)))
        t = f"   t = {times[i]:.1f} s" if times is not None else ""
        ax.set_title(f"{title}{t}")

    anim = FuncAnimation(fig, draw, frames=n_t, blit=False)
    out = _prepare(path)
    anim.save(str(out), writer=PillowWriter(fps=fps), dpi=dpi)
    plt.close(fig)
    return out


def animate_torus_3d(edge_value, *, path, R0=3.0, a=1.0, n_u=80, n_v=40,
                     cmap="inferno", title="", fps=20, dpi=90, rotate=True,
                     vmin=None, vmax=None):
    """Rotating 3-D torus surface colored by a per-frame scalar (showcase view).

    `edge_value` is a 1-D array (one value per frame) used to color the whole
    surface — a simple stand-in until later rungs color by a poloidally-resolved
    field. The camera azimuth sweeps if `rotate` is True. Returns the saved Path.
    """
    edge_value = np.asarray(edge_value, dtype=float)
    n_t = edge_value.size
    X, Y, Z = torus_surface(R0, a, n_u, n_v)
    vmin = float(edge_value.min()) if vmin is None else vmin
    vmax = float(edge_value.max()) if vmax is None else vmax
    if vmax <= vmin:
        vmax = vmin + 1.0
    norm = plt.Normalize(vmin=vmin, vmax=vmax)

    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(111, projection="3d")
    cmap_obj = matplotlib.colormaps[cmap]

    def draw(i):
        ax.clear()
        ax.set_axis_off()
        facecolor = cmap_obj(norm(edge_value[i]))
        ax.plot_surface(X, Y, Z, color=facecolor, rstride=2, cstride=2,
                        linewidth=0, antialiased=True, shade=True)
        rng = R0 + a
        ax.set_xlim(-rng, rng); ax.set_ylim(-rng, rng); ax.set_zlim(-rng, rng)
        ax.set_box_aspect((1, 1, 1))
        if rotate:
            ax.view_init(elev=35, azim=(360.0 * i / max(n_t, 1)))
        ax.set_title(title)

    anim = FuncAnimation(fig, draw, frames=n_t, blit=False)
    out = _prepare(path)
    anim.save(str(out), writer=PillowWriter(fps=fps), dpi=dpi)
    plt.close(fig)
    return out
