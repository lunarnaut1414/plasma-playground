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


HOUSE_BG = "#0e1116"          # the shared "plasma" dark background
HOUSE_FG = "#f0f0f0"          # light foreground text on that background


def apply_house_style(fig, axes=(), *, dark=True):
    """Apply the gallery's shared 'plasma' look so all gifs read as one designed set.

    Dark figure + axes background with light text (hot colormaps pop against it).
    Pass the figure and any axes (one, or a list/tuple) to restyle. Returns the
    text color to use for titles/labels so callers stay consistent. With dark=False
    it's a no-op that returns black, so a caller can keep one code path.
    """
    if not dark:
        return "black"
    fig.patch.set_facecolor(HOUSE_BG)
    if not isinstance(axes, (list, tuple)):
        axes = [axes]
    for ax in axes:
        ax.set_facecolor(HOUSE_BG)
        # style 2-D axis furniture so light text/ticks read on the dark background
        # (harmless for 3-D axes, which carry no 2-D spines). Re-call after ax.clear().
        for spine in getattr(ax, "spines", {}).values():
            spine.set_color(HOUSE_FG)
            spine.set_alpha(0.4)
        ax.tick_params(colors=HOUSE_FG, labelcolor=HOUSE_FG)
        ax.xaxis.label.set_color(HOUSE_FG)
        ax.yaxis.label.set_color(HOUSE_FG)
        ax.title.set_color(HOUSE_FG)
    return HOUSE_FG


def animate_profiles(x, frames, times=None, *, path, labels=None, xlabel="r/a",
                     ylabel="", title="", fps=20, dpi=90, ylim=None, colors=None,
                     dark=True, shade_between=None, shade_label=""):
    """Animate one or several radial line profiles over time, save to `path` (gif).

    `frames` is (n_t, n_x) for a single curve, or (n_t, n_series, n_x) for several
    curves drawn together (e.g. Te and Ti). `colors` is an optional per-series color
    list. `shade_between=(i, j)` shades the band between series i and j every frame
    (e.g. the Ti-Te equipartition gap) and `shade_label` annotates it. Returns the
    saved Path.
    """
    frames = np.asarray(frames, dtype=float)
    x = np.asarray(x, dtype=float)
    if frames.ndim == 2:
        frames = frames[:, None, :]            # -> (n_t, 1, n_x)
    n_t, n_series, _ = frames.shape
    labels = labels or [None] * n_series
    colors = colors or [None] * n_series
    if ylim is None:
        ylim = (min(0.0, frames.min()), frames.max() * 1.08 + 1e-30)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    txt = apply_house_style(fig, [ax], dark=dark)
    lines = [ax.plot([], [], label=labels[s], color=colors[s], lw=2.1)[0]
             for s in range(n_series)]
    ax.set(xlim=(x.min(), x.max()), ylim=ylim, xlabel=xlabel, ylabel=ylabel)
    if any(lbl is not None for lbl in labels):
        leg = ax.legend(loc="upper right", facecolor=HOUSE_BG if dark else "white",
                        edgecolor=txt, labelcolor=txt if dark else "black")
        leg.get_frame().set_alpha(0.6)
    if shade_between is not None and shade_label:
        ax.text(0.04, 0.06, shade_label, transform=ax.transAxes, color=txt,
                fontsize=8, alpha=0.75)
    ttl = ax.set_title(title, color=txt)
    holder = {"poly": None}

    def draw(i):
        for s, ln in enumerate(lines):
            ln.set_data(x, frames[i, s])
        if shade_between is not None:
            a, b = shade_between
            if holder["poly"] is not None:
                holder["poly"].remove()
            shade_c = colors[b] or txt
            holder["poly"] = ax.fill_between(x, frames[i, a], frames[i, b],
                                             color=shade_c, alpha=0.18, lw=0)
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
                        fps=20, dpi=90, logx=False, dark=True, band=None,
                        band_label="burning band"):
    """Animate a 2-D trajectory (x[i], y[i]) as a growing, fading tail with a head.

    The reusable "phase-space movie" used for burn (n,T) ignition tracks: a point
    sweeps the plane over time, leaving a tail; if `color` (one scalar per frame) is
    given the head is colored by it with a colorbar (e.g. ash fraction). `band` is an
    optional (y_lo, y_hi) shaded horizontal strip — e.g. the efficient-burn (Lawson)
    temperature window the track ignites INTO. Returns the saved Path. Pure rendering
    — the *physics* lives in the (x,y,color) arrays.
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
    txt = apply_house_style(fig, [ax], dark=dark)
    ax.set(xlim=xlim, ylim=ylim, xlabel=xlabel, ylabel=ylabel)
    if logx:
        ax.set_xscale("log")
    if band is not None:
        ax.axhspan(band[0], band[1], color="#ffd166" if dark else "0.85",
                   alpha=0.12 if dark else 0.6, zorder=0)
        ax.text(xlim[0] + 0.03 * (xlim[1] - xlim[0]), band[1], band_label,
                color="#ffd166" if dark else "0.4", fontsize=8, va="top")
    (tail,) = ax.plot([], [], color="#8a8f99", lw=1.4, alpha=0.9)
    head = ax.scatter([x[0]], [y[0]], s=90, zorder=3,
                      c=([color[0]] if color is not None else "crimson"),
                      cmap=cmap, norm=norm, edgecolors=txt, linewidths=0.7)
    if color is not None:
        cb = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), ax=ax,
                          label=clabel, shrink=0.85)
        cb.ax.yaxis.label.set_color(txt); cb.ax.tick_params(colors=txt)
        cb.outline.set_edgecolor(txt)
    ttl = ax.set_title(title, color=txt)

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
                            fps=20, dpi=90, logx=False, dark=True):
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
    txt = apply_house_style(fig, [ax], dark=dark)
    ax.set(xlim=xlim, ylim=ylim, xlabel=xlabel, ylabel=ylabel)
    if logx:
        ax.set_xscale("log")
    if band is not None:
        band_color = "#ffd166" if dark else "0.85"
        ax.axhspan(band[0], band[1], color=band_color, alpha=0.12 if dark else 0.6,
                   zorder=0, label="burning band")
    for xv, lbl in vlines:
        ax.axvline(xv, ls="--", color="#ff5a5a", lw=1.1)
        ax.text(xv, ylim[1] * 0.96, lbl, color="#ff5a5a", ha="right", va="top",
                rotation=90, fontsize=8)
    lines, heads = [], []
    edge = txt if dark else "k"
    for tk in tracks:
        (ln,) = ax.plot([], [], color=tk.get("color"), lw=1.8, label=tk.get("label"))
        hd = ax.scatter([], [], s=60, color=tk.get("color"), edgecolors=edge,
                        linewidths=0.6, zorder=3)
        lines.append(ln); heads.append(hd)
    leg = ax.legend(loc="upper left", fontsize=8,
                    facecolor=HOUSE_BG if dark else "white",
                    edgecolor=txt, labelcolor=txt if dark else "black")
    leg.get_frame().set_alpha(0.6)
    ttl = ax.set_title(title, color=txt)

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
                           vmin=0.0, vmax=None, levels=40, dark=True,
                           show_boundary=True):
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

    bnd = mask.astype(float) if (mask is not None and show_boundary) else None

    fig, ax = plt.subplots(figsize=(4.9, 5.6))
    txt = apply_house_style(fig, [ax], dark=dark)
    ax.set_aspect("equal")
    ax.set_xlabel("R [m]"); ax.set_ylabel("Z [m]")
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    cb = fig.colorbar(sm, ax=ax, label=clabel, shrink=0.85)
    cb.ax.yaxis.label.set_color(txt); cb.ax.tick_params(colors=txt)
    cb.outline.set_edgecolor(txt)

    def draw(i):
        ax.clear()
        ax.set_aspect("equal")
        if dark:
            ax.set_facecolor(HOUSE_BG)
        ax.set_xlabel("R [m]", color=txt); ax.set_ylabel("Z [m]", color=txt)
        ax.tick_params(colors=txt)
        ax.contourf(RR, ZZ, frames[i], levels=lv, cmap=cmap, extend="both")
        if bnd is not None:                                       # crisp plasma boundary
            ax.contour(RR, ZZ, bnd, levels=[0.5], colors=[txt], linewidths=1.1,
                       alpha=0.6)
        t = f"   t = {times[i]:.1f} s" if times is not None else ""
        ax.set_title(f"{title}{t}", color=txt, fontsize=10)

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


def torus_field_lines(R0, a, iota, n_lines=3, *, shape=None, n_tor=3.0, npts=600,
                      rscale=1.0, phase0=0.0):
    """Helical magnetic field lines wrapping a torus surface (the rotational transform).

    Returns a list of (X, Y, Z) curves. Each advances poloidally by `iota` radians per
    radian of toroidal angle u (ι = 1/q is the field-line twist). `shape(v, u)` is an
    optional radial multiplier for a shaped (e.g. l=2 stellarator) tube; default is a
    circular tube of radius `a`. `n_tor` toroidal turns are traced (more -> a denser
    visible helix); `rscale` lifts the line just off the surface so it stays visible.
    """
    u = np.linspace(0.0, 2.0 * np.pi * n_tor, npts)
    lines = []
    for k in range(n_lines):
        v = phase0 + 2.0 * np.pi * k / n_lines + iota * u
        rr = a * rscale * (shape(v, u) if shape is not None else np.ones_like(u))
        rl = R0 + rr * np.cos(v)
        lines.append((rl * np.cos(u), rl * np.sin(u), rr * np.sin(v)))
    return lines


def poloidal_bp_quiver(rings=(0.42, 0.78), n_ang=13, scale=0.17, shape_r=None):
    """Arrows for the poloidal field B_p circulating around the magnetic axis.

    Returns (X, Y, U, V). Arrows sit on the flux surfaces at normalized radii `rings`,
    point in the circulating (poloidal) direction, and grow with minor radius (B_p
    increases outward). `shape_r(v)` optionally scales the surface radius at poloidal
    angle v for a shaped cross-section.
    """
    xs, ys, us, vs = [], [], [], []
    for s in rings:
        for v in np.linspace(0.0, 2.0 * np.pi, n_ang, endpoint=False):
            rb = shape_r(v) if shape_r is not None else 1.0
            r = s * rb
            xs.append(r * np.cos(v)); ys.append(r * np.sin(v))
            us.append(-scale * s * np.sin(v)); vs.append(scale * s * np.cos(v))
    return np.array(xs), np.array(ys), np.array(us), np.array(vs)


def tf_coils(R0, coil_r, n_coils=12, npts=80):
    """Planar toroidal-field coils: circular rings encircling the plasma tube at evenly
    spaced toroidal angles (the iconic tokamak 'rings around the donut'). Each ring's
    field threads through it the long way -> the toroidal field. Returns a list of
    (X, Y, Z) loops."""
    th = np.linspace(0.0, 2.0 * np.pi, npts)
    coils = []
    for k in range(n_coils):
        phi = 2.0 * np.pi * k / n_coils
        rp = R0 + coil_r * np.cos(th)
        coils.append((rp * np.cos(phi), rp * np.sin(phi), coil_r * np.sin(th)))
    return coils


def central_solenoid(R_sol, z_half, n_rings=15, npts=60):
    """Central solenoid: a vertical stack of rings on the machine axis — the transformer
    that drives the tokamak's toroidal plasma current (and hence its poloidal field). A
    stellarator has none. Returns a list of (X, Y, Z) loops."""
    ph = np.linspace(0.0, 2.0 * np.pi, npts)
    return [(R_sol * np.cos(ph), R_sol * np.sin(ph), np.full_like(ph, z))
            for z in np.linspace(-z_half, z_half, n_rings)]


def helical_coils(R0, coil_r, n_coils=2, n_wind=5, npts=700, shape=None):
    """Helical stellarator coils: conductors that wind poloidally `n_wind` times per
    toroidal turn as they go around — the twisted external winding that builds the
    rotational transform from geometry (no plasma current). `shape(theta, phi)` is an
    optional radial multiplier so the coils hug a shaped (l=2) plasma. Returns a list of
    (X, Y, Z) curves."""
    phi = np.linspace(0.0, 2.0 * np.pi, npts)
    coils = []
    for k in range(n_coils):
        th = 2.0 * np.pi * k / n_coils + n_wind * phi
        rr = coil_r * (shape(th, phi) if shape is not None else np.ones_like(phi))
        rp = R0 + rr * np.cos(th)
        coils.append((rp * np.cos(phi), rp * np.sin(phi), rr * np.sin(th)))
    return coils


COIL_COLOR = "#e2e8f0"          # bright silver — the external coil hardware


def _coil_legend(ax, has_field, has_coils, coil_name):
    """Small color-coded legend (top-left of a 3-D panel): field lines + coils."""
    y = 0.97
    if has_field:
        ax.text2D(0.02, y, "— confinement field", transform=ax.transAxes,
                  color="#67e8f9", fontsize=8)
        y -= 0.05
    if has_coils:
        ax.text2D(0.02, y, f"— {coil_name}", transform=ax.transAxes,
                  color=COIL_COLOR, fontsize=8)


def animate_discharge_3d(rho, T_rt, times=None, *, path, R0=3.0, a=1.0, n_u=80, n_v=40,
                         cmap="inferno", title="", fps=16, dpi=120, vmin=0.0, vmax=None,
                         crashes=None, dark=True, field_iota=1.0 / 2.2, n_field=3,
                         field_tor=3.0, show_bp=True, coils=()):
    """Two-panel 3-D discharge: a glowing rotating torus beside its poloidal bullseye.

    `rho` is the normalized minor-radius grid (0..1); `T_rt` is (n_t, n_rho), the
    temperature profile over time. LEFT: the whole torus surface colored by the CORE
    temperature `T_rt[:, 0]`, so the donut brightens through the burn and dims on a
    sawtooth crash (a single honest scalar — the radial structure lives in the right
    panel, not faked onto the 3-D surface). RIGHT: the face-on poloidal cross-section
    T(rho) as a filled 'bullseye' (hot core -> cold edge), where each crash visibly
    flattens the core. `crashes`, if given, is a per-frame count of sawtooth crashes
    since the previous frame -> a flash + running counter (so the fast crashes read as
    activity, not aliasing). The torus sweeps a full turn (seamless loop). Returns the
    saved Path.
    """
    rho = np.asarray(rho, dtype=float)
    T_rt = np.asarray(T_rt, dtype=float)
    n_t = T_rt.shape[0]
    vmax = float(T_rt.max()) if vmax is None else vmax
    norm = plt.Normalize(vmin=vmin, vmax=max(vmax, vmin + 1e-9))
    cmap_obj = matplotlib.colormaps[cmap]
    core = T_rt[:, 0]
    saw_cum = np.cumsum(crashes) if crashes is not None else None

    X, Y, Z = torus_surface(R0, a, n_u, n_v)                 # left: full torus
    ang = np.linspace(0.0, 2.0 * np.pi, 80)                  # right: polar bullseye mesh
    Rg, Ag = np.meshgrid(rho, ang, indexing="ij")
    Xc, Yc = Rg * np.cos(Ag), Rg * np.sin(Ag)
    levels = np.linspace(vmin, max(vmax, vmin + 1e-9), 41)
    # confinement field lines ON two nested flux surfaces (boundary + core); the core
    # twists tighter (q->1 on axis -> larger ι) — the sheared rotational transform
    flines = []
    if field_iota:
        flines += [(c, 1.6, 0.95) for c in
                   torus_field_lines(R0, a, field_iota, n_field, n_tor=field_tor)]
        flines += [(c, 1.2, 0.8) for c in
                   torus_field_lines(R0, 0.55 * a, field_iota * 1.9, n_field,
                                     n_tor=field_tor, phase0=0.6)]
    bpq = poloidal_bp_quiver() if show_bp else None

    fig = plt.figure(figsize=(10.4, 5.0))
    axL = fig.add_subplot(1, 2, 1, projection="3d")
    axR = fig.add_subplot(1, 2, 2)
    txt = apply_house_style(fig, [axR], dark=dark)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cb = fig.colorbar(sm, ax=axR, label="T [keV]", shrink=0.85, pad=0.02)
    cb.ax.yaxis.label.set_color(txt)
    cb.ax.tick_params(colors=txt)
    cb.outline.set_edgecolor(txt)

    if dark:                                          # 3-D panes -> dark, once
        for axis in (axL.xaxis, axL.yaxis, axL.zaxis):
            axis.set_pane_color((0.055, 0.067, 0.086, 1.0))

    def draw(i):
        axL.clear(); axL.set_axis_off()
        if dark:
            axL.patch.set_alpha(0.0)                   # cleared each frame -> re-hide
        for cx, cy, cz in coils:                        # external coils (the hardware)
            axL.plot(cx, cy, cz, color=COIL_COLOR, lw=1.5, alpha=0.85)
        axL.plot_surface(X, Y, Z, color=cmap_obj(norm(core[i])), rstride=2, cstride=2,
                         linewidth=0, antialiased=True, shade=True, alpha=0.5)
        for (fx, fy, fz), lw, al in flines:             # confinement field lines (ι=1/q)
            axL.plot(fx, fy, fz, color="#67e8f9", lw=lw, alpha=al)
        rng = R0 + a * 1.45
        axL.set_xlim(-rng, rng); axL.set_ylim(-rng, rng); axL.set_zlim(-rng, rng)
        axL.set_box_aspect((1, 1, 1))
        axL.view_init(elev=34, azim=360.0 * i / max(n_t, 1))
        axL.set_title("3-D torus (plasma glow = core T)", color=txt, fontsize=10, pad=0)
        _coil_legend(axL, bool(flines), bool(len(coils)), "TF coils + solenoid")
        if crashes is not None and crashes[i] > 0:     # flash over the torus panel
            axL.text2D(0.5, 0.04, "⚡ sawtooth crash", transform=axL.transAxes,
                       color="#ffd166", fontsize=13, fontweight="bold", ha="center")

        axR.clear(); axR.set_aspect("equal"); axR.set_xticks([]); axR.set_yticks([])
        if dark:
            axR.set_facecolor(HOUSE_BG)
        Tc = np.broadcast_to(T_rt[i][:, None], Rg.shape)
        axR.contourf(Xc, Yc, Tc, levels=levels, cmap=cmap, norm=norm, extend="max")
        axR.plot(np.cos(ang), np.sin(ang), color=txt, lw=0.8, alpha=0.5)
        if bpq is not None:                             # poloidal field B_p (circulating)
            axR.quiver(bpq[0], bpq[1], bpq[2], bpq[3], color="#67e8f9", alpha=0.8,
                       width=0.006, scale=2.2, zorder=4)
        axR.set_xlim(-1.12, 1.12); axR.set_ylim(-1.18, 1.12)
        bptxt = "   ·  arrows = $B_p$" if bpq is not None else ""
        axR.set_title(f"poloidal cross-section  T(ρ){bptxt}", color=txt, fontsize=10)
        parts = []
        if times is not None:
            parts.append(f"t = {times[i]:5.1f} s")
        if saw_cum is not None:
            parts.append(f"sawteeth: {int(saw_cum[i]):3d}")
        if parts:
            axR.text(0.0, -1.12, "   ".join(parts), color=txt, fontsize=10,
                     family="monospace", ha="center")

    fig.suptitle(title, color=txt, fontsize=12, y=0.98)
    anim = FuncAnimation(fig, draw, frames=n_t, blit=False)
    out = _prepare(path)
    anim.save(str(out), writer=PillowWriter(fps=fps), dpi=dpi)
    plt.close(fig)
    return out


def animate_stellarator_3d(rho, T_rt, times=None, *, path, R0=3.0, a=1.0, delta=0.30,
                           n_periods=5, n_u=150, n_v=46, n_s=46, cmap="inferno",
                           title="", fps=16, dpi=120, vmin=0.0, vmax=None, dark=True,
                           field_iota=0.45, n_field=3, field_tor=3.0, show_bp=True,
                           coils=()):
    """Two-panel STELLARATOR burn: a glowing twisty torus beside its elliptical bullseye.

    The stellarator analog of `animate_discharge_3d`. `rho` is the normalized minor-radius
    grid (0..1); `T_rt` is (n_t, n_rho), the temperature profile over time. LEFT: a torus
    whose l=2 elliptical cross-section ROTATES helically around the machine (the
    current-free stellarator twist), the whole surface colored by the CORE temperature so
    it brightens through the burn. RIGHT: the shaped (elliptical, l=2) poloidal
    cross-section T(rho) as a filled bullseye. Unlike the tokamak there are NO sawtooth
    crashes (no plasma current -> no q=1 kink): the burn is steady-state, flagged as such.
    The helical shaping geometry is ILLUSTRATIVE (an l=2 envelope mapped onto a torus for
    the render, as in `stellarator_flux_surfaces`); the transport T(rho, t) is the real
    exp-09 burn. The torus sweeps a full turn (seamless loop). Returns the saved Path.
    """
    rho = np.asarray(rho, dtype=float)
    T_rt = np.asarray(T_rt, dtype=float)
    n_t = T_rt.shape[0]
    vmax = float(T_rt.max()) if vmax is None else vmax
    norm = plt.Normalize(vmin=vmin, vmax=max(vmax, vmin + 1e-9))
    cmap_obj = matplotlib.colormaps[cmap]
    core = T_rt[:, 0]

    # left: twisty torus — an l=2 elliptical tube cross-section rotating helically with u
    u = np.linspace(0.0, 2.0 * np.pi, n_u)
    v = np.linspace(0.0, 2.0 * np.pi, n_v)
    U, V = np.meshgrid(u, v, indexing="ij")
    rshape = a * (1.0 + delta * np.cos(2.0 * V - n_periods * U))
    Rl = R0 + rshape * np.cos(V)
    Xs, Ys, Zs = Rl * np.cos(U), Rl * np.sin(U), rshape * np.sin(V)

    # right: l=2 shaped (elliptical) poloidal bullseye mesh
    s = np.linspace(0.0, 1.0, n_s)
    vv = np.linspace(0.0, 2.0 * np.pi, 90)
    Sg, Vg = np.meshgrid(s, vv, indexing="ij")
    rad = Sg * a * (1.0 - delta * np.cos(2.0 * Vg))        # vertically elongated ellipse
    Xc, Yc = rad * np.cos(Vg), rad * np.sin(Vg)
    vb = np.linspace(0.0, 2.0 * np.pi, 220)
    rb = a * (1.0 - delta * np.cos(2.0 * vb))              # plasma boundary ellipse
    levels = np.linspace(vmin, max(vmax, vmin + 1e-9), 41)
    lim = a * (1.0 + delta) * 1.14
    # confinement field lines on the SHAPED nested surfaces (twist ι from geometry, not
    # current); boundary + a core surface, both following the l=2 lobes
    def _shape(vv_, uu_):
        return 1.0 + delta * np.cos(2.0 * vv_ - n_periods * uu_)

    flines = []
    if field_iota:
        flines += [(c, 1.6, 0.95) for c in
                   torus_field_lines(R0, a, field_iota, n_field, n_tor=field_tor,
                                     shape=_shape)]
        flines += [(c, 1.2, 0.8) for c in
                   torus_field_lines(R0, 0.55 * a, field_iota * 1.3, n_field,
                                     n_tor=field_tor, shape=_shape, phase0=0.6)]
    bpq = poloidal_bp_quiver(shape_r=lambda vq: 1.0 - delta * np.cos(2.0 * vq)) \
        if show_bp else None

    fig = plt.figure(figsize=(10.4, 5.0))
    axL = fig.add_subplot(1, 2, 1, projection="3d")
    axR = fig.add_subplot(1, 2, 2)
    txt = apply_house_style(fig, [axR], dark=dark)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cb = fig.colorbar(sm, ax=axR, label="T [keV]", shrink=0.85, pad=0.02)
    cb.ax.yaxis.label.set_color(txt); cb.ax.tick_params(colors=txt)
    cb.outline.set_edgecolor(txt)
    if dark:
        for axis in (axL.xaxis, axL.yaxis, axL.zaxis):
            axis.set_pane_color((0.055, 0.067, 0.086, 1.0))
    rng = R0 + a * 1.7

    def draw(i):
        axL.clear(); axL.set_axis_off()
        if dark:
            axL.patch.set_alpha(0.0)
        for cx, cy, cz in coils:                        # twisted helical coils (hardware)
            axL.plot(cx, cy, cz, color=COIL_COLOR, lw=1.5, alpha=0.85)
        axL.plot_surface(Xs, Ys, Zs, color=cmap_obj(norm(core[i])), rstride=2, cstride=2,
                         linewidth=0, antialiased=True, shade=True, alpha=0.5)
        for (fx, fy, fz), lw, al in flines:             # confinement field lines (geometry)
            axL.plot(fx, fy, fz, color="#67e8f9", lw=lw, alpha=al)
        axL.set_xlim(-rng, rng); axL.set_ylim(-rng, rng); axL.set_zlim(-rng, rng)
        axL.set_box_aspect((1, 1, 1))
        axL.view_init(elev=34, azim=360.0 * i / max(n_t, 1))
        axL.set_title("3-D stellarator (plasma glow = core T)", color=txt, fontsize=10,
                      pad=0)
        _coil_legend(axL, bool(flines), bool(len(coils)), "helical coils")
        axL.text2D(0.5, 0.04, "steady · no sawteeth", transform=axL.transAxes,
                   color="#34d399", fontsize=11, fontweight="bold", ha="center")

        axR.clear(); axR.set_aspect("equal"); axR.set_xticks([]); axR.set_yticks([])
        if dark:
            axR.set_facecolor(HOUSE_BG)
        Tc = np.broadcast_to(np.interp(s, rho, T_rt[i])[:, None], Sg.shape)
        axR.contourf(Xc, Yc, Tc, levels=levels, cmap=cmap, norm=norm, extend="max")
        axR.plot(rb * np.cos(vb), rb * np.sin(vb), color=txt, lw=0.9, alpha=0.5)
        if bpq is not None:                             # poloidal field B_p (circulating)
            axR.quiver(bpq[0], bpq[1], bpq[2], bpq[3], color="#67e8f9", alpha=0.8,
                       width=0.006, scale=2.2, zorder=4)
        axR.set_xlim(-lim, lim); axR.set_ylim(-lim * 1.05, lim)
        bptxt = "  ·  arrows = $B_p$" if bpq is not None else ""
        axR.set_title(f"poloidal cross-section  T(ρ){bptxt}", color=txt, fontsize=9)
        parts = []
        if times is not None:
            parts.append(f"t = {times[i]:5.1f} s")
        parts.append("sawteeth: 0")
        axR.text(0.0, -lim * 0.98, "    ".join(parts), color=txt, fontsize=10,
                 family="monospace", ha="center")

    fig.suptitle(title, color=txt, fontsize=12, y=0.98)
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
