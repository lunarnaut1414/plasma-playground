"""Plotting helpers shared across experiments."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt


def ensure_outputs_dir(experiment_file: str) -> Path:
    """Return (and create) an `outputs/` dir next to the calling experiment."""
    out = Path(experiment_file).resolve().parent / "outputs"
    out.mkdir(exist_ok=True)
    return out


def plot_trajectory_3d(positions, title="Particle trajectory", ax=None):
    """Plot a 3D particle path. `positions` is (N, 3) in metres."""
    if ax is None:
        fig = plt.figure(figsize=(7, 6))
        ax = fig.add_subplot(111, projection="3d")
    ax.plot(positions[:, 0], positions[:, 1], positions[:, 2], lw=0.8)
    ax.scatter(*positions[0], color="green", s=30, label="start")
    ax.scatter(*positions[-1], color="red", s=30, label="end")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_zlabel("z [m]")
    ax.set_title(title)
    ax.legend()
    return ax
