"""Smoke tests: every experiment's run.py executes end-to-end without error.

These don't check physics (that's the V-suite) — they guard against the
experiment *scripts* breaking: an API rename, a numpy-2 gotcha, a bad index. We
load each run.py, stub out the GUI (Agg backend + no-op plt.show), and exercise
every case with tiny particle counts so the whole file stays fast.
"""

import importlib.util
import pathlib

import matplotlib
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

EXPERIMENTS = pathlib.Path(__file__).resolve().parent.parent / "experiments"


def _load(exp_dir):
    path = EXPERIMENTS / exp_dir / "run.py"
    spec = importlib.util.spec_from_file_location(f"run_{exp_dir}", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(autouse=True)
def _no_gui(monkeypatch):
    """Make plt.show() a no-op and close figures so nothing pops up or leaks."""
    monkeypatch.setattr(plt, "show", lambda *a, **k: None)
    yield
    plt.close("all")


def test_exp01_single_particle_runs():
    _load("01_single_particle_motion").main(save=False)


def test_exp04_tokamak_runs():
    # small grid keeps the sparse GS solve quick
    _load("04_tokamak_equilibrium").main(save=False, n=31)


def test_exp05_stellarator_runs():
    # reduced crossing counts keep the field-line tracing quick
    _load("05_stellarator_field_lines").main(save=False, n_crossings=5, profile_n=2)


def test_exp03_cold_runs():
    _load("03_pic_1d").run_cold(save=False, n_particles=2000, n_steps=200)


def test_exp03_landau_runs():
    _load("03_pic_1d").run_landau(save=False, n_particles=3000, n_steps=300)


def test_exp03_twostream_runs():
    _load("03_pic_1d").run_twostream(save=False, n_particles=4000, n_steps=320)


def test_exp06_mhd_briowu_runs():
    _load("06_ideal_mhd").run_briowu(save=False, n=120, t=0.02)


def test_exp06_mhd_alfven_runs():
    _load("06_ideal_mhd").run_alfven(save=False, n=64, t=0.05)


def test_exp07_drive_scaling_runs():
    _load("07_mhd_space_drive").run_scaling(save=False)


def test_exp07_drive_channel_runs():
    _load("07_mhd_space_drive").run_channel(save=False)
