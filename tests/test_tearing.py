"""T4 — the linear resistive tearing mode (`plasmaplay.tearing`).

The two falsifiable rung-T4 checks from `docs/3D_TOKAMAK_GUIDE.md`:

  * the ideal outer region reproduces the **exact analytic Δ'** of the tanh sheet,
    Δ' = (2/a)(1/(ka) − ka) — and so the mode is unstable only for k a < 1;
  * the resistive growth rate obeys the Furth–Killeen–Rosenbluth scaling
    **γ ∝ S^(−3/5)** (S = Lundquist number).
"""

import numpy as np

from plasmaplay.tearing import (
    delta_prime_analytic,
    delta_prime_slab,
    tearing_growth_rate,
)


# --- Δ': the ideal outer region ------------------------------------------

def test_delta_prime_matches_analytic_tanh_sheet():
    # Newcomb integration must reproduce Δ' = (2/a)(1/(ka) − ka) across k.
    for ka in (0.3, 0.5, 0.7, 0.9):
        num = delta_prime_slab(ka)
        ana = delta_prime_analytic(ka)
        assert abs(num - ana) / abs(ana) < 0.01


def test_tearing_threshold_at_ka_equals_one():
    # Δ' > 0 (tearing-unstable) for k a < 1, Δ' < 0 (stable) for k a > 1, ≈ 0 at 1.
    assert delta_prime_slab(0.6) > 0
    assert delta_prime_slab(1.4) < 0
    assert abs(delta_prime_slab(1.0)) < 0.02


# --- γ: the resistive layer and the S^(−3/5) law -------------------------

def test_unstable_mode_has_positive_growth_rate():
    # k a = 0.8 (Δ' > 0) is tearing-unstable: a real, positive growth rate exists.
    g = tearing_growth_rate(0.8, 3e4, N=2000)
    assert g > 0
    assert g < 0.1                     # and it is slow (resistive), not ideal-fast


def test_growth_rate_scales_as_S_minus_three_fifths():
    # The headline T4 result: γ ∝ S^(−3/5). Fit the slope of ln γ vs ln S in the
    # constant-ψ regime (small Δ', k a = 0.8) over a decade-resolved S range.
    Ss = np.array([1e4, 3e4, 1e5, 3e5])
    gam = np.array([tearing_growth_rate(0.8, S, N=3000) for S in Ss])
    assert np.all(gam > 0)
    slope = np.polyfit(np.log(Ss), np.log(gam), 1)[0]
    assert -0.68 < slope < -0.52       # FKR: −0.60
