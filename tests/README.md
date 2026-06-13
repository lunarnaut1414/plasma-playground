# Validation suite

Falsifiable tests for the shared `plasmaplay` kernels — each pins a kernel
against an analytic result or a published reference. Test IDs (V0, V1, ...) map
to the suite in [`../docs/FUNDAMENTALS.md`](../docs/FUNDAMENTALS.md) §2.

```bash
pytest                 # run everything
pytest -v              # show each V-test by name
pytest -k constants    # just the V0 group
```

A kernel is not "done" until its V-test is green. Tests use only numpy/scipy by
default; cross-checks against PlasmaPy are skipped automatically if it isn't
installed.

| ID | Proves | Status |
|----|--------|--------|
| V0 | constants / formulary (ω_pe, λ_D, v_A, ω_c, r_L) | ✅ |
| V1 | Boris pusher energy conservation + gyroradius | ✅ |
| V2 | E×B drift = E/B | ✅ |
| V3 | RK4 integrator 4th-order convergence | ✅ |
| V4 | FFT Poisson 1-D vs analytic (spectral) | ✅ |
| V5 | cold PIC oscillation at ω_pe | ✅ |
| V6 | Landau damping rate | ✅ |
| V7 | two-stream growth rate | ✅ |
| V10 | FFT Poisson 2-D vs analytic | ✅ |
| V11 | Biot–Savart loop vs analytic on-axis field | ✅ |
| V13 | field-line tracing → ι; flux surfaces close | ✅ |
| others | see FUNDAMENTALS.md (built alongside their experiments) | ☐ |

Beyond the physics-validation V-tests, each step-1 module also has **unit-level
tests** of its contract and edge cases:

- `test_constants.py` — aliases, scaling laws (ω_c∝B, λ_D∝√T, …), and
  cross-relations (e.g. λ_D = v_th / (√2 ω_pe)).
- `test_fields.py` — field factories return constant (3,) vectors; the magnetic
  mirror's on-axis profile, symmetry, inward radial field, and **∇·B = 0**.
- `test_pushers.py` — output shapes, zero-field straight-line motion, pure-E
  uniform acceleration, charge-sign reversal, one-period return, list/array parity.
- `test_integrators.py` — RK4 convergence order + exactness on a constant RHS.
- `test_diagnostics.py` — field-line tracing stays on flux surfaces; ι is
  radius-independent for a shearless profile.
- `test_solvers.py` — spectral Poisson exactness (1-D & 2-D), E = -∇φ, zero-mean gauge.
- `test_pic.py` — CIC charge conservation, partition of unity, linear exactness; loaders.
