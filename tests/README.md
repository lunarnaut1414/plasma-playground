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
| V3–V15 | see FUNDAMENTALS.md (built alongside their experiments) | ☐ |
