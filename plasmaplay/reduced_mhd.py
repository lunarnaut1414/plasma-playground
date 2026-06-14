"""Nonlinear 2-D reduced MHD — a tearing mode growing into a magnetic island (B2).

B1 (`cylinder_mhd.py`) gave the *linear* stability of a current profile — does a
tearing mode grow, and how fast. B2 follows the mode into the **nonlinear** regime,
where it stops growing exponentially and **saturates** into a finite magnetic island
(the Rutherford regime). This is the laptop-scale cousin of JOREK/NIMROD.

We solve the **Strauss reduced-MHD equations** in a 2-D slab — the periodic-y
"straightened" tokamak around one rational surface — for the poloidal flux psi and
the vorticity U = laplacian(phi) (phi the velocity stream function, v = z x grad phi):

    d psi/dt = -[phi, psi] + eta (lap psi - J_eq)
    d U/dt   = -[phi, U]   + [psi, lap psi] + nu lap U

with the Poisson bracket [a,b] = a_x b_y - a_y b_x (advection / Lorentz drive) and
J = lap psi the current. The equilibrium is the Harris sheet B_y0(x) = tanh(x),
psi0(x) = ln cosh(x), J_eq = sech^2 x — the same sheet whose linear tearing growth
and Delta' are validated analytically in `tearing.py` (T4). Reconnection at the
neutral line x = 0 grows an island of width W = 4 sqrt(psi_island / B_y0'(0)).

Numerics: finite differences in x (Dirichlet walls), spectral (FFT) in the periodic
y, an FFT + tridiagonal elliptic solve for phi from U, and an SSP-RK2 time step.
Normalisation: lengths in the sheet width a=1, B in its asymptotic value, time in the
Alfven time tau_A, resistivity eta = 1/S (S the Lundquist number).
"""

from __future__ import annotations

import numpy as np


def _thomas_vec(lower, diag, upper, rhs):
    """Thomas algorithm solving many tridiagonal systems at once.

    `diag` is (N, M) — N rows, M independent systems (the ky columns); `lower`,
    `upper` are (N, M) off-diagonals; `rhs` is (N, M). The two O(N) sweeps loop over
    the N rows but are vectorised across the M systems, so the per-ky Python loop is
    gone. Complex-capable. Returns (N, M)."""
    n = diag.shape[0]
    c = upper.astype(complex).copy()
    d = rhs.astype(complex).copy()
    b = diag.astype(complex).copy()
    c[0] /= b[0]
    d[0] /= b[0]
    for i in range(1, n):
        m = b[i] - lower[i] * c[i - 1]
        c[i] = upper[i] / m
        d[i] = (d[i] - lower[i] * d[i - 1]) / m
    x = d
    for i in range(n - 2, -1, -1):
        x[i] -= c[i] * x[i + 1]
    return x


class ReducedMHD:
    """Strauss reduced-MHD on a 2-D slab (x finite-difference, y spectral).

    A Harris current sheet B_y0 = tanh(x) on x in [-Lx, Lx] (Dirichlet walls),
    periodic in y over one wavelength Ly = 2 pi / k so the box holds exactly the
    fundamental tearing harmonic k. Construct, `seed` a small perturbation, then
    `step`/`run`; read the island via `mode_amplitude` / `island_width`.
    """

    def __init__(self, k, *, S=1e3, Pm=1.0, Lx=6.0, nx=192, ny=32):
        self.k = float(k)
        self.eta = 1.0 / S
        self.nu = Pm / S                      # viscosity = Pm * eta (magnetic Prandtl)
        self.Lx = Lx
        self.nx, self.ny = nx, ny
        self.x = np.linspace(-Lx, Lx, nx)
        self.dx = self.x[1] - self.x[0]
        self.Ly = 2.0 * np.pi / self.k
        self.y = np.linspace(0.0, self.Ly, ny, endpoint=False)
        self.dy = self.y[1] - self.y[0]
        self.ky = 2.0 * np.pi * np.fft.rfftfreq(ny, d=self.dy)   # (ny//2+1,)
        self.ix0 = int(np.argmin(np.abs(self.x)))               # neutral-line index

        X = self.x[:, None]
        self.psi0 = np.log(np.cosh(X)) * np.ones((nx, ny))      # equilibrium flux
        self.By0 = np.tanh(self.x)                              # B_y0(x) = psi0'
        self.J_eq = (1.0 / np.cosh(X) ** 2) * np.ones((nx, ny))  # sech^2 x
        self.psi = self.psi0.copy()
        self.U = np.zeros((nx, ny))
        self.t = 0.0

    # --- spectral y-derivatives / finite-difference x-derivatives -----------
    def _dy(self, f):
        return np.fft.irfft(1j * self.ky[None, :] * np.fft.rfft(f, axis=1), n=self.ny, axis=1)

    def _dyy(self, f):
        return np.fft.irfft(-self.ky[None, :] ** 2 * np.fft.rfft(f, axis=1), n=self.ny, axis=1)

    def _dx(self, f):
        g = np.zeros_like(f)
        g[1:-1] = (f[2:] - f[:-2]) / (2 * self.dx)
        g[0] = (f[1] - f[0]) / self.dx
        g[-1] = (f[-1] - f[-2]) / self.dx
        return g

    def _dxx(self, f):
        g = np.zeros_like(f)
        g[1:-1] = (f[2:] - 2 * f[1:-1] + f[:-2]) / self.dx ** 2
        g[0] = g[1]
        g[-1] = g[-2]
        return g

    def _lap(self, f):
        return self._dxx(f) + self._dyy(f)

    def _bracket(self, a, b):
        return self._dx(a) * self._dy(b) - self._dy(a) * self._dx(b)

    def invert_phi(self, U):
        """Solve lap(phi) = U with phi = 0 on the x-walls, periodic in y.

        FFT in y decouples the ky modes; each is a 1-D Helmholtz solve
        (d^2/dx^2 - ky^2) phi_hat = U_hat by the tridiagonal Thomas algorithm.
        """
        Uhat = np.fft.rfft(U, axis=1)                 # (nx, nky)
        nx, dx = self.nx, self.dx
        nky = self.ky.size
        off = 1.0 / dx ** 2
        diag = np.empty((nx, nky), complex)
        diag[:] = -2.0 / dx ** 2 - self.ky[None, :] ** 2
        diag[0] = diag[-1] = 1.0                       # Dirichlet phi=0 at walls
        lower = np.full((nx, nky), off, complex); lower[-1] = 0.0
        upper = np.full((nx, nky), off, complex); upper[0] = 0.0
        rhs = Uhat.copy(); rhs[0] = 0.0; rhs[-1] = 0.0
        phihat = _thomas_vec(lower, diag, upper, rhs)
        return np.fft.irfft(phihat, n=self.ny, axis=1)

    # --- seeding, RHS, stepping ---------------------------------------------
    def seed(self, amp=1e-4):
        """Seed a small m=1 flux perturbation ~ amp cos(k y) sech(x) (vanishes at walls)."""
        env = 1.0 / np.cosh(self.x)[:, None]
        self.psi = self.psi0 + amp * env * np.cos(self.k * self.y)[None, :]
        self.U = np.zeros((self.nx, self.ny))
        self.t = 0.0
        return self

    def _rhs(self, psi, U):
        phi = self.invert_phi(U)
        J = self._lap(psi)
        dpsi = -self._bracket(phi, psi) + self.eta * (J - self.J_eq)
        dU = -self._bracket(phi, U) + self._bracket(psi, J) + self.nu * self._lap(U)
        # Dirichlet walls: hold the equilibrium flux, no vorticity at the walls
        dpsi[0] = dpsi[-1] = 0.0
        dU[0] = dU[-1] = 0.0
        return dpsi, dU

    def step(self, dt):
        """One SSP-RK2 (Heun) step of dt Alfven times."""
        k1p, k1u = self._rhs(self.psi, self.U)
        p1 = self.psi + dt * k1p
        u1 = self.U + dt * k1u
        k2p, k2u = self._rhs(p1, u1)
        self.psi = self.psi + 0.5 * dt * (k1p + k2p)
        self.U = self.U + 0.5 * dt * (k1u + k2u)
        self.t += dt
        return self

    def run(self, t_end, dt):
        n = int(round(t_end / dt))
        for _ in range(n):
            self.step(dt)
        return self

    # --- diagnostics ---------------------------------------------------------
    def mode_amplitude(self):
        """Amplitude of the reconnected (m=1) flux at the neutral line x=0."""
        row = self.psi[self.ix0] - self.psi0[self.ix0]
        return 2.0 / self.ny * np.abs(np.fft.rfft(row)[1])

    def island_width(self):
        """Magnetic-island full width W = 4 sqrt(psi_rec / B_y0'(0)), B_y0'(0)=1."""
        return 4.0 * np.sqrt(max(self.mode_amplitude(), 0.0))

    def flux_function(self):
        """The total poloidal flux psi(x, y) (its contours are the field lines /
        island separatrix) — what the island gif renders."""
        return self.psi
