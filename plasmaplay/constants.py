"""Physical constants in SI units.

Thin wrapper over scipy.constants so experiments read clearly. For anything
serious, prefer astropy.units / plasmapy which carry units through the math.
"""

from scipy import constants as _c

# Fundamental
ELEMENTARY_CHARGE = _c.elementary_charge        # C
SPEED_OF_LIGHT = _c.speed_of_light              # m/s
BOLTZMANN = _c.Boltzmann                        # J/K
EPSILON_0 = _c.epsilon_0                         # F/m
MU_0 = _c.mu_0                                   # H/m

# Particle masses (kg)
ELECTRON_MASS = _c.electron_mass
PROTON_MASS = _c.proton_mass

# Convenience aliases used in the experiments
e = ELEMENTARY_CHARGE
c = SPEED_OF_LIGHT
m_e = ELECTRON_MASS
m_p = PROTON_MASS


def thermal_velocity(temperature_eV: float, mass: float) -> float:
    """Thermal speed v_th = sqrt(2 k T / m), with temperature given in eV."""
    kT = temperature_eV * ELEMENTARY_CHARGE  # eV -> Joules
    return (2.0 * kT / mass) ** 0.5


def gyrofrequency(charge: float, B: float, mass: float) -> float:
    """Cyclotron (gyro) angular frequency omega_c = qB/m  [rad/s]."""
    return abs(charge) * B / mass


def gyroradius(velocity_perp: float, charge: float, B: float, mass: float) -> float:
    """Larmor radius r_L = m v_perp / (|q| B)  [m]."""
    return mass * velocity_perp / (abs(charge) * B)
