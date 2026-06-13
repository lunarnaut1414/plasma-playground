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


def plasma_frequency(number_density: float, charge: float = e,
                     mass: float = ELECTRON_MASS) -> float:
    """Plasma (Langmuir) angular frequency omega_p = sqrt(n q^2 / (eps0 m)) [rad/s].

    Defaults to electrons. This is the fundamental oscillation frequency of a
    perturbed plasma — the rate at which it springs back toward neutrality.
    """
    return (number_density * charge**2 / (EPSILON_0 * mass)) ** 0.5


def debye_length(temperature_eV: float, number_density: float,
                 charge: float = e) -> float:
    """Debye screening length lambda_D = sqrt(eps0 k T / (n q^2)) [m].

    Temperature given in eV. The distance over which a charge is screened out;
    a plasma must be much larger than lambda_D to behave collectively.
    """
    kT = temperature_eV * ELEMENTARY_CHARGE  # eV -> Joules
    return (EPSILON_0 * kT / (number_density * charge**2)) ** 0.5


def alfven_speed(B: float, number_density: float, mass: float = PROTON_MASS) -> float:
    """Alfven speed v_A = B / sqrt(mu0 rho)  [m/s], with rho = n * mass.

    Defaults to a proton (hydrogen) plasma. The speed at which magnetic-tension
    waves travel along field lines — the characteristic speed of ideal MHD.
    """
    mass_density = number_density * mass
    return B / (MU_0 * mass_density) ** 0.5
