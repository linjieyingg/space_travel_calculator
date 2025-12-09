"""
celestial_data.py

Purpose:
This module serves as the centralized and definitive source for astrophysical data within the
space_travel_calculator application. It stores detailed parameters for celestial bodies in our
solar system and provides functions to retrieve this data.
It also re-exports universal physical constants from the `constants` module for convenience
and consistency, resolving redundancy with previous definitions.

This file consolidates the functionality previously found in `celestial_bodies.py` and
`celestial_data.py` into a single, comprehensive module, ensuring all necessary physical
and orbital properties (mass, radius, gravitational_parameter_mu, semi_major_axis_from_sun,
semi_major_axis_au, and orbital_period_days) are present for calculations.

Crucially, it now includes classical orbital elements (eccentricity, inclination, longitude of
ascending node, argument of periapsis, and mean anomaly at J2000 epoch) for each body. These
elements are primarily heliocentric (relative to the ecliptic J2000.0) for planets and geocentric
(relative to Earth's ecliptic J2000.0) for the Moon, facilitating more accurate 3D orbital
calculations.

Dependencies:
    - constants: For universal physical constants like GRAVITATIONAL_CONSTANT and SPEED_OF_LIGHT.
"""

import math # Used for fabs() in distance calculation
from constants import G_GRAVITATIONAL, C_LIGHT_MPS

# Re-export Universal Constants, maintaining original names for backward compatibility
# if other modules explicitly imported these from an earlier version of celestial_data.py.
# Gravitational Constant (G) in cubic meters per kilogram per square second (m^3 kg^-1 s^-2)
GRAVITATIONAL_CONSTANT = G_GRAVITATIONAL

# Speed of Light in vacuum (c) in meters per second (m/s)
SPEED_OF_LIGHT = C_LIGHT_MPS

# Astronomical Unit (AU) in meters, for consistent conversion
AU_TO_METERS = 149_597_870_700.0

# Astrophysical data for various celestial bodies
# semi_major_axis_from_sun in meters
# semi_major_axis_au in Astronomical Units (calculated from semi_major_axis_from_sun)
# gravitational_parameter_mu = G * mass (in m^3 kg^-1 s^-2 * kg = m^3 s^-2)
# orbital_period_days: Sidereal orbital period around the Sun (or Earth for the Moon)
# New Orbital Elements (J2000.0 Epoch):
# - eccentricity (dimensionless)
# - inclination_deg (degrees, relative to ecliptic J2000.0 for planets, to Earth's ecliptic for Moon)
# - longitude_of_ascending_node_deg (degrees, relative to ecliptic J2000.0 for planets, to Earth's ecliptic for Moon)
# - argument_of_periapsis_deg (degrees, relative to ecliptic J2000.0 for planets, to Earth's ecliptic for Moon)
# - mean_anomaly_at_j2000_deg (degrees, relative to ecliptic J2000.0 for planets, to Earth's ecliptic for Moon)
CELESTIAL_BODIES_DATA = {
    "SUN": {
        "name": "Sun",
        "mass": 1.989e30, # kg
        "radius": 6.957e8, # meters
        "gravitational_parameter_mu": 1.32712440018e20, # m^3/s^2 (G*Mass_Sun)
        "semi_major_axis_from_sun": 0.0, # Sun is the center
        "semi_major_axis_au": 0.0, # Sun is the center
        "orbital_period_days": 0.0, # Not orbiting the sun
        # Classical orbital elements are not applicable for the central body in a heliocentric frame
        "eccentricity": None,
        "inclination_deg": None,
        "longitude_of_ascending_node_deg": None,
        "argument_of_periapsis_deg": None,
        "mean_anomaly_at_j2000_deg": None
    },
    "MERCURY": {
        "name": "Mercury",
        "mass": 3.3011e23,
        "radius": 2.4397e6,
        "gravitational_parameter_mu": 2.2032e13,
        "semi_major_axis_from_sun": 5.790905e10, # meters
        "semi_major_axis_au": 5.790905e10 / AU_TO_METERS,
        "orbital_period_days": 87.969, # Earth days
        "eccentricity": 0.20563593,
        "inclination_deg": 7.0049790, # J2000.0 Ecliptic
        "longitude_of_ascending_node_deg": 48.3307659, # J2000.0 Ecliptic
        "argument_of_periapsis_deg": 29.1253534, # J2000.0 Ecliptic (Longitude of Perihelion - LoAN)
        "mean_anomaly_at_j2000_deg": 174.791572 # J2000.0 Ecliptic
    },
    "VENUS": {
        "name": "Venus",
        "mass": 4.8675e24,
        "radius": 6.0518e6,
        "gravitational_parameter_mu": 3.24859e14,
        "semi_major_axis_from_sun": 1.0820893e11, # meters
        "semi_major_axis_au": 1.0820893e11 / AU_TO_METERS,
        "orbital_period_days": 224.701, # Earth days
        "eccentricity": 0.00677192,
        "inclination_deg": 3.39467605, # J2000.0 Ecliptic
        "longitude_of_ascending_node_deg": 76.67984, # J2000.0 Ecliptic
        "argument_of_periapsis_deg": 54.85314, # J2000.0 Ecliptic (Longitude of Perihelion - LoAN)
        "mean_anomaly_at_j2000_deg": 50.115204 # J2000.0 Ecliptic
    },
    "EARTH": {
        "name": "Earth",
        "mass": 5.972e24,
        "radius": 6.371e6,
        "gravitational_parameter_mu": 3.986004418e14, # m^3/s^2
        "semi_major_axis_from_sun": 1.49598023e11, # meters (1 AU)
        "semi_major_axis_au": 1.49598023e11 / AU_TO_METERS, # ~1.0 AU
        "orbital_period_days": 365.256, # Earth days (sidereal year)
        "eccentricity": 0.01671022,
        "inclination_deg": 0.00005, # J2000.0 Ecliptic (nearly zero by definition of ecliptic)
        "longitude_of_ascending_node_deg": 0.0, # J2000.0 Ecliptic (by definition of ecliptic plane)
        "argument_of_periapsis_deg": 102.94719, # J2000.0 Ecliptic (Longitude of Perihelion - LoAN)
        "mean_anomaly_at_j2000_deg": 100.464488 # J2000.0 Ecliptic
    },
    "MOON": {
        "name": "Moon",
        "mass": 7.342e22,
        "radius": 1.7374e6,
        "gravitational_parameter_mu": 4.9048695e12, # G * Mass_Moon
        "semi_major_axis_from_sun": 1.496e11, # Approximately Earth's distance from Sun (for heliocentric reference)
        "semi_major_axis_au": 1.496e11 / AU_TO_METERS, # Approximately 1.0 AU
        "orbital_period_days": 27.3217, # Earth days (sidereal period around Earth)
        # These are GEOCENTRIC orbital elements for the Moon's orbit around Earth
        "eccentricity": 0.054900,
        "inclination_deg": 5.145, # Relative to Earth's Ecliptic J2000.0
        "longitude_of_ascending_node_deg": 125.08, # Relative to Earth's Ecliptic J2000.0
        "argument_of_periapsis_deg": 318.15, # Relative to Earth's Ecliptic J2000.0
        "mean_anomaly_at_j2000_deg": 115.36 # Relative to Earth's Ecliptic J2000.0
    },
    "MARS": {
        "name": "Mars",
        "mass": 6.4171e23,
        "radius": 3.3895e6,
        "gravitational_parameter_mu": 4.282837e13,
        "semi_major_axis_from_sun": 2.2793911e11, # meters
        "semi_major_axis_au": 2.2793911e11 / AU_TO_METERS,
        "orbital_period_days": 686.971, # Earth days
        "eccentricity": 0.09340056,
        "inclination_deg": 1.849726, # J2000.0 Ecliptic
        "longitude_of_ascending_node_deg": 49.5785367, # J2000.0 Ecliptic
        "argument_of_periapsis_deg": 286.4623033, # J2000.0 Ecliptic (Longitude of Perihelion - LoAN)
        "mean_anomaly_at_j2000_deg": 19.41846 # J2000.0 Ecliptic
    },
    "JUPITER": {
        "name": "Jupiter",
        "mass": 1.8982e27,
        "radius": 6.9911e7,
        "gravitational_parameter_mu": 1.26686534e17,
        "semi_major_axis_from_sun": 7.7841202e11, # meters
        "semi_major_axis_au": 7.7841202e11 / AU_TO_METERS,
        "orbital_period_days": 4332.59, # Earth days
        "eccentricity": 0.048498,
        "inclination_deg": 1.30456, # J2000.0 Ecliptic
        "longitude_of_ascending_node_deg": 100.55615, # J2000.0 Ecliptic
        "argument_of_periapsis_deg": 274.17232, # J2000.0 Ecliptic (Longitude of Perihelion - LoAN)
        "mean_anomaly_at_j2000_deg": 20.02020 # J2000.0 Ecliptic
    },
    "SATURN": {
        "name": "Saturn",
        "mass": 5.6834e26,
        "radius": 5.8232e7,
        "gravitational_parameter_mu": 3.7931207e16,
        "semi_major_axis_from_sun": 1.4293946e12, # meters
        "semi_major_axis_au": 1.4293946e12 / AU_TO_METERS,
        "orbital_period_days": 10759.22, # Earth days
        "eccentricity": 0.0541506,
        "inclination_deg": 2.48446, # J2000.0 Ecliptic
        "longitude_of_ascending_node_deg": 113.6655, # J2000.0 Ecliptic
        "argument_of_periapsis_deg": 338.76644, # J2000.0 Ecliptic (Longitude of Perihelion - LoAN)
        "mean_anomaly_at_j2000_deg": 320.3256 # J2000.0 Ecliptic
    },
    "URANUS": {
        "name": "Uranus",
        "mass": 8.6810e25,
        "radius": 2.5362e7,
        "gravitational_parameter_mu": 5.793939e15,
        "semi_major_axis_from_sun": 2.8709722e12, # meters
        "semi_major_axis_au": 2.8709722e12 / AU_TO_METERS,
        "orbital_period_days": 30687.15, # Earth days
        "eccentricity": 0.04716771,
        "inclination_deg": 0.76986, # J2000.0 Ecliptic
        "longitude_of_ascending_node_deg": 74.22988, # J2000.0 Ecliptic
        "argument_of_periapsis_deg": 96.73436, # J2000.0 Ecliptic (Longitude of Perihelion - LoAN)
        "mean_anomaly_at_j2000_deg": 314.2023 # J2000.0 Ecliptic
    },
    "NEPTUNE": {
        "name": "Neptune",
        "mass": 1.02413e26,
        "radius": 2.4622e7,
        "gravitational_parameter_mu": 6.836529e15,
        "semi_major_axis_from_sun": 4.4984284e12, # meters
        "semi_major_axis_au": 4.4984284e12 / AU_TO_METERS,
        "orbital_period_days": 60190.03, # Earth days
        "eccentricity": 0.00858587,
        "inclination_deg": 1.76917, # J2000.0 Ecliptic
        "longitude_of_ascending_node_deg": 131.7820, # J2000.0 Ecliptic
        "argument_of_periapsis_deg": 273.18935, # J2000.0 Ecliptic (Longitude of Perihelion - LoAN)
        "mean_anomaly_at_j2000_deg": 304.34866 # J2000.0 Ecliptic
    }
}

def get_celestial_body_data(body_name: str) -> dict | None:
    """
    Retrieves astrophysical data for a given celestial body.

    Args:
        body_name (str): The case-insensitive name of the celestial body.

    Returns:
        dict: A dictionary containing the body's full data if found, or None.
    """
    if not isinstance(body_name, str):
        return None
    return CELESTIAL_BODIES_DATA.get(body_name.upper())