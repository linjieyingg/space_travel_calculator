"""
celestial_data.py

Purpose:
This module serves as the centralized and definitive source for astrophysical data within the
space_travel_calculator application. It stores detailed parameters for celestial bodies in our
solar system and provides functions to retrieve this data and list potential destinations.
It also re-exports universal physical constants from the `constants` module for convenience
and consistency, resolving redundancy with previous definitions.

This file consolidates the functionality previously found in `celestial_bodies.py` and
`celestial_data.py` into a single, comprehensive module, ensuring all necessary physical
and orbital properties (mass, radius, gravitational_parameter_mu, semi_major_axis_from_sun,
average_distance_from_earth for the Moon, and orbital_period_days) are present for calculations.

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

# Astrophysical data for various celestial bodies
# semi_major_axis_from_sun in meters
# gravitational_parameter_mu = G * mass (in m^3 kg^-1 s^-2 * kg = m^3 s^-2)
# orbital_period_days: Sidereal orbital period around the Sun (or Earth for the Moon)
CELESTIAL_BODIES_DATA = {
    "SUN": {
        "name": "Sun",
        "mass": 1.989e30, # kg
        "radius": 6.957e8, # meters
        "gravitational_parameter_mu": 1.32712440018e20, # m^3/s^2 (G*Mass_Sun)
        "semi_major_axis_from_sun": 0.0, # Sun is the center
        "average_distance_from_earth": 1.496e11, # For reference, not direct orbit
        "orbital_period_days": 0.0 # Not orbiting the sun
    },
    "MERCURY": {
        "name": "Mercury",
        "mass": 3.3011e23,
        "radius": 2.4397e6,
        "gravitational_parameter_mu": 2.2032e13,
        "semi_major_axis_from_sun": 5.7909e10, # meters
        "orbital_period_days": 87.97 # Earth days
    },
    "VENUS": {
        "name": "Venus",
        "mass": 4.8675e24,
        "radius": 6.0518e6,
        "gravitational_parameter_mu": 3.24859e14,
        "semi_major_axis_from_sun": 1.0821e11, # meters
        "orbital_period_days": 224.70 # Earth days
    },
    "EARTH": {
        "name": "Earth",
        "mass": 5.972e24,
        "radius": 6.371e6,
        "gravitational_parameter_mu": 3.986004418e14, # m^3/s^2
        "semi_major_axis_from_sun": 1.496e11, # meters (1 AU)
        "orbital_period_days": 365.256 # Earth days (sidereal year)
    },
    "MOON": {
        "name": "Moon",
        "mass": 7.342e22,
        "radius": 1.7374e6,
        "gravitational_parameter_mu": 4.9048695e12, # G * Mass_Moon
        "semi_major_axis_from_sun": 1.496e11, # Approximately Earth's distance from Sun
        "average_distance_from_earth": 3.844e8, # Average distance from Earth in meters
        "orbital_period_days": 27.3217 # Earth days (sidereal period around Earth)
    },
    "MARS": {
        "name": "Mars",
        "mass": 6.4171e23,
        "radius": 3.3895e6,
        "gravitational_parameter_mu": 4.282837e13,
        "semi_major_axis_from_sun": 2.2792e11, # meters
        "orbital_period_days": 686.98 # Earth days
    },
    "JUPITER": {
        "name": "Jupiter",
        "mass": 1.8982e27,
        "radius": 6.9911e7,
        "gravitational_parameter_mu": 1.26686534e17,
        "semi_major_axis_from_sun": 7.7857e11, # meters
        "orbital_period_days": 4332.59 # Earth days
    },
    "SATURN": {
        "name": "Saturn",
        "mass": 5.6834e26,
        "radius": 5.8232e7,
        "gravitational_parameter_mu": 3.7931207e16,
        "semi_major_axis_from_sun": 1.4335e12, # meters
        "orbital_period_days": 10759.22 # Earth days
    },
    "URANUS": {
        "name": "Uranus",
        "mass": 8.6810e25,
        "radius": 2.5362e7,
        "gravitational_parameter_mu": 5.793939e15,
        "semi_major_axis_from_sun": 2.8709e12, # meters
        "orbital_period_days": 30687.15 # Earth days
    },
    "NEPTUNE": {
        "name": "Neptune",
        "mass": 1.02413e26,
        "radius": 2.4622e7,
        "gravitational_parameter_mu": 6.836529e15,
        "semi_major_axis_from_sun": 4.4983e12, # meters
        "orbital_period_days": 60190.03 # Earth days
    }
}

def get_celestial_body_data(body_name: str) -> dict | None:
    """
    Retrieves astrophysical data for a given celestial body.

    Args:
        body_name (str): The case-insensitive name of the celestial body.

    Returns:
        dict: A dictionary containing the body's data if found, otherwise None.
    """
    if not isinstance(body_name, str):
        return None
    return CELESTIAL_BODIES_DATA.get(body_name.upper())

def get_all_solar_system_destinations() -> list[dict] | None:
    """
    Generates a list of all solar system bodies as potential destinations from Earth,
    sorted by their average distance from Earth.

    Returns:
        list: A list of dictionaries, each with 'name' and 'distance' from Earth in meters.
              Returns None if Earth's data is missing.
    """
    earth_data = get_celestial_body_data('Earth')
    if not earth_data or 'semi_major_axis_from_sun' not in earth_data:
        # Earth data is critical for distance calculations, return None if unavailable
        return None

    earth_semi_major_axis = earth_data['semi_major_axis_from_sun']
    destinations = []

    for body_name_key, data in CELESTIAL_BODIES_DATA.items():
        if body_name_key == 'EARTH':
            continue

        distance = 0.0
        # Check for 'name' key, as it's used in the output dictionary
        if 'name' not in data:
            continue # Skip bodies that don't have a 'name' defined

        if body_name_key == 'MOON' and 'average_distance_from_earth' in data:
            distance = data['average_distance_from_earth']
        elif 'semi_major_axis_from_sun' in data:
            distance = math.fabs(data['semi_major_axis_from_sun'] - earth_semi_major_axis)
        else:
            # Skip bodies without sufficient data to calculate distance
            continue

        destinations.append({'name': data['name'], 'distance': distance})

    # Sort destinations by distance
    destinations.sort(key=lambda x: x['distance'])

    return destinations