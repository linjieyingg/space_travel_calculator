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
and average_distance_from_earth for the Moon) are present for calculations.

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

# Dictionary containing astrophysical data for solar system bodies.
# Keys are the celestial body names (uppercase for consistent lookup).
# Values are dictionaries with the following properties:
#   'mass': Mass in kilograms (kg).
#   'radius': Equatorial radius in meters (m).
#   'gravitational_parameter_mu': Standard Gravitational Parameter (GM) in m^3/s^2.
#                                 Often more precisely known than G and M independently.
#   'semi_major_axis_from_sun': Average orbital radius from the Sun in meters (m) for planets,
#                               or for the Sun itself (0.0). For the Moon, this is approximated
#                               as Earth's average orbital radius from the Sun for interplanetary
#                               Hohmann transfer contexts.
#   'average_distance_from_earth': Explicitly stored for the Moon as its average distance
#                                  from Earth. For other planets, this is calculated dynamically.
CELESTIAL_BODIES_DATA = {
    "SUN": {
        "mass": 1.9885e30,
        "radius": 6.9634e8,
        "gravitational_parameter_mu": 1.32712440018e20,
        "semi_major_axis_from_sun": 0.0  # The Sun is the primary for other planets
    },
    "MERCURY": {
        "mass": 3.3011e23,
        "radius": 2.4397e6,
        "gravitational_parameter_mu": 2.2032e13,
        "semi_major_axis_from_sun": 5.7909050e10
    },
    "VENUS": {
        "mass": 4.8675e24,
        "radius": 6.0518e6,
        "gravitational_parameter_mu": 3.24859e14,
        "semi_major_axis_from_sun": 1.08208000e11
    },
    "EARTH": {
        "mass": 5.9722e24,
        "radius": 6.3781e6,
        "gravitational_parameter_mu": 3.986004418e14,
        "semi_major_axis_from_sun": 1.49598023e11
    },
    "MOON": {  # Earth's Moon
        "mass": 7.342e22,
        "radius": 1.7374e6,
        "gravitational_parameter_mu": 4.9048695e12,
        "semi_major_axis_from_sun": 1.49598023e11, # Approximated as Earth's for interplanetary context
        "average_distance_from_earth": 3.84400e8 # Average distance from Earth
    },
    "MARS": {
        "mass": 6.4169e23,
        "radius": 3.3895e6,
        "gravitational_parameter_mu": 4.282837e13,
        "semi_major_axis_from_sun": 2.279392e11
    },
    "JUPITER": {
        "mass": 1.8982e27,
        "radius": 7.1492e7,
        "gravitational_parameter_mu": 1.26686534e17,
        "semi_major_axis_from_sun": 7.78412010e11
    },
    "SATURN": {
        "mass": 5.6834e26,
        "radius": 6.0268e7,
        "gravitational_parameter_mu": 3.7931207e16,
        "semi_major_axis_from_sun": 1.426725400e12
    },
    "URANUS": {
        "mass": 8.6810e25,
        "radius": 2.5559e7,
        "gravitational_parameter_mu": 5.793939e15,
        "semi_major_axis_from_sun": 2.870972220e12
    },
    "NEPTUNE": {
        "mass": 1.0241e26,
        "radius": 2.4764e7,
        "gravitational_parameter_mu": 6.83652e15,
        "semi_major_axis_from_sun": 4.498252900e12
    },
    "PLUTO": {  # Dwarf Planet
        "mass": 1.303e22,
        "radius": 1.1883e6,
        "gravitational_parameter_mu": 8.71e11,
        "semi_major_axis_from_sun": 5.906376272e12
    }
}


def get_celestial_body_data(body_name: str) -> dict | None:
    """
    Retrieves the astrophysical data for a given celestial body from the CELESTIAL_BODIES_DATA.

    Args:
        body_name (str): The name of the celestial body (e.g., "EARTH", "MARS").
                         The lookup is case-insensitive.

    Returns:
        dict | None: A dictionary containing the body's data if found, otherwise None.
    """
    if not isinstance(body_name, str):
        return None
    return CELESTIAL_BODIES_DATA.get(body_name.upper())


def get_all_solar_system_destinations() -> list[dict] | None:
    """
    Retrieves a list of all celestial bodies in the solar system suitable as travel destinations,
    excluding Earth, and including their estimated average distance from Earth.

    Returns:
        list[dict] | None: A list of dictionaries, where each dictionary has 'name' (str)
                           and 'distance' (float, average distance from Earth in meters).
                           Returns None if Earth's data is critically missing from
                           CELESTIAL_BODIES_DATA or if essential orbital data is missing.
                           The list is sorted by distance in ascending order.
    """
    destinations = []
    earth_data = get_celestial_body_data("EARTH")
    if earth_data is None:
        print("Error: Earth's data not found in CELESTIAL_BODIES_DATA. Cannot calculate destinations.")
        return None

    earth_orbital_radius_from_sun = earth_data.get("semi_major_axis_from_sun")
    if earth_orbital_radius_from_sun is None:
        print("Error: Earth's 'semi_major_axis_from_sun' data not found. Cannot calculate destinations.")
        return None

    for body_name, data in CELESTIAL_BODIES_DATA.items():
        # Exclude Earth itself as a destination
        if body_name == "EARTH":
            continue

        distance_from_earth = 0.0
        if body_name == "MOON":
            # For the Moon, use its explicitly stored average distance from Earth
            moon_distance_from_earth = data.get("average_distance_from_earth")
            if moon_distance_from_earth is None:
                print(f"Warning: Moon's 'average_distance_from_earth' data not found. Skipping Moon.")
                continue
            distance_from_earth = moon_distance_from_earth
        else:
            # For other planets, calculate the approximate average distance from Earth
            # by taking the absolute difference of their average orbital radii from the Sun.
            # This is a simplification and does not account for specific orbital phases
            # or elliptical orbits for precise travel planning, but serves for average distance.
            body_orbital_radius_from_sun = data.get("semi_major_axis_from_sun")
            if body_orbital_radius_from_sun is None:
                print(f"Warning: {body_name}'s 'semi_major_axis_from_sun' data not found. Skipping {body_name}.")
                continue
            distance_from_earth = math.fabs(body_orbital_radius_from_sun - earth_orbital_radius_from_sun)

        destinations.append({
            "name": body_name.capitalize(), # Capitalize for consistent display
            "distance": distance_from_earth
        })

    # Sort destinations by distance from Earth
    destinations.sort(key=lambda x: x["distance"])

    return destinations