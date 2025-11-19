"""
This module stores structured data for solar system celestial bodies (planets, sun, moon).
It includes their physical properties (mass, radius), standard gravitational parameter (mu),
and orbital characteristics (e.g., semi-major axis) relevant for orbital mechanics
calculations and simplified space travel planning.
"""

# Gravitational Constant (G) in m^3 kg^-1 s^-2.
# This constant is provided for reference, though not directly used in the functions below,
# it is fundamental for many orbital mechanics calculations.
GRAVITATIONAL_CONSTANT = 6.67430e-11

# SOLAR_SYSTEM_BODIES_DATA
# A dictionary containing data for various solar system bodies.
#
# Each entry is a dictionary with the following keys:
# - 'mass': The mass of the celestial body in kilograms (kg).
# - 'radius': The mean equatorial radius of the celestial body in meters (m).
# - 'gravitational_parameter_mu': The standard gravitational parameter (G * mass) in m^3/s^2.
#   This is often used in orbital mechanics calculations as it combines G and M into one precise value.
# - 'semi_major_axis_from_sun': The average orbital radius (semi-major axis) from the Sun in meters (m).
#   For the Sun itself, this value is 0.0. For the Moon, it's approximately Earth's semi-major axis.
# - 'average_distance_from_earth': An approximate average distance from Earth in meters (m).
#   This value is a simplification for basic travel calculations, as actual distances vary significantly
#   due to orbital mechanics. For Earth, this value is 0.0 (as it's the assumed origin).
#   For the Moon, it's the average Earth-Moon distance.
SOLAR_SYSTEM_BODIES_DATA = {
    "Sun": {
        "mass": 1.9885e+30,
        "radius": 6.957e+8,
        "gravitational_parameter_mu": 1.32712440018e+20,
        "semi_major_axis_from_sun": 0.0,
        "average_distance_from_earth": 1.496e+11,  # Approximately 1 AU
    },
    "Mercury": {
        "mass": 3.3011e+23,
        "radius": 2.4397e+6,
        "gravitational_parameter_mu": 2.2032e+13,
        "semi_major_axis_from_sun": 5.7909227e+10,
        "average_distance_from_earth": 1.22e+11,  # Simplified average from Earth
    },
    "Venus": {
        "mass": 4.8675e+24,
        "radius": 6.0518e+6,
        "gravitational_parameter_mu": 3.24859e+14,
        "semi_major_axis_from_sun": 1.08208e+11,
        "average_distance_from_earth": 4.14e+10,  # Simplified average from Earth
    },
    "Earth": {
        "mass": 5.972168e+24,
        "radius": 6.371e+6,
        "gravitational_parameter_mu": 3.986004418e+14,
        "semi_major_axis_from_sun": 1.49598023e+11,
        "average_distance_from_earth": 0.0,  # Origin planet
    },
    "Moon": {  # Earth's Moon
        "mass": 7.342e+22,
        "radius": 1.7374e+6,
        "gravitational_parameter_mu": 4.9048695e+12,
        "semi_major_axis_from_sun": 1.49598023e+11,  # Approximately Earth's orbit
        "average_distance_from_earth": 3.844e+8,  # Average Earth-Moon distance
    },
    "Mars": {
        "mass": 6.4171e+23,
        "radius": 3.3895e+6,
        "gravitational_parameter_mu": 4.282837e+13,
        "semi_major_axis_from_sun": 2.279391e+11,
        "average_distance_from_earth": 7.83e+10,  # Simplified average from Earth
    },
    "Jupiter": {
        "mass": 1.89819e+27,
        "radius": 6.9911e+7,
        "gravitational_parameter_mu": 1.26686534e+17,
        "semi_major_axis_from_sun": 7.7857e+11,
        "average_distance_from_earth": 6.287e+11,  # Simplified average from Earth
    },
    "Saturn": {
        "mass": 5.6834e+26,
        "radius": 5.8232e+7,
        "gravitational_parameter_mu": 3.7931187e+16,
        "semi_major_axis_from_sun": 1.43353e+12,
        "average_distance_from_earth": 1.277e+12,  # Simplified average from Earth
    },
    "Uranus": {
        "mass": 8.6810e+25,
        "radius": 2.5362e+7,
        "gravitational_parameter_mu": 5.793939e+15,
        "semi_major_axis_from_sun": 2.87097e+12,
        "average_distance_from_earth": 2.723e+12,  # Simplified average from Earth
    },
    "Neptune": {
        "mass": 1.02413e+26,
        "radius": 2.4622e+7,
        "gravitational_parameter_mu": 6.83652e+15,
        "semi_major_axis_from_sun": 4.49839e+12,
        "average_distance_from_earth": 4.347e+12,  # Simplified average from Earth
    },
}

def get_celestial_body_data(name: str) -> dict | None:
    """
    Retrieves the structured data for a specific celestial body.

    Args:
        name (str): The name of the celestial body (e.g., "Earth", "Mars", "Moon").
                    The lookup is case-insensitive (e.g., "mars" will find "Mars").

    Returns:
        dict or None: A dictionary containing the body's mass, radius, gravitational
                      parameter (mu), semi-major axis from the Sun, and average
                      distance from Earth. Returns None if the body is not found
                      or if the input `name` is not a string.
    """
    if not isinstance(name, str):
        # Handle non-string input gracefully
        return None
    return SOLAR_SYSTEM_BODIES_DATA.get(name.capitalize())

def get_all_solar_system_destinations() -> list[dict]:
    """
    Retrieves a list of solar system bodies formatted as potential travel destinations.
    Each dictionary contains 'name' and 'distance' (average distance from Earth).
    This format is designed to be compatible with the `planets` list used in `main.py`
    for the space travel calculator. Earth is excluded as it is typically the origin.

    Returns:
        list[dict]: A list of dictionaries, where each dictionary has:
                    - 'name' (str): The name of the celestial body.
                    - 'distance' (float): The 'average_distance_from_earth' for that body in meters.
                    The list is sorted by distance for consistent output.
    """
    destinations = []
    for name, data in SOLAR_SYSTEM_BODIES_DATA.items():
        if name != "Earth":  # Exclude Earth as it's the presumed origin
            destinations.append({
                "name": name,
                "distance": data["average_distance_from_earth"]
            })
    # Sort the destinations by their distance from Earth for a logical order
    destinations.sort(key=lambda x: x["distance"])
    return destinations

if __name__ == "__main__":
    print("--- Solar System Celestial Bodies Data Test ---")

    # Test retrieval of specific body data
    earth_data = get_celestial_body_data("Earth")
    if earth_data:
        print("\nEarth Data:")
        for key, value in earth_data.items():
            print(f"  {key}: {value}")
    else:
        print("\nEarth data not found.")

    mars_data = get_celestial_body_data("mars")  # Test case-insensitivity
    if mars_data:
        print("\nMars Data:")
        for key, value in mars_data.items():
            print(f"  {key}: {value}")
    else:
        print("\nMars data not found.")

    moon_data = get_celestial_body_data("Moon")
    if moon_data:
        print("\nMoon Data:")
        for key, value in moon_data.items():
            print(f"  {key}: {value}")
    else:
        print("\nMoon data not found.")

    # Test for a non-existent body
    unknown_data = get_celestial_body_data("Pluto")
    if unknown_data:
        print("\nPluto Data (unexpectedly found):")
        for key, value in unknown_data.items():
            print(f"  {key}: {value}")
    else:
        print("\nPluto data not found (as expected for this dataset).")

    # Test getting all solar system destinations in the format compatible with main.py
    solar_system_destinations_list = get_all_solar_system_destinations()
    print("\n--- Solar System Destinations (for main.py integration) ---")
    if solar_system_destinations_list:
        for dest in solar_system_destinations_list:
            # Format distance for readability in kilometers
            print(f"  Name: {dest['name']}, Distance from Earth: {dest['distance']/1e3:.2f} km")
    else:
        print("No solar system destinations found.")

    # Test error handling for invalid input type
    invalid_input_data = get_celestial_body_data(123)
    if invalid_input_data is None:
        print("\nSuccessfully handled invalid input type (expected None).")
    else:
        print("\nFailed to handle invalid input type.")