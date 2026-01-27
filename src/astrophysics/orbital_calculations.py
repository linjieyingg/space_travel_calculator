import math
from typing import Dict, Any, Optional

# Relative import to celestial_data.py located at the repository root.
# This assumes 'src' is a package and 'orbital_calculations.py' is within a subpackage 'astrophysics'.
# For this to work correctly, the top-level directory (where celestial_data.py and src are)
# must be treated as a Python package root, or the 'src' directory (and thus 'src/astrophysics')
# must be run in a way that its parent directory is in sys.path.
from ..celestial_data import get_celestial_body_data, GRAVITATIONAL_CONSTANT


def solve_kepler_equation(
    mean_anomaly: float,
    eccentricity: float,
    tolerance: float = 1e-8,
    max_iterations: int = 100
) -> float:
    """
    Solves Kepler's Equation (M = E - e * sin(E)) for the eccentric anomaly (E)
    using the Newton-Raphson method.

    This function is crucial for determining a celestial body's position in an
    elliptical orbit at a given time.

    Args:
        mean_anomaly (float): The mean anomaly in radians (M). It represents
                              the angle a fictitious body would have if it
                              moved in a circular orbit with the same period
                              as the actual body.
        eccentricity (float): The orbital eccentricity (e). For elliptical
                              orbits, 0 <= e < 1.
        tolerance (float, optional): The convergence criterion for the
                                     Newton-Raphson method. Defaults to 1e-8.
        max_iterations (int, optional): The maximum number of iterations
                                        allowed for convergence. Defaults to 100.

    Returns:
        float: The eccentric anomaly (E) in radians.

    Raises:
        ValueError: If eccentricity is out of the valid range for elliptical
                    orbits (0 <= e < 1), or if the equation does not converge
                    within the specified maximum iterations.
        TypeError: If inputs are not of the expected numeric types.
    """
    if not isinstance(mean_anomaly, (int, float)) or \
       not isinstance(eccentricity, (int, float)):
        raise TypeError("Mean anomaly and eccentricity must be numeric types.")
    if not isinstance(tolerance, (int, float)) or \
       not isinstance(max_iterations, int):
        raise TypeError("Tolerance must be numeric and max_iterations an integer.")

    if not (0 <= eccentricity < 1):
        raise ValueError("Eccentricity must be between 0 (inclusive) and 1 (exclusive) for elliptical orbits.")

    # Normalize mean_anomaly to be within [0, 2*pi]
    mean_anomaly = math.fmod(mean_anomaly, 2 * math.pi)
    if mean_anomaly < 0:
        mean_anomaly += 2 * math.pi

    # Initial guess for Eccentric Anomaly (E)
    # A common and simple initial guess is M.
    eccentric_anomaly = mean_anomaly

    for _ in range(max_iterations):
        f = eccentric_anomaly - eccentricity * math.sin(eccentric_anomaly) - mean_anomaly
        f_prime = 1 - eccentricity * math.cos(eccentric_anomaly)

        # f_prime should generally be positive for e < 1, preventing division by zero.
        # min value of f_prime is 1-e (at E=0 or E=2pi).
        # Since e < 1, 1-e > 0.
        delta_e = f / f_prime
        eccentric_anomaly -= delta_e

        if abs(delta_e) < tolerance:
            return eccentric_anomaly

    raise ValueError(
        f"Kepler's Equation failed to converge within {max_iterations} iterations."
    )


def eccentric_to_true_anomaly(eccentric_anomaly: float, eccentricity: float) -> float:
    """
    Converts the eccentric anomaly to the true anomaly.

    The true anomaly is the angle between the periapsis direction and the
    current position of the orbiting body, as seen from the focus of the ellipse.

    Args:
        eccentric_anomaly (float): The eccentric anomaly in radians (E).
        eccentricity (float): The orbital eccentricity (e). For elliptical
                              orbits, 0 <= e < 1.

    Returns:
        float: The true anomaly (nu) in radians, normalized to the range [0, 2*pi).

    Raises:
        ValueError: If eccentricity is out of the valid range for elliptical
                    orbits (0 <= e < 1).
        TypeError: If inputs are not of the expected numeric types.
    """
    if not isinstance(eccentric_anomaly, (int, float)) or \
       not isinstance(eccentricity, (int, float)):
        raise TypeError("Eccentric anomaly and eccentricity must be numeric types.")

    if not (0 <= eccentricity < 1):
        raise ValueError("Eccentricity must be between 0 (inclusive) and 1 (exclusive) for elliptical orbits.")

    # Using math.atan2 is robust as it handles all quadrants correctly.
    # sin(nu) = (sqrt(1 - e^2) * sin(E)) / (1 - e * cos(E))
    # cos(nu) = (cos(E) - e) / (1 - e * cos(E))
    # The denominator `(1 - e * cos(E))` is always positive for e < 1,
    # so no division by zero issues are expected.
    y_component = math.sqrt(1 - eccentricity**2) * math.sin(eccentric_anomaly)
    x_component = math.cos(eccentric_anomaly) - eccentricity

    true_anomaly = math.atan2(y_component, x_component)

    # Normalize true_anomaly to [0, 2*pi)
    if true_anomaly < 0:
        true_anomaly += 2 * math.pi

    return true_anomaly


def calculate_orbital_distance(
    body_name: str, time_since_periapsis: float
) -> Optional[float]:
    """
    Calculates the current orbital distance of a celestial body from its primary
    (assumed to be the Sun for this application) at a given time since periapsis.

    This function integrates data retrieval, Kepler's equation solving, and
    orbital geometry to provide a position estimate.

    Args:
        body_name (str): The case-insensitive name of the celestial body.
        time_since_periapsis (float): The time in seconds that has elapsed
                                      since the body last passed its periapsis.
                                      Must be non-negative.

    Returns:
        Optional[float]: The orbital distance in meters if calculation is successful,
                         otherwise None if the body is not found.

    Raises:
        ValueError: If `time_since_periapsis` is negative, or if required
                    orbital data (semi-major axis, eccentricity, orbital period)
                    is missing or invalid for the specified body, or if Kepler's
                    equation fails to converge.
        TypeError: If inputs are not of the expected types.
    """
    if not isinstance(body_name, str):
        raise TypeError("Body name must be a string.")
    if not isinstance(time_since_periapsis, (int, float)):
        raise TypeError("Time since periapsis must be a numeric type.")
    if time_since_periapsis < 0:
        raise ValueError("Time since periapsis cannot be negative.")

    body_data = get_celestial_body_data(body_name)

    if body_data is None:
        return None  # Body not found

    # Special case for the Sun: it is the central body, its distance from itself is 0.
    if body_data.get('name', '').lower() == 'sun':
        return 0.0

    # Extract required orbital elements
    required_keys = ['semi_major_axis_m', 'eccentricity', 'orbital_period_s']
    for key in required_keys:
        if key not in body_data or not isinstance(body_data[key], (int, float)):
            raise ValueError(
                f"Missing or invalid orbital data '{key}' for {body_name}. "
                "Cannot calculate orbital distance."
            )

    semi_major_axis = body_data['semi_major_axis_m']
    eccentricity = body_data['eccentricity']
    orbital_period = body_data['orbital_period_s']

    if semi_major_axis <= 0:
        raise ValueError(
            f"Semi-major axis for {body_name} must be positive for an orbit."
        )
    if not (0 <= eccentricity < 1):
        raise ValueError(
            f"Eccentricity for {body_name} must be between 0 (inclusive) and 1 (exclusive) for elliptical orbits."
        )
    if orbital_period <= 0:
        raise ValueError(
            f"Orbital period for {body_name} must be positive."
        )

    # Calculate Mean Motion (n) in radians/second
    mean_motion = (2 * math.pi) / orbital_period

    # Calculate Mean Anomaly (M) at the given time
    mean_anomaly = mean_motion * time_since_periapsis

    # Solve Kepler's Equation for Eccentric Anomaly (E)
    try:
        eccentric_anomaly = solve_kepler_equation(mean_anomaly, eccentricity)
    except ValueError as e:
        raise ValueError(f"Failed to solve Kepler's equation for {body_name}: {e}")
    except TypeError as e:
        raise TypeError(f"Invalid input type for Kepler's equation for {body_name}: {e}")

    # Calculate the orbital distance (r) from the central body (Sun)
    # The formula is r = a * (1 - e * cos(E))
    orbital_distance = semi_major_axis * (1 - eccentricity * math.cos(eccentric_anomaly))

    return orbital_distance


def calculate_hohmann_delta_v(
    initial_body_name: str, final_body_name: str, central_body_name: str
) -> Dict[str, float]:
    """
    Calculates the two delta-v impulses and the total delta-v required for a
    Hohmann transfer maneuver between two circular (or near-circular) orbits
    around a central body.

    A Hohmann transfer is an elliptical orbit used to transfer between two
    circular orbits. It involves two impulsive burns: one to enter the
    transfer ellipse (delta_v1) and another to circularize the orbit at the
    destination (delta_v2).

    Assumptions:
    - Initial and final orbits are circular or approximated as such by their semi-major axes.
    - Orbits are coplanar.
    - Impulsive burns.

    Args:
        initial_body_name (str): The name of the celestial body in the initial orbit.
        final_body_name (str): The name of the celestial body in the final (target) orbit.
        central_body_name (str): The name of the central celestial body around which
                                 the transfer occurs (e.g., "Sun").

    Returns:
        Dict[str, float]: A dictionary containing:
            - 'delta_v1': The magnitude of the first delta-v impulse (m/s).
            - 'delta_v2': The magnitude of the second delta-v impulse (m/s).
            - 'total_delta_v': The sum of delta_v1 and delta_v2 (m/s).

    Raises:
        TypeError: If any input `body_name` is not a string.
        ValueError: If:
            - Any of the specified celestial bodies are not found in the data.
            - Required orbital parameters (semi-major axis, gravitational parameter)
              are missing or invalid for the specified bodies.
            - The initial or final body is the same as the central body.
            - Semi-major axes or central body gravitational parameter are non-positive.
            - Division by zero occurs in intermediate calculations (e.g., due to invalid semi-major axis).
    """
    if not all(isinstance(arg, str) for arg in [initial_body_name, final_body_name, central_body_name]):
        raise TypeError("All body names must be strings.")

    # If initial and final bodies are the same, no transfer is needed, delta-v is zero.
    if initial_body_name.lower() == final_body_name.lower():
        return {"delta_v1": 0.0, "delta_v2": 0.0, "total_delta_v": 0.0}

    # Fetch data for all three bodies
    initial_body_data = get_celestial_body_data(initial_body_name)
    final_body_data = get_celestial_body_data(final_body_name)
    central_body_data = get_celestial_body_data(central_body_name)

    if initial_body_data is None:
        raise ValueError(f"Celestial body data not found for initial body: {initial_body_name}")
    if final_body_data is None:
        raise ValueError(f"Celestial body data not found for final body: {final_body_name}")
    if central_body_data is None:
        raise ValueError(f"Celestial body data not found for central body: {central_body_name}")

    # Ensure initial and final bodies are not the central body
    if initial_body_data.get('name', '').lower() == central_body_name.lower():
        raise ValueError(f"Initial body '{initial_body_name}' cannot be the central body '{central_body_name}'.")
    if final_body_data.get('name', '').lower() == central_body_name.lower():
        raise ValueError(f"Final body '{final_body_name}' cannot be the central body '{central_body_name}'.")

    # Extract gravitational parameter of the central body (mu = G * M_central)
    mu_central = central_body_data.get('gravitational_parameter_mu')
    if mu_central is None or not isinstance(mu_central, (int, float)) or mu_central <= 0:
        raise ValueError(
            f"Missing or invalid 'gravitational_parameter_mu' for central body {central_body_name}. "
            "Cannot calculate Hohmann transfer (must be positive)."
        )

    # Extract semi-major axes for initial and final orbits.
    # For Hohmann transfer, these are treated as radii of the initial/final circular orbits.
    r1 = initial_body_data.get('semi_major_axis_m')
    r2 = final_body_data.get('semi_major_axis_m')

    if r1 is None or not isinstance(r1, (int, float)) or r1 <= 0:
        raise ValueError(
            f"Missing or invalid 'semi_major_axis_m' for initial body {initial_body_name}. "
            "Cannot calculate Hohmann transfer (r1 must be positive)."
        )
    if r2 is None or not isinstance(r2, (int, float)) or r2 <= 0:
        raise ValueError(
            f"Missing or invalid 'semi_major_axis_m' for final body {final_body_name}. "
            "Cannot calculate Hohmann transfer (r2 must be positive)."
        )

    # Calculate velocities for circular orbits (v = sqrt(mu/r))
    v1_circular = math.sqrt(mu_central / r1)
    v2_circular = math.sqrt(mu_central / r2)

    # Calculate the semi-major axis of the transfer ellipse (a_transfer = (r1 + r2) / 2)
    a_transfer = (r1 + r2) / 2

    if a_transfer <= 0:  # This check ensures validity for the vis-viva equation
        raise ValueError("Calculated semi-major axis of transfer ellipse is non-positive, indicating invalid orbital radii.")

    # Calculate velocities at periapsis (r1) and apoapsis (r2) of the transfer ellipse
    # using the Vis-Viva equation: v = sqrt(mu * ((2/r) - (1/a)))
    try:
        v_periapsis_transfer = math.sqrt(mu_central * ((2 / r1) - (1 / a_transfer)))
        v_apoapsis_transfer = math.sqrt(mu_central * ((2 / r2) - (1 / a_transfer)))
    except ValueError: # math.sqrt raises ValueError for negative input
        raise ValueError(
            "Orbital parameters lead to physically impossible velocities "
            "(e.g., negative value under square root). Check input radii."
        )

    # Calculate the delta-v impulses
    # delta_v1: difference between transfer ellipse periapsis velocity and initial circular orbit velocity
    delta_v1 = abs(v_periapsis_transfer - v1_circular)
    # delta_v2: difference between final circular orbit velocity and transfer ellipse apoapsis velocity
    delta_v2 = abs(v2_circular - v_apoapsis_transfer)

    total_delta_v = delta_v1 + delta_v2

    return {
        "delta_v1": delta_v1,
        "delta_v2": delta_v2,
        "total_delta_v": total_delta_v
    }