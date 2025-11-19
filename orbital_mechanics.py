import math
# Assuming constants.py and celestial_bodies.py are available in the Python path
# or within the same directory as this file, following the project's import patterns.
import constants
import celestial_bodies

# Gravitational Constant (G) in m^3 kg^-1 s^-2, assumed to be defined in constants.py
# e.g., constants.GRAVITATIONAL_CONSTANT = 6.67430e-11
G = constants.GRAVITATIONAL_CONSTANT

def calculate_gravitational_parameter(central_body_name: str) -> float:
    """
    Calculates the standard gravitational parameter (mu) for a given central celestial body.
    The gravitational parameter is a constant for a given body, equal to the product of
    the gravitational constant (G) and the body's mass (M).

    Args:
        central_body_name (str): The name of the central body (e.g., "Sun", "Earth", "Mars").
                                 The mass of this body is retrieved from the `celestial_bodies` module.

    Returns:
        float: The gravitational parameter (mu) in m^3/s^2.

    Raises:
        ValueError: If the central body's mass cannot be found in the `celestial_bodies` data,
                    or if the mass is non-positive.
    """
    body_mass = celestial_bodies.get_body_mass(central_body_name)
    if body_mass is None:
        raise ValueError(f"Mass for central body '{central_body_name}' not found in celestial_bodies data.")
    if body_mass <= 0:
        raise ValueError(f"Mass for central body '{central_body_name}' must be positive.")
    return G * body_mass

def calculate_circular_orbital_velocity(mu: float, radius: float) -> float:
    """
    Calculates the orbital velocity required for a stable circular orbit at a given radius
    around a central body with a specified gravitational parameter.

    Args:
        mu (float): The gravitational parameter of the central body (m^3/s^2).
        radius (float): The radius of the circular orbit from the center of the central body (m).

    Returns:
        float: The orbital velocity in m/s.

    Raises:
        ValueError: If the gravitational parameter (mu) or orbital radius is non-positive.
    """
    if mu <= 0:
        raise ValueError("Gravitational parameter (mu) must be positive.")
    if radius <= 0:
        raise ValueError("Orbital radius must be positive.")
    return math.sqrt(mu / radius)

def calculate_hohmann_transfer_delta_v(mu: float, r1: float, r2: float) -> tuple[float, float, float]:
    """
    Calculates the Delta-V (change in velocity) required for a Hohmann transfer
    between two coplanar circular orbits. A Hohmann transfer is a two-impulse
    maneuver used to transfer a spacecraft between two orbits.

    Args:
        mu (float): Gravitational parameter of the central body (m^3/s^2).
        r1 (float): Radius of the initial circular orbit (m).
        r2 (float): Radius of the final circular orbit (m).

    Returns:
        tuple[float, float, float]: A tuple containing:
            - delta_v1 (float): Delta-V for the first burn to enter the transfer orbit (m/s).
                                Positive for acceleration, negative for deceleration.
            - delta_v2 (float): Delta-V for the second burn to enter the final orbit (m/s).
                                Positive for acceleration, negative for deceleration.
            - total_delta_v (float): Absolute sum of Delta-V for both burns (m/s),
                                     representing the total fuel cost.

    Raises:
        ValueError: If mu, r1, or r2 are non-positive, or if r1 equals r2 (no transfer needed).
    """
    if mu <= 0:
        raise ValueError("Gravitational parameter (mu) must be positive.")
    if r1 <= 0 or r2 <= 0:
        raise ValueError("Orbital radii (r1, r2) must be positive.")
    if r1 == r2:
        raise ValueError("Initial and final orbital radii cannot be equal for a Hohmann transfer.")

    # Semi-major axis of the elliptical transfer orbit
    a_transfer = (r1 + r2) / 2

    # Velocities in the initial and final circular orbits
    v_circular_r1 = calculate_circular_orbital_velocity(mu, r1)
    v_circular_r2 = calculate_circular_orbital_velocity(mu, r2)

    # Velocities at r1 and r2 within the transfer ellipse
    v_transfer_r1 = math.sqrt(mu * ((2 / r1) - (1 / a_transfer)))
    v_transfer_r2 = math.sqrt(mu * ((2 / r2) - (1 / a_transfer)))

    # Delta-V for the first burn (at r1)
    delta_v1 = v_transfer_r1 - v_circular_r1

    # Delta-V for the second burn (at r2)
    delta_v2 = v_circular_r2 - v_transfer_r2

    # Total Delta-V is the sum of the absolute values of the two burns
    total_delta_v = abs(delta_v1) + abs(delta_v2)

    return delta_v1, delta_v2, total_delta_v

def calculate_hohmann_transfer_time_of_flight(mu: float, r1: float, r2: float) -> float:
    """
    Calculates the time of flight for a Hohmann transfer. This is half the orbital period
    of the elliptical transfer orbit.

    Args:
        mu (float): Gravitational parameter of the central body (m^3/s^2).
        r1 (float): Radius of the initial circular orbit (m).
        r2 (float): Radius of the final circular orbit (m).

    Returns:
        float: The time of flight in seconds.

    Raises:
        ValueError: If mu, r1, or r2 are non-positive, or if r1 equals r2.
    """
    if mu <= 0:
        raise ValueError("Gravitational parameter (mu) must be positive.")
    if r1 <= 0 or r2 <= 0:
        raise ValueError("Orbital radii (r1, r2) must be positive.")
    if r1 == r2:
        raise ValueError("Initial and final orbital radii cannot be equal for a Hohmann transfer.")

    # Semi-major axis of the elliptical transfer orbit
    a_transfer = (r1 + r2) / 2

    # Calculate the orbital period of the transfer ellipse using Kepler's Third Law
    period_transfer = 2 * math.pi * math.sqrt(a_transfer**3 / mu)

    # Time of flight is half the period of the transfer ellipse
    time_of_flight = period_transfer / 2

    return time_of_flight

def calculate_elliptical_orbital_period(mu: float, semi_major_axis: float) -> float:
    """
    Calculates the orbital period of an elliptical (or circular) orbit using Kepler's Third Law.

    Args:
        mu (float): Gravitational parameter of the central body (m^3/s^2).
        semi_major_axis (float): The semi-major axis of the elliptical orbit (m).
                                 For a circular orbit, this is simply the orbit's radius.

    Returns:
        float: The orbital period in seconds.

    Raises:
        ValueError: If the gravitational parameter (mu) or semi-major axis is non-positive.
    """
    if mu <= 0:
        raise ValueError("Gravitational parameter (mu) must be positive.")
    if semi_major_axis <= 0:
        raise ValueError("Semi-major axis must be positive.")

    return 2 * math.pi * math.sqrt(semi_major_axis**3 / mu)

def calculate_escape_velocity(mu: float, radius: float) -> float:
    """
    Calculates the escape velocity from a central body at a given radius.
    Escape velocity is the minimum speed needed for a free, non-propelled object
    to escape from the gravitational influence of a massive body.

    Args:
        mu (float): The gravitational parameter of the central body (m^3/s^2).
        radius (float): The distance from the center of the central body (m).

    Returns:
        float: The escape velocity in m/s.

    Raises:
        ValueError: If the gravitational parameter (mu) or radius is non-positive.
    """
    if mu <= 0:
        raise ValueError("Gravitational parameter (mu) must be positive.")
    if radius <= 0:
        raise ValueError("Radius must be positive.")

    return math.sqrt(2 * mu / radius)