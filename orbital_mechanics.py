import math

# Universal Gravitational Constant (m^3 kg^-1 s^-2)
G = 6.67430e-11

# The `celestial_data` module is listed as a dependency but is not provided in the
# repository context. In a real application, this would be an actual module
# providing planetary masses, radii, or pre-calculated gravitational parameters.
# For this file, functions will assume that necessary parameters (like mass
# or gravitational parameter) are passed directly to them.
# Example if celestial_data existed:
# from celestial_data import get_body_mass, get_standard_gravitational_parameter


def calculate_gravitational_parameter(mass_central_body: float) -> float:
    """
    Calculates the standard gravitational parameter (mu) for a central body.

    The standard gravitational parameter (mu) is the product of the universal
    gravitational constant (G) and the mass (M) of the central body.
    It is fundamental for many orbital mechanics calculations.

    Args:
        mass_central_body (float): The mass of the central body in kilograms (kg).

    Returns:
        float: The standard gravitational parameter (mu) in cubic meters per
               second squared (m^3/s^2).

    Raises:
        TypeError: If `mass_central_body` is not a numeric type.
        ValueError: If `mass_central_body` is not a positive value.
    """
    if not isinstance(mass_central_body, (float, int)):
        raise TypeError("mass_central_body must be a numeric value.")
    if mass_central_body <= 0:
        raise ValueError("mass_central_body must be a positive value.")

    return G * mass_central_body


def calculate_orbital_period(semimajor_axis: float, gravitational_parameter: float) -> float:
    """
    Calculates the orbital period of an object in an elliptical orbit using Kepler's Third Law.

    This function is applicable for elliptical and circular orbits where the
    semi-major axis is positive.

    Args:
        semimajor_axis (float): The semi-major axis of the orbit in meters (m).
        gravitational_parameter (float): The standard gravitational parameter (mu)
                                         of the central body in m^3/s^2.

    Returns:
        float: The orbital period in seconds (s).

    Raises:
        TypeError: If `semimajor_axis` or `gravitational_parameter` are not numeric types.
        ValueError: If `semimajor_axis` or `gravitational_parameter` are not positive,
                    or if the calculation results in an invalid mathematical operation
                    (e.g., square root of a negative number, though handled by positive checks).
    """
    if not isinstance(semimajor_axis, (float, int)):
        raise TypeError("semimajor_axis must be a numeric value.")
    if not isinstance(gravitational_parameter, (float, int)):
        raise TypeError("gravitational_parameter must be a numeric value.")
    if semimajor_axis <= 0:
        raise ValueError("semimajor_axis must be a positive value for orbital period calculation.")
    if gravitational_parameter <= 0:
        raise ValueError("gravitational_parameter must be a positive value.")

    try:
        # Kepler's Third Law: T = 2 * pi * sqrt(a^3 / mu)
        period = 2 * math.pi * math.sqrt(semimajor_axis**3 / gravitational_parameter)
        return period
    except Exception as e:
        # Catching potential math errors like sqrt of negative if inputs bypassed checks
        raise ValueError(f"Error calculating orbital period: {e}")


def calculate_specific_orbital_energy(semimajor_axis: float, gravitational_parameter: float) -> float:
    """
    Calculates the specific orbital energy (epsilon) of an object in orbit.

    Specific orbital energy is a constant for a given orbit and is used to classify
    orbit types:
    - Epsilon < 0: Elliptical orbit (including circular)
    - Epsilon = 0: Parabolic orbit
    - Epsilon > 0: Hyperbolic orbit

    Args:
        semimajor_axis (float): The semi-major axis of the orbit in meters (m).
                                Note: For hyperbolic orbits, `semimajor_axis` is negative.
                                For parabolic orbits, `semimajor_axis` is considered infinite.
        gravitational_parameter (float): The standard gravitational parameter (mu)
                                         of the central body in m^3/s^2.

    Returns:
        float: The specific orbital energy in joules per kilogram (J/kg) or m^2/s^2.

    Raises:
        TypeError: If `semimajor_axis` or `gravitational_parameter` are not numeric types.
        ValueError: If `gravitational_parameter` is not positive, or if `semimajor_axis`
                    is zero (which would lead to division by zero and is physically
                    unrealistic for a stable orbit definition in this context).
    """
    if not isinstance(semimajor_axis, (float, int)):
        raise TypeError("semimajor_axis must be a numeric value.")
    if not isinstance(gravitational_parameter, (float, int)):
        raise TypeError("gravitational_parameter must be a numeric value.")
    if gravitational_parameter <= 0:
        raise ValueError("gravitational_parameter must be a positive value.")
    if semimajor_axis == 0:
        raise ValueError("Semi-major axis cannot be zero for specific orbital energy calculation (implies collision).")

    # Epsilon = -mu / (2 * a)
    return -gravitational_parameter / (2 * semimajor_axis)


def calculate_orbital_velocity_circular(radius: float, gravitational_parameter: float) -> float:
    """
    Calculates the orbital velocity required for a stable circular orbit at a given radius.

    Args:
        radius (float): The radius of the circular orbit in meters (m). This is
                        measured from the center of the central body.
        gravitational_parameter (float): The standard gravitational parameter (mu)
                                         of the central body in m^3/s^2.

    Returns:
        float: The orbital velocity in meters per second (m/s).

    Raises:
        TypeError: If `radius` or `gravitational_parameter` are not numeric types.
        ValueError: If `radius` or `gravitational_parameter` are not positive.
    """
    if not isinstance(radius, (float, int)):
        raise TypeError("radius must be a numeric value.")
    if not isinstance(gravitational_parameter, (float, int)):
        raise TypeError("gravitational_parameter must be a numeric value.")
    if radius <= 0:
        raise ValueError("radius must be a positive value.")
    if gravitational_parameter <= 0:
        raise ValueError("gravitational_parameter must be a positive value.")

    try:
        # v = sqrt(mu / r)
        velocity = math.sqrt(gravitational_parameter / radius)
        return velocity
    except Exception as e:
        raise ValueError(f"Error calculating circular orbital velocity: {e}")


def calculate_delta_v_hohmann_transfer(radius_initial: float, radius_final: float, gravitational_parameter: float) -> tuple[float, float, float]:
    """
    Calculates the delta-V required for a Hohmann transfer between two circular orbits.

    A Hohmann transfer is an elliptical orbit used as an energy-efficient maneuver
    to transfer a spacecraft between two co-planar circular orbits of different
    altitudes around a central body.

    Args:
        radius_initial (float): The radius of the initial circular orbit in meters (m).
        radius_final (float): The radius of the final circular orbit in meters (m).
        gravitational_parameter (float): The standard gravitational parameter (mu)
                                         of the central body in m^3/s^2.

    Returns:
        tuple[float, float, float]: A tuple containing (delta_v1, delta_v2, total_delta_v)
                                    in meters per second (m/s).
                                    - `delta_v1`: The magnitude of the first burn, required
                                                  to inject the spacecraft into the transfer orbit.
                                    - `delta_v2`: The magnitude of the second burn, required
                                                  to circularize the orbit at the final radius.
                                    - `total_delta_v`: The sum of `delta_v1` and `delta_v2`,
                                                       representing the total propellant cost.

    Raises:
        TypeError: If any input (`radius_initial`, `radius_final`, `gravitational_parameter`)
                   is not a numeric type.
        ValueError: If any radius or `gravitational_parameter` are not positive,
                    or if the initial and final radii are identical (no transfer needed).
    """
    if not isinstance(radius_initial, (float, int)):
        raise TypeError("radius_initial must be a numeric value.")
    if not isinstance(radius_final, (float, int)):
        raise TypeError("radius_final must be a numeric value.")
    if not isinstance(gravitational_parameter, (float, int)):
        raise TypeError("gravitational_parameter must be a numeric value.")

    if radius_initial <= 0 or radius_final <= 0:
        raise ValueError("Orbit radii must be positive values.")
    if gravitational_parameter <= 0:
        raise ValueError("gravitational_parameter must be a positive value.")
    if radius_initial == radius_final:
        raise ValueError("Initial and final radii cannot be the same for a Hohmann transfer; no transfer is needed.")

    try:
        # 1. Calculate velocities for initial and final circular orbits
        v_initial_circular = calculate_orbital_velocity_circular(radius_initial, gravitational_parameter)
        v_final_circular = calculate_orbital_velocity_circular(radius_final, gravitational_parameter)

        # 2. Calculate semi-major axis of the Hohmann transfer ellipse
        semimajor_axis_transfer = (radius_initial + radius_final) / 2

        # 3. Calculate velocities at periapsis and apoapsis of the transfer orbit
        # Velocity in an elliptical orbit: v = sqrt(mu * ((2/r) - (1/a)))
        v_periapsis_transfer = math.sqrt(gravitational_parameter * ((2 / radius_initial) - (1 / semimajor_axis_transfer)))
        v_apoapsis_transfer = math.sqrt(gravitational_parameter * ((2 / radius_final) - (1 / semimajor_axis_transfer)))

        # 4. Calculate delta-V for the first burn (injection into transfer orbit)
        delta_v1 = abs(v_periapsis_transfer - v_initial_circular)

        # 5. Calculate delta-V for the second burn (circularization at final orbit)
        delta_v2 = abs(v_final_circular - v_apoapsis_transfer)
        
        total_delta_v = delta_v1 + delta_v2

        return delta_v1, delta_v2, total_delta_v
    except ValueError as ve:
        # Re-raise ValueErrors from helper functions with additional context
        raise ValueError(f"Error during Hohmann transfer calculation: {ve}")
    except Exception as e:
        # Catch any other unexpected errors
        raise ValueError(f"An unexpected error occurred during Hohmann transfer calculation: {e}")