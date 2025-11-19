# orbital_mechanics.py

import math
# from datetime import datetime, timedelta # Not directly used in this module after refactoring
from constants import G_GRAVITATIONAL  # Importing G_GRAVITATIONAL directly
from celestial_data import get_celestial_body_data  # Importing the specific function from celestial_data

# --- 3D Vector Utility Functions ---

def _validate_vector(vec: tuple, name: str):
    """Internal helper to validate a 3D vector."""
    if not isinstance(vec, (tuple, list)) or len(vec) != 3:
        raise ValueError(f"{name} must be a 3-element tuple or list of numbers.")
    if not all(isinstance(x, (int, float)) for x in vec):
        raise ValueError(f"All elements in {name} must be numeric.")

def vector_magnitude(vec: tuple) -> float:
    """
    Calculates the Euclidean magnitude of a 3D vector.

    Args:
        vec (tuple): A 3-element tuple representing the vector (x, y, z).

    Returns:
        float: The magnitude of the vector.

    Raises:
        ValueError: If the input is not a 3-element tuple of numbers.
    """
    _validate_vector(vec, "vector")
    return math.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)

def vector_subtract(vec1: tuple, vec2: tuple) -> tuple:
    """
    Subtracts one 3D vector from another (vec1 - vec2).

    Args:
        vec1 (tuple): The first 3-element vector.
        vec2 (tuple): The second 3-element vector.

    Returns:
        tuple: A new 3-element vector (x, y, z) representing the result.

    Raises:
        ValueError: If inputs are not 3-element tuples of numbers.
    """
    _validate_vector(vec1, "vec1")
    _validate_vector(vec2, "vec2")
    return (vec1[0] - vec2[0], vec1[1] - vec2[1], vec1[2] - vec2[2])

def vector_add(vec1: tuple, vec2: tuple) -> tuple:
    """
    Adds two 3D vectors (vec1 + vec2).

    Args:
        vec1 (tuple): The first 3-element vector.
        vec2 (tuple): The second 3-element vector.

    Returns:
        tuple: A new 3-element vector (x, y, z) representing the result.

    Raises:
        ValueError: If inputs are not 3-element tuples of numbers.
    """
    _validate_vector(vec1, "vec1")
    _validate_vector(vec2, "vec2")
    return (vec1[0] + vec2[0], vec1[1] + vec2[1], vec1[2] + vec2[2])

def vector_scale(vec: tuple, scalar: float) -> tuple:
    """
    Scales a 3D vector by a scalar.

    Args:
        vec (tuple): A 3-element vector.
        scalar (float): The scalar value to multiply the vector by.

    Returns:
        tuple: A new 3-element vector (x, y, z) representing the scaled result.

    Raises:
        ValueError: If the vector is not a 3-element tuple of numbers, or scalar is not numeric.
    """
    _validate_vector(vec, "vector")
    if not isinstance(scalar, (int, float)):
        raise ValueError("Scalar must be a numeric value.")
    return (vec[0] * scalar, vec[1] * scalar, vec[2] * scalar)

def vector_dot(vec1: tuple, vec2: tuple) -> float:
    """
    Calculates the dot product of two 3D vectors.

    Args:
        vec1 (tuple): The first 3-element vector.
        vec2 (tuple): The second 3-element vector.

    Returns:
        float: The dot product.

    Raises:
        ValueError: If inputs are not 3-element tuples of numbers.
    """
    _validate_vector(vec1, "vec1")
    _validate_vector(vec2, "vec2")
    return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2]

def vector_cross(vec1: tuple, vec2: tuple) -> tuple:
    """
    Calculates the cross product of two 3D vectors.

    Args:
        vec1 (tuple): The first 3-element vector.
        vec2 (tuple): The second 3-element vector.

    Returns:
        tuple: A new 3-element vector (x, y, z) representing the cross product.

    Raises:
        ValueError: If inputs are not 3-element tuples of numbers.
    """
    _validate_vector(vec1, "vec1")
    _validate_vector(vec2, "vec2")
    x = vec1[1] * vec2[2] - vec1[2] * vec2[1]
    y = vec1[2] * vec2[0] - vec1[0] * vec2[2]
    z = vec1[0] * vec2[1] - vec1[1] * vec2[0]
    return (x, y, z)

# --- Existing Orbital Mechanics Functions ---

def calculate_gravitational_parameter(central_body_name: str) -> float:
    """
    Retrieves the standard gravitational parameter (mu) for a given central celestial body.
    This function first attempts to fetch 'gravitational_parameter_mu' directly from the
    celestial_data module. If not found, it falls back to calculating mu using the body's
    mass and the universal gravitational constant G_GRAVITATIONAL.

    Args:
        central_body_name (str): The name of the central body (e.g., "Sun", "Earth", "Mars").

    Returns:
        float: The gravitational parameter (mu) in m^3/s^2.

    Raises:
        ValueError: If the central body's data cannot be found, if 'gravitational_parameter_mu'
                    is missing or non-positive, or if 'mass' is missing or non-positive
                    when 'gravitational_parameter_mu' is not available.
    """
    body_data = get_celestial_body_data(central_body_name)
    if body_data is None:
        raise ValueError(f"Celestial body '{central_body_name}' not found in celestial data.")

    # Prioritize 'gravitational_parameter_mu' if available
    mu = body_data.get('gravitational_parameter_mu')
    if mu is not None:
        if not isinstance(mu, (int, float)) or mu <= 0:
            raise ValueError(
                f"Gravitational parameter 'gravitational_parameter_mu' for '{central_body_name}' "
                f"must be a positive number, got {mu}."
            )
        return float(mu)
    
    # Fallback: calculate mu from mass if 'gravitational_parameter_mu' is not directly provided
    mass = body_data.get('mass')
    if mass is None:
        raise ValueError(
            f"Gravitational parameter ('gravitational_parameter_mu') or mass data "
            f"missing for '{central_body_name}' in celestial data."
        )
    
    if not isinstance(mass, (int, float)) or mass <= 0:
        raise ValueError(f"Mass for '{central_body_name}' must be a positive number, got {mass}.")
        
    return G_GRAVITATIONAL * mass

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

# --- Universal Variable Lambert Solver Helper Functions ---

def _stumpff_c(z: float) -> float:
    """Calculates the Stumpff function C(z)."""
    if z > 1e-6:
        return (1 - math.cos(math.sqrt(z))) / z
    elif z < -1e-6:
        return (math.cosh(math.sqrt(-z)) - 1) / (-z)
    else: # z is close to 0
        return 0.5 # Taylor expansion for z=0

def _stumpff_s(z: float) -> float:
    """Calculates the Stumpff function S(z)."""
    if z > 1e-6:
        return (math.sqrt(z) - math.sin(math.sqrt(z))) / (z * math.sqrt(z))
    elif z < -1e-6:
        return (math.sinh(math.sqrt(-z)) - math.sqrt(-z)) / (-z * math.sqrt(-z))
    else: # z is close to 0
        return 1/6 # Taylor expansion for z=0

def _calc_time_from_psi(psi: float, A_param: float, mu: float) -> float:
    """
    Calculates the time of flight for Lambert's problem given psi (universal anomaly squared).
    This function implements the time equation (F(psi) in a Newton-Raphson context).
    Based on Vallado's Algorithm 50, `tof_func`.
    """
    sqrt_mu = math.sqrt(mu)
    x_val = math.sqrt(abs(psi)) # Universal anomaly (x in Vallado)

    C = _stumpff_c(psi)
    S = _stumpff_s(psi)

    # Vallado's Algorithm 50, tof_func (p. 288)
    time_val = (x_val**3 * S + A_param * x_val * (1 - C)) / sqrt_mu
    
    return time_val

def _calc_d_time_d_psi(psi: float, A_param: float, mu: float) -> float:
    """
    Calculates the derivative of time of flight with respect to psi (universal anomaly squared).
    Used in Newton-Raphson iteration.
    Based on Vallado's Algorithm 50, `d_tof_d_x` transformed to `d_tof_d_psi`.
    """
    sqrt_mu = math.sqrt(mu)
    x_val = math.sqrt(abs(psi))

    C = _stumpff_c(psi)
    S = _stumpff_s(psi)

    # Vallado's Algorithm 50, d_tof_d_x (derivative w.r.t. `x`, universal anomaly)
    d_tof_d_x = (x_val**2 * S + A_param * (1 - C)) / sqrt_mu

    # Using the chain rule: dt/d_psi = (dt/d_x) * (dx/d_psi)
    # where x = sqrt(psi), so dx/d_psi = 1 / (2 * sqrt(psi)) = 1 / (2*x_val)
    if x_val < 1e-10: # If x_val is practically zero, derivative is infinite.
        return float('inf')
        
    return d_tof_d_x / (2 * x_val)


def calculate_lambert_transfer(mu: float, r1_vec: tuple, r2_vec: tuple, delta_t: float) -> tuple[float, float, float, float]:
    """
    Calculates the Delta-V for a Lambert transfer between two position vectors in a given time of flight.
    This function implements a full iterative Lambert solver using the universal variable formulation
    and Newton-Raphson iteration to determine the transfer orbit.

    The solver finds the initial and final velocity vectors (v1_transfer_vec, v2_transfer_vec)
    required to traverse from r1_vec to r2_vec in delta_t seconds. The returned Delta-V magnitudes
    are the magnitudes of these transfer velocities relative to the central body.

    Args:
        mu (float): Gravitational parameter of the central body (m^3/s^2).
        r1_vec (tuple): Initial position vector (x, y, z) in meters from the central body.
        r2_vec (tuple): Final position vector (x, y, z) in meters from the central body.
        delta_t (float): Time of flight in seconds.

    Returns:
        tuple[float, float, float, float]: A tuple containing:
            - delta_v1_mag (float): Magnitude of Delta-V for the first burn (m/s).
            - delta_v2_mag (float): Magnitude of Delta-V for the second burn (m/s).
            - total_delta_v (float): Total Delta-V (sum of magnitudes) for both burns (m/s).
            - calculated_time_of_flight (float): The actual delta_t used for calculations (will be the input delta_t if successful).

    Raises:
        ValueError: If inputs are invalid (e.g., non-positive mu/delta_t, invalid vectors,
                    transfer impossible, or solver fails to converge).
    """
    if mu <= 0:
        raise ValueError("Gravitational parameter (mu) must be positive.")
    if delta_t <= 0:
        raise ValueError("Time of flight (delta_t) must be positive.")
    _validate_vector(r1_vec, "r1_vec")
    _validate_vector(r2_vec, "r2_vec")

    r1_mag = vector_magnitude(r1_vec)
    r2_mag = vector_magnitude(r2_vec)

    if r1_mag <= 0 or r2_mag <= 0:
        raise ValueError("Position vector magnitudes must be positive.")
    if r1_vec == r2_vec:
        return 0.0, 0.0, 0.0, delta_t # No transfer needed

    # Calculate geometrical parameters
    cos_nu = vector_dot(r1_vec, r2_vec) / (r1_mag * r2_mag)
    
    # Handle floating point inaccuracies
    cos_nu = max(-1.0, min(1.0, cos_nu))
    
    nu = math.acos(cos_nu) # Angle between r1 and r2
    
    # A parameter (constant for a given geometry, for prograde short-way transfer)
    # This `A_param` is Vallado's `A` (Alg 50, p 287), which is `sqrt(r1*r2)*cos(nu/2)`
    A_param = math.sqrt(r1_mag * r2_mag) * math.cos(nu / 2)
    
    # --- Newton-Raphson Iteration to find psi (universal anomaly squared) ---
    # `psi` (universal anomaly squared, `z` in Vallado)
    
    # Initial guess for psi (z). A small positive value often works for typical transfers.
    psi_k = 0.1 

    max_iterations = 100
    tolerance = 1e-8 # Tolerance for time difference

    for _ in range(max_iterations):
        time_k = _calc_time_from_psi(psi_k, A_param, mu)
        
        if time_k < 0 and psi_k > 0: # If time is negative for elliptic solutions, means psi_k is too high
            psi_k = psi_k * 0.5 # Try a smaller psi
            continue
        elif time_k > delta_t and psi_k < 0: # If time is too large for hyperbolic, means psi_k is too low
            psi_k = psi_k * 2.0 # Try a larger psi
            continue

        time_prime_k = _calc_d_time_d_psi(psi_k, A_param, mu)
        
        if time_prime_k == 0 or math.isinf(time_prime_k) or math.isnan(time_prime_k):
            raise ValueError("Lambert solver: Derivative of time is invalid, cannot converge with Newton-Raphson.")

        # Check for convergence
        if abs(time_k - delta_t) < tolerance:
            break

        # Newton-Raphson step
        psi_k_new = psi_k - (time_k - delta_t) / time_prime_k
        
        # Simple clamping/adjustment to prevent wild jumps and ensure progress towards convergence
        if abs(psi_k_new) > 1e10: 
            psi_k = psi_k * 0.5 # Reduce step if it's too large
        else:
            psi_k = psi_k_new
    else:
        raise ValueError(f"Lambert solver failed to converge within {max_iterations} iterations for delta_t={delta_t:.2f}s.")

    final_psi = psi_k
    chi_final = math.sqrt(abs(final_psi)) # Universal anomaly (x in Vallado)
    
    calculated_delta_t = _calc_time_from_psi(final_psi, A_param, mu)
    if abs(calculated_delta_t - delta_t) > tolerance:
         raise ValueError(f"Lambert solver converged but time error {abs(calculated_delta_t - delta_t):.2e}s exceeds tolerance.")

    # --- Calculate transfer velocities (v1_transfer_vec, v2_transfer_vec) ---
    sqrt_mu = math.sqrt(mu)
    C_final = _stumpff_c(final_psi)
    S_final = _stumpff_s(final_psi)
    
    # Using the standard F, G, F_dot, G_dot coefficients (from Bate, Mueller, White Eq 5.4.10 and 5.4.11, p 195)
    # and relation: v = (r_final - F*r_initial) / G, v_final = (G_dot*r_final - r_initial) / G
    
    f_val = 1 - (final_psi * C_final) / r1_mag
    g_val = delta_t - (chi_final**3 * S_final) / sqrt_mu
    g_dot_val = 1 - (final_psi * C_final) / r2_mag
    
    if g_val == 0:
        raise ValueError("Denominator `g` for Lambert transfer velocity calculation is zero. Transfer is physically impossible (e.g., singular case).")
    
    v1_transfer_vec = vector_scale(vector_subtract(r2_vec, vector_scale(r1_vec, f_val)), 1/g_val)
    v2_transfer_vec = vector_scale(vector_subtract(vector_scale(r2_vec, g_dot_val), r1_vec), 1/g_val)

    delta_v1_mag = vector_magnitude(v1_transfer_vec)
    delta_v2_mag = vector_magnitude(v2_transfer_vec)
    total_delta_v = delta_v1_mag + delta_v2_mag

    return delta_v1_mag, delta_v2_mag, total_delta_v, delta_t