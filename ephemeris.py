import math
import datetime
from celestial_data import get_celestial_body_data
from orbital_mechanics import calculate_gravitational_parameter # Only mu is needed from orbital_mechanics

# --- Constants ---
# J2000.0 epoch (January 1, 2000, 12:00 TT). Using UTC for simplicity in calculations.
J2000_EPOCH = datetime.datetime(2000, 1, 1, 12, 0, 0, tzinfo=datetime.timezone.utc)

# --- Internal Helper for Sun's Gravitational Parameter ---
_MU_SUN_CACHE = None # Initialized at module level as requested

def _get_mu_sun() -> float:
    """
    Retrieves the Sun's gravitational parameter (mu) from celestial data,
    caching the result for subsequent calls.

    Returns:
        float: The gravitational parameter of the Sun in m^3/s^2.

    Raises:
        ValueError: If the Sun's gravitational parameter cannot be determined.
    """
    global _MU_SUN_CACHE # Placed as the very first executable statement as requested
    if _MU_SUN_CACHE is None:
        try:
            # Prefer using orbital_mechanics function which handles retrieval logic
            _MU_SUN_CACHE = calculate_gravitational_parameter('Sun')
        except ValueError as e:
            # Fallback directly to celestial_data if orbital_mechanics fails for some reason
            sun_data = get_celestial_body_data('Sun')
            if sun_data and 'gravitational_parameter_mu' in sun_data:
                _MU_SUN_CACHE = sun_data['gravitational_parameter_mu']
            else:
                raise ValueError(f"Failed to retrieve Sun's gravitational parameter: {e}. "
                                 "Check 'Sun' entry in celestial data.")
        if _MU_SUN_CACHE is None or not isinstance(_MU_SUN_CACHE, (int, float)) or _MU_SUN_CACHE <= 0:
            raise ValueError("Sun's gravitational parameter could not be determined or is invalid.")
    return _MU_SUN_CACHE

# --- NEW HELPER FUNCTIONS FOR 3D ORBITAL MECHANICS ---

def _normalize_angle(angle_rad: float) -> float:
    """Normalizes an angle to be within [0, 2*pi)."""
    return angle_rad % (2 * math.pi)

def _calculate_mean_anomaly(M0_rad: float, mean_motion_rad_s: float, delta_t_s: float) -> float:
    """
    Calculates the mean anomaly at a given time relative to the J2000 epoch.

    Args:
        M0_rad (float): Mean anomaly at epoch (J2000) in radians.
        mean_motion_rad_s (float): Mean motion in radians per second.
        delta_t_s (float): Time difference from epoch in seconds.

    Returns:
        float: Mean anomaly in radians, normalized to [0, 2*pi).
    """
    M = M0_rad + mean_motion_rad_s * delta_t_s
    return _normalize_angle(M)

def _solve_keplers_equation(M_rad: float, e: float, tolerance: float = 1e-10, max_iterations: int = 100) -> float:
    """
    Solves Kepler's Equation E - e * sin(E) = M for Eccentric Anomaly (E) using Newton-Raphson iteration.

    Args:
        M_rad (float): Mean anomaly in radians.
        e (float): Eccentricity. Must be 0 <= e < 1 for elliptical orbits.
        tolerance (float): Desired accuracy for E.
        max_iterations (int): Maximum number of iterations for the solver.

    Returns:
        float: Eccentric Anomaly (E) in radians.

    Raises:
        ValueError: If the eccentricity is invalid or the solver fails to converge.
    """
    if not (0 <= e < 1):
        # Allow e=0 for circular orbits where E=M, but raise for other invalid eccentricities
        if not math.isclose(e, 0.0, abs_tol=tolerance):
             raise ValueError(f"Eccentricity {e} is not valid for elliptical orbits (must be 0 <= e < 1).")

    # Initial guess for E, often M itself is a good starting point.
    E = M_rad

    for _ in range(max_iterations):
        f = E - e * math.sin(E) - M_rad
        fp = 1 - e * math.cos(E) # Derivative of f with respect to E
        
        # Avoid division by zero if derivative is too small
        if math.isclose(fp, 0.0, abs_tol=tolerance):
            raise ValueError("Kepler's Equation solver encountered a near-zero derivative, indicating potential "
                             "convergence issues or invalid orbital parameters (e.g., e close to 1 near periapsis).")
        
        E_new = E - f / fp
        if abs(E_new - E) < tolerance:
            return E_new
        E = E_new
    raise ValueError(f"Kepler's Equation solver failed to converge after {max_iterations} iterations "
                     f"for Mean Anomaly={M_rad:.4f} rad and Eccentricity={e:.4f}.")

def _calculate_true_anomaly_and_radius(a: float, e: float, E_rad: float) -> tuple[float, float]:
    """
    Calculates the true anomaly and orbital radius from semi-major axis, eccentricity,
    and eccentric anomaly.

    Args:
        a (float): Semi-major axis in meters. Must be positive.
        e (float): Eccentricity.
        E_rad (float): Eccentric Anomaly in radians.

    Returns:
        tuple[float, float]: (true_anomaly_rad, radius_m)
    
    Raises:
        ValueError: If semi-major axis is non-positive.
    """
    if a <= 0:
        raise ValueError(f"Semi-major axis '{a}' must be positive.")

    # Orbital radius
    r = a * (1 - e * math.cos(E_rad))

    # True anomaly using atan2 for correct quadrant handling
    nu = math.atan2(math.sqrt(1 - e*e) * math.sin(E_rad), math.cos(E_rad) - e)
    return _normalize_angle(nu), r

def _elements_to_cartesian(mu: float, a: float, e: float, i_rad: float, Omega_rad: float, omega_rad: float, nu_rad: float, r_m: float) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    """
    Converts classical orbital elements to heliocentric Cartesian position and velocity vectors
    in the J2000 equatorial inertial frame.

    Args:
        mu (float): Gravitational parameter of the central body (Sun) in m^3/s^2.
        a (float): Semi-major axis in meters.
        e (float): Eccentricity.
        i_rad (float): Inclination in radians.
        Omega_rad (float): Longitude of ascending node (RAAN) in radians.
        omega_rad (float): Argument of periapsis in radians.
        nu_rad (float): True anomaly in radians.
        r_m (float): Orbital radius in meters.

    Returns:
        tuple[tuple[float, float, float], tuple[float, float, float]]:
            (position_vector_m, velocity_vector_mps)
            - position_vector_m (tuple): (x, y, z) position in meters relative to the Sun.
            - velocity_vector_mps (tuple): (vx, vy, vz) velocity in meters/second relative to the Sun.

    Raises:
        ValueError: If semi-latus rectum is non-positive due to invalid orbital elements.
    """
    # Position in perifocal frame (P, Q)
    x_perifocal = r_m * math.cos(nu_rad)
    y_perifocal = r_m * math.sin(nu_rad)
    # z_perifocal is 0 in the perifocal frame as it's the orbital plane.

    # Velocity in perifocal frame
    p = a * (1 - e*e) # Semi-latus rectum
    if p <= 0:
        raise ValueError("Calculated semi-latus rectum is non-positive. Invalid orbital elements (a or e).")
    
    sqrt_mu_div_p = math.sqrt(mu / p)
    vx_perifocal = -sqrt_mu_div_p * math.sin(nu_rad)
    vy_perifocal = sqrt_mu_div_p * (e + math.cos(nu_rad))
    # vz_perifocal is 0 in the perifocal frame.

    # Pre-calculate trigonometric values for rotation
    cos_Omega = math.cos(Omega_rad)
    sin_Omega = math.sin(Omega_rad)
    cos_i = math.cos(i_rad)
    sin_i = math.sin(i_rad)
    cos_omega = math.cos(omega_rad)
    sin_omega = math.sin(omega_rad)

    # Transformation from perifocal (P, Q, W) to J2000 equatorial inertial (I, J, K) frame
    # This matrix is R_Z(Omega) * R_X(i) * R_Z(omega)
    # Reference: Vallado, Fundamentals of Astrodynamics and Applications, 4th ed., page 110 (Eq. 2-71)
    
    # Position components
    r_I = x_perifocal * (cos_omega * cos_Omega - sin_omega * sin_Omega * cos_i) \
        - y_perifocal * (sin_omega * cos_Omega + cos_omega * sin_Omega * cos_i)

    r_J = x_perifocal * (cos_omega * sin_Omega + sin_omega * cos_Omega * cos_i) \
        + y_perifocal * (cos_omega * cos_Omega * cos_i - sin_omega * sin_Omega)

    r_K = x_perifocal * (sin_omega * sin_i) \
        + y_perifocal * (cos_omega * sin_i)

    position_vector_m = (r_I, r_J, r_K)

    # Velocity components (same rotation matrix applied)
    v_I = vx_perifocal * (cos_omega * cos_Omega - sin_omega * sin_Omega * cos_i) \
        - vy_perifocal * (sin_omega * cos_Omega + cos_omega * sin_Omega * cos_i)

    v_J = vx_perifocal * (cos_omega * sin_Omega + sin_omega * cos_Omega * cos_i) \
        + vy_perifocal * (cos_omega * cos_Omega * cos_i - sin_omega * sin_Omega)

    v_K = vx_perifocal * (sin_omega * sin_i) \
        + vy_perifocal * (cos_omega * sin_i)

    velocity_vector_mps = (v_I, v_J, v_K)

    return position_vector_m, velocity_vector_mps


def get_heliocentric_state(body_name: str, date: datetime.datetime) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    """
    Calculates the approximate heliocentric orbital state (position and velocity vectors)
    for a given celestial body at a specified date and time using a 3D elliptical orbit model.

    This function explicitly does NOT support:
    - The Sun itself, as it is the central body.
    - Earth's Moon, as its heliocentric orbit is more complex (a barycentric orbit
      around the Sun with Earth) and requires a full N-body or patched conic approach,
      which is beyond this model's scope.

    Args:
        body_name (str): The case-insensitive name of the celestial body (e.g., 'Mars', 'Jupiter').
                         Must be a non-empty string.
        date (datetime.datetime): The specific date and time for which to calculate
                                  the orbital state. It is highly recommended to provide
                                  timezone-aware datetime objects, preferably in UTC,
                                  for consistent results. If the datetime object is
                                  timezone-naive, it will be assumed to be in local time
                                  and converted to UTC before calculation.

    Returns:
        tuple[tuple[float, float, float], tuple[float, float, float]]:
            A tuple containing two 3-element tuples:
            - position_vector_m (tuple): (x, y, z) position in meters relative to the Sun.
            - velocity_vector_mps (tuple): (vx, vy, vz) velocity in meters/second relative to the Sun.

    Raises:
        ValueError: If inputs are invalid, required celestial body data is missing or invalid,
                    or the calculation cannot be performed for the specified body.
    """
    if not isinstance(body_name, str) or not body_name.strip():
        raise ValueError("`body_name` must be a non-empty string.")
    if not isinstance(date, datetime.datetime):
        raise ValueError("`date` must be a datetime.datetime object.")

    body_name_lower = body_name.lower()
    if body_name_lower == 'sun':
        raise ValueError("Cannot calculate heliocentric state for the Sun itself (it is the central body).")
    if body_name_lower == 'moon' or body_name_lower == "earth's moon":
        raise ValueError("Heliocentric state for Earth's Moon is not directly supported by this simplified model. "
                         "It requires calculating Earth's position first.")

    # Get Sun's gravitational parameter
    mu_sun = _get_mu_sun()

    # Retrieve celestial body data
    body_data = get_celestial_body_data(body_name)
    if not body_data:
        raise ValueError(f"Data for celestial body '{body_name}' not found.")

    # Validate and extract required orbital elements
    # ASSUMPTION: The 'celestial_data.py' module now provides these extended orbital elements.
    required_elements = [
        'semi_major_axis_from_sun', 'eccentricity', 'inclination_rad',
        'longitude_of_ascending_node_rad', 'argument_of_periapsis_rad',
        'mean_anomaly_at_J2000_rad', 'orbital_period'
    ]
    
    for element in required_elements:
        if body_data.get(element) is None:
            raise ValueError(f"Missing orbital element '{element}' for '{body_name}' in celestial data. "
                             "Cannot perform 3D elliptical orbit calculation.")
        if not isinstance(body_data.get(element), (int, float)):
            raise ValueError(f"Invalid type for orbital element '{element}' for '{body_name}'. Must be numeric.")

    a = body_data['semi_major_axis_from_sun']
    e = body_data['eccentricity']
    i_rad = body_data['inclination_rad']
    Omega_rad = body_data['longitude_of_ascending_node_rad']
    omega_rad = body_data['argument_of_periapsis_rad']
    M0_rad = body_data['mean_anomaly_at_J2000_rad']
    orbital_period_s = body_data['orbital_period']

    # Basic validation of retrieved elements
    if a <= 0:
        raise ValueError(f"Invalid semi-major axis '{a}' for '{body_name}'. Must be positive.")
    if not (0 <= e < 1): # For elliptical orbits, e must be less than 1
         if not math.isclose(e, 0.0, abs_tol=1e-10): # Allow e=0 for circular
            raise ValueError(f"Invalid eccentricity '{e}' for '{body_name}'. Must be 0 <= e < 1 for elliptical orbits.")
    if orbital_period_s <= 0:
        raise ValueError(f"Invalid orbital period '{orbital_period_s}' for '{body_name}'. Must be positive.")

    # Standardize input date to UTC for consistent time difference calculation
    if date.tzinfo is None:
        # If naive, assume it's local time and convert to UTC
        date_utc = date.astimezone(datetime.timezone.utc)
    else:
        # If timezone-aware, convert to UTC
        date_utc = date.astimezone(datetime.timezone.utc)

    # Calculate time difference from J2000 epoch in seconds
    delta_t_s = (date_utc - J2000_EPOCH).total_seconds()

    # 1. Calculate Mean Motion
    mean_motion_rad_s = (2 * math.pi) / orbital_period_s

    # 2. Calculate Mean Anomaly at `date`
    M_rad = _calculate_mean_anomaly(M0_rad, mean_motion_rad_s, delta_t_s)

    # 3. Solve Kepler's Equation for Eccentric Anomaly
    E_rad = _solve_keplers_equation(M_rad, e)

    # 4. Calculate True Anomaly and Orbital Radius
    nu_rad, r_m = _calculate_true_anomaly_and_radius(a, e, E_rad)

    # 5. Convert Orbital Elements to Cartesian Position and Velocity Vectors
    position_vector_m, velocity_vector_mps = _elements_to_cartesian(mu_sun, a, e, i_rad, Omega_rad, omega_rad, nu_rad, r_m)

    return position_vector_m, velocity_vector_mps

if __name__ == '__main__':
    print("--- Ephemeris Module Self-Test (3D Elliptical Model) ---")

    # --- MOCKING CELESTIAL DATA FOR SELF-TEST ---
    # This class temporarily replaces the celestial_data module's get_celestial_body_data
    # to provide the necessary extended orbital elements for testing.
    class MockCelestialData:
        def get_celestial_body_data(self, name: str):
            # Values are illustrative and approximate for J2000 epoch.
            # RAAN and Arg of Peri often precess, so these are snapshots.
            data = {
                'sun': {
                    'gravitational_parameter_mu': 1.32712440018e20, # m^3/s^2
                },
                'earth': {
                    'mass': 5.972e24, 'radius': 6.371e6,
                    'gravitational_parameter_mu': 3.986004418e14,
                    'semi_major_axis_from_sun': 149.598e9, # m (approx 1 AU)
                    'orbital_period': 31557600.0, # s (approx 1 sidereal year)
                    'eccentricity': 0.0167086,
                    'inclination_rad': math.radians(0.00005), # Close to 0 for ecliptic, J2000 eqatorial is 0.0
                    'longitude_of_ascending_node_rad': math.radians(-11.26064), # J2000 RAAN
                    'argument_of_periapsis_rad': math.radians(102.93735), # J2000 Arg of Periapsis
                    'mean_anomaly_at_J2000_rad': math.radians(357.51716), # J2000 Mean Anomaly
                },
                'mars': {
                    'mass': 6.4171e23, 'radius': 3.3895e6,
                    'gravitational_parameter_mu': 4.2828375e13,
                    'semi_major_axis_from_sun': 227.9392e9, # m
                    'orbital_period': 59354000.0, # s (approx 687 Earth days)
                    'eccentricity': 0.0934,
                    'inclination_rad': math.radians(1.85040), # Relative to ecliptic for Mars
                    'longitude_of_ascending_node_rad': math.radians(49.55743), # J2000 RAAN
                    'argument_of_periapsis_rad': math.radians(286.5016), # J2000 Arg of Periapsis
                    'mean_anomaly_at_J2000_rad': math.radians(18.6021), # J2000 Mean Anomaly
                },
                'jupiter': {
                    'mass': 1.898e27, 'radius': 69.911e6,
                    'gravitational_parameter_mu': 1.26686534e17,
                    'semi_major_axis_from_sun': 778.57e9, # m
                    'orbital_period': 374330000.0, # s (approx 11.86 Earth years)
                    'eccentricity': 0.0489,
                    'inclination_rad': math.radians(1.3033), # Relative to ecliptic for Jupiter
                    'longitude_of_ascending_node_rad': math.radians(100.46444), # J2000 RAAN
                    'argument_of_periapsis_rad': math.radians(273.8643), # J2000 Arg of Periapsis
                    'mean_anomaly_at_J2000_rad': math.radians(20.020), # J2000 Mean Anomaly
                }
            }
            return data.get(name.lower())

    # Temporarily patch the celestial_data.get_celestial_body_data function for testing
    # This is a basic patch for self-contained testing; in a real pytest environment, use unittest.mock.patch.
    import sys
    sys_module_backup = sys.modules.get('celestial_data')
    mock_celestial_data_instance = MockCelestialData()
    sys.modules['celestial_data'] = mock_celestial_data_instance
    # Re-import get_celestial_body_data to ensure it uses the mocked module
    from celestial_data import get_celestial_body_data as _get_celestial_body_data_mocked
    # Overwrite the function reference in the current module's scope
    globals()['get_celestial_body_data'] = _get_celestial_body_data_mocked
    # Also clear the cached Sun mu to force recalculation with mock data
    global _MU_SUN_CACHE
    _MU_SUN_CACHE = None
    print("Celestial data temporarily mocked for self-test.")
    # --- END MOCKING ---


    test_date_now = datetime.datetime.now(datetime.timezone.utc)
    print(f"\nCurrent Test Date (UTC): {test_date_now}")

    # --- Test valid celestial bodies ---
    print("\nTesting valid celestial bodies (Earth, Mars, Jupiter):")
    bodies_to_test_valid = ['Earth', 'Mars', 'Jupiter']
    for body in bodies_to_test_valid:
        try:
            pos_vec, vel_vec = get_heliocentric_state(body, test_date_now)
            print(f"  {body}:")
            print(f"    Position Vector (m): ({pos_vec[0]:.2e}, {pos_vec[1]:.2e}, {pos_vec[2]:.2e})")
            print(f"    Velocity Vector (m/s): ({vel_vec[0]:.2e}, {vel_vec[1]:.2e}, {vel_vec[2]:.2e})")
            print(f"    Distance from Sun (m): {math.sqrt(sum(x**2 for x in pos_vec)):.2e}")
            print(f"    Speed (m/s): {math.sqrt(sum(x**2 for x in vel_vec)):.2e}")
        except ValueError as e:
            print(f"  Error calculating state for {body}: {e}")
        print("-" * 30)

    # --- Test invalid/edge case inputs ---
    print("\nTesting invalid/edge case inputs:")

    # 1. Invalid body_name type
    try:
        get_heliocentric_state(123, test_date_now)
    except ValueError as e:
        print(f"  Caught expected error for invalid body_name type (123): {e}")

    # 2. Empty body_name string
    try:
        get_heliocentric_state("", test_date_now)
    except ValueError as e:
        print(f"  Caught expected error for empty body_name string (''): {e}")

    # 3. Non-existent celestial body (will now raise error due to missing orbital elements)
    try:
        get_heliocentric_state("Krypton", test_date_now)
    except ValueError as e:
        print(f"  Caught expected error for non-existent body ('Krypton'): {e}")

    # 4. Sun itself (central body)
    try:
        get_heliocentric_state("Sun", test_date_now)
    except ValueError as e:
        print(f"  Caught expected error for 'Sun' as the body: {e}")

    # 5. Earth's Moon (explicitly not supported by this simplified model)
    try:
        get_heliocentric_state("Moon", test_date_now)
    except ValueError as e:
        print(f"  Caught expected error for 'Moon' as the body: {e}")
    try:
        get_heliocentric_state("Earth's Moon", test_date_now)
    except ValueError as e:
        print(f"  Caught expected error for 'Earth's Moon' as the body: {e}")

    # 6. Invalid date type
    try:
        get_heliocentric_state("Earth", "this is not a date string")
    except ValueError as e:
        print(f"  Caught expected error for invalid date type ('not a date string'): {e}")

    # --- Test with a past date (before J2000) ---
    past_date = datetime.datetime(1998, 5, 15, 6, 0, 0, tzinfo=datetime.timezone.utc)
    print(f"\nTesting with a past date (UTC): {past_date}")
    try:
        pos_vec_past, vel_vec_past = get_heliocentric_state('Earth', past_date)
        print(f"  Earth (past date):")
        print(f"    Position Vector (m): ({pos_vec_past[0]:.2e}, {pos_vec_past[1]:.2e}, {pos_vec_past[2]:.2e})")
        print(f"    Velocity Vector (m/s): ({vel_vec_past[0]:.2e}, {vel_vec_past[1]:.2e}, {vel_vec_past[2]:.2e})")
    except ValueError as e:
        print(f"  Error calculating state for Earth (past date): {e}")

    # --- Test with a timezone-naive datetime (should be handled by converting to UTC) ---
    unaware_date = datetime.datetime(2023, 10, 26, 10, 30, 0) # No tzinfo
    print(f"\nTesting with a timezone-naive datetime (local time assumed then converted to UTC): {unaware_date}")
    try:
        pos_vec_unaware, vel_vec_unaware = get_heliocentric_state('Mars', unaware_date)
        print(f"  Mars (from unaware date):")
        print(f"    Position Vector (m): ({pos_vec_unaware[0]:.2e}, {pos_vec_unaware[1]:.2e}, {pos_vec_unaware[2]:.2e})")
        print(f"    Velocity Vector (m/s): ({vel_vec_unaware[0]:.2e}, {vel_vec_unaware[1]:.2e}, {vel_vec_unaware[2]:.2e})")
    except ValueError as e:
        print(f"  Error calculating state for Mars (unaware date): {e}")
    
    # --- Cleanup Mocking ---
    if sys_module_backup:
        sys.modules['celestial_data'] = sys_module_backup
    else:
        del sys.modules['celestial_data']
    print("\nCleaned up celestial data mock.")
    # --- End Cleanup ---

    print("\n--- Ephemeris Module Self-Test Complete ---")