import math
import datetime
from celestial_data import get_celestial_body_data
from orbital_mechanics import calculate_gravitational_parameter, calculate_elliptical_orbital_period

# --- Constants ---
# J2000.0 epoch (January 1, 2000, 12:00 TT). Using UTC for simplicity in calculations.
J2000_EPOCH = datetime.datetime(2000, 1, 1, 12, 0, 0, tzinfo=datetime.timezone.utc)

# --- Internal Helper for Sun's Gravitational Parameter ---
_MU_SUN_CACHE = None

def _get_mu_sun() -> float:
    """
    Retrieves the Sun's gravitational parameter (mu) from celestial data,
    caching the result for subsequent calls.

    Returns:
        float: The gravitational parameter of the Sun in m^3/s^2.

    Raises:
        ValueError: If the Sun's gravitational parameter cannot be determined.
    """
    global _MU_SUN_CACHE
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

def get_heliocentric_state(body_name: str, date: datetime.datetime) -> dict:
    """
    Calculates the approximate heliocentric orbital state (distance from the Sun and
    angular position) for a given celestial body at a specified date and time.
    This function uses a simplified circular orbit model based on the body's
    semi-major axis from the Sun.

    This function explicitly does NOT support:
    - The Sun itself, as it is the central body.
    - Earth's Moon, as its heliocentric orbit is more complex (a barycentric orbit
      around the Sun with Earth) and requires a two-body calculation (Earth's
      position + Moon's geocentric position), which is beyond this simplified model's scope.

    Args:
        body_name (str): The case-insensitive name of the celestial body (e.g., 'Mars', 'Jupiter').
        date (datetime.datetime): The specific date and time for which to calculate
                                  the orbital state. It is highly recommended to provide
                                  timezone-aware datetime objects, preferably in UTC,
                                  for consistent results. If timezone-naive, it will be
                                  converted to UTC assuming local time.

    Returns:
        dict: A dictionary containing the approximate heliocentric state:
            - "distance_from_sun_m" (float): The approximate distance from the Sun in meters.
            - "angular_position_rad" (float): The approximate angular position in radians
              (ranging from 0 to 2*pi), measured counter-clockwise from an arbitrary
              reference direction at the J2000.0 epoch.

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

    # Extract semi-major axis (used as approximate orbital radius for a circular orbit)
    semi_major_axis_from_sun = body_data.get('semi_major_axis_from_sun')
    if semi_major_axis_from_sun is None:
        raise ValueError(f"Missing 'semi_major_axis_from_sun' for '{body_name}' in celestial data. "
                         "Cannot calculate heliocentric state with this simplified model.")
    if not isinstance(semi_major_axis_from_sun, (int, float)) or semi_major_axis_from_sun <= 0:
        raise ValueError(f"Invalid 'semi_major_axis_from_sun' value for '{body_name}'. "
                         "Must be a positive numeric value.")

    # Calculate orbital period using Kepler's Third Law
    try:
        orbital_period_s = calculate_elliptical_orbital_period(mu_sun, semi_major_axis_from_sun)
    except ValueError as e:
        raise ValueError(f"Failed to calculate orbital period for '{body_name}': {e}")

    if orbital_period_s <= 0:
        raise ValueError(f"Calculated orbital period for '{body_name}' is non-positive ({orbital_period_s} s). "
                         "Cannot determine angular position.")

    # Calculate angular speed (radians per second)
    angular_speed_rad_s = (2 * math.pi) / orbital_period_s

    # Standardize input date to UTC for consistent time difference calculation
    if date.tzinfo is None:
        # If naive, assume it's local time and convert to UTC
        date_utc = date.astimezone(datetime.timezone.utc)
    else:
        # If timezone-aware, convert to UTC
        date_utc = date.astimezone(datetime.timezone.utc)

    # Calculate time difference from J2000 epoch in seconds
    time_difference_s = (date_utc - J2000_EPOCH).total_seconds()

    # Calculate approximate angular position.
    # We assume an initial angular position of 0 radians at J2000.0 for all bodies
    # for this simplified model. This means the result is relative to an arbitrary
    # reference direction at the epoch.
    angular_position_rad = (time_difference_s * angular_speed_rad_s) % (2 * math.pi)

    # Ensure angular position is within [0, 2*pi)
    if angular_position_rad < 0:
        angular_position_rad += 2 * math.pi

    return {
        "distance_from_sun_m": semi_major_axis_from_sun,
        "angular_position_rad": angular_position_rad,
    }

if __name__ == '__main__':
    print("--- Ephemeris Module Self-Test ---")

    # Current time in UTC for testing
    test_date_now = datetime.datetime.now(datetime.timezone.utc)
    print(f"\nCurrent Test Date (UTC): {test_date_now}")

    # --- Test valid celestial bodies ---
    print("\nTesting valid celestial bodies (Earth, Mars, Jupiter):")
    bodies_to_test_valid = ['Earth', 'Mars', 'Jupiter']
    for body in bodies_to_test_valid:
        try:
            state = get_heliocentric_state(body, test_date_now)
            print(f"  {body}:")
            print(f"    Distance from Sun: {state['distance_from_sun_m']:.2e} m")
            print(f"    Angular Position: {math.degrees(state['angular_position_rad']):.2f}°")
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

    # 3. Non-existent celestial body
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
        state_past = get_heliocentric_state('Earth', past_date)
        print(f"  Earth (past date):")
        print(f"    Distance from Sun: {state_past['distance_from_sun_m']:.2e} m")
        print(f"    Angular Position: {math.degrees(state_past['angular_position_rad']):.2f}°")
    except ValueError as e:
        print(f"  Error calculating state for Earth (past date): {e}")

    # --- Test with a timezone-naive datetime (should be handled by converting to UTC) ---
    unaware_date = datetime.datetime(2023, 10, 26, 10, 30, 0) # No tzinfo
    print(f"\nTesting with a timezone-naive datetime (local time assumed then converted to UTC): {unaware_date}")
    try:
        state_unaware = get_heliocentric_state('Mars', unaware_date)
        print(f"  Mars (from unaware date):")
        print(f"    Distance from Sun: {state_unaware['distance_from_sun_m']:.2e} m")
        print(f"    Angular Position: {math.degrees(state_unaware['angular_position_rad']):.2f}°")
    except ValueError as e:
        print(f"  Error calculating state for Mars (unaware date): {e}")

    print("\n--- Ephemeris Module Self-Test Complete ---")