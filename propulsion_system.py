import math

# Constants
STANDARD_GRAVITY = 9.80665  # m/s^2, standard acceleration due to gravity


def calculate_required_fuel_mass(
    delta_v: float,
    dry_mass: float,
    specific_impulse: float = None,
    exhaust_velocity: float = None
) -> float:
    """
    Calculates the required fuel mass using the Tsiolkovsky rocket equation.

    This function determines the amount of fuel needed for a spacecraft to achieve
    a specific change in velocity (delta-V), given its dry mass and engine
    characteristics (either specific impulse or exhaust velocity).

    Args:
        delta_v (float): The required change in velocity (delta-V) in meters/second (m/s).
                         Must be non-negative.
        dry_mass (float): The mass of the spacecraft without fuel (payload + structure)
                          in kilograms (kg). Must be positive.
        specific_impulse (float, optional): The engine's specific impulse in seconds (s).
                                            Must be positive if exhaust_velocity is not provided.
        exhaust_velocity (float, optional): The effective exhaust velocity of the engine
                                            in meters/second (m/s). Must be positive if
                                            specific_impulse is not provided.

    Returns:
        float: The required fuel mass in kilograms (kg).

    Raises:
        ValueError: If `delta_v` is negative.
        ValueError: If `dry_mass` is not positive.
        ValueError: If neither `specific_impulse` nor `exhaust_velocity` is provided.
        ValueError: If `specific_impulse` is provided but not positive.
        ValueError: If `exhaust_velocity` is provided but not positive.
        ValueError: If calculated effective exhaust velocity is zero or negative, which would
                    lead to mathematical errors.
    """

    # Input validation
    if not isinstance(delta_v, (int, float)) or delta_v < 0:
        raise ValueError("Delta-V must be a non-negative numeric value.")
    if not isinstance(dry_mass, (int, float)) or dry_mass <= 0:
        raise ValueError("Dry mass must be a positive numeric value.")

    effective_exhaust_velocity: float = 0.0

    if exhaust_velocity is not None:
        if not isinstance(exhaust_velocity, (int, float)) or exhaust_velocity <= 0:
            raise ValueError("Exhaust velocity must be a positive numeric value.")
        effective_exhaust_velocity = exhaust_velocity
    elif specific_impulse is not None:
        if not isinstance(specific_impulse, (int, float)) or specific_impulse <= 0:
            raise ValueError("Specific impulse must be a positive numeric value.")
        effective_exhaust_velocity = specific_impulse * STANDARD_GRAVITY
    else:
        raise ValueError("Either specific_impulse or exhaust_velocity must be provided.")

    if effective_exhaust_velocity <= 0:
        raise ValueError(
            "Calculated effective exhaust velocity cannot be zero or negative. "
            "Check specific_impulse or exhaust_velocity input."
        )

    # If delta_v is 0, no change in velocity is required, thus no fuel is needed.
    if delta_v == 0:
        return 0.0

    # Tsiolkovsky rocket equation: delta_v = Ve * ln(m0 / mf)
    # Where:
    #   Ve = effective_exhaust_velocity
    #   m0 = initial_mass (dry_mass + fuel_mass)
    #   mf = dry_mass (final mass after fuel is expended)
    #
    # We want to find fuel_mass = m0 - mf
    # Rearranging the equation to solve for the mass ratio (m0 / mf):
    # m0 / mf = exp(delta_v / Ve)
    #
    # Then, m0 = mf * exp(delta_v / Ve)
    # And fuel_mass = (mf * exp(delta_v / Ve)) - mf

    mass_ratio = math.exp(delta_v / effective_exhaust_velocity)
    initial_mass = dry_mass * mass_ratio
    fuel_mass = initial_mass - dry_mass

    return fuel_mass