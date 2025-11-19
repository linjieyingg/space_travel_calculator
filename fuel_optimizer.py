import math
import sys
import os

# Assume 'propulsion_system.py' and 'trajectory_planner.py' are peer modules
# in the same project structure.
# If they are in a different package or directory, relative imports might change.
# For now, assuming direct imports if they are at the root level like this file.

# Placeholder imports - these files are expected to be created and contain the described functions
# For this file to be executable, these modules must exist in the Python path.
try:
    from propulsion_system import calculate_fuel_mass_from_delta_v
except ImportError:
    # This block provides a mock implementation for development/testing if the file doesn't exist yet
    print("Warning: propulsion_system.py not found. Using mock function.", file=sys.stderr)
    G0_STANDARD = 9.80665  # Standard gravity in m/s^2

    def calculate_fuel_mass_from_delta_v(delta_v: float, spacecraft_dry_mass: float, specific_impulse: float) -> float:
        """
        MOCK: Calculates the fuel mass required for a given delta-V.
        In a real implementation, this would be in propulsion_system.py.
        """
        if delta_v < 0:
            raise ValueError("Delta-V cannot be negative.")
        if spacecraft_dry_mass <= 0:
            raise ValueError("Spacecraft dry mass must be positive.")
        if specific_impulse <= 0:
            raise ValueError("Specific impulse must be positive.")

        effective_exhaust_velocity = specific_impulse * G0_STANDARD
        if effective_exhaust_velocity == 0:
            raise ValueError("Effective exhaust velocity cannot be zero.")

        mass_ratio_exponent = delta_v / effective_exhaust_velocity
        fuel_mass = spacecraft_dry_mass * (math.exp(mass_ratio_exponent) - 1)
        return fuel_mass

try:
    from trajectory_planner import calculate_delta_v
except ImportError:
    # This block provides a mock implementation for development/testing if the file doesn't exist yet
    print("Warning: trajectory_planner.py not found. Using mock function.", file=sys.stderr)

    def calculate_delta_v(start_body: str, end_body: str, trajectory_type: str, **kwargs) -> float:
        """
        MOCK: Calculates the required delta-V for a given trajectory.
        In a real implementation, this would be in trajectory_planner.py.
        """
        start_body = start_body.lower()
        end_body = end_body.lower()
        trajectory_type = trajectory_type.lower()

        if start_body == 'earth' and end_body == 'mars':
            if trajectory_type == 'hohmann':
                return 5.7 * 1000  # Example m/s
            elif trajectory_type == 'direct':
                return 8.0 * 1000  # Example m/s
        elif start_body == 'earth' and end_body == 'moon':
            if trajectory_type == 'direct':
                return 3.1 * 1000  # Example m/s
        elif start_body == 'moon' and end_body == 'earth':
            if trajectory_type == 'direct':
                return 1.2 * 1000  # Example m/s
        elif start_body == 'earth' and end_body == 'jupiter':
            if trajectory_type == 'hohmann':
                return 15.0 * 1000 # Example m/s
        
        raise ValueError(f"MOCK ERROR: Unsupported trajectory: {start_body} to {end_body} via {trajectory_type}")


def optimize_fuel_for_trajectory(
    start_body: str,
    end_body: str,
    trajectory_type: str,
    spacecraft_dry_mass: float,
    engine_specific_impulse: float,
    **trajectory_args
) -> float:
    """
    Calculates the optimal fuel mass required for a planned trajectory.

    This function integrates trajectory data (required delta-V) with propulsion system
    models to determine the necessary fuel mass using the Tsiolkovsky rocket equation.

    Args:
        start_body (str): The starting celestial body (e.g., 'Earth').
        end_body (str): The target celestial body (e.g., 'Mars').
        trajectory_type (str): The type of trajectory (e.g., 'Hohmann', 'direct').
        spacecraft_dry_mass (float): The mass of the spacecraft without fuel, in kg.
        engine_specific_impulse (float): The specific impulse (Isp) of the engine, in seconds.
                                        Isp is a measure of the efficiency of a rocket.
        **trajectory_args: Additional keyword arguments to pass to the trajectory planner
                         (e.g., launch_date, arrival_date, intermediate_maneuvers,
                         parking_orbit_altitude).

    Returns:
        float: The calculated optimal fuel mass in kilograms.

    Raises:
        ValueError: If `spacecraft_dry_mass` or `engine_specific_impulse` are non-positive,
                    or if the trajectory planner or propulsion system return invalid values.
        RuntimeError: If an unexpected error occurs during the calculation,
                      e.g., a dependency function fails.
    """
    if not isinstance(start_body, str) or not start_body:
        raise ValueError("Start body must be a non-empty string.")
    if not isinstance(end_body, str) or not end_body:
        raise ValueError("End body must be a non-empty string.")
    if not isinstance(trajectory_type, str) or not trajectory_type:
        raise ValueError("Trajectory type must be a non-empty string.")
    if not isinstance(spacecraft_dry_mass, (int, float)) or spacecraft_dry_mass <= 0:
        raise ValueError("Spacecraft dry mass must be a positive numeric value.")
    if not isinstance(engine_specific_impulse, (int, float)) or engine_specific_impulse <= 0:
        raise ValueError("Engine specific impulse must be a positive numeric value.")

    try:
        # 1. Calculate the required Delta-V (change in velocity) for the trajectory
        # The trajectory_planner function will perform complex orbital mechanics calculations
        # and return the total delta-V needed for the mission.
        required_delta_v = calculate_delta_v(
            start_body=start_body,
            end_body=end_body,
            trajectory_type=trajectory_type,
            **trajectory_args
        )

        if not isinstance(required_delta_v, (int, float)) or required_delta_v < 0:
            raise ValueError(
                f"Trajectory planner returned an invalid or negative delta-V: {required_delta_v}. "
                "Delta-V must be a non-negative numeric value."
            )

        # 2. Calculate the fuel mass based on the required Delta-V and engine parameters
        # The propulsion_system function applies the Tsiolkovsky rocket equation.
        fuel_mass = calculate_fuel_mass_from_delta_v(
            delta_v=required_delta_v,
            spacecraft_dry_mass=spacecraft_dry_mass,
            specific_impulse=engine_specific_impulse
        )

        if not isinstance(fuel_mass, (int, float)) or fuel_mass < 0:
            raise ValueError(
                f"Propulsion system returned an invalid or negative fuel mass: {fuel_mass}. "
                "Fuel mass must be a non-negative numeric value."
            )

        return fuel_mass

    except ValueError as ve:
        # Re-raise ValueErrors with context
        raise ValueError(f"Fuel optimization input or calculation error: {ve}") from ve
    except Exception as e:
        # Catch any other unexpected errors from dependency calls or internal logic
        raise RuntimeError(f"An unexpected error occurred during fuel optimization: {e}") from e