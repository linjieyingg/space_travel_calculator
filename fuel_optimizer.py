import math

# Assume 'propulsion_system.py' and 'trajectory_planner.py' are peer modules
# in the same project structure.

# Import actual functions/classes from propulsion_system
# The mock implementation is being removed, so we directly attempt to import
# the calculate_required_fuel_mass function and alias it for consistent usage
# within this module.
from propulsion_system import calculate_required_fuel_mass as calculate_fuel_mass_from_delta_v

# Import TrajectoryPlanner class from trajectory_planner
try:
    from trajectory_planner import TrajectoryPlanner
except ImportError as e:
    # If TrajectoryPlanner cannot be imported, it's a critical error for this file's operation
    raise ImportError(f"Error: Could not import TrajectoryPlanner from trajectory_planner.py. "
                      f"Please ensure trajectory_planner.py is in the Python path. Details: {e}") from e


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

    This function integrates trajectory data (required delta-V) from the TrajectoryPlanner
    with propulsion system models to determine the necessary fuel mass using the
    Tsiolkovsky rocket equation.

    Args:
        start_body (str): The starting celestial body (e.g., 'Earth').
        end_body (str): The target celestial body (e.g., 'Mars').
        trajectory_type (str): The type of trajectory (e.g., 'Hohmann').
                               Currently, only 'Hohmann' is explicitly supported for planning.
        spacecraft_dry_mass (float): The mass of the spacecraft without fuel, in kg.
        engine_specific_impulse (float): The specific impulse (Isp) of the engine, in seconds.
                                        Isp is a measure of the efficiency of a rocket.
        **trajectory_args: Additional keyword arguments to pass to the trajectory planner.
                         (e.g., launch_date, arrival_date, intermediate_maneuvers,
                         parking_orbit_altitude). Note: TrajectoryPlanner.plan_hohmann_trajectory's
                         actual signature might not support all arbitrary kwargs, which will raise
                         a ValueError if passed.

    Returns:
        float: The calculated optimal fuel mass in kilograms.

    Raises:
        ValueError: If `spacecraft_dry_mass` or `engine_specific_impulse` are non-positive,
                    if the trajectory type is unsupported, or if the trajectory planner
                    or propulsion system return invalid values.
        RuntimeError: If an unexpected error occurs during the calculation,
                      e.g., a dependency function fails.
        ImportError: If 'propulsion_system.py' or 'trajectory_planner.py' cannot be found/imported.
    """
    if not isinstance(start_body, str) or not start_body:
        raise ValueError("Start body must be a non-empty string.")
    if not isinstance(end_body, str) or not end_body:
        raise ValueError("End body must be a non-empty string.")
    if not isinstance(trajectory_type, str) or not trajectory_type:
        raise ValueError("Trajectory type must be a non-empty string.")

    # Enforce 'Hohmann' as plan_hohmann_trajectory is the method used
    if trajectory_type.lower() != 'hohmann':
        raise ValueError(f"Unsupported trajectory type: '{trajectory_type}'. "
                         "Only 'Hohmann' is supported for planning at this time via TrajectoryPlanner.")

    if not isinstance(spacecraft_dry_mass, (int, float)) or spacecraft_dry_mass <= 0:
        raise ValueError("Spacecraft dry mass must be a positive numeric value.")
    if not isinstance(engine_specific_impulse, (int, float)) or engine_specific_impulse <= 0:
        raise ValueError("Engine specific impulse must be a positive numeric value.")

    try:
        # Instantiate the TrajectoryPlanner
        planner = TrajectoryPlanner()

        # 1. Plan the trajectory to get the required Delta-V
        # The prompt explicitly states to pass **trajectory_args to plan_hohmann_trajectory.
        # This might cause a TypeError if TrajectoryPlanner's signature doesn't match the kwargs.
        try:
            trajectory_result = planner.plan_hohmann_trajectory(
                departure_body_name=start_body,
                arrival_body_name=end_body,
                **trajectory_args
            )
        except TypeError as te:
            # Catch TypeError if plan_hohmann_trajectory's signature doesn't accept all **trajectory_args
            raise ValueError(
                f"Trajectory planner method signature mismatch: "
                f"plan_hohmann_trajectory may not accept all provided trajectory arguments. "
                f"Details: {te}"
            ) from te
        except Exception as e:
            # Catch other potential errors during trajectory planning
            raise RuntimeError(f"An unexpected error occurred during trajectory planning: {e}") from e

        if not isinstance(trajectory_result, dict) or not trajectory_result:
            raise ValueError(
                f"Trajectory planner returned an invalid or empty result: {trajectory_result}. "
                "Expected a dictionary with trajectory details."
            )

        if not trajectory_result.get('success', False):
            error_message = trajectory_result.get('error_message', 'Unknown error during trajectory planning.')
            raise ValueError(f"Trajectory planning failed: {error_message}")

        required_delta_v = trajectory_result.get('total_delta_v')

        if not isinstance(required_delta_v, (int, float)) or required_delta_v < 0:
            raise ValueError(
                f"Trajectory planner returned an invalid or negative 'total_delta_v': {required_delta_v}. "
                "The 'total_delta_v' must be a non-negative numeric value."
            )

        # 2. Calculate the fuel mass based on the required Delta-V and engine parameters
        # Calling the aliased function which now points to the real implementation from propulsion_system.
        fuel_mass = calculate_fuel_mass_from_delta_v(
            delta_v=required_delta_v,
            dry_mass=spacecraft_dry_mass,  # Use 'dry_mass' for clarity as per propulsion_system signature
            specific_impulse=engine_specific_impulse
        )

        if not isinstance(fuel_mass, (int, float)) or fuel_mass < 0:
            raise ValueError(
                f"Propulsion system returned an invalid or negative fuel mass: {fuel_mass}. "
                "Fuel mass must be a non-negative numeric value."
            )

        return fuel_mass

    except ValueError as ve:
        # Re-raise ValueErrors with added context
        raise ValueError(f"Fuel optimization input or calculation error: {ve}") from ve
    except Exception as e:
        # Catch any other unexpected errors from dependency calls or internal logic
        raise RuntimeError(f"An unexpected error occurred during fuel optimization: {e}") from e