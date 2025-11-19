import math
from datetime import timedelta

# --- MOCK IMPORTS (for testing/standalone execution) ---
# In a full project, these would be separate, fully-fleshed out modules.
# We include try-except blocks to allow the file to run even if these are not
# yet fully implemented or correctly structured in the user's project,
# by providing minimal mock implementations.

try:
    import orbital_mechanics
except ImportError:
    class MockOrbitalMechanics:
        """
        Mocks orbital_mechanics module for trajectory calculations.
        Provides a simplified Hohmann transfer calculation.
        Assumes circular, coplanar orbits around a central body.
        """
        def calculate_hohmann_transfer(self, mu_central: float, r1: float, r2: float) -> tuple[float, float]:
            """
            Mocks the calculation of delta-v and travel time for a Hohmann transfer.
            Args:
                mu_central (float): Standard gravitational parameter of the central body (m^3/s^2).
                r1 (float): Radius of the initial circular orbit (m).
                r2 (float): Radius of the final circular orbit (m).
            Returns:
                tuple[float, float]: (total_delta_v_mps, travel_time_s)
            """
            if not all(isinstance(x, (int, float)) and x > 0 for x in [mu_central, r1, r2]):
                raise ValueError("Inputs for Hohmann transfer must be positive numbers.")
            
            # Orbital velocities of the circular orbits
            v1 = math.sqrt(mu_central / r1)
            v2 = math.sqrt(mu_central / r2)

            # Semi-major axis of the transfer ellipse
            a_transfer = (r1 + r2) / 2.0

            # Velocities at periapsis and apoapsis of the transfer ellipse
            vp_transfer = math.sqrt(mu_central * (2/r1 - 1/a_transfer))
            va_transfer = math.sqrt(mu_central * (2/r2 - 1/a_transfer))

            # Delta-v for departure and arrival
            delta_v1 = abs(vp_transfer - v1)
            delta_v2 = abs(v2 - va_transfer)

            total_delta_v = delta_v1 + delta_v2

            # Travel time (half the period of the transfer ellipse)
            period_transfer = 2 * math.pi * math.sqrt(a_transfer**3 / mu_central)
            travel_time_s = period_transfer / 2.0

            return total_delta_v, travel_time_s
    orbital_mechanics = MockOrbitalMechanics()


try:
    import celestial_bodies
except ImportError:
    class MockCelestialBodies:
        """
        Mocks a celestial_bodies module to provide data for various celestial bodies.
        """
        _bodies = {
            "mercury": {"name": "Mercury", "mass_kg": 3.3011e23, "radius_m": 2.4397e6, "orbital_radius_m": 5.7909e10, "std_grav_param_m3_s2": 2.2032e13},
            "venus": {"name": "Venus", "mass_kg": 4.8675e24, "radius_m": 6.0518e6, "orbital_radius_m": 1.0821e11, "std_grav_param_m3_s2": 3.24859e14},
            "earth": {"name": "Earth", "mass_kg": 5.972e24, "radius_m": 6.371e6, "orbital_radius_m": 1.496e11, "std_grav_param_m3_s2": 3.986004418e14},
            "mars": {"name": "Mars", "mass_kg": 6.39e23, "radius_m": 3.3895e6, "orbital_radius_m": 2.279e11, "std_grav_param_m3_s2": 4.282837e13},
            "jupiter": {"name": "Jupiter", "mass_kg": 1.898e27, "radius_m": 6.9911e7, "orbital_radius_m": 7.785e11, "std_grav_param_m3_s2": 1.26686534e17},
            "saturn": {"name": "Saturn", "mass_kg": 5.683e26, "radius_m": 5.8232e7, "orbital_radius_m": 1.434e12, "std_grav_param_m3_s2": 3.79312e16},
            "uranus": {"name": "Uranus", "mass_kg": 8.681e25, "radius_m": 2.5362e7, "orbital_radius_m": 2.871e12, "std_grav_param_m3_s2": 5.793939e15},
            "neptune": {"name": "Neptune", "mass_kg": 1.024e26, "radius_m": 2.4622e7, "orbital_radius_m": 4.495e12, "std_grav_param_m3_s2": 6.836534e15},
            "moon": {"name": "Moon", "mass_kg": 7.342e22, "radius_m": 1.7374e6, "orbital_radius_m": 3.844e8, "std_grav_param_m3_s2": 4.9048695e12} # Relative to Earth
        }

        def get_celestial_body_data(self, body_name: str) -> dict | None:
            """
            Retrieves data for a given celestial body.
            """
            return self._bodies.get(body_name.lower())

        def get_all_celestial_bodies(self) -> list[dict]:
            """
            Retrieves data for all known celestial bodies.
            """
            return list(self._bodies.values())
    celestial_bodies = MockCelestialBodies()


try:
    import constants
except ImportError:
    class MockConstants:
        """
        Mocks a constants module to provide universal physical constants.
        """
        GRAVITATIONAL_CONSTANT = 6.67430e-11  # N(m/kg)^2
        SUN_STANDARD_GRAVITATIONAL_PARAMETER = 1.32712440018e20 # m^3/s^2 (GM_sun)
    constants = MockConstants()
# --- END MOCK IMPORTS ---


class TrajectoryPlanner:
    """
    A class to provide a higher-level interface for planning interplanetary trajectories.
    It orchestrates calls to orbital mechanics functions to determine the optimal
    (e.g., minimum energy) path between two celestial bodies, returning details
    like total delta-v and total travel time.
    """

    def __init__(self):
        """
        Initializes the TrajectoryPlanner.
        Sets the standard gravitational parameter of the central body (assumed to be the Sun
        for interplanetary travel) from the constants module.
        """
        try:
            self.mu_sun = constants.SUN_STANDARD_GRAVITATIONAL_PARAMETER
        except AttributeError:
            raise AttributeError("SUN_STANDARD_GRAVITATIONAL_PARAMETER not found in constants module.")

    def plan_hohmann_trajectory(self, departure_body_name: str, arrival_body_name: str) -> dict:
        """
        Plans a Hohmann transfer trajectory between two celestial bodies orbiting the Sun.
        This is typically a minimum-energy transfer for coplanar, circular orbits.

        Args:
            departure_body_name (str): The name of the celestial body to depart from (e.g., "Earth").
            arrival_body_name (str): The name of the celestial body to arrive at (e.g., "Mars").

        Returns:
            dict: A dictionary containing trajectory details, including:
                  - 'departure_body': Name of the departure body.
                  - 'arrival_body': Name of the arrival body.
                  - 'transfer_type': "Hohmann Transfer".
                  - 'total_delta_v_mps': Total change in velocity required in meters per second.
                  - 'travel_time_days': Total travel time in days.
                  - 'travel_time_h_m_s': Formatted string of travel time (days, hours, minutes, seconds).
                  - 'success': True if planning was successful, False otherwise.
                  - 'error': Error message if planning failed.

        Raises:
            ValueError: If celestial body data cannot be retrieved or input is invalid.
        """
        if not isinstance(departure_body_name, str) or not departure_body_name.strip():
            return {"success": False, "error": "Departure body name must be a non-empty string."}
        if not isinstance(arrival_body_name, str) or not arrival_body_name.strip():
            return {"success": False, "error": "Arrival body name must be a non-empty string."}
        if departure_body_name.lower() == arrival_body_name.lower():
            return {"success": False, "error": "Departure and arrival bodies cannot be the same."}

        departure_body_data = celestial_bodies.get_celestial_body_data(departure_body_name)
        arrival_body_data = celestial_bodies.get_celestial_body_data(arrival_body_name)

        if not departure_body_data:
            return {
                "success": False,
                "error": f"Could not find data for departure body: '{departure_body_name}'"
            }
        if not arrival_body_data:
            return {
                "success": False,
                "error": f"Could not find data for arrival body: '{arrival_body_name}'"
            }

        r1 = departure_body_data.get("orbital_radius_m")
        r2 = arrival_body_data.get("orbital_radius_m")

        if r1 is None or not isinstance(r1, (int, float)) or r1 <= 0:
            return {
                "success": False,
                "error": f"Invalid or missing 'orbital_radius_m' for {departure_body_data['name']} in celestial_bodies data."
            }
        if r2 is None or not isinstance(r2, (int, float)) or r2 <= 0:
            return {
                "success": False,
                "error": f"Invalid or missing 'orbital_radius_m' for {arrival_body_data['name']} in celestial_bodies data."
            }

        try:
            total_delta_v, travel_time_s = orbital_mechanics.calculate_hohmann_transfer(
                self.mu_sun, r1, r2
            )
        except ValueError as ve:
            return {
                "success": False,
                "error": f"Orbital mechanics calculation error: {ve}"
            }
        except Exception as e:
            # Catch any unexpected errors from orbital_mechanics module
            return {
                "success": False,
                "error": f"An unexpected error occurred during Hohmann transfer calculation: {e}"
            }

        # Convert travel time to more readable formats
        travel_time_days = travel_time_s / (24 * 3600)
        
        td = timedelta(seconds=int(travel_time_s)) # Use int to avoid float precision issues with timedelta
        days = td.days
        hours, remainder = divmod(td.seconds, 3600)
        minutes, seconds = divmod(remainder, 60)
        
        travel_time_h_m_s = f"{days} days, {hours:02d} hours, {minutes:02d} minutes, {seconds:02d} seconds"


        return {
            "departure_body": departure_body_data["name"],
            "arrival_body": arrival_body_data["name"],
            "transfer_type": "Hohmann Transfer",
            "total_delta_v_mps": total_delta_v,
            "travel_time_days": travel_time_days,
            "travel_time_h_m_s": travel_time_h_m_s,
            "success": True,
            "error": None
        }

    # Future methods could be added here for other trajectory types, e.g.:
    # def plan_gravity_assist_trajectory(self, departure_body: str, arrival_body: str, flyby_body: str) -> dict:
    #     """
    #     Plans a trajectory utilizing a gravity assist.
    #     """
    #     pass


if __name__ == '__main__':
    # Example usage for testing
    planner = TrajectoryPlanner()

    print("--- Planning Earth to Mars Hohmann Transfer ---")
    earth_to_mars = planner.plan_hohmann_trajectory("Earth", "Mars")
    if earth_to_mars["success"]:
        print(f"Departure: {earth_to_mars['departure_body']}")
        print(f"Arrival: {earth_to_mars['arrival_body']}")
        print(f"Transfer Type: {earth_to_mars['transfer_type']}")
        print(f"Total Delta-V: {earth_to_mars['total_delta_v_mps']:.2f} m/s")
        print(f"Travel Time (days): {earth_to_mars['travel_time_days']:.2f}")
        print(f"Travel Time (H:M:S): {earth_to_mars['travel_time_h_m_s']}")
    else:
        print(f"Error: {earth_to_mars['error']}")

    print("\n--- Planning Mars to Jupiter Hohmann Transfer ---")
    mars_to_jupiter = planner.plan_hohmann_trajectory("Mars", "Jupiter")
    if mars_to_jupiter["success"]:
        print(f"Departure: {mars_to_jupiter['departure_body']}")
        print(f"Arrival: {mars_to_jupiter['arrival_body']}")
        print(f"Transfer Type: {mars_to_jupiter['transfer_type']}")
        print(f"Total Delta-V: {mars_to_jupiter['total_delta_v_mps']:.2f} m/s")
        print(f"Travel Time (days): {mars_to_jupiter['travel_time_days']:.2f}")
        print(f"Travel Time (H:M:S): {mars_to_jupiter['travel_time_h_m_s']}")
    else:
        print(f"Error: {mars_to_jupiter['error']}")

    print("\n--- Planning Invalid Trajectory (same bodies) ---")
    invalid_plan = planner.plan_hohmann_trajectory("Earth", "earth")
    print(f"Result: {invalid_plan['success']}, Error: {invalid_plan['error']}")

    print("\n--- Planning Invalid Trajectory (unknown body) ---")
    unknown_plan = planner.plan_hohmann_trajectory("Earth", "Pluto")
    print(f"Result: {unknown_plan['success']}, Error: {unknown_plan['error']}")

    print("\n--- Planning Invalid Trajectory (empty body name) ---")
    empty_name_plan = planner.plan_hohmann_trajectory("  ", "Mars")
    print(f"Result: {empty_name_plan['success']}, Error: {empty_name_plan['error']}")