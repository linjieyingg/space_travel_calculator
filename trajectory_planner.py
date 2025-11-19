import math
from datetime import timedelta
import orbital_mechanics
import celestial_data
from constants import G_GRAVITATIONAL


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
        for interplanetary travel) from the celestial_data module.
        """
        sun_data = celestial_data.get_celestial_body_data('Sun')
        if not sun_data:
            raise ValueError("Data for 'Sun' not found in celestial_data module. Cannot initialize TrajectoryPlanner.")
        
        self.mu_sun = sun_data.get('gravitational_parameter_mu')
        if self.mu_sun is None:
            raise AttributeError("'gravitational_parameter_mu' not found for 'Sun' in celestial_data. Please check celestial_data.py.")

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

        departure_body_data = celestial_data.get_celestial_body_data(departure_body_name)
        arrival_body_data = celestial_data.get_celestial_body_data(arrival_body_name)

        if not departure_body_data:
            return {
                "success": False,
                "error": f"Could not find data for departure body: '{departure_body_name}' in celestial_data."
            }
        if not arrival_body_data:
            return {
                "success": False,
                "error": f"Could not find data for arrival body: '{arrival_body_name}' in celestial_data."
            }

        # Use 'average_orbital_radius_m' from celestial_data
        r1 = departure_body_data.get("average_orbital_radius_m")
        r2 = arrival_body_data.get("average_orbital_radius_m")

        if r1 is None or not isinstance(r1, (int, float)) or r1 <= 0:
            return {
                "success": False,
                "error": f"Invalid or missing 'average_orbital_radius_m' for {departure_body_data['name']} in celestial_data."
            }
        if r2 is None or not isinstance(r2, (int, float)) or r2 <= 0:
            return {
                "success": False,
                "error": f"Invalid or missing 'average_orbital_radius_m' for {arrival_body_data['name']} in celestial_data."
            }

        try:
            # Use specific orbital_mechanics functions for delta-v and time of flight
            # calculate_hohmann_transfer_delta_v returns (delta_v1, delta_v2, total_delta_v)
            _, _, total_delta_v = orbital_mechanics.calculate_hohmann_transfer_delta_v(
                self.mu_sun, r1, r2
            )
            travel_time_s = orbital_mechanics.calculate_hohmann_transfer_time_of_flight(
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
    try:
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
    except (ValueError, AttributeError) as e:
        print(f"Initialization Error: {e}")