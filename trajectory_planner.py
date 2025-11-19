import math
from datetime import timedelta, datetime
import orbital_mechanics
import celestial_data
import ephemeris

class TrajectoryPlanner:
    """
    A class to provide a higher-level interface for planning interplanetary trajectories.
    It orchestrates calls to orbital mechanics functions to determine paths between
    celestial bodies, returning details like total delta-v and total travel time.
    This version integrates dynamic orbital positions (3D position and velocity vectors)
    based on a departure date using the ephemeris module and supports multiple
    (and extensible) trajectory types, including Hohmann transfers and direct transfers
    via Lambert's problem.
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
        if not isinstance(self.mu_sun, (int, float)) or self.mu_sun <= 0:
            raise ValueError(
                f"Invalid or non-positive 'gravitational_parameter_mu' ({self.mu_sun}) found for 'Sun' in celestial_data. "
                "Cannot initialize TrajectoryPlanner."
            )

    def _plan_hohmann_trajectory(self, departure_body_name: str, arrival_body_name: str, departure_date: datetime) -> dict:
        """
        Plans a Hohmann transfer trajectory between two celestial bodies orbiting the Sun.
        This is typically a minimum-energy transfer for coplanar, circular orbits.
        It considers the dynamic 3D orbital positions of the bodies on the specified departure date
        obtained from the ephemeris module.

        Args:
            departure_body_name (str): The name of the celestial body to depart from (e.g., "Earth").
            arrival_body_name (str): The name of the celestial body to arrive at (e.g., "Mars").
            departure_date (datetime.datetime): The date of departure, used to determine the
                                                dynamic orbital positions of the bodies.

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
                  - 'departure_distance_from_sun_m': Dynamic distance (magnitude of position vector) of departure body from Sun.
                  - 'arrival_distance_from_sun_m': Dynamic distance (magnitude of position vector) of arrival body from Sun.
        """
        if not isinstance(departure_body_name, str) or not departure_body_name.strip():
            return {"success": False, "error": "Departure body name must be a non-empty string."}
        if not isinstance(arrival_body_name, str) or not arrival_body_name.strip():
            return {"success": False, "error": "Arrival body name must be a non-empty string."}
        if not isinstance(departure_date, datetime):
            return {"success": False, "error": "Departure date must be a datetime object."}
        
        # Normalize names for comparison (case-insensitive and strip whitespace)
        if departure_body_name.strip().lower() == arrival_body_name.strip().lower():
            return {"success": False, "error": "Departure and arrival bodies cannot be the same."}

        # Use ephemeris to get dynamic orbital states (position and velocity vectors)
        try:
            # ephemeris.get_heliocentric_state is now expected to return (position_vector_m, velocity_vector_mps)
            # For Hohmann, we only need the position vector to derive radial distances.
            dep_position_vector_m, _ = ephemeris.get_heliocentric_state(departure_body_name, departure_date)
            arr_position_vector_m, _ = ephemeris.get_heliocentric_state(arrival_body_name, departure_date)
        except Exception as e:
            # Catch errors from ephemeris module (e.g., unsupported body, internal error)
            return {"success": False, "error": f"Error retrieving orbital states via ephemeris for date '{departure_date.strftime('%Y-%m-%d')}': {e}"}

        # Validate the retrieved position vectors
        if not dep_position_vector_m or not isinstance(dep_position_vector_m, tuple) or len(dep_position_vector_m) != 3:
            return {
                "success": False,
                "error": f"Invalid or missing 3D position vector for departure body: '{departure_body_name}' on {departure_date.strftime('%Y-%m-%d')} from ephemeris. Ensure body name is correct and ephemeris provides valid 3D vectors."
            }
        if not arr_position_vector_m or not isinstance(arr_position_vector_m, tuple) or len(arr_position_vector_m) != 3:
            return {
                "success": False,
                "error": f"Invalid or missing 3D position vector for arrival body: '{arrival_body_name}' on {departure_date.strftime('%Y-%m-%d')} from ephemeris. Ensure body name is correct and ephemeris provides valid 3D vectors."
            }

        try:
            # Calculate magnitudes from the 3D position vectors for Hohmann transfer radii (r1, r2)
            r1 = orbital_mechanics.vector_magnitude(dep_position_vector_m)
            r2 = orbital_mechanics.vector_magnitude(arr_position_vector_m)
        except ValueError as ve:
            return {
                "success": False,
                "error": f"Error calculating magnitude of position vector for Hohmann transfer: {ve}. Ensure ephemeris returns valid 3D vectors."
            }

        if r1 <= 0 or r2 <= 0:
             return {
                "success": False,
                "error": f"Calculated radial distance from Sun is zero or negative for a celestial body. Ensure ephemeris data is valid and provides non-zero position vectors."
            }
        
        # Ensure that the central body's gravitational parameter is valid before proceeding
        if self.mu_sun is None or not isinstance(self.mu_sun, (int, float)) or self.mu_sun <= 0:
             return {
                "success": False,
                "error": "Central body's (Sun's) gravitational parameter (mu_sun) is invalid. Cannot perform orbital calculations."
            }

        try:
            # orbital_mechanics.calculate_hohmann_transfer_delta_v returns (delta_v1, delta_v2, total_delta_v)
            _, _, total_delta_v = orbital_mechanics.calculate_hohmann_transfer_delta_v(
                self.mu_sun, r1, r2
            )
            travel_time_s = orbital_mechanics.calculate_hohmann_transfer_time_of_flight(
                self.mu_sun, r1, r2
            )
        except ValueError as ve:
            # Specific error from orbital_mechanics module functions
            return {
                "success": False,
                "error": f"Orbital mechanics calculation error for Hohmann transfer: {ve}"
            }
        except Exception as e:
            # Catch any unexpected errors from orbital_mechanics module or other calculations
            return {
                "success": False,
                "error": f"An unexpected error occurred during Hohmann transfer calculation: {e}"
            }

        # Convert travel time to more readable formats
        travel_time_days = travel_time_s / (24 * 3600)
        
        td = timedelta(seconds=int(travel_time_s)) 
        days = td.days
        hours, remainder = divmod(td.seconds, 3600)
        minutes, seconds = divmod(remainder, 60)
        
        travel_time_h_m_s = f"{days} days, {hours:02d} hours, {minutes:02d} minutes, {seconds:02d} seconds"

        return {
            "departure_body": departure_body_name,
            "arrival_body": arrival_body_name,
            "transfer_type": "Hohmann Transfer",
            "total_delta_v_mps": total_delta_v,
            "travel_time_days": travel_time_days,
            "travel_time_h_m_s": travel_time_h_m_s,
            "success": True,
            "error": None,
            "departure_distance_from_sun_m": r1,
            "arrival_distance_from_sun_m": r2,
        }

    def _plan_direct_transfer(self, departure_body_name: str, arrival_body_name: str, departure_date: datetime) -> dict:
        """
        Plans a direct transfer trajectory using Lambert's problem, optimizing for minimum total Delta-V
        across a range of flight times. It utilizes dynamic 3D orbital position and velocity vectors
        from the ephemeris module.

        Args:
            departure_body_name (str): The name of the celestial body to depart from (e.g., "Earth").
            arrival_body_name (str): The name of the celestial body to arrive at (e.g., "Mars").
            departure_date (datetime.datetime): The date of departure.

        Returns:
            dict: A dictionary containing trajectory details, including:
                  - 'departure_body': Name of the departure body.
                  - 'arrival_body': Name of the arrival body.
                  - 'transfer_type': "Direct Transfer".
                  - 'total_delta_v_mps': Total change in velocity required in meters per second.
                  - 'travel_time_days': Total travel time in days.
                  - 'travel_time_h_m_s': Formatted string of travel time (days, hours, minutes, seconds).
                  - 'success': True if planning was successful, False otherwise.
                  - 'error': Error message if planning failed.
                  - 'departure_distance_from_sun_m': Dynamic distance (magnitude of position vector) of departure body from Sun.
                  - 'arrival_distance_from_sun_m': Dynamic distance (magnitude of position vector) of arrival body from Sun.
        """
        if not isinstance(departure_body_name, str) or not departure_body_name.strip():
            return {"success": False, "error": "Departure body name must be a non-empty string."}
        if not isinstance(arrival_body_name, str) or not arrival_body_name.strip():
            return {"success": False, "error": "Arrival body name must be a non-empty string."}
        if not isinstance(departure_date, datetime):
            return {"success": False, "error": "Departure date must be a datetime object."}
        
        if departure_body_name.strip().lower() == arrival_body_name.strip().lower():
            return {"success": False, "error": "Departure and arrival bodies cannot be the same."}

        if self.mu_sun is None or not isinstance(self.mu_sun, (int, float)) or self.mu_sun <= 0:
             return {
                "success": False,
                "error": "Central body's (Sun's) gravitational parameter (mu_sun) is invalid. Cannot perform orbital calculations."
            }

        try:
            # 1. Get initial state (3D position and velocity vectors) for departure body at departure_date
            # The required changes state: "new (position_vector_m, velocity_vector_mps) output from ephemeris.get_heliocentric_state"
            dep_position_vector_m, dep_velocity_vector_mps = ephemeris.get_heliocentric_state(departure_body_name, departure_date)
            
            # Validate the retrieved position vector
            if not dep_position_vector_m or not isinstance(dep_position_vector_m, tuple) or len(dep_position_vector_m) != 3:
                return {"success": False, "error": f"Invalid or missing 3D position vector for departure body '{departure_body_name}' on {departure_date.strftime('%Y-%m-%d')} from ephemeris."}
            
            r1_vec = dep_position_vector_m
            r1_mag = orbital_mechanics.vector_magnitude(r1_vec)
            if r1_mag <= 0:
                return {"success": False, "error": f"Calculated radial distance from Sun is zero or negative for departure body '{departure_body_name}'. Ensure ephemeris data is valid and provides non-zero position vectors."}

        except ValueError as ve: # Catching specific validation errors from orbital_mechanics.vector_magnitude etc.
            return {"success": False, "error": f"Data processing error for departure body in Direct Transfer: {ve}"}
        except Exception as e: # Catching potential errors from ephemeris.get_heliocentric_state
            return {"success": False, "error": f"Error retrieving departure body data via ephemeris for Direct Transfer: {e}"}

        # 2. Iterate through possible times of flight (TOF) to find minimum Delta-V
        # Define TOF range and step (in days)
        min_tof_days = 50 
        max_tof_days = 700 
        tof_step_days = 5 

        min_total_delta_v = float('inf')
        best_trajectory_data = None
        
        # Iteration through possible TOFs
        for current_tof_days in range(int(min_tof_days), int(max_tof_days) + int(tof_step_days), int(tof_step_days)):
            tof_s = current_tof_days * 24 * 3600
            current_arrival_date = departure_date + timedelta(seconds=tof_s)

            try:
                # Get final state (3D position and velocity vectors) for arrival body at current_arrival_date
                arr_position_vector_m, arr_velocity_vector_mps = ephemeris.get_heliocentric_state(arrival_body_name, current_arrival_date)
                
                # Validate the retrieved position vector
                if not arr_position_vector_m or not isinstance(arr_position_vector_m, tuple) or len(arr_position_vector_m) != 3:
                    # Skip this TOF if arrival data is invalid (e.g., position vector is None or malformed)
                    continue 

                r2_vec = arr_position_vector_m
                r2_mag = orbital_mechanics.vector_magnitude(r2_vec)
                if r2_mag <= 0:
                    continue # Skip if radial distance is invalid
                
                # 3. Call Lambert Solver with actual 3D position vectors
                # orbital_mechanics.calculate_lambert_transfer returns:
                # (delta_v1_mag, delta_v2_mag, total_delta_v, calculated_time_of_flight)
                # as per the orbital_mechanics.py summary.
                dv1_mag_burn, dv2_mag_burn, calculated_total_delta_v, _ = orbital_mechanics.calculate_lambert_transfer(
                    self.mu_sun, r1_vec, r2_vec, tof_s
                )

                # 4. Check for valid delta-V results and update total Delta-V
                if not isinstance(calculated_total_delta_v, (int, float)) or calculated_total_delta_v < 0:
                    # Negative delta-V or non-numeric indicates an invalid solution from the simplified solver
                    # This check is crucial for the iterative Lambert solver, as some TOFs might not yield valid solutions.
                    raise ValueError("Lambert solver returned an invalid total Delta-V (negative or non-numeric).")
                
                # 5. Select the TOF that yields the minimum total delta-V
                if calculated_total_delta_v < min_total_delta_v:
                    min_total_delta_v = calculated_total_delta_v
                    
                    td = timedelta(seconds=int(tof_s)) 
                    days = td.days
                    hours, remainder = divmod(td.seconds, 3600)
                    minutes, seconds = divmod(remainder, 60)
                    travel_time_h_m_s = f"{days} days, {hours:02d} hours, {minutes:02d} minutes, {seconds:02d} seconds"

                    best_trajectory_data = {
                        "departure_body": departure_body_name,
                        "arrival_body": arrival_body_name,
                        "transfer_type": "Direct Transfer",
                        "total_delta_v_mps": min_total_delta_v,
                        "travel_time_days": current_tof_days,
                        "travel_time_h_m_s": travel_time_h_m_s,
                        "success": True,
                        "error": None,
                        "departure_distance_from_sun_m": r1_mag,
                        "arrival_distance_from_sun_m": r2_mag, # Store the arrival distance for *this* optimal TOF
                        # Optionally, add individual burn delta-Vs if needed by other modules
                        # "delta_v_burn1_mps": dv1_mag_burn, 
                        # "delta_v_burn2_mps": dv2_mag_burn,
                    }
            except ValueError as ve:
                # Specific errors from orbital_mechanics.calculate_lambert_transfer (e.g., no real solution found)
                # or from our check for negative delta-V
                # print(f"DEBUG: Lambert solver failed for TOF {current_tof_days} days: {ve}") # For debugging purposes
                continue # Try next TOF
            except Exception as e:
                # General errors during iteration (e.g., ephemeris failure for a specific date or vector magnitude issues)
                # print(f"DEBUG: General error for TOF {current_tof_days} days: {e}") # For debugging purposes
                continue # Try next TOF

        if best_trajectory_data:
            return best_trajectory_data
        else:
            return {"success": False, "error": "No suitable Direct Transfer trajectory found within the specified time of flight range. Consider adjusting the date or TOF range, or checking ephemeris data availability and correctness."}

    def plan_trajectory(self, departure_body_name: str, arrival_body_name: str, trajectory_type: str = 'Hohmann', **kwargs) -> dict:
        """
        Plans a space trajectory between two celestial bodies using various algorithms.
        This method acts as a dispatcher for different trajectory planning algorithms.

        Args:
            departure_body_name (str): The name of the celestial body to depart from (e.g., "Earth").
            arrival_body_name (str): The name of the celestial body to arrive at (e.g., "Mars").
            trajectory_type (str, optional): The type of trajectory to plan. Defaults to 'Hohmann'.
                                             Supported types: 'Hohmann', 'Direct'.
            **kwargs: Additional arguments specific to the chosen trajectory type.
                      For 'Hohmann' and 'Direct' transfers, 'departure_date' (str, 'YYYY-MM-DD') is required.

        Returns:
            dict: A dictionary containing trajectory details, including:
                  - 'departure_body': Name of the departure body.
                  - 'arrival_body': Name of the arrival body.
                  - 'transfer_type': The type of transfer.
                  - 'total_delta_v_mps': Total change in velocity required in meters per second.
                  - 'travel_time_days': Total travel time in days.
                  - 'travel_time_h_m_s': Formatted string of travel time (days, hours, minutes, seconds).
                  - 'success': True if planning was successful, False otherwise.
                  - 'error': Error message if planning failed.
                  - 'departure_distance_from_sun_m' (optional): Dynamic distance (magnitude of position vector) of departure body from Sun.
                  - 'arrival_distance_from_sun_m' (optional): Dynamic distance (magnitude of position vector) of arrival body from Sun.
        """
        if not isinstance(departure_body_name, str) or not departure_body_name.strip():
            return {"success": False, "error": "Departure body name must be a non-empty string."}
        if not isinstance(arrival_body_name, str) or not arrival_body_name.strip():
            return {"success": False, "error": "Arrival body name must be a non-empty string."}
        if not isinstance(trajectory_type, str) or not trajectory_type.strip():
            return {"success": False, "error": "Trajectory type must be a non-empty string."}
        
        normalized_type = trajectory_type.strip().lower()

        # Extract and validate departure_date from kwargs, convert to datetime object
        departure_date_str = kwargs.get('departure_date')
        if not departure_date_str:
            return {"success": False, "error": f"'{trajectory_type}' trajectory planning requires a 'departure_date' in kwargs."}
        if not isinstance(departure_date_str, str) or not departure_date_str.strip():
            return {"success": False, "error": "Departure date in kwargs must be a non-empty string in 'YYYY-MM-DD' format."}
        
        try:
            # Assuming 'YYYY-MM-DD' format for parsing from user input
            departure_date_obj = datetime.strptime(departure_date_str, '%Y-%m-%d')
        except ValueError:
            return {"success": False, "error": f"Invalid 'departure_date' format: '{departure_date_str}'. Expected 'YYYY-MM-DD'."}

        if normalized_type == 'hohmann':
            return self._plan_hohmann_trajectory(departure_body_name, arrival_body_name, departure_date_obj)
        elif normalized_type == 'direct':
            return self._plan_direct_transfer(departure_body_name, arrival_body_name, departure_date_obj)
        else:
            return {"success": False, "error": f"Unsupported trajectory type: '{trajectory_type}'. Supported types: 'Hohmann', 'Direct'."}

if __name__ == '__main__':
    # Example usage for testing
    try:
        planner = TrajectoryPlanner()
        test_date_str = "2024-07-20" # Define a test date string
        
        print("--- Planning Earth to Mars Hohmann Transfer ---")
        earth_to_mars = planner.plan_trajectory("Earth", "Mars", trajectory_type='Hohmann', departure_date=test_date_str)
        if earth_to_mars["success"]:
            print(f"Departure: {earth_to_mars['departure_body']}")
            print(f"Arrival: {earth_to_mars['arrival_body']}")
            print(f"Transfer Type: {earth_to_mars['transfer_type']}")
            print(f"Departure Distance from Sun: {earth_to_mars['departure_distance_from_sun_m']:.2e} m")
            print(f"Arrival Distance from Sun: {earth_to_mars['arrival_distance_from_sun_m']:.2e} m")
            print(f"Total Delta-V: {earth_to_mars['total_delta_v_mps']:.2f} m/s")
            print(f"Travel Time (days): {earth_to_mars['travel_time_days']:.2f}")
            print(f"Travel Time (H:M:S): {earth_to_mars['travel_time_h_m_s']}")
        else:
            print(f"Error: {earth_to_mars['error']}")

        print("\n--- Planning Mars to Jupiter Hohmann Transfer (using default type) ---")
        mars_to_jupiter = planner.plan_trajectory("Mars", "Jupiter", departure_date=test_date_str)
        if mars_to_jupiter["success"]:
            print(f"Departure: {mars_to_jupiter['departure_body']}")
            print(f"Arrival: {mars_to_jupiter['arrival_body']}")
            print(f"Transfer Type: {mars_to_jupiter['transfer_type']}")
            print(f"Departure Distance from Sun: {mars_to_jupiter['departure_distance_from_sun_m']:.2e} m")
            print(f"Arrival Distance from Sun: {mars_to_jupiter['arrival_distance_from_sun_m']:.2e} m")
            print(f"Total Delta-V: {mars_to_jupiter['total_delta_v_mps']:.2f} m/s")
            print(f"Travel Time (days): {mars_to_jupiter['travel_time_days']:.2f}")
            print(f"Travel Time (H:M:S): {mars_to_jupiter['travel_time_h_m_s']}")
        else:
            print(f"Error: {mars_to_jupiter['error']}")

        print("\n--- Planning Earth to Venus Direct Transfer (using implemented Lambert solver logic) ---")
        # Note: Actual Lambert solver in orbital_mechanics.py must be properly implemented for accurate results.
        # Placeholder implementation might yield illustrative but not physically accurate results.
        earth_to_venus_direct = planner.plan_trajectory("Earth", "Venus", trajectory_type='Direct', departure_date=test_date_str)
        if earth_to_venus_direct["success"]:
            print(f"Departure: {earth_to_venus_direct['departure_body']}")
            print(f"Arrival: {earth_to_venus_direct['arrival_body']}")
            print(f"Transfer Type: {earth_to_venus_direct['transfer_type']}")
            print(f"Departure Distance from Sun: {earth_to_venus_direct['departure_distance_from_sun_m']:.2e} m")
            print(f"Arrival Distance from Sun: {earth_to_venus_direct['arrival_distance_from_sun_m']:.2e} m")
            print(f"Total Delta-V: {earth_to_venus_direct['total_delta_v_mps']:.2f} m/s")
            print(f"Travel Time (days): {earth_to_venus_direct['travel_time_days']:.2f}")
            print(f"Travel Time (H:M:S): {earth_to_venus_direct['travel_time_h_m_s']}")
        else:
            print(f"Result: {earth_to_venus_direct['success']}, Error: {earth_to_venus_direct['error']}")

        print("\n--- Planning Invalid Trajectory (same bodies) ---")
        invalid_plan = planner.plan_trajectory("Earth", "earth", departure_date=test_date_str) # Default type
        print(f"Result: {invalid_plan['success']}, Error: {invalid_plan['error']}")

        print("\n--- Planning Invalid Trajectory (unknown body for ephemeris) ---")
        # Pluto is in celestial_data, but ephemeris.get_heliocentric_state might raise an error or return None for it.
        unknown_plan = planner.plan_trajectory("Earth", "Pluto", departure_date=test_date_str)
        print(f"Result: {unknown_plan['success']}, Error: {unknown_plan['error']}")

        print("\n--- Planning Invalid Trajectory (empty body name) ---")
        empty_name_plan = planner.plan_trajectory("  ", "Mars", departure_date=test_date_str)
        print(f"Result: {empty_name_plan['success']}, Error: {empty_name_plan['error']}")
        
        print("\n--- Planning Invalid Trajectory (non-string body name) ---")
        # This will now fail due to explicit type validation within plan_trajectory
        non_string_plan = planner.plan_trajectory(123, "Mars", departure_date=test_date_str)
        print(f"Result: {non_string_plan['success']}, Error: {non_string_plan['error']}")

        print("\n--- Planning Invalid Trajectory (unsupported type) ---")
        unsupported_plan = planner.plan_trajectory("Earth", "Mars", trajectory_type='Bi-elliptic', departure_date=test_date_str)
        print(f"Result: {unsupported_plan['success']}, Error: {unsupported_plan['error']}")

        print("\n--- Planning Invalid Trajectory (missing departure_date) ---")
        missing_date_plan = planner.plan_trajectory("Earth", "Mars", trajectory_type='Hohmann')
        print(f"Result: {missing_date_plan['success']}, Error: {missing_date_plan['error']}")

        print("\n--- Planning Invalid Trajectory (invalid departure_date format) ---")
        invalid_date_format_plan = planner.plan_trajectory("Earth", "Mars", trajectory_type='Hohmann', departure_date="2024/07/20")
        print(f"Result: {invalid_date_format_plan['success']}, Error: {invalid_date_format_plan['error']}")


    except (ValueError, AttributeError) as e:
        print(f"Initialization Error: {e}")