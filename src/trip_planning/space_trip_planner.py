```python
import math
from typing import Optional, Dict, Any

# Relative imports from the project structure
from src.astrophysics.orbital_calculations import calculate_hohmann_delta_v
from celestial_data import get_celestial_body_data, AU_TO_METERS, CELESTIAL_BODIES_DATA

class SpaceTripPlanner:
    """
    Implements the core logic for planning Hohmann transfer space trips,
    displaying trip summaries, and providing an interactive user interface
    for trip planning.
    """

    SECONDS_PER_DAY = 86400.0
    DAYS_PER_YEAR = 365.25

    def __init__(self):
        """
        Initializes the SpaceTripPlanner by loading necessary celestial data.
        """
        # Direct access to CELESTIAL_BODIES_DATA for listing and specific parameters
        self._celestial_body_raw_data = CELESTIAL_BODIES_DATA 
        self._mu_sun = self._celestial_body_raw_data.get("SUN", {}).get("gravitational_parameter_mu")
        if self._mu_sun is None:
            raise ValueError("Sun's gravitational parameter (mu) not found in celestial data. Cannot initialize planner.")

    def _get_body_orbital_radius_meters(self, body_name: str) -> Optional[float]:
        """
        Retrieves the semi-major axis (approximated as orbital radius for Hohmann transfers)
        of a celestial body in meters.

        Args:
            body_name (str): The name of the celestial body.

        Returns:
            Optional[float]: The semi-major axis in meters, or None if data is not found.
        """
        body_data = get_celestial_body_data(body_name)
        if body_data and "semi_major_axis" in body_data:
            # semi_major_axis in celestial_data is in AU, convert to meters
            return body_data["semi_major_axis"] * AU_TO_METERS
        return None

    def _calculate_hohmann_transfer_time(self, r1: float, r2: float, mu_central: float) -> float:
        """
        Calculates the time of flight for a Hohmann transfer orbit.

        Args:
            r1 (float): Initial orbital radius (meters).
            r2 (float): Final orbital radius (meters).
            mu_central (float): Standard gravitational parameter of the central body (m^3/s^2).

        Returns:
            float: Transfer time in seconds.

        Raises:
            ValueError: If orbital radii or central body mu are non-positive.
        """
        if r1 <= 0 or r2 <= 0 or mu_central <= 0:
            raise ValueError("Orbital radii and central body mu must be positive for transfer time calculation.")

        # Semi-major axis of the transfer ellipse
        a_transfer = (r1 + r2) / 2
        
        # Period of the transfer ellipse (Kepler's Third Law: P = 2*pi*sqrt(a^3/mu))
        p_transfer = 2 * math.pi * math.sqrt(a_transfer**3 / mu_central)
        
        # Hohmann transfer time is half the period of the transfer ellipse
        return p_transfer / 2

    def plan_trip(self, source_body_name: str, target_body_name: str) -> Optional[Dict[str, Any]]:
        """
        Plans a space trip between two celestial bodies using a Hohmann transfer.

        Args:
            source_body_name (str): The name of the departure celestial body.
            target_body_name (str): The name of the destination celestial body.

        Returns:
            Optional[Dict[str, Any]]: A dictionary containing trip details
            (delta_v, transfer time) or None if planning fails.
        """
        if source_body_name.lower() == target_body_name.lower():
            print(f"Cannot plan a trip from {source_body_name} to itself. Please choose different bodies.")
            return None

        # Retrieve orbital radii in meters for transfer time calculation
        r_source = self._get_body_orbital_radius_meters(source_body_name)
        r_target = self._get_body_orbital_radius_meters(target_body_name)

        if r_source is None:
            print(f"Error: Could not retrieve orbital data for '{source_body_name}'.")
            return None
        if r_target is None:
            print(f"Error: Could not retrieve orbital data for '{target_body_name}'.")
            return None

        # Calculate Delta-V for the Hohmann transfer
        # Assumes transfers are heliocentric ('Sun' as central body)
        hohmann_dv_results = calculate_hohmann_delta_v(source_body_name, target_body_name, "Sun")

        if hohmann_dv_results is None:
            print("Error: Failed to calculate Hohmann Delta-V. Check body names or data integrity.")
            return None

        total_delta_v = hohmann_dv_results.get("total_delta_v")
        if total_delta_v is None:
            print("Error: Total Delta-V was not found in Hohmann transfer results.")
            return None

        # Calculate transfer time
        try:
            transfer_time_seconds = self._calculate_hohmann_transfer_time(
                r_source, r_target, self._mu_sun
            )
        except ValueError as e:
            print(f"Error calculating transfer time: {e}")
            return None
        
        transfer_time_days = transfer_time_seconds / self.SECONDS_PER_DAY
        transfer_time_years = transfer_time_days / self.DAYS_PER_YEAR

        return {
            "source": source_body_name,
            "target": target_body_name,
            "delta_v1": hohmann_dv_results.get("delta_v1"), # Delta-V at departure
            "delta_v2": hohmann_dv_results.get("delta_v2"), # Delta-V at arrival
            "total_delta_v": total_delta_v,
            "transfer_time_seconds": transfer_time_seconds,
            "transfer_time_days": transfer_time_days,
            "transfer_time_years": transfer_time_years,
        }

    def display_trip_summary(self, trip_data: Dict[str, Any]):
        """
        Displays a formatted summary of the planned Hohmann transfer trip.

        Args:
            trip_data (Dict[str, Any]): The dictionary containing trip details.
        """
        print("\n" + "="*50)
        print(f"  HOHMANN TRANSFER SUMMARY: {trip_data['source']} to {trip_data['target']}")
        print("="*50)
        print(f"  Delta-V at Departure (Δv₁): {trip_data['delta_v1']:.2f} m/s")
        print(f"  Delta-V at Arrival (Δv₂):   {trip_data['delta_v2']:.2f} m/s")
        print(f"  Total Delta-V Required:    {trip_data['total_delta_v']:.2f} m/s")
        print("-" * 50)
        print(f"  Hohmann Transfer Time of Flight:")
        print(f"    {trip_data['transfer_time_seconds'] / 3600:.2f} hours")
        print(f"    {trip_data['transfer_time_days']:.2f} days")
        print(f"    {trip_data['transfer_time_years']:.2f} years")
        print("="*50 + "\n")

    def run_interactive_planner(self):
        """
        Runs an interactive command-line interface for planning space trips.
        """
        print("\n" + "#" * 50)
        print("  WELCOME TO THE SPACE TRAVEL CALCULATOR!")
        print("  Plan your Hohmann transfer missions within the Solar System.")
        print("#" * 50)
        print("\nCommands:")
        print("  - Type 'list' to see all available celestial bodies.")
        print("  - Type 'exit' to quit the planner at any prompt.")

        while True:
            source_input = input("\nEnter source celestial body (e.g., Earth): ").strip()
            if source_input.lower() == 'exit':
                break
            if source_input.lower() == 'list':
                self._list_available_bodies()
                continue

            target_input = input("Enter target celestial body (e.g., Mars): ").strip()
            if target_input.lower() == 'exit':
                break
            if target_input.lower() == 'list':
                self._list_available_bodies()
                continue

            # Validate input bodies existence
            source_body_data = get_celestial_body_data(source_input)
            target_body_data = get_celestial_body_data(target_input)

            if source_body_data is None:
                print(f"Error: '{source_input}' is not a recognized celestial body. Please try again.")
                continue
            if target_body_data is None:
                print(f"Error: '{target_input}' is not a recognized celestial body. Please try again.")
                continue

            # Use the actual capitalized names from data for consistency in display
            source_body_name = source_body_data['name']
            target_body_name = target_body_data['name']

            print(f"\nInitiating plan for transfer from {source_body_name} to {target_body_name}...")
            trip_details = self.plan_trip(source_body_name, target_body_name)

            if trip_details:
                self.display_trip_summary(trip_details)
            else:
                print("Failed to plan trip. Please review your input bodies and try again.")

        print("\nThank you for using the Space Travel Calculator! Farewell, future astronaut.")

    def _list_available_bodies(self):
        """
        Prints a sorted list of all available celestial bodies in the planner's data.
        """
        print("\n" + "-" * 30)
        print("  Available Celestial Bodies:")
        print("-" * 30)
        
        # Get actual names (case-corrected) and sort them
        body_names = sorted([body_data['name'] for body_data in self._celestial_body_raw_data.values() if 'name' in body_data])
        
        if not body_names:
            print("No celestial bodies found in data.")
            print("-" * 30)
            return

        for i, name in enumerate(body_names):
            print(f"- {name}")
            # Add a pause for long lists, e.g., every 10 items
            if (i + 1) % 10 == 0 and i != len(body_names) - 1:
                input("Press Enter to see more bodies...")
        print("-" * 30)

if __name__ == "__main__":
    # Example usage when this file is run directly
    try:
        planner = SpaceTripPlanner()
        planner.run_interactive_planner()
    except ValueError as e:
        print(f"Initialization Error: {e}")
        print("Please ensure celestial_data.py is correctly configured.")

```