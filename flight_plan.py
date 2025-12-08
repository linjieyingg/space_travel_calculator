import datetime
from dataclasses import dataclass, field
from typing import Union

@dataclass(frozen=True)
class FlightPlan:
    """
    Encapsulates all details of a planned space journey, providing a single,
    structured representation of a complete flight plan.

    This dataclass is designed to be immutable (frozen=True) as it represents
    a finalized plan, ensuring data integrity once an instance is created.

    Attributes:
        # User Inputs
        departure_body_name (str): The name of the celestial body from which the journey starts.
        arrival_body_name (str): The name of the celestial body at which the journey ends.
        trajectory_type (str): The type of trajectory planned (e.g., 'Hohmann', 'Direct').
        spacecraft_dry_mass_kg (float): The dry mass of the spacecraft in kilograms.
        engine_specific_impulse_s (float): The specific impulse of the spacecraft engine in seconds.
        departure_date (datetime.datetime): The Earth-frame date and time of departure.
        traveler_starting_age_years (int): The starting age of the traveler in years.
        fuel_price_per_unit (float): The cost of fuel per unit mass (e.g., per kg).

        # Trajectory Specifics
        total_delta_v_mps (float): The total change in velocity required for the trajectory in meters per second.
        travel_time_traveler_frame_years (float): The duration of travel as experienced by the traveler in years.
        travel_time_earth_frame_years (float): The duration of travel as observed from Earth in years.
        departure_distance_from_sun_m (float): The distance of the departure body from the Sun at departure time in meters.
        arrival_distance_from_sun_m (float): The distance of the arrival body from the Sun at arrival time in meters.

        # Fuel Calculations
        total_fuel_mass_needed_kg (float): The total mass of fuel required for the journey in kilograms.
        total_fuel_cost (float): The total cost of the required fuel.

        # Relativistic Effects
        traveler_arrival_age_years (float): The traveler's age upon arrival at the destination (traveler's frame).
        arrival_date_earth_frame (datetime.datetime): The Earth-frame date and time of arrival.
    """
    # User Inputs
    departure_body_name: str
    arrival_body_name: str
    trajectory_type: str
    spacecraft_dry_mass_kg: float
    engine_specific_impulse_s: float
    departure_date: datetime.datetime
    traveler_starting_age_years: int
    fuel_price_per_unit: float

    # Trajectory Specifics
    total_delta_v_mps: float
    travel_time_traveler_frame_years: float
    travel_time_earth_frame_years: float
    departure_distance_from_sun_m: float
    arrival_distance_from_sun_m: float

    # Fuel Calculations
    total_fuel_mass_needed_kg: float
    total_fuel_cost: float

    # Relativistic Effects
    traveler_arrival_age_years: float
    arrival_date_earth_frame: datetime.datetime

    def __post_init__(self):
        """
        Performs post-initialization validation to ensure the logical integrity
        of the flight plan data.
        """
        # --- Validate User Inputs ---
        if not isinstance(self.departure_body_name, str) or not self.departure_body_name:
            raise ValueError("Departure body name must be a non-empty string.")
        if not isinstance(self.arrival_body_name, str) or not self.arrival_body_name:
            raise ValueError("Arrival body name must be a non-empty string.")
        if not isinstance(self.trajectory_type, str) or not self.trajectory_type:
            raise ValueError("Trajectory type must be a non-empty string.")

        if not isinstance(self.spacecraft_dry_mass_kg, (int, float)) or self.spacecraft_dry_mass_kg <= 0:
            raise ValueError("Spacecraft dry mass must be a positive number.")
        if not isinstance(self.engine_specific_impulse_s, (int, float)) or self.engine_specific_impulse_s <= 0:
            raise ValueError("Engine specific impulse must be a positive number.")
        if not isinstance(self.departure_date, datetime.datetime):
            raise TypeError("Departure date must be a datetime.datetime object.")
        if not isinstance(self.traveler_starting_age_years, int) or self.traveler_starting_age_years < 0:
            raise ValueError("Traveler starting age must be a non-negative integer.")
        if not isinstance(self.fuel_price_per_unit, (int, float)) or self.fuel_price_per_unit < 0:
            raise ValueError("Fuel price per unit must be a non-negative number.")

        # --- Validate Trajectory Specifics ---
        if not isinstance(self.total_delta_v_mps, (int, float)) or self.total_delta_v_mps < 0:
            raise ValueError("Total Delta-V must be a non-negative number.")
        if not isinstance(self.travel_time_traveler_frame_years, (int, float)) or self.travel_time_traveler_frame_years < 0:
            raise ValueError("Travel time (traveler frame) must be a non-negative number.")
        if not isinstance(self.travel_time_earth_frame_years, (int, float)) or self.travel_time_earth_frame_years < 0:
            raise ValueError("Travel time (Earth frame) must be a non-negative number.")
        if not isinstance(self.departure_distance_from_sun_m, (int, float)) or self.departure_distance_from_sun_m < 0:
            raise ValueError("Departure distance from Sun must be a non-negative number.")
        if not isinstance(self.arrival_distance_from_sun_m, (int, float)) or self.arrival_distance_from_sun_m < 0:
            raise ValueError("Arrival distance from Sun must be a non-negative number.")

        # --- Validate Fuel Calculations ---
        if not isinstance(self.total_fuel_mass_needed_kg, (int, float)) or self.total_fuel_mass_needed_kg < 0:
            raise ValueError("Total fuel mass needed must be a non-negative number.")
        if not isinstance(self.total_fuel_cost, (int, float)) or self.total_fuel_cost < 0:
            raise ValueError("Total fuel cost must be a non-negative number.")

        # --- Validate Relativistic Effects ---
        if not isinstance(self.traveler_arrival_age_years, (int, float)) or self.traveler_arrival_age_years < self.traveler_starting_age_years:
            # Arrival age can be same as starting age for 0 travel time.
            raise ValueError("Traveler arrival age cannot be less than starting age.")
        if not isinstance(self.arrival_date_earth_frame, datetime.datetime):
            raise TypeError("Arrival date (Earth frame) must be a datetime.datetime object.")
        if self.arrival_date_earth_frame < self.departure_date:
             raise ValueError("Arrival date cannot be before departure date.")

    def __str__(self) -> str:
        """
        Returns a user-friendly string representation of the flight plan.
        """
        return (
            f"Flight Plan from {self.departure_body_name} to {self.arrival_body_name} "
            f"({self.trajectory_type} Trajectory)\n"
            f"  Departure Date: {self.departure_date.strftime('%Y-%m-%d')}\n"
            f"  Delta-V Required: {self.total_delta_v_mps:.2f} m/s\n"
            f"  Travel Time (Traveler): {self.travel_time_traveler_frame_years:.2f} years\n"
            f"  Travel Time (Earth): {self.travel_time_earth_frame_years:.2f} years\n"
            f"  Fuel Mass: {self.total_fuel_mass_needed_kg:.2f} kg\n"
            f"  Total Cost: ${self.total_fuel_cost:.2f}\n"
            f"  Traveler's Age on Arrival: {self.traveler_arrival_age_years:.2f} years\n"
            f"  Arrival Date (Earth Frame): {self.arrival_date_earth_frame.strftime('%Y-%m-%d')}"
        )

# Example Usage (for demonstration, not part of the production code itself)
if __name__ == '__main__':
    try:
        # Valid Flight Plan
        plan1 = FlightPlan(
            departure_body_name="Earth",
            arrival_body_name="Mars",
            trajectory_type="Hohmann",
            spacecraft_dry_mass_kg=1000.0,
            engine_specific_impulse_s=450.0,
            departure_date=datetime.datetime(2025, 7, 15, 10, 0, 0),
            traveler_starting_age_years=30,
            fuel_price_per_unit=1500.0,
            total_delta_v_mps=15000.0,
            travel_time_traveler_frame_years=0.7,
            travel_time_earth_frame_years=0.7000000001, # Slightly different due to time dilation
            departure_distance_from_sun_m=1.496e11, # ~1 AU
            arrival_distance_from_sun_m=2.279e11, # ~1.52 AU
            total_fuel_mass_needed_kg=25000.0,
            total_fuel_cost=37500000.0,
            traveler_arrival_age_years=30.7,
            arrival_date_earth_frame=datetime.datetime(2026, 3, 23, 10, 0, 0)
        )
        print("--- Valid Flight Plan 1 ---")
        print(plan1)

        # Another Valid Flight Plan (Direct Transfer)
        plan2 = FlightPlan(
            departure_body_name="Earth",
            arrival_body_name="Jupiter",
            trajectory_type="Direct",
            spacecraft_dry_mass_kg=5000.0,
            engine_specific_impulse_s=800.0,
            departure_date=datetime.datetime(2030, 1, 1),
            traveler_starting_age_years=40,
            fuel_price_per_unit=2000.0,
            total_delta_v_mps=30000.0,
            travel_time_traveler_frame_years=2.5,
            travel_time_earth_frame_years=2.5000001,
            departure_distance_from_sun_m=1.496e11,
            arrival_distance_from_sun_m=7.785e11,
            total_fuel_mass_needed_kg=80000.0,
            total_fuel_cost=160000000.0,
            traveler_arrival_age_years=42.5,
            arrival_date_earth_frame=datetime.datetime(2032, 7, 1, 0, 0, 0)
        )
        print("\n--- Valid Flight Plan 2 ---")
        print(plan2)


        # Invalid Flight Plan (negative dry mass)
        print("\n--- Attempting Invalid Flight Plan (negative dry mass) ---")
        try:
            invalid_plan = FlightPlan(
                departure_body_name="Earth",
                arrival_body_name="Mars",
                trajectory_type="Hohmann",
                spacecraft_dry_mass_kg=-100.0, # Invalid!
                engine_specific_impulse_s=450.0,
                departure_date=datetime.datetime(2025, 7, 15),
                traveler_starting_age_years=30,
                fuel_price_per_unit=1500.0,
                total_delta_v_mps=15000.0,
                travel_time_traveler_frame_years=0.7,
                travel_time_earth_frame_years=0.7000000001,
                departure_distance_from_sun_m=1.496e11,
                arrival_distance_from_sun_m=2.279e11,
                total_fuel_mass_needed_kg=25000.0,
                total_fuel_cost=37500000.0,
                traveler_arrival_age_years=30.7,
                arrival_date_earth_frame=datetime.datetime(2026, 3, 23)
            )
            print(invalid_plan)
        except ValueError as e:
            print(f"Caught expected error: {e}")

        # Invalid Flight Plan (arrival date before departure)
        print("\n--- Attempting Invalid Flight Plan (arrival date before departure) ---")
        try:
            invalid_plan_date = FlightPlan(
                departure_body_name="Earth",
                arrival_body_name="Mars",
                trajectory_type="Hohmann",
                spacecraft_dry_mass_kg=1000.0,
                engine_specific_impulse_s=450.0,
                departure_date=datetime.datetime(2025, 7, 15),
                traveler_starting_age_years=30,
                fuel_price_per_unit=1500.0,
                total_delta_v_mps=15000.0,
                travel_time_traveler_frame_years=0.7,
                travel_time_earth_frame_years=0.7000000001,
                departure_distance_from_sun_m=1.496e11,
                arrival_distance_from_sun_m=2.279e11,
                total_fuel_mass_needed_kg=25000.0,
                total_fuel_cost=37500000.0,
                traveler_arrival_age_years=30.7,
                arrival_date_earth_frame=datetime.datetime(2024, 3, 23) # Invalid!
            )
            print(invalid_plan_date)
        except ValueError as e:
            print(f"Caught expected error: {e}")

    except Exception as e:
        print(f"An unexpected error occurred during testing: {e}")