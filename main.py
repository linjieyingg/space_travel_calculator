"""
Course Number: ENGR 13300
Semester: Fall 2024

Description:
    Space travel calculator that calculates time it takes for a user to
    reach a chosen celestial destination, accounting for relativistic effects,
    orbital mechanics, and fuel requirements.

Assignment Information:
    Assignment:     Individual Project
    Team ID:        LC5 - 21
    Author:         Jieying Lin, lin1914@purdue.edu
    Date:           12/01/2024

Contributors:
    Name, login@purdue [repeat for each]

    My contributor(s) helped me:
    [ ] understand the assignment expectations without
        telling me how they will approach it.
    [ ] understand different ways to think about a solution
        without helping me plan my solution.
    [ ] think through the meaning of a specific error or
        bug present in my code without looking at my code.
    Note that if you helped somebody else with their code, you
    have to list that person as a contributor here as well.

Academic Integrity Statement:
    I have not used source code obtained from any unauthorized
    source, either modified or unmodified; nor have I provided
    another student access to my code.  The project I am
    submitting is my own original work.
"""

import math
from datetime import datetime, timedelta
import checks as c
import celestial_data
from orbital_mechanics import calculate_circular_orbital_velocity # Not directly used in main, but good to have if needed
from trajectory_planner import TrajectoryPlanner
import propulsion_system
import fuel_calc
import travel_logger
from constants import C_LIGHT_MPS # Import directly from constants for clarity
import fuel_optimizer # Added import for fuel_optimizer

# Constants
# Conversion factor from years to seconds
YEAR_TO_SECONDS = 365.25 * 24 * 3600

def calc_time_earth(years_traveler_frame, average_speed_ms):
    """
    Calculates travel time in Earth's reference frame using relativistic effects.

    Args:
        years_traveler_frame (float): Travel time in the traveler's (proper) frame in years.
        average_speed_ms (float): The average effective speed of the spacecraft relative to Earth in m/s.
                                  Used to calculate the Lorentz factor for relativistic time dilation.

    Returns:
        float: Travel time in Earth's reference frame in years.
    """
    if average_speed_ms < 0:
        print("Warning: Negative speed provided. Returning traveler's time.")
        return years_traveler_frame
    
    if average_speed_ms == 0:
        # If speed is zero, travel conceptually takes infinite time from Earth's perspective for a journey
        return float('inf')

    if average_speed_ms >= C_LIGHT_MPS:
        # Speed at or above light speed means Lorentz factor is undefined or imaginary.
        # For practical purposes, return traveler's time as a fallback or indicate error.
        print("Warning: Speed is at or above the speed of light. Relativistic time dilation calculation may be invalid or undefined. Returning traveler's time.")
        return years_traveler_frame 

    try:
        lorentz_factor = 1 / math.sqrt(1 - (average_speed_ms / C_LIGHT_MPS)**2)
        time_earth_frame = years_traveler_frame * lorentz_factor
        return time_earth_frame
    except Exception as e:
        print(f"Warning: Error during relativistic time calculation: {e}. Falling back to traveler's time.")
        return years_traveler_frame

def calc_age(years_travel_time: float, starting_age: int) -> float:
    """
    Calculates the user's age upon arrival.

    Args:
        years_travel_time (float): Total travel time in years (traveler's frame).
        starting_age (int): User's starting age.

    Returns:
        float: User's age upon arrival.
    """
    if not isinstance(starting_age, (int, float)) or starting_age < 0:
        raise ValueError("Starting age must be a non-negative number.")
    if not isinstance(years_travel_time, (int, float)) or years_travel_time < 0:
        raise ValueError("Travel time must be a non-negative number.")

    return starting_age + years_travel_time

def calc_arrival(departure_date: datetime, years_earth_frame: float) -> datetime | None:
    """
    Calculates the estimated arrival date based on departure date and travel time
    in Earth's reference frame.

    Args:
        departure_date (datetime): The date of departure.
        years_earth_frame (float): Total travel years in Earth's reference frame.

    Returns:
        datetime.datetime: Estimated arrival date, or None if calculation fails.
    """
    try:
        # Convert years to days for timedelta
        days_earth_frame = years_earth_frame * 365.25 # Account for leap years roughly
        arrival_date = departure_date + timedelta(days=days_earth_frame)
        return arrival_date
    except OverflowError:
        print("Warning: Arrival date calculation resulted in a date out of range due to extremely long travel time.")
        return None
    except TypeError:
        print("Error: Invalid date or time input for arrival calculation.")
        return None
    except Exception as e:
        print(f"Warning: An unexpected error occurred calculating arrival date: {e}")
        return None

def convert_date(date_str: str) -> datetime:
    """
    Converts a date string in YYYY-MM-DD format to a datetime object.

    Args:
        date_str (str): Date string in 'YYYY-MM-DD' format.

    Returns:
        datetime.datetime: Converted datetime object.
    """
    return datetime.strptime(date_str, '%Y-%m-%d')

def main():
    print("Welcome to the Space Travel Calculator!")

    # Initialize TrajectoryPlanner
    planner = TrajectoryPlanner()

    # --- Input Collection ---
    solar_system_destinations = celestial_data.get_all_solar_system_destinations()
    if not solar_system_destinations:
        print("Error: Could not retrieve solar system destinations. Exiting.")
        return

    print("\nAvailable destinations (average distance from Earth):")
    for i, dest in enumerate(solar_system_destinations):
        # Display distances in Astronomical Units or millions of km for better readability
        distance_km = dest['distance'] / 1000 # Convert meters to km
        if distance_km > 1e9: # If distance is greater than a billion km, use billion km
             print(f"{i+1}. {dest['name']} ({distance_km / 1e9:.2f} billion km)")
        elif distance_km > 1e6: # If distance is greater than a million km, use million km
            print(f"{i+1}. {dest['name']} ({distance_km / 1e6:.2f} million km)")
        else: # Otherwise, use km
            print(f"{i+1}. {dest['name']} ({distance_km:.2f} km)")
    
    source_planet_data = None
    while True:
        source_planet_input = input("Enter your source celestial body (e.g., Earth): ").strip()
        source_planet_data = celestial_data.get_celestial_body_data(source_planet_input)
        if source_planet_data:
            source_planet = source_planet_data['name']
            break
        else:
            print("Invalid source planet. Please choose from available bodies (case-insensitive).")
    
    destination_planet_data = None
    while True:
        destination_planet_input = input("Enter your destination celestial body (e.g., Mars): ").strip()
        destination_planet_data = celestial_data.get_celestial_body_data(destination_planet_input)
        if destination_planet_data and destination_planet_data['name'].lower() != source_planet.lower():
            destination_planet = destination_planet_data['name']
            break
        else:
            print("Invalid destination planet or same as source. Please choose a different body (case-insensitive).")

    spacecraft_dry_mass = 0.0
    while True:
        try:
            spacecraft_dry_mass = float(input("Enter spacecraft dry mass (kg): "))
            if spacecraft_dry_mass <= 0:
                raise ValueError
            break
        except ValueError:
            print("Invalid input. Please enter a positive number for spacecraft dry mass.")
        except Exception:
            print("An unexpected error occurred with spacecraft dry mass input. Please re-enter.")


    engine_specific_impulse = 0.0
    while True:
        try:
            engine_specific_impulse = float(input("Enter engine specific impulse (s): "))
            if engine_specific_impulse <= 0:
                raise ValueError
            break
        except ValueError:
            print("Invalid input. Please enter a positive number for engine specific impulse.")
        except Exception:
            print("An unexpected error occurred with engine specific impulse input. Please re-enter.")

    departure_date = None
    while True:
        departure_date_str = input("Enter departure date (YYYY-MM-DD): ")
        if c.is_valid_date(departure_date_str):
            try:
                departure_date = convert_date(departure_date_str)
                break
            except ValueError:
                print("Error: Date format is superficially correct but could not be parsed. Please re-enter a valid date.")
            except Exception:
                print("An unexpected error occurred while processing the date. Please re-enter a valid date.")
        else:
            print("Invalid date format or date. Please use YYYY-MM-DD.")

    user_age = 0
    while True:
        try:
            user_age = int(input("Enter your current age (years): "))
            if c.is_valid_age(user_age):
                break
            else:
                print("Invalid age. Age must be between 18 (inclusive) and 75 (exclusive).")
        except ValueError:
            print("Invalid input. Please enter a whole number for your age.")
        except Exception:
            print("An unexpected error occurred with age input. Please re-enter age.")

    fuel_price_per_unit = 0.0
    while True:
        try:
            fuel_price_per_unit = float(input("Enter fuel price per unit mass (e.g., cost per kg): "))
            if fuel_price_per_unit < 0:
                raise ValueError
            break
        except ValueError:
            print("Invalid input. Please enter a non-negative number for fuel price.")
        except Exception:
            print("An unexpected error occurred with fuel price input. Please re-enter.")

    print("\nCalculating trajectory...")

    # --- Trajectory Planning ---
    # Assume Hohmann Transfer for simplicity as described in trajectory_planner
    trajectory_result = None
    try:
        trajectory_result = planner.plan_hohmann_trajectory(
            departure_body_name=source_planet, 
            arrival_body_name=destination_planet
        )
    except ValueError as e:
        print(f"Error during trajectory planning: {e}.")
        return
    except KeyError as e:
        print(f"Error: Missing expected data key during trajectory planning: {e}.")
        return
    except Exception as e:
        print(f"An unexpected error occurred during trajectory planning: {e}.")
        return

    if not trajectory_result or not trajectory_result.get('success', False):
        print(f"Error planning trajectory: {trajectory_result.get('error_message', 'Unknown error')}")
        return

    # Extract delta-V and travel time from the trajectory planning result
    delta_v_required = trajectory_result['total_delta_v_mps']
    travel_time_seconds = trajectory_result['total_travel_time_seconds']
    transfer_type = trajectory_result['transfer_type']

    # --- Fuel Calculation (using fuel_optimizer) ---
    fuel_mass_needed = 0.0
    try:
        fuel_mass_needed = fuel_optimizer.optimize_fuel_for_trajectory(
            start_body=source_planet,
            end_body=destination_planet,
            trajectory_type=transfer_type, # Use the transfer type from trajectory planning
            spacecraft_dry_mass=spacecraft_dry_mass,
            engine_specific_impulse=engine_specific_impulse
        )
    except ValueError as e:
        print(f"Error optimizing fuel mass: {e}. Fuel calculation skipped.")
        return
    except RuntimeError as e: # Catch RuntimeError as per fuel_optimizer summary
        print(f"Runtime error during fuel optimization: {e}. Fuel calculation skipped.")
        return
    except Exception as e:
        print(f"An unexpected error occurred during fuel optimization: {e}.")
        return

    if fuel_mass_needed is None or fuel_mass_needed < 0:
        print("Error calculating fuel mass or received invalid value from optimizer. Check input parameters.")
        return

    fuel_cost_data = None
    try:
        fuel_cost_data = fuel_calc.calculate_fuel_cost(
            total_fuel_mass_needed=fuel_mass_needed, 
            fuel_price_per_unit=fuel_price_per_unit
        )
    except Exception as e:
        print(f"An unexpected error occurred during fuel cost calculation: {e}.")
        return
    
    total_fuel_cost = 0.0
    if fuel_cost_data and 'total_cost' in fuel_cost_data:
        total_fuel_cost = fuel_cost_data['total_cost']
    else:
        print("Error calculating fuel cost. Check input parameters.")
        return

    # --- Time and Age Calculations ---
    # Convert travel time from seconds to years (traveler's frame is assumed to be the Hohmann transfer time initially)
    travel_time_years_traveler_frame = travel_time_seconds / YEAR_TO_SECONDS
    
    # REVISIT: The relativistic time dilation calculation `calc_time_earth` requires an `average_speed_ms`.
    # For Hohmann transfers, the speed is not constant. Using an 'average' speed for relativistic 
    # effects on a Hohmann trajectory is a significant simplification. 
    # If the intent is to show *some* relativistic effect, we need a representative speed.
    # The `trajectory_planner` provides `average_hohmann_speed_ms` in its output.
    # Let's use the provided average speed from the planner if available, otherwise fallback to a placeholder.
    average_speed_for_relativistic_calc = trajectory_result.get('average_hohmann_speed_ms', 0.0)

    # Fallback to a placeholder if the average speed from planner is not valid or zero.
    if average_speed_for_relativistic_calc <= 0:
        print("Warning: Average Hohmann speed is zero or negative. Using a placeholder speed for relativistic calculation (0.1% of C).")
        average_speed_for_relativistic_calc = 0.001 * C_LIGHT_MPS 

    travel_time_years_earth_frame = calc_time_earth(travel_time_years_traveler_frame, average_speed_for_relativistic_calc)
    
    arrival_age = 0.0
    try:
        arrival_age = calc_age(travel_time_years_traveler_frame, user_age)
    except ValueError as e:
        print(f"Error calculating arrival age: {e}. Age calculation skipped.")
        # arrival_age remains 0.0 or could be set to user_age

    estimated_arrival_date = calc_arrival(departure_date, travel_time_years_earth_frame)

    # --- Output Summary ---
    print("\n--- Travel Summary ---")
    print(f"Source: {source_planet}")
    print(f"Destination: {destination_planet}")
    print(f"Transfer Type: {transfer_type}")
    print(f"Delta-V Required: {delta_v_required:.2f} m/s")
    print(f"Estimated Fuel Mass Needed: {fuel_mass_needed:.2f} kg")
    print(f"Total Fuel Cost: ${total_fuel_cost:,.2f}")
    print(f"Travel Time (Traveler's Frame): {travel_time_years_traveler_frame:.2f} years ({travel_time_seconds / (3600 * 24):.2f} days)")
    if travel_time_years_earth_frame == float('inf'):
        print("Travel Time (Earth's Frame): Effectively infinite due to zero average speed.")
    else:
        print(f"Travel Time (Earth's Frame): {travel_time_years_earth_frame:.2f} years ({travel_time_years_earth_frame * 365.25:.2f} days)")
    print(f"Your Age Upon Arrival: {arrival_age:.0f} years")
    if estimated_arrival_date:
        print(f"Estimated Arrival Date: {estimated_arrival_date.strftime('%Y-%m-%d')}")
    else:
        print("Estimated Arrival Date: Could not be determined")

    # --- Log Travel Details ---
    try:
        travel_logger.save_travel_log(
            source_planet=source_planet,
            destination_planet=destination_planet,
            speed=average_speed_for_relativistic_calc, # Using the speed used for relativistic calculation
            travel_time=travel_time_seconds, # Log original seconds from planner
            delta_v_required=delta_v_required,
            fuel_mass_needed=fuel_mass_needed,
            transfer_type=transfer_type
        )
        print("\nTravel details logged successfully.")
    except Exception as e:
        print(f"\nError saving travel log: {e}")

if __name__ == "__main__":
    main()