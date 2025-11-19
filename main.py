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

from datetime import datetime, timedelta
import checks as c
import math
import celestial_data # Provides celestial body data
import orbital_mechanics # Foundational for trajectory_planner
import trajectory_planner # Calculates orbital trajectories
import propulsion_system # Calculates fuel requirements
import fuel_calc # For fuel cost calculation
import travel_logger # For logging travel details
import constants # Provides universal physical constants

# Constants
FUEL_PRICE_PER_UNIT = 1000000.0 # Example high price for space fuel per kg/unit

# List of solar system bodies for user selection (simplified for example)
SOLAR_SYSTEM_BODIES = [
    "Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"
]

def calc_time_earth(years_traveler_frame, average_speed_ms):
    """
    Calculates the travel time in Earth's reference frame, considering relativistic
    time dilation.

    Args:
        years_traveler_frame (float): Travel time in the traveler's (proper) frame in years.
        average_speed_ms (float): The average effective speed of the spacecraft relative to Earth in m/s.
                                  This speed is used to calculate the Lorentz factor.

    Returns:
        float: Travel time in Earth's reference frame in years.
    """
    if average_speed_ms <= 0:
        # If speed is non-positive, no relativistic effect, or conceptually impossible travel.
        # Return traveler's time as a fallback.
        return years_traveler_frame
    
    # Avoid division by zero or negative square root if speed is too close or exceeds light speed
    if average_speed_ms >= constants.C_LIGHT_MPS:
        print("Warning: Speed is at or above the speed of light. Relativistic time dilation calculation may be invalid or undefined.")
        # In a theoretical scenario, if v=c, Lorentz factor is undefined. If v>c, sqrt is imaginary.
        # For practical purposes, we might cap or return a large value, or the traveler's time.
        return years_traveler_frame # Fallback

    try:
        # Calculate the Lorentz factor term (1 - v^2/c^2)
        v_c_squared = (average_speed_ms / constants.C_LIGHT_MPS)**2
        sqrt_term_arg = 1 - v_c_squared

        if sqrt_term_arg <= 0: # Should be caught by the average_speed_ms >= C_LIGHT_MPS check, but good defensive programming
            print("Warning: Relativistic term argument is non-positive. Returning traveler's time.")
            return years_traveler_frame
        
        # Calculate time in Earth's reference frame (coordinate time)
        # Delta_t_prime = Delta_t / sqrt(1 - v^2/c^2)
        years_earth_frame = years_traveler_frame / math.sqrt(sqrt_term_arg)
        
        return years_earth_frame
    except Exception as e:
        print(f"Warning: Error during relativistic time calculation: {e}. Falling back to traveler's time.")
        return years_traveler_frame

def calc_age(years, age):
    """
    Calculates the user's age upon arrival.

    Args:
        years (float): Travel time in the traveler's frame.
        age (int): User's starting age.

    Returns:
        int: User's age upon arrival, floored to a whole number.
    """
    new_age = math.floor(age + years)
    return new_age

def calc_arrival(date, years_earth_frame):
    """
    Calculates the estimated arrival date on Earth's calendar.

    Args:
        date (datetime.datetime): Departure date.
        years_earth_frame (float): Total travel years in Earth's reference frame.

    Returns:
        datetime.datetime: Estimated arrival date, or None if calculation fails.
    """
    try:
        days = years_earth_frame * 365.25
        td = timedelta(days=days)
        new_date = date + td
        return new_date
    except OverflowError:
        print("Warning: Arrival date calculation resulted in a date out of range due to extremely long travel time.")
        return None
    except Exception as e:
        print(f"Warning: An unexpected error occurred calculating arrival date: {e}")
        return None

def convert_date(date_str):
    """
    Converts a date string in 'YYYY-MM-DD' format to a datetime object.

    Args:
        date_str (str): Date string.

    Returns:
        datetime.datetime: Converted datetime object.
    """
    date = datetime.strptime(date_str, '%Y-%m-%d')
    return date

def main():
    print("Welcome to the Interplanetary Travel Calculator!")

    # --- User Input for Trajectory and Spacecraft ---

    # Source Planet Input
    source_data = None
    while True:
        source_planet = input(f"Enter your source planet ({', '.join(SOLAR_SYSTEM_BODIES)}): ").strip().title()
        if source_planet in SOLAR_SYSTEM_BODIES:
            try:
                source_data = celestial_data.get_celestial_body_data(source_planet)
                if source_data: # Ensure data was successfully retrieved
                    break
                else:
                    print(f"Error: Data for {source_planet} not found in celestial_data. Please try again.")
            except Exception as e:
                print(f"An unexpected error occurred retrieving data for {source_planet}: {e}. Please try again.")
        else:
            print(f"Invalid source planet. Please choose from: {', '.join(SOLAR_SYSTEM_BODIES)}")

    # Destination Planet Input
    destination_data = None
    while True:
        destination_planet = input(f"Enter your destination planet ({', '.join(SOLAR_SYSTEM_BODIES)}, excluding {source_planet}): ").strip().title()
        if destination_planet == source_planet:
            print("Destination planet cannot be the same as the source planet. Please choose a different one.")
            continue
        if destination_planet in SOLAR_SYSTEM_BODIES:
            try:
                destination_data = celestial_data.get_celestial_body_data(destination_planet)
                if destination_data: # Ensure data was successfully retrieved
                    break
                else:
                    print(f"Error: Data for {destination_planet} not found in celestial_data. Please try again.")
            except Exception as e:
                print(f"An unexpected error occurred retrieving data for {destination_planet}: {e}. Please try again.")
        else:
            print(f"Invalid destination planet. Please choose from: {', '.join(SOLAR_SYSTEM_BODIES)}")

    # Spacecraft Dry Mass Input
    while True:
        try:
            spacecraft_dry_mass = float(input("Enter spacecraft dry mass in kilograms (kg): "))
            if spacecraft_dry_mass > 0:
                break
            else:
                print("Spacecraft dry mass must be a positive number. Please re-enter.")
        except ValueError:
            print("Invalid input. Spacecraft dry mass must be a numeric value. Please re-enter.")
        except Exception:
            print("An unexpected error occurred with spacecraft mass input. Please re-enter.")

    # Engine Specific Impulse (ISP) Input
    while True:
        try:
            engine_isp = float(input("Enter engine specific impulse (ISP) in seconds (s): "))
            if engine_isp > 0:
                break
            else:
                print("Engine specific impulse must be a positive number. Please re-enter.")
        except ValueError:
            print("Invalid input. Engine ISP must be a numeric value. Please re-enter.")
        except Exception:
            print("An unexpected error occurred with engine ISP input. Please re-enter.")

    # Departure Date Input
    while True:
        date_str = input("Date of departure (YYYY-MM-DD): ")
        if c.is_valid_date(date_str):
            try:
                departure_date = convert_date(date_str)
                break
            except ValueError:
                print("Error: Date format is superficially correct but could not be parsed. Please re-enter a valid date.")
            except Exception:
                print("An unexpected error occurred while processing the date. Please re-enter a valid date.")
        else:
            print("Error with date format. Re-enter valid date.")

    # User Age Input
    while True:
        try:
            user_age = int(input("Enter your age in years: "))
            if c.is_valid_age(user_age):
                break
            else:
                print("Enter a valid age between 18 and 74.")
        except ValueError:
            print("Invalid input. Age must be a whole number. Re-enter age.")
        except Exception:
            print("An unexpected error occurred with age input. Please re-enter age.")
            
    print("\n--- Calculating Trajectory and Fuel Requirements ---")

    # --- Trajectory Planning ---
    delta_v = 0.0
    travel_time_traveler_frame_years = 0.0
    average_hohmann_speed_ms = 0.0
    transfer_type = "Hohmann Transfer" # Default transfer type

    try:
        planner = trajectory_planner.TrajectoryPlanner()
        hohmann_result = planner.plan_hohmann_trajectory(source_planet, destination_planet)
        
        # Check if planning was successful
        if hohmann_result and hohmann_result.get('success', False):
            delta_v = hohmann_result.get('total_delta_v_mps', 0.0)
            # These are expected to be returned by the trajectory planner as per main.py's existing logic
            travel_time_traveler_frame_years = hohmann_result.get('travel_time_traveler_frame_years', 0.0)
            average_hohmann_speed_ms = hohmann_result.get('average_hohmann_speed_ms', 0.0)

            print(f"Calculated Delta-V for {source_planet} to {destination_planet} {transfer_type}: {delta_v:.3f} m/s")
        else:
            error_message = hohmann_result.get('error', 'Unknown error during trajectory planning.') if hohmann_result else 'Trajectory planner returned no result.'
            print(f"Error during trajectory planning: {error_message}. Cannot proceed with calculations.")
            return # Exit if trajectory planning fails

    except ValueError as e:
        print(f"Error during trajectory planning: {e}. Cannot proceed with calculations.")
        return # Exit if trajectory planning fails
    except KeyError as e:
        print(f"Error: Missing expected data key from trajectory planner: {e}. Cannot proceed.")
        return
    except Exception as e:
        print(f"An unexpected error occurred during trajectory planning: {e}. Cannot proceed.")
        return

    # --- Fuel Calculation ---
    fuel_mass_kg = 0.0
    total_fuel_cost = 0.0
    if delta_v > 0 and spacecraft_dry_mass > 0 and engine_isp > 0:
        try:
            fuel_mass_kg = propulsion_system.calculate_required_fuel_mass(delta_v, spacecraft_dry_mass, engine_isp)
            
            fuel_cost_result = fuel_calc.calculate_fuel_cost(total_fuel_mass_needed=fuel_mass_kg, fuel_price_per_unit=FUEL_PRICE_PER_UNIT)
            if fuel_cost_result and 'total_cost' in fuel_cost_result:
                total_fuel_cost = fuel_cost_result['total_cost']
            else:
                print("Warning: Fuel cost calculation failed or returned invalid result.")
                total_fuel_cost = 0.0 # Default to 0 if calculation fails

            print(f"Required Fuel Mass: {fuel_mass_kg:.3f} kg")
            print(f"Estimated Fuel Cost: ${total_fuel_cost:,.2f}")
        except ValueError as e:
            print(f"Error calculating fuel mass or cost: {e}.")
        except Exception as e:
            print(f"An unexpected error occurred during fuel calculation: {e}.")
    else:
        print("Skipping fuel calculation due to invalid Delta-V, spacecraft mass, or engine ISP.")

    # --- Time and Age Calculations ---
    years_earth_frame = calc_time_earth(travel_time_traveler_frame_years, average_hohmann_speed_ms)
    arrival_age = calc_age(travel_time_traveler_frame_years, user_age)
    arrival_date = calc_arrival(departure_date, years_earth_frame)

    # --- Output Results ---
    print("\n--- Travel Summary ---")
    print(f"Source Planet: {source_planet}")
    print(f"Destination: {destination_planet}")
    print(f"Transfer Type: {transfer_type}")
    print(f"Required Delta-V: {delta_v:.3f} m/s")
    print(f"Required Fuel Mass: {fuel_mass_kg:.3f} kg")
    print(f"Estimated Fuel Cost: ${total_fuel_cost:,.2f}")
    print(f"Travel Time (Traveler's Frame): {travel_time_traveler_frame_years:.3g} years")
    print(f"Travel Time (Earth's Frame): {years_earth_frame:.5g} years")
    print(f"Your Age Upon Arrival: {arrival_age} years")
    if arrival_date is not None:
        print(f"Estimated Arrival Date (Earth Calendar): {arrival_date:%m/%d/%Y}")
    else:
        print("Estimated Arrival Date: Could not be determined due to extreme duration or date limits.")

    # --- Log Travel Details ---
    try:
        travel_logger.save_travel_log(
            source_planet=source_planet,
            destination_planet=destination_planet,
            speed=average_hohmann_speed_ms, # Using average speed for logging
            travel_time=travel_time_traveler_frame_years, # Travel time in traveler's frame (years) for logging
            delta_v_required=delta_v,
            fuel_mass_needed=fuel_mass_kg,
            transfer_type=transfer_type
        )
        print("\nTravel details successfully logged.")
    except TypeError as e: # Catch TypeError if signature is truly mismatched
        print(f"\nWarning: travel_logger.save_travel_log signature mismatch: {e}. Logging failed.")
    except Exception as e:
        print(f"\nError saving travel log: {e}")

if __name__ == "__main__":
    main()