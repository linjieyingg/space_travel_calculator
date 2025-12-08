"""
Course Number: ENGR 13300
Semester: Fall 2024

Description:
    Space travel calculator that calculates time it takes for a user to
    reach a chosen celestial destination, accounting for relativistic effects,
    orbital mechanics, and fuel requirements. This refactored version
    leverages a structured FlightPlan data model and delegates core
    calculation logic to dedicated modules.

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
import propulsion_system # Imported but not directly called in main after fuel_optimizer was added
import fuel_calc
import travel_logger
from constants import C_LIGHT_MPS # Import directly from constants for clarity
import fuel_optimizer # Added import for fuel_optimizer
from flight_plan import FlightPlan # New: Import FlightPlan dataclass
import flight_data_processor as fdp # New: Import flight data processing functions

# Constants
# Conversion factor from years to seconds
YEAR_TO_SECONDS = 365.25 * 24 * 3600

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
        source_planet_input = input("Enter your Departure Body (e.g., Earth): ").strip()
        source_planet_data = celestial_data.get_celestial_body_data(source_planet_input)
        if source_planet_data:
            source_planet = source_planet_data['name']
            break
        else:
            print("Invalid departure body. Please choose from available bodies (case-insensitive).")
    
    destination_planet_data = None
    while True:
        destination_planet_input = input("Enter your Arrival Body (e.g., Mars): ").strip()
        destination_planet_data = celestial_data.get_celestial_body_data(destination_planet_input)
        if destination_planet_data and destination_planet_data['name'].lower() != source_planet.lower():
            destination_planet = destination_planet_data['name']
            break
        else:
            print("Invalid arrival body or same as departure. Please choose a different body (case-insensitive).")

    # New: Trajectory Type Selection
    supported_trajectory_types = ['Hohmann', 'Direct Transfer'] # 'Direct Transfer' is for future expansion
    selected_trajectory_type = 'Hohmann' # Default to Hohmann
    while True:
        type_input = input(f"Enter trajectory type (e.g., {', '.join(supported_trajectory_types)}). Default is Hohmann: ").strip()
        if not type_input:
            print(f"No trajectory type entered. Defaulting to '{selected_trajectory_type}'.")
            break
        
        # Normalize input for comparison
        normalized_input = type_input.replace(' ', '').lower() 
        
        if normalized_input == 'hohmann':
            selected_trajectory_type = 'Hohmann'
            break
        elif normalized_input == 'directtransfer': # For future proofing 'Direct Transfer'
            selected_trajectory_type = 'Direct Transfer'
            print(f"Note: '{selected_trajectory_type}' is selected, but may not be fully implemented and will default to Hohmann behavior if not supported by the planner.")
            break
        else:
            print(f"Invalid trajectory type. Please choose from {', '.join(supported_trajectory_types)}.")


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
                departure_date = datetime.strptime(departure_date_str, '%Y-%m-%d')
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
    trajectory_result = None
    try:
        trajectory_result = planner.plan_trajectory(
            departure_body_name=source_planet, 
            arrival_body_name=destination_planet,
            trajectory_type=selected_trajectory_type, # Pass the selected trajectory type
            departure_date=departure_date 
        )
    except ValueError as e:
        print(f"Error during trajectory planning: {e}. Please ensure valid departure and arrival bodies are selected and their data is complete.")
        return
    except KeyError as e:
        print(f"Error: Essential data (like Sun's gravitational parameter) is missing for trajectory planning: {e}. Please check the celestial data configuration.")
        return
    except Exception as e:
        print(f"An unexpected error occurred during trajectory planning: {e}. Please review your inputs or contact support if the issue persists.")
        return

    # Ensure trajectory_result is valid before parsing
    if not trajectory_result or not trajectory_result.get('success', False):
        print(f"Error planning trajectory: {trajectory_result.get('error_message', 'Unknown error')}. Please verify your selected celestial bodies.")
        return

    # Extract delta-V and travel time from the trajectory planning result using .get() for robustness
    delta_v_required = trajectory_result.get('total_delta_v_mps', 0.0)
    travel_time_seconds = trajectory_result.get('total_travel_time_seconds', 0.0)
    transfer_type = trajectory_result.get('transfer_type', 'Unknown Transfer')

    # Add explicit check for critical values being non-positive, as they lead to issues downstream
    if delta_v_required <= 0 and travel_time_seconds <= 0:
        print("Critical error: Trajectory planning failed to provide valid Delta-V and travel time. Exiting.")
        return
    elif delta_v_required <= 0:
        print("Warning: Calculated Delta-V is non-positive. This might indicate an error in trajectory planning or no required burn.")
        # Continue as fuel_optimizer might handle zero delta-v
    elif travel_time_seconds <= 0:
        print("Warning: Calculated travel time is non-positive. This might indicate an error in trajectory planning. Using a default small positive time for further calculations.")
        travel_time_seconds = 1.0 # Avoid division by zero or other issues downstream if 0 or negative.


    # --- Fuel Calculation (using fuel_optimizer) ---
    fuel_mass_needed = 0.0
    try:
        fuel_mass_needed = fuel_optimizer.optimize_fuel_for_trajectory(
            start_body=source_planet,
            end_body=destination_planet,
            trajectory_type=selected_trajectory_type, # Pass the user's selected trajectory type
            spacecraft_dry_mass=spacecraft_dry_mass,
            engine_specific_impulse=engine_specific_impulse,
            departure_date=departure_date # Pass the departure date to fuel optimizer for trajectory planning
        )
    except ValueError as e:
        print(f"Input error during fuel optimization: {e}. Please ensure spacecraft dry mass and engine specific impulse are positive numbers. Fuel calculation skipped.")
        return
    except RuntimeError as e: # Catch RuntimeError as per fuel_optimizer summary
        print(f"Operational error during fuel optimization: {e}. This might indicate an issue with internal calculations or mock components. Fuel calculation skipped.")
        return
    except Exception as e:
        print(f"An unexpected error occurred during fuel optimization: {e}. Fuel calculation could not be completed. Please review your inputs.")
        return

    if fuel_mass_needed is None or fuel_mass_needed < 0:
        print("Fuel mass calculation failed or resulted in an invalid value. This might occur if the required Delta-V is non-positive or other parameters are out of range. Please review your spacecraft and engine inputs.")
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
        print("Error calculating fuel cost. Check input parameters or the result from fuel mass calculation.")
        return

    # --- Time and Age Calculations (using flight_data_processor) ---
    travel_time_years_traveler_frame = travel_time_seconds / YEAR_TO_SECONDS
    
    average_speed_for_relativistic_calc_mps = trajectory_result.get('average_hohmann_speed_ms', 0.0)

    # Fallback to a placeholder if the average speed from planner is not valid or zero.
    if average_speed_for_relativistic_calc_mps <= 0:
        print("Warning: Average Hohmann speed is zero or negative. Using a placeholder speed for relativistic calculation (0.1% of C).")
        average_speed_for_relativistic_calc_mps = 0.001 * C_LIGHT_MPS 

    travel_time_years_earth_frame = 0.0
    try:
        travel_time_years_earth_frame = fdp.calculate_earth_frame_time(
            travel_time_years_traveler_frame, average_speed_for_relativistic_calc_mps
        )
    except Exception as e:
        print(f"Error calculating Earth frame time: {e}. Defaulting to traveler's frame time.")
        travel_time_years_earth_frame = travel_time_years_traveler_frame


    arrival_age_years = 0.0
    try:
        arrival_age_years = fdp.calculate_arrival_age(travel_time_years_traveler_frame, user_age)
    except ValueError as e:
        print(f"Error calculating arrival age: {e}. Arrival age calculation skipped.")
        arrival_age_years = float(user_age) # Set to starting age as a fallback

    estimated_arrival_date = None
    try:
        estimated_arrival_date = fdp.calculate_arrival_date(departure_date, travel_time_years_earth_frame)
    except Exception as e:
        print(f"Error calculating estimated arrival date: {e}.")
        # estimated_arrival_date remains None

    # --- Assemble FlightPlan Object ---
    flight_plan = FlightPlan(
        source_planet=source_planet,
        destination_planet=destination_planet,
        selected_trajectory_type=selected_trajectory_type,
        calculated_transfer_type=transfer_type,
        departure_date=departure_date,
        spacecraft_dry_mass_kg=spacecraft_dry_mass,
        engine_specific_impulse_s=engine_specific_impulse,
        user_starting_age_years=user_age,
        fuel_price_per_unit=fuel_price_per_unit,
        delta_v_required_mps=delta_v_required,
        fuel_mass_needed_kg=fuel_mass_needed,
        total_fuel_cost=total_fuel_cost,
        average_speed_for_relativistic_calc_mps=average_speed_for_relativistic_calc_mps,
        travel_time_traveler_frame_years=travel_time_years_traveler_frame,
        travel_time_earth_frame_years=travel_time_years_earth_frame,
        arrival_age_years=arrival_age_years,
        estimated_arrival_date=estimated_arrival_date
    )

    # --- Output Summary ---
    print("\n--- Travel Summary ---")
    print(f"Departure Body: {flight_plan.source_planet}")
    print(f"Arrival Body: {flight_plan.destination_planet}")
    print(f"Trajectory Type Selected: {flight_plan.selected_trajectory_type}") 
    print(f"Calculated Transfer Type: {flight_plan.calculated_transfer_type}")
    print(f"Delta-V Required: {flight_plan.delta_v_required_mps:.2f} m/s")
    print(f"Estimated Fuel Mass Needed: {flight_plan.fuel_mass_needed_kg:.2f} kg")
    print(f"Total Fuel Cost: ${flight_plan.total_fuel_cost:,.2f}")
    print(f"Travel Time (Traveler's Frame): {flight_plan.travel_time_traveler_frame_years:.2f} years ({flight_plan.travel_time_traveler_frame_years * YEAR_TO_SECONDS / (3600 * 24):.2f} days)")
    if flight_plan.travel_time_earth_frame_years == float('inf'):
        print("Travel Time (Earth's Frame): Effectively infinite due to zero average speed.")
    else:
        print(f"Travel Time (Earth's Frame): {flight_plan.travel_time_earth_frame_years:.2f} years ({flight_plan.travel_time_earth_frame_years * 365.25:.2f} days)")
    print(f"Your Age Upon Arrival: {flight_plan.arrival_age_years:.0f} years")
    if flight_plan.estimated_arrival_date:
        print(f"Estimated Arrival Date: {flight_plan.estimated_arrival_date.strftime('%Y-%m-%d')}")
    else:
        print("Estimated Arrival Date: Could not be determined")

    # --- Log Travel Details ---
    try:
        # Assuming travel_logger.save_travel_log has been updated to accept a FlightPlan object
        travel_logger.save_travel_log(flight_plan)
        print("\nTravel details logged successfully.")
    except Exception as e:
        print(f"\nError saving travel log: {e}")

if __name__ == "__main__":
    main()