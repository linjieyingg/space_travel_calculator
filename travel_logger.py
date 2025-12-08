import json
import datetime
import os
# Assuming FlightPlan class is defined in flight_plan.py
# and is accessible in the same directory or via a configured Python path.
from flight_plan import FlightPlan 

def save_travel_log(flight_plan: FlightPlan):
    """
    Saves detailed space travel calculation details from a FlightPlan object to a JSON log file.

    The log file 'travel_log.json' will store a list of dictionaries,
    each containing a timestamp, origin body, destination body, spacecraft speed,
    travel time (traveler's frame), required delta-v, estimated fuel mass, and transfer type.

    If the log file does not exist, it will be created.
    If it exists but is empty or malformed, it will be initialized or overwritten
    in a way that new data can be appended.

    Args:
        flight_plan (FlightPlan): An object containing all details of the planned space travel.
                                  Expected attributes include:
                                  - origin_body (str)
                                  - destination_body (str)
                                  - speed_Mm_s (float)
                                  - travel_time_years (float)
                                  - delta_v_required_km_s (float)
                                  - fuel_mass_needed_kg (float)
                                  - transfer_type (str)
    """
    log_file = 'travel_log.json'
    log_data = []

    # Get current timestamp
    timestamp = datetime.datetime.now().isoformat()

    # Create the new entry with detailed orbital trajectory planning information
    # Extract data directly from the FlightPlan object
    new_entry = {
        'timestamp': timestamp,
        'origin_body': flight_plan.origin_body,
        'destination_body': flight_plan.destination_body,
        'speed_Mm_s': flight_plan.speed_Mm_s,
        'travel_time_years': flight_plan.travel_time_years,
        'delta_v_required_km_s': flight_plan.delta_v_required_km_s,
        'fuel_mass_needed_kg': flight_plan.fuel_mass_needed_kg,
        'transfer_type': flight_plan.transfer_type
    }

    # Load existing data if the file exists and is valid JSON
    if os.path.exists(log_file):
        try:
            with open(log_file, 'r') as f:
                content = f.read().strip()
                if content: # Check if file is not empty
                    log_data = json.loads(content)
                else:
                    log_data = [] # File is empty, start with an empty list
        except json.JSONDecodeError:
            # Handle cases where the file might be corrupted or malformed JSON
            print(f"Warning: {log_file} is corrupted or not valid JSON. Starting a new log.")
            log_data = []
        except Exception as e:
            print(f"An unexpected error occurred while reading {log_file}: {e}")
            log_data = [] # Fallback to empty list on other errors

    # Ensure log_data is always a list for appending
    if not isinstance(log_data, list):
        print(f"Warning: {log_file} content is not a list. Resetting log.")
        log_data = []

    # Append the new entry
    log_data.append(new_entry)

    # Write the updated data back to the file
    with open(log_file, 'w') as f:
        json.dump(log_data, f, indent=4)