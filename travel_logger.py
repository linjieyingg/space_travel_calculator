import json
import datetime
import os

def save_travel_log(source_planet: str, destination_planet: str, speed: float, travel_time: float,
                    delta_v_required: float, fuel_mass_needed: float, transfer_type: str):
    """
    Saves detailed space travel calculation details to a JSON log file.

    The log file 'travel_log.json' will store a list of dictionaries,
    each containing a timestamp, origin body, destination body, spacecraft speed,
    travel time (traveler's frame), required delta-v, estimated fuel mass, and transfer type.

    If the log file does not exist, it will be created.
    If it exists but is empty or malformed, it will be initialized or overwritten
    in a way that new data can be appended.

    Args:
        source_planet (str): The name of the celestial body of origin (e.g., 'Earth').
        destination_planet (str): The name of the target celestial body (destination).
        speed (float): The calculated or estimated travel speed of the spacecraft (e.g., in Mm/s).
        travel_time (float): The calculated or estimated travel time from the traveler's perspective (e.g., in years).
        delta_v_required (float): The total change in velocity (delta-v) required for the mission (e.g., in km/s).
        fuel_mass_needed (float): The estimated mass of fuel required for the mission (e.g., in kg).
        transfer_type (str): The type of orbital transfer used (e.g., 'Hohmann Transfer', 'Ballistic Capture').
    """
    log_file = 'travel_log.json'
    log_data = []

    # Get current timestamp
    timestamp = datetime.datetime.now().isoformat()

    # Create the new entry with detailed orbital trajectory planning information
    new_entry = {
        'timestamp': timestamp,
        'origin_body': source_planet, # Renamed key for clarity as requested
        'destination_body': destination_planet, # Renamed key for clarity as requested
        'speed_Mm_s': speed,
        'travel_time_years': travel_time,
        'delta_v_required_km_s': delta_v_required,
        'fuel_mass_needed_kg': fuel_mass_needed,
        'transfer_type': transfer_type
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