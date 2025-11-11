import json
import datetime
import os

def save_travel_log(destination: str, speed: float, travel_time: float):
    """
    Saves travel calculation details to a JSON log file.

    The log file 'travel_log.json' will store a list of dictionaries,
    each containing a timestamp, destination, speed, and travel time.

    If the log file does not exist, it will be created.
    If it exists but is empty or malformed, it will be initialized or overwritten
    in a way that new data can be appended.

    Args:
        destination (str): The name of the travel destination.
        speed (float): The calculated or estimated travel speed.
        travel_time (float): The calculated or estimated travel time.
    """
    log_file = 'travel_log.json'
    log_data = []

    # Get current timestamp
    timestamp = datetime.datetime.now().isoformat()

    # Create the new entry
    new_entry = {
        'timestamp': timestamp,
        'destination': destination,
        'speed': speed,
        'travel_time': travel_time
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
            # We can log a warning here if needed, but for simplicity,
            # we'll just start with an empty list for new data.
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
