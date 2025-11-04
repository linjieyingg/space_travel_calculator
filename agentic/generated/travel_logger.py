import json
from datetime import datetime
import os

def save_travel_log(destination: str, speed: float, travel_time: float, log_file: str = 'travel_log.json'):
    """
    Saves travel calculations to a JSON log file.

    Args:
        destination (str): The destination of the travel.
        speed (float): The average speed during travel (e.g., in km/h or mph).
        travel_time (float): The total travel time (e.g., in hours).
        log_file (str, optional): The name of the JSON log file. Defaults to 'travel_log.json'.
    """
    
    timestamp = datetime.now().isoformat()
    
    new_entry = {
        "timestamp": timestamp,
        "destination": destination,
        "speed": speed,
        "travel_time": travel_time
    }
    
    log_data = []
    
    # Load existing data if the file exists
    if os.path.exists(log_file):
        try:
            with open(log_file, 'r') as f:
                log_data = json.load(f)
        except json.JSONDecodeError:
            # Handle case where file is empty or malformed JSON
            print(f"Warning: {log_file} is empty or contains invalid JSON. Starting a new log.")
            log_data = []
    
    log_data.append(new_entry)
    
    # Save all data back to the file
    with open(log_file, 'w') as f:
        json.dump(log_data, f, indent=4)