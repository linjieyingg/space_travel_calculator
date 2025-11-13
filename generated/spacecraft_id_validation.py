import re

def validate_spacecraft_id(spacecraft_id: str) -> bool:
    """
    Validates a spacecraft ID based on a set of criteria.

    A valid spacecraft ID must:
    1. Be a string.
    2. Not be empty.
    3. Be alphanumeric (contain only letters and numbers).
    4. Have a length between 5 and 20 characters, inclusive.

    Args:
        spacecraft_id (str): The ID string to validate.

    Returns:
        bool: True if the spacecraft ID is valid, False otherwise.
    """
    if not isinstance(spacecraft_id, str):
        return False
    if not spacecraft_id:
        return False
    if not (5 <= len(spacecraft_id) <= 20):
        return False
    # Check if the ID is alphanumeric (letters and numbers only)
    if not re.fullmatch(r'^[a-zA-Z0-9]*$', spacecraft_id):
        return False
    return True