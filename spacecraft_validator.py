def validate_spacecraft_name(name):
    """
    Validates a spacecraft name based on the following rules:
    - Must be between 3 and 20 characters long (inclusive).
    - Must consist only of alphanumeric characters (letters and numbers).

    Args:
        name (str): The spacecraft name to validate.

    Returns:
        bool: True if the name is valid, False otherwise.
    """
    if not isinstance(name, str):
        return False
    
    name_length = len(name)
    
    if not (3 <= name_length <= 20):
        return False
        
    if not name.isalnum():
        return False
        
    return True
