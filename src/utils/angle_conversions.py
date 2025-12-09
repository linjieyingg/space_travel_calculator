import math

def _is_number(value):
    """Helper to check if a value is a number."""
    return isinstance(value, (int, float))

def deg_to_rad(degrees: float) -> float:
    """
    Converts an angle from degrees to radians.

    Args:
        degrees (float): The angle in degrees.

    Returns:
        float: The angle in radians.

    Raises:
        ValueError: If the input is not a number.
    """
    if not _is_number(degrees):
        raise ValueError("Input 'degrees' must be a number.")
    return math.radians(degrees)

def rad_to_deg(radians: float) -> float:
    """
    Converts an angle from radians to degrees.

    Args:
        radians (float): The angle in radians.

    Returns:
        float: The angle in degrees.

    Raises:
        ValueError: If the input is not a number.
    """
    if not _is_number(radians):
        raise ValueError("Input 'radians' must be a number.")
    return math.degrees(radians)

def dms_to_deg(degrees: float, minutes: float, seconds: float) -> float:
    """
    Converts an angle from Degrees, Minutes, Seconds (DMS) format to decimal degrees.

    The sign of the angle is determined solely by the 'degrees' component.
    Minutes and seconds should always be non-negative in DMS representation.

    Args:
        degrees (float): The degrees component (can be negative).
        minutes (float): The minutes component (0 <= minutes < 60).
        seconds (float): The seconds component (0 <= seconds < 60).

    Returns:
        float: The angle in decimal degrees.

    Raises:
        ValueError: If any input is not a number, or if minutes/seconds are out of range.
    """
    if not (_is_number(degrees) and _is_number(minutes) and _is_number(seconds)):
        raise ValueError("All inputs (degrees, minutes, seconds) must be numbers.")
    if not (0 <= minutes < 60):
        raise ValueError("Minutes must be between 0 (inclusive) and 60 (exclusive).")
    if not (0 <= seconds < 60):
        raise ValueError("Seconds must be between 0 (inclusive) and 60 (exclusive).")

    sign = -1 if degrees < 0 else 1
    abs_degrees = abs(degrees)

    decimal_degrees = abs_degrees + (minutes / 60.0) + (seconds / 3600.0)
    return sign * decimal_degrees

def deg_to_dms(decimal_degrees: float) -> tuple[int, int, float]:
    """
    Converts an angle from decimal degrees to Degrees, Minutes, Seconds (DMS) format.

    Args:
        decimal_degrees (float): The angle in decimal degrees.

    Returns:
        tuple[int, int, float]: A tuple containing (degrees, minutes, seconds).
                                Degrees will have the original sign, minutes and seconds are non-negative.

    Raises:
        ValueError: If the input is not a number.
    """
    if not _is_number(decimal_degrees):
        raise ValueError("Input 'decimal_degrees' must be a number.")

    sign = -1 if decimal_degrees < 0 else 1
    abs_decimal_degrees = abs(decimal_degrees)

    degrees_part = int(abs_decimal_degrees)
    remaining_degrees = abs_decimal_degrees - degrees_part

    minutes_float = remaining_degrees * 60
    minutes_part = int(minutes_float)
    remaining_minutes = minutes_float - minutes_part

    seconds_part = remaining_minutes * 60

    return (sign * degrees_part, minutes_part, seconds_part)

def rad_to_dms(radians: float) -> tuple[int, int, float]:
    """
    Converts an angle from radians to Degrees, Minutes, Seconds (DMS) format.

    Args:
        radians (float): The angle in radians.

    Returns:
        tuple[int, int, float]: A tuple containing (degrees, minutes, seconds).
                                Degrees will have the original sign, minutes and seconds are non-negative.

    Raises:
        ValueError: If the input is not a number.
    """
    if not _is_number(radians):
        raise ValueError("Input 'radians' must be a number.")
    
    decimal_degrees = rad_to_deg(radians)
    return deg_to_dms(decimal_degrees)

def dms_to_rad(degrees: float, minutes: float, seconds: float) -> float:
    """
    Converts an angle from Degrees, Minutes, Seconds (DMS) format to radians.

    Args:
        degrees (float): The degrees component.
        minutes (float): The minutes component (0 <= minutes < 60).
        seconds (float): The seconds component (0 <= seconds < 60).

    Returns:
        float: The angle in radians.

    Raises:
        ValueError: If any input is not a number, or if minutes/seconds are out of range.
    """
    # The validation for minutes/seconds and numerical types is handled by dms_to_deg
    decimal_degrees = dms_to_deg(degrees, minutes, seconds)
    return deg_to_rad(decimal_degrees)