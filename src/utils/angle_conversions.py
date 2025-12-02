import math
from typing import Union

def _validate_numeric_input(value: Union[int, float], name: str):
    """
    Internal helper to validate if a value is an integer or float.

    Args:
        value (Union[int, float]): The value to validate.
        name (str): The name of the value for error messages.

    Raises:
        ValueError: If the input value is not a numeric type.
    """
    if not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a numeric value (int or float), got {type(value).__name__}.")

def degrees_to_radians(degrees: Union[int, float]) -> float:
    """
    Converts an angle from degrees to radians.

    Args:
        degrees (Union[int, float]): The angle in degrees.

    Returns:
        float: The angle in radians.

    Raises:
        ValueError: If the input `degrees` is not a numeric type.
    """
    _validate_numeric_input(degrees, "Degrees")
    return math.radians(degrees)

def radians_to_degrees(radians: Union[int, float]) -> float:
    """
    Converts an angle from radians to degrees.

    Args:
        radians (Union[int, float]): The angle in radians.

    Returns:
        float: The angle in degrees.

    Raises:
        ValueError: If the input `radians` is not a numeric type.
    """
    _validate_numeric_input(radians, "Radians")
    return math.degrees(radians)

def dms_to_degrees(degrees: Union[int, float], minutes: Union[int, float], seconds: Union[int, float]) -> float:
    """
    Converts an angle from Degrees, Minutes, Seconds (DMS) format to decimal degrees.
    The sign of the `degrees` component determines the sign of the overall angle.
    Minutes and seconds are treated as non-negative magnitudes.

    Args:
        degrees (Union[int, float]): The degrees component (can be negative).
        minutes (Union[int, float]): The minutes component (0 <= minutes < 60).
        seconds (Union[int, float]): The seconds component (0 <= seconds < 60).

    Returns:
        float: The angle in decimal degrees.

    Raises:
        ValueError: If any input is not numeric, or if minutes/seconds are out of range.
    """
    _validate_numeric_input(degrees, "Degrees")
    _validate_numeric_input(minutes, "Minutes")
    _validate_numeric_input(seconds, "Seconds")

    if not (0 <= minutes < 60):
        raise ValueError(f"Minutes must be between 0 (inclusive) and 60 (exclusive), got {minutes}.")
    if not (0 <= seconds < 60):
        raise ValueError(f"Seconds must be between 0 (inclusive) and 60 (exclusive), got {seconds}.")

    sign = -1 if degrees < 0 else 1
    abs_degrees = abs(degrees)
    
    total_abs_degrees = abs_degrees + (minutes / 60.0) + (seconds / 3600.0)
    return sign * total_abs_degrees

def degrees_to_dms(degrees: Union[int, float]) -> tuple[int, int, float]:
    """
    Converts an angle from decimal degrees to Degrees, Minutes, Seconds (DMS) format.
    The sign is applied to the degrees component, and minutes/seconds are always positive.

    Args:
        degrees (Union[int, float]): The angle in decimal degrees.

    Returns:
        tuple[int, int, float]: A tuple containing (degrees, minutes, seconds).
                                Degrees is an int, minutes is an int, seconds is a float.

    Raises:
        ValueError: If the input `degrees` is not a numeric type.
    """
    _validate_numeric_input(degrees, "Degrees")

    is_negative = degrees < 0
    abs_degrees = abs(degrees)

    deg_abs = int(abs_degrees)
    fractional_degrees = abs_degrees - deg_abs

    minutes_float = fractional_degrees * 60
    min_val = int(minutes_float)
    
    seconds_float = (minutes_float - min_val) * 60

    # Handle floating point inaccuracies that might push seconds up to 60.0
    # Use a small tolerance for comparison
    if math.isclose(seconds_float, 60.0, rel_tol=1e-9, abs_tol=1e-12):
        seconds_float = 0.0
        min_val += 1
    
    # Handle floating point inaccuracies that might push minutes up to 60.0
    if math.isclose(min_val, 60.0, rel_tol=1e-9, abs_tol=1e-12):
        min_val = 0
        deg_abs += 1
    
    # Apply the original sign to the degree component
    if is_negative:
        final_deg = -deg_abs
    else:
        final_deg = deg_abs

    return (final_deg, min_val, seconds_float)

def dms_to_radians(degrees: Union[int, float], minutes: Union[int, float], seconds: Union[int, float]) -> float:
    """
    Converts an angle from Degrees, Minutes, Seconds (DMS) format to radians.

    Args:
        degrees (Union[int, float]): The degrees component.
        minutes (Union[int, float]): The minutes component.
        seconds (Union[int, float]): The seconds component.

    Returns:
        float: The angle in radians.

    Raises:
        ValueError: Propagates ValueError from `dms_to_degrees` if inputs are invalid.
    """
    decimal_degrees = dms_to_degrees(degrees, minutes, seconds)
    return degrees_to_radians(decimal_degrees)

def radians_to_dms(radians: Union[int, float]) -> tuple[int, int, float]:
    """
    Converts an angle from radians to Degrees, Minutes, Seconds (DMS) format.

    Args:
        radians (Union[int, float]): The angle in radians.

    Returns:
        tuple[int, int, float]: A tuple containing (degrees, minutes, seconds).
                                Degrees is an int, minutes is an int, seconds is a float.

    Raises:
        ValueError: Propagates ValueError from `radians_to_degrees` if input is not numeric.
    """
    decimal_degrees = radians_to_degrees(radians)
    return degrees_to_dms(decimal_degrees)

if __name__ == "__main__":
    print("--- Angle Conversion Self-Tests ---")

    # Test degrees_to_radians
    print(f"90 degrees to radians: {degrees_to_radians(90):.6f}") # Expected: 1.570796
    print(f"180 degrees to radians: {degrees_to_radians(180):.6f}") # Expected: 3.141593
    try:
        degrees_to_radians("invalid")
    except ValueError as e:
        print(f"Error (degrees_to_radians): {e}")

    print("\n--------------------\n")

    # Test radians_to_degrees
    print(f"math.pi/2 radians to degrees: {radians_to_degrees(math.pi / 2):.6f}") # Expected: 90.000000
    print(f"math.pi radians to degrees: {radians_to_degrees(math.pi):.6f}") # Expected: 180.000000
    try:
        radians_to_degrees([1, 2])
    except ValueError as e:
        print(f"Error (radians_to_degrees): {e}")

    print("\n--------------------\n")

    # Test dms_to_degrees
    print(f"30 deg 30 min 0 sec to degrees: {dms_to_degrees(30, 30, 0):.6f}") # Expected: 30.500000
    print(f"-10 deg 30 min 0 sec to degrees: {dms_to_degrees(-10, 30, 0):.6f}") # Expected: -10.500000
    print(f"0 deg 0 min 30 sec to degrees: {dms_to_degrees(0, 0, 30):.6f}") # Expected: 0.008333
    try:
        dms_to_degrees(10, 60, 0)
    except ValueError as e:
        print(f"Error (dms_to_degrees - invalid minutes): {e}")
    try:
        dms_to_degrees(10, 30, 60)
    except ValueError as e:
        print(f"Error (dms_to_degrees - invalid seconds): {e}")

    print("\n--------------------\n")

    # Test degrees_to_dms
    print(f"30.5 degrees to DMS: {degrees_to_dms(30.5)}") # Expected: (30, 30, 0.0)
    print(f"-10.5 degrees to DMS: {degrees_to_dms(-10.5)}") # Expected: (-10, 30, 0.0)
    print(f"0.008333333333333333 degrees to DMS: {degrees_to_dms(0.008333333333333333)}") # Expected: (0, 0, ~30.0)
    print(f"10.0 degrees to DMS: {degrees_to_dms(10.0)}") # Expected: (10, 0, 0.0)
    print(f"0.0 degrees to DMS: {degrees_to_dms(0.0)}") # Expected: (0, 0, 0.0)
    print(f"72.999 degrees to DMS: {degrees_to_dms(72.999)}") # Expected: (72, 59, 56.4)
    # Test for floating point inaccuracies near 60 seconds/minutes boundary
    # Value that is extremely close to 1 degree
    close_to_one_degree = 0.9999999999999999
    print(f"{close_to_one_degree} degrees to DMS: {degrees_to_dms(close_to_one_degree)}") # Expected: (0, 59, ~59.999)

    print("\n--------------------\n")

    # Test dms_to_radians
    deg_in, min_in, sec_in = 45, 15, 30
    dms_rad = dms_to_radians(deg_in, min_in, sec_in)
    print(f"{deg_in}d {min_in}m {sec_in}s to radians: {dms_rad:.6f}") # Expected: 0.795415
    print(f"Convert back to DMS: {radians_to_dms(dms_rad)}") # Expected: (45, 15, ~30.0)

    print("\n--------------------\n")

    # Test radians_to_dms
    rad_in = math.pi / 4
    rad_dms = radians_to_dms(rad_in)
    print(f"math.pi/4 radians to DMS: {rad_dms}") # Expected: (45, 0, 0.0)
    print(f"Convert back to radians: {dms_to_radians(*rad_dms):.6f}") # Expected: 0.785398

    print("\n--- All Tests Complete ---")