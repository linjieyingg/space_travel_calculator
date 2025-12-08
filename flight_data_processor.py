import math
from datetime import datetime, timedelta
from typing import Union

import constants

"""
This module provides functions for post-trajectory calculations related to
relativistic time dilation, estimated arrival age, and estimated arrival date.
It centralizes logic previously found in main.py, improving modularity and testability.
"""


def calc_time_earth(years_traveler_frame: Union[int, float], average_speed_ms: Union[int, float]) -> float:
    """
    Calculates the travel time in Earth's reference frame, accounting for
    relativistic time dilation.

    Args:
        years_traveler_frame (Union[int, float]): The travel time experienced by the traveler, in years.
        average_speed_ms (Union[int, float]): The average speed of the spacecraft in meters per second.

    Returns:
        float: The equivalent travel time in Earth's reference frame, in years.
               Returns float('inf') if average_speed_ms is greater than or equal to the speed of light.

    Raises:
        ValueError: If `years_traveler_frame` is negative, `average_speed_ms` is negative,
                    or if inputs are not numeric.
    """
    if not isinstance(years_traveler_frame, (int, float)) or not isinstance(average_speed_ms, (int, float)):
        raise ValueError("Travel time in traveler's frame and average speed must be numeric.")
    if years_traveler_frame < 0:
        raise ValueError("Travel time in traveler's frame cannot be negative.")
    if average_speed_ms < 0:
        raise ValueError("Average speed cannot be negative.")

    if average_speed_ms >= constants.C_LIGHT_MPS:
        # As per relativistic physics, if speed is at or above 'c', time in Earth's
        # reference frame becomes infinite or undefined.
        return float('inf')
    
    # If speed is zero, there is no time dilation. Avoid potential division by zero
    # or floating point issues when speed_ratio_squared is extremely small.
    if math.isclose(average_speed_ms, 0.0, abs_tol=1e-9):
        return float(years_traveler_frame)

    speed_ratio_squared = (average_speed_ms / constants.C_LIGHT_MPS) ** 2
    
    # Calculate the denominator for the Lorentz factor.
    # It should be > 0 because average_speed_ms is strictly less than C_LIGHT_MPS.
    denominator = math.sqrt(1 - speed_ratio_squared)
    
    # A safety check for floating point inaccuracies that might make denominator very close to zero.
    if math.isclose(denominator, 0.0, abs_tol=1e-9):
        return float('inf') # Treat as infinite due to extreme relativistic effects

    years_earth_frame = years_traveler_frame / denominator
    return years_earth_frame


def calc_age(years_travel_time: Union[int, float], starting_age: Union[int, float]) -> float:
    """
    Calculates the traveler's age upon arrival, given their starting age and
    the total travel time.

    Args:
        years_travel_time (Union[int, float]): The total travel time in years (e.g., in Earth's reference frame).
        starting_age (Union[int, float]): The traveler's age at the beginning of the journey.

    Returns:
        float: The traveler's estimated age upon arrival.

    Raises:
        ValueError: If `years_travel_time` or `starting_age` is negative, or if inputs are not numeric.
    """
    if not isinstance(years_travel_time, (int, float)) or not isinstance(starting_age, (int, float)):
        raise ValueError("Travel time and starting age must be numeric.")
    if years_travel_time < 0:
        raise ValueError("Travel time cannot be negative.")
    if starting_age < 0:
        raise ValueError("Starting age cannot be negative.")

    arrival_age = starting_age + years_travel_time
    return float(arrival_age) # Return as float for consistency with years_travel_time


def calc_arrival(departure_date: datetime, years_earth_frame: Union[int, float]) -> datetime:
    """
    Calculates the estimated arrival date based on the departure date and
    the total travel time in Earth's reference frame.

    Args:
        departure_date (datetime): The exact date and time of departure.
        years_earth_frame (Union[int, float]): The total travel time in Earth's reference frame, in years.

    Returns:
        datetime: The estimated date and time of arrival.

    Raises:
        ValueError: If `years_earth_frame` is negative or not numeric, or if `departure_date`
                    is not a datetime object.
    """
    if not isinstance(departure_date, datetime):
        raise ValueError("Departure date must be a datetime object.")
    if not isinstance(years_earth_frame, (int, float)):
        raise ValueError("Years in Earth frame must be numeric.")
    if years_earth_frame < 0:
        raise ValueError("Travel time in Earth's frame cannot be negative.")
    if years_earth_frame == float('inf'):
        return datetime.max # Represents an infinitely far future date

    # Using 365.25 days per year to account for leap years on average.
    total_days = years_earth_frame * 365.25
    
    # Guard against excessively large timedelta, which can raise OverflowError
    # datetime.timedelta constructor takes an integer number of days.
    # Max days for timedelta is limited, approximately 999999999 days.
    # 999999999 / 365.25 years ~ 2.7 million years.
    # While float('inf') is handled, large but finite years_earth_frame could still exceed timedelta capacity.
    if total_days > timedelta.max.days:
        return datetime.max
        
    arrival_date = departure_date + timedelta(days=total_days)
    return arrival_date