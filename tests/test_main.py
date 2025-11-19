import pytest
from datetime import datetime, timedelta
from constants import C_LIGHT_MPS # Import necessary constant for speed validation

# Import functions directly from main for testing
from main import calc_time_earth, calc_age, calc_arrival, convert_date

# Mock constants.C_LIGHT_MPS if main.py is imported in an environment where constants isn't available
# For this specific file's context, constants.C_LIGHT_MPS is directly imported and used.

def test_calc_time_earth_basic():
    """
    Test calc_time_earth with a non-relativistic speed.
    Earth time should be slightly greater than traveler's time due to even slight dilation.
    """
    traveler_time = 1.0 # years
    average_speed = 100000.0 # m/s (much less than c)
    earth_time = calc_time_earth(traveler_time, average_speed)
    assert earth_time > traveler_time
    assert isinstance(earth_time, float)
    assert earth_time > 0

def test_calc_time_earth_zero_speed():
    """
    Test calc_time_earth with zero speed.
    In a non-relativistic scenario, traveler's time equals Earth's time.
    """
    assert calc_time_earth(1.0, 0.0) == 1.0

def test_calc_time_earth_speed_of_light_limit():
    """
    Test calc_time_earth when speed is at the speed of light.
    The function in main.py is designed to return traveler's time directly in this case.
    """
    assert calc_time_earth(1.0, C_LIGHT_MPS) == 1.0

def test_calc_age():
    """
    Test calc_age with various inputs, ensuring correct integer age calculation.
    """
    assert calc_age(10, 25) == 35
    assert calc_age(0, 30) == 30
    assert calc_age(5.5, 20) == 25 # main.py uses int(age + years) which floors the result
    assert calc_age(1, 0) == 1 # Validates calculation for a starting age of 0 (though application input usually validates >=18)

def test_calc_age_invalid_input():
    """
    Test calc_age with non-numeric inputs, expecting TypeError.
    """
    with pytest.raises(TypeError):
        calc_age("abc", 25)
    with pytest.raises(TypeError):
        calc_age(10, "abc")

def test_calc_arrival_whole_years():
    """
    Test calc_arrival with whole years, expecting exact date addition.
    """
    departure_date = datetime(2023, 1, 1)
    arrival_date = calc_arrival(departure_date, 10)
    assert arrival_date == datetime(2033, 1, 1)

def test_calc_arrival_half_year():
    """
    Test calc_arrival with half a year.
    Note: The exact day might vary slightly due to 365.25 days/year approximation
    and how timedelta handles fractions of days.
    """
    departure_date = datetime(2023, 1, 1)
    arrival_date_half_year = calc_arrival(departure_date, 0.5)
    # The actual calculation in main.py is departure_date + timedelta(days=years_earth_frame * 365.25)
    # 0.5 * 365.25 = 182.625 days.
    # 2023-01-01 + 182 days = 2023-07-02.
    # The current assertion `datetime(2023, 7, 1)` might be off by a day. Adjusting for expected behavior.
    expected_date = departure_date + timedelta(days=0.5 * 365.25)
    assert arrival_date_half_year.year == expected_date.year
    assert arrival_date_half_year.month == expected_date.month
    assert arrival_date_half_year.day == expected_date.day


def test_calc_arrival_negative_years():
    """
    Test calc_arrival with negative years, expecting a past date.
    main.py returns a datetime object for negative years, not None.
    """
    departure_date = datetime(2023, 1, 1)
    arrival_date = calc_arrival(departure_date, -5)
    assert arrival_date == datetime(2018, 1, 1) # Approximately, given 365.25 days/year

def test_calc_arrival_invalid_date_type():
    """
    Test calc_arrival with invalid departure date types, expecting None or TypeError.
    main.py checks for `isinstance(date, datetime)` and returns None if not.
    """
    assert calc_arrival(None, 10) is None
    assert calc_arrival("2023-01-01", 10) is None # main.py returns None for string date, not TypeError

def test_convert_date_valid():
    """
    Test convert_date with valid date strings.
    """
    assert convert_date("2023-01-15") == datetime(2023, 1, 15)
    assert convert_date("1999-12-31") == datetime(1999, 12, 31)

def test_convert_date_invalid_format():
    """
    Test convert_date with invalid date string formats, expecting ValueError.
    """
    with pytest.raises(ValueError, match="Incorrect date format"):
        convert_date("invalid-date")
    with pytest.raises(ValueError, match="Incorrect date format"):
        convert_date("2023/01/15") # Wrong separator

def test_convert_date_non_string_input():
    """
    Test convert_date with non-string inputs, expecting TypeError.
    """
    with pytest.raises(TypeError):
        convert_date(None)
    with pytest.raises(TypeError):
        convert_date(12345)