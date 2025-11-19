import pytest
import math
from datetime import datetime, timedelta
from unittest.mock import patch
import main # Import the main module
import constants # Import the constants module

# Mocking external functions for isolated testing of main's components
@patch('main.trajectory_planner.TrajectoryPlanner')
@patch('main.fuel_calc.calculate_fuel_cost')
@patch('main.propulsion_system.calculate_required_fuel_mass')
@patch('main.travel_logger.save_travel_log')
def test_main_function_smoke_test(mock_save_log, mock_calculate_fuel_mass, mock_calculate_fuel_cost, MockTrajectoryPlanner):
    """
    Smoke test for the main() function to ensure it runs without crashing
    and orchestrates calls to its dependencies.
    """
    # Mock user inputs for main function
    mock_inputs = [
        'Earth', # source_planet
        'Mars',  # destination_planet
        '10000', # spacecraft_dry_mass (kg)
        '450',   # engine_specific_impulse (s)
        '2024-01-01', # departure_date_str
        '30',    # user_age
    ]

    with patch('builtins.input', side_effect=mock_inputs):
        # Mock the TrajectoryPlanner instance and its plan_hohmann_trajectory method
        mock_planner_instance = MockTrajectoryPlanner.return_value
        mock_planner_instance.plan_hohmann_trajectory.return_value = {
            'success': True,
            'total_delta_v': 10000.0, # Example delta-v
            'total_travel_time_days': 250.0, # Example travel time
            'transfer_type': 'Hohmann Transfer',
            'error': None
        }

        mock_calculate_fuel_mass.return_value = 50000.0 # Example fuel mass
        # Assuming FUEL_PRICE_PER_UNIT is 100.0 as seen in main.py context
        mock_calculate_fuel_cost.return_value = {'total_cost': 5000000.0} # Example fuel cost

        # Call the main function
        main.main()

        # Assert that key functions were called
        mock_planner_instance.plan_hohmann_trajectory.assert_called_once_with('Earth', 'Mars')
        mock_calculate_fuel_mass.assert_called_once_with(
            delta_v=pytest.approx(10000.0),
            dry_mass=pytest.approx(10000.0),
            specific_impulse=pytest.approx(450.0)
        )
        mock_calculate_fuel_cost.assert_called_once_with(
            total_fuel_mass_needed=pytest.approx(50000.0),
            fuel_price_per_unit=pytest.approx(100.0) # From main.py's FUEL_PRICE_PER_UNIT
        )
        mock_save_log.assert_called_once()
        # Further assertions could check the exact arguments passed to mock_save_log if needed

# --- Start of unit tests for main's helper functions ---

def test_calc_time_earth_no_relativity():
    """
    Test calc_time_earth with a very low speed where relativistic effects are negligible.
    Expected Earth time should be approximately equal to traveler's time.
    """
    years_traveler = 1.0
    average_speed_ms = 1.0  # A very low speed to ensure gamma is ~1
    expected_years_earth = 1.0
    calculated_years_earth = main.calc_time_earth(years_traveler, average_speed_ms)
    assert calculated_years_earth == pytest.approx(expected_years_earth, abs=1e-9)

def test_calc_time_earth_high_speed():
    """
    Test calc_time_earth with a speed close to c to observe significant relativistic time dilation.
    """
    years_traveler = 1.0
    average_speed_ms = 0.8 * constants.C_LIGHT_MPS # 80% speed of light
    # Lorentz factor gamma = 1 / sqrt(1 - (v/c)^2)
    # Earth time = traveler time * gamma
    expected_years_earth = years_traveler / math.sqrt(1 - (0.8**2))
    calculated_years_earth = main.calc_time_earth(years_traveler, average_speed_ms)
    assert calculated_years_earth == pytest.approx(expected_years_earth, rel=1e-9)

def test_calc_time_earth_speed_of_light_limit():
    """
    Test calc_time_earth when speed is exactly the speed of light.
    The function in main.py is designed to return traveler's time directly in this case
    to avoid division by zero or complex numbers due to floating point inaccuracies
    at v=c, as time dilation would theoretically be infinite, but in practice,
    v cannot equal c for massive objects.
    """
    assert main.calc_time_earth(1.0, constants.C_LIGHT_MPS) == pytest.approx(1.0)

def test_calc_time_earth_zero_speed():
    """
    Test calc_time_earth with zero speed.
    In this scenario, traveler's time should equal Earth's time (no relative motion).
    """
    assert main.calc_time_earth(1.0, 0.0) == pytest.approx(1.0)

def test_calc_age_whole_years():
    """
    Test calc_age with whole number inputs for years.
    """
    assert main.calc_age(10, 25) == 35
    assert main.calc_age(0, 30) == 30
    assert main.calc_age(1, 0) == 1

def test_calc_age_fractional_years():
    """
    Test calc_age with fractional years, ensuring the result is floored (integer).
    """
    assert main.calc_age(5.5, 20) == 25 # int(20 + 5.5) = int(25.5) = 25
    assert main.calc_age(0.5, 20) == 20 # int(20 + 0.5) = 20
    assert main.calc_age(0.9, 20) == 20 # int(20 + 0.9) = 20
    assert main.calc_age(1.0, 20) == 21 # int(20 + 1.0) = 21

def test_calc_age_invalid_input_types():
    """
    Test calc_age with non-numeric inputs, expecting TypeError.
    """
    with pytest.raises(TypeError):
        main.calc_age("abc", 25)
    with pytest.raises(TypeError):
        main.calc_age(10, "abc")
    with pytest.raises(TypeError):
        main.calc_age(None, 25)
    with pytest.raises(TypeError):
        main.calc_age(10, None)

def test_calc_arrival_whole_years():
    """
    Test calc_arrival with whole years, expecting exact date addition.
    """
    departure_date = datetime(2024, 1, 1) # Starting on a leap year to test robustness
    arrival_date = main.calc_arrival(departure_date, 10.0)
    assert arrival_date == datetime(2034, 1, 1) # 10 years later

def test_calc_arrival_fractional_years():
    """
    Test calc_arrival with fractional years, comparing against exact timedelta calculation.
    The function uses 365.25 days per year, and timedelta truncates fractional days.
    """
    departure_date = datetime(2024, 1, 1) # Start on a leap year
    
    # 0.5 year travel time: 0.5 * 365.25 = 182.625 days. Timedelta adds 182 days.
    arrival_date_half_year = main.calc_arrival(departure_date, 0.5)
    expected_date_half_year = departure_date + timedelta(days=math.floor(0.5 * 365.25))
    assert arrival_date_half_year == expected_date_half_year
    assert arrival_date_half_year == datetime(2024, 7, 1) # 2024-01-01 + 182 days = 2024-07-01

    # Small fraction of a year: 0.01 * 365.25 = 3.6525 days. Timedelta adds 3 days.
    arrival_date_small_fraction = main.calc_arrival(departure_date, 0.01)
    expected_date_small_fraction = departure_date + timedelta(days=math.floor(0.01 * 365.25))
    assert arrival_date_small_fraction == expected_date_small_fraction
    assert arrival_date_small_fraction == datetime(2024, 1, 4)

def test_calc_arrival_negative_years():
    """
    Test calc_arrival with negative years, expecting a past date.
    main.py returns a datetime object for negative years, and timedelta correctly
    calculates the past date.
    """
    departure_date = datetime(2023, 1, 1)
    # -5 * 365.25 = -1826.25 days. Timedelta correctly applies negative days.
    arrival_date = main.calc_arrival(departure_date, -5.0)
    expected_date = departure_date + timedelta(days=-5 * 365.25)
    assert arrival_date == expected_date
    # Specifically assert the resulting date and time components
    assert arrival_date == datetime(2017, 12, 31, 15, 0) # 2023-01-01 minus 1826 days 6 hours

def test_calc_arrival_invalid_date_type():
    """
    Test calc_arrival with invalid departure date types, expecting None as per main.py's
    internal validation (`isinstance(date, datetime)`).
    """
    assert main.calc_arrival(None, 10) is None
    assert main.calc_arrival("2023-01-01", 10) is None # main.py returns None for string date, not TypeError

def test_convert_date_valid():
    """
    Test convert_date with valid date strings in 'YYYY-MM-DD' format.
    """
    assert main.convert_date("2023-01-15") == datetime(2023, 1, 15)
    assert main.convert_date("1999-12-31") == datetime(1999, 12, 31)
    assert main.convert_date("2024-02-29") == datetime(2024, 2, 29) # Leap year

def test_convert_date_invalid_format():
    """
    Test convert_date with invalid date string formats, expecting ValueError.
    """
    with pytest.raises(ValueError, match="Incorrect date format"):
        main.convert_date("invalid-date")
    with pytest.raises(ValueError, match="Incorrect date format"):
        main.convert_date("2023/01/15") # Wrong separator
    with pytest.raises(ValueError, match="Incorrect date format"):
        main.convert_date("2023-1-1") # Month/day not zero-padded

def test_convert_date_non_string_input():
    """
    Test convert_date with non-string inputs, expecting TypeError.
    """
    with pytest.raises(TypeError):
        main.convert_date(None)
    with pytest.raises(TypeError):
        main.convert_date(12345)
    with pytest.raises(TypeError):
        main.convert_date(['2023-01-01'])