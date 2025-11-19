import pytest
import math
from datetime import datetime, timedelta
from unittest.mock import patch, MagicMock
import main # Import the main module
import constants # Import the constants module

# Mock external dependencies for isolated testing of main.py
@pytest.fixture
def mock_dependencies():
    with patch('main.celestial_data.get_celestial_body_data') as mock_get_celestial, \
         patch('main.celestial_data.get_all_solar_system_destinations') as mock_get_all_dest, \
         patch('main.TrajectoryPlanner') as MockTrajectoryPlanner, \
         patch('main.propulsion_system.calculate_required_fuel_mass') as mock_calc_fuel_mass, \
         patch('main.fuel_calc.calculate_fuel_cost') as mock_calc_fuel_cost, \
         patch('main.travel_logger.save_travel_log') as mock_save_log, \
         patch('builtins.input') as mock_input, \
         patch('builtins.print') as mock_print:

        # Configure mock celestial data
        mock_get_celestial.side_effect = lambda name: {
            'earth': {'semi_major_axis_from_sun': 1.496e11, 'gravitational_parameter_mu': 3.986004418e14, 'radius': 6.371e6},
            'mars': {'semi_major_axis_from_sun': 2.279e11, 'gravitational_parameter_mu': 4.282837e13, 'radius': 3.389e6}
        }.get(name.lower())

        mock_get_all_dest.return_value = [{'name': 'Mars', 'distance': 7.83e10}]

        # Configure mock TrajectoryPlanner
        mock_planner_instance = MockTrajectoryPlanner.return_value
        mock_planner_instance.plan_hohmann_trajectory.return_value = {
            'success': True,
            'total_delta_v': 15000.0,
            'total_travel_time_seconds': 250 * 24 * 3600, # ~250 days in seconds
            'transfer_type': 'Hohmann Transfer'
        }

        # Configure mock fuel calculations
        mock_calc_fuel_mass.return_value = 50000.0
        mock_calc_fuel_cost.return_value = {'total_fuel_mass_needed': 50000.0, 'total_cost': 5000000.0}

        # Simulate user input
        mock_input.side_effect = [
            'Earth', # source_planet
            'Mars',  # destination_planet
            '100000', # spacecraft_dry_mass (kg)
            '450',   # engine_specific_impulse (s)
            '2024-01-01', # departure_date_str
            '30',    # user_age
            # main.py uses a constant for FUEL_PRICE_PER_UNIT, so it's not an interactive input.
        ]

        yield {
            'mock_get_celestial': mock_get_celestial,
            'mock_get_all_dest': mock_get_all_dest,
            'MockTrajectoryPlanner': MockTrajectoryPlanner,
            'mock_planner_instance': mock_planner_instance, # provide instance for direct checks
            'mock_calc_fuel_mass': mock_calc_fuel_mass,
            'mock_calc_fuel_cost': mock_calc_fuel_cost,
            'mock_save_log': mock_save_log,
            'mock_input': mock_input,
            'mock_print': mock_print
        }

def test_main_function_integrates_trajectory_planner_and_calculates_correctly(mock_dependencies):
    """
    Test the main function's orchestration logic, ensuring it correctly:
    1. Prompts for departure/arrival bodies and other inputs.
    2. Calls TrajectoryPlanner.plan_hohmann_trajectory with the correct arguments.
    3. Processes the planner's output for fuel and time calculations.
    4. Logs travel details and prints a comprehensive summary.
    """
    main.main()

    # Verify that input prompts were presented as expected
    input_calls = [call.args[0] for call in mock_dependencies['mock_input'].call_args_list]
    assert "Enter the source celestial body (e.g., Earth): " in input_calls
    assert "Enter the destination celestial body (e.g., Mars): " in input_calls
    assert "Enter spacecraft dry mass in kg (e.g., 100000): " in input_calls
    assert "Enter engine specific impulse in seconds (e.g., 450): " in input_calls
    assert "Enter departure date (YYYY-MM-DD): " in input_calls
    assert "Enter your current age in years: " in input_calls


    # Verify TrajectoryPlanner was instantiated and its method called correctly
    mock_dependencies['MockTrajectoryPlanner'].assert_called_once()
    mock_dependencies['mock_planner_instance'].plan_hohmann_trajectory.assert_called_once_with(
        departure_body_name='Earth', arrival_body_name='Mars'
    )

    # Verify fuel mass calculation uses delta-V from planner and dry mass input
    mock_dependencies['mock_calc_fuel_mass'].assert_called_once_with(
        delta_v=pytest.approx(15000.0), dry_mass=pytest.approx(100000.0), specific_impulse=pytest.approx(450.0)
    )

    # Verify fuel cost calculation uses fuel mass and the constant fuel price from main.py
    # NOTE: The value 100.0 is inferred from the `main.py` context in the repo summary,
    # as FUEL_PRICE_PER_UNIT is a constant in `main.py`.
    mock_dependencies['mock_calc_fuel_cost'].assert_called_once_with(
        total_fuel_mass_needed=pytest.approx(50000.0), fuel_price_per_unit=pytest.approx(100.0)
    )

    # Verify travel logger was called with correct data
    mock_dependencies['mock_save_log'].assert_called_once()
    save_log_args = mock_dependencies['mock_save_log'].call_args.kwargs
    assert save_log_args['source_planet'] == 'Earth'
    assert save_log_args['destination_planet'] == 'Mars'

    # The 'speed' for logging is derived from (abs(SMA_dest - SMA_src)) / total_travel_time_seconds
    # distance = abs(2.279e11 - 1.496e11) = 7.83e10 m
    # total_travel_time_seconds = 250 * 24 * 3600 = 21600000 seconds
    # Expected avg_speed_ms = 7.83e10 / 21600000 = 3625000 m/s
    assert save_log_args['speed'] == pytest.approx(3625000.0)

    # The 'travel_time' for logging is total_travel_time_seconds converted to days.
    assert save_log_args['travel_time'] == pytest.approx(250.0) # 250 days
    assert save_log_args['delta_v_required'] == pytest.approx(15000.0)
    assert save_log_args['fuel_mass_needed'] == pytest.approx(50000.0)
    assert save_log_args['transfer_type'] == 'Hohmann Transfer'

    # Collect all print calls (args[0] because print passes one string argument normally)
    printed_output = [call.args[0] for call in mock_dependencies['mock_print'].call_args_list if call.args]

    # Verify key information is printed using 'in' operator for substring matching
    assert any("Travel Summary" in s for s in printed_output)
    assert any("Source Planet: Earth" in s for s in printed_output)
    assert any("Destination Planet: Mars" in s for s in printed_output)
    assert any("Delta-V Required: 15000.0 m/s" in s for s in printed_output)
    assert any("Fuel Mass Needed: 50000.0 kg" in s for s in printed_output)
    assert any("Total Fuel Cost: 5000000.0" in s for s in printed_output)
    assert any("Transfer Type: Hohmann Transfer" in s for s in printed_output)
    # Check time outputs
    assert any("Travel Time (Traveler's Frame):" in s for s in printed_output)
    assert any("Travel Time (Earth's Frame):" in s for s in printed_output)
    assert any("Traveler's Age Upon Arrival:" in s for s in printed_output)
    assert any("Estimated Arrival Date:" in s for s in printed_output)


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
    # timedelta implicitly handles fractions by rounding to seconds.
    expected_date_half_year = departure_date + timedelta(days=0.5 * 365.25)
    assert arrival_date_half_year == expected_date_half_year
    assert arrival_date_half_year == datetime(2024, 7, 1, 3, 0) # 2024-01-01 + 182 days 6 hours

    # Small fraction of a year: 0.01 * 365.25 = 3.6525 days. Timedelta adds 3 days + 15.66 minutes.
    arrival_date_small_fraction = main.calc_arrival(departure_date, 0.01)
    expected_date_small_fraction = departure_date + timedelta(days=0.01 * 365.25)
    assert arrival_date_small_fraction == expected_date_small_fraction
    assert arrival_date_small_fraction == datetime(2024, 1, 4, 15, 39, 36) # 2024-01-01 + 3 days 15 hours 39 minutes 36 seconds

def test_calc_arrival_negative_years():
    """
    Test calc_arrival with negative years, expecting a past date.
    main.py returns a datetime object for negative years, and timedelta correctly
    calculates the past date.
    """
    departure_date = datetime(2023, 1, 1)
    # -5 * 365.25 = -1826.25 days. timedelta(days=-1826.25) is -1826 days and -6 hours.
    arrival_date = main.calc_arrival(departure_date, -5.0)
    expected_date = departure_date + timedelta(days=-5.0 * 365.25)
    assert arrival_date == expected_date
    assert arrival_date == datetime(2017, 12, 30, 18, 0) # 2023-01-01 minus 1826 days 6 hours

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