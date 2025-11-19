import pytest
import math
from datetime import datetime, timedelta
from unittest.mock import patch, MagicMock, call
import main # Import the main module
import constants # Import the constants module
import celestial_data # Required for some mocks

# Mock external dependencies for isolated testing of main.py
@pytest.fixture
def mock_dependencies():
    # Use patch context managers for cleaner setup and teardown
    with patch('main.celestial_data.get_celestial_body_data') as mock_get_celestial, \
         patch('main.celestial_data.get_all_solar_system_destinations') as mock_get_all_dest, \
         patch('main.TrajectoryPlanner') as MockTrajectoryPlanner, \
         patch('main.fuel_optimizer.optimize_fuel_for_trajectory') as mock_optimize_fuel_for_trajectory, \
         patch('main.fuel_calc.calculate_fuel_cost') as mock_calc_fuel_cost, \
         patch('main.travel_logger.save_travel_log') as mock_save_log, \
         patch('builtins.input') as mock_input, \
         patch('builtins.print') as mock_print:

        # Configure mock celestial data
        mock_get_celestial.side_effect = lambda name: {
            'sun': {'name': 'Sun', 'gravitational_parameter_mu': 1.327e20},
            'earth': {'name': 'Earth', 'semi_major_axis_from_sun': 1.496e11, 'gravitational_parameter_mu': 3.986004418e14, 'radius': 6.371e6},
            'mars': {'name': 'Mars', 'semi_major_axis_from_sun': 2.279e11, 'gravitational_parameter_mu': 4.282837e13, 'radius': 3.389e6}
        }.get(name.lower())

        mock_get_all_dest.return_value = [{'name': 'Mars', 'distance': 7.83e10}]

        # Configure mock TrajectoryPlanner
        mock_planner_instance = MockTrajectoryPlanner.return_value
        # This output structure must match what `main.py` expects from TrajectoryPlanner
        mock_planner_instance.plan_hohmann_trajectory.return_value = {
            'success': True,
            'departure_body': 'Earth',
            'arrival_body': 'Mars',
            'transfer_type': 'Hohmann Transfer',
            'total_delta_v_mps': 15000.0,
            'travel_time_days': 250.0,
            'travel_time_h_m_s': '250 days, 00 hours, 00 minutes, 00 seconds',
            'total_travel_time_seconds': 250 * 24 * 3600, # Added for main.py to use directly
            # No 'average_hohmann_speed_ms' in actual trajectory_planner output, so main.py will use fallback
        }

        # Configure mock fuel calculations
        mock_optimize_fuel_for_trajectory.return_value = 50000.0 # kg
        mock_calc_fuel_cost.return_value = {'total_fuel_mass_needed': 50000.0, 'total_cost': 5000000.0} # $

        # Simulate user input
        # Ensure this matches the sequence of input() calls in main.py exactly
        mock_input.side_effect = [
            'Earth',        # source_planet_input
            'Mars',         # destination_planet_input
            '100000',       # spacecraft_dry_mass (kg)
            '450',          # engine_specific_impulse (s)
            '2024-01-01',   # departure_date_str
            '30',           # user_age
            '100'           # fuel_price_per_unit (cost per kg)
        ]

        # Use spec=True for better mocking behavior, ensuring mock has methods of original
        MockTrajectoryPlanner.return_value.plan_hohmann_trajectory.spec = True

        yield {
            'mock_get_celestial': mock_get_celestial,
            'mock_get_all_dest': mock_get_all_dest,
            'MockTrajectoryPlanner': MockTrajectoryPlanner,
            'mock_planner_instance': mock_planner_instance,
            'mock_optimize_fuel_for_trajectory': mock_optimize_fuel_for_trajectory,
            'mock_calc_fuel_cost': mock_calc_fuel_cost,
            'mock_save_log': mock_save_log,
            'mock_input': mock_input,
            'mock_print': mock_print
        }

def test_main_function_integrates_all_components_correctly(mock_dependencies):
    """
    Tests the main function's end-to-end integration and orchestration logic,
    ensuring correct calls to mocked dependencies and proper output.
    """
    main.main()

    # Verify input prompts were presented (checking actual prompts, not just "Enter...")
    mock_dependencies['mock_input'].assert_has_calls([
        call("Enter your Departure Body (e.g., Earth): "),
        call("Enter your Arrival Body (e.g., Mars): "),
        call("Enter spacecraft dry mass (kg): "),
        call("Enter engine specific impulse (s): "),
        call("Enter departure date (YYYY-MM-DD): "),
        call("Enter your current age (years): "),
        call("Enter fuel price per unit mass (e.g., cost per kg): "),
    ])

    # Verify TrajectoryPlanner was instantiated and its method called correctly
    mock_dependencies['MockTrajectoryPlanner'].assert_called_once()
    mock_dependencies['mock_planner_instance'].plan_hohmann_trajectory.assert_called_once_with(
        departure_body_name='Earth', arrival_body_name='Mars'
    )

    # Verify fuel optimizer was called with correct arguments from main's parsed inputs and planner result
    mock_dependencies['mock_optimize_fuel_for_trajectory'].assert_called_once_with(
        start_body='Earth',
        end_body='Mars',
        trajectory_type='Hohmann Transfer',
        spacecraft_dry_mass=pytest.approx(100000.0),
        engine_specific_impulse=pytest.approx(450.0)
    )

    # Verify fuel cost calculation
    mock_dependencies['mock_calc_fuel_cost'].assert_called_once_with(
        total_fuel_mass_needed=pytest.approx(50000.0), fuel_price_per_unit=pytest.approx(100.0)
    )

    # Verify travel logger was called with correct data
    mock_dependencies['mock_save_log'].assert_called_once()
    save_log_args = mock_dependencies['mock_save_log'].call_args.kwargs
    assert save_log_args['source_planet'] == 'Earth'
    assert save_log_args['destination_planet'] == 'Mars'
    assert save_log_args['delta_v_required'] == pytest.approx(15000.0)
    assert save_log_args['fuel_mass_needed'] == pytest.approx(50000.0)
    assert save_log_args['transfer_type'] == 'Hohmann Transfer'
    
    # Check average_speed_for_relativistic_calc which falls back to 0.001 * C_LIGHT_MPS
    assert save_log_args['speed'] == pytest.approx(0.001 * constants.C_LIGHT_MPS)
    # Check travel_time for logger, which is total_travel_time_seconds from planner mock
    assert save_log_args['travel_time'] == pytest.approx(250 * 24 * 3600)

    # Collect all print calls (args[0] because print passes one string argument normally)
    # Filter out empty calls or non-string calls
    printed_output = [call.args[0] for call in mock_dependencies['mock_print'].call_args_list if call.args and isinstance(call.args[0], str)]

    # Verify key summary information is printed
    assert any("Travel Summary" in s for s in printed_output)
    assert any("Departure Body: Earth" in s for s in printed_output)
    assert any("Arrival Body: Mars" in s for s in printed_output)
    assert any("Transfer Type: Hohmann Transfer" in s for s in printed_output)
    assert any("Delta-V Required: 15000.00 m/s" in s for s in printed_output)
    assert any("Estimated Fuel Mass Needed: 50000.00 kg" in s for s in printed_output)
    assert any("Total Fuel Cost: $5,000,000.00" in s for s in printed_output)
    
    # Verify travel times and age, accounting for the fallback relativistic speed
    # Traveler's time is 250 days = (250 / 365.25) years
    travel_time_traveler_years = (250 * 24 * 3600) / main.YEAR_TO_SECONDS
    assert any(f"Travel Time (Traveler's Frame): {travel_time_traveler_years:.2f} years" in s for s in printed_output)
    
    # Earth time based on fallback speed
    lorentz_factor = 1 / math.sqrt(1 - (0.001 * constants.C_LIGHT_MPS / constants.C_LIGHT_MPS)**2)
    expected_earth_time_years = travel_time_traveler_years * lorentz_factor
    assert any(f"Travel Time (Earth's Frame): {expected_earth_time_years:.2f} years" in s for s in printed_output)

    # Arrival age
    # main.calc_age returns float, so exact float comparison for sum, then .0f for print assertion
    assert any(f"Your Age Upon Arrival: {30 + travel_time_traveler_years:.0f} years" in s for s in printed_output)
    
    # Arrival date
    departure_date = datetime(2024, 1, 1)
    expected_arrival_date = departure_date + timedelta(days=expected_earth_time_years * 365.25)
    assert any(f"Estimated Arrival Date: {expected_arrival_date.strftime('%Y-%m-%d')}" in s for s in printed_output)

    assert any("Travel details logged successfully." in s for s in printed_output)


def test_main_function_trajectory_planner_fails(mock_dependencies):
    """
    Test main function's behavior when TrajectoryPlanner returns a failure.
    """
    mock_dependencies['mock_planner_instance'].plan_hohmann_trajectory.return_value = {
        'success': False,
        'error': 'Cannot plan trajectory between these bodies.'
    }
    mock_dependencies['mock_input'].side_effect = [
        'Earth', 'Pluto', '100000', '450', '2024-01-01', '30', '100'
    ]
    
    main.main()

    printed_output = [call.args[0] for call in mock_dependencies['mock_print'].call_args_list if call.args and isinstance(call.args[0], str)]
    assert any("Error planning trajectory: Cannot plan trajectory between these bodies." in s for s in printed_output)
    # Ensure subsequent calculations (fuel, log) are NOT called
    mock_dependencies['mock_optimize_fuel_for_trajectory'].assert_not_called()
    mock_dependencies['mock_calc_fuel_cost'].assert_not_called()
    mock_dependencies['mock_save_log'].assert_not_called()

def test_main_function_fuel_optimizer_raises_value_error(mock_dependencies):
    """
    Test main function's behavior when fuel_optimizer.optimize_fuel_for_trajectory raises ValueError.
    """
    mock_dependencies['mock_optimize_fuel_for_trajectory'].side_effect = ValueError("Invalid fuel parameters.")
    mock_dependencies['mock_input'].side_effect = [
        'Earth', 'Mars', '100000', '450', '2024-01-01', '30', '100'
    ]

    main.main()

    printed_output = [call.args[0] for call in mock_dependencies['mock_print'].call_args_list if call.args and isinstance(call.args[0], str)]
    assert any("Error optimizing fuel mass: Invalid fuel parameters." in s for s in printed_output)
    # Ensure fuel cost and logging are NOT called
    mock_dependencies['mock_calc_fuel_cost'].assert_not_called()
    mock_dependencies['mock_save_log'].assert_not_called()

def test_main_function_fuel_cost_calculation_fails(mock_dependencies):
    """
    Test main function's behavior when fuel_calc.calculate_fuel_cost returns None or incomplete data.
    """
    mock_dependencies['mock_calc_fuel_cost'].return_value = None # Simulate failure
    mock_dependencies['mock_input'].side_effect = [
        'Earth', 'Mars', '100000', '450', '2024-01-01', '30', '100'
    ]

    main.main()

    printed_output = [call.args[0] for call in mock_dependencies['mock_print'].call_args_list if call.args and isinstance(call.args[0], str)]
    assert any("Error calculating fuel cost. Check input parameters." in s for s in printed_output)
    # Ensure logging is NOT called
    mock_dependencies['mock_save_log'].assert_not_called()

def test_main_function_get_all_solar_system_destinations_fails(mock_dependencies):
    """
    Test main function's behavior when celestial_data.get_all_solar_system_destinations returns an empty list.
    """
    mock_dependencies['mock_get_all_dest'].return_value = []
    
    main.main()

    printed_output = [call.args[0] for call in mock_dependencies['mock_print'].call_args_list if call.args and isinstance(call.args[0], str)]
    assert any("Error: Could not retrieve solar system destinations. Exiting." in s for s in printed_output)
    # No further inputs or calculations should occur
    mock_dependencies['mock_input'].assert_not_called()
    mock_dependencies['MockTrajectoryPlanner'].assert_not_called()

# --- Unit tests for main's helper functions (updated) ---

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
    expected_years_earth = years_traveler / math.sqrt(1 - (0.8**2))
    calculated_years_earth = main.calc_time_earth(years_traveler, average_speed_ms)
    assert calculated_years_earth == pytest.approx(expected_years_earth, rel=1e-9)

def test_calc_time_earth_at_speed_of_light():
    """
    Test calc_time_earth when speed is exactly the speed of light.
    main.py returns traveler's time directly in this case.
    """
    assert main.calc_time_earth(1.0, constants.C_LIGHT_MPS) == pytest.approx(1.0)
    # Also ensure a warning is printed
    with patch('builtins.print') as mock_print:
        main.calc_time_earth(1.0, constants.C_LIGHT_MPS)
        mock_print.assert_called_with("Warning: Speed is at or above the speed of light. Relativistic time dilation calculation may be invalid or undefined. Returning traveler's time.")

def test_calc_time_earth_above_speed_of_light():
    """
    Test calc_time_earth when speed is above the speed of light.
    main.py returns traveler's time directly in this case.
    """
    assert main.calc_time_earth(1.0, 1.1 * constants.C_LIGHT_MPS) == pytest.approx(1.0)
    with patch('builtins.print') as mock_print:
        main.calc_time_earth(1.0, 1.1 * constants.C_LIGHT_MPS)
        mock_print.assert_called_with("Warning: Speed is at or above the speed of light. Relativistic time dilation calculation may be invalid or undefined. Returning traveler's time.")


def test_calc_time_earth_zero_speed():
    """
    Test calc_time_earth with zero speed.
    main.py now returns float('inf') for average_speed_ms == 0.
    """
    assert main.calc_time_earth(1.0, 0.0) == float('inf')

def test_calc_time_earth_negative_speed():
    """
    Test calc_time_earth with negative speed.
    main.py should return traveler's time and print a warning.
    """
    with patch('builtins.print') as mock_print:
        calculated_years = main.calc_time_earth(5.0, -100.0)
        assert calculated_years == 5.0
        mock_print.assert_called_with("Warning: Negative speed provided. Returning traveler's time.")

def test_calc_age_whole_numbers():
    """
    Test calc_age with whole number inputs.
    """
    assert main.calc_age(10.0, 25) == 35.0
    assert main.calc_age(0.0, 30) == 30.0

def test_calc_age_fractional_numbers():
    """
    Test calc_age with fractional inputs.
    """
    assert main.calc_age(5.5, 20) == 25.5
    assert main.calc_age(0.5, 20) == 20.5

def test_calc_age_invalid_inputs():
    """
    Test calc_age with invalid input values (negative).
    (Note: type errors for non-numeric input are primarily handled by `int(input())` in `main()` before reaching `calc_age`)
    """
    with pytest.raises(ValueError, match="Starting age must be a non-negative number."):
        main.calc_age(10, -5)
    with pytest.raises(ValueError, match="Travel time must be a non-negative number."):
        main.calc_age(-10.0, 25)
    
    # Test non-numeric types directly on calc_age for robustness, even if main() usually handles this
    with pytest.raises(ValueError): 
        main.calc_age("abc", 25)
    with pytest.raises(ValueError):
        main.calc_age(10, "abc")


def test_calc_arrival_whole_years():
    """
    Test calc_arrival with whole years, expecting exact date addition.
    """
    departure_date = datetime(2024, 1, 1) # Starting on a leap year to test robustness
    arrival_date = main.calc_arrival(departure_date, 10.0)
    # 10 years * 365.25 days/year = 3652.5 days
    # datetime.timedelta handles fractional days, adding 0.5 * 24 = 12 hours
    expected_date = datetime(2034, 1, 1, 12, 0, 0) # Account for the 0.5 day
    assert arrival_date == expected_date

def test_calc_arrival_fractional_years():
    """
    Test calc_arrival with fractional years, comparing against exact timedelta calculation.
    """
    departure_date = datetime(2024, 1, 1) # Start on a leap year
    
    # 0.5 year travel time: 0.5 * 365.25 = 182.625 days.
    arrival_date_half_year = main.calc_arrival(departure_date, 0.5)
    # 182 days, 0.625 * 24 hours = 15 hours
    expected_date_half_year = datetime(2024, 7, 1, 15, 0) 
    assert arrival_date_half_year == expected_date_half_year

def test_calc_arrival_negative_years():
    """
    Test calc_arrival with negative years, expecting a past date.
    """
    departure_date = datetime(2023, 1, 1)
    # -5 years * 365.25 = -1826.25 days.
    # So, departure_date - 1826 days - 0.25 days (6 hours)
    arrival_date = main.calc_arrival(departure_date, -5.0)
    expected_date = datetime(2017, 12, 30, 18, 0) 
    assert arrival_date == expected_date

def test_calc_arrival_invalid_date_type():
    """
    Test calc_arrival with invalid departure date types, expecting None.
    main.py handles TypeError and returns None.
    """
    assert main.calc_arrival(None, 10) is None
    assert main.calc_arrival("2023-01-01", 10) is None # Passes non-datetime type

def test_calc_arrival_overflow_error():
    """
    Test calc_arrival with extremely large years to trigger OverflowError.
    """
    departure_date = datetime(2000, 1, 1)
    # Large enough to cause overflow, e.g., max timedelta allows ~27000 years, go beyond
    very_large_years = 1_000_000.0 
    
    with patch('builtins.print') as mock_print:
        arrival_date = main.calc_arrival(departure_date, very_large_years)
        assert arrival_date is None
        mock_print.assert_called_with("Warning: Arrival date calculation resulted in a date out of range due to extremely long travel time.")


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
    with pytest.raises(ValueError): # main.py directly calls datetime.strptime
        main.convert_date("invalid-date")
    with pytest.raises(ValueError):
        main.convert_date("2023/01/15") # Wrong separator
    with pytest.raises(ValueError):
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
