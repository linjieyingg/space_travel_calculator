import pytest
from unittest.mock import patch, MagicMock, call
from datetime import datetime, timedelta
import math
import sys

# Assume main.py is in the same directory as tests/
# We need to import main to test it, and its dependencies for mocking
import main
from constants import C_LIGHT_MPS # Need this for specific constant assertions

# --- Fixtures for main.py tests ---
@pytest.fixture
def mock_dependencies():
    with patch('main.celestial_data.get_celestial_body_data') as mock_get_celestial_body_data, \
         patch('main.celestial_data.get_all_solar_system_destinations') as mock_get_all_solar_system_destinations, \
         patch('main.trajectory_planner.TrajectoryPlanner') as MockTrajectoryPlanner, \
         patch('main.fuel_optimizer.optimize_fuel_for_trajectory') as mock_optimize_fuel_for_trajectory, \
         patch('main.fuel_calc.calculate_fuel_cost') as mock_calculate_fuel_cost, \
         patch('main.travel_logger.save_travel_log') as mock_save_travel_log, \
         patch('builtins.input') as mock_input, \
         patch('builtins.print') as mock_print, \
         patch('main.c.is_valid_date') as mock_is_valid_date, \
         patch('main.c.is_valid_age') as mock_is_valid_age:

        # Mock celestial data
        mock_get_celestial_body_data.side_effect = lambda name: {
            'earth': {'name': 'Earth', 'mass': 5.972e24, 'radius': 6.371e6, 'gravitational_parameter_mu': 3.986e14, 'semi_major_axis_from_sun': 1.496e11, 'average_distance_from_earth': 0},
            'mars': {'name': 'Mars', 'mass': 6.39e23, 'radius': 3.3895e6, 'gravitational_parameter_mu': 4.2828e13, 'semi_major_axis_from_sun': 2.279e11, 'average_distance_from_earth': 7.834e10},
            'sun': {'name': 'Sun', 'mass': 1.989e30, 'gravitational_parameter_mu': 1.327e20}
        }.get(name.lower())

        mock_get_all_solar_system_destinations.return_value = [
            {'name': 'Mars', 'distance': 7.834e10},
            {'name': 'Jupiter', 'distance': 6.287e11}
        ]

        # Mock TrajectoryPlanner
        mock_planner_instance = MockTrajectoryPlanner.return_value
        # Updated to use plan_trajectory with trajectory_type
        mock_planner_instance.plan_trajectory.return_value = {
            'success': True,
            'total_delta_v_mps': 10000.0,
            'total_travel_time_seconds': 100 * 24 * 3600, # 100 days
            'transfer_type': 'Hohmann',
            'average_hohmann_speed_ms': 500000.0 # Example average speed
        }

        # Mock fuel_optimizer
        mock_optimize_fuel_for_trajectory.return_value = 50000.0 # kg

        # Mock fuel_calc
        mock_calculate_fuel_cost.return_value = {
            'total_fuel_mass_needed': 50000.0,
            'total_cost': 500000.0 # $
        }

        # Mock checks module functions
        mock_is_valid_date.return_value = True
        mock_is_valid_age.return_value = True

        mocks = {
            "mock_get_celestial_body_data": mock_get_celestial_body_data,
            "mock_get_all_solar_system_destinations": mock_get_all_solar_system_destinations,
            "MockTrajectoryPlanner": MockTrajectoryPlanner,
            "mock_planner_instance": mock_planner_instance,
            "mock_optimize_fuel_for_trajectory": mock_optimize_fuel_for_trajectory,
            "mock_calculate_fuel_cost": mock_calculate_fuel_cost,
            "mock_save_travel_log": mock_save_travel_log,
            "mock_input": mock_input,
            "mock_print": mock_print,
            "mock_is_valid_date": mock_is_valid_date,
            "mock_is_valid_age": mock_is_valid_age
        }
        yield mocks

# --- Tests for main.py ---

def test_main_function_integrates_all_components_correctly(mock_dependencies):
    mocks = mock_dependencies

    # Simulate user input
    mocks["mock_input"].side_effect = [
        "Earth",  # Departure Body
        "Mars",   # Arrival Body
        "10000",  # Spacecraft dry mass (kg)
        "3000",   # Engine specific impulse (s)
        "2024-12-01", # Departure date
        "30",     # User age
        "10"      # Fuel price per unit mass ($/kg)
    ]

    main.main()

    # Assertions on mock calls
    mocks["mock_get_all_solar_system_destinations"].assert_called_once()
    mocks["mock_get_celestial_body_data"].assert_any_call("Earth")
    mocks["mock_get_celestial_body_data"].assert_any_call("Mars")
    mocks["MockTrajectoryPlanner"].assert_called_once()
    # Updated to assert plan_trajectory call
    mocks["mock_planner_instance"].plan_trajectory.assert_called_once_with(
        departure_body_name="Earth",
        arrival_body_name="Mars",
        trajectory_type="Hohmann" # Added argument
    )
    mocks["mock_optimize_fuel_for_trajectory"].assert_called_once_with(
        start_body="Earth",
        end_body="Mars",
        trajectory_type="Hohmann",
        spacecraft_dry_mass=10000.0,
        engine_specific_impulse=3000.0
    )
    mocks["mock_calculate_fuel_cost"].assert_called_once_with(
        total_fuel_mass_needed=50000.0,
        fuel_price_per_unit=10.0
    )
    mocks["mock_save_travel_log"].assert_called_once()
    # Check specific arguments for save_travel_log. travel_time is in seconds.
    # average_speed_for_relativistic_calc comes from `trajectory_result.get('average_hohmann_speed_ms', 0.0)`
    # The fixture sets `average_hohmann_speed_ms` to 500000.0
    mocks["mock_save_travel_log"].assert_called_once_with(
        source_planet='Earth',
        destination_planet='Mars',
        speed=500000.0, # From average_hohmann_speed_ms
        travel_time=float(100 * 24 * 3600), # 100 days in seconds
        delta_v_required=10000.0,
        fuel_mass_needed=50000.0,
        transfer_type='Hohmann'
    )

    # Assert printed output (checking for key phrases, not exact string due to formatting)
    output_calls = [c.args[0] for c in mocks["mock_print"].call_args_list if c.args]
    assert "Welcome to the Space Travel Calculator!" in output_calls[0]
    assert "Travel Summary" in "\n".join(output_calls)
    assert "Departure Body: Earth" in "\n".join(output_calls)
    assert "Estimated Fuel Mass Needed: 50000.00 kg" in "\n".join(output_calls)
    assert "Total Fuel Cost: $500,000.00" in "\n".join(output_calls)
    # Check relativistic time calculation result (100 days traveler, 100 days * lorentz_factor earth)
    # traveler_time_years = 100 days / 365.25 days/year
    traveler_time_years = 100 / 365.25
    average_speed_ms = 500000.0
    lorentz_factor = 1 / math.sqrt(1 - (average_speed_ms / C_LIGHT_MPS)**2)
    expected_earth_frame_years = traveler_time_years * lorentz_factor
    # The printed output is formatted to 2 decimal places, so check for that.
    assert f"Travel Time (Traveler's Frame): {traveler_time_years:.2f} years" in "\n".join(output_calls)
    assert f"Travel Time (Earth's Frame): {expected_earth_frame_years:.2f} years" in "\n".join(output_calls)
    assert "Your Age Upon Arrival: 30 years" in "\n".join(output_calls) # Age is rounded to 0f
    # Calculate expected arrival date: 2024-12-01 + (100 * 24 * 3600 / (365.25 * 24 * 3600)) * 365.25 days
    # Which is effectively 100 days from 2024-12-01 = 2025-03-11
    assert "Estimated Arrival Date: 2025-03-11" in "\n".join(output_calls)


def test_main_function_trajectory_planner_fails(mock_dependencies):
    mocks = mock_dependencies
    # Updated to mock plan_trajectory
    mocks["mock_planner_instance"].plan_trajectory.return_value = {
        'success': False,
        'error_message': 'Departure body data missing'
    }
    mocks["mock_input"].side_effect = [
        "Earth", "Mars", "10000", "3000", "2024-12-01", "30", "10"
    ]

    main.main()

    output_calls = [c.args[0] for c in mocks["mock_print"].call_args_list if c.args]
    assert "Error planning trajectory: Departure body data missing" in "\n".join(output_calls)
    mocks["mock_optimize_fuel_for_trajectory"].assert_not_called()
    mocks["mock_calculate_fuel_cost"].assert_not_called()
    mocks["mock_save_travel_log"].assert_not_called()

def test_main_function_fuel_optimizer_raises_value_error(mock_dependencies):
    mocks = mock_dependencies
    mocks["mock_optimize_fuel_for_trajectory"].side_effect = ValueError("Dry mass must be positive.")
    mocks["mock_input"].side_effect = [
        "Earth", "Mars", "10000", "3000", "2024-12-01", "30", "10"
    ]

    main.main()

    output_calls = [c.args[0] for c in mocks["mock_print"].call_args_list if c.args]
    assert "Input error during fuel optimization: Dry mass must be positive." in "\n".join(output_calls)
    mocks["mock_calculate_fuel_cost"].assert_not_called()
    mocks["mock_save_travel_log"].assert_not_called()

def test_main_function_fuel_optimizer_raises_runtime_error(mock_dependencies):
    mocks = mock_dependencies
    mocks["mock_optimize_fuel_for_trajectory"].side_effect = RuntimeError("Internal fuel calculation error.")
    mocks["mock_input"].side_effect = [
        "Earth", "Mars", "10000", "3000", "2024-12-01", "30", "10"
    ]

    main.main()

    output_calls = [c.args[0] for c in mocks["mock_print"].call_args_list if c.args]
    assert "Operational error during fuel optimization: Internal fuel calculation error." in "\n".join(output_calls)
    mocks["mock_calculate_fuel_cost"].assert_not_called()
    mocks["mock_save_travel_log"].assert_not_called()

def test_main_function_fuel_cost_calculation_fails(mock_dependencies):
    mocks = mock_dependencies
    mocks["mock_calculate_fuel_cost"].return_value = None # Simulate failure
    mocks["mock_input"].side_effect = [
        "Earth", "Mars", "10000", "3000", "2024-12-01", "30", "10"
    ]

    main.main()

    output_calls = [c.args[0] for c in mocks["mock_print"].call_args_list if c.args]
    assert "Error calculating fuel cost. Check input parameters or the result from fuel mass calculation." in "\n".join(output_calls)
    mocks["mock_save_travel_log"].assert_not_called()

def test_main_function_get_all_solar_system_destinations_fails(mock_dependencies):
    mocks = mock_dependencies
    mocks["mock_get_all_solar_system_destinations"].return_value = [] # Simulate failure

    main.main()

    output_calls = [c.args[0] for c in mocks["mock_print"].call_args_list if c.args]
    assert "Error: Could not retrieve solar system destinations. Exiting." in "\n".join(output_calls)
    mocks["mock_input"].assert_not_called() # No inputs should be requested
    mocks["MockTrajectoryPlanner"].assert_not_called()
    mocks["mock_optimize_fuel_for_trajectory"].assert_not_called()
    mocks["mock_calculate_fuel_cost"].assert_not_called()
    mocks["mock_save_travel_log"].assert_not_called()

# Test cases for main's helper functions
def test_calc_time_earth_relativistic_zero_speed():
    result = main.calc_time_earth(1.0, 0.0)
    assert result == float('inf')

def test_calc_time_earth_relativistic_speed_of_light():
    # Should print warning and return traveler's time
    with patch('builtins.print') as mock_print:
        result = main.calc_time_earth(10.0, C_LIGHT_MPS)
        assert result == 10.0
        mock_print.assert_called_once_with(
            "Warning: Speed is at or above the speed of light. Relativistic time dilation calculation may be invalid or undefined. Returning traveler's time."
        )

def test_calc_time_earth_relativistic_negative_speed():
    # Should print warning and return traveler's time
    with patch('builtins.print') as mock_print:
        result = main.calc_time_earth(5.0, -100.0)
        assert result == 5.0
        mock_print.assert_called_once_with(
            "Warning: Negative speed provided. Returning traveler's time."
        )

def test_calc_age_valid_inputs():
    assert main.calc_age(5.0, 30) == 35.0
    assert main.calc_age(0.5, 25) == 25.5

def test_calc_age_invalid_inputs():
    with pytest.raises(ValueError, match="Starting age must be a non-negative number."):
        main.calc_age(5.0, -10)
    with pytest.raises(ValueError, match="Travel time must be a non-negative number."):
        main.calc_age(-2.0, 30)
    with pytest.raises(ValueError, match="Starting age must be a non-negative number."):
        main.calc_age(5.0, "invalid")

def test_calc_arrival_valid_inputs():
    departure = datetime(2024, 12, 1)
    years_earth_frame = 0.25 # approx 3 months
    # The actual main.py code uses float days directly: timedelta(days=days_earth_frame)
    days_for_calc = years_earth_frame * 365.25
    expected_arrival_exact = departure + timedelta(days=days_for_calc)
    result = main.calc_arrival(departure, years_earth_frame)
    assert result == expected_arrival_exact

def test_calc_arrival_long_travel_time_overflow():
    departure = datetime(2024, 1, 1)
    years_earth_frame = 1000000000000.0 # Extremely long
    with patch('builtins.print') as mock_print:
        result = main.calc_arrival(departure, years_earth_frame)
        assert result is None
        mock_print.assert_called_once_with(
            "Warning: Arrival date calculation resulted in a date out of range due to extremely long travel time."
        )

def test_convert_date_valid():
    date_str = "2023-01-15"
    expected_date = datetime(2023, 1, 15)
    assert main.convert_date(date_str) == expected_date

def test_convert_date_invalid_format():
    date_str = "15/01/2023"
    with pytest.raises(ValueError):
        main.convert_date(date_str)

def test_main_function_average_speed_fallback(mock_dependencies):
    mocks = mock_dependencies
    # Simulate trajectory planner returning 0 for average_hohmann_speed_ms
    # Updated to mock plan_trajectory
    mocks["mock_planner_instance"].plan_trajectory.return_value = {
        'success': True,
        'total_delta_v_mps': 10000.0,
        'total_travel_time_seconds': 100 * 24 * 3600, # 100 days
        'transfer_type': 'Hohmann',
        'average_hohmann_speed_ms': 0.0 # This will trigger the fallback
    }

    mocks["mock_input"].side_effect = [
        "Earth", "Mars", "10000", "3000", "2024-12-01", "30", "10"
    ]

    main.main()

    output_calls = [c.args[0] for c in mocks["mock_print"].call_args_list if c.args]
    assert "Warning: Average Hohmann speed is zero or negative. Using a placeholder speed for relativistic calculation (0.1% of C)." in "\n".join(output_calls)
    
    # Assert save_travel_log received the fallback speed
    # Fallback speed is 0.001 * C_LIGHT_MPS
    mocks["mock_save_travel_log"].assert_called_once()
    logged_speed = mocks["mock_save_travel_log"].call_args[1]['speed']
    assert logged_speed == pytest.approx(0.001 * C_LIGHT_MPS)