import pytest
from datetime import datetime, timedelta
import math
from unittest.mock import patch, MagicMock

# Import functions directly from main for easier testing
from main import (
    calc_time_earth, calc_age, calc_arrival, convert_date, main,
    FUEL_PRICE_PER_UNIT, SOLAR_SYSTEM_BODIES
)

# Mock constants from the constants module as main depends on it
class MockConstants:
    C_LIGHT_MPS = 299792458.0
    G_GRAVITATIONAL = 6.67430e-11

# Mock data for celestial bodies that main.py's main() function will request
MOCK_EARTH_DATA = {
    "mass": 5.9722e24, "radius": 6.3781e6, "gravitational_parameter_mu": 3.986004418e14,
    "semi_major_axis_from_sun": 1.49598023e11
}
MOCK_MARS_DATA = {
    "mass": 6.4169e23, "radius": 3.3895e6, "gravitational_parameter_mu": 4.282837e13,
    "semi_major_axis_from_sun": 2.279392e11
}
MOCK_SUN_DATA = {
    "mass": 1.9885e30, "radius": 6.9634e8, "gravitational_parameter_mu": 1.32712440018e20,
    "semi_major_axis_from_sun": 0.0
}

# Mock return for successful trajectory planning
MOCK_TRAJECTORY_SUCCESS = {
    'success': True,
    'total_delta_v_mps': 10000.0,
    'travel_time_traveler_frame_years': 0.7, # Roughly 255 days
    'average_hohmann_speed_ms': 50000.0, # ~0.000167c
    'transfer_type': 'Hohmann Transfer',
    'departure_body_name': 'Earth',
    'arrival_body_name': 'Mars'
}

# Mock return for fuel calculation
MOCK_FUEL_MASS = 500000.0 # kg
MOCK_FUEL_COST_RESULT = {'total_fuel_mass_needed': MOCK_FUEL_MASS, 'total_cost': MOCK_FUEL_MASS * FUEL_PRICE_PER_UNIT}


@pytest.fixture(autouse=True)
def mock_constants_in_main(mocker):
    """Automatically mock constants in main.py for consistent testing."""
    mocker.patch('main.constants.C_LIGHT_MPS', MockConstants.C_LIGHT_MPS)
    mocker.patch('main.constants.G_GRAVITATIONAL', MockConstants.G_GRAVITATIONAL)


# --- Test calc_time_earth function ---
@pytest.mark.parametrize("traveler_years, speed_ms, expected_earth_years", [
    (1.0, 0.0, 1.0), # No speed, no dilation
    (1.0, -100.0, 1.0), # Negative speed, no dilation
    (1.0, MockConstants.C_LIGHT_MPS / 2, 1.1547005383792515), # v=0.5c, gamma = 1/sqrt(1-0.25) = 1/sqrt(0.75) ~ 1.1547
    (1.0, MockConstants.C_LIGHT_MPS * 0.99, 7.088812037946979), # v=0.99c, high dilation
    (1.0, MockConstants.C_LIGHT_MPS * 0.999999, 707.1068800560875), # v=0.999999c, very high dilation
    (5.0, MockConstants.C_LIGHT_MPS * 0.8, 8.333333333333334), # 5 years @ 0.8c
])
def test_calc_time_earth_relativistic(traveler_years, speed_ms, expected_earth_years):
    """Test calc_time_earth with various speeds, including relativistic effects."""
    result = calc_time_earth(traveler_years, speed_ms)
    assert math.isclose(result, expected_earth_years, rel_tol=1e-9)


def test_calc_time_earth_at_light_speed():
    """Test calc_time_earth when speed is exactly the speed of light (should return traveler's time)."""
    result = calc_time_earth(1.0, MockConstants.C_LIGHT_MPS)
    assert result == 1.0

def test_calc_time_earth_above_light_speed(capsys):
    """Test calc_time_earth when speed exceeds the speed of light (should return traveler's time and print warning)."""
    result = calc_time_earth(1.0, MockConstants.C_LIGHT_MPS * 1.1)
    assert result == 1.0
    captured = capsys.readouterr()
    assert "Warning: Speed is at or above the speed of light." in captured.out


# --- Test calc_age function ---
@pytest.mark.parametrize("travel_years, start_age, expected_arrival_age", [
    (10.5, 20, 30), # Basic addition, floor
    (0.0, 25, 25),  # Zero travel time
    (5.99, 18, 23), # Floor test
    (6.01, 18, 24), # Floor test
    (0.001, 74, 74), # Small travel time
])
def test_calc_age(travel_years, start_age, expected_arrival_age):
    """Test calc_age with various travel times and starting ages."""
    assert calc_age(travel_years, start_age) == expected_arrival_age


# --- Test calc_arrival function ---
@pytest.mark.parametrize("start_date_str, earth_years, expected_arrival_date_str", [
    ("2024-01-01", 1.0, "2025-01-01"),
    ("2024-01-01", 0.5, "2024-07-01"), # ~182.625 days
    ("2024-01-01", 10.0, "2034-01-01"),
    ("2023-03-15", 2.75, "2025-12-15"), # Accounting for leap years (~365.25 days/year)
])
def test_calc_arrival_basic(start_date_str, earth_years, expected_arrival_date_str):
    """Test calc_arrival with valid departure dates and earth frame travel times."""
    departure_date = datetime.strptime(start_date_str, '%Y-%m-%d')
    expected_arrival_date = datetime.strptime(expected_arrival_date_str, '%Y-%m-%d')
    
    # The expected date calculation needs to be precise with 365.25
    calculated_arrival_date = departure_date + timedelta(days=earth_years * 365.25)
    
    result = calc_arrival(departure_date, earth_years)
    # Allow for minor differences due to float precision, compare year, month, day
    assert result.year == calculated_arrival_date.year
    assert result.month == calculated_arrival_date.month
    assert result.day == calculated_arrival_date.day


def test_calc_arrival_large_years_overflow(capsys):
    """Test calc_arrival with extremely large years to trigger OverflowError."""
    departure_date = datetime(2000, 1, 1)
    # A number large enough to cause OverflowError in timedelta (e.g., > 10^10 days)
    extremely_large_years = 1e18 # This will definitely overflow if multiplied by 365.25

    # We need to mock timedelta or the datetime.__add__ method to simulate OverflowError
    # Patch datetime.timedelta's __new__ to raise OverflowError when created with large days
    original_timedelta = datetime.__add__

    def mock_add(self, other):
        if isinstance(other, timedelta) and other.days > 1e10:
            raise OverflowError("simulated timedelta overflow")
        return original_timedelta(self, other)

    with patch('datetime.datetime.__add__', new=mock_add):
        result = calc_arrival(departure_date, extremely_large_years)
        assert result is None
        captured = capsys.readouterr()
        assert "Warning: Arrival date calculation resulted in a date out of range" in captured.out


def test_calc_arrival_unexpected_error(mocker, capsys):
    """Test calc_arrival with an unexpected error during date addition."""
    mocker.patch('datetime.timedelta', side_effect=Exception("Mocked general error"))
    result = calc_arrival(datetime(2024, 1, 1), 1.0)
    assert result is None
    captured = capsys.readouterr()
    assert "Warning: An unexpected error occurred calculating arrival date" in captured.out


# --- Test convert_date function ---
@pytest.mark.parametrize("date_str, expected_date", [
    ("2024-12-01", datetime(2024, 12, 1)),
    ("1999-07-15", datetime(1999, 7, 15)),
])
def test_convert_date_valid(date_str, expected_date):
    """Test convert_date with valid date strings."""
    assert convert_date(date_str) == expected_date


@pytest.mark.parametrize("invalid_date_str", [
    "01-12-2024", # Wrong format
    "2024/12/01", # Wrong separator
    "2024-13-01", # Invalid month
    "2024-02-30", # Invalid day for month
    "not-a-date", # Not a date
    "",           # Empty string
])
def test_convert_date_invalid_format(invalid_date_str):
    """Test convert_date with invalid date string formats, expecting ValueError."""
    with pytest.raises(ValueError):
        convert_date(invalid_date_str)


# --- Test main function (integration and orchestration) ---

@pytest.fixture
def mock_main_dependencies(mocker):
    """Fixture to mock external dependencies for the main function."""
    # Mock internal check functions from 'checks' module
    mocker.patch('main.c.is_valid_date', return_value=True)
    mocker.patch('main.c.is_valid_age', return_value=True)

    # Mock celestial_data lookups
    mocker.patch('main.celestial_data.get_celestial_body_data', side_effect=[
        MOCK_EARTH_DATA, # For source planet input validation
        MOCK_MARS_DATA   # For destination planet input validation
    ])
    
    # Mock TrajectoryPlanner and its method
    mock_planner_instance = MagicMock()
    mock_planner_instance.plan_hohmann_trajectory.return_value = MOCK_TRAJECTORY_SUCCESS
    mocker.patch('main.trajectory_planner.TrajectoryPlanner', return_value=mock_planner_instance)

    # Mock propulsion_system
    mocker.patch('main.propulsion_system.calculate_required_fuel_mass', return_value=MOCK_FUEL_MASS)

    # Mock fuel_calc
    mocker.patch('main.fuel_calc.calculate_fuel_cost', return_value=MOCK_FUEL_COST_RESULT)

    # Mock travel_logger
    mock_save_log = mocker.patch('main.travel_logger.save_travel_log')
    return {
        'mock_planner_instance': mock_planner_instance,
        'mock_save_log': mock_save_log
    }


def test_main_successful_run(mock_main_dependencies, capsys, mocker):
    """Test a complete successful execution of the main function with valid inputs."""
    # Prepare inputs for the main function's interactive prompts
    mock_inputs = [
        "Earth",
        "Mars",
        "10000", # spacecraft_dry_mass
        "3000",  # engine_isp
        "2025-01-01", # date_str
        "30"     # user_age
    ]
    mocker.patch('builtins.input', side_effect=mock_inputs)

    main()

    captured = capsys.readouterr()
    # Assertions for output prints
    assert "Welcome to the Interplanetary Travel Calculator!" in captured.out
    assert f"Source Planet: Earth" in captured.out
    assert f"Destination: Mars" in captured.out
    assert f"Required Delta-V: {MOCK_TRAJECTORY_SUCCESS['total_delta_v_mps']:.3f} m/s" in captured.out
    assert f"Required Fuel Mass: {MOCK_FUEL_MASS:.3f} kg" in captured.out
    assert f"Estimated Fuel Cost: ${MOCK_FUEL_COST_RESULT['total_cost']:,.2f}" in captured.out
    assert "Travel details successfully logged." in captured.out

    # Assert that key functions from other modules were called
    mock_main_dependencies['mock_planner_instance'].plan_hohmann_trajectory.assert_called_once_with(
        "Earth", "Mars"
    )
    mock_main_dependencies['mock_save_log'].assert_called_once_with(
        source_planet="Earth",
        destination_planet="Mars",
        speed=MOCK_TRAJECTORY_SUCCESS['average_hohmann_speed_ms'],
        travel_time=MOCK_TRAJECTORY_SUCCESS['travel_time_traveler_frame_years'],
        delta_v_required=MOCK_TRAJECTORY_SUCCESS['total_delta_v_mps'],
        fuel_mass_needed=MOCK_FUEL_MASS,
        transfer_type=MOCK_TRAJECTORY_SUCCESS['transfer_type']
    )
    assert 'celestial_data.get_celestial_body_data' in mocker.patch.get_original('main.celestial_data.get_celestial_body_data').call_args_list[0][0][0].lower()
    assert 'celestial_data.get_celestial_body_data' in mocker.patch.get_original('main.celestial_data.get_celestial_body_data').call_args_list[1][0][0].lower()


def test_main_trajectory_planning_failure(mock_main_dependencies, capsys, mocker):
    """Test main function when trajectory planning fails."""
    mock_main_dependencies['mock_planner_instance'].plan_hohmann_trajectory.return_value = {
        'success': False, 'error': 'Cannot reach that destination.'
    }
    mock_inputs = [
        "Earth",
        "Sun", # An impossible destination for Hohmann
        "10000",
        "3000",
        "2025-01-01",
        "30"
    ]
    mocker.patch('builtins.input', side_effect=mock_inputs)
    
    # Mock celestial_data get_celestial_body_data to correctly return for Earth and Sun
    mocker.patch('main.celestial_data.get_celestial_body_data', side_effect=[
        MOCK_EARTH_DATA, # For source planet
        MOCK_SUN_DATA    # For destination planet
    ])

    main()

    captured = capsys.readouterr()
    assert "Error during trajectory planning: Cannot reach that destination." in captured.out
    assert "Cannot proceed with calculations." in captured.out
    mock_main_dependencies['mock_save_log'].assert_not_called()


def test_main_fuel_calculation_failure(mock_main_dependencies, capsys, mocker):
    """Test main function when fuel calculation raises an error."""
    mocker.patch('main.propulsion_system.calculate_required_fuel_mass', side_effect=ValueError("Invalid dry mass"))

    mock_inputs = [
        "Earth",
        "Mars",
        "10000",
        "3000",
        "2025-01-01",
        "30"
    ]
    mocker.patch('builtins.input', side_effect=mock_inputs)

    main()

    captured = capsys.readouterr()
    assert "Error calculating fuel mass or cost: Invalid dry mass." in captured.out
    # Fuel cost should be 0.0 if fuel mass calculation failed.
    assert "Estimated Fuel Cost: $0.00" in captured.out
    # Logger should still be called, but with 0 fuel mass
    mock_main_dependencies['mock_save_log'].assert_called_once()
    args, _ = mock_main_dependencies['mock_save_log'].call_args
    assert args['fuel_mass_needed'] == 0.0


def test_main_travel_logging_failure(mock_main_dependencies, capsys, mocker):
    """Test main function when travel logging fails."""
    mock_main_dependencies['mock_save_log'].side_effect = Exception("Logger error")

    mock_inputs = [
        "Earth",
        "Mars",
        "10000",
        "3000",
        "2025-01-01",
        "30"
    ]
    mocker.patch('builtins.input', side_effect=mock_inputs)

    main()

    captured = capsys.readouterr()
    assert "Error saving travel log: Logger error" in captured.out
    # Ensure other parts of the summary are still printed
    assert "--- Travel Summary ---" in captured.out


def test_main_invalid_input_reprompts(mock_main_dependencies, capsys, mocker):
    """Test that main function reprompts for invalid user inputs."""
    # Scenario: Invalid dry mass -> Valid dry mass
    mock_inputs = [
        "Earth",
        "Mars",
        "-100",    # Invalid: negative dry mass
        "invalid", # Invalid: non-numeric dry mass
        "10000",   # Valid dry mass
        "3000",
        "2025-01-01",
        "30"
    ]
    mocker.patch('builtins.input', side_effect=mock_inputs)

    main()

    captured = capsys.readouterr()
    assert "Spacecraft dry mass must be a positive number. Please re-enter." in captured.out
    assert "Invalid input. Spacecraft dry mass must be a numeric value. Please re-enter." in captured.out
    # Assert a successful log, implying that eventually valid input was provided and processed
    mock_main_dependencies['mock_save_log'].assert_called_once()
