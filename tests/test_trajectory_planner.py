import pytest
from unittest.mock import patch, MagicMock, call
from datetime import timedelta, datetime
import math

# Use relative imports as per instructions
from trajectory_planner import TrajectoryPlanner
import celestial_data
import orbital_mechanics
import constants
# Assuming ephemeris module exists and is imported by trajectory_planner for dynamic data
# import ephemeris # Not directly imported here, but patched via trajectory_planner.ephemeris

# Fixture for mocking celestial_data.get_celestial_body_data
@pytest.fixture
def mock_celestial_data_for_planner():
    with patch('trajectory_planner.celestial_data.get_celestial_body_data') as mock_get_data:
        # Added 'orbital_period_days' as requested
        mock_get_data.side_effect = lambda name: {
            'sun': {'name': 'Sun', 'mass': 1.989e30, 'radius': 6.957e8, 'gravitational_parameter_mu': 1.327e20},
            'earth': {'name': 'Earth', 'mass': 5.972e24, 'radius': 6.371e6, 'gravitational_parameter_mu': 3.986e14, 'semi_major_axis_from_sun': 1.496e11, 'orbital_period_days': 365.25},
            'mars': {'name': 'Mars', 'mass': 6.39e23, 'radius': 3.389e6, 'gravitational_parameter_mu': 4.283e13, 'semi_major_axis_from_sun': 2.279e11, 'orbital_period_days': 687},
            'jupiter': {'name': 'Jupiter', 'mass': 1.898e27, 'radius': 6.9911e7, 'gravitational_parameter_mu': 1.267e17, 'semi_major_axis_from_sun': 7.785e11, 'orbital_period_days': 4333}
        }.get(name.lower())
        yield mock_get_data

# Fixture for mocking orbital_mechanics functions
@pytest.fixture
def mock_orbital_mechanics():
    with patch('trajectory_planner.orbital_mechanics.calculate_hohmann_transfer_delta_v') as mock_delta_v, \
         patch('trajectory_planner.orbital_mechanics.calculate_hohmann_transfer_time_of_flight') as mock_time_of_flight:
        mock_delta_v.return_value = (1000.0, 2000.0, 3000.0) # Delta-V1, Delta-V2, Total Delta-V
        mock_time_of_flight.return_value = 100 * 24 * 3600 # 100 days in seconds
        yield mock_delta_v, mock_time_of_flight

# New Fixture for mocking ephemeris.get_body_orbital_state
@pytest.fixture
def mock_ephemeris_get_body_orbital_state():
    with patch('trajectory_planner.ephemeris.get_body_orbital_state') as mock_get_state:
        # Default behavior: return distinct, valid distances for Earth and Mars
        def side_effect_ephemeris(body_name: str, departure_date: datetime):
            # Simulate slight variation for ephemeris data compared to static semi-major axis
            if body_name.lower() == 'earth':
                return {'distance_from_sun_m': 1.500e11} # Slightly different from 1.496e11
            elif body_name.lower() == 'mars':
                return {'distance_from_sun_m': 2.300e11} # Slightly different from 2.279e11
            return None # For other bodies or if no specific distance is defined
        mock_get_state.side_effect = side_effect_ephemeris
        yield mock_get_state

# --- Tests for TrajectoryPlanner.__init__ ---

def test_trajectory_planner_init_success(mock_celestial_data_for_planner):
    """
    Test successful initialization of TrajectoryPlanner when Sun data is valid.
    """
    planner = TrajectoryPlanner()
    assert planner.mu_sun == pytest.approx(1.327e20)
    mock_celestial_data_for_planner.assert_called_with('Sun')

def test_trajectory_planner_init_no_sun_data(mock_celestial_data_for_planner):
    """
    Test TrajectoryPlanner initialization when Sun data is not found.
    """
    mock_celestial_data_for_planner.side_effect = lambda name: None if name.lower() == 'sun' else {}
    with pytest.raises(ValueError, match="Data for 'Sun' not found"):
        TrajectoryPlanner()

def test_trajectory_planner_init_missing_mu_sun(mock_celestial_data_for_planner):
    """
    Test TrajectoryPlanner initialization when 'gravitational_parameter_mu' is missing for Sun.
    """
    mock_celestial_data_for_planner.side_effect = lambda name: {'name': 'Sun', 'mass': 1.989e30} if name.lower() == 'sun' else {}
    with pytest.raises(AttributeError, match="'gravitational_parameter_mu' not found for 'Sun'"):
        TrajectoryPlanner()

def test_trajectory_planner_init_invalid_mu_sun_type(mock_celestial_data_for_planner):
    """
    Test TrajectoryPlanner initialization when 'gravitational_parameter_mu' for Sun is of an invalid type.
    """
    mock_celestial_data_for_planner.side_effect = lambda name: {'name': 'Sun', 'gravitational_parameter_mu': 'invalid'} if name.lower() == 'sun' else {}
    with pytest.raises(ValueError, match="Invalid or non-positive 'gravitational_parameter_mu'"):
        TrajectoryPlanner()

def test_trajectory_planner_init_non_positive_mu_sun(mock_celestial_data_for_planner):
    """
    Test TrajectoryPlanner initialization when 'gravitational_parameter_mu' for Sun is non-positive.
    """
    mock_celestial_data_for_planner.side_effect = lambda name: {'name': 'Sun', 'gravitational_parameter_mu': -1.0} if name.lower() == 'sun' else {}
    with pytest.raises(ValueError, match="Invalid or non-positive 'gravitational_parameter_mu'"):
        TrajectoryPlanner()

# --- Tests for TrajectoryPlanner.plan_trajectory (Hohmann type) ---

def test_plan_trajectory_hohmann_success_earth_mars_with_ephemeris(
    mock_celestial_data_for_planner, mock_orbital_mechanics, mock_ephemeris_get_body_orbital_state
):
    """
    Test successful planning of an Earth to Mars Hohmann transfer using plan_trajectory
    with a departure_date, leveraging ephemeris data for distances.
    """
    planner = TrajectoryPlanner()
    departure_date = datetime(2025, 1, 1)
    result = planner.plan_trajectory("Earth", "Mars", trajectory_type="Hohmann", departure_date=departure_date)

    assert result["success"] is True
    assert result["departure_body"] == "Earth"
    assert result["arrival_body"] == "Mars"
    assert result["transfer_type"] == "Hohmann Transfer"
    assert result["total_delta_v_mps"] == pytest.approx(3000.0)
    assert result["travel_time_days"] == pytest.approx(100.0)
    assert result["travel_time_h_m_s"] == "100 days, 00 hours, 00 minutes, 00 seconds"
    assert result["error"] is None

    # Verify ephemeris was called for both bodies with the correct date
    mock_ephemeris_get_body_orbital_state.assert_has_calls([
        call('Earth', departure_date),
        call('Mars', departure_date)
    ], any_order=True)

    # Verify orbital mechanics functions were called with the distances from ephemeris
    ephemeris_earth_distance = mock_ephemeris_get_body_orbital_state.side_effect('Earth', departure_date)['distance_from_sun_m']
    ephemeris_mars_distance = mock_ephemeris_get_body_orbital_state.side_effect('Mars', departure_date)['distance_from_sun_m']

    mock_orbital_mechanics[0].assert_called_once_with(
        planner.mu_sun,
        ephemeris_earth_distance,
        ephemeris_mars_distance
    )
    mock_orbital_mechanics[1].assert_called_once_with(
        planner.mu_sun,
        ephemeris_earth_distance,
        ephemeris_mars_distance
    )

def test_plan_trajectory_hohmann_success_earth_mars_no_departure_date_fallback_to_static(
    mock_celestial_data_for_planner, mock_orbital_mechanics, mock_ephemeris_get_body_orbital_state
):
    """
    Test successful planning of an Earth to Mars Hohmann transfer when no departure_date is provided,
    falling back to celestial_data's semi_major_axis_from_sun.
    """
    planner = TrajectoryPlanner()
    result = planner.plan_trajectory("Earth", "Mars", trajectory_type="Hohmann") # No departure_date

    assert result["success"] is True
    assert result["departure_body"] == "Earth"
    assert result["arrival_body"] == "Mars"
    assert result["transfer_type"] == "Hohmann Transfer"
    assert result["total_delta_v_mps"] == pytest.approx(3000.0)
    assert result["travel_time_days"] == pytest.approx(100.0)
    assert result["travel_time_h_m_s"] == "100 days, 00 hours, 00 minutes, 00 seconds"
    assert result["error"] is None

    # Verify ephemeris was NOT called
    mock_ephemeris_get_body_orbital_state.assert_not_called()

    # Verify orbital mechanics functions were called with semi_major_axis_from_sun from celestial_data
    static_earth_distance = mock_celestial_data_for_planner('Earth')['semi_major_axis_from_sun']
    static_mars_distance = mock_celestial_data_for_planner('Mars')['semi_major_axis_from_sun']

    mock_orbital_mechanics[0].assert_called_once_with(
        planner.mu_sun,
        static_earth_distance,
        static_mars_distance
    )
    mock_orbital_mechanics[1].assert_called_once_with(
        planner.mu_sun,
        static_earth_distance,
        static_mars_distance
    )

def test_plan_trajectory_hohmann_input_validation_empty_departure_with_date(mock_celestial_data_for_planner):
    """
    Test planning with an empty departure body name for Hohmann transfer, with a date.
    """
    planner = TrajectoryPlanner()
    departure_date = datetime(2025, 1, 1)
    result = planner.plan_trajectory("", "Mars", trajectory_type="Hohmann", departure_date=departure_date)
    assert result["success"] is False
    assert "non-empty string" in result["error"]

def test_plan_trajectory_hohmann_input_validation_whitespace_arrival_with_date(mock_celestial_data_for_planner):
    """
    Test planning with a whitespace-only arrival body name for Hohmann transfer, with a date.
    """
    planner = TrajectoryPlanner()
    departure_date = datetime(2025, 1, 1)
    result = planner.plan_trajectory("Earth", "  ", trajectory_type="Hohmann", departure_date=departure_date)
    assert result["success"] is False
    assert "non-empty string" in result["error"]

def test_plan_trajectory_hohmann_input_validation_non_string_departure_with_date(mock_celestial_data_for_planner):
    """
    Test planning with a non-string departure body name for Hohmann transfer, with a date.
    """
    planner = TrajectoryPlanner()
    departure_date = datetime(2025, 1, 1)
    result = planner.plan_trajectory(123, "Mars", trajectory_type="Hohmann", departure_date=departure_date)
    assert result["success"] is False
    assert "non-empty string" in result["error"]

def test_plan_trajectory_hohmann_same_departure_arrival_bodies_with_date(mock_celestial_data_for_planner, mock_ephemeris_get_body_orbital_state):
    """
    Test planning where departure and arrival bodies are the same (case-insensitive) for Hohmann transfer, with a date.
    """
    planner = TrajectoryPlanner()
    departure_date = datetime(2025, 1, 1)
    result = planner.plan_trajectory("Earth", "earth", trajectory_type="Hohmann", departure_date=departure_date)
    assert result["success"] is False
    assert "cannot be the same" in result["error"]
    # Ephemeris will be called for the departure body before the same-body check
    mock_ephemeris_get_body_orbital_state.assert_called_once_with('Earth', departure_date)

def test_plan_trajectory_hohmann_ephemeris_fails_departure_fallback_succeeds(
    mock_celestial_data_for_planner, mock_orbital_mechanics, mock_ephemeris_get_body_orbital_state
):
    """
    Test planning when ephemeris fails for departure body, but static data successfully provides fallback.
    Expected: Success, with warning. Orbital mechanics uses static for departure, ephemeris for arrival.
    """
    planner = TrajectoryPlanner()
    departure_date = datetime(2025, 1, 1)
    # Configure ephemeris mock to return None for Earth, but success for Mars
    mock_ephemeris_get_body_orbital_state.side_effect = lambda name, date: \
        None if name.lower() == 'earth' else {'distance_from_sun_m': 2.30e11} if name.lower() == 'mars' else None

    with patch('builtins.print') as mock_print:
        result = planner.plan_trajectory("Earth", "Mars", trajectory_type="Hohmann", departure_date=departure_date)

    assert result["success"] is True
    assert result["departure_body"] == "Earth"
    assert result["arrival_body"] == "Mars"
    assert result["transfer_type"] == "Hohmann Transfer"
    assert result["total_delta_v_mps"] == pytest.approx(3000.0) # Still expects default mock return
    
    # Check if a warning was printed for Earth
    mock_print.assert_any_call(
        "Warning: Could not retrieve dynamic distance for Earth from ephemeris on 2025-01-01. Falling back to average semi-major axis from static data. This might impact accuracy."
    )

    # Verify ephemeris was called for both, and orbital mechanics used mixed data
    mock_ephemeris_get_body_orbital_state.assert_has_calls([
        call('Earth', departure_date),
        call('Mars', departure_date)
    ], any_order=True)

    static_earth_distance = mock_celestial_data_for_planner('Earth')['semi_major_axis_from_sun']
    ephemeris_mars_distance = mock_ephemeris_get_body_orbital_state.side_effect('Mars', departure_date)['distance_from_sun_m']

    mock_orbital_mechanics[0].assert_called_once_with(
        planner.mu_sun,
        static_earth_distance, # Static for Earth
        ephemeris_mars_distance # Ephemeris for Mars
    )
    mock_orbital_mechanics[1].assert_called_once_with(
        planner.mu_sun,
        static_earth_distance,
        ephemeris_mars_distance
    )

def test_plan_trajectory_hohmann_ephemeris_fails_arrival_fallback_succeeds(
    mock_celestial_data_for_planner, mock_orbital_mechanics, mock_ephemeris_get_body_orbital_state
):
    """
    Test planning when ephemeris fails for arrival body, but static data successfully provides fallback.
    Expected: Success, with warning. Orbital mechanics uses ephemeris for departure, static for arrival.
    """
    planner = TrajectoryPlanner()
    departure_date = datetime(2025, 1, 1)
    # Configure ephemeris mock to return success for Earth, but None for Mars
    mock_ephemeris_get_body_orbital_state.side_effect = lambda name, date: \
        {'distance_from_sun_m': 1.50e11} if name.lower() == 'earth' else None if name.lower() == 'mars' else None

    with patch('builtins.print') as mock_print:
        result = planner.plan_trajectory("Earth", "Mars", trajectory_type="Hohmann", departure_date=departure_date)

    assert result["success"] is True
    # Check if a warning was printed for Mars
    mock_print.assert_any_call(
        "Warning: Could not retrieve dynamic distance for Mars from ephemeris on 2025-01-01. Falling back to average semi-major axis from static data. This might impact accuracy."
    )
    
    # Verify ephemeris was called for both, and orbital mechanics used mixed data
    mock_ephemeris_get_body_orbital_state.assert_has_calls([
        call('Earth', departure_date),
        call('Mars', departure_date)
    ], any_order=True)

    ephemeris_earth_distance = mock_ephemeris_get_body_orbital_state.side_effect('Earth', departure_date)['distance_from_sun_m']
    static_mars_distance = mock_celestial_data_for_planner('Mars')['semi_major_axis_from_sun']

    mock_orbital_mechanics[0].assert_called_once_with(
        planner.mu_sun,
        ephemeris_earth_distance, # Ephemeris for Earth
        static_mars_distance # Static for Mars
    )
    mock_orbital_mechanics[1].assert_called_once_with(
        planner.mu_sun,
        ephemeris_earth_distance,
        static_mars_distance
    )

def test_plan_trajectory_hohmann_ephemeris_and_static_fail_departure(
    mock_celestial_data_for_planner, mock_ephemeris_get_body_orbital_state
):
    """
    Test planning when ephemeris fails for departure, AND static data is also missing 'semi_major_axis_from_sun'.
    Expected: Failure.
    """
    planner = TrajectoryPlanner()
    departure_date = datetime(2025, 1, 1)

    # Make ephemeris fail for Earth, succeed for Mars
    mock_ephemeris_get_body_orbital_state.side_effect = lambda name, date: \
        None if name.lower() == 'earth' else {'distance_from_sun_m': 2.30e11} if name.lower() == 'mars' else None

    # Temporarily configure celestial data mock to be missing semi_major_axis_from_sun for Earth
    with patch('trajectory_planner.celestial_data.get_celestial_body_data') as temp_mock_celestial_data, \
         patch('builtins.print') as mock_print:
        temp_mock_celestial_data.side_effect = lambda name: {
            'sun': {'name': 'Sun', 'gravitational_parameter_mu': 1.327e20},
            'earth': {'name': 'Earth', 'gravitational_parameter_mu': 3.986e14}, # Missing semi_major_axis_from_sun
            'mars': {'name': 'Mars', 'gravitational_parameter_mu': 4.283e13, 'semi_major_axis_from_sun': 2.279e11}
        }.get(name.lower())

        result = planner.plan_trajectory("Earth", "Mars", trajectory_type="Hohmann", departure_date=departure_date)

    assert result["success"] is False
    assert "Could not find sufficient ephemeris data for departure body 'Earth' on 2025-01-01. Falling back to average semi-major axis from static data. This might impact accuracy. Invalid or missing 'semi_major_axis_from_sun' for Earth." in result["error"]
    mock_ephemeris_get_body_orbital_state.assert_called_once_with('Earth', departure_date)
    mock_print.assert_any_call(
        "Warning: Could not retrieve dynamic distance for Earth from ephemeris on 2025-01-01. Falling back to average semi-major axis from static data. This might impact accuracy."
    )


def test_plan_trajectory_hohmann_ephemeris_and_static_fail_arrival(
    mock_celestial_data_for_planner, mock_ephemeris_get_body_orbital_state
):
    """
    Test planning when ephemeris fails for arrival, AND static data is also missing 'semi_major_axis_from_sun'.
    Expected: Failure.
    """
    planner = TrajectoryPlanner()
    departure_date = datetime(2025, 1, 1)

    # Make ephemeris succeed for Earth, fail for Mars
    mock_ephemeris_get_body_orbital_state.side_effect = lambda name, date: \
        {'distance_from_sun_m': 1.50e11} if name.lower() == 'earth' else None if name.lower() == 'mars' else None

    # Temporarily configure celestial data mock to be missing semi_major_axis_from_sun for Mars
    with patch('trajectory_planner.celestial_data.get_celestial_body_data') as temp_mock_celestial_data, \
         patch('builtins.print') as mock_print:
        temp_mock_celestial_data.side_effect = lambda name: {
            'sun': {'name': 'Sun', 'gravitational_parameter_mu': 1.327e20},
            'earth': {'name': 'Earth', 'semi_major_axis_from_sun': 1.496e11},
            'mars': {'name': 'Mars', 'gravitational_parameter_mu': 4.283e13} # Missing semi_major_axis_from_sun
        }.get(name.lower())

        result = planner.plan_trajectory("Earth", "Mars", trajectory_type="Hohmann", departure_date=departure_date)

    assert result["success"] is False
    assert "Could not find sufficient ephemeris data for arrival body 'Mars' on 2025-01-01. Falling back to average semi-major axis from static data. This might impact accuracy. Invalid or missing 'semi_major_axis_from_sun' for Mars." in result["error"]
    mock_ephemeris_get_body_orbital_state.assert_has_calls([
        call('Earth', departure_date),
        call('Mars', departure_date)
    ], any_order=True)
    mock_print.assert_any_call(
        "Warning: Could not retrieve dynamic distance for Mars from ephemeris on 2025-01-01. Falling back to average semi-major axis from static data. This might impact accuracy."
    )


def test_plan_trajectory_hohmann_ephemeris_missing_distance_key_fallback_succeeds(
    mock_celestial_data_for_planner, mock_orbital_mechanics, mock_ephemeris_get_body_orbital_state
):
    """
    Test planning when ephemeris returns data but without 'distance_from_sun_m' for departure,
    but static data provides fallback.
    """
    planner = TrajectoryPlanner()
    departure_date = datetime(2025, 1, 1)
    # Configure ephemeris mock to return partial data for Earth
    mock_ephemeris_get_body_orbital_state.side_effect = lambda name, date: \
        {'some_other_key': 123} if name.lower() == 'earth' else {'distance_from_sun_m': 2.30e11} if name.lower() == 'mars' else None

    with patch('builtins.print') as mock_print:
        result = planner.plan_trajectory("Earth", "Mars", trajectory_type="Hohmann", departure_date=departure_date)

    assert result["success"] is True
    # Check if a warning was printed for Earth
    mock_print.assert_any_call(
        "Warning: Could not retrieve dynamic distance for Earth from ephemeris on 2025-01-01. Falling back to average semi-major axis from static data. This might impact accuracy."
    )
    mock_ephemeris_get_body_orbital_state.assert_has_calls([
        call('Earth', departure_date),
        call('Mars', departure_date)
    ], any_order=True)
    
    static_earth_distance = mock_celestial_data_for_planner('Earth')['semi_major_axis_from_sun']
    ephemeris_mars_distance = mock_ephemeris_get_body_orbital_state.side_effect('Mars', departure_date)['distance_from_sun_m']

    mock_orbital_mechanics[0].assert_called_once_with(planner.mu_sun, static_earth_distance, ephemeris_mars_distance)
    mock_orbital_mechanics[1].assert_called_once_with(planner.mu_sun, static_earth_distance, ephemeris_mars_distance)

def test_plan_trajectory_hohmann_ephemeris_invalid_distance_value_fallback_succeeds(
    mock_celestial_data_for_planner, mock_orbital_mechanics, mock_ephemeris_get_body_orbital_state
):
    """
    Test planning when ephemeris returns invalid 'distance_from_sun_m' (e.g., non-positive) for arrival,
    but static data provides fallback.
    """
    planner = TrajectoryPlanner()
    departure_date = datetime(2025, 1, 1)
    # Configure ephemeris mock to return invalid distance for Mars
    mock_ephemeris_get_body_orbital_state.side_effect = lambda name, date: \
        {'distance_from_sun_m': 1.50e11} if name.lower() == 'earth' else {'distance_from_sun_m': -10.0} if name.lower() == 'mars' else None

    with patch('builtins.print') as mock_print:
        result = planner.plan_trajectory("Earth", "Mars", trajectory_type="Hohmann", departure_date=departure_date)

    assert result["success"] is True
    # Check if a warning was printed for Mars
    mock_print.assert_any_call(
        "Warning: Could not retrieve dynamic distance for Mars from ephemeris on 2025-01-01. Falling back to average semi-major axis from static data. This might impact accuracy."
    )

    mock_ephemeris_get_body_orbital_state.assert_has_calls([
        call('Earth', departure_date),
        call('Mars', departure_date)
    ], any_order=True)

    ephemeris_earth_distance = mock_ephemeris_get_body_orbital_state.side_effect('Earth', departure_date)['distance_from_sun_m']
    static_mars_distance = mock_celestial_data_for_planner('Mars')['semi_major_axis_from_sun']

    mock_orbital_mechanics[0].assert_called_once_with(planner.mu_sun, ephemeris_earth_distance, static_mars_distance)
    mock_orbital_mechanics[1].assert_called_once_with(planner.mu_sun, ephemeris_earth_distance, static_mars_distance)


def test_plan_trajectory_hohmann_non_existent_departure_body_with_date(
    mock_celestial_data_for_planner, mock_ephemeris_get_body_orbital_state
):
    """
    Test planning with a non-existent departure body for Hohmann transfer, with a date.
    Ephemeris will return None, celestial_data will also return None for this body.
    """
    planner = TrajectoryPlanner()
    departure_date = datetime(2025, 1, 1)
    # Ensure ephemeris doesn't have data for 'NonExistentPlanet'
    mock_ephemeris_get_body_orbital_state.side_effect = lambda name, date: \
        {'distance_from_sun_m': 2.30e11} if name.lower() == 'mars' else None

    # Ensure celestial data also doesn't have data for 'NonExistentPlanet' (fixture default)
    
    with patch('builtins.print'): # Suppress warning for this specific test as the overall failure is what we check
        result = planner.plan_trajectory("NonExistentPlanet", "Mars", trajectory_type="Hohmann", departure_date=departure_date)
    assert result["success"] is False
    assert "Could not find data for departure body 'NonExistentPlanet'." in result["error"]
    mock_ephemeris_get_body_orbital_state.assert_called_once_with('NonExistentPlanet', departure_date) # ephemeris tried first

def test_plan_trajectory_hohmann_non_existent_arrival_body_with_date(
    mock_celestial_data_for_planner, mock_ephemeris_get_body_orbital_state
):
    """
    Test planning with a non-existent arrival body for Hohmann transfer, with a date.
    Ephemeris will return None, celestial_data will also return None for this body.
    """
    planner = TrajectoryPlanner()
    departure_date = datetime(2025, 1, 1)
    # Ensure ephemeris doesn't have data for 'NonExistentPlanet'
    mock_ephemeris_get_body_orbital_state.side_effect = lambda name, date: \
        {'distance_from_sun_m': 1.50e11} if name.lower() == 'earth' else None

    # Ensure celestial data also doesn't have data for 'NonExistentPlanet' (fixture default)

    with patch('builtins.print'): # Suppress warning for this specific test
        result = planner.plan_trajectory("Earth", "NonExistentPlanet", trajectory_type="Hohmann", departure_date=departure_date)
    assert result["success"] is False
    assert "Could not find data for arrival body 'NonExistentPlanet'." in result["error"]
    mock_ephemeris_get_body_orbital_state.assert_has_calls([
        call('Earth', departure_date),
        call('NonExistentPlanet', departure_date)
    ], any_order=True)

def test_plan_trajectory_hohmann_orbital_mechanics_value_error_with_date(mock_celestial_data_for_planner, mock_orbital_mechanics, mock_ephemeris_get_body_orbital_state):
    """
    Test planning when orbital_mechanics functions raise a ValueError for Hohmann transfer, with date.
    """
    mock_orbital_mechanics[0].side_effect = ValueError("Test orbital calculation error")
    planner = TrajectoryPlanner()
    departure_date = datetime(2025, 1, 1)
    result = planner.plan_trajectory("Earth", "Mars", trajectory_type="Hohmann", departure_date=departure_date)
    assert result["success"] is False
    assert "Orbital mechanics calculation error: Test orbital calculation error" in result["error"]
    mock_ephemeris_get_body_orbital_state.assert_has_calls([
        call('Earth', departure_date),
        call('Mars', departure_date)
    ], any_order=True)


def test_plan_trajectory_hohmann_orbital_mechanics_unexpected_error_with_date(mock_celestial_data_for_planner, mock_orbital_mechanics, mock_ephemeris_get_body_orbital_state):
    """
    Test planning when orbital_mechanics functions raise an unexpected error for Hohmann transfer, with date.
    """
    mock_orbital_mechanics[1].side_effect = TypeError("Unexpected type error in orbital mechanics")
    planner = TrajectoryPlanner()
    departure_date = datetime(2025, 1, 1)
    result = planner.plan_trajectory("Earth", "Mars", trajectory_type="Hohmann", departure_date=departure_date)
    assert result["success"] is False
    assert "An unexpected error occurred during Hohmann transfer calculation: Unexpected type error in orbital mechanics" in result["error"]
    mock_ephemeris_get_body_orbital_state.assert_has_calls([
        call('Earth', departure_date),
        call('Mars', departure_date)
    ], any_order=True)


def test_plan_trajectory_hohmann_with_fractional_travel_time_with_date(mock_celestial_data_for_planner, mock_orbital_mechanics, mock_ephemeris_get_body_orbital_state):
    """
    Test planning with fractional travel time to ensure correct H:M:S formatting for Hohmann transfer, with date.
    """
    mock_orbital_mechanics[1].return_value = 100 * 24 * 3600 + 3 * 3600 + 30 * 60 + 15 # 100 days, 3h 30m 15s
    planner = TrajectoryPlanner()
    departure_date = datetime(2025, 1, 1)
    result = planner.plan_trajectory("Earth", "Mars", trajectory_type="Hohmann", departure_date=departure_date)

    assert result["success"] is True
    assert result["travel_time_days"] == pytest.approx(100 + (3*3600 + 30*60 + 15) / (24*3600))
    assert result["travel_time_h_m_s"] == "100 days, 03 hours, 30 minutes, 15 seconds"
    mock_ephemeris_get_body_orbital_state.assert_has_calls([
        call('Earth', departure_date),
        call('Mars', departure_date)
    ], any_order=True)


# --- New Tests for TrajectoryPlanner.plan_trajectory specific dispatching and error handling ---

def test_plan_trajectory_default_type_hohmann_with_date(
    mock_celestial_data_for_planner, mock_orbital_mechanics, mock_ephemeris_get_body_orbital_state
):
    """
    Test plan_trajectory with no specified trajectory type, expecting it to default to Hohmann, with a date.
    """
    planner = TrajectoryPlanner()
    departure_date = datetime(2025, 1, 1)
    result = planner.plan_trajectory("Earth", "Mars", departure_date=departure_date) # No trajectory_type specified

    assert result["success"] is True
    assert result["departure_body"] == "Earth"
    assert result["arrival_body"] == "Mars"
    assert result["transfer_type"] == "Hohmann Transfer" # Should still be Hohmann
    assert result["total_delta_v_mps"] == pytest.approx(3000.0)
    assert result["travel_time_days"] == pytest.approx(100.0)
    assert result["travel_time_h_m_s"] == "100 days, 00 hours, 00 minutes, 00 seconds"
    assert result["error"] is None

    # Verify ephemeris was called for both bodies with the correct date
    mock_ephemeris_get_body_orbital_state.assert_has_calls([
        call('Earth', departure_date),
        call('Mars', departure_date)
    ], any_order=True)

    # Verify orbital mechanics functions were called with the distances from ephemeris
    ephemeris_earth_distance = mock_ephemeris_get_body_orbital_state.side_effect('Earth', departure_date)['distance_from_sun_m']
    ephemeris_mars_distance = mock_ephemeris_get_body_orbital_state.side_effect('Mars', departure_date)['distance_from_sun_m']

    mock_orbital_mechanics[0].assert_called_once_with(
        planner.mu_sun,
        ephemeris_earth_distance,
        ephemeris_mars_distance
    )
    mock_orbital_mechanics[1].assert_called_once_with(
        planner.mu_sun,
        ephemeris_earth_distance,
        ephemeris_mars_distance
    )

def test_plan_trajectory_unsupported_type_with_date(mock_celestial_data_for_planner, mock_ephemeris_get_body_orbital_state):
    """
    Test plan_trajectory with an unsupported trajectory type, with a date.
    """
    planner = TrajectoryPlanner()
    departure_date = datetime(2025, 1, 1)
    result = planner.plan_trajectory("Earth", "Mars", trajectory_type="Unsupported", departure_date=departure_date)
    assert result["success"] is False
    assert "Unsupported trajectory type: 'Unsupported'" in result["error"]
    assert "transfer_type" not in result # Ensure it doesn't set an invalid type
    mock_ephemeris_get_body_orbital_state.assert_not_called() # No ephemeris calls for unsupported types

def test_plan_trajectory_invalid_type_non_string_with_date(mock_celestial_data_for_planner, mock_ephemeris_get_body_orbital_state):
    """
    Test plan_trajectory with a non-string trajectory type, with a date.
    """
    planner = TrajectoryPlanner()
    departure_date = datetime(2025, 1, 1)
    result = planner.plan_trajectory("Earth", "Mars", trajectory_type=123, departure_date=departure_date)
    assert result["success"] is False
    assert "Trajectory type must be a non-empty string." in result["error"]
    mock_ephemeris_get_body_orbital_state.assert_not_called() # No ephemeris calls for invalid types

def test_plan_trajectory_direct_transfer_not_implemented(
    mock_celestial_data_for_planner, mock_orbital_mechanics, mock_ephemeris_get_body_orbital_state
):
    """
    Test plan_trajectory with a 'Direct Transfer' trajectory type, which is currently not implemented.
    Expected: Failure with an "Unsupported" error message.
    """
    planner = TrajectoryPlanner()
    departure_date = datetime(2025, 1, 1)
    result = planner.plan_trajectory("Earth", "Mars", trajectory_type="Direct Transfer", departure_date=departure_date)

    assert result["success"] is False
    assert "Unsupported trajectory type: 'Direct Transfer'" in result["error"]
    assert "transfer_type" not in result

    # Ensure no orbital mechanics or ephemeris functions were called as it should fail early
    mock_ephemeris_get_body_orbital_state.assert_not_called()
    mock_orbital_mechanics[0].assert_not_called()
    mock_orbital_mechanics[1].assert_not_called()