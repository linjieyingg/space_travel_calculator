import pytest
from unittest.mock import patch, MagicMock
from datetime import timedelta
import math

# Use relative imports as per instructions
from trajectory_planner import TrajectoryPlanner
import celestial_data
import orbital_mechanics
import constants

# Fixture for mocking celestial_data.get_celestial_body_data
@pytest.fixture
def mock_celestial_data_for_planner():
    with patch('trajectory_planner.celestial_data.get_celestial_body_data') as mock_get_data:
        mock_get_data.side_effect = lambda name: {
            'sun': {'name': 'Sun', 'mass': 1.989e30, 'radius': 6.957e8, 'gravitational_parameter_mu': 1.327e20},
            'earth': {'name': 'Earth', 'mass': 5.972e24, 'radius': 6.371e6, 'gravitational_parameter_mu': 3.986e14, 'semi_major_axis_from_sun': 1.496e11},
            'mars': {'name': 'Mars', 'mass': 6.39e23, 'radius': 3.389e6, 'gravitational_parameter_mu': 4.283e13, 'semi_major_axis_from_sun': 2.279e11},
            'jupiter': {'name': 'Jupiter', 'mass': 1.898e27, 'radius': 6.9911e7, 'gravitational_parameter_mu': 1.267e17, 'semi_major_axis_from_sun': 7.785e11}
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

def test_plan_trajectory_hohmann_success_earth_mars(mock_celestial_data_for_planner, mock_orbital_mechanics):
    """
    Test successful planning of an Earth to Mars Hohmann transfer using plan_trajectory.
    """
    planner = TrajectoryPlanner()
    result = planner.plan_trajectory("Earth", "Mars", trajectory_type="Hohmann")

    assert result["success"] is True
    assert result["departure_body"] == "Earth"
    assert result["arrival_body"] == "Mars"
    assert result["transfer_type"] == "Hohmann Transfer"
    assert result["total_delta_v_mps"] == pytest.approx(3000.0)
    assert result["travel_time_days"] == pytest.approx(100.0)
    assert result["travel_time_h_m_s"] == "100 days, 00 hours, 00 minutes, 00 seconds"
    assert result["error"] is None

    mock_orbital_mechanics[0].assert_called_once_with(
        planner.mu_sun,
        mock_celestial_data_for_planner('Earth')['semi_major_axis_from_sun'],
        mock_celestial_data_for_planner('Mars')['semi_major_axis_from_sun']
    )
    mock_orbital_mechanics[1].assert_called_once_with(
        planner.mu_sun,
        mock_celestial_data_for_planner('Earth')['semi_major_axis_from_sun'],
        mock_celestial_data_for_planner('Mars')['semi_major_axis_from_sun']
    )

def test_plan_trajectory_hohmann_input_validation_empty_departure(mock_celestial_data_for_planner):
    """
    Test planning with an empty departure body name for Hohmann transfer.
    """
    planner = TrajectoryPlanner()
    result = planner.plan_trajectory("", "Mars", trajectory_type="Hohmann")
    assert result["success"] is False
    assert "non-empty string" in result["error"]

def test_plan_trajectory_hohmann_input_validation_whitespace_arrival(mock_celestial_data_for_planner):
    """
    Test planning with a whitespace-only arrival body name for Hohmann transfer.
    """
    planner = TrajectoryPlanner()
    result = planner.plan_trajectory("Earth", "  ", trajectory_type="Hohmann")
    assert result["success"] is False
    assert "non-empty string" in result["error"]

def test_plan_trajectory_hohmann_input_validation_non_string_departure(mock_celestial_data_for_planner):
    """
    Test planning with a non-string departure body name for Hohmann transfer.
    """
    planner = TrajectoryPlanner()
    result = planner.plan_trajectory(123, "Mars", trajectory_type="Hohmann")
    assert result["success"] is False
    assert "non-empty string" in result["error"]

def test_plan_trajectory_hohmann_same_departure_arrival_bodies(mock_celestial_data_for_planner):
    """
    Test planning where departure and arrival bodies are the same (case-insensitive) for Hohmann transfer.
    """
    planner = TrajectoryPlanner()
    result = planner.plan_trajectory("Earth", "earth", trajectory_type="Hohmann")
    assert result["success"] is False
    assert "cannot be the same" in result["error"]

def test_plan_trajectory_hohmann_non_existent_departure_body(mock_celestial_data_for_planner):
    """
    Test planning with a non-existent departure body for Hohmann transfer.
    """
    planner = TrajectoryPlanner()
    result = planner.plan_trajectory("NonExistentPlanet", "Mars", trajectory_type="Hohmann")
    assert result["success"] is False
    assert "Could not find data for departure body" in result["error"]

def test_plan_trajectory_hohmann_non_existent_arrival_body(mock_celestial_data_for_planner):
    """
    Test planning with a non-existent arrival body for Hohmann transfer.
    """
    planner = TrajectoryPlanner()
    result = planner.plan_trajectory("Earth", "NonExistentPlanet", trajectory_type="Hohmann")
    assert result["success"] is False
    assert "Could not find data for arrival body" in result["error"]

def test_plan_trajectory_hohmann_missing_semi_major_axis_departure(mock_celestial_data_for_planner):
    """
    Test planning when 'semi_major_axis_from_sun' is missing for the departure body for Hohmann transfer.
    """
    mock_celestial_data_for_planner.side_effect = lambda name: {
        'sun': {'name': 'Sun', 'gravitational_parameter_mu': 1.327e20},
        'earth': {'name': 'Earth', 'gravitational_parameter_mu': 3.986e14}, # Missing semi_major_axis_from_sun
        'mars': {'name': 'Mars', 'gravitational_parameter_mu': 4.283e13, 'semi_major_axis_from_sun': 2.279e11}
    }.get(name.lower())
    planner = TrajectoryPlanner()
    result = planner.plan_trajectory("Earth", "Mars", trajectory_type="Hohmann")
    assert result["success"] is False
    assert "Invalid or missing 'semi_major_axis_from_sun' for Earth" in result["error"]

def test_plan_trajectory_hohmann_invalid_semi_major_axis_arrival(mock_celestial_data_for_planner):
    """
    Test planning when 'semi_major_axis_from_sun' is invalid for the arrival body for Hohmann transfer.
    """
    mock_celestial_data_for_planner.side_effect = lambda name: {
        'sun': {'name': 'Sun', 'gravitational_parameter_mu': 1.327e20},
        'earth': {'name': 'Earth', 'semi_major_axis_from_sun': 1.496e11},
        'mars': {'name': 'Mars', 'semi_major_axis_from_sun': -1.0} # Invalid semi_major_axis_from_sun
    }.get(name.lower())
    planner = TrajectoryPlanner()
    result = planner.plan_trajectory("Earth", "Mars", trajectory_type="Hohmann")
    assert result["success"] is False
    assert "Invalid or missing 'semi_major_axis_from_sun' for Mars" in result["error"]

def test_plan_trajectory_hohmann_orbital_mechanics_value_error(mock_celestial_data_for_planner, mock_orbital_mechanics):
    """
    Test planning when orbital_mechanics functions raise a ValueError for Hohmann transfer.
    """
    mock_orbital_mechanics[0].side_effect = ValueError("Test orbital calculation error")
    planner = TrajectoryPlanner()
    result = planner.plan_trajectory("Earth", "Mars", trajectory_type="Hohmann")
    assert result["success"] is False
    assert "Orbital mechanics calculation error: Test orbital calculation error" in result["error"]

def test_plan_trajectory_hohmann_orbital_mechanics_unexpected_error(mock_celestial_data_for_planner, mock_orbital_mechanics):
    """
    Test planning when orbital_mechanics functions raise an unexpected error for Hohmann transfer.
    """
    mock_orbital_mechanics[1].side_effect = TypeError("Unexpected type error in orbital mechanics")
    planner = TrajectoryPlanner()
    result = planner.plan_trajectory("Earth", "Mars", trajectory_type="Hohmann")
    assert result["success"] is False
    assert "An unexpected error occurred during Hohmann transfer calculation: Unexpected type error in orbital mechanics" in result["error"]

def test_plan_trajectory_hohmann_with_fractional_travel_time(mock_celestial_data_for_planner, mock_orbital_mechanics):
    """
    Test planning with fractional travel time to ensure correct H:M:S formatting for Hohmann transfer.
    """
    mock_orbital_mechanics[1].return_value = 100 * 24 * 3600 + 3 * 3600 + 30 * 60 + 15 # 100 days, 3h 30m 15s
    planner = TrajectoryPlanner()
    result = planner.plan_trajectory("Earth", "Mars", trajectory_type="Hohmann")

    assert result["success"] is True
    assert result["travel_time_days"] == pytest.approx(100 + (3*3600 + 30*60 + 15) / (24*3600))
    assert result["travel_time_h_m_s"] == "100 days, 03 hours, 30 minutes, 15 seconds"


# --- New Tests for TrajectoryPlanner.plan_trajectory specific dispatching and error handling ---

def test_plan_trajectory_default_type_hohmann(mock_celestial_data_for_planner, mock_orbital_mechanics):
    """
    Test plan_trajectory with no specified trajectory type, expecting it to default to Hohmann.
    """
    planner = TrajectoryPlanner()
    result = planner.plan_trajectory("Earth", "Mars") # No trajectory_type specified

    assert result["success"] is True
    assert result["departure_body"] == "Earth"
    assert result["arrival_body"] == "Mars"
    assert result["transfer_type"] == "Hohmann Transfer" # Should still be Hohmann
    assert result["total_delta_v_mps"] == pytest.approx(3000.0)
    assert result["travel_time_days"] == pytest.approx(100.0)
    assert result["travel_time_h_m_s"] == "100 days, 00 hours, 00 minutes, 00 seconds"
    assert result["error"] is None

    # Verify orbital mechanics functions were called
    mock_orbital_mechanics[0].assert_called_once_with(
        planner.mu_sun,
        mock_celestial_data_for_planner('Earth')['semi_major_axis_from_sun'],
        mock_celestial_data_for_planner('Mars')['semi_major_axis_from_sun']
    )
    mock_orbital_mechanics[1].assert_called_once_with(
        planner.mu_sun,
        mock_celestial_data_for_planner('Earth')['semi_major_axis_from_sun'],
        mock_celestial_data_for_planner('Mars')['semi_major_axis_from_sun']
    )

def test_plan_trajectory_unsupported_type(mock_celestial_data_for_planner):
    """
    Test plan_trajectory with an unsupported trajectory type.
    """
    planner = TrajectoryPlanner()
    result = planner.plan_trajectory("Earth", "Mars", trajectory_type="Unsupported")
    assert result["success"] is False
    assert "Unsupported trajectory type: 'Unsupported'" in result["error"]
    assert "transfer_type" not in result # Ensure it doesn't set an invalid type

def test_plan_trajectory_invalid_type_non_string(mock_celestial_data_for_planner):
    """
    Test plan_trajectory with a non-string trajectory type.
    """
    planner = TrajectoryPlanner()
    result = planner.plan_trajectory("Earth", "Mars", trajectory_type=123)
    assert result["success"] is False
    assert "Trajectory type must be a non-empty string." in result["error"]