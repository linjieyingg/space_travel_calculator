import pytest
from unittest.mock import MagicMock, patch
from fuel_optimizer import optimize_fuel_for_trajectory, TrajectoryPlanner
# Import the actual function to check its signature, though not strictly needed for patching
# from propulsion_system import calculate_required_fuel_mass 

# Mock TrajectoryPlanner and propulsion_system for fuel_optimizer tests
@pytest.fixture
def mock_dependencies_for_fuel_optimizer():
    with patch('fuel_optimizer.TrajectoryPlanner') as MockTrajectoryPlanner, \
         patch('fuel_optimizer.calculate_fuel_mass_from_delta_v') as MockCalculateFuelMass: # This is the aliased function
        
        # Configure TrajectoryPlanner mock
        mock_planner_instance = MockTrajectoryPlanner.return_value
        mock_planner_instance.plan_hohmann_trajectory.return_value = {
            'success': True,
            'total_delta_v': 10000.0, # meters/second
            'travel_time_days': 100.0,
            'transfer_type': 'Hohmann',
            'average_hohmann_speed_ms': 500000.0 # m/s (0.5 million m/s)
        }

        # Configure calculate_fuel_mass_from_delta_v mock (aliased from propulsion_system)
        MockCalculateFuelMass.return_value = 50000.0 # kg

        yield {
            "MockTrajectoryPlanner": MockTrajectoryPlanner,
            "MockCalculateFuelMass": MockCalculateFuelMass,
            "mock_planner_instance": mock_planner_instance
        }

def test_optimize_fuel_for_trajectory_success(mock_dependencies_for_fuel_optimizer):
    mocks = mock_dependencies_for_fuel_optimizer
    
    start_body = "Earth"
    end_body = "Mars"
    trajectory_type = "Hohmann"
    spacecraft_dry_mass = 10000.0
    engine_specific_impulse = 3000.0

    fuel_mass = optimize_fuel_for_trajectory(
        start_body=start_body,
        end_body=end_body,
        trajectory_type=trajectory_type,
        spacecraft_dry_mass=spacecraft_dry_mass,
        engine_specific_impulse=engine_specific_impulse
    )

    assert fuel_mass == 50000.0
    mocks["MockTrajectoryPlanner"].assert_called_once()
    mocks["mock_planner_instance"].plan_hohmann_trajectory.assert_called_once_with(
        departure_body_name=start_body,
        arrival_body_name=end_body
    )
    mocks["MockCalculateFuelMass"].assert_called_once_with(
        delta_v=10000.0,
        dry_mass=spacecraft_dry_mass,
        specific_impulse=engine_specific_impulse
    )

def test_optimize_fuel_for_trajectory_unsupported_trajectory_type(mock_dependencies_for_fuel_optimizer):
    mocks = mock_dependencies_for_fuel_optimizer # Access mocks even if the test fails before they are yielded
    with pytest.raises(ValueError, match="Unsupported trajectory type: 'Ballistic'"):
        optimize_fuel_for_trajectory(
            start_body="Earth",
            end_body="Mars",
            trajectory_type="Ballistic",
            spacecraft_dry_mass=10000.0,
            engine_specific_impulse=3000.0
        )
    mocks["MockTrajectoryPlanner"].assert_not_called()
    mocks["MockCalculateFuelMass"].assert_not_called()

@pytest.mark.parametrize("dry_mass, specific_impulse, expected_error", [
    (0, 3000, "Spacecraft dry mass must be a positive numeric value."),
    (-100, 3000, "Spacecraft dry mass must be a positive numeric value."),
    (10000, 0, "Engine specific impulse must be a positive numeric value."),
    (10000, -50, "Engine specific impulse must be a positive numeric value."),
    ("invalid", 3000, "Spacecraft dry mass must be a positive numeric value."),
    (10000, "invalid", "Engine specific impulse must be a positive numeric value."),
])
def test_optimize_fuel_for_trajectory_invalid_inputs(dry_mass, specific_impulse, expected_error):
    with pytest.raises(ValueError, match=expected_error):
        optimize_fuel_for_trajectory(
            start_body="Earth",
            end_body="Mars",
            trajectory_type="Hohmann",
            spacecraft_dry_mass=dry_mass,
            engine_specific_impulse=specific_impulse
        )

def test_optimize_fuel_for_trajectory_planner_returns_unsuccessful(mock_dependencies_for_fuel_optimizer):
    mocks = mock_dependencies_for_fuel_optimizer
    mocks["mock_planner_instance"].plan_hohmann_trajectory.return_value = {
        'success': False,
        'error_message': 'Departure body data missing'
    }

    with pytest.raises(ValueError, match="Trajectory planning failed: Departure body data missing"):
        optimize_fuel_for_trajectory(
            start_body="Earth",
            end_body="Mars",
            trajectory_type="Hohmann",
            spacecraft_dry_mass=10000.0,
            engine_specific_impulse=3000.0
        )
    mocks["MockCalculateFuelMass"].assert_not_called()

def test_optimize_fuel_for_trajectory_planner_returns_invalid_delta_v(mock_dependencies_for_fuel_optimizer):
    mocks = mock_dependencies_for_fuel_optimizer
    mocks["mock_planner_instance"].plan_hohmann_trajectory.return_value = {
        'success': True,
        'total_delta_v': -500.0, # Invalid delta_v
        'transfer_type': 'Hohmann'
    }

    with pytest.raises(ValueError, match="Trajectory planner returned an invalid or negative 'total_delta_v': -500.0"):
        optimize_fuel_for_trajectory(
            start_body="Earth",
            end_body="Mars",
            trajectory_type="Hohmann",
            spacecraft_dry_mass=10000.0,
            engine_specific_impulse=3000.0
        )
    mocks["MockCalculateFuelMass"].assert_not_called()

def test_optimize_fuel_for_trajectory_propulsion_system_raises_error(mock_dependencies_for_fuel_optimizer):
    mocks = mock_dependencies_for_fuel_optimizer
    mocks["MockCalculateFuelMass"].side_effect = ValueError("Engine parameters invalid.")

    with pytest.raises(ValueError, match="Fuel optimization input or calculation error: Engine parameters invalid."):
        optimize_fuel_for_trajectory(
            start_body="Earth",
            end_body="Mars",
            trajectory_type="Hohmann",
            spacecraft_dry_mass=10000.0,
            engine_specific_impulse=3000.0
        )

def test_optimize_fuel_for_trajectory_propulsion_returns_negative_fuel(mock_dependencies_for_fuel_optimizer):
    mocks = mock_dependencies_for_fuel_optimizer
    mocks["MockCalculateFuelMass"].return_value = -100.0

    with pytest.raises(ValueError, match="Propulsion system returned an invalid or negative fuel mass: -100.0"):
        optimize_fuel_for_trajectory(
            start_body="Earth",
            end_body="Mars",
            trajectory_type="Hohmann",
            spacecraft_dry_mass=10000.0,
            engine_specific_impulse=3000.0
        )

def test_optimize_fuel_for_trajectory_planner_raises_runtime_error(mock_dependencies_for_fuel_optimizer):
    mocks = mock_dependencies_for_fuel_optimizer
    mocks["mock_planner_instance"].plan_hohmann_trajectory.side_effect = RuntimeError("Something went wrong internally")

    with pytest.raises(RuntimeError, match="An unexpected error occurred during trajectory planning: Something went wrong internally"):
        optimize_fuel_for_trajectory(
            start_body="Earth",
            end_body="Mars",
            trajectory_type="Hohmann",
            spacecraft_dry_mass=10000.0,
            engine_specific_impulse=3000.0
        )
    mocks["MockCalculateFuelMass"].assert_not_called()

def test_optimize_fuel_for_trajectory_planner_raises_type_error_with_kwargs(mock_dependencies_for_fuel_optimizer):
    mocks = mock_dependencies_for_fuel_optimizer
    # Simulate TrajectoryPlanner.plan_hohmann_trajectory not accepting extra kwargs
    mocks["mock_planner_instance"].plan_hohmann_trajectory.side_effect = TypeError("plan_hohmann_trajectory() got an unexpected keyword argument 'extra_arg'")

    with pytest.raises(ValueError, match="Trajectory planner method signature mismatch"):
        optimize_fuel_for_trajectory(
            start_body="Earth",
            end_body="Mars",
            trajectory_type="Hohmann",
            spacecraft_dry_mass=10000.0,
            engine_specific_impulse=3000.0,
            extra_arg="some_value" # This extra arg should cause a TypeError in the mock
        )
    mocks["MockCalculateFuelMass"].assert_not_called()