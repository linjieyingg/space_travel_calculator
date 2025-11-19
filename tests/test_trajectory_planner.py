import pytest
from unittest.mock import MagicMock, patch, call
from datetime import datetime, timedelta
import math

# Import the module under test and its dependencies for patching
from trajectory_planner import TrajectoryPlanner
import celestial_data
import orbital_mechanics
import ephemeris
import constants

# --- Helper Function ---
def _parse_date_string(date_str: str) -> datetime:
    """Helper to parse date string for tests."""
    return datetime.strptime(date_str, '%Y-%m-%d')

# --- Fixtures ---

@pytest.fixture
def mock_celestial_data_for_planner():
    """
    Fixture to mock celestial_data.get_celestial_body_data for TrajectoryPlanner tests.
    It provides predefined data for Sun, Earth, Mars, and Moon.
    """
    with patch('trajectory_planner.celestial_data.get_celestial_body_data') as mock_get_data:
        mock_get_data.side_effect = lambda name: {
            'sun': {'name': 'Sun', 'mass': 1.989e30, 'radius': 6.957e8, 'gravitational_parameter_mu': 1.327e20, 'semi_major_axis_from_sun': 0, 'orbital_period_seconds': 0},
            'earth': {'name': 'Earth', 'mass': 5.972e24, 'radius': 6.371e6, 'gravitational_parameter_mu': 3.986e14, 'semi_major_axis_from_sun': 1.496e11, 'orbital_period_seconds': 31557600},
            'mars': {'name': 'Mars', 'mass': 6.417e23, 'radius': 3.389e6, 'gravitational_parameter_mu': 4.283e13, 'semi_major_axis_from_sun': 2.279e11, 'orbital_period_seconds': 59354240},
            'moon': {'name': 'Moon', 'mass': 7.342e22, 'radius': 1.737e6, 'gravitational_parameter_mu': 4.905e12, 'semi_major_axis_from_sun': 1.496e11, 'orbital_period_seconds': 2360584},
        }.get(name.lower())
        yield mock_get_data

@pytest.fixture
def mock_orbital_mechanics():
    """
    Fixture to mock functions from the orbital_mechanics module, including a refined
    mock for `calculate_lambert_transfer` to simulate an iterative solver's behavior.
    """
    with patch('trajectory_planner.orbital_mechanics.calculate_hohmann_transfer_delta_v') as mock_hohmann_dv, \
         patch('trajectory_planner.orbital_mechanics.calculate_hohmann_transfer_time_of_flight') as mock_hohmann_tof, \
         patch('trajectory_planner.orbital_mechanics.calculate_lambert_transfer') as mock_lambert_transfer, \
         patch('trajectory_planner.orbital_mechanics.calculate_gravitational_parameter') as mock_grav_param:

        # Default Hohmann mock values
        mock_hohmann_dv.return_value = (1000.0, 500.0, 1500.0)  # dv1, dv2, total_dv
        mock_hohmann_tof.return_value = 100 * 24 * 3600  # 100 days in seconds

        # Mock for the new iterative Lambert solver behavior
        def mock_lambert_solver_effect(mu: float, r1_vec: tuple, r2_vec: tuple, delta_t: float):
            """
            Simulates a Lambert solver that returns different delta_v based on delta_t
            and the relative positions (implicitly assumed by r1_vec and r2_vec).
            It has an 'optimal' time window for lower delta_v.
            """
            if delta_t <= 0:
                raise ValueError("Time of flight must be positive for Lambert transfer.")

            # Assume r1_vec and r2_vec are (r_mag, angle_rad, z_pos) or (x,y,z)
            # For simplicity in mock, we primarily vary total_dv with delta_t.
            # However, for Earth-to-Earth (r1_vec == r2_vec), return 0 delta_v.
            if r1_vec == r2_vec:
                return 0.0, 0.0, 0.0, delta_t

            # Optimal time of flight for a plausible Earth-Mars transfer (approx 200 days)
            optimal_tof_seconds = 17_280_000  # Roughly 200 days
            min_total_dv = 15_000.0  # m/s for a plausible direct transfer

            # Simulate a "U-shaped" delta_v profile around the optimal_tof
            if 16_000_000 <= delta_t <= 19_000_000:  # ~185 to ~220 days
                # Optimal window: dv close to min_total_dv
                deviation = abs(delta_t - optimal_tof_seconds)
                total_dv = min_total_dv + (deviation * 0.0005)  # Small penalty for deviation
                dv1 = total_dv * 0.45
                dv2 = total_dv * 0.55
            elif 10_000_000 < delta_t < 25_000_000:  # ~115 days to ~290 days
                # Wider window, higher dv
                deviation = abs(delta_t - optimal_tof_seconds)
                total_dv = min_total_dv + 5000 + (deviation * 0.001)  # Larger penalty
                dv1 = total_dv * 0.4
                dv2 = total_dv * 0.6
            else:
                # Very short or very long transfers are very high cost or impossible
                if delta_t < 5_000_000 or delta_t > 40_000_000:  # < ~58 days or > ~460 days
                    raise ValueError(f"Lambert solver failed to find solution for delta_t: {delta_t/86400:.1f} days")
                total_dv = min_total_dv + 15000 + (abs(delta_t - optimal_tof_seconds) * 0.002)
                dv1 = total_dv * 0.3
                dv2 = total_dv * 0.7

            return dv1, dv2, total_dv, delta_t # The calculated TOF is usually close to input delta_t

        mock_lambert_transfer.side_effect = mock_lambert_solver_effect
        mock_grav_param.return_value = constants.G_GRAVITATIONAL

        yield {
            'hohmann_dv': mock_hohmann_dv,
            'hohmann_tof': mock_hohmann_tof,
            'lambert_transfer': mock_lambert_transfer,
            'grav_param': mock_grav_param,
        }

@pytest.fixture
def mock_ephemeris_get_heliocentric_state():
    """
    Fixture to mock ephemeris.get_heliocentric_state.
    Provides illustrative heliocentric positions (distance and angular position) for
    Earth and Mars on specific dates.
    """
    with patch('trajectory_planner.ephemeris.get_heliocentric_state') as mock_get_state:
        # Define some illustrative positions. Planner converts these to 3D vectors.
        # We assume a simple 2D conversion: (r*cos(theta), r*sin(theta), 0)
        earth_pos_2024_01_01 = {'distance_from_sun_m': 1.496e11, 'angular_position_rad': 0.0}
        mars_pos_2024_01_01 = {'distance_from_sun_m': 2.279e11, 'angular_position_rad': math.pi} # Opposition

        mock_get_state.side_effect = lambda body_name, date: {
            ('earth', _parse_date_string('2024-01-01')): earth_pos_2024_01_01,
            ('mars', _parse_date_string('2024-01-01')): mars_pos_2024_01_01,
            ('earth', _parse_date_string('2025-01-01')): {'distance_from_sun_m': 1.496e11, 'angular_position_rad': 0.5}, # Example different angle
            ('mars', _parse_date_string('2025-01-01')): {'distance_from_sun_m': 2.279e11, 'angular_position_rad': 0.1}, # Example different angle
        }.get((body_name.lower(), date))
        yield mock_get_state

# --- Test TrajectoryPlanner.__init__ ---

def test_trajectory_planner_init_success(mock_celestial_data_for_planner: MagicMock):
    """
    Tests successful initialization of TrajectoryPlanner, ensuring Sun's gravitational
    parameter is correctly loaded.
    """
    planner = TrajectoryPlanner()
    assert planner.mu_sun == 1.327e20
    mock_celestial_data_for_planner.assert_called_once_with('sun')

def test_trajectory_planner_init_missing_sun_data(mock_celestial_data_for_planner: MagicMock):
    """
    Tests initialization failure when Sun data cannot be retrieved from `celestial_data`.
    """
    mock_celestial_data_for_planner.return_value = None
    with pytest.raises(ValueError, match="Sun data not found"):
        TrajectoryPlanner()

def test_trajectory_planner_init_invalid_sun_mu(mock_celestial_data_for_planner: MagicMock):
    """
    Tests initialization failure when Sun's gravitational parameter is invalid or missing.
    """
    mock_celestial_data_for_planner.return_value = {'name': 'Sun', 'mass': 1.989e30, 'gravitational_parameter_mu': None}
    with pytest.raises(ValueError, match="Sun's gravitational parameter is invalid"):
        TrajectoryPlanner()

# --- Test TrajectoryPlanner.plan_trajectory (Hohmann Type) ---

def test_plan_hohmann_trajectory_success_ephemeris(
    mock_celestial_data_for_planner: MagicMock,
    mock_orbital_mechanics: dict,
    mock_ephemeris_get_heliocentric_state: MagicMock
):
    """
    Tests successful Hohmann transfer planning, verifying that ephemeris data is
    prioritized for calculating radii and that orbital mechanics functions are called.
    """
    planner = TrajectoryPlanner()
    departure_date = _parse_date_string('2024-01-01')

    # Configure ephemeris mock specifically for this test to return distinct values
    earth_pos = {'distance_from_sun_m': 1.496e11, 'angular_position_rad': 0.0}
    mars_pos = {'distance_from_sun_m': 2.279e11, 'angular_position_rad': math.pi}
    mock_ephemeris_get_heliocentric_state.side_effect = lambda body, date: {
        ('earth', departure_date): earth_pos,
        ('mars', departure_date): mars_pos
    }.get((body.lower(), date))

    # Re-initialize planner for fresh mocks and ensure specific calls
    planner = TrajectoryPlanner()
    result = planner.plan_trajectory(
        departure_body_name='Earth',
        arrival_body_name='Mars',
        trajectory_type='Hohmann',
        departure_date=departure_date
    )

    assert result['success'] is True
    assert result['total_delta_v'] == 1500.0  # From mock_hohmann_dv
    assert result['travel_time_seconds'] == 100 * 24 * 3600  # From mock_hohmann_tof
    assert result['trajectory_type'] == 'Hohmann'

    mock_ephemeris_get_heliocentric_state.assert_has_calls([
        call('earth', departure_date),
        call('mars', departure_date)
    ], any_order=True)

    mock_orbital_mechanics['hohmann_dv'].assert_called_once_with(
        planner.mu_sun, earth_pos['distance_from_sun_m'], mars_pos['distance_from_sun_m']
    )
    mock_orbital_mechanics['hohmann_tof'].assert_called_once_with(
        planner.mu_sun, earth_pos['distance_from_sun_m'], mars_pos['distance_from_sun_m']
    )

def test_plan_hohmann_trajectory_missing_ephemeris_data_fallback_to_static(
    mock_celestial_data_for_planner: MagicMock,
    mock_orbital_mechanics: dict,
    mock_ephemeris_get_heliocentric_state: MagicMock
):
    """
    Tests Hohmann transfer planning gracefully falls back to static celestial data
    if ephemeris data retrieval fails.
    """
    mock_ephemeris_get_heliocentric_state.side_effect = ValueError("Ephemeris system offline")

    planner = TrajectoryPlanner()
    departure_date = _parse_date_string('2024-01-01')
    result = planner.plan_trajectory(
        departure_body_name='Earth',
        arrival_body_name='Mars',
        trajectory_type='Hohmann',
        departure_date=departure_date
    )

    assert result['success'] is True
    assert result['total_delta_v'] == 1500.0
    assert result['travel_time_seconds'] == 100 * 24 * 3600
    assert result['trajectory_type'] == 'Hohmann'
    mock_ephemeris_get_heliocentric_state.assert_called() # Should still attempt to call ephemeris
    # Ensure static data was used as fallback
    mock_celestial_data_for_planner.assert_has_calls([call('earth'), call('mars')], any_order=True)

def test_plan_hohmann_trajectory_unsupported_body(mock_celestial_data_for_planner: MagicMock):
    """
    Tests Hohmann planning for a non-existent celestial body.
    """
    planner = TrajectoryPlanner()
    departure_date = _parse_date_string('2024-01-01')
    result = planner.plan_trajectory(
        departure_body_name='Earth',
        arrival_body_name='NonExistentPlanet',
        trajectory_type='Hohmann',
        departure_date=departure_date
    )

    assert result['success'] is False
    assert "Data not found for NonExistentPlanet" in result['message']
    mock_celestial_data_for_planner.assert_called_with('nonexistentplanet')

def test_plan_hohmann_trajectory_orbital_mechanics_error(
    mock_celestial_data_for_planner: MagicMock,
    mock_orbital_mechanics: dict,
    mock_ephemeris_get_heliocentric_state: MagicMock
):
    """
    Tests error handling when underlying orbital mechanics functions for Hohmann transfer fail.
    """
    mock_orbital_mechanics['hohmann_dv'].side_effect = ValueError("Hohmann DV calculation logic error")
    planner = TrajectoryPlanner()
    departure_date = _parse_date_string('2024-01-01')
    result = planner.plan_trajectory(
        departure_body_name='Earth',
        arrival_body_name='Mars',
        trajectory_type='Hohmann',
        departure_date=departure_date
    )

    assert result['success'] is False
    assert "Hohmann DV calculation logic error" in result['message']

# --- Test TrajectoryPlanner.plan_trajectory (Direct Transfer Type - Lambert) ---

def test_plan_direct_transfer_success(
    mock_celestial_data_for_planner: MagicMock,
    mock_orbital_mechanics: dict,
    mock_ephemeris_get_heliocentric_state: MagicMock
):
    """
    Tests successful direct transfer (Lambert) planning, verifying interaction with
    ephemeris data and the iterative search for optimal transfer time via the Lambert solver mock.
    """
    planner = TrajectoryPlanner()
    departure_date = _parse_date_string('2024-01-01')

    # Configure ephemeris mock to return distinct positions for Earth and Mars
    earth_pos = {'distance_from_sun_m': 1.496e11, 'angular_position_rad': 0.0}
    mars_pos = {'distance_from_sun_m': 2.279e11, 'angular_position_rad': math.pi}
    mock_ephemeris_get_heliocentric_state.side_effect = lambda body, date: {
        ('earth', departure_date): earth_pos,
        ('mars', departure_date): mars_pos
    }.get((body.lower(), date))

    # Expected optimal values from the mock_lambert_solver_effect (around 200 days / 17.28M seconds)
    expected_travel_time_seconds = 17_280_000 # Optimal point in the mock
    expected_total_dv = 15_000.0 # Minimum total_dv from the mock

    result = planner.plan_trajectory(
        departure_body_name='Earth',
        arrival_body_name='Mars',
        trajectory_type='Direct',
        departure_date=departure_date
    )

    assert result['success'] is True
    assert math.isclose(result['total_delta_v'], expected_total_dv, rel_tol=1e-3)
    assert math.isclose(result['travel_time_seconds'], expected_travel_time_seconds, rel_tol=1e-3)
    assert result['trajectory_type'] == 'Direct'

    # Verify that calculate_lambert_transfer was called multiple times, indicating iteration
    assert mock_orbital_mechanics['lambert_transfer'].call_count > 1
    # Check that calls were made with relevant arguments (e.g., mu, 3D position vectors, delta_t)
    first_call_args = mock_orbital_mechanics['lambert_transfer'].call_args_list[0].args
    assert first_call_args[0] == planner.mu_sun # mu
    # Planner converts ephemeris (distance, angle) to (x,y,z) 3-tuples
    assert isinstance(first_call_args[1], tuple) and len(first_call_args[1]) == 3 # r1_vec
    assert isinstance(first_call_args[2], tuple) and len(first_call_args[2]) == 3 # r2_vec
    assert isinstance(first_call_args[3], (int, float)) # delta_t

def test_plan_direct_transfer_no_viable_solution_found(
    mock_celestial_data_for_planner: MagicMock,
    mock_orbital_mechanics: dict,
    mock_ephemeris_get_heliocentric_state: MagicMock
):
    """
    Tests direct transfer when the Lambert solver (mock) consistently returns very high delta-v
    or raises errors for most/all probed delta_t, implying no viable optimal solution can be found
    within the planner's search parameters.
    """
    planner = TrajectoryPlanner()
    departure_date = _parse_date_string('2024-01-01')

    # Configure ephemeris mock
    earth_pos = {'distance_from_sun_m': 1.496e11, 'angular_position_rad': 0.0}
    mars_pos = {'distance_from_sun_m': 2.279e11, 'angular_position_rad': 0.1} # Positions making transfer difficult
    mock_ephemeris_get_heliocentric_state.side_effect = lambda body, date: {
        ('earth', departure_date): earth_pos,
        ('mars', departure_date): mars_pos
    }.get((body.lower(), date))

    # Make mock Lambert solver always return a very high delta_v for any valid input,
    # to simulate a situation where no "optimal" solution is found by the planner.
    def consistently_high_dv_solver(mu, r1_vec, r2_vec, delta_t):
        if delta_t <= 0:
            raise ValueError("Time of flight must be positive for Lambert transfer.")
        # Returns a dv that's always considered 'suboptimal' by the planner's thresholds
        return 100_000.0, 100_000.0, 200_000.0, delta_t
    mock_orbital_mechanics['lambert_transfer'].side_effect = consistently_high_dv_solver

    result = planner.plan_trajectory(
        departure_body_name='Earth',
        arrival_body_name='Mars',
        trajectory_type='Direct',
        departure_date=departure_date
    )

    assert result['success'] is False
    assert "Could not find a viable direct transfer trajectory" in result['message']
    assert mock_orbital_mechanics['lambert_transfer'].call_count > 1 # Still tries multiple times

def test_plan_direct_transfer_lambert_solver_error_propagation(
    mock_celestial_data_for_planner: MagicMock,
    mock_orbital_mechanics: dict,
    mock_ephemeris_get_heliocentric_state: MagicMock
):
    """
    Tests direct transfer when the underlying Lambert solver raises a critical error
    for all attempted `delta_t` values, leading to overall plan failure.
    """
    planner = TrajectoryPlanner()
    departure_date = _parse_date_string('2024-01-01')

    # Configure ephemeris mock
    earth_pos = {'distance_from_sun_m': 1.496e11, 'angular_position_rad': 0.0}
    mars_pos = {'distance_from_sun_m': 2.279e11, 'angular_position_rad': math.pi}
    mock_ephemeris_get_heliocentric_state.side_effect = lambda body, date: {
        ('earth', departure_date): earth_pos,
        ('mars', departure_date): mars_pos
    }.get((body.lower(), date))

    # Make mock Lambert solver always raise an error
    mock_orbital_mechanics['lambert_transfer'].side_effect = ValueError("Lambert calculation failed catastrophically due to bad geometry.")

    result = planner.plan_trajectory(
        departure_body_name='Earth',
        arrival_body_name='Mars',
        trajectory_type='Direct',
        departure_date=departure_date
    )

    assert result['success'] is False
    assert "Lambert calculation failed catastrophically due to bad geometry." in result['message']
    assert mock_orbital_mechanics['lambert_transfer'].call_count >= 1 # At least one attempt before error is fatal or caught

def test_plan_direct_transfer_identical_bodies_earth_to_earth(
    mock_celestial_data_for_planner: MagicMock,
    mock_orbital_mechanics: dict,
    mock_ephemeris_get_heliocentric_state: MagicMock
):
    """
    Tests planning a direct transfer from a body to itself (e.g., Earth to Earth).
    The planner should detect this special case and return zero delta-v and travel time.
    """
    planner = TrajectoryPlanner()
    departure_date = _parse_date_string('2024-01-01')

    # Ephemeris will be called for 'Earth' twice, returning identical positions.
    earth_pos = {'distance_from_sun_m': 1.496e11, 'angular_position_rad': 0.0}
    mock_ephemeris_get_heliocentric_state.side_effect = lambda body, date: \
        earth_pos if body.lower() == 'earth' and date == departure_date else None

    result = planner.plan_trajectory(
        departure_body_name='Earth',
        arrival_body_name='Earth', # Plan from Earth to Earth
        trajectory_type='Direct',
        departure_date=departure_date
    )

    assert result['success'] is True
    assert math.isclose(result['total_delta_v'], 0.0, abs_tol=1e-9)
    assert math.isclose(result['travel_time_seconds'], 0.0, abs_tol=1e-9)
    assert result['trajectory_type'] == 'Direct'
    assert "already at destination" in result['message'].lower()

    # The planner should recognize identical bodies before calling ephemeris or orbital mechanics
    mock_ephemeris_get_heliocentric_state.assert_not_called()
    mock_orbital_mechanics['lambert_transfer'].assert_not_called()

# --- Test TrajectoryPlanner.plan_trajectory (Dispatching and Validation) ---

def test_plan_trajectory_unsupported_type(mock_celestial_data_for_planner: MagicMock):
    """
    Tests `plan_trajectory` with an unsupported trajectory type string.
    """
    planner = TrajectoryPlanner()
    result = planner.plan_trajectory(
        departure_body_name='Earth',
        arrival_body_name='Mars',
        trajectory_type='UnsupportedType',
        departure_date=_parse_date_string('2024-01-01')
    )
    assert result['success'] is False
    assert "Unsupported trajectory type: UnsupportedType" in result['message']

def test_plan_trajectory_invalid_date_format(mock_celestial_data_for_planner: MagicMock):
    """
    Tests `plan_trajectory` when `departure_date` is provided in an invalid string format.
    """
    planner = TrajectoryPlanner()
    result = planner.plan_trajectory(
        departure_body_name='Earth',
        arrival_body_name='Mars',
        trajectory_type='Hohmann',
        departure_date='2024/01/01' # Invalid format
    )
    assert result['success'] is False
    assert "Invalid departure date format" in result['message']

def test_plan_trajectory_missing_departure_date_for_dynamic_plan(mock_celestial_data_for_planner: MagicMock):
    """
    Tests `plan_trajectory` when `departure_date` is omitted for a trajectory type
    that dynamically requires it (e.g., 'Direct' or 'Hohmann' using ephemeris).
    """
    planner = TrajectoryPlanner()
    result = planner.plan_trajectory(
        departure_body_name='Earth',
        arrival_body_name='Mars',
        trajectory_type='Direct' # Requires departure_date
    )
    assert result['success'] is False
    assert "departure_date is required for this trajectory type" in result['message']