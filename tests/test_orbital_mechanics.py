import pytest
import math
from datetime import datetime, timedelta
from unittest.mock import patch, MagicMock

# Assuming these imports are correct relative to the test file's location
# For example, if tests/test_orbital_mechanics.py, then `from .. import orbital_mechanics`
# For this structure, direct import from the root will be assumed, or `from . import ...` might be needed
# if `orbital_mechanics.py` is in the same directory as `tests/`.
# Given the `REPOSITORY CONTEXT` it means `orbital_mechanics.py` is at the root level.
# So, `from orbital_mechanics import ...` is appropriate if the project root is in PYTHONPATH.
# For pytest, this usually means running from the root directory.
from orbital_mechanics import (
    vector_magnitude, vector_subtract, vector_add, vector_scale, vector_dot, vector_cross,
    calculate_gravitational_parameter, calculate_circular_orbital_velocity,
    calculate_hohmann_transfer_delta_v, calculate_hohmann_transfer_time_of_flight,
    calculate_elliptical_orbital_period, calculate_escape_velocity,
    calculate_lambert_transfer,
    _get_body_orbital_velocity_vector,
)
# We need to import celestial_data to patch it, not necessarily to use its functions directly
# If celestial_data.py is in the same directory as orbital_mechanics.py, then it's a sibling import.
# The original analysis said 'local custom module' which indicates sibling or sub-package.
# For simplicity, assuming direct import from root.
import celestial_data
import constants

# --- Fixtures ---

@pytest.fixture
def mock_celestial_data():
    """
    Fixture to mock celestial_data.get_celestial_body_data.
    Provides data for Sun, Earth, Mars, Moon for tests.
    """
    mock_data = {
        'sun': {
            'mass': 1.989e30,
            'radius': 6.957e8,
            'gravitational_parameter_mu': 1.32712440018e20,  # m^3/s^2
            'semi_major_axis_from_sun': 0,  # N/A for Sun itself
            'orbital_period': 0,  # N/A
        },
        'earth': {
            'mass': 5.972e24,
            'radius': 6.371e6,
            'gravitational_parameter_mu': 3.986004418e14,
            'semi_major_axis_from_sun': 1.496e11,  # 1 AU
            'orbital_period': 3.15576e7  # 1 sidereal year in seconds
        },
        'mars': {
            'mass': 6.417e23,
            'radius': 3.3895e6,
            'gravitational_parameter_mu': 4.282837e13,
            'semi_major_axis_from_sun': 2.279e11,  # 1.52 AU
            'orbital_period': 5.93542944e7  # ~687 Earth days in seconds
        },
        'moon': {
            'mass': 7.342e22,
            'radius': 1.7374e6,
            'gravitational_parameter_mu': 4.9048695e12,
            'semi_major_axis_from_sun': 1.496e11,  # Approx same as Earth for heliocentric
            'orbital_period': 2.36058e6  # ~27.3 Earth days in seconds (sidereal)
        }
    }

    # Patch celestial_data.get_celestial_body_data in the orbital_mechanics module's context
    with patch('orbital_mechanics.celestial_data.get_celestial_body_data') as mock_get_data:
        mock_get_data.side_effect = lambda name: mock_data.get(name.lower())
        yield mock_get_data

@pytest.fixture
def mock_ephemeris():
    """
    Fixture to mock ephemeris.get_heliocentric_state.
    Provides 3-element position and velocity tuples for Earth and Mars at a specific date.
    This fixture assumes that ephemeris.get_heliocentric_state will return a dictionary
    containing 'position_vector_m' and 'velocity_vector_mps'.
    """
    # Using Sun's mu from mock_celestial_data for consistency in calculating these hardcoded values
    mu_sun = 1.32712440018e20 # m^3/s^2

    # Earth's data for 2023-01-01 (simplified circular, at +X axis, theta = 0)
    r_earth = 1.496e11 # m
    v_mag_earth = math.sqrt(mu_sun / r_earth) # m/s
    # Position vector (r*cos(theta), r*sin(theta), 0)
    pos_earth_vec = (r_earth, 0.0, 0.0)
    # Velocity vector (-v_mag*sin(theta), v_mag*cos(theta), 0)
    vel_earth_vec = (0.0, v_mag_earth, 0.0) # ~29784.85 m/s

    # Mars' data for 2023-01-01 (simplified circular, at -X axis, theta = pi)
    r_mars = 2.279e11 # m
    v_mag_mars = math.sqrt(mu_sun / r_mars) # m/s
    # Position vector (r*cos(theta), r*sin(theta), 0)
    pos_mars_vec = (-r_mars, 0.0, 0.0)
    # Velocity vector (-v_mag*sin(theta), v_mag*cos(theta), 0)
    vel_mars_vec = (0.0, -v_mag_mars, 0.0) # ~-24128.53 m/s

    mock_states = {
        ('earth', datetime(2023, 1, 1, 0, 0, 0, tzinfo=None)): {
            'position_vector_m': pos_earth_vec,
            'velocity_vector_mps': vel_earth_vec,
        },
        ('mars', datetime(2023, 1, 1, 0, 0, 0, tzinfo=None)): {
            'position_vector_m': pos_mars_vec,
            'velocity_vector_mps': vel_mars_vec,
        }
    }

    # Patch ephemeris.get_heliocentric_state in the orbital_mechanics module's context
    with patch('orbital_mechanics.ephemeris.get_heliocentric_state') as mock_get_state:
        def side_effect(body_name, date):
            # Normalizing date to be naive for lookup
            naive_date = date.replace(tzinfo=None) if date.tzinfo else date
            return mock_states.get((body_name.lower(), naive_date))
        mock_get_state.side_effect = side_effect
        yield mock_get_state

# --- Vector Operations Tests ---

def test_vector_magnitude_valid_input():
    """Tests vector_magnitude with valid 3D vectors including zero and non-zero components."""
    assert vector_magnitude((3, 4, 0)) == pytest.approx(5.0)
    assert vector_magnitude((1, 1, 1)) == pytest.approx(math.sqrt(3))
    assert vector_magnitude((0, 0, 0)) == pytest.approx(0.0)

def test_vector_magnitude_invalid_input_type():
    """Tests vector_magnitude with invalid input types or incorrect dimensions."""
    with pytest.raises(ValueError, match="Vector must be a tuple or list of 3 numeric elements"):
        vector_magnitude([1, 2])  # Incorrect length
    with pytest.raises(ValueError, match="Vector must be a tuple or list of 3 numeric elements"):
        vector_magnitude((1, 2, 'a'))  # Non-numeric element
    with pytest.raises(ValueError, match="Vector must be a tuple or list of 3 numeric elements"):
        vector_magnitude(None)  # Incorrect type

def test_vector_subtract_valid_inputs():
    """Tests vector_subtract with valid 3D vectors."""
    assert vector_subtract((5, 5, 5), (1, 2, 3)) == (4, 3, 2)
    assert vector_subtract((1, 1, 1), (1, 1, 1)) == (0, 0, 0)
    assert vector_subtract((0, 0, 0), (1, 2, 3)) == (-1, -2, -3)

def test_vector_subtract_invalid_inputs():
    """Tests vector_subtract with invalid vector inputs (type, length, or content)."""
    with pytest.raises(ValueError, match="Vector must be a tuple or list of 3 numeric elements"):
        vector_subtract((1, 2), (3, 4, 5))
    with pytest.raises(ValueError, match="Vector must be a tuple or list of 3 numeric elements"):
        vector_subtract((1, 2, 3), (3, 4))
    with pytest.raises(ValueError, match="Vector must be a tuple or list of 3 numeric elements"):
        vector_subtract((1, 2, 'a'), (3, 4, 5))

def test_vector_add_valid_inputs():
    """Tests vector_add with valid 3D vectors."""
    assert vector_add((1, 2, 3), (4, 5, 6)) == (5, 7, 9)
    assert vector_add((0, 0, 0), (1, 2, 3)) == (1, 2, 3)

def test_vector_add_invalid_inputs():
    """Tests vector_add with invalid vector inputs."""
    with pytest.raises(ValueError, match="Vector must be a tuple or list of 3 numeric elements"):
        vector_add((1, 2), (3, 4, 5))

def test_vector_scale_valid_inputs():
    """Tests vector_scale with valid 3D vector and scalar."""
    assert vector_scale((1, 2, 3), 2) == (2, 4, 6)
    assert vector_scale((1, 2, 3), 0) == (0, 0, 0)
    assert vector_scale((1, 2, 3), -1) == (-1, -2, -3)
    assert vector_scale((0, 0, 0), 5) == (0, 0, 0)

def test_vector_scale_invalid_vector():
    """Tests vector_scale with invalid vector input."""
    with pytest.raises(ValueError, match="Vector must be a tuple or list of 3 numeric elements"):
        vector_scale((1, 2), 2)
    with pytest.raises(ValueError, match="Vector must be a tuple or list of 3 numeric elements"):
        vector_scale((1, 2, 'a'), 2)

def test_vector_scale_invalid_scalar():
    """Tests vector_scale with invalid scalar input."""
    with pytest.raises(ValueError, match="Scalar must be a numeric value"):
        vector_scale((1, 2, 3), 'a')

def test_vector_dot_valid_inputs():
    """Tests vector_dot with valid 3D vectors."""
    assert vector_dot((1, 0, 0), (0, 1, 0)) == pytest.approx(0.0)
    assert vector_dot((1, 2, 3), (4, 5, 6)) == pytest.approx(4 + 10 + 18)  # 32.0
    assert vector_dot((1, 1, 1), (1, 1, 1)) == pytest.approx(3.0)

def test_vector_dot_invalid_inputs():
    """Tests vector_dot with invalid vector inputs."""
    with pytest.raises(ValueError, match="Vector must be a tuple or list of 3 numeric elements"):
        vector_dot((1, 2), (3, 4, 5))

def test_vector_cross_valid_inputs():
    """Tests vector_cross with valid 3D vectors, including parallel vectors."""
    assert vector_cross((1, 0, 0), (0, 1, 0)) == (0, 0, 1)
    assert vector_cross((0, 1, 0), (1, 0, 0)) == (0, 0, -1)
    assert vector_cross((1, 2, 3), (4, 5, 6)) == (-3, 6, -3)
    assert vector_cross((1, 1, 1), (1, 1, 1)) == (0, 0, 0)  # Parallel vectors

def test_vector_cross_invalid_inputs():
    """Tests vector_cross with invalid vector inputs."""
    with pytest.raises(ValueError, match="Vector must be a tuple or list of 3 numeric elements"):
        vector_cross((1, 2), (3, 4, 5))

# --- Gravitational Parameter Test ---

def test_calculate_gravitational_parameter_valid(mock_celestial_data):
    """Tests calculate_gravitational_parameter with a valid celestial body name."""
    mu = calculate_gravitational_parameter('Earth')
    assert mu == pytest.approx(mock_celestial_data.side_effect('Earth')['gravitational_parameter_mu'])

def test_calculate_gravitational_parameter_non_existent(mock_celestial_data):
    """Tests calculate_gravitational_parameter with a non-existent celestial body."""
    with pytest.raises(ValueError, match="Celestial body data not found for NonExistentBody"):
        calculate_gravitational_parameter('NonExistentBody')

def test_calculate_gravitational_parameter_missing_data(mock_celestial_data):
    """
    Tests calculate_gravitational_parameter when the celestial data is missing the
    'gravitational_parameter_mu' key.
    """
    # Configure mock to return data without 'gravitational_parameter_mu' for a specific body
    mock_celestial_data.side_effect = lambda name: {'mass': 1, 'radius': 1} if name.lower() == 'badbody' else None
    with pytest.raises(ValueError, match="Gravitational parameter 'gravitational_parameter_mu' missing"):
        calculate_gravitational_parameter('BadBody')

def test_calculate_gravitational_parameter_invalid_mu_value(mock_celestial_data):
    """
    Tests calculate_gravitational_parameter when the retrieved gravitational parameter
    is not a positive value.
    """
    mock_celestial_data.side_effect = lambda name: {'gravitational_parameter_mu': -1} if name.lower() == 'badbody' else None
    with pytest.raises(ValueError, match="Gravitational parameter must be positive"):
        calculate_gravitational_parameter('BadBody')

# --- Circular Orbital Velocity Test ---

def test_calculate_circular_orbital_velocity_valid():
    """Tests calculate_circular_orbital_velocity with valid inputs (e.g., LEO)."""
    mu = 3.986004418e14  # Earth's standard gravitational parameter
    radius = 6.371e6 + 400e3  # Earth's radius + 400km LEO altitude
    expected_v = math.sqrt(mu / radius)
    assert calculate_circular_orbital_velocity(mu, radius) == pytest.approx(expected_v)
    assert calculate_circular_orbital_velocity(mu, radius) > 0

def test_calculate_circular_orbital_velocity_invalid_mu():
    """Tests calculate_circular_orbital_velocity with non-positive gravitational parameter."""
    with pytest.raises(ValueError, match="Gravitational parameter \(mu\) must be positive"):
        calculate_circular_orbital_velocity(0, 7e6)
    with pytest.raises(ValueError, match="Gravitational parameter \(mu\) must be positive"):
        calculate_circular_orbital_velocity(-1, 7e6)

def test_calculate_circular_orbital_velocity_invalid_radius():
    """Tests calculate_circular_orbital_velocity with non-positive radius."""
    with pytest.raises(ValueError, match="Radius must be positive"):
        calculate_circular_orbital_velocity(3.98e14, 0)
    with pytest.raises(ValueError, match="Radius must be positive"):
        calculate_circular_orbital_velocity(3.98e14, -100)

# --- Hohmann Transfer Tests ---

def test_calculate_hohmann_transfer_delta_v_valid(mock_celestial_data):
    """
    Tests calculate_hohmann_transfer_delta_v for a representative Earth-Mars transfer.
    Calculates expected values using known formulas for verification.
    """
    mu_sun = mock_celestial_data.side_effect('sun')['gravitational_parameter_mu']
    r1 = mock_celestial_data.side_effect('earth')['semi_major_axis_from_sun']
    r2 = mock_celestial_data.side_effect('mars')['semi_major_axis_from_sun']

    v_circ1 = math.sqrt(mu_sun / r1)
    v_circ2 = math.sqrt(mu_sun / r2)

    a_transfer = (r1 + r2) / 2
    v_periapsis_transfer = math.sqrt(mu_sun * ((2 / r1) - (1 / a_transfer)))
    v_apoapsis_transfer = math.sqrt(mu_sun * ((2 / r2) - (1 / a_transfer)))

    delta_v1_expected = abs(v_periapsis_transfer - v_circ1)
    delta_v2_expected = abs(v_circ2 - v_apoapsis_transfer)
    total_delta_v_expected = delta_v1_expected + delta_v2_expected

    delta_v1, delta_v2, total_delta_v = calculate_hohmann_transfer_delta_v(mu_sun, r1, r2)

    assert delta_v1 == pytest.approx(delta_v1_expected)
    assert delta_v2 == pytest.approx(delta_v2_expected)
    assert total_delta_v == pytest.approx(total_delta_v_expected)
    assert total_delta_v > 0

def test_calculate_hohmann_transfer_delta_v_invalid_inputs():
    """Tests calculate_hohmann_transfer_delta_v with various invalid input values."""
    mu_sun = 1.327e20
    r1 = 1.496e11
    r2 = 2.279e11

    with pytest.raises(ValueError, match="Gravitational parameter \(mu\) must be positive"):
        calculate_hohmann_transfer_delta_v(0, r1, r2)
    with pytest.raises(ValueError, match="Radii \(r1, r2\) must be positive"):
        calculate_hohmann_transfer_delta_v(mu_sun, 0, r2)
    with pytest.raises(ValueError, match="Radii \(r1, r2\) must be positive"):
        calculate_hohmann_transfer_delta_v(mu_sun, r1, -r2)
    with pytest.raises(ValueError, match="Radii \(r1, r2\) cannot be equal for a Hohmann transfer"):
        calculate_hohmann_transfer_delta_v(mu_sun, r1, r1)

def test_calculate_hohmann_transfer_time_of_flight_valid(mock_celestial_data):
    """
    Tests calculate_hohmann_transfer_time_of_flight for a representative Earth-Mars transfer.
    Verifies against Kepler's Third Law for half an elliptical period.
    """
    mu_sun = mock_celestial_data.side_effect('sun')['gravitational_parameter_mu']
    r1 = mock_celestial_data.side_effect('earth')['semi_major_axis_from_sun']
    r2 = mock_celestial_data.side_effect('mars')['semi_major_axis_from_sun']

    a_transfer = (r1 + r2) / 2
    tof_expected = math.pi * math.sqrt(a_transfer**3 / mu_sun)

    tof = calculate_hohmann_transfer_time_of_flight(mu_sun, r1, r2)
    assert tof == pytest.approx(tof_expected)
    assert tof > 0

def test_calculate_hohmann_transfer_time_of_flight_invalid_inputs():
    """Tests calculate_hohmann_transfer_time_of_flight with various invalid input values."""
    mu_sun = 1.327e20
    r1 = 1.496e11
    r2 = 2.279e11

    with pytest.raises(ValueError, match="Gravitational parameter \(mu\) must be positive"):
        calculate_hohmann_transfer_time_of_flight(0, r1, r2)
    with pytest.raises(ValueError, match="Radii \(r1, r2\) must be positive"):
        calculate_hohmann_transfer_time_of_flight(mu_sun, 0, r2)
    with pytest.raises(ValueError, match="Radii \(r1, r2\) must be positive"):
        calculate_hohmann_transfer_time_of_flight(mu_sun, r1, -r2)
    with pytest.raises(ValueError, match="Radii \(r1, r2\) cannot be equal for a Hohmann transfer"):
        calculate_hohmann_transfer_time_of_flight(mu_sun, r1, r1)

# --- Elliptical Orbital Period Test ---

def test_calculate_elliptical_orbital_period_valid(mock_celestial_data):
    """
    Tests calculate_elliptical_orbital_period for Earth's orbit around the Sun.
    Verifies against Kepler's Third Law.
    """
    mu_sun = mock_celestial_data.side_effect('sun')['gravitational_parameter_mu']
    semi_major_axis_earth = mock_celestial_data.side_effect('earth')['semi_major_axis_from_sun']

    period_expected = 2 * math.pi * math.sqrt(semi_major_axis_earth**3 / mu_sun)
    assert calculate_elliptical_orbital_period(mu_sun, semi_major_axis_earth) == pytest.approx(period_expected)
    assert calculate_elliptical_orbital_period(mu_sun, semi_major_axis_earth) > 0

def test_calculate_elliptical_orbital_period_invalid_inputs():
    """Tests calculate_elliptical_orbital_period with non-positive inputs."""
    with pytest.raises(ValueError, match="Gravitational parameter \(mu\) must be positive"):
        calculate_elliptical_orbital_period(0, 1e11)
    with pytest.raises(ValueError, match="Gravitational parameter \(mu\) must be positive"):
        calculate_elliptical_orbital_period(-1, 1e11)
    with pytest.raises(ValueError, match="Semi-major axis must be positive"):
        calculate_elliptical_orbital_period(1e20, 0)
    with pytest.raises(ValueError, match="Semi-major axis must be positive"):
        calculate_elliptical_orbital_period(1e20, -1e11)

# --- Escape Velocity Test ---

def test_calculate_escape_velocity_valid(mock_celestial_data):
    """Tests calculate_escape_velocity for Earth's surface."""
    mu_earth = mock_celestial_data.side_effect('earth')['gravitational_parameter_mu']
    radius_earth = mock_celestial_data.side_effect('earth')['radius']

    escape_v_expected = math.sqrt(2 * mu_earth / radius_earth)
    assert calculate_escape_velocity(mu_earth, radius_earth) == pytest.approx(escape_v_expected)
    assert calculate_escape_velocity(mu_earth, radius_earth) > 0

def test_calculate_escape_velocity_invalid_inputs():
    """Tests calculate_escape_velocity with non-positive inputs."""
    with pytest.raises(ValueError, match="Gravitational parameter \(mu\) must be positive"):
        calculate_escape_velocity(0, 6.371e6)
    with pytest.raises(ValueError, match="Gravitational parameter \(mu\) must be positive"):
        calculate_escape_velocity(-1, 6.371e6)
    with pytest.raises(ValueError, match="Radius must be positive"):
        calculate_escape_velocity(3.98e14, 0)
    with pytest.raises(ValueError, match="Radius must be positive"):
        calculate_escape_velocity(3.98e14, -1e6)

# --- _get_body_orbital_velocity_vector Test ---

def test_get_body_orbital_velocity_vector_valid(mock_celestial_data, mock_ephemeris):
    """
    Tests _get_body_orbital_velocity_vector with mocked ephemeris data.
    Verifies that the function correctly returns the velocity vector provided by the mock.
    This assumes _get_body_orbital_velocity_vector is updated to extract 'velocity_vector_mps' directly.
    """
    body_name = 'Earth'
    test_date = datetime(2023, 1, 1, 0, 0, 0, tzinfo=None)
    mu_sun = mock_celestial_data.side_effect('sun')['gravitational_parameter_mu'] # Mu is still a required input

    # Get the expected velocity vector from the mock's internal state
    ephemeris_state = mock_ephemeris.side_effect(body_name, test_date)
    expected_velocity_vector = ephemeris_state['velocity_vector_mps']
    expected_v_mag = vector_magnitude(expected_velocity_vector)

    velocity_vector = _get_body_orbital_velocity_vector(body_name, test_date, mu_sun)

    assert len(velocity_vector) == 3
    assert velocity_vector[0] == pytest.approx(expected_velocity_vector[0])
    assert velocity_vector[1] == pytest.approx(expected_velocity_vector[1])
    assert velocity_vector[2] == pytest.approx(expected_velocity_vector[2])
    assert vector_magnitude(velocity_vector) == pytest.approx(expected_v_mag)

def test_get_body_orbital_velocity_vector_ephemeris_not_found(mock_celestial_data, mock_ephemeris):
    """
    Tests _get_body_orbital_velocity_vector when no ephemeris data is available for the given body/date.
    Assumes _get_body_orbital_velocity_vector has been updated to expect structured data.
    """
    mock_ephemeris.side_effect = lambda body, date: None  # Simulate no ephemeris data

    body_name = 'Jupiter'
    test_date = datetime(2023, 1, 1, 0, 0, 0, tzinfo=None)
    mu_sun = mock_celestial_data.side_effect('sun')['gravitational_parameter_mu']

    with pytest.raises(ValueError, match=f"Ephemeris data not available for {body_name}"):
        _get_body_orbital_velocity_vector(body_name, test_date, mu_sun)

def test_get_body_orbital_velocity_vector_missing_velocity_data(mock_celestial_data, mock_ephemeris):
    """
    Tests _get_body_orbital_velocity_vector when ephemeris returns data but without 'velocity_vector_mps'.
    Assumes _get_body_orbital_velocity_vector has been updated to expect this key.
    """
    mock_ephemeris.side_effect = lambda name, date: {'position_vector_m': (1,2,3)} if name.lower() == 'badbody' else None
    body_name = 'BadBody'
    test_date = datetime(2023, 1, 1, 0, 0, 0, tzinfo=None)
    mu_sun = mock_celestial_data.side_effect('sun')['gravitational_parameter_mu']

    with pytest.raises(ValueError, match=f"Required 'velocity_vector_mps' missing for {body_name} in ephemeris data."):
        _get_body_orbital_velocity_vector(body_name, test_date, mu_sun)

def test_get_body_orbital_velocity_vector_invalid_velocity_data_format(mock_celestial_data, mock_ephemeris):
    """
    Tests _get_body_orbital_velocity_vector when ephemeris returns a malformed 'velocity_vector_mps'.
    Assumes _get_body_orbital_velocity_vector has been updated to validate the format.
    """
    mock_ephemeris.side_effect = lambda name, date: {'velocity_vector_mps': [1,2]} if name.lower() == 'badbody' else None
    body_name = 'BadBody'
    test_date = datetime(2023, 1, 1, 0, 0, 0, tzinfo=None)
    mu_sun = mock_celestial_data.side_effect('sun')['gravitational_parameter_mu']

    with pytest.raises(ValueError, match=f"Velocity vector must be a 3-element numeric tuple/list for {body_name}."):
        _get_body_orbital_velocity_vector(body_name, test_date, mu_sun)


def test_get_body_orbital_velocity_vector_invalid_mu(mock_celestial_data, mock_ephemeris):
    """
    Tests _get_body_orbital_velocity_vector with an invalid gravitational parameter.
    """
    body_name = 'Earth'
    test_date = datetime(2023, 1, 1, 0, 0, 0, tzinfo=None)
    invalid_mu = 0

    # Ensure ephemeris returns valid data, as the mu check happens before accessing it.
    mock_ephemeris.side_effect = lambda body, date: {
        'position_vector_m': (1.496e11, 0.0, 0.0),
        'velocity_vector_mps': (0.0, 29784.85, 0.0)
    } if body.lower() == 'earth' else None

    with pytest.raises(ValueError, match="Gravitational parameter \(mu\) must be positive"):
        _get_body_orbital_velocity_vector(body_name, test_date, invalid_mu)


# --- Simplified Lambert Transfer Test ---

def test_calculate_lambert_transfer_hohmann_like_collinear(mock_celestial_data):
    """
    Tests calculate_lambert_transfer for a simplified, collinear Hohmann-like transfer.
    This test adheres to the function's description of being a 'simplified analytical approximation',
    where delta_t is assumed to be half the transfer orbit period.
    We set up a scenario where r1 and r2 are collinear and the time of flight matches a half-ellipse.
    """
    mu_sun = mock_celestial_data.side_effect('sun')['gravitational_parameter_mu']
    r1_mag = 1.0 * 1.496e11  # Earth-like distance
    r2_mag = 1.52 * 1.496e11  # Mars-like distance

    # Position vectors for a collinear transfer (e.g., from +X to -X)
    r1_vec = (r1_mag, 0.0, 0.0)
    r2_vec = (-r2_mag, 0.0, 0.0)

    # Calculate the exact Hohmann transfer time for these radii, which is expected by the function
    a_transfer = (r1_mag + r2_mag) / 2
    delta_t_hohmann = math.pi * math.sqrt(a_transfer**3 / mu_sun)

    # The function calculates the magnitudes of the velocity vectors on the transfer ellipse
    # at the start and end points. These are *not* delta-V burns from circular orbits.
    expected_v_transfer_at_r1_mag = math.sqrt(mu_sun * ((2 / r1_mag) - (1 / a_transfer)))
    expected_v_transfer_at_r2_mag = math.sqrt(mu_sun * ((2 / r2_mag) - (1 / a_transfer)))
    expected_total_dv = expected_v_transfer_at_r1_mag + expected_v_transfer_at_r2_mag

    actual_v1_mag, actual_v2_mag, actual_total_dv, actual_tof = calculate_lambert_transfer(
        mu_sun, r1_vec, r2_vec, delta_t_hohmann
    )

    assert actual_v1_mag == pytest.approx(expected_v_transfer_at_r1_mag)
    assert actual_v2_mag == pytest.approx(expected_v_transfer_at_r2_mag)
    assert actual_total_dv == pytest.approx(expected_total_dv)
    assert actual_tof == pytest.approx(delta_t_hohmann)

def test_calculate_lambert_transfer_invalid_inputs():
    """Tests calculate_lambert_transfer with invalid input types or values."""
    mu_sun = 1.327e20
    r1_vec = (1.496e11, 0, 0)
    r2_vec = (2.279e11, 0, 0)
    delta_t = 1e7

    with pytest.raises(ValueError, match="Gravitational parameter \(mu\) must be positive"):
        calculate_lambert_transfer(0, r1_vec, r2_vec, delta_t)
    with pytest.raises(ValueError, match="Position vectors must be tuples or lists of 3 numeric elements"):
        calculate_lambert_transfer(mu_sun, (1, 2), r2_vec, delta_t)
    with pytest.raises(ValueError, match="Position vectors must be tuples or lists of 3 numeric elements"):
        calculate_lambert_transfer(mu_sun, r1_vec, (1, 2, 'a'), delta_t)
    with pytest.raises(ValueError, match="Time of flight \(delta_t\) must be positive"):
        calculate_lambert_transfer(mu_sun, r1_vec, r2_vec, 0)
    with pytest.raises(ValueError, match="Position vectors cannot be zero vectors"):
        calculate_lambert_transfer(mu_sun, (0, 0, 0), r2_vec, delta_t)
    with pytest.raises(ValueError, match="Position vectors cannot be zero vectors"):
        calculate_lambert_transfer(mu_sun, r1_vec, (0, 0, 0), delta_t)

def test_calculate_lambert_transfer_same_point_distinct_error():
    """
    Tests calculate_lambert_transfer when initial and final position vectors are the same.
    This should typically result in an error for a non-zero time of flight, as it's an ill-defined problem
    for a simple transfer.
    """
    mu_sun = 1.327e20
    r1_vec = (1.496e11, 0, 0)
    r2_vec = (1.496e11, 0, 0)  # Same as r1_vec
    delta_t = 1000  # Non-zero time

    with pytest.raises(ValueError, match="Position vectors must be distinct"):
        calculate_lambert_transfer(mu_sun, r1_vec, r2_vec, delta_t)

def test_calculate_lambert_transfer_impossible_simplified_transfer():
    """
    Tests calculate_lambert_transfer when the given time of flight is too short or too long
    for the simplified half-ellipse approximation between the two points.
    This should trigger the error condition mentioned in the function's description.
    """
    mu_sun = 1.327e20
    r1_vec = (1.496e11, 0, 0)
    r2_vec = (1.496e11 + 1e11, 0, 0)  # A bit further out
    delta_t = 100  # A very short time for this distance

    # The internal calculation for semi-major axis 'a' of the transfer ellipse
    # based on the simplified assumption (delta_t = pi * sqrt(a^3 / mu)) might lead
    # to an `a` value that is too small or other issues if the transfer is not a half-ellipse.
    # The error message should indicate this failure.
    with pytest.raises(ValueError, match="Transfer impossible or requires more advanced Lambert solver. Time of flight too short or semi-major axis calculation issue."):
        calculate_lambert_transfer(mu_sun, r1_vec, r2_vec, delta_t)