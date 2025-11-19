# tests/test_celestial_data.py
import pytest
import math
from celestial_data import (
    GRAVITATIONAL_CONSTANT,
    SPEED_OF_LIGHT,
    CELESTIAL_BODIES_DATA,
    get_celestial_body_data,
    get_all_solar_system_destinations
)
from constants import G_GRAVITATIONAL, C_LIGHT_MPS

# --- Test Constants ---

def test_gravitational_constant_value():
    """Verify GRAVITATIONAL_CONSTANT is correctly re-exported from constants."""
    assert GRAVITATIONAL_CONSTANT == pytest.approx(G_GRAVITATIONAL)
    assert isinstance(GRAVITATIONAL_CONSTANT, float)

def test_speed_of_light_value():
    """Verify SPEED_OF_LIGHT is correctly re-exported from constants."""
    assert SPEED_OF_LIGHT == pytest.approx(C_LIGHT_MPS)
    assert isinstance(SPEED_OF_LIGHT, float)

# --- Test get_celestial_body_data ---

def test_get_celestial_body_data_valid_body():
    """Test retrieval of data for a valid celestial body (Earth), including new orbital elements."""
    earth_data = get_celestial_body_data('Earth')
    assert earth_data is not None
    assert earth_data['name'] == 'Earth'
    assert earth_data['mass'] == pytest.approx(5.972e24)
    assert 'eccentricity' in earth_data
    assert 'inclination_deg' in earth_data
    assert 'longitude_of_ascending_node_deg' in earth_data
    assert 'argument_of_periapsis_deg' in earth_data
    assert 'mean_anomaly_at_j2000_deg' in earth_data
    assert isinstance(earth_data['eccentricity'], float)

def test_get_celestial_body_data_case_insensitivity():
    """Test retrieval with case-insensitive input."""
    mars_data_lower = get_celestial_body_data('mars')
    mars_data_upper = get_celestial_body_data('MARS')
    assert mars_data_lower is not None
    assert mars_data_upper is not None
    assert mars_data_lower == mars_data_upper
    assert mars_data_lower['name'] == 'Mars'

def test_get_celestial_body_data_non_existent_body():
    """Test retrieval for a celestial body that does not exist."""
    non_existent_data = get_celestial_body_data('Pluto')
    assert non_existent_data is None

def test_get_celestial_body_data_invalid_input_type():
    """Test retrieval with an invalid input type (e.g., non-string)."""
    invalid_data = get_celestial_body_data(123)
    assert invalid_data is None
    invalid_data_none = get_celestial_body_data(None)
    assert invalid_data_none is None

def test_get_celestial_body_data_sun_orbital_elements():
    """Test that the Sun's orbital elements are None as expected."""
    sun_data = get_celestial_body_data('Sun')
    assert sun_data is not None
    assert sun_data['name'] == 'Sun'
    assert sun_data['eccentricity'] is None
    assert sun_data['inclination_deg'] is None

def test_get_celestial_body_data_moon_geocentric_elements():
    """Test that the Moon's elements are present and noted as geocentric."""
    moon_data = get_celestial_body_data('Moon')
    assert moon_data is not None
    assert moon_data['name'] == 'Moon'
    assert 'average_distance_from_earth' in moon_data
    assert moon_data['eccentricity'] == pytest.approx(0.054900)
    assert moon_data['inclination_deg'] == pytest.approx(5.145)
    # The docstring specifies "These are GEOCENTRIC orbital elements for the Moon's orbit around Earth"

# --- Test get_all_solar_system_destinations ---

def test_get_all_solar_system_destinations_valid_output():
    """Test get_all_solar_system_destinations for correct format and sorting."""
    destinations = get_all_solar_system_destinations()
    assert isinstance(destinations, list)
    assert len(destinations) > 0

    # Earth should not be in the list
    assert not any(d['name'] == 'Earth' for d in destinations)

    # All items should be dictionaries with 'name' and 'distance'
    for dest in destinations:
        assert isinstance(dest, dict)
        assert 'name' in dest
        assert 'distance' in dest
        assert isinstance(dest['name'], str)
        assert isinstance(dest['distance'], float)
        assert dest['distance'] >= 0

    # Verify sorting by distance (ascending)
    for i in range(len(destinations) - 1):
        assert destinations[i]['distance'] <= destinations[i+1]['distance']

    # Check first and last destinations for plausibility (Moon closest, Neptune furthest)
    # The list is sorted, so 'Moon' should be first.
    # Note: If other bodies are added closer than Moon, this assert might fail.
    # But for solar system, Moon is the closest major body to Earth.
    assert destinations[0]['name'] == 'Moon'
    assert destinations[0]['distance'] == pytest.approx(CELESTIAL_BODIES_DATA['MOON']['average_distance_from_earth'])

    # Neptune should be the furthest among planets listed
    assert destinations[-1]['name'] == 'Neptune'
    expected_neptune_dist = math.fabs(CELESTIAL_BODIES_DATA['NEPTUNE']['semi_major_axis_from_sun'] - CELESTIAL_BODIES_DATA['EARTH']['semi_major_axis_from_sun'])
    assert destinations[-1]['distance'] == pytest.approx(expected_neptune_dist)

def test_get_all_solar_system_destinations_earth_data_missing(monkeypatch):
    """
    Test get_all_solar_system_destinations when Earth's data is missing.
    Should return None.
    """
    # Temporarily remove Earth from CELESTIAL_BODIES_DATA
    original_earth_data = CELESTIAL_BODIES_DATA.pop('EARTH', None)
    
    try:
        destinations = get_all_solar_system_destinations()
        assert destinations is None
    finally:
        # Restore Earth data
        if original_earth_data:
            CELESTIAL_BODIES_DATA['EARTH'] = original_earth_data

# tests/test_ephemeris.py
import pytest
import math
import datetime
from unittest.mock import patch, MagicMock
from ephemeris import (
    J2000_EPOCH,
    _get_mu_sun,
    _normalize_angle,
    _calculate_mean_anomaly,
    _solve_keplers_equation,
    _calculate_true_anomaly_and_radius,
    _elements_to_cartesian,
    get_heliocentric_state,
    _MU_SUN_CACHE # For clearing cache in tests
)

# Constants for test data (simplified, but consistent for mocks)
TEST_MU_SUN = 1.32712440018e20 # m^3/s^2
EARTH_SEMI_MAJOR_AXIS = 1.49598023e11 # m
EARTH_ORBITAL_PERIOD_S = 31557600.0 # s (~365.25 days)
EARTH_E = 0.01671022
EARTH_I_RAD = math.radians(0.00005)
EARTH_OMEGA_RAD = math.radians(-11.26064) # J2000 RAAN
EARTH_OMEGA_ARG_PERI_RAD = math.radians(102.93735) # J2000 Arg of Periapsis
EARTH_M0_RAD = math.radians(357.51716) # J2000 Mean Anomaly

MARS_SEMI_MAJOR_AXIS = 2.2793911e11 # m
MARS_ORBITAL_PERIOD_S = 59354000.0 # s (~687 days)
MARS_E = 0.09340056
MARS_I_RAD = math.radians(1.849726)
MARS_OMEGA_RAD = math.radians(49.5785367)
MARS_OMEGA_ARG_PERI_RAD = math.radians(286.4623033)
MARS_M0_RAD = math.radians(19.41846)

# --- Fixtures ---

@pytest.fixture
def mock_celestial_data_for_ephemeris():
    """
    Mocks `ephemeris.celestial_data.get_celestial_body_data` to provide
    extended orbital elements in radians and seconds.
    """
    with patch('ephemeris.celestial_data.get_celestial_body_data') as mock_get_data:
        mock_get_data.side_effect = lambda name: {
            'sun': {'name': 'Sun', 'gravitational_parameter_mu': TEST_MU_SUN},
            'earth': {
                'name': 'Earth',
                'semi_major_axis_from_sun': EARTH_SEMI_MAJOR_AXIS,
                'orbital_period': EARTH_ORBITAL_PERIOD_S,
                'eccentricity': EARTH_E,
                'inclination_rad': EARTH_I_RAD,
                'longitude_of_ascending_node_rad': EARTH_OMEGA_RAD,
                'argument_of_periapsis_rad': EARTH_OMEGA_ARG_PERI_RAD,
                'mean_anomaly_at_J2000_rad': EARTH_M0_RAD,
            },
            'mars': {
                'name': 'Mars',
                'semi_major_axis_from_sun': MARS_SEMI_MAJOR_AXIS,
                'orbital_period': MARS_ORBITAL_PERIOD_S,
                'eccentricity': MARS_E,
                'inclination_rad': MARS_I_RAD,
                'longitude_of_ascending_node_rad': MARS_OMEGA_RAD,
                'argument_of_periapsis_rad': MARS_OMEGA_ARG_PERI_RAD,
                'mean_anomaly_at_J2000_rad': MARS_M0_RAD,
            },
        }.get(name.lower())
        yield mock_get_data

@pytest.fixture
def mock_orbital_mechanics_for_ephemeris():
    """
    Mocks `ephemeris.orbital_mechanics.calculate_gravitational_parameter`.
    """
    with patch('ephemeris.orbital_mechanics.calculate_gravitational_parameter') as mock_calc_grav_param:
        mock_calc_grav_param.return_value = TEST_MU_SUN
        yield mock_calc_grav_param

@pytest.fixture(autouse=True)
def clear_mu_sun_cache():
    """Fixture to clear the _MU_SUN_CACHE before each test."""
    global _MU_SUN_CACHE
    _MU_SUN_CACHE = None
    yield
    _MU_SUN_CACHE = None # Clear after test too

# --- Test _get_mu_sun ---

def test_get_mu_sun_success(mock_orbital_mechanics_for_ephemeris: MagicMock):
    """Test successful retrieval and caching of Sun's gravitational parameter."""
    mu = _get_mu_sun()
    assert mu == pytest.approx(TEST_MU_SUN)
    mock_orbital_mechanics_for_ephemeris.assert_called_once_with('Sun')
    # Second call should use cache, not call mock again
    _get_mu_sun()
    mock_orbital_mechanics_for_ephemeris.assert_called_once() # Still only one call

def test_get_mu_sun_failure_no_data(mock_orbital_mechanics_for_ephemeris: MagicMock):
    """Test _get_mu_sun raises ValueError if Sun's data is missing."""
    mock_orbital_mechanics_for_ephemeris.side_effect = ValueError("Sun data not found in OM.")
    with patch('ephemeris.celestial_data.get_celestial_body_data', return_value=None):
        with pytest.raises(ValueError, match="Failed to retrieve Sun's gravitational parameter"):
            _get_mu_sun()

def test_get_mu_sun_failure_invalid_mu(mock_orbital_mechanics_for_ephemeris: MagicMock):
    """Test _get_mu_sun raises ValueError if Sun's mu is non-positive."""
    mock_orbital_mechanics_for_ephemeris.return_value = -1.0
    with pytest.raises(ValueError, match="Sun's gravitational parameter could not be determined or is invalid."):
        _get_mu_sun()

# --- Test Internal Helper Functions ---

def test_normalize_angle():
    """Test _normalize_angle for various inputs."""
    assert _normalize_angle(0) == pytest.approx(0)
    assert _normalize_angle(math.pi / 2) == pytest.approx(math.pi / 2)
    assert _normalize_angle(2 * math.pi) == pytest.approx(0)
    assert _normalize_angle(3 * math.pi) == pytest.approx(math.pi)
    assert _normalize_angle(-math.pi / 2) == pytest.approx(3 * math.pi / 2)
    assert _normalize_angle(-3 * math.pi) == pytest.approx(math.pi)

def test_calculate_mean_anomaly():
    """Test _calculate_mean_anomaly."""
    M0 = math.pi / 4
    mean_motion = (2 * math.pi) / 100 # 1 revolution in 100 seconds
    delta_t = 25 # Quarter of a revolution
    expected_M = M0 + mean_motion * delta_t # math.pi/4 + math.pi/2 = 3*math.pi/4
    assert _calculate_mean_anomaly(M0, mean_motion, delta_t) == pytest.approx(expected_M)

    # Test with normalization (more than one revolution)
    delta_t_long = 125 # 1.25 revolutions
    expected_M_normalized = M0 + mean_motion * delta_t_long # math.pi/4 + 1.25*2*math.pi = 2.5*math.pi + math.pi/4
    assert _calculate_mean_anomaly(M0, mean_motion, delta_t_long) == pytest.approx(math.pi / 4 + math.pi / 2)

def test_solve_keplers_equation_circular():
    """Test _solve_keplers_equation for a circular orbit (e=0)."""
    M = math.pi / 3
    E = _solve_keplers_equation(M, 0.0)
    assert E == pytest.approx(M)

def test_solve_keplers_equation_elliptical():
    """Test _solve_keplers_equation for an elliptical orbit."""
    M = math.radians(90) # 90 degrees mean anomaly
    e = 0.1
    # Expected E for M=90 deg, e=0.1 is approximately 95.77 degrees or 1.6715 rad
    # Using an online calculator (e.g., https://www.omnicalculator.com/physics/kepler-equation)
    expected_E = 1.6715396594 # radians for M=90deg, e=0.1
    E = _solve_keplers_equation(M, e)
    assert E == pytest.approx(expected_E, rel=1e-6)
    # Check Kepler's equation: E - e*sin(E) should be close to M
    assert (E - e * math.sin(E)) == pytest.approx(M, rel=1e-10)

def test_solve_keplers_equation_invalid_eccentricity():
    """Test _solve_keplers_equation with invalid eccentricity."""
    with pytest.raises(ValueError, match="Eccentricity .* not valid"):
        _solve_keplers_equation(math.pi/2, 1.0) # e >= 1
    with pytest.raises(ValueError, match="Eccentricity .* not valid"):
        _solve_keplers_equation(math.pi/2, -0.1) # e < 0

def test_calculate_true_anomaly_and_radius_circular():
    """Test _calculate_true_anomaly_and_radius for circular orbit."""
    a = 1e11
    e = 0.0
    E_rad = math.pi / 2 # E = M = nu for circular
    nu, r = _calculate_true_anomaly_and_radius(a, e, E_rad)
    assert nu == pytest.approx(math.pi / 2)
    assert r == pytest.approx(a) # Radius should be semi-major axis for circular

def test_calculate_true_anomaly_and_radius_elliptical():
    """Test _calculate_true_anomaly_and_radius for elliptical orbit."""
    a = 1.5e11 # semi-major axis
    e = 0.5
    E_rad = math.pi / 3 # Eccentric Anomaly
    # Expected radius: r = a * (1 - e * cos(E))
    expected_r = a * (1 - e * math.cos(E_rad))
    # Expected true anomaly: tan(nu/2) = sqrt((1+e)/(1-e)) * tan(E/2)
    expected_nu_half = math.sqrt((1 + e) / (1 - e)) * math.tan(E_rad / 2)
    expected_nu = 2 * math.atan(expected_nu_half)

    nu, r = _calculate_true_anomaly_and_radius(a, e, E_rad)
    assert r == pytest.approx(expected_r)
    assert nu == pytest.approx(_normalize_angle(expected_nu))

def test_calculate_true_anomaly_and_radius_invalid_semi_major_axis():
    """Test _calculate_true_anomaly_and_radius with invalid semi-major axis."""
    with pytest.raises(ValueError, match="Semi-major axis .* must be positive"):
        _calculate_true_anomaly_and_radius(0, 0.1, math.pi/2)
    with pytest.raises(ValueError, match="Semi-major axis .* must be positive"):
        _calculate_true_anomaly_and_radius(-1e11, 0.1, math.pi/2)

def test_elements_to_cartesian_circular_planar_orbit():
    """Test _elements_to_cartesian for a simple circular, planar orbit."""
    mu = TEST_MU_SUN
    a = EARTH_SEMI_MAJOR_AXIS
    e = 0.0
    i_rad = 0.0 # Planar on ecliptic
    Omega_rad = 0.0
    omega_rad = 0.0
    nu_rad = math.pi / 2 # At (0, a, 0)
    r_m = a # For circular orbit, r = a

    pos, vel = _elements_to_cartesian(mu, a, e, i_rad, Omega_rad, omega_rad, nu_rad, r_m)

    v_mag_circ = math.sqrt(mu / a)

    # Expected position at (0, a, 0)
    assert pos[0] == pytest.approx(0.0, abs=1e-6)
    assert pos[1] == pytest.approx(a, rel=1e-6)
    assert pos[2] == pytest.approx(0.0, abs=1e-6)

    # Expected velocity tangential (-v_mag, 0, 0)
    assert vel[0] == pytest.approx(-v_mag_circ, rel=1e-6)
    assert vel[1] == pytest.approx(0.0, abs=1e-6)
    assert vel[2] == pytest.approx(0.0, abs=1e-6)
    assert math.sqrt(sum(x**2 for x in vel)) == pytest.approx(v_mag_circ)


def test_elements_to_cartesian_inclined_orbit():
    """Test _elements_to_cartesian with inclination and RAAN."""
    mu = TEST_MU_SUN
    a = EARTH_SEMI_MAJOR_AXIS
    e = 0.0
    i_rad = math.radians(30) # 30 degrees inclination
    Omega_rad = math.radians(90) # Ascending node along Y-axis
    omega_rad = 0.0
    nu_rad = 0.0 # At periapsis, which is on ascending node
    r_m = a # For circular orbit

    pos, vel = _elements_to_cartesian(mu, a, e, i_rad, Omega_rad, omega_rad, nu_rad, r_m)

    # At nu=0 and omega=0, spacecraft is at the ascending node (r_I, r_J, r_K)
    # The periapsis is at RAAN, so it's a *rotation around Z by Omega_rad*, then X by i_rad
    # Since nu=0, pos in perifocal is (a, 0, 0)
    # R_Z(Omega) * R_X(i) * R_Z(omega=0) * (a, 0, 0)^T
    # This simplifies to R_Z(Omega) * R_X(i) * (a, 0, 0)^T
    # Which gives (a*cos(Omega), a*sin(Omega)*cos(i), a*sin(Omega)*sin(i)) -- wait, this is incorrect.
    # The coordinate transformation is standard, per Vallado:
    # x = r_PQW * (cos_omega * cos_Omega - sin_omega * sin_Omega * cos_i) ...
    # At nu = 0, x_perifocal = r_m, y_perifocal = 0
    # r_I = r_m * (cos_omega * cos_Omega - sin_omega * sin_Omega * cos_i)
    # r_J = r_m * (cos_omega * sin_Omega + sin_omega * cos_Omega * cos_i)
    # r_K = r_m * (sin_omega * sin_i)
    # With omega_rad = 0:
    # r_I = r_m * cos_Omega
    # r_J = r_m * sin_Omega
    # r_K = 0

    # With Omega_rad = 90 deg (pi/2)
    assert pos[0] == pytest.approx(0.0, abs=1e-6)
    assert pos[1] == pytest.approx(a, rel=1e-6)
    assert pos[2] == pytest.approx(0.0, abs=1e-6)

    # Velocity components (same rotation matrix applied)
    # At nu = 0, vx_perifocal = 0, vy_perifocal = sqrt(mu/p) (for e=0, p=a, so vy_perifocal = sqrt(mu/a))
    # v_I = - vy_perifocal * (sin_omega * cos_Omega + cos_omega * sin_Omega * cos_i)
    # v_J = + vy_perifocal * (cos_omega * cos_Omega * cos_i - sin_omega * sin_Omega)
    # v_K = + vy_perifocal * (cos_omega * sin_i)
    # With omega_rad = 0:
    # v_I = - vy_perifocal * sin_Omega * cos_i
    # v_J = + vy_perifocal * cos_Omega * cos_i
    # v_K = + vy_perifocal * sin_i

    # With Omega_rad = 90 deg (pi/2) and i_rad = 30 deg (pi/6)
    v_mag_circ = math.sqrt(mu / a)
    assert vel[0] == pytest.approx(-v_mag_circ * math.cos(i_rad), rel=1e-6) # -v_circ * cos(30)
    assert vel[1] == pytest.approx(0.0, abs=1e-6)
    assert vel[2] == pytest.approx(v_mag_circ * math.sin(i_rad), rel=1e-6) # v_circ * sin(30)
    assert math.sqrt(sum(x**2 for x in vel)) == pytest.approx(v_mag_circ)


def test_elements_to_cartesian_invalid_semi_latus_rectum():
    """Test _elements_to_cartesian when semi-latus rectum is non-positive."""
    mu = TEST_MU_SUN
    a = -1e11 # Invalid semi-major axis
    e = 0.5
    i_rad = 0.0
    Omega_rad = 0.0
    omega_rad = 0.0
    nu_rad = 0.0
    r_m = 1e11

    with pytest.raises(ValueError, match="Calculated semi-latus rectum is non-positive"):
        _elements_to_cartesian(mu, a, e, i_rad, Omega_rad, omega_rad, nu_rad, r_m)

# --- Test get_heliocentric_state ---

def test_get_heliocentric_state_earth_valid(
    mock_celestial_data_for_ephemeris: MagicMock,
    mock_orbital_mechanics_for_ephemeris: MagicMock
):
    """
    Test get_heliocentric_state for Earth at a specific date (J2000 epoch).
    This should put Earth close to its mean anomaly at J2000.
    """
    date = J2000_EPOCH # Exact J2000 epoch
    pos, vel = get_heliocentric_state('Earth', date)

    assert isinstance(pos, tuple) and len(pos) == 3
    assert isinstance(vel, tuple) and len(vel) == 3

    # Check magnitudes for plausibility (Earth's orbital radius and speed)
    pos_mag = math.sqrt(sum(x**2 for x in pos))
    vel_mag = math.sqrt(sum(x**2 for x in vel))

    # At J2000 epoch, mean anomaly is ~357.5 degrees (close to perihelion)
    # So distance should be slightly less than semi_major_axis, speed slightly more than circular velocity
    assert pos_mag == pytest.approx(EARTH_SEMI_MAJOR_AXIS * (1 - EARTH_E), rel=0.1) # Approx at perihelion
    assert vel_mag == pytest.approx(math.sqrt(TEST_MU_SUN * ((2 / pos_mag) - (1 / EARTH_SEMI_MAJOR_AXIS))), rel=0.1) # Vis-viva equation

    # Specifically check the first mock call for 'Earth' at 'J2000_EPOCH'
    mock_celestial_data_for_ephemeris.assert_called_with('Earth')
    mock_orbital_mechanics_for_ephemeris.assert_called_once() # For Sun's mu

def test_get_heliocentric_state_mars_valid_different_date(
    mock_celestial_data_for_ephemeris: MagicMock,
    mock_orbital_mechanics_for_ephemeris: MagicMock
):
    """Test get_heliocentric_state for Mars at a later date."""
    date = J2000_EPOCH + datetime.timedelta(days=100) # 100 days after J2000
    pos, vel = get_heliocentric_state('Mars', date)

    assert isinstance(pos, tuple) and len(pos) == 3
    assert isinstance(vel, tuple) and len(vel) == 3

    pos_mag = math.sqrt(sum(x**2 for x in pos))
    vel_mag = math.sqrt(sum(x**2 for x in vel))

    # Check magnitudes for plausibility (Mars' orbital radius and speed)
    # Semi-major axis is 2.279e11 m, eccentricity is 0.0934
    assert pos_mag == pytest.approx(MARS_SEMI_MAJOR_AXIS, rel=0.1) # Within 10% of SMA
    assert vel_mag == pytest.approx(math.sqrt(TEST_MU_SUN / MARS_SEMI_MAJOR_AXIS), rel=0.1) # Within 10% of circular V

    mock_celestial_data_for_ephemeris.assert_called_with('Mars') # Should be called for Mars

def test_get_heliocentric_state_invalid_body_name_type():
    """Test get_heliocentric_state with invalid body_name type."""
    with pytest.raises(ValueError, match="`body_name` must be a non-empty string."):
        get_heliocentric_state(123, datetime.datetime.now())

def test_get_heliocentric_state_empty_body_name():
    """Test get_heliocentric_state with an empty body_name string."""
    with pytest.raises(ValueError, match="`body_name` must be a non-empty string."):
        get_heliocentric_state("", datetime.datetime.now())

def test_get_heliocentric_state_non_existent_body(mock_celestial_data_for_ephemeris: MagicMock):
    """Test get_heliocentric_state for a non-existent celestial body."""
    mock_celestial_data_for_ephemeris.return_value = None # Ensure it returns None
    with pytest.raises(ValueError, match="Data for celestial body 'Krypton' not found."):
        get_heliocentric_state("Krypton", datetime.datetime.now())

def test_get_heliocentric_state_unsupported_sun():
    """Test get_heliocentric_state for the Sun (explicitly unsupported)."""
    with pytest.raises(ValueError, match="Cannot calculate heliocentric state for the Sun itself"):
        get_heliocentric_state("Sun", datetime.datetime.now())

def test_get_heliocentric_state_unsupported_moon():
    """Test get_heliocentric_state for the Moon (explicitly unsupported)."""
    with pytest.raises(ValueError, match="Heliocentric state for Earth's Moon is not directly supported"):
        get_heliocentric_state("Moon", datetime.datetime.now())
    with pytest.raises(ValueError, match="Heliocentric state for Earth's Moon is not directly supported"):
        get_heliocentric_state("Earth's Moon", datetime.datetime.now())


def test_get_heliocentric_state_invalid_date_type():
    """Test get_heliocentric_state with invalid date type."""
    with pytest.raises(ValueError, match="`date` must be a datetime.datetime object."):
        get_heliocentric_state("Earth", "not a date")

def test_get_heliocentric_state_missing_orbital_element(mock_celestial_data_for_ephemeris: MagicMock):
    """Test get_heliocentric_state when celestial data is missing an orbital element."""
    # Configure mock to return data without 'eccentricity'
    mock_celestial_data_for_ephemeris.side_effect = lambda name: \
        {'name': 'BadBody', 'semi_major_axis_from_sun': 1e11, 'orbital_period': 1e7,
         'eccentricity': None, # This makes it fail later in validation
         'inclination_rad': 0.0, 'longitude_of_ascending_node_rad': 0.0,
         'argument_of_periapsis_rad': 0.0, 'mean_anomaly_at_J2000_rad': 0.0} \
        if name.lower() == 'badbody' else None

    with pytest.raises(ValueError, match="Missing orbital element 'eccentricity' for 'BadBody'"):
        get_heliocentric_state("BadBody", datetime.datetime.now())

def test_get_heliocentric_state_invalid_orbital_element_value(mock_celestial_data_for_ephemeris: MagicMock):
    """Test get_heliocentric_state when celestial data has an invalid orbital element value."""
    # Configure mock to return data with invalid semi_major_axis_from_sun
    mock_celestial_data_for_ephemeris.side_effect = lambda name: \
        {'name': 'BadBody', 'semi_major_axis_from_sun': -1e11, 'orbital_period': 1e7,
         'eccentricity': 0.1, 'inclination_rad': 0.0, 'longitude_of_ascending_node_rad': 0.0,
         'argument_of_periapsis_rad': 0.0, 'mean_anomaly_at_J2000_rad': 0.0} \
        if name.lower() == 'badbody' else None

    with pytest.raises(ValueError, match="Invalid semi-major axis '-1.0e\\+11' for 'BadBody'"):
        get_heliocentric_state("BadBody", datetime.datetime.now())

def test_get_heliocentric_state_naive_datetime_conversion(
    mock_celestial_data_for_ephemeris: MagicMock,
    mock_orbital_mechanics_for_ephemeris: MagicMock
):
    """
    Test that a timezone-naive datetime is handled by converting it to UTC.
    The internal calculations should still yield correct results relative to J2000_EPOCH UTC.
    """
    naive_date = datetime.datetime(2023, 10, 26, 10, 30, 0) # No tzinfo
    aware_date = naive_date.astimezone(datetime.timezone.utc) # What it should convert to internally

    # We mock _get_mu_sun to avoid issues if other tests left it in a bad state
    mock_orbital_mechanics_for_ephemeris.return_value = TEST_MU_SUN

    # The calculation is complex, so we mostly ensure it runs without error and returns vectors
    pos, vel = get_heliocentric_state('Earth', naive_date)

    assert pos is not None
    assert vel is not None
    assert isinstance(pos, tuple) and len(pos) == 3
    assert isinstance(vel, tuple) and len(vel) == 3

    # Hard to assert exact values without a reference, but we ensure it processes
    # and the logic around J2000_EPOCH (which is UTC) would use this converted UTC time.
    # The crucial part is that `(date_utc - J2000_EPOCH).total_seconds()` gets called correctly.
    # We can't directly check the intermediate `date_utc` within `get_heliocentric_state` from here.
    # We primarily ensure no `ValueError` related to timezone handling is raised.
    # If the _elements_to_cartesian and prior functions work, this should too.

# tests/test_trajectory_planner.py
import pytest
from unittest.mock import MagicMock, patch, call
from datetime import datetime, timedelta
import math

# Import the module under test and its dependencies for patching
from trajectory_planner import TrajectoryPlanner
import celestial_data
import orbital_mechanics
import ephemeris
import constants # For G_GRAVITATIONAL if needed by orbital_mechanics mock

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
    Crucially, it now includes the new orbital elements in radians and orbital period in seconds,
    which `ephemeris` module expects.
    """
    # Using constants from ephemeris test for consistency
    sun_data = {'name': 'Sun', 'mass': 1.989e30, 'radius': 6.957e8, 'gravitational_parameter_mu': 1.32712440018e20, 'semi_major_axis_from_sun': 0, 'orbital_period': 0}
    earth_data = {
        'name': 'Earth', 'mass': 5.972e24, 'radius': 6.371e6, 'gravitational_parameter_mu': 3.986004418e14,
        'semi_major_axis_from_sun': 1.49598023e11, 'orbital_period': 31557600.0,
        'eccentricity': 0.01671022, 'inclination_rad': math.radians(0.00005),
        'longitude_of_ascending_node_rad': math.radians(-11.26064),
        'argument_of_periapsis_rad': math.radians(102.93735),
        'mean_anomaly_at_J2000_rad': math.radians(357.51716),
    }
    mars_data = {
        'name': 'Mars', 'mass': 6.417e23, 'radius': 3.389e6, 'gravitational_parameter_mu': 4.282837e13,
        'semi_major_axis_from_sun': 2.2793911e11, 'orbital_period': 59354000.0,
        'eccentricity': 0.09340056, 'inclination_rad': math.radians(1.849726),
        'longitude_of_ascending_node_rad': math.radians(49.5785367),
        'argument_of_periapsis_rad': math.radians(286.4623033),
        'mean_anomaly_at_J2000_rad': math.radians(19.41846),
    }
    moon_data = {
        'name': 'Moon', 'mass': 7.342e22, 'radius': 1.737e6, 'gravitational_parameter_mu': 4.905e12,
        'semi_major_axis_from_sun': 1.496e11, 'orbital_period': 2360584, # Geocentric period
        'eccentricity': 0.054900, 'inclination_rad': math.radians(5.145),
        'longitude_of_ascending_node_rad': math.radians(125.08),
        'argument_of_periapsis_rad': math.radians(318.15),
        'mean_anomaly_at_J2000_rad': math.radians(115.36),
    }

    with patch('trajectory_planner.celestial_data.get_celestial_body_data') as mock_get_data:
        mock_get_data.side_effect = lambda name: {
            'sun': sun_data,
            'earth': earth_data,
            'mars': mars_data,
            'moon': moon_data,
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
         patch('trajectory_planner.orbital_mechanics.vector_magnitude') as mock_vector_magnitude:

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

            # Check for identical position vectors (e.g., Earth to Earth)
            if all(math.isclose(a, b, abs_tol=1e-9) for a, b in zip(r1_vec, r2_vec)):
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
        mock_vector_magnitude.side_effect = lambda vec: math.sqrt(sum(x**2 for x in vec)) # Real implementation

        yield {
            'hohmann_dv': mock_hohmann_dv,
            'hohmann_tof': mock_hohmann_tof,
            'lambert_transfer': mock_lambert_transfer,
            'vector_magnitude': mock_vector_magnitude,
        }

@pytest.fixture
def mock_ephemeris_get_heliocentric_state():
    """
    Fixture to mock ephemeris.get_heliocentric_state.
    Provides illustrative heliocentric states (position and velocity vectors) for
    Earth and Mars on specific dates, representing a simplified 3D elliptical orbit model.
    """
    with patch('trajectory_planner.ephemeris.get_heliocentric_state') as mock_get_state:
        # Define some illustrative states (position_vector_m, velocity_vector_mps)
        # Assuming simplified 2D planar orbits for mock, z components are 0.
        # Earth at ~1 AU, velocity ~29.8 km/s
        earth_r = 1.49598023e11 # m (Earth's semi-major axis)
        earth_v = 29800.0 # m/s (approx circular orbital velocity)

        # Mars at ~1.52 AU, velocity ~24.1 km/s
        mars_r = 2.2793911e11 # m (Mars' semi-major axis)
        mars_v = 24100.0 # m/s (approx circular orbital velocity)

        def get_mock_state_effect(body_name, date):
            # Define specific positions/velocities for certain dates
            # These are simplified for testing; actual ephemeris would calculate them dynamically.
            # J2000_EPOCH is datetime.datetime(2000, 1, 1, 12, 0, 0, tzinfo=datetime.timezone.utc)
            # We'll use a fixed reference date for mocks for simplicity.
            ref_date = _parse_date_string('2024-01-01') # Reference point for our mock state

            # Calculate time difference from ref_date
            delta_days = (date - ref_date).total_seconds() / (24 * 3600)

            # Earth: ~365.25 days period
            earth_period_days = 365.25
            earth_angle_offset = (delta_days % earth_period_days / earth_period_days) * 2 * math.pi
            # Start Earth at (r, 0, 0) on ref_date
            earth_pos_x = earth_r * math.cos(earth_angle_offset)
            earth_pos_y = earth_r * math.sin(earth_angle_offset)
            earth_vel_x = -earth_v * math.sin(earth_angle_offset)
            earth_vel_y = earth_v * math.cos(earth_angle_offset)

            # Mars: ~687 days period
            mars_period_days = 687.0
            mars_angle_offset = (delta_days % mars_period_days / mars_period_days) * 2 * math.pi
            # Start Mars at (-r, 0, 0) (opposition) on ref_date
            mars_pos_x = mars_r * math.cos(math.pi + mars_angle_offset)
            mars_pos_y = mars_r * math.sin(math.pi + mars_angle_offset)
            mars_vel_x = -mars_v * math.sin(math.pi + mars_angle_offset)
            mars_vel_y = mars_v * math.cos(math.pi + mars_angle_offset)

            if body_name.lower() == 'earth':
                return {
                    'position_vector_m': (earth_pos_x, earth_pos_y, 0.0),
                    'velocity_vector_mps': (earth_vel_x, earth_vel_y, 0.0)
                }
            elif body_name.lower() == 'mars':
                return {
                    'position_vector_m': (mars_pos_x, mars_pos_y, 0.0),
                    'velocity_vector_mps': (mars_vel_x, mars_vel_y, 0.0)
                }
            elif body_name.lower() == 'venus':
                # Just some mock data for Venus to pass direct transfer test
                venus_r = 1.0820893e11 # m
                venus_v = 35000.0 # m/s
                venus_angle_offset = (delta_days % 224.7 / 224.7) * 2 * math.pi
                venus_pos_x = venus_r * math.cos(venus_angle_offset + math.pi/2) # Offset from Earth
                venus_pos_y = venus_r * math.sin(venus_angle_offset + math.pi/2)
                venus_vel_x = -venus_v * math.sin(venus_angle_offset + math.pi/2)
                venus_vel_y = venus_v * math.cos(venus_angle_offset + math.pi/2)
                return {
                    'position_vector_m': (venus_pos_x, venus_pos_y, 0.0),
                    'velocity_vector_mps': (venus_vel_x, venus_vel_y, 0.0)
                }
            return None # For unsupported bodies or dates out of scope if needed

        mock_get_state.side_effect = get_mock_state_effect
        yield mock_get_state

# --- Test TrajectoryPlanner.__init__ ---

def test_trajectory_planner_init_success(mock_celestial_data_for_planner: MagicMock):
    """
    Tests successful initialization of TrajectoryPlanner, ensuring Sun's gravitational
    parameter is correctly loaded.
    """
    planner = TrajectoryPlanner()
    # Using the exact value from mock_celestial_data_for_planner for 'Sun'
    assert planner.mu_sun == pytest.approx(1.32712440018e20)
    mock_celestial_data_for_planner.assert_called_once_with('Sun')

def test_trajectory_planner_init_missing_sun_data(mock_celestial_data_for_planner: MagicMock):
    """
    Tests initialization failure when Sun data cannot be retrieved from `celestial_data`.
    """
    mock_celestial_data_for_planner.return_value = None
    with pytest.raises(ValueError, match="Data for 'Sun' not found"):
        TrajectoryPlanner()

def test_trajectory_planner_init_invalid_sun_mu(mock_celestial_data_for_planner: MagicMock):
    """
    Tests initialization failure when Sun's gravitational parameter is invalid or missing.
    """
    # Simulate valid sun data but with a problematic mu value
    mock_celestial_data_for_planner.side_effect = lambda name: {'name': 'Sun', 'gravitational_parameter_mu': -1.0} if name.lower() == 'sun' else None
    with pytest.raises(ValueError, match="Invalid or non-positive 'gravitational_parameter_mu'"):
        TrajectoryPlanner()

# --- Test TrajectoryPlanner.plan_trajectory (Hohmann Type) ---

def test_plan_hohmann_trajectory_success_ephemeris(
    mock_celestial_data_for_planner: MagicMock, # Not directly called for radii, but for mu_sun in init
    mock_orbital_mechanics: dict,
    mock_ephemeris_get_heliocentric_state: MagicMock
):
    """
    Tests successful Hohmann transfer planning, verifying that ephemeris data (position vector magnitude) is
    prioritized for calculating radii and that orbital mechanics functions are called.
    """
    planner = TrajectoryPlanner()
    departure_date = _parse_date_string('2024-01-01')

    # Get specific mock states for assertion. The side_effect function is defined
    # so we can call it here to get the expected outputs.
    earth_state = mock_ephemeris_get_heliocentric_state.side_effect('Earth', departure_date)
    mars_state = mock_ephemeris_get_heliocentric_state.side_effect('Mars', departure_date)
    
    earth_r_mag = orbital_mechanics.vector_magnitude(earth_state['position_vector_m'])
    mars_r_mag = orbital_mechanics.vector_magnitude(mars_state['position_vector_m'])

    result = planner.plan_trajectory(
        departure_body_name='Earth',
        arrival_body_name='Mars',
        trajectory_type='Hohmann',
        departure_date=departure_date.strftime('%Y-%m-%d') # Pass as string as per plan_trajectory
    )

    assert result['success'] is True
    assert result['total_delta_v_mps'] == 1500.0  # From mock_hohmann_dv
    assert result['travel_time_days'] == pytest.approx(100.0) # From mock_hohmann_tof (in days)
    assert result['transfer_type'] == 'Hohmann Transfer'
    assert result['departure_distance_from_sun_m'] == pytest.approx(earth_r_mag)
    assert result['arrival_distance_from_sun_m'] == pytest.approx(mars_r_mag)

    mock_ephemeris_get_heliocentric_state.assert_has_calls([
        call('earth', departure_date),
        call('mars', departure_date)
    ], any_order=True)

    mock_orbital_mechanics['hohmann_dv'].assert_called_once_with(
        planner.mu_sun, earth_r_mag, mars_r_mag
    )
    mock_orbital_mechanics['hohmann_tof'].assert_called_once_with(
        planner.mu_sun, earth_r_mag, mars_r_mag
    )

def test_plan_hohmann_trajectory_ephemeris_error(
    mock_celestial_data_for_planner: MagicMock,
    mock_ephemeris_get_heliocentric_state: MagicMock
):
    """
    Tests that Hohmann planning returns an error if ephemeris data retrieval fails
    for either body (no fallback to static data, as per current implementation).
    """
    mock_ephemeris_get_heliocentric_state.side_effect = ValueError("Ephemeris data unavailable for Mars.")

    planner = TrajectoryPlanner()
    departure_date = _parse_date_string('2024-01-01')
    result = planner.plan_trajectory(
        departure_body_name='Earth',
        arrival_body_name='Mars',
        trajectory_type='Hohmann',
        departure_date=departure_date.strftime('%Y-%m-%d')
    )

    assert result['success'] is False
    assert "Error retrieving orbital states via ephemeris for date" in result['error']
    assert "Ephemeris data unavailable for Mars." in result['error']
    mock_ephemeris_get_heliocentric_state.assert_has_calls([
        call('earth', departure_date),
        call('mars', departure_date) # It tries for both before the error is raised by the side_effect
    ], any_order=True)

def test_plan_hohmann_trajectory_orbital_mechanics_error(
    mock_celestial_data_for_planner: MagicMock,
    mock_orbital_mechanics: dict,
    mock_ephemeris_get_heliocentric_state: MagicMock
):
    """
    Tests error handling when underlying orbital mechanics functions for Hohmann transfer fail.
    """
    mock_orbital_mechanics['hohmann_dv'].side_effect = ValueError("Hohmann DV calculation logic error.")
    planner = TrajectoryPlanner()
    departure_date = _parse_date_string('2024-01-01')
    result = planner.plan_trajectory(
        departure_body_name='Earth',
        arrival_body_name='Mars',
        trajectory_type='Hohmann',
        departure_date=departure_date.strftime('%Y-%m-%d')
    )

    assert result['success'] is False
    assert "Orbital mechanics calculation error for Hohmann transfer: Hohmann DV calculation logic error." in result['error']

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

    # Expected optimal values from the mock_lambert_solver_effect (around 200 days / 17.28M seconds)
    expected_travel_time_days = 200.0 # From mock (17_280_000 seconds)
    expected_total_dv = 15_000.0 # Minimum total_dv from the mock

    result = planner.plan_trajectory(
        departure_body_name='Earth',
        arrival_body_name='Mars',
        trajectory_type='Direct',
        departure_date=departure_date.strftime('%Y-%m-%d')
    )

    assert result['success'] is True
    assert math.isclose(result['total_delta_v_mps'], expected_total_dv, rel_tol=1e-3)
    assert math.isclose(result['travel_time_days'], expected_travel_time_days, rel_tol=1e-3)
    assert result['transfer_type'] == 'Direct Transfer'
    assert result['departure_distance_from_sun_m'] > 0
    assert result['arrival_distance_from_sun_m'] > 0

    # Verify that calculate_lambert_transfer was called multiple times, indicating iteration
    assert mock_orbital_mechanics['lambert_transfer'].call_count > 1
    
    # Verify that ephemeris was called for both bodies, including for the iterated target date
    mock_ephemeris_get_heliocentric_state.assert_any_call('earth', departure_date)
    # It should call Mars multiple times for different dates (departure_date + delta_t)
    # Let's check a call for Mars at the optimal TOF from the mock
    optimal_arrival_date = departure_date + timedelta(seconds=expected_travel_time_days * 24 * 3600)
    mock_ephemeris_get_heliocentric_state.assert_any_call('mars', optimal_arrival_date)


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

    # Make mock Lambert solver always return a very high delta_v for any valid input,
    # to simulate a situation where no "optimal" solution is found by the planner.
    def consistently_high_dv_solver(mu, r1_vec, r2_vec, delta_t):
        if delta_t <= 0:
            raise ValueError("Time of flight must be positive for Lambert transfer.")
        if all(math.isclose(a, b, abs_tol=1e-9) for a, b in zip(r1_vec, r2_vec)):
            return 0.0, 0.0, 0.0, delta_t
        # Returns a dv that's always considered 'suboptimal' by the planner's thresholds
        return 100_000.0, 100_000.0, 200_000.0, delta_t # Force high DV
    mock_orbital_mechanics['lambert_transfer'].side_effect = consistently_high_dv_solver

    result = planner.plan_trajectory(
        departure_body_name='Earth',
        arrival_body_name='Mars',
        trajectory_type='Direct',
        departure_date=departure_date.strftime('%Y-%m-%d')
    )

    assert result['success'] is False
    assert "No suitable Direct Transfer trajectory found" in result['error']
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

    # Make mock Lambert solver always raise an error
    mock_orbital_mechanics['lambert_transfer'].side_effect = ValueError("Lambert calculation failed catastrophically due to bad geometry.")

    result = planner.plan_trajectory(
        departure_body_name='Earth',
        arrival_body_name='Mars',
        trajectory_type='Direct',
        departure_date=departure_date.strftime('%Y-%m-%d')
    )

    assert result['success'] is False
    assert "Lambert calculation failed catastrophically due to bad geometry." in result['error']
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

    result = planner.plan_trajectory(
        departure_body_name='Earth',
        arrival_body_name='Earth', # Plan from Earth to Earth
        trajectory_type='Direct',
        departure_date=departure_date.strftime('%Y-%m-%d')
    )

    assert result['success'] is False # The new plan_trajectory will return false as "Departure and arrival bodies cannot be the same."
    assert "Departure and arrival bodies cannot be the same." in result['error']

    # The planner *should* recognize identical bodies BEFORE calling ephemeris or orbital mechanics
    # for positions or Lambert calculations.
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
        departure_date=_parse_date_string('2024-01-01').strftime('%Y-%m-%d')
    )
    assert result['success'] is False
    assert "Unsupported trajectory type: 'UnsupportedType'" in result['error']

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
    assert "Invalid 'departure_date' format: '2024/01/01'. Expected 'YYYY-MM-DD'." in result['error']

def test_plan_trajectory_missing_departure_date_for_dynamic_plan(mock_celestial_data_for_planner: MagicMock):
    """
    Tests `plan_trajectory` when `departure_date` is omitted for a trajectory type
    that dynamically requires it (e.g., 'Direct' or 'Hohmann' using ephemeris).
    """
    planner = TrajectoryPlanner()
    result = planner.plan_trajectory(
        departure_body_name='Earth',
        arrival_body_name='Mars',
        trajectory_type='Direct' # Requires departure_date, but not provided in kwargs
    )
    assert result['success'] is False
    assert "'Direct' trajectory planning requires a 'departure_date' in kwargs." in result['error']

def test_plan_trajectory_invalid_departure_body_name(mock_celestial_data_for_planner: MagicMock):
    """Tests `plan_trajectory` with an invalid (non-string) departure body name."""
    planner = TrajectoryPlanner()
    result = planner.plan_trajectory(
        departure_body_name=123,
        arrival_body_name='Mars',
        trajectory_type='Hohmann',
        departure_date='2024-01-01'
    )
    assert result['success'] is False
    assert "Departure body name must be a non-empty string." in result['error']

def test_plan_trajectory_empty_arrival_body_name(mock_celestial_data_for_planner: MagicMock):
    """Tests `plan_trajectory` with an empty arrival body name."""
    planner = TrajectoryPlanner()
    result = planner.plan_trajectory(
        departure_body_name='Earth',
        arrival_body_name='  ',
        trajectory_type='Hohmann',
        departure_date='2024-01-01'
    )
    assert result['success'] is False
    assert "Arrival body name must be a non-empty string." in result['error']

# tests/test_main_helpers.py
import pytest
import math
from datetime import datetime, timedelta
from main import calc_time_earth, calc_age, calc_arrival, convert_date, YEAR_TO_SECONDS
from constants import C_LIGHT_MPS

# --- Test calc_time_earth ---

def test_calc_time_earth_relativistic_effect():
    """Test calc_time_earth with a speed less than c, expecting time dilation."""
    traveler_time = 1.0  # year
    average_speed = 0.5 * C_LIGHT_MPS # Half the speed of light

    expected_lorentz_factor = 1 / math.sqrt(1 - (0.5)**2)
    expected_earth_time = traveler_time * expected_lorentz_factor

    earth_time = calc_time_earth(traveler_time, average_speed)
    assert earth_time == pytest.approx(expected_earth_time)
    assert earth_time > traveler_time # Earth time should be greater

def test_calc_time_earth_zero_speed():
    """Test calc_time_earth with zero speed, expecting infinite time."""
    traveler_time = 1.0
    average_speed = 0.0
    earth_time = calc_time_earth(traveler_time, average_speed)
    assert earth_time == float('inf')

def test_calc_time_earth_speed_of_light():
    """Test calc_time_earth with speed equal to C_LIGHT_MPS."""
    traveler_time = 1.0
    average_speed = C_LIGHT_MPS
    earth_time = calc_time_earth(traveler_time, average_speed)
    assert earth_time == pytest.approx(traveler_time) # Returns traveler's time as fallback/warning

def test_calc_time_earth_greater_than_speed_of_light():
    """Test calc_time_earth with speed greater than C_LIGHT_MPS."""
    traveler_time = 1.0
    average_speed = 1.1 * C_LIGHT_MPS
    earth_time = calc_time_earth(traveler_time, average_speed)
    assert earth_time == pytest.approx(traveler_time) # Returns traveler's time as fallback/warning

def test_calc_time_earth_negative_speed():
    """Test calc_time_earth with a negative speed, expecting traveler's time."""
    traveler_time = 1.0
    average_speed = -1000.0 # m/s
    earth_time = calc_time_earth(traveler_time, average_speed)
    assert earth_time == pytest.approx(traveler_time) # Returns traveler's time as fallback/warning

# --- Test calc_age ---

def test_calc_age_normal():
    """Test calc_age with positive inputs."""
    travel_time = 5.0 # years
    starting_age = 30 # years
    expected_age = 35.0
    assert calc_age(travel_time, starting_age) == pytest.approx(expected_age)

def test_calc_age_zero_travel_time():
    """Test calc_age with zero travel time."""
    travel_time = 0.0
    starting_age = 25
    expected_age = 25.0
    assert calc_age(travel_time, starting_age) == pytest.approx(expected_age)

def test_calc_age_invalid_travel_time():
    """Test calc_age with negative travel time, expecting ValueError."""
    with pytest.raises(ValueError, match="Travel time must be a non-negative number."):
        calc_age(-1.0, 30)

def test_calc_age_invalid_starting_age():
    """Test calc_age with negative starting age, expecting ValueError."""
    with pytest.raises(ValueError, match="Starting age must be a non-negative number."):
        calc_age(5.0, -10)
    with pytest.raises(ValueError, match="Starting age must be a non-negative number."):
        calc_age(5.0, "invalid") # Non-numeric

# --- Test calc_arrival ---

def test_calc_arrival_normal():
    """Test calc_arrival with positive travel time."""
    departure = datetime(2024, 1, 1)
    travel_years = 1.0 # 1 Earth year
    expected_arrival = datetime(2025, 1, 1, 6, 0, 0) # 365.25 days = 1 year and 6 hours
    
    arrival = calc_arrival(departure, travel_years)
    assert arrival is not None
    # Compare year, month, day, hour, minute for approximate match
    # Since 365.25 days, expecting to be around 6 AM on Jan 1, 2025
    assert arrival.year == 2025
    assert arrival.month == 1
    assert arrival.day == 1
    assert arrival.hour == 6 # Due to 0.25 day = 6 hours
    
    # Check total seconds difference
    total_seconds_expected = travel_years * YEAR_TO_SECONDS
    total_seconds_actual = (arrival - departure).total_seconds()
    assert total_seconds_actual == pytest.approx(total_seconds_expected, rel=1e-3) # A bit of tolerance

def test_calc_arrival_zero_travel_time():
    """Test calc_arrival with zero travel time."""
    departure = datetime(2024, 1, 1)
    travel_years = 0.0
    expected_arrival = datetime(2024, 1, 1)
    assert calc_arrival(departure, travel_years) == expected_arrival

def test_calc_arrival_invalid_departure_date():
    """Test calc_arrival with an invalid departure date type, expecting None."""
    travel_years = 1.0
    arrival = calc_arrival("not a date", travel_years)
    assert arrival is None

def test_calc_arrival_extremely_long_travel_time():
    """Test calc_arrival with extremely long travel time leading to OverflowError."""
    departure = datetime(2024, 1, 1)
    # A travel time that would exceed the maximum representable datetime
    very_long_years = 1_000_000_000_000_000_000.0 # 1 quintillion years
    arrival = calc_arrival(departure, very_long_years)
    assert arrival is None # Should return None due to OverflowError handling

# --- Test convert_date ---

def test_convert_date_valid_format():
    """Test convert_date with a valid YYYY-MM-DD string."""
    date_str = "2023-10-26"
    expected_date = datetime(2023, 10, 26)
    assert convert_date(date_str) == expected_date

def test_convert_date_invalid_format():
    """Test convert_date with an invalid date string format, expecting ValueError."""
    with pytest.raises(ValueError):
        convert_date("26/10/2023")
    with pytest.raises(ValueError):
        convert_date("20231026")
    with pytest.raises(ValueError):
        convert_date("not-a-date")