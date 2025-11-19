import pytest
import math
from unittest.mock import patch, MagicMock

# Assuming celestial_data.py and constants.py are in the same parent directory
# or that celestial_data correctly imports from constants.
from celestial_data import (
    GRAVITATIONAL_CONSTANT, SPEED_OF_LIGHT, CELESTIAL_BODIES_DATA,
    get_celestial_body_data, get_all_solar_system_destinations
)

# Mock constants.py to ensure consistent testing environment
class MockConstants:
    G_GRAVITATIONAL = 6.67430e-11
    C_LIGHT_MPS = 299792458.0

@pytest.fixture(autouse=True)
def mock_constants_in_celestial_data(mocker):
    """Automatically mock constants in celestial_data for consistent testing."""
    mocker.patch('celestial_data.G_GRAVITATIONAL', MockConstants.G_GRAVITATIONAL)
    mocker.patch('celestial_data.C_LIGHT_MPS', MockConstants.C_LIGHT_MPS)
    mocker.patch('celestial_data.constants.G_GRAVITATIONAL', MockConstants.G_GRAVITATIONAL)
    mocker.patch('celestial_data.constants.C_LIGHT_MPS', MockConstants.C_LIGHT_MPS)


def test_gravitational_constant_re_export():
    """Test that GRAVITATIONAL_CONSTANT is correctly re-exported."""
    assert GRAVITATIONAL_CONSTANT == MockConstants.G_GRAVITATIONAL


def test_speed_of_light_re_export():
    """Test that SPEED_OF_LIGHT is correctly re-exported."""
    assert SPEED_OF_LIGHT == MockConstants.C_LIGHT_MPS


@pytest.mark.parametrize("body_name, expected_name_in_data", [
    ("Earth", "EARTH"),
    ("earth", "EARTH"),
    ("EARTH", "EARTH"),
    ("Moon", "MOON"),
    ("mArS", "MARS"),
])
def test_get_celestial_body_data_valid_names(body_name, expected_name_in_data):
    """Test get_celestial_body_data with valid, case-insensitive body names."""
    data = get_celestial_body_data(body_name)
    assert data is not None
    assert isinstance(data, dict)
    assert data == CELESTIAL_BODIES_DATA[expected_name_in_data]


def test_get_celestial_body_data_non_existent():
    """Test get_celestial_body_data with a non-existent body name."""
    data = get_celestial_body_data("Krypton")
    assert data is None


@pytest.mark.parametrize("invalid_input", [
    None,
    123,
    [],
    {}
])
def test_get_celestial_body_data_invalid_input_type(invalid_input):
    """Test get_celestial_body_data with invalid input types."""
    data = get_celestial_body_data(invalid_input)
    assert data is None


def test_get_celestial_body_data_data_integrity():
    """Test that essential keys are present for a known celestial body."""
    earth_data = get_celestial_body_data("Earth")
    assert "mass" in earth_data
    assert "radius" in earth_data
    assert "gravitational_parameter_mu" in earth_data
    assert "semi_major_axis_from_sun" in earth_data


def test_get_all_solar_system_destinations_basic():
    """Test basic functionality of get_all_solar_system_destinations.
    Ensures Earth is excluded, list is sorted, and Moon's explicit distance is used.
    """
    destinations = get_all_solar_system_destinations()
    assert isinstance(destinations, list)
    assert len(destinations) > 0

    # Check Earth is not in destinations
    assert not any(d['name'] == 'Earth' for d in destinations)

    # Check sorting by distance
    for i in range(len(destinations) - 1):
        assert destinations[i]['distance'] <= destinations[i+1]['distance']

    # Check Moon's presence and distance (using actual data)
    moon_entry = next((d for d in destinations if d['name'] == 'Moon'), None)
    assert moon_entry is not None
    assert math.isclose(moon_entry['distance'], CELESTIAL_BODIES_DATA["MOON"]["average_distance_from_earth"])

    # Check a few other known planets
    mars_entry = next((d for d in destinations if d['name'] == 'Mars'), None)
    assert mars_entry is not None
    # The calculated distance for Mars should be |Mars_SMA - Earth_SMA|
    earth_sma = CELESTIAL_BODIES_DATA["EARTH"]["semi_major_axis_from_sun"]
    mars_sma = CELESTIAL_BODIES_DATA["MARS"]["semi_major_axis_from_sun"]
    assert math.isclose(mars_entry['distance'], abs(mars_sma - earth_sma))


def test_get_all_solar_system_destinations_earth_data_missing(mocker, capsys):
    """Test get_all_solar_system_destinations when Earth's data is critically missing."""
    # Mock get_celestial_body_data to return None for 'Earth'
    original_get_data = get_celestial_body_data
    def mock_get_celestial_body_data_for_earth_missing(body_name):
        if body_name.upper() == "EARTH":
            return None
        return original_get_data(body_name) # Call original for other bodies

    mocker.patch('celestial_data.get_celestial_body_data', side_effect=mock_get_celestial_body_data_for_earth_missing)

    destinations = get_all_solar_system_destinations()
    assert destinations is None
    captured = capsys.readouterr()
    assert "Error: Earth's data not found in CELESTIAL_BODIES_DATA." in captured.out


def test_get_all_solar_system_destinations_earth_semi_major_axis_missing(mocker, capsys):
    """Test get_all_solar_system_destinations when Earth's semi_major_axis_from_sun is missing."""
    # Mock get_celestial_body_data to return Earth data without 'semi_major_axis_from_sun'
    mock_earth_data_incomplete = CELESTIAL_BODIES_DATA["EARTH"].copy()
    del mock_earth_data_incomplete["semi_major_axis_from_sun"]

    original_get_data = get_celestial_body_data
    def mock_get_celestial_body_data_for_earth_incomplete(body_name):
        if body_name.upper() == "EARTH":
            return mock_earth_data_incomplete
        return original_get_data(body_name) # Call original for other bodies

    mocker.patch('celestial_data.get_celestial_body_data', side_effect=mock_get_celestial_body_data_for_earth_incomplete)

    destinations = get_all_solar_system_destinations()
    assert destinations is None
    captured = capsys.readouterr()
    assert "Error: Earth's 'semi_major_axis_from_sun' data not found." in captured.out


def test_get_all_solar_system_destinations_moon_distance_missing(mocker, capsys):
    """Test get_all_solar_system_destinations when Moon's average_distance_from_earth is missing."""
    mock_moon_data_incomplete = CELESTIAL_BODIES_DATA["MOON"].copy()
    del mock_moon_data_incomplete["average_distance_from_earth"]

    original_get_data = get_celestial_body_data
    def mock_get_celestial_body_data_for_moon_incomplete(body_name):
        if body_name.upper() == "MOON":
            return mock_moon_data_incomplete
        return original_get_data(body_name) # Call original for other bodies

    mocker.patch('celestial_data.get_celestial_body_data', side_effect=mock_get_celestial_body_data_for_moon_incomplete)

    destinations = get_all_solar_system_destinations()
    assert destinations is not None
    assert not any(d['name'] == 'Moon' for d in destinations)
    captured = capsys.readouterr()
    assert "Warning: Moon's 'average_distance_from_earth' data not found. Skipping Moon." in captured.err


def test_get_all_solar_system_destinations_planet_semi_major_axis_missing(mocker, capsys):
    """Test get_all_solar_system_destinations when a planet's semi_major_axis_from_sun is missing."""
    # Mock Mars data to be incomplete
    mock_mars_data_incomplete = CELESTIAL_BODIES_DATA["MARS"].copy()
    del mock_mars_data_incomplete["semi_major_axis_from_sun"]

    original_get_data = get_celestial_body_data
    def mock_get_celestial_body_data_for_mars_incomplete(body_name):
        if body_name.upper() == "MARS":
            return mock_mars_data_incomplete
        return original_get_data(body_name) # Call original for other bodies

    mocker.patch('celestial_data.get_celestial_body_data', side_effect=mock_get_celestial_body_data_for_mars_incomplete)

    destinations = get_all_solar_system_destinations()
    assert destinations is not None
    assert not any(d['name'] == 'Mars' for d in destinations)
    captured = capsys.readouterr()
    assert "Warning: MARS's 'semi_major_axis_from_sun' data not found. Skipping MARS." in captured.err
