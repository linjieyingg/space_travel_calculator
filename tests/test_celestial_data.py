# tests/test_celestial_data.py
import pytest
import math
from unittest.mock import patch, MagicMock
from celestial_data import (
    GRAVITATIONAL_CONSTANT,
    SPEED_OF_LIGHT,
    get_celestial_body_data,
    get_all_solar_system_destinations,
    CELESTIAL_BODIES_DATA
)
from constants import G_GRAVITATIONAL, C_LIGHT_MPS # To compare against re-exported values

# --- Test Constants Re-exported from celestial_data ---
def test_gravitational_constant_re_export_value():
    """
    Verify that GRAVITATIONAL_CONSTANT re-exported by celestial_data
    matches the definitive G_GRAVITATIONAL from constants.
    """
    assert isinstance(GRAVITATIONAL_CONSTANT, float)
    assert GRAVITATIONAL_CONSTANT == pytest.approx(G_GRAVITATIONAL, rel=1e-10)

def test_speed_of_light_re_export_value():
    """
    Verify that SPEED_OF_LIGHT re-exported by celestial_data
    matches the definitive C_LIGHT_MPS from constants.
    """
    assert isinstance(SPEED_OF_LIGHT, float)
    assert SPEED_OF_LIGHT == pytest.approx(C_LIGHT_MPS, rel=1e-9)

# --- Test get_celestial_body_data function ---
def test_get_celestial_body_data_valid_body():
    """
    Test retrieval of data for a valid celestial body like 'Earth'.
    Ensures correct structure and type of returned data.
    """
    earth_data = get_celestial_body_data('Earth')
    assert earth_data is not None
    assert isinstance(earth_data, dict)
    assert earth_data['name'] == 'Earth'
    assert 'mass' in earth_data and isinstance(earth_data['mass'], float)
    assert 'radius' in earth_data and isinstance(earth_data['radius'], float)
    assert 'gravitational_parameter_mu' in earth_data and isinstance(earth_data['gravitational_parameter_mu'], float)
    assert earth_data['semi_major_axis_from_sun'] > 0 # Earth has a semi-major axis from the Sun

def test_get_celestial_body_data_case_insensitivity():
    """
    Verify that get_celestial_body_data handles case-insensitive body names.
    """
    mars_data_lower = get_celestial_body_data('mars')
    mars_data_upper = get_celestial_body_data('MARS')
    mars_data_mixed = get_celestial_body_data('MaRs')
    assert mars_data_lower is not None
    assert mars_data_lower == mars_data_upper
    assert mars_data_lower == mars_data_mixed
    assert mars_data_lower['name'] == 'Mars' # Check the canonical name

def test_get_celestial_body_data_non_existent_body():
    """
    Test that `None` is returned for a celestial body that does not exist in the data.
    """
    assert get_celestial_body_data('Pluto') is None
    assert get_celestial_body_data('Vulcan') is None

def test_get_celestial_body_data_invalid_input_type():
    """
    Test that `None` is returned for invalid input types (non-string).
    """
    assert get_celestial_body_data(123) is None
    assert get_celestial_body_data(None) is None
    assert get_celestial_body_data([]) is None

def test_get_celestial_body_data_sun_orbital_elements_expected_values():
    """
    Test that Sun's data is retrieved correctly and its orbital elements
    relative to itself (the central body) are zero or None as expected.
    """
    sun_data = get_celestial_body_data('Sun')
    assert sun_data is not None
    assert sun_data['name'] == 'Sun'
    assert sun_data.get('semi_major_axis_from_sun', 0) == 0 # Sun's semi-major axis from itself is 0
    # Other orbital elements relative to itself should ideally be absent or zero
    assert sun_data.get('orbital_period_days') is None or sun_data.get('orbital_period_days') == 0

# --- Test get_all_solar_system_destinations function ---
def test_get_all_solar_system_destinations_valid_output():
    """
    Test that `get_all_solar_system_destinations` returns a sorted list of destinations,
    excluding Earth, with correct structure.
    """
    destinations = get_all_solar_system_destinations()
    assert isinstance(destinations, list)
    assert len(destinations) > 0
    
    # Earth should not be in the list of destinations
    assert 'Earth' not in [d['name'] for d in destinations]
    
    # Check that the list is sorted by distance
    for i in range(len(destinations) - 1):
        assert destinations[i]['distance'] <= destinations[i+1]['distance']
        assert isinstance(destinations[i]['name'], str)
        assert isinstance(destinations[i]['distance'], (float, int))
        assert destinations[i]['distance'] >= 0 # Distances should be non-negative

    # Assuming Moon is closest and Neptune furthest based on typical solar system data
    # This test relies on the existing data in CELESTIAL_BODIES_DATA
    assert destinations[0]['name'] == 'Moon'
    assert destinations[-1]['name'] == 'Neptune'

@patch('celestial_data.get_celestial_body_data')
def test_get_all_solar_system_destinations_earth_data_missing(mock_get_celestial_body_data):
    """
    Test that `get_all_solar_system_destinations` returns `None`
    if Earth's data cannot be retrieved (e.g., it's missing or corrupted).
    """
    # Configure the mock to return valid data for others but None for Earth
    def side_effect(name):
        if name.lower() == 'earth':
            return None
        # Provide minimal data for other bodies for the sorting logic to proceed
        # These mock values need to be consistent with what CELESTIAL_BODIES_DATA would provide
        # for other bodies in a real scenario.
        # This setup needs to simulate the actual iteration process if CELESTIAL_BODIES_DATA is still used
        # or if get_celestial_body_data is called for each.
        # A simpler way is to mock CELESTIAL_BODIES_DATA directly to remove Earth, but
        # the summary specifically mentioned `get_celestial_body_data` for Earth.
        # Let's assume the function calls get_celestial_body_data for Earth explicitly
        # and also internally iterates CELESTIAL_BODIES_DATA.
        # For simplicity, if get_celestial_body_data('Earth') returns None,
        # the function should fail as per the summary.
        return CELESTIAL_BODIES_DATA.get(name.capitalize())

    mock_get_celestial_body_data.side_effect = side_effect

    destinations = get_all_solar_system_destinations()
    assert destinations is None
    mock_get_celestial_body_data.assert_any_call('Earth') # Ensure Earth's data was attempted to be fetched

def test_celestial_bodies_data_integrity():
    """
    Perform a general check on the structure and basic types/values of all
    entries in `CELESTIAL_BODIES_DATA`.
    """
    for body_name, data in CELESTIAL_BODIES_DATA.items():
        assert isinstance(body_name, str)
        assert isinstance(data, dict)
        assert 'name' in data and isinstance(data['name'], str)
        assert 'mass' in data and isinstance(data['mass'], (float, int)) and data['mass'] >= 0
        assert 'radius' in data and isinstance(data['radius'], (float, int)) and data['radius'] >= 0
        assert 'gravitational_parameter_mu' in data and isinstance(data['gravitational_parameter_mu'], (float, int)) and data['gravitational_parameter_mu'] >= 0

        # Check orbital elements if they exist (not all bodies have them relative to Sun)
        if 'semi_major_axis_from_sun' in data:
            assert isinstance(data['semi_major_axis_from_sun'], (float, int))
            assert data['semi_major_axis_from_sun'] >= 0
        if 'orbital_period_days' in data:
            assert isinstance(data['orbital_period_days'], (float, int))
            assert data['orbital_period_days'] >= 0

        # Average distance from Earth is crucial for destinations list
        if 'average_distance_from_earth' in data:
            assert isinstance(data['average_distance_from_earth'], (float, int))
            assert data['average_distance_from_earth'] >= 0