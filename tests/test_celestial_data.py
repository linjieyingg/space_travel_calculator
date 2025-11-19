import pytest
from celestial_data import get_celestial_body_data, get_all_solar_system_destinations, CELESTIAL_BODIES_DATA

def test_get_celestial_body_data_valid_input():
    data = get_celestial_body_data('earth')
    assert data is not None
    assert isinstance(data, dict)
    assert 'mass' in data
    assert 'gravitational_parameter_mu' in data

def test_get_celestial_body_data_invalid_input():
    assert get_celestial_body_data(123) is None
    assert get_celestial_body_data(None) is None

def test_get_celestial_body_data_non_existent():
    assert get_celestial_body_data('plutonot') is None

def test_get_all_solar_system_destinations():
    destinations = get_all_solar_system_destinations()
    assert isinstance(destinations, list)
    assert len(destinations) > 0
    assert all(isinstance(d, dict) for d in destinations)
    assert all('name' in d and 'distance' in d for d in destinations)
    assert 'Earth' not in [d['name'] for d in destinations]
    # Check if sorted by distance
    assert all(destinations[i]['distance'] <= destinations[i+1]['distance'] for i in range(len(destinations) - 1))
    # Check for specific body (e.g., Mars distance calculation sanity)
    mars_data = next((d for d in destinations if d['name'] == 'Mars'), None)
    assert mars_data is not None
    # Add more specific checks if CELESTIAL_BODIES_DATA becomes more stable in test environment