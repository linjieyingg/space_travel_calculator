import pytest
from celestial_data import get_celestial_body_data, get_all_solar_system_destinations, GRAVITATIONAL_CONSTANT, SPEED_OF_LIGHT

def test_gravitational_constant_exists():
    """Test that GRAVITATIONAL_CONSTANT is defined, is a float, and is positive."""
    assert GRAVITATIONAL_CONSTANT is not None
    assert isinstance(GRAVITATIONAL_CONSTANT, float)
    assert GRAVITATIONAL_CONSTANT > 0

def test_speed_of_light_exists():
    """Test that SPEED_OF_LIGHT is defined, is a float, and is positive."""
    assert SPEED_OF_LIGHT is not None
    assert isinstance(SPEED_OF_LIGHT, float)
    assert SPEED_OF_LIGHT > 0

def test_get_celestial_body_data_earth():
    """Test retrieving data for 'Earth' and verify its basic structure."""
    data = get_celestial_body_data("Earth")
    assert data is not None
    assert isinstance(data, dict)
    assert data.get("name") == "Earth"
    assert "mass" in data
    assert "radius" in data
    assert "gravitational_parameter_mu" in data
    assert "semi_major_axis_from_sun" in data
    assert data["gravitational_parameter_mu"] > 0

def test_get_celestial_body_data_sun():
    """Test retrieving data for 'Sun' and verify its basic structure."""
    data = get_celestial_body_data("Sun")
    assert data is not None
    assert isinstance(data, dict)
    assert data.get("name") == "Sun"
    assert "mass" in data
    assert "radius" in data
    assert "gravitational_parameter_mu" in data
    assert data["gravitational_parameter_mu"] > 0

def test_get_celestial_body_data_case_insensitivity():
    """Test that body names are handled case-insensitively."""
    data_lower = get_celestial_body_data("mars")
    data_upper = get_celestial_body_data("MARS")
    assert data_lower is not None
    assert data_lower == data_upper
    assert data_lower.get("name") == "Mars"

def test_get_celestial_body_data_non_existent():
    """Test retrieving data for a non-existent celestial body returns None."""
    data = get_celestial_body_data("NonExistentBody")
    assert data is None

@pytest.mark.parametrize("invalid_input", [
    123,
    None,
    [],
    {}
])
def test_get_celestial_body_data_invalid_input(invalid_input):
    """Test retrieving data with invalid input types returns None."""
    assert get_celestial_body_data(invalid_input) is None

def test_get_all_solar_system_destinations():
    """Test get_all_solar_system_destinations for proper output format, exclusion of Earth, and sorting."""
    destinations = get_all_solar_system_destinations()
    assert isinstance(destinations, list)
    assert len(destinations) > 0
    assert all(isinstance(d, dict) for d in destinations)
    assert all("name" in d and "distance" in d for d in destinations)
    
    # Earth should not be in the destinations list
    assert "Earth" not in [d["name"] for d in destinations]
    
    # Check if sorted by distance (ascending)
    distances = [d["distance"] for d in destinations]
    assert distances == sorted(distances)

    # Check for specific known bodies (e.g., Mars, Jupiter, Moon)
    body_names = [d["name"] for d in destinations]
    assert "Mars" in body_names
    assert "Jupiter" in body_names
    assert "Moon" in body_names
    assert "Sun" not in body_names # Sun is not a destination from Earth, but the central body for SMA calculations

def test_get_all_solar_system_destinations_with_earth_data_integrity():
    """
    Test that the function correctly calculates distances relative to Earth's semi-major axis,
    and includes Moon's average distance from Earth if available.
    """
    destinations = get_all_solar_system_destinations()
    earth_data = get_celestial_body_data("Earth")
    earth_sma = earth_data["semi_major_axis_from_sun"]

    # Verify Mars distance calculation
    mars_data = get_celestial_body_data("Mars")
    mars_sma = mars_data["semi_major_axis_from_sun"]
    expected_mars_distance = abs(mars_sma - earth_sma)
    
    mars_destination = next((d for d in destinations if d["name"] == "Mars"), None)
    assert mars_destination is not None
    assert pytest.approx(mars_destination["distance"], rel=1e-6) == expected_mars_distance

    # Verify Moon's distance (which is explicitly defined, not calculated from SMA)
    moon_data = get_celestial_body_data("Moon")
    expected_moon_distance = moon_data["average_distance_from_earth"]
    
    moon_destination = next((d for d in destinations if d["name"] == "Moon"), None)
    assert moon_destination is not None
    assert pytest.approx(moon_destination["distance"], rel=1e-6) == expected_moon_distance