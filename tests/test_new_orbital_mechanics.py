import pytest
import math
from unittest.mock import patch, MagicMock
from typing import Optional, Dict, Any

# Import the functions to be tested from the actual module
from src.astrophysics.orbital_calculations import (
    calculate_orbital_velocity,
    calculate_escape_velocity,
    calculate_orbital_period,
)

# Import celestial data for mocking
from src.celestial_data import CELESTIAL_BODIES_DATA

# Define some test data based on CELESTIAL_BODIES_DATA for consistency
# Using approximations for simplicity in tests
SUN_MU = CELESTIAL_BODIES_DATA['SUN']['gravitational_parameter_mu']
EARTH_RADIUS = CELESTIAL_BODIES_DATA['EARTH']['radius_m']
EARTH_MU = CELESTIAL_BODIES_DATA['EARTH']['gravitational_parameter_mu']
EARTH_SEMI_MAJOR_AXIS = CELESTIAL_BODIES_DATA['EARTH']['semi_major_axis_m']
EARTH_ORBITAL_PERIOD_S = CELESTIAL_BODIES_DATA['EARTH']['orbital_period_s']

MARS_SEMI_MAJOR_AXIS = CELESTIAL_BODIES_DATA['MARS']['semi_major_axis_m']
MARS_ORBITAL_PERIOD_S = CELESTIAL_BODIES_DATA['MARS']['orbital_period_s']


# Mock function to simulate get_celestial_body_data behavior
def mock_get_celestial_body_data(body_name: str) -> Optional[Dict[str, Any]]:
    return CELESTIAL_BODIES_DATA.get(body_name.upper())


@patch('src.astrophysics.orbital_calculations.get_celestial_body_data', side_effect=mock_get_celestial_body_data)
class TestOrbitalVelocity:
    """
    Test suite for calculate_orbital_velocity.
    """

    @pytest.mark.parametrize(
        "body_name, distance, central_body_name, expected_velocity_m_s, tolerance",
        [
            # Earth around Sun (approximate circular orbit velocity at its semi-major axis)
            ("Earth", EARTH_SEMI_MAJOR_AXIS, "Sun", 29780, 100), # Actual ~29.78 km/s
            # Mars around Sun (approximate circular orbit velocity at its semi-major axis)
            ("Mars", MARS_SEMI_MAJOR_AXIS, "Sun", 24130, 100), # Actual ~24.13 km/s
            # A hypothetical body in circular orbit around Earth at its surface (very low orbit)
            ("Moon", EARTH_RADIUS + 1, "Earth", math.sqrt(EARTH_MU / (EARTH_RADIUS + 1)), 0.1),
        ]
    )
    def test_calculate_orbital_velocity_valid(self, mock_get_data, body_name, distance, central_body_name, expected_velocity_m_s, tolerance):
        """Tests calculate_orbital_velocity with valid inputs."""
        velocity = calculate_orbital_velocity(body_name, distance, central_body_name)
        assert math.isclose(velocity, expected_velocity_m_s, abs_tol=tolerance)

    @pytest.mark.parametrize(
        "body_name, distance, central_body_name",
        [
            (123, EARTH_SEMI_MAJOR_AXIS, "Sun"),
            ("Earth", "invalid", "Sun"),
            ("Earth", EARTH_SEMI_MAJOR_AXIS, 456),
        ]
    )
    def test_calculate_orbital_velocity_type_errors(
        self, mock_get_data, body_name, distance, central_body_name
    ):
        """Tests calculate_orbital_velocity with invalid input types."""
        with pytest.raises(TypeError):
            calculate_orbital_velocity(body_name, distance, central_body_name)

    @pytest.mark.parametrize(
        "body_name, distance, central_body_name, expected_error, error_message_part",
        [
            ("NonExistent", EARTH_SEMI_MAJOR_AXIS, "Sun", ValueError, "Celestial body data not found for orbiting body: NonExistent"),
            ("Earth", EARTH_SEMI_MAJOR_AXIS, "NonExistent", ValueError, "Celestial body data not found for central body: NonExistent"),
            ("Earth", 0.0, "Sun", ValueError, "Distance must be positive."),
            ("Earth", -100.0, "Sun", ValueError, "Distance must be positive."),
            # Test missing gravitational_parameter_mu for central body (mock it for Sun)
            ("Earth", EARTH_SEMI_MAJOR_AXIS, "Sun", ValueError, "gravitational_parameter_mu.*must be positive"),
            # Test missing semi_major_axis_m for orbiting body (mock it for Jupiter)
            ("Jupiter", EARTH_SEMI_MAJOR_AXIS, "Sun", ValueError, "semi_major_axis_m.*must be positive"),
            # Test for a scenario where `(2/r - 1/a)` becomes negative, leading to math.sqrt(negative)
            (
                "Earth",
                EARTH_SEMI_MAJOR_AXIS * 3,
                "Sun",
                ValueError,
                "physically impossible velocity",
            ),
        ],
    )
    def test_calculate_orbital_velocity_error_handling(
        self, mock_get_data, body_name, distance, central_body_name, expected_error, error_message_part
    ):
        """Tests calculate_orbital_velocity with invalid values or missing data."""
        if body_name.upper() == "JUPITER" and "semi_major_axis_m" in error_message_part:
            modified_jupiter_data = CELESTIAL_BODIES_DATA.get('JUPITER', {}).copy()
            modified_jupiter_data.pop('semi_major_axis_m', None)
            with patch.dict(CELESTIAL_BODIES_DATA, {'JUPITER': modified_jupiter_data}):
                with pytest.raises(expected_error, match=error_message_part):
                    calculate_orbital_velocity(body_name, distance, central_body_name)
        elif central_body_name.upper() == "SUN" and "gravitational_parameter_mu" in error_message_part:
            modified_sun_data = CELESTIAL_BODIES_DATA.get('SUN', {}).copy()
            modified_sun_data.pop('gravitational_parameter_mu', None)
            with patch.dict(CELESTIAL_BODIES_DATA, {'SUN': modified_sun_data}):
                with pytest.raises(expected_error, match=error_message_part):
                    calculate_orbital_velocity(body_name, distance, central_body_name)
        else:
            with pytest.raises(expected_error, match=error_message_part):
                calculate_orbital_velocity(body_name, distance, central_body_name)


@patch('src.astrophysics.orbital_calculations.get_celestial_body_data', side_effect=mock_get_celestial_body_data)
class TestEscapeVelocity:
    """
    Test suite for calculate_escape_velocity.
    """

    @pytest.mark.parametrize(
        "body_name, distance, expected_velocity_m_s, tolerance",
        [
            # Earth escape velocity at its surface
            ("Earth", EARTH_RADIUS, 11186, 10), # Actual ~11.186 km/s
            # Moon escape velocity at its surface
            ("Moon", CELESTIAL_BODIES_DATA['MOON']['radius_m'], math.sqrt(2 * CELESTIAL_BODIES_DATA['MOON']['gravitational_parameter_mu'] / CELESTIAL_BODIES_DATA['MOON']['radius_m']), 0.1),
        ]
    )
    def test_calculate_escape_velocity_valid(self, mock_get_data, body_name, distance, expected_velocity_m_s, tolerance):
        """Tests calculate_escape_velocity with valid inputs."""
        velocity = calculate_escape_velocity(body_name, distance)
        assert math.isclose(velocity, expected_velocity_m_s, abs_tol=tolerance)

    @pytest.mark.parametrize(
        "body_name, distance",
        [
            (123, EARTH_RADIUS),
            ("Earth", "invalid"),
        ]
    )
    def test_calculate_escape_velocity_type_errors(
        self, mock_get_data, body_name, distance
    ):
        """Tests calculate_escape_velocity with invalid input types."""
        with pytest.raises(TypeError):
            calculate_escape_velocity(body_name, distance)

    @pytest.mark.parametrize(
        "body_name, distance, expected_error, error_message_part",
        [
            ("NonExistent", EARTH_RADIUS, ValueError, "Celestial body data not found for: NonExistent"),
            ("Earth", 0.0, ValueError, "Distance must be positive."),
            ("Earth", -100.0, ValueError, "Distance must be positive."),
            # Test missing gravitational_parameter_mu for the body (mock it for Venus)
            ("Venus", CELESTIAL_BODIES_DATA['VENUS']['radius_m'], ValueError, "gravitational_parameter_mu.*must be positive"),
        ]
    )
    def test_calculate_escape_velocity_error_handling(
        self, mock_get_data, body_name, distance, expected_error, error_message_part
    ):
        """Tests calculate_escape_velocity with invalid values or missing data."""
        if body_name.upper() == "VENUS" and "gravitational_parameter_mu" in error_message_part:
            modified_venus_data = CELESTIAL_BODIES_DATA.get('VENUS', {}).copy()
            modified_venus_data.pop('gravitational_parameter_mu', None)
            with patch.dict(CELESTIAL_BODIES_DATA, {'VENUS': modified_venus_data}):
                with pytest.raises(expected_error, match=error_message_part):
                    calculate_escape_velocity(body_name, distance)
        else:
            with pytest.raises(expected_error, match=error_message_part):
                calculate_escape_velocity(body_name, distance)


@patch('src.astrophysics.orbital_calculations.get_celestial_body_data', side_effect=mock_get_celestial_body_data)
class TestOrbitalPeriod:
    """
    Test suite for calculate_orbital_period.
    """

    @pytest.mark.parametrize(
        "body_name, central_body_name, expected_period_s, tolerance",
        [
            # Earth around Sun
            ("Earth", "Sun", EARTH_ORBITAL_PERIOD_S, 10000), # ~31,557,600 s (1 year)
            # Mars around Sun
            ("Mars", "Sun", MARS_ORBITAL_PERIOD_S, 10000), # ~59,354,294 s (1.88 years)
        ]
    )
    def test_calculate_orbital_period_valid(
        self, mock_get_data, body_name, central_body_name, expected_period_s, tolerance
    ):
        """Tests calculate_orbital_period with valid inputs."""
        period = calculate_orbital_period(body_name, central_body_name)
        assert math.isclose(period, expected_period_s, abs_tol=tolerance)

    @pytest.mark.parametrize(
        "body_name, central_body_name",
        [
            (123, "Sun"),
            ("Earth", 456),
        ]
    )
    def test_calculate_orbital_period_type_errors(
        self, mock_get_data, body_name, central_body_name
    ):
        """Tests calculate_orbital_period with invalid input types."""
        with pytest.raises(TypeError):
            calculate_orbital_period(body_name, central_body_name)

    @pytest.mark.parametrize(
        "body_name, central_body_name, expected_error, error_message_part",
        [
            ("NonExistent", "Sun", ValueError, "Celestial body data not found for orbiting body: NonExistent"),
            ("Earth", "NonExistent", ValueError, "Celestial body data not found for central body: NonExistent"),
            # Test missing gravitational_parameter_mu for central body (mock it for Sun)
            ("Earth", "Sun", ValueError, "gravitational_parameter_mu.*must be positive"),
            # Test missing semi_major_axis_m for orbiting body (mock it for Jupiter)
            ("Jupiter", "Sun", ValueError, "semi_major_axis_m.*must be positive"),
        ]
    )
    def test_calculate_orbital_period_error_handling(
        self, mock_get_data, body_name, central_body_name, expected_error, error_message_part
    ):
        """Tests calculate_orbital_period with invalid values or missing data."""
        if body_name.upper() == "JUPITER" and "semi_major_axis_m" in error_message_part:
            modified_jupiter_data = CELESTIAL_BODIES_DATA.get('JUPITER', {}).copy()
            modified_jupiter_data.pop('semi_major_axis_m', None)
            with patch.dict(CELESTIAL_BODIES_DATA, {'JUPITER': modified_jupiter_data}):
                with pytest.raises(expected_error, match=error_message_part):
                    calculate_orbital_period(body_name, central_body_name)
        elif central_body_name.upper() == "SUN" and "gravitational_parameter_mu" in error_message_part:
            modified_sun_data = CELESTIAL_BODIES_DATA.get('SUN', {}).copy()
            modified_sun_data.pop('gravitational_parameter_mu', None)
            with patch.dict(CELESTIAL_BODIES_DATA, {'SUN': modified_sun_data}):
                with pytest.raises(expected_error, match=error_message_part):
                    calculate_orbital_period(body_name, central_body_name)
        else:
            with pytest.raises(expected_error, match=error_message_part):
                calculate_orbital_period(body_name, central_body_name)