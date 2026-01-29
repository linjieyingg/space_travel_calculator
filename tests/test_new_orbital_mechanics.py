import pytest
import math
from unittest.mock import patch

# Import the functions to be tested from the actual module
from src.astrophysics.orbital_calculations import (
    calculate_orbital_velocity,
    calculate_escape_velocity,
    calculate_orbital_period,
)

# Mock data for celestial bodies for testing purposes
MOCK_CELESTIAL_DATA = {
    "SUN": {
        "name": "Sun",
        "mass_kg": 1.989e30,
        "radius_m": 6.957e8,
        "gravitational_parameter_mu": 1.32712440018e20,  # G * M_sun
        "semi_major_axis_m": 0.0,  # Sun is the central body, its 'semi_major_axis_m' is not applicable for orbiting itself
        "orbital_period_s": 0.0,
    },
    "EARTH": {
        "name": "Earth",
        "mass_kg": 5.972e24,
        "radius_m": 6.371e6,
        "gravitational_parameter_mu": 3.986004418e14,  # G * M_earth
        "semi_major_axis_m": 1.495978707e11,  # ~1 AU, Earth's avg distance from Sun
        "orbital_period_s": 31557600.0,  # ~1 year, Earth's orbital period around Sun
    },
    "MOON": {
        "name": "Moon",
        "mass_kg": 7.342e22,
        "radius_m": 1.7374e6,
        "gravitational_parameter_mu": 4.9048695e12,  # G * M_moon
        "semi_major_axis_m": 3.844e8,  # Avg distance of Moon from Earth
        "orbital_period_s": 2360592.0,  # ~27.3 days, Moon's orbital period around Earth
    },
    "MARS": {
        "name": "Mars",
        "mass_kg": 6.4171e23,
        "radius_m": 3.3895e6,
        "gravitational_parameter_mu": 4.282837e13,
        "semi_major_axis_m": 2.279392e11,  # ~1.52 AU, Mars's avg distance from Sun
        "orbital_period_s": 59354000.0,  # ~687 days, Mars's orbital period around Sun
    },
    # Test bodies for specific error cases (missing/invalid data)
    "BROKEN_BODY_MISSING_MU": {
        "name": "Broken Body Missing Mu",
        "semi_major_axis_m": 1.0e11,
    },  # Missing 'gravitational_parameter_mu'
    "BROKEN_BODY_MISSING_A": {
        "name": "Broken Body Missing A",
        "gravitational_parameter_mu": 1.0e10,
    },  # Missing 'semi_major_axis_m'
    "BODY_WITH_ZERO_MU": {
        "name": "Body with Zero Mu",
        "gravitational_parameter_mu": 0.0,
        "semi_major_axis_m": 1.0e11,
    },
    "BODY_WITH_ZERO_A": {
        "name": "Body with Zero A",
        "gravitational_parameter_mu": 1.0e10,
        "semi_major_axis_m": 0.0,
    },
    "BODY_WITH_NEGATIVE_MU": {
        "name": "Body with Negative Mu",
        "gravitational_parameter_mu": -1.0e10,
        "semi_major_axis_m": 1.0e11,
    },
    "BODY_WITH_NEGATIVE_A": {
        "name": "Body with Negative A",
        "gravitational_parameter_mu": 1.0e10,
        "semi_major_axis_m": -1.0e11,
    },
}


# Fixture to mock get_celestial_body_data for all tests in this module
@pytest.fixture
def mock_get_celestial_body_data():
    with patch(
        "src.astrophysics.orbital_calculations.get_celestial_body_data"
    ) as mock_func:
        mock_func.side_effect = lambda body_name: MOCK_CELESTIAL_DATA.get(
            body_name.upper()
        )
        yield mock_func


class TestNewOrbitalMechanics:
    """
    Test suite for new orbital mechanics functions:
    calculate_orbital_velocity, calculate_escape_velocity, calculate_orbital_period.
    """

    # --- Tests for calculate_orbital_velocity ---
    @pytest.mark.parametrize(
        "body, distance, central_body, expected_velocity, tolerance",
        [
            # Earth around Sun (approximate circular velocity at 1 AU)
            (
                "Earth",
                MOCK_CELESTIAL_DATA["EARTH"]["semi_major_axis_m"],
                "Sun",
                29782.87,
                0.01,
            ),
            # Moon around Earth (approximate circular velocity at its semi-major axis)
            (
                "Moon",
                MOCK_CELESTIAL_DATA["MOON"]["semi_major_axis_m"],
                "Earth",
                1018.33,
                0.01,
            ),
        ],
    )
    def test_calculate_orbital_velocity_success(
        self, mock_get_celestial_body_data, body, distance, central_body, expected_velocity, tolerance
    ):
        """Tests calculate_orbital_velocity with valid inputs."""
        velocity = calculate_orbital_velocity(body, distance, central_body)
        assert math.isclose(velocity, expected_velocity, rel_tol=tolerance)

    @pytest.mark.parametrize(
        "body, distance, central_body, expected_error, error_message_part",
        [
            ("Earth", "invalid", "Sun", TypeError, "Distance must be a numeric type."),
            (123, MOCK_CELESTIAL_DATA["EARTH"]["semi_major_axis_m"], "Sun", TypeError, "Body names must be strings."),
            ("Earth", MOCK_CELESTIAL_DATA["EARTH"]["semi_major_axis_m"], 456, TypeError, "Body names must be strings."),
        ],
    )
    def test_calculate_orbital_velocity_type_errors(
        self, mock_get_celestial_body_data, body, distance, central_body, expected_error, error_message_part
    ):
        """Tests calculate_orbital_velocity with invalid input types."""
        with pytest.raises(expected_error, match=error_message_part):
            calculate_orbital_velocity(body, distance, central_body)

    @pytest.mark.parametrize(
        "body, distance, central_body, expected_error, error_message_part",
        [
            ("Earth", 0.0, "Sun", ValueError, "Distance must be positive."),
            ("Earth", -100.0, "Sun", ValueError, "Distance must be positive."),
            (
                "NonExistent",
                MOCK_CELESTIAL_DATA["EARTH"]["semi_major_axis_m"],
                "Sun",
                ValueError,
                "Celestial body data not found for orbiting body: NonExistent",
            ),
            (
                "Earth",
                MOCK_CELESTIAL_DATA["EARTH"]["semi_major_axis_m"],
                "NonExistent",
                ValueError,
                "Celestial body data not found for central body: NonExistent",
            ),
            (
                "Earth",
                MOCK_CELESTIAL_DATA["EARTH"]["semi_major_axis_m"],
                "BROKEN_BODY_MISSING_MU",
                ValueError,
                "Missing or invalid 'gravitational_parameter_mu' for central body BROKEN_BODY_MISSING_MU",
            ),
            (
                "BROKEN_BODY_MISSING_A",
                MOCK_CELESTIAL_DATA["EARTH"]["semi_major_axis_m"],
                "Sun",
                ValueError,
                "Missing or invalid 'semi_major_axis_m' for orbiting body BROKEN_BODY_MISSING_A",
            ),
            (
                "Earth",
                MOCK_CELESTIAL_DATA["EARTH"]["semi_major_axis_m"],
                "BODY_WITH_ZERO_MU",
                ValueError,
                "must be positive",
            ),
            (
                "BODY_WITH_ZERO_A",
                MOCK_CELESTIAL_DATA["EARTH"]["semi_major_axis_m"],
                "Sun",
                ValueError,
                "must be positive",
            ),
            # Test for a scenario where `(2/r - 1/a)` becomes negative, leading to math.sqrt(negative)
            (
                "Earth",
                MOCK_CELESTIAL_DATA["EARTH"]["semi_major_axis_m"] * 3,
                "Sun",
                ValueError,
                "physically impossible velocity",
            ),
            (
                "BODY_WITH_NEGATIVE_A",
                MOCK_CELESTIAL_DATA["EARTH"]["semi_major_axis_m"],
                "Sun",
                ValueError,
                "must be positive",
            ),
            (
                "Earth",
                MOCK_CELESTIAL_DATA["EARTH"]["semi_major_axis_m"],
                "BODY_WITH_NEGATIVE_MU",
                ValueError,
                "must be positive",
            ),
        ],
    )
    def test_calculate_orbital_velocity_value_errors(
        self, mock_get_celestial_body_data, body, distance, central_body, expected_error, error_message_part
    ):
        """Tests calculate_orbital_velocity with invalid values or missing data."""
        with pytest.raises(expected_error, match=error_message_part):
            calculate_orbital_velocity(body, distance, central_body)

    # --- Tests for calculate_escape_velocity ---
    @pytest.mark.parametrize(
        "body, distance, expected_velocity, tolerance",
        [
            # From Earth surface
            ("Earth", MOCK_CELESTIAL_DATA["EARTH"]["radius_m"], 11184.8, 0.01),
            # From Sun surface
            ("Sun", MOCK_CELESTIAL_DATA["SUN"]["radius_m"], 617743.0, 0.01),
            # From Moon surface
            ("Moon", MOCK_CELESTIAL_DATA["MOON"]["radius_m"], 2379.7, 0.01),
        ],
    )
    def test_calculate_escape_velocity_success(
        self, mock_get_celestial_body_data, body, distance, expected_velocity, tolerance
    ):
        """Tests calculate_escape_velocity with valid inputs."""
        velocity = calculate_escape_velocity(body, distance)
        assert math.isclose(velocity, expected_velocity, rel_tol=tolerance)

    @pytest.mark.parametrize(
        "body, distance, expected_error, error_message_part",
        [
            ("Earth", "invalid", TypeError, "Distance must be a numeric type."),
            (123, MOCK_CELESTIAL_DATA["EARTH"]["radius_m"], TypeError, "Body name must be a string."),
        ],
    )
    def test_calculate_escape_velocity_type_errors(
        self, mock_get_celestial_body_data, body, distance, expected_error, error_message_part
    ):
        """Tests calculate_escape_velocity with invalid input types."""
        with pytest.raises(expected_error, match=error_message_part):
            calculate_escape_velocity(body, distance)

    @pytest.mark.parametrize(
        "body, distance, expected_error, error_message_part",
        [
            ("Earth", 0.0, ValueError, "Distance must be positive."),
            ("Earth", -100.0, ValueError, "Distance must be positive."),
            ("NonExistent", MOCK_CELESTIAL_DATA["EARTH"]["radius_m"], ValueError, "Celestial body data not found for: NonExistent"),
            (
                "BROKEN_BODY_MISSING_MU",
                MOCK_CELESTIAL_DATA["EARTH"]["radius_m"],
                ValueError,
                "Missing or invalid 'gravitational_parameter_mu' for BROKEN_BODY_MISSING_MU",
            ),
            ("BODY_WITH_ZERO_MU", MOCK_CELESTIAL_DATA["EARTH"]["radius_m"], ValueError, "must be positive"),
            ("BODY_WITH_NEGATIVE_MU", MOCK_CELESTIAL_DATA["EARTH"]["radius_m"], ValueError, "must be positive"),
        ],
    )
    def test_calculate_escape_velocity_value_errors(
        self, mock_get_celestial_body_data, body, distance, expected_error, error_message_part
    ):
        """Tests calculate_escape_velocity with invalid values or missing data."""
        with pytest.raises(expected_error, match=error_message_part):
            calculate_escape_velocity(body, distance)

    # --- Tests for calculate_orbital_period ---
    @pytest.mark.parametrize(
        "body, central_body, expected_period, tolerance",
        [
            # Earth around Sun (~1 year)
            ("Earth", "Sun", 31557600.0, 0.01),
            # Moon around Earth (~27.3 days)
            ("Moon", "Earth", 2360592.0, 0.01),
        ],
    )
    def test_calculate_orbital_period_success(
        self, mock_get_celestial_body_data, body, central_body, expected_period, tolerance
    ):
        """Tests calculate_orbital_period with valid inputs."""
        period = calculate_orbital_period(body, central_body)
        assert math.isclose(period, expected_period, rel_tol=tolerance)

    @pytest.mark.parametrize(
        "body, central_body, expected_error, error_message_part",
        [
            (123, "Sun", TypeError, "Body names must be strings."),
            ("Earth", 456, TypeError, "Body names must be strings."),
        ],
    )
    def test_calculate_orbital_period_type_errors(
        self, mock_get_celestial_body_data, body, central_body, expected_error, error_message_part
    ):
        """Tests calculate_orbital_period with invalid input types."""
        with pytest.raises(expected_error, match=error_message_part):
            calculate_orbital_period(body, central_body)

    @pytest.mark.parametrize(
        "body, central_body, expected_error, error_message_part",
        [
            (
                "NonExistent",
                "Sun",
                ValueError,
                "Celestial body data not found for orbiting body: NonExistent",
            ),
            (
                "Earth",
                "NonExistent",
                ValueError,
                "Celestial body data not found for central body: NonExistent",
            ),
            (
                "Earth",
                "BROKEN_BODY_MISSING_MU",
                ValueError,
                "Missing or invalid 'gravitational_parameter_mu' for central body BROKEN_BODY_MISSING_MU",
            ),
            (
                "BROKEN_BODY_MISSING_A",
                "Sun",
                ValueError,
                "Missing or invalid 'semi_major_axis_m' for orbiting body BROKEN_BODY_MISSING_A",
            ),
            ("Earth", "BODY_WITH_ZERO_MU", ValueError, "must be positive"),
            ("BODY_WITH_ZERO_A", "Sun", ValueError, "must be positive"),
            ("Earth", "BODY_WITH_NEGATIVE_MU", ValueError, "must be positive"),
            ("BODY_WITH_NEGATIVE_A", "Sun", ValueError, "must be positive"),
        ],
    )
    def test_calculate_orbital_period_value_errors(
        self, mock_get_celestial_body_data, body, central_body, expected_error, error_message_part
    ):
        """Tests calculate_orbital_period with invalid values or missing data."""
        with pytest.raises(expected_error, match=error_message_part):
            calculate_orbital_period(body, central_body)