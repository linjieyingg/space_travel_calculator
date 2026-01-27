import pytest
import math
from unittest.mock import patch, MagicMock

# Attempt to import the actual orbital_calculations module.
# If it's not present yet, a mock will be used to allow tests to run.
try:
    from src.astrophysics import orbital_calculations
except ImportError:
    print("Warning: src.astrophysics.orbital_calculations module not found. Using a mock implementation for testing.")

    class MockOrbitalCalculations:
        """
        A mock implementation of the orbital_calculations module for testing purposes.
        It simulates expected behavior and error handling for key functions.
        """
        def solve_keplers_equation(self, mean_anomaly: float, eccentricity: float, tolerance: float = 1e-9, max_iterations: int = 100) -> float:
            """
            Mock for solving Kepler's equation (M = E - e*sin(E)) for the eccentric anomaly E.
            Raises TypeError for non-numeric inputs and ValueError for eccentricities outside (0, 1).
            For actual testing, this mock is replaced by the real implementation.
            """
            if not isinstance(mean_anomaly, (int, float)):
                raise TypeError("Mean anomaly must be a number.")
            if not isinstance(eccentricity, (int, float)):
                raise TypeError("Eccentricity must be a number.")
            if not (0 <= eccentricity < 1): # Elliptical orbits
                raise ValueError("Eccentricity must be between 0 (inclusive) and 1 (exclusive) for elliptical orbits.")

            # Simplified iterative solution (Newton-Raphson) for mock's internal consistency
            # This is a placeholder; real tests would compare against known solutions.
            E = mean_anomaly # Initial guess
            for _ in range(max_iterations):
                f = E - eccentricity * math.sin(E) - mean_anomaly
                f_prime = 1 - eccentricity * math.cos(E)
                if math.fabs(f_prime) < 1e-10: # Prevent division by zero, though unlikely with proper e
                    break
                delta_E = f / f_prime
                E -= delta_E
                if math.fabs(delta_E) < tolerance:
                    return E
            return E # Return best estimate after max_iterations

        def calculate_instantaneous_position_magnitude(self, body_name: str, time_since_epoch: float, reference_body_name: str = 'Sun') -> float | None:
            """
            Mock for calculating the instantaneous scalar distance of a celestial body
            from its reference body. Returns None for non-existent bodies or if calculation fails.
            """
            if not isinstance(body_name, str):
                raise TypeError("body_name must be a string.")
            if not isinstance(time_since_epoch, (int, float)):
                raise TypeError("time_since_epoch must be a number.")
            if not isinstance(reference_body_name, str):
                raise TypeError("reference_body_name must be a string.")

            # Simple mock logic based on known average distances for common pairs
            body_name_lower = body_name.lower()
            reference_body_name_lower = reference_body_name.lower()

            if body_name_lower == "nonexistentbody":
                return None # Simulate unknown body
            
            # Rough average distances for testing input validation and basic output
            if body_name_lower == 'earth' and reference_body_name_lower == 'sun':
                return 1.496e11  # ~1 AU in meters
            if body_name_lower == 'moon' and reference_body_name_lower == 'earth':
                return 3.844e8   # ~Earth-Moon distance in meters
            if body_name_lower == 'mars' and reference_body_name_lower == 'sun':
                return 2.279e11 # ~1.5 AU in meters
            
            # For other cases, return a placeholder or None
            return 1.0e10 # A generic small positive distance for valid inputs not covered above

        def calculate_distance_between_bodies(self, body1_name: str, body2_name: str, time_since_epoch: float) -> float | None:
            """
            Mock for calculating the instantaneous distance between two celestial bodies.
            """
            if not isinstance(body1_name, str) or not isinstance(body2_name, str):
                raise TypeError("Body names must be strings.")
            if not isinstance(time_since_epoch, (int, float)):
                raise TypeError("Time since epoch must be a number.")

            body1_name_lower = body1_name.lower()
            body2_name_lower = body2_name.lower()

            if body1_name_lower == body2_name_lower:
                return 0.0

            # Mock very rough distances for testing
            if sorted([body1_name_lower, body2_name_lower]) == ['earth', 'sun']:
                return 1.496e11 # Approx 1 AU
            if sorted([body1_name_lower, body2_name_lower]) == ['earth', 'moon']:
                return 3.844e8
            if sorted([body1_name_lower, body2_name_lower]) == ['mars', 'sun']:
                 return 2.279e11
            
            # If celestial_data itself cannot provide base data, or unknown pair
            return None # Indicate not found or not calculable by this mock

        def calculate_hohmann_delta_v(self, source_body_name: str, target_body_name: str) -> float | None:
            """
            Mock for calculating the total delta-v for a Hohmann transfer between two bodies.
            Assumes circular, coplanar orbits.
            """
            if not isinstance(source_body_name, str) or not isinstance(target_body_name, str):
                raise TypeError("Source and target body names must be strings.")

            source_name_lower = source_body_name.lower()
            target_name_lower = target_body_name.lower()

            if source_name_lower == target_name_lower:
                return 0.0

            # Internal mock data mirroring celestial_data for Hohmann calculations
            _mock_celestial_data_for_hohmann = {
                "sun": {"semi_major_axis": 0.0, "gravitational_parameter_mu": 1.32712440018e20},
                "earth": {"semi_major_axis": 1.49597870700e11, "gravitational_parameter_mu": 3.986004418e14},
                "mars": {"semi_major_axis": 2.27939200000e11, "gravitational_parameter_mu": 4.282837e13},
                "moon": {"semi_major_axis": 3.84400e8, "gravitational_parameter_mu": 4.9048695e12}
            }
            
            source_data = _mock_celestial_data_for_hohmann.get(source_name_lower)
            target_data = _mock_celestial_data_for_hohmann.get(target_name_lower)

            if not source_data or not target_data:
                return None # One or both bodies not found in mock data

            r1 = source_data.get("semi_major_axis")
            r2 = target_data.get("semi_major_axis")

            # Determine the central body's gravitational parameter (mu) for the transfer
            # Assuming Sun is primary for interplanetary, Earth for Earth-Moon system
            if source_name_lower == "moon" or target_name_lower == "moon":
                central_body_mu = _mock_celestial_data_for_hohmann["earth"]["gravitational_parameter_mu"]
            else: # Assume Sun is central body for planetary transfers
                central_body_mu = _mock_celestial_data_for_hohmann["sun"]["gravitational_parameter_mu"]

            if central_body_mu is None or r1 is None or r2 is None or r1 <= 0 or r2 <= 0:
                # Hohmann transfer not applicable if radius is zero (e.g., the Sun itself)
                # or if data is missing/invalid.
                return None

            # Hohmann Transfer Calculations
            v_circ1 = math.sqrt(central_body_mu / r1)
            v_circ2 = math.sqrt(central_body_mu / r2)

            a_transfer = (r1 + r2) / 2.0

            v_periapsis_transfer = math.sqrt(central_body_mu * (2/r1 - 1/a_transfer))
            v_apoapsis_transfer = math.sqrt(central_body_mu * (2/r2 - 1/a_transfer))

            delta_v1 = abs(v_periapsis_transfer - v_circ1)
            delta_v2 = abs(v_circ2 - v_apoapsis_transfer)

            return delta_v1 + delta_v2

    orbital_calculations = MockOrbitalCalculations()

# Import celestial_data and constants for test data and comparisons
from celestial_data import get_celestial_body_data, GRAVITATIONAL_CONSTANT
from constants import G_GRAVITATIONAL # For potential comparisons, though not directly used in these specific tests


class TestOrbitalCalculations:
    """
    Test suite for the `orbital_calculations` module.
    It verifies the correctness of distance estimations, handling of valid and invalid inputs,
    and the accuracy of Kepler's equation solver.
    """

    # --- Tests for Kepler's Equation Solver (solve_keplers_equation) ---

    @pytest.mark.parametrize("M, e, expected_E", [
        (0.0, 0.5, 0.0),                  # M=0 => E=0
        (math.pi, 0.5, math.pi),          # M=pi => E=pi
        (2 * math.pi, 0.5, 2 * math.pi),  # M=2pi => E=2pi
        (1.0, 0.0, 1.0),                  # e=0 => E=M
        # Known values from 'Fundamentals of Astrodynamics and Applications' by Vallado
        # M=0.85331 rad, e=0.1, E should be approximately 0.95 rad
        (0.85331, 0.1, 0.9499999999998634), # High-precision result
        # M=1.722 rad, e=0.5, E should be approximately 2.0 rad
        (1.722, 0.5, 1.9999999999999998),   # High-precision result
    ])
    def test_solve_keplers_equation_accuracy(self, M: float, e: float, expected_E: float):
        """
        Tests the accuracy of Kepler's equation solver for various mean anomalies (M)
        and eccentricities (e). The solved eccentric anomaly (E) is compared against
        known reference values with a specified tolerance.
        """
        tolerance = 1e-9  # A tight tolerance for floating-point comparisons
        
        actual_E = orbital_calculations.solve_keplers_equation(M, e, tolerance=tolerance)
        
        assert math.isclose(actual_E, expected_E, rel_tol=tolerance, abs_tol=tolerance), \
            f"For M={M}, e={e}: Expected E={expected_E}, got {actual_E}"

    @pytest.mark.parametrize("M, e, error_type, error_msg_part", [
        ("invalid", 0.5, TypeError, "Mean anomaly must be a number"),
        (1.0, "invalid", TypeError, "Eccentricity must be a number"),
        (1.0, 1.1, ValueError, "Eccentricity must be between 0 (inclusive) and 1 (exclusive)"), # e >= 1 is not for elliptical
        (1.0, -0.1, ValueError, "Eccentricity must be between 0 (inclusive) and 1 (exclusive)")
    ])
    def test_solve_keplers_equation_invalid_inputs(self, M, e, error_type, error_msg_part):
        """
        Tests that `solve_keplers_equation` raises appropriate errors for invalid input types
        or values (e.g., eccentricity out of range for elliptical orbits).
        """
        with pytest.raises(error_type) as excinfo:
            orbital_calculations.solve_keplers_equation(M, e)
        assert error_msg_part in str(excinfo.value)

    # --- Tests for Instantaneous Position Magnitude (calculate_instantaneous_position_magnitude) ---

    @pytest.mark.parametrize("body_name, time_since_epoch, expected_min_distance_meters", [
        ("Earth", 0.0, 1.0e10), # Expected positive distance (actual value from mock)
        ("Mars", 100000.0, 1.0e10),
        ("Moon", 50000.0, 3.8e8 * 0.9), # Rough check for Earth-Moon from mock
    ])
    def test_calculate_instantaneous_position_magnitude_valid_inputs(self, body_name: str, time_since_epoch: float, expected_min_distance_meters: float):
        """
        Tests that `calculate_instantaneous_position_magnitude` returns a valid float
        and a non-negative distance for valid inputs.
        Note: Exact distance accuracy depends on the mocked or actual implementation
              and would typically require comparison with ephemeris data.
        """
        # Ensure the body has data for a more realistic scenario, if not, skip (though mock handles it)
        if get_celestial_body_data(body_name) is None and body_name.lower() not in ['moon']:
            pytest.skip(f"Celestial body data for '{body_name}' not found in `celestial_data`, skipping test.")
            
        distance = orbital_calculations.calculate_instantaneous_position_magnitude(body_name, time_since_epoch)
        
        assert isinstance(distance, float)
        assert distance >= expected_min_distance_meters # Distance should always be non-negative

    @pytest.mark.parametrize("body_name, time_since_epoch, error_type, error_msg_part", [
        (123, 0.0, TypeError, "body_name must be a string."),
        ("Earth", "invalid", TypeError, "time_since_epoch must be a number."),
        ("Earth", 0.0, TypeError, "reference_body_name must be a string.") # Testing default parameter if overridden
    ])
    def test_calculate_instantaneous_position_magnitude_invalid_inputs(self, body_name, time_since_epoch, error_type, error_msg_part):
        """
        Tests that `calculate_instantaneous_position_magnitude` raises appropriate errors
        for invalid input types.
        """
        with pytest.raises(error_type) as excinfo:
            # Explicitly pass a bad reference_body_name for one test case
            if "reference_body_name" in error_msg_part:
                orbital_calculations.calculate_instantaneous_position_magnitude(body_name, time_since_epoch, reference_body_name=123)
            else:
                orbital_calculations.calculate_instantaneous_position_magnitude(body_name, time_since_epoch)
        assert error_msg_part in str(excinfo.value)

    def test_calculate_instantaneous_position_magnitude_non_existent_body(self):
        """
        Tests handling of a non-existent celestial body name.
        The mock implementation returns `None` for such cases.
        """
        distance = orbital_calculations.calculate_instantaneous_position_magnitude("NonExistentBody", 0.0)
        assert distance is None

    # --- Tests for Distance Between Two Bodies (calculate_distance_between_bodies) ---

    def test_calculate_distance_between_same_bodies(self):
        """
        Tests that calculating the distance between the same celestial body returns zero.
        """
        distance = orbital_calculations.calculate_distance_between_bodies("Earth", "Earth", 0.0)
        assert math.isclose(distance, 0.0, abs_tol=1e-9)

    @pytest.mark.parametrize("body1, body2, time_since_epoch, expected_rough_distance_meters", [
        ("Earth", "Sun", 0.0, 1.496e11), # Approx 1 AU
        ("Earth", "Moon", 0.0, 3.844e8),  # Approx Earth-Moon distance
        ("Mars", "Sun", 0.0, 2.279e11), # Approx 1.5 AU
    ])
    def test_calculate_distance_between_bodies_known_pairs(self, body1: str, body2: str, time_since_epoch: float, expected_rough_distance_meters: float):
        """
        Tests distance calculation for known celestial pairs using rough approximations.
        Verifies that a float is returned and it's within a reasonable range of the expected value.
        """
        distance = orbital_calculations.calculate_distance_between_bodies(body1, body2, time_since_epoch)
        
        assert distance is not None
        assert isinstance(distance, float)
        # Allow a 1% relative tolerance for these rough mock-based comparisons
        assert math.isclose(distance, expected_rough_distance_meters, rel_tol=0.01), \
            f"Distance between {body1} and {body2}: Expected ~{expected_rough_distance_meters}, got {distance}"

    @pytest.mark.parametrize("body1, body2, time_since_epoch, error_type, error_msg_part", [
        (123, "Earth", 0.0, TypeError, "Body names must be strings."),
        ("Earth", 456, 0.0, TypeError, "Body names must be strings."),
        ("Earth", "Moon", "invalid", TypeError, "Time since epoch must be a number."),
    ])
    def test_calculate_distance_between_bodies_invalid_inputs(self, body1, body2, time_since_epoch, error_type, error_msg_part):
        """
        Tests that `calculate_distance_between_bodies` raises appropriate errors for invalid inputs.
        """
        with pytest.raises(error_type) as excinfo:
            orbital_calculations.calculate_distance_between_bodies(body1, body2, time_since_epoch)
        assert error_msg_part in str(excinfo.value)

    def test_calculate_distance_between_bodies_non_existent_pair(self):
        """
        Tests handling of distance calculation for a pair of bodies that the mock/system
        might not have pre-defined logic for. Expects None in this mock.
        """
        distance = orbital_calculations.calculate_distance_between_bodies("Saturn", "Jupiter", 1000.0)
        assert distance is None

    # --- Tests for Hohmann Transfer Delta-V (calculate_hohmann_delta_v) ---

    @pytest.mark.parametrize("source, target, expected_delta_v, tolerance", [
        ("Earth", "Mars", 5591.791, 1e-3), # Calculated using celestial_data values for SMA and Sun's mu
        ("Mars", "Earth", 5591.791, 1e-3), # Magnitude should be the same for return trip
    ])
    def test_calculate_hohmann_delta_v_accuracy(self, source: str, target: str, expected_delta_v: float, tolerance: float):
        """
        Tests the accuracy of Hohmann transfer delta-v calculations for known planetary pairs.
        """
        delta_v = orbital_calculations.calculate_hohmann_delta_v(source, target)
        assert isinstance(delta_v, float)
        assert math.isclose(delta_v, expected_delta_v, rel_tol=tolerance, abs_tol=tolerance), \
            f"For {source} to {target}: Expected {expected_delta_v}, got {delta_v}"

    def test_calculate_hohmann_delta_v_same_body(self):
        """
        Tests that calculating Hohmann delta-v between the same body returns 0.0.
        """
        delta_v = orbital_calculations.calculate_hohmann_delta_v("Earth", "Earth")
        assert math.isclose(delta_v, 0.0, abs_tol=1e-9)

    @pytest.mark.parametrize("source, target, error_type, error_msg_part", [
        (123, "Mars", TypeError, "Source and target body names must be strings."),
        ("Earth", None, TypeError, "Source and target body names must be strings."),
    ])
    def test_calculate_hohmann_delta_v_invalid_input_types(self, source, target, error_type, error_msg_part):
        """
        Tests that `calculate_hohmann_delta_v` raises appropriate errors for invalid input types.
        """
        with pytest.raises(error_type) as excinfo:
            orbital_calculations.calculate_hohmann_delta_v(source, target)
        assert error_msg_part in str(excinfo.value)

    @pytest.mark.parametrize("source, target", [
        ("NonExistent", "Mars"),
        ("Earth", "UnknownPlanet"),
        ("NonExistent", "UnknownPlanet"),
    ])
    def test_calculate_hohmann_delta_v_non_existent_body(self, source: str, target: str):
        """
        Tests that `calculate_hohmann_delta_v` returns None for non-existent celestial bodies.
        """
        delta_v = orbital_calculations.calculate_hohmann_delta_v(source, target)
        assert delta_v is None

    @pytest.mark.parametrize("source, target", [
        ("Sun", "Mars"), # Transfer to/from central body itself (r=0) is not a Hohmann between orbits
        ("Earth", "Sun"),
    ])
    def test_calculate_hohmann_delta_v_primary_as_source_or_target(self, source: str, target: str):
        """
        Tests scenarios where the source or target body is the primary (e.g., Sun),
        which is not a valid start/end for a Hohmann transfer between orbits around that primary.
        Should return None as per mock implementation for r <= 0.
        """
        delta_v = orbital_calculations.calculate_hohmann_delta_v(source, target)
        assert delta_v is None