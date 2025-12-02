import pytest
import math
from propulsion_system import calculate_required_fuel_mass, STANDARD_GRAVITY, calculate_exhaust_velocity

# Define a reasonable tolerance for floating-point comparisons
# Absolute tolerance for values close to zero, relative tolerance for larger values.
ABS_TOL = 1e-9
REL_TOL = 1e-6

# --- Tests for calculate_exhaust_velocity (new function in propulsion_system.py) ---

def test_calculate_exhaust_velocity_valid_input():
    """
    Test `calculate_exhaust_velocity` with a valid positive specific impulse.
    Verifies the output is correct based on STANDARD_GRAVITY.
    """
    specific_impulse = 300.0  # seconds
    expected_exhaust_velocity = specific_impulse * STANDARD_GRAVITY
    actual_exhaust_velocity = calculate_exhaust_velocity(specific_impulse)
    assert math.isclose(actual_exhaust_velocity, expected_exhaust_velocity, rel_tol=REL_TOL, abs_tol=ABS_TOL)

def test_calculate_exhaust_velocity_float_input():
    """
    Test `calculate_exhaust_velocity` with a float specific impulse.
    Ensures floating-point values are handled correctly.
    """
    specific_impulse = 450.5  # seconds
    expected_exhaust_velocity = specific_impulse * STANDARD_GRAVITY
    actual_exhaust_velocity = calculate_exhaust_velocity(specific_impulse)
    assert math.isclose(actual_exhaust_velocity, expected_exhaust_velocity, rel_tol=REL_TOL, abs_tol=ABS_TOL)

def test_calculate_exhaust_velocity_zero_specific_impulse_raises_error():
    """
    Test `calculate_exhaust_velocity` with zero specific impulse.
    Should raise a ValueError as specific impulse must be positive.
    """
    with pytest.raises(ValueError, match="Specific impulse must be a positive value."):
        calculate_exhaust_velocity(0.0)

def test_calculate_exhaust_velocity_negative_specific_impulse_raises_error():
    """
    Test `calculate_exhaust_velocity` with a negative specific impulse.
    Should raise a ValueError as specific impulse must be positive.
    """
    with pytest.raises(ValueError, match="Specific impulse must be a positive value."):
        calculate_exhaust_velocity(-100.0)

def test_calculate_exhaust_velocity_invalid_type_raises_error():
    """
    Test `calculate_exhaust_velocity` with non-numeric specific impulse.
    Should raise a ValueError for invalid input types.
    """
    with pytest.raises(ValueError, match="Specific impulse must be a numeric value."):
        calculate_exhaust_velocity("not_a_number")
    with pytest.raises(ValueError, match="Specific impulse must be a numeric value."):
        calculate_exhaust_velocity(None)
    with pytest.raises(ValueError, match="Specific impulse must be a numeric value."):
        calculate_exhaust_velocity([300])

# --- Tests for calculate_required_fuel_mass ---

@pytest.mark.parametrize(
    "delta_v, dry_mass, specific_impulse, expected_fuel_mass",
    [
        (1000.0, 1000.0, 300.0, 381.82855),  # Common LEO burn example
        (5000.0, 5000.0, 450.0, 4897.43321), # Interplanetary example
        (0.0, 1000.0, 300.0, 0.0),           # No delta_v, no fuel needed
        (2000.0, 500.0, 250.0, 1373.19323),
    ]
)
def test_calculate_required_fuel_mass_with_specific_impulse(
    delta_v, dry_mass, specific_impulse, expected_fuel_mass
):
    """
    Test `calculate_required_fuel_mass` with valid inputs, providing `specific_impulse`.
    """
    actual_fuel_mass = calculate_required_fuel_mass(
        delta_v=delta_v, dry_mass=dry_mass, specific_impulse=specific_impulse
    )
    assert math.isclose(actual_fuel_mass, expected_fuel_mass, rel_tol=REL_TOL, abs_tol=ABS_TOL)

@pytest.mark.parametrize(
    "delta_v, dry_mass, specific_impulse, expected_fuel_mass",
    [
        (1000.0, 1000.0, 300.0, 381.82855),
        (5000.0, 5000.0, 450.0, 4897.43321),
        (0.0, 1000.0, 300.0, 0.0),
        (2000.0, 500.0, 250.0, 1373.19323),
    ]
)
def test_calculate_required_fuel_mass_with_exhaust_velocity(
    delta_v, dry_mass, specific_impulse, expected_fuel_mass
):
    """
    Test `calculate_required_fuel_mass` with valid inputs, providing `exhaust_velocity`.
    The `specific_impulse` is converted to `exhaust_velocity` for testing.
    """
    exhaust_velocity = specific_impulse * STANDARD_GRAVITY
    actual_fuel_mass = calculate_required_fuel_mass(
        delta_v=delta_v, dry_mass=dry_mass, exhaust_velocity=exhaust_velocity
    )
    assert math.isclose(actual_fuel_mass, expected_fuel_mass, rel_tol=REL_TOL, abs_tol=ABS_TOL)

@pytest.mark.parametrize("delta_v", [-10.0, -0.001, "invalid", None])
def test_calculate_required_fuel_mass_invalid_delta_v_raises_error(delta_v):
    """
    Test `calculate_required_fuel_mass` with invalid (negative or non-numeric) `delta_v`.
    Should raise a ValueError.
    """
    with pytest.raises(ValueError, match="Delta-V must be a non-negative numeric value."):
        calculate_required_fuel_mass(delta_v=delta_v, dry_mass=1000.0, specific_impulse=300.0)

@pytest.mark.parametrize("dry_mass", [-10.0, 0.0, "invalid", None])
def test_calculate_required_fuel_mass_invalid_dry_mass_raises_error(dry_mass):
    """
    Test `calculate_required_fuel_mass` with invalid (non-positive or non-numeric) `dry_mass`.
    Should raise a ValueError.
    """
    with pytest.raises(ValueError, match="Dry mass must be a positive numeric value."):
        calculate_required_fuel_mass(delta_v=1000.0, dry_mass=dry_mass, specific_impulse=300.0)

def test_calculate_required_fuel_mass_missing_engine_characteristics_raises_error():
    """
    Test `calculate_required_fuel_mass` when neither `specific_impulse` nor `exhaust_velocity` is provided.
    Should raise a ValueError.
    """
    with pytest.raises(ValueError, match="Either specific_impulse or exhaust_velocity must be provided."):
        calculate_required_fuel_mass(delta_v=1000.0, dry_mass=1000.0)

@pytest.mark.parametrize("invalid_isp", [-10.0, 0.0, "invalid", None])
def test_calculate_required_fuel_mass_invalid_specific_impulse_raises_error(invalid_isp):
    """
    Test `calculate_required_fuel_mass` with invalid `specific_impulse` (non-positive or non-numeric).
    `None` case is handled by the missing engine characteristics check.
    """
    if invalid_isp is None:
        # This case is caught by the "missing engine characteristics" validation earlier.
        with pytest.raises(ValueError, match="Either specific_impulse or exhaust_velocity must be provided."):
            calculate_required_fuel_mass(delta_v=1000.0, dry_mass=1000.0, specific_impulse=invalid_isp)
    else:
        with pytest.raises(ValueError, match="Specific impulse must be a positive numeric value."):
            calculate_required_fuel_mass(delta_v=1000.0, dry_mass=1000.0, specific_impulse=invalid_isp)

@pytest.mark.parametrize("invalid_exhaust_vel", [-10.0, 0.0, "invalid", None])
def test_calculate_required_fuel_mass_invalid_exhaust_velocity_raises_error(invalid_exhaust_vel):
    """
    Test `calculate_required_fuel_mass` with invalid `exhaust_velocity` (non-positive or non-numeric).
    `None` case is handled by the missing engine characteristics check.
    """
    if invalid_exhaust_vel is None:
        # This case is caught by the "missing engine characteristics" validation earlier.
        with pytest.raises(ValueError, match="Either specific_impulse or exhaust_velocity must be provided."):
            calculate_required_fuel_mass(delta_v=1000.0, dry_mass=1000.0, exhaust_velocity=invalid_exhaust_vel)
    else:
        with pytest.raises(ValueError, match="Exhaust velocity must be a positive numeric value."):
            calculate_required_fuel_mass(delta_v=1000.0, dry_mass=1000.0, exhaust_velocity=invalid_exhaust_vel)

def test_calculate_required_fuel_mass_both_engine_characteristics_consistent():
    """
    Test `calculate_required_fuel_mass` when both `specific_impulse` and `exhaust_velocity` are provided
    and are consistent with each other. The function should proceed without error.
    """
    specific_impulse = 300.0
    exhaust_velocity = specific_impulse * STANDARD_GRAVITY
    expected_fuel_mass = 381.82855
    actual_fuel_mass = calculate_required_fuel_mass(
        delta_v=1000.0, dry_mass=1000.0,
        specific_impulse=specific_impulse, exhaust_velocity=exhaust_velocity
    )
    assert math.isclose(actual_fuel_mass, expected_fuel_mass, rel_tol=REL_TOL, abs_tol=ABS_TOL)

def test_calculate_required_fuel_mass_both_engine_characteristics_inconsistent_raises_error():
    """
    Test `calculate_required_fuel_mass` when both `specific_impulse` and `exhaust_velocity` are provided
    but are inconsistent with each other. Should raise a ValueError.
    """
    specific_impulse = 300.0
    # Make exhaust_velocity inconsistent with specific_impulse
    exhaust_velocity = specific_impulse * STANDARD_GRAVITY * 1.1
    with pytest.raises(ValueError, match="Provided specific_impulse and exhaust_velocity are inconsistent."):
        calculate_required_fuel_mass(
            delta_v=1000.0, dry_mass=1000.0,
            specific_impulse=specific_impulse, exhaust_velocity=exhaust_velocity
        )

def test_calculate_required_fuel_mass_extreme_delta_v_overflow_error():
    """
    Test `calculate_required_fuel_mass` with an extremely high delta_v
    that could lead to an `OverflowError` during `math.exp` calculation.
    The function should catch this and re-raise as a `ValueError`.
    """
    # A delta_v value that will make delta_v / eff_exhaust_velocity very large
    # For Isp=300, eff_exhaust_velocity is approx 2941.995 m/s.
    # If delta_v is 1_000_000, then exp(1_000_000 / 2941.995) = exp(340.57) which will overflow.
    extreme_delta_v = 1_000_000.0 # 1000 km/s, far beyond practical for a single stage
    dry_mass = 1000.0
    specific_impulse = 300.0
    with pytest.raises(ValueError, match="Calculated fuel mass is too large"):
        calculate_required_fuel_mass(
            delta_v=extreme_delta_v, dry_mass=dry_mass, specific_impulse=specific_impulse
        )