import math
from typing import Union

# Standard gravity acceleration (g0) used in specific impulse calculation (m/s^2)
STANDARD_GRAVITY_ACCELERATION = 9.80665

def calculate_effective_exhaust_velocity(specific_impulse: float) -> float:
    """
    Calculates the effective exhaust velocity of a rocket engine.

    The effective exhaust velocity (ve) is derived from the specific impulse (Isp)
    and the standard gravity acceleration (g0).

    Args:
        specific_impulse (float): The specific impulse of the engine in seconds (s).

    Returns:
        float: The effective exhaust velocity in meters per second (m/s).

    Raises:
        TypeError: If `specific_impulse` is not a number.
        ValueError: If `specific_impulse` is non-positive.
    """
    if not isinstance(specific_impulse, (int, float)):
        raise TypeError("Specific impulse must be a number.")
    if specific_impulse <= 0:
        raise ValueError("Specific impulse must be a positive value.")

    # ve = Isp * g0
    return specific_impulse * STANDARD_GRAVITY_ACCELERATION

def calculate_mass_ratio(delta_v: float, effective_exhaust_velocity: float) -> float:
    """
    Calculates the mass ratio (initial mass / final mass) required for a given delta-v
    using the Tsiolkovsky rocket equation.

    The Tsiolkovsky rocket equation is: delta_v = ve * ln(m0 / mf),
    where m0 is initial mass, mf is final mass, and ve is effective exhaust velocity.
    Rearranging for mass ratio (m0/mf): m0/mf = exp(delta_v / ve).

    Args:
        delta_v (float): The change in velocity required in meters per second (m/s).
        effective_exhaust_velocity (float): The effective exhaust velocity
                                           of the engine in meters per second (m/s).

    Returns:
        float: The mass ratio (m0 / mf), which is unitless.

    Raises:
        TypeError: If `delta_v` or `effective_exhaust_velocity` are not numbers.
        ValueError: If `effective_exhaust_velocity` is non-positive.
                    If `delta_v` is negative (though 0 is allowed, meaning mass_ratio = 1).
    """
    if not isinstance(delta_v, (int, float)):
        raise TypeError("Delta-v must be a number.")
    if not isinstance(effective_exhaust_velocity, (int, float)):
        raise TypeError("Effective exhaust velocity must be a number.")

    if effective_exhaust_velocity <= 0:
        raise ValueError("Effective exhaust velocity must be a positive value.")
    if delta_v < 0:
        raise ValueError("Delta-v cannot be negative.")

    return math.exp(delta_v / effective_exhaust_velocity)

def calculate_fuel_mass(mass_ratio: float, dry_mass: float) -> float:
    """
    Calculates the required fuel mass for a mission, given the mass ratio and dry mass.

    The relationship is: Initial Mass (m0) = Final Mass (mf) + Fuel Mass (m_fuel).
    Since Mass Ratio (MR) = m0 / mf, then m0 = MR * mf.
    Substituting mf as dry_mass: m0 = MR * dry_mass.
    So, (MR * dry_mass) = dry_mass + m_fuel.
    Therefore, m_fuel = dry_mass * (MR - 1).

    Args:
        mass_ratio (float): The mass ratio (initial mass / final mass) required (unitless).
                            This is typically the output from `calculate_mass_ratio`.
        dry_mass (float): The final mass of the spacecraft (mf) without fuel,
                          often referred to as 'dry mass', in kilograms (kg).

    Returns:
        float: The required fuel mass in kilograms (kg).

    Raises:
        TypeError: If `mass_ratio` or `dry_mass` are not numbers.
        ValueError: If `mass_ratio` is less than 1.0 (implying negative or zero fuel).
                    If `dry_mass` is non-positive.
    """
    if not isinstance(mass_ratio, (int, float)):
        raise TypeError("Mass ratio must be a number.")
    if not isinstance(dry_mass, (int, float)):
        raise TypeError("Dry mass must be a number.")

    if mass_ratio < 1.0:
        raise ValueError("Mass ratio must be 1.0 or greater (1.0 means no fuel needed).")
    if dry_mass <= 0:
        raise ValueError("Dry mass must be a positive value.")

    return dry_mass * (mass_ratio - 1)

if __name__ == '__main__':
    # --- Test Cases ---
    print("--- Rocket Equation Calculations ---")

    # Test 1: Calculate effective exhaust velocity
    print("\n--- Effective Exhaust Velocity ---")
    try:
        isp_test = 450.0  # seconds
        ve = calculate_effective_exhaust_velocity(isp_test)
        print(f"Specific Impulse: {isp_test} s")
        print(f"Effective Exhaust Velocity: {ve:.2f} m/s")
        assert math.isclose(ve, 450.0 * STANDARD_GRAVITY_ACCELERATION), "VE calculation failed"
        print("  Test Passed: VE for positive Isp")
    except Exception as e:
        print(f"  Test Failed: VE for positive Isp - {e}")

    try:
        calculate_effective_exhaust_velocity(-100)
    except ValueError as e:
        print(f"  Test Passed: VE for negative Isp (expected error: {e})")
    except Exception as e:
        print(f"  Test Failed: VE for negative Isp - {e}")

    try:
        calculate_effective_exhaust_velocity("abc")
    except TypeError as e:
        print(f"  Test Passed: VE for non-numeric Isp (expected error: {e})")
    except Exception as e:
        print(f"  Test Failed: VE for non-numeric Isp - {e}")

    # Test 2: Calculate mass ratio
    print("\n--- Mass Ratio ---")
    delta_v_test = 5000.0  # m/s
    ve_test = calculate_effective_exhaust_velocity(450.0) # Use previously calculated VE
    try:
        mass_ratio = calculate_mass_ratio(delta_v_test, ve_test)
        print(f"Delta-V: {delta_v_test} m/s")
        print(f"Effective Exhaust Velocity: {ve_test:.2f} m/s")
        print(f"Mass Ratio (m0/mf): {mass_ratio:.4f}")
        assert math.isclose(mass_ratio, math.exp(delta_v_test / ve_test)), "Mass Ratio calculation failed"
        print("  Test Passed: Mass ratio for valid inputs")
    except Exception as e:
        print(f"  Test Failed: Mass ratio for valid inputs - {e}")

    try:
        mass_ratio_zero_dv = calculate_mass_ratio(0.0, ve_test)
        print(f"Mass Ratio for 0 Delta-V: {mass_ratio_zero_dv:.4f}")
        assert math.isclose(mass_ratio_zero_dv, 1.0), "Mass ratio for 0 delta-v should be 1.0"
        print("  Test Passed: Mass ratio for 0 delta-v")
    except Exception as e:
        print(f"  Test Failed: Mass ratio for 0 delta-v - {e}")

    try:
        calculate_mass_ratio(delta_v_test, -100.0)
    except ValueError as e:
        print(f"  Test Passed: Mass ratio for non-positive VE (expected error: {e})")
    except Exception as e:
        print(f"  Test Failed: Mass ratio for non-positive VE - {e}")
    
    try:
        calculate_mass_ratio(-10.0, ve_test)
    except ValueError as e:
        print(f"  Test Passed: Mass ratio for negative Delta-V (expected error: {e})")
    except Exception as e:
        print(f"  Test Failed: Mass ratio for negative Delta-V - {e}")

    # Test 3: Calculate fuel mass
    print("\n--- Fuel Mass ---")
    dry_mass_test = 1000.0  # kg
    try:
        fuel_mass = calculate_fuel_mass(mass_ratio, dry_mass_test)
        print(f"Mass Ratio: {mass_ratio:.4f}")
        print(f"Dry Mass: {dry_mass_test} kg")
        print(f"Required Fuel Mass: {fuel_mass:.2f} kg")
        assert math.isclose(fuel_mass, dry_mass_test * (mass_ratio - 1)), "Fuel Mass calculation failed"
        print("  Test Passed: Fuel mass for valid inputs")
    except Exception as e:
        print(f"  Test Failed: Fuel mass for valid inputs - {e}")

    try:
        fuel_mass_no_fuel = calculate_fuel_mass(1.0, dry_mass_test)
        print(f"Fuel Mass for Mass Ratio 1.0: {fuel_mass_no_fuel:.2f} kg")
        assert math.isclose(fuel_mass_no_fuel, 0.0), "Fuel mass for MR 1.0 should be 0.0"
        print("  Test Passed: Fuel mass for mass ratio 1.0")
    except Exception as e:
        print(f"  Test Failed: Fuel mass for mass ratio 1.0 - {e}")

    try:
        calculate_fuel_mass(0.5, dry_mass_test)
    except ValueError as e:
        print(f"  Test Passed: Fuel mass for mass ratio < 1.0 (expected error: {e})")
    except Exception as e:
        print(f"  Test Failed: Fuel mass for mass ratio < 1.0 - {e}")

    try:
        calculate_fuel_mass(mass_ratio, -500.0)
    except ValueError as e:
        print(f"  Test Passed: Fuel mass for non-positive dry mass (expected error: {e})")
    except Exception as e:
        print(f"  Test Failed: Fuel mass for non-positive dry mass - {e}")

    print("\n--- End of Tests ---")