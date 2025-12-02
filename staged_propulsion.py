import math
from typing import List, Union
from . import propulsion_system


def calculate_multistage_fuel_mass(
    payload_mass: float,
    stage_dry_masses: List[float],
    stage_specific_impulses: List[float],
    stage_delta_vs: List[float]
) -> List[float]:
    """
    Calculates the fuel mass required for each stage of a multi-stage rocket.

    The calculation proceeds from the last stage to the first, as the "effective payload"
    of an earlier stage includes the dry mass and fuel mass of all subsequent stages,
    plus the actual payload.

    Args:
        payload_mass (float): The mass of the payload (e.g., spacecraft) in kilograms,
                              that needs to be delivered to the final orbit. Must be positive.
        stage_dry_masses (List[float]): A list of dry masses for each stage in kilograms.
                                       The list order corresponds to stages from first to last.
                                       Each mass must be positive.
        stage_specific_impulses (List[float]): A list of specific impulses (in seconds)
                                               for the engine of each stage.
                                               The list order corresponds to stages from first to last.
                                               Each impulse must be positive.
        stage_delta_vs (List[float]): A list of delta-V requirements (in m/s) for each stage.
                                      The list order corresponds to stages from first to last.
                                      Each delta-V must be non-negative.

    Returns:
        List[float]: A list of fuel masses (in kilograms) required for each stage,
                     ordered from first stage to last stage.

    Raises:
        ValueError: If inputs are invalid (e.g., non-positive masses/ISPs, negative delta-Vs,
                    mismatched list lengths, empty lists).
    """
    if not isinstance(payload_mass, (int, float)) or payload_mass <= 0:
        raise ValueError("Payload mass must be a positive number.")

    num_stages = len(stage_dry_masses)

    if num_stages == 0:
        raise ValueError("Cannot calculate fuel for zero stages. Provide at least one stage.")

    if not (len(stage_specific_impulses) == num_stages and len(stage_delta_vs) == num_stages):
        raise ValueError("The length of 'stage_dry_masses', 'stage_specific_impulses', "
                         "and 'stage_delta_vs' lists must be identical.")

    for i, dry_mass in enumerate(stage_dry_masses):
        if not isinstance(dry_mass, (int, float)) or dry_mass <= 0:
            raise ValueError(f"Stage dry mass at index {i} must be a positive number.")
    for i, isp in enumerate(stage_specific_impulses):
        if not isinstance(isp, (int, float)) or isp <= 0:
            raise ValueError(f"Specific impulse at index {i} must be a positive number.")
    for i, dv in enumerate(stage_delta_vs):
        if not isinstance(dv, (int, float)) or dv < 0:
            raise ValueError(f"Delta-V at index {i} must be a non-negative number.")

    calculated_fuel_masses = [0.0] * num_stages

    # This represents the total mass (payload + subsequent stages' dry mass + subsequent stages' fuel mass)
    # that the current stage needs to accelerate after its own fuel is expended and itself is jettisoned.
    mass_to_be_lifted_by_current_stage_dry_mass = payload_mass

    # Iterate backwards from the last stage to the first
    for i in range(num_stages - 1, -1, -1):
        current_stage_dry_mass = stage_dry_masses[i]
        current_stage_isp = stage_specific_impulses[i]
        current_stage_delta_v = stage_delta_vs[i]

        # The 'dry mass' input for the Tsiolkovsky equation for this stage calculation
        # is its own dry mass PLUS all the mass that it needs to propel.
        effective_dry_mass_for_calc = current_stage_dry_mass + mass_to_be_lifted_by_current_stage_dry_mass

        try:
            fuel_needed_for_this_stage = propulsion_system.calculate_required_fuel_mass(
                delta_v=current_stage_delta_v,
                dry_mass=effective_dry_mass_for_calc,
                specific_impulse=current_stage_isp
            )
        except ValueError as e:
            # Re-raise with more context if propulsion_system validation fails
            raise ValueError(f"Error calculating fuel for stage {i+1}: {e}")

        calculated_fuel_masses[i] = fuel_needed_for_this_stage

        # Update the total mass that the *previous* stage (i-1) would need to lift.
        # This new mass includes the dry mass of the current stage, its calculated fuel,
        # and the mass that this stage itself had to lift (which was the old
        # mass_to_be_lifted_by_current_stage_dry_mass).
        mass_to_be_lifted_by_current_stage_dry_mass += current_stage_dry_mass + fuel_needed_for_this_stage

    return calculated_fuel_masses


if __name__ == "__main__":
    print("--- Staged Propulsion Fuel Calculator Demo ---")

    # Mock propulsion_system.calculate_required_fuel_mass for standalone testing
    # In a real application, this module would import the actual propulsion_system.
    class MockPropulsionSystem:
        STANDARD_GRAVITY = 9.80665

        @staticmethod
        def calculate_required_fuel_mass(delta_v: float, dry_mass: float, specific_impulse: float) -> float:
            if not isinstance(delta_v, (int, float)) or delta_v < 0:
                raise ValueError("delta_v must be a non-negative number.")
            if not isinstance(dry_mass, (int, float)) or dry_mass <= 0:
                raise ValueError("dry_mass must be a positive number.")
            if not isinstance(specific_impulse, (int, float)) or specific_impulse <= 0:
                raise ValueError("specific_impulse must be a positive number.")

            exhaust_velocity = specific_impulse * MockPropulsionSystem.STANDARD_GRAVITY

            if delta_v == 0:
                return 0.0

            mass_ratio = math.exp(delta_v / exhaust_velocity)
            fuel_mass = dry_mass * (mass_ratio - 1)
            return fuel_mass

    # Replace the actual import for this self-test block
    propulsion_system = MockPropulsionSystem()

    # Test Case 1: Simple two-stage rocket
    print("\nTest Case 1: Two-stage rocket to orbit")
    payload = 1000.0  # kg
    stage_dry_masses_t1 = [20000.0, 5000.0]  # kg (Stage 1, Stage 2)
    stage_specific_impulses_t1 = [300.0, 450.0]  # s (Stage 1, Stage 2)
    stage_delta_vs_t1 = [4000.0, 3000.0]  # m/s (Stage 1, Stage 2)

    try:
        fuel_masses_t1 = calculate_multistage_fuel_mass(
            payload, stage_dry_masses_t1, stage_specific_impulses_t1, stage_delta_vs_t1
        )
        print(f"Payload Mass: {payload:.2f} kg")
        print(f"Stage Dry Masses: {stage_dry_masses_t1} kg")
        print(f"Stage Specific Impulses: {stage_specific_impulses_t1} s")
        print(f"Stage Delta-Vs: {stage_delta_vs_t1} m/s")
        for i, fuel in enumerate(fuel_masses_t1):
            print(f"  Fuel Mass for Stage {i+1}: {fuel:.2f} kg")
        print(f"  Total Fuel Mass: {sum(fuel_masses_t1):.2f} kg")
    except ValueError as e:
        print(f"Error in Test Case 1: {e}")

    # Test Case 2: Single-stage rocket
    print("\nTest Case 2: Single-stage rocket")
    payload_t2 = 500.0  # kg
    stage_dry_masses_t2 = [1000.0]  # kg
    stage_specific_impulses_t2 = [350.0]  # s
    stage_delta_vs_t2 = [7000.0]  # m/s

    try:
        fuel_masses_t2 = calculate_multistage_fuel_mass(
            payload_t2, stage_dry_masses_t2, stage_specific_impulses_t2, stage_delta_vs_t2
        )
        print(f"Payload Mass: {payload_t2:.2f} kg")
        print(f"Stage Dry Masses: {stage_dry_masses_t2} kg")
        print(f"Stage Specific Impulses: {stage_specific_impulses_t2} s")
        print(f"Stage Delta-Vs: {stage_delta_vs_t2} m/s")
        for i, fuel in enumerate(fuel_masses_t2):
            print(f"  Fuel Mass for Stage {i+1}: {fuel:.2f} kg")
        print(f"  Total Fuel Mass: {sum(fuel_masses_t2):.2f} kg")
    except ValueError as e:
        print(f"Error in Test Case 2: {e}")

    # Test Case 3: Error handling - Mismatched list lengths
    print("\nTest Case 3: Error - Mismatched list lengths")
    try:
        calculate_multistage_fuel_mass(100, [1000, 500], [300], [4000, 3000, 100])
    except ValueError as e:
        print(f"Caught expected error: {e}")

    # Test Case 4: Error handling - Zero payload mass
    print("\nTest Case 4: Error - Zero payload mass")
    try:
        calculate_multistage_fuel_mass(0, [1000], [300], [4000])
    except ValueError as e:
        print(f"Caught expected error: {e}")

    # Test Case 5: Error handling - Negative specific impulse
    print("\nTest Case 5: Error - Negative specific impulse")
    try:
        calculate_multistage_fuel_mass(100, [1000], [-300], [4000])
    except ValueError as e:
        print(f"Caught expected error: {e}")

    # Test Case 6: Error handling - Zero stages
    print("\nTest Case 6: Error - Zero stages")
    try:
        calculate_multistage_fuel_mass(100, [], [], [])
    except ValueError as e:
        print(f"Caught expected error: {e}")

    # Test Case 7: Three-stage rocket with some 0 delta-v (e.g. for coasting or small adjustments)
    print("\nTest Case 7: Three-stage rocket with 0 delta-v for a middle stage")
    payload_t3 = 50.0  # kg
    stage_dry_masses_t3 = [1000.0, 200.0, 50.0]  # kg (Stage 1, Stage 2, Stage 3)
    stage_specific_impulses_t3 = [400.0, 300.0, 250.0]  # s
    stage_delta_vs_t3 = [5000.0, 0.0, 1000.0]  # m/s

    try:
        fuel_masses_t3 = calculate_multistage_fuel_mass(
            payload_t3, stage_dry_masses_t3, stage_specific_impulses_t3, stage_delta_vs_t3
        )
        print(f"Payload Mass: {payload_t3:.2f} kg")
        print(f"Stage Dry Masses: {stage_dry_masses_t3} kg")
        print(f"Stage Specific Impulses: {stage_specific_impulses_t3} s")
        print(f"Stage Delta-Vs: {stage_delta_vs_t3} m/s")
        for i, fuel in enumerate(fuel_masses_t3):
            print(f"  Fuel Mass for Stage {i+1}: {fuel:.2f} kg")
        print(f"  Total Fuel Mass: {sum(fuel_masses_t3):.2f} kg")
    except ValueError as e:
        print(f"Error in Test Case 7: {e}")