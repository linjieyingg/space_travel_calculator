"""
trajectory_planner.py

Purpose:
This module provides functionality to calculate specific interplanetary trajectories,
such as Hohmann transfer orbits, between two celestial bodies orbiting a central star
(assumed to be the Sun). It determines the required delta-V for maneuvers and the
transfer time, leveraging functions from `orbital_mechanics.py` and data from
`celestial_data.py`.
"""

import math
import orbital_mechanics
import celestial_data

def calculate_hohmann_transfer(start_body_name: str, target_body_name: str) -> dict:
    """
    Calculates the parameters for a Hohmann transfer orbit between two celestial bodies
    orbiting a central body (assumed to be the Sun).

    Args:
        start_body_name (str): The name of the starting celestial body (e.g., "Earth").
        target_body_name (str): The name of the target celestial body (e.g., "Mars").

    Returns:
        dict: A dictionary containing the transfer parameters.
              Keys include:
              - 'delta_v1_ms' (float): Delta-V for the first burn at the starting orbit (m/s).
              - 'delta_v2_ms' (float): Delta-V for the second burn at the target orbit (m/s).
              - 'total_delta_v_ms' (float): Total Delta-V required (m/s).
              - 'transfer_time_seconds' (float): Time taken for the transfer in seconds.
              - 'transfer_time_days' (float): Time taken for the transfer in days.
              - 'error' (str, optional): An error message if the calculation fails.
              - 'message' (str, optional): A descriptive message for certain conditions (e.g., same planet).

    Raises:
        This function handles internal ValueErrors and other Exceptions and returns them
        in the 'error' key of the dictionary, rather than raising them directly, for easier
        integration into a broader application flow.
    """
    try:
        # Retrieve celestial data
        if start_body_name not in celestial_data.CELESTIAL_BODIES:
            raise ValueError(f"Starting body '{start_body_name}' not found in celestial data.")
        if target_body_name not in celestial_data.CELESTIAL_BODIES:
            raise ValueError(f"Target body '{target_body_name}' not found in celestial data.")

        r_start_original = celestial_data.CELESTIAL_BODIES[start_body_name]["orbital_radius_m"]
        r_target_original = celestial_data.CELESTIAL_BODIES[target_body_name]["orbital_radius_m"]
        mu_sun = celestial_data.MU_SUN

        if r_start_original <= 0 or r_target_original <= 0:
            raise ValueError("Orbital radii must be positive for both start and target bodies.")
        
        if r_start_original == r_target_original:
            return {
                "delta_v1_ms": 0.0,
                "delta_v2_ms": 0.0,
                "total_delta_v_ms": 0.0,
                "transfer_time_seconds": 0.0,
                "transfer_time_days": 0.0,
                "message": "Start and target bodies are at the same orbital radius; no transfer needed."
            }

        # Determine inner and outer radii for Hohmann transfer ellipse semi-major axis calculation
        r_inner = min(r_start_original, r_target_original)
        r_outer = max(r_start_original, r_target_original)

        # 1. Calculate semi-major axis of the Hohmann transfer orbit
        a_transfer = orbital_mechanics.calculate_transfer_orbit_semimajor_axis(r_inner, r_outer)

        # 2. Calculate velocities in the circular orbits
        v_circular_start = orbital_mechanics.calculate_orbital_velocity(mu_sun, r_start_original)
        v_circular_target = orbital_mechanics.calculate_orbital_velocity(mu_sun, r_target_original)

        # 3. Calculate velocities in the transfer orbit at the start and target radii
        v_transfer_at_start = orbital_mechanics.calculate_elliptical_velocity(mu_sun, r_start_original, a_transfer)
        v_transfer_at_target = orbital_mechanics.calculate_elliptical_velocity(mu_sun, r_target_original, a_transfer)

        # 4. Calculate Delta-Vs required for the two burns
        # Delta-V1: Burn to enter transfer orbit from starting circular orbit
        delta_v1 = v_transfer_at_start - v_circular_start

        # Delta-V2: Burn to enter target circular orbit from transfer orbit
        delta_v2 = v_circular_target - v_transfer_at_target

        # Total Delta-V is the sum of the magnitudes of the two burns
        total_delta_v = abs(delta_v1) + abs(delta_v2)

        # 5. Calculate transfer time
        # Period of the transfer ellipse: T = 2 * pi * sqrt(a_transfer^3 / mu)
        transfer_orbit_period_seconds = 2 * math.pi * math.sqrt(a_transfer**3 / mu_sun)
        transfer_time_seconds = transfer_orbit_period_seconds / 2
        transfer_time_days = transfer_time_seconds / (24 * 3600)

        return {
            "delta_v1_ms": delta_v1,
            "delta_v2_ms": delta_v2,
            "total_delta_v_ms": total_delta_v,
            "transfer_time_seconds": transfer_time_seconds,
            "transfer_time_days": transfer_time_days,
        }

    except ValueError as ve:
        return {"error": str(ve)}
    except Exception as e:
        # Catch any other unexpected errors during calculation
        return {"error": f"An unexpected error occurred during trajectory calculation: {e}"}