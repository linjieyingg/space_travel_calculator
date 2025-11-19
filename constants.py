"""
constants.py

This module defines universal physical constants and other widely used
scientific constants relevant to space travel calculations within the
space_travel_calculator application.

Constants defined here ensure a single source of truth for critical values,
improving consistency and maintainability across the project.
"""

# Universal Physical Constants

# Speed of light in a vacuum (meters per second)
# This value is based on the exact definition of the meter.
# It is used for relativistic calculations and speed validation.
C_LIGHT_MPS = 299_792_458.0  # Speed of light in a vacuum in m/s

# Gravitational Constant (m^3 kg^-1 s^-2)
# This is the CODATA 2018 recommended value.
# While not directly used in all current calculations, it's a fundamental
# constant for space mechanics and included for completeness and future expansion
# related to orbital mechanics or gravitational effects.
G_GRAVITATIONAL = 6.67430e-11  # Gravitational constant in m^3 kg^-1 s^-2

# Alias for gravitational constant to ensure compatibility with modules expecting 'GRAVITATIONAL_CONSTANT'
# Some modules in the repository, such as orbital_mechanics.py, may expect
# to import 'GRAVITATIONAL_CONSTANT' directly. This alias provides backward
# compatibility or consistency with other parts of the codebase.
GRAVITATIONAL_CONSTANT = G_GRAVITATIONAL

# Other potentially useful constants can be added here as needed, such as:
# # Standard Gravitational Parameter of the Sun (m^3 s^-2)
# MU_SUN = 1.32712440018e20
#
# # Astronomical Unit (meters)
# AU_TO_METERS = 149_597_870_700.0
#
# # Earth's average radius (meters)
# EARTH_RADIUS_METERS = 6_371_000.0
#
# # Earth's standard gravitational parameter (m^3 s^-2)
# MU_EARTH = 3.986004418e14