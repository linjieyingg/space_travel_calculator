```python
# constants.py

"""
This module defines universal physical constants and other widely used
scientific constants relevant to space travel calculations within the
space_travel_calculator application.

Constants defined here ensure a single source of truth for critical values,
improving consistency and maintainability across the project.
"""

# Universal physical constants
C_LIGHT_MPS = 299792458.0  # Speed of light in a vacuum in m/s
G_GRAVITATIONAL = 6.67430e-11  # Gravitational Constant in m^3 kg^-1 s^-2

# Alias for backward compatibility
# This provides an alternative name for the gravitational constant,
# ensuring modules that might expect this specific variable name continue to function.
GRAVITATIONAL_CONSTANT = G_GRAVITATIONAL

# Commented-out examples for future expansion
# Other potentially useful constants can be added here as needed.
# MU_SUN = 1.32712440018e20  # Sun's standard gravitational parameter in m^3/s^2
# AU_TO_METERS = 149597870700  # Astronomical Unit in meters
# EARTH_RADIUS_METERS = 6371000 # Average Earth radius in meters
```