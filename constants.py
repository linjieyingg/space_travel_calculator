```python
"""
constants.py

This module defines universal physical constants and other widely used
scientific constants relevant to space travel calculations within the
space_travel_calculator application.

Constants defined here ensure a single source of truth for critical values,
improving consistency and maintainability across the project.
"""

# Universal Physical Constants
C_LIGHT_MPS = 299_792_458.0  # Speed of light in a vacuum (m/s)
G_GRAVITATIONAL = 6.67430e-11 # Gravitational Constant (m^3 kg^-1 s^-2)

# For backward compatibility / alias
# This provides an alternative name for the gravitational constant,
# ensuring modules that might expect this specific variable name continue to function.
GRAVITATIONAL_CONSTANT = G_GRAVITATIONAL

# Commented-out Examples for future expansion:
# Other potentially useful constants can be added here as needed.
# MU_SUN = 1.32712440018e20  # Sun's gravitational parameter (m^3/s^2)
# AU_TO_METERS = 149597870700.0 # Astronomical Unit to meters
# EARTH_RADIUS_METERS = 6371000.0 # Average radius of Earth in meters
```