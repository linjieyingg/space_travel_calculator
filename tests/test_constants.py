# tests/test_constants.py
import pytest
from constants import G_GRAVITATIONAL, C_LIGHT_MPS

def test_gravitational_constant_value_and_type():
    """
    Test that G_GRAVITATIONAL from constants.py is a float
    and has the expected value (CODATA 2018).
    """
    assert isinstance(G_GRAVITATIONAL, float)
    assert G_GRAVITATIONAL == pytest.approx(6.67430e-11, rel=1e-10)

def test_speed_of_light_value_and_type():
    """
    Test that C_LIGHT_MPS from constants.py is a float
    and has the exact defined value.
    """
    assert isinstance(C_LIGHT_MPS, float)
    assert C_LIGHT_MPS == pytest.approx(299_792_458.0, rel=1e-9)

def test_no_deprecated_gravitational_constant_alias():
    """
    Verify that the old `GRAVITATIONAL_CONSTANT` alias is no longer present
    directly in the `constants` module, reflecting the consolidation change.
    """
    # Attempt to access the old alias, which should now raise an AttributeError
    # because it has been removed from constants.py
    with pytest.raises(AttributeError):
        # We need to explicitly check the module's attributes, not just import
        # as a direct import would fail before even reaching the test body.
        # This workaround ensures the module is loaded and we check its dict.
        import sys
        if 'constants' in sys.modules:
            del sys.modules['constants'] # Reload to ensure fresh state if needed
        import constants as test_constants
        _ = test_constants.GRAVITATIONAL_CONSTANT