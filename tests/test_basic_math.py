```python
# tests/test_basic_math.py

import pytest
from src.utils.basic_math import add_numbers

def test_add_two_positive_integers():
    """Test addition of two positive integers."""
    assert add_numbers(2, 3) == 5

def test_add_positive_and_negative_integer():
    """Test addition of a positive and a negative integer."""
    assert add_numbers(5, -2) == 3

def test_add_two_negative_integers():
    """Test addition of two negative integers."""
    assert add_numbers(-7, -3) == -10

def test_add_zero_to_number():
    """Test addition involving zero."""
    assert add_numbers(10, 0) == 10
    assert add_numbers(0, -5) == -5
    assert add_numbers(0, 0) == 0

def test_add_two_floats():
    """Test addition of two floating-point numbers."""
    assert add_numbers(2.5, 3.5) == 6.0
    # Use pytest.approx for floating-point comparisons
    assert pytest.approx(add_numbers(0.1, 0.2)) == 0.3

def test_add_integer_and_float():
    """Test addition of an integer and a float."""
    assert add_numbers(2, 3.5) == 5.5
    assert add_numbers(2.5, 3) == 5.5

def test_add_large_numbers():
    """Test addition with large integer values."""
    assert add_numbers(1_000_000, 2_000_000) == 3_000_000
    assert add_numbers(999999999999, 1) == 1000000000000

# Assuming add_numbers might raise an error for non-numeric types
# If add_numbers relies directly on Python's `+`, these tests might change to check for Python's TypeError.
def test_add_non_numeric_inputs_raises_type_error():
    """Test that non-numeric inputs raise a TypeError."""
    with pytest.raises(TypeError):
        add_numbers("hello", 5)
    with pytest.raises(TypeError):
        add_numbers(5, [1, 2])
    with pytest.raises(TypeError):
        add_numbers("a", "b")
```