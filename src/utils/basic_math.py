def add_numbers(a, b):
    """
    Sums two numbers and returns the result.

    Args:
        a (int or float): The first number.
        b (int or float): The second number.

    Returns:
        (int or float): The sum of the two numbers.

    Raises:
        TypeError: If either 'a' or 'b' is not a number (int or float).
    """
    if not isinstance(a, (int, float)) or not isinstance(b, (int, float)):
        raise TypeError("Inputs 'a' and 'b' must be numbers (int or float).")
    return a + b

def multiply_numbers(a, b):
    """
    Multiplies two numbers and returns the result.

    Args:
        a (int or float): The first number.
        b (int or float): The second number.

    Returns:
        (int or float): The product of the two numbers.

    Raises:
        TypeError: If either 'a' or 'b' is not a number (int or float).
    """
    if not isinstance(a, (int, float)) or not isinstance(b, (int, float)):
        raise TypeError("Inputs 'a' and 'b' must be numbers (int or float).")
    return a * b