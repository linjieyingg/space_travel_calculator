def count_vowels(text):
    """
    Counts the number of vowels in a given string.

    Vowels are considered to be 'a', 'e', 'i', 'o', 'u', case-insensitive.

    Args:
        text (str): The input string.

    Returns:
        int: The total count of vowels in the string.
    """
    vowels = "aeiouAEIOU"
    count = 0
    for char in text:
        if char in vowels:
            count += 1
    return count
