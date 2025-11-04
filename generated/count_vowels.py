def count_vowels(text):
    """
    Counts the number of vowels (a, e, i, o, u, case-insensitive) in a given string.

    Args:
        text (str): The input string.

    Returns:
        int: The total count of vowels in the string.
    """
    vowels = "aeiou"
    vowel_count = 0
    for char in text:
        if char.lower() in vowels:
            vowel_count += 1
    return vowel_count