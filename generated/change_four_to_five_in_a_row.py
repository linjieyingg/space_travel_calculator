def update_game_reference(text: str) -> str:
    """
    Replaces all occurrences of the phrase "four in a row" with "five in a row"
    within the provided text. This function is typically used to update game titles
    or descriptions.
    """
    return text.replace("four in a row", "five in a row")