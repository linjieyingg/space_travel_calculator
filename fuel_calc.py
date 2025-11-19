def calculate_fuel_cost(total_fuel_mass_needed, fuel_price_per_unit):
    """
    Calculates the total cost of fuel given the total mass of fuel needed and its price per unit mass.

    Inputs:
        total_fuel_mass_needed (numeric): The total mass of fuel required for the journey
                                         (e.g., in kilograms, tons, or any consistent mass unit).
        fuel_price_per_unit (numeric): The cost of one unit of fuel mass
                                       (e.g., cost per kilogram, cost per ton).

    Outputs:
        Returns a dictionary containing 'total_fuel_mass_needed' (the input fuel mass)
        and 'total_cost' (the calculated total cost) if inputs are valid.
        Returns None if `total_fuel_mass_needed` or `fuel_price_per_unit` is
        negative or not a valid numeric type.
    """
    # Validate inputs to ensure they are numeric and non-negative
    if not isinstance(total_fuel_mass_needed, (int, float)) or total_fuel_mass_needed < 0:
        return None
    if not isinstance(fuel_price_per_unit, (int, float)) or fuel_price_per_unit < 0:
        return None

    # Calculate total cost based on fuel mass and price per unit
    total_cost = total_fuel_mass_needed * fuel_price_per_unit

    return {
        'total_fuel_mass_needed': total_fuel_mass_needed,
        'total_cost': total_cost
    }