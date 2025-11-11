def calculate_fuel_cost(distance, fuel_efficiency, fuel_price_per_unit):
    if fuel_efficiency <= 0:
        return None

    total_fuel_needed = distance / fuel_efficiency
    total_cost = total_fuel_needed * fuel_price_per_unit

    return {
        'total_fuel_needed': total_fuel_needed,
        'total_cost': total_cost
    }
