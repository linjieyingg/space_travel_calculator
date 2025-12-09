"""
Course Number: ENGR 13300
Semester: Fall 2024

Description:
    This file serves as the main entry point for the space travel calculator application.
    (Note: The core logic for trajectory planning, calculations, and user interaction
    has been externalized to other modules as part of refactoring. This file
    now primarily orchestrates or provides a placeholder for the application start).

Assignment Information:
    Assignment:     Individual Project
    Team ID:        LC5 - 21
    Author:         Jieying Lin, lin1914@purdue.edu
    Date:           12/01/2024

Contributors:
    Name, login@purdue [repeat for each]

    My contributor(s) helped me:
    [ ] understand the assignment expectations without
        telling me how they will approach it.
    [ ] understand different ways to think about a solution
        without helping me plan my solution.
    [ ] think through the meaning of a specific error or
        bug present in my code without looking at my code.
    Note that if you helped somebody else with their code, you
    have to list that person as a contributor here as well.

Academic Integrity Statement:
    I have not used source code obtained from any unauthorized
    source, either modified or unmodified; nor have I provided
    another student access to my code.  The project I am
    submitting is my own original work.
"""

# Import the new angle conversion functions
from src.utils.angle_conversions import degrees_to_radians, radians_to_degrees

def main():
    """
    Placeholder for the main application entry point.
    All core space travel calculation logic has been moved to other modules.
    This function now also demonstrates the usage of angle conversion utilities.
    """
    print("Welcome to the Space Travel Calculator!")
    print("Application logic will be orchestrated from other modules.")
    print("For now, this is a placeholder demonstrating utility functions.")

    # --- Demonstrate angle conversion utilities ---
    print("\n--- Demonstrating Angle Conversion Utilities ---")

    # Example 1: Degrees to Radians
    angle_deg = 90.0
    try:
        angle_rad = degrees_to_radians(angle_deg)
        print(f"Converting {angle_deg} degrees to radians: {angle_rad:.4f} radians")
    except TypeError as e:
        print(f"Error converting {angle_deg} degrees: {e}")

    angle_deg_neg = -45.0
    try:
        angle_rad_neg = degrees_to_radians(angle_deg_neg)
        print(f"Converting {angle_deg_neg} degrees to radians: {angle_rad_neg:.4f} radians")
    except TypeError as e:
        print(f"Error converting {angle_deg_neg} degrees: {e}")

    # Example 2: Radians to Degrees
    angle_rad_pi_half = 1.57079632679
    try:
        angle_deg_converted = radians_to_degrees(angle_rad_pi_half)
        print(f"Converting {angle_rad_pi_half:.4f} radians to degrees: {angle_deg_converted:.4f} degrees")
    except TypeError as e:
        print(f"Error converting {angle_rad_pi_half:.4f} radians: {e}")

    angle_rad_neg_pi = -3.14159265359
    try:
        angle_deg_neg_converted = radians_to_degrees(angle_rad_neg_pi)
        print(f"Converting {angle_rad_neg_pi:.4f} radians to degrees: {angle_deg_neg_converted:.4f} degrees")
    except TypeError as e:
        print(f"Error converting {angle_rad_neg_pi:.4f} radians: {e}")

    # Example 3: Invalid input type for degrees_to_radians
    invalid_input = "forty-five"
    try:
        degrees_to_radians(invalid_input)
    except TypeError as e:
        print(f"Attempted to convert invalid input '{invalid_input}' (string) to radians, caught expected error: {e}")

    print("\n--- End of Utility Demonstration ---")


if __name__ == "__main__":
    main()