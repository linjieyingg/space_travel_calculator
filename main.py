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

from src.user_management.auth import register_user, login_user
from src.trip_planning.space_trip_planner import interactive_trip_planner

def display_auth_menu():
    """
    Displays the authentication menu (register/login) and handles user choice.
    Returns True if a user successfully logs in, False otherwise.
    """
    while True:
        print("\n--- Authentication Menu ---")
        print("1. Register")
        print("2. Login")
        print("3. Exit")
        choice = input("Enter your choice: ").strip()

        if choice == '1':
            if register_user():
                print("Registration successful! Please log in.")
            else:
                print("Registration failed. Please try again.")
        elif choice == '2':
            if login_user():
                print("Login successful!")
                return True
            else:
                print("Login failed. Please try again.")
        elif choice == '3':
            print("Exiting application. Goodbye!")
            return False
        else:
            print("Invalid choice. Please enter 1, 2, or 3.")

def main():
    """
    Primary entry point for the Space Travel Calculator application.
    Handles user authentication before proceeding to the main application logic.
    """
    print("Welcome to the Space Travel Calculator!")

    if display_auth_menu():
        print("\nLogin successful! Welcome to the Space Travel Calculator.")
        while True:
            print("\n--- Space Travel Calculator Menu ---")
            print("1. Plan a Space Trip")
            # Add other main application features here as they are developed
            print("2. Exit Application")
            
            main_app_choice = input("Please choose an option: ").strip()

            if main_app_choice == '1':
                interactive_trip_planner()
            elif main_app_choice == '2':
                print("Exiting application. Goodbye!")
                break # Exit the while loop, which will then exit the main function
            else:
                print("Invalid option. Please try again.")
    else:
        print("Authentication not completed. Application will exit.")


if __name__ == "__main__":
    main()