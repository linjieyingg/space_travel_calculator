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
        # User has successfully logged in, proceed to main application logic
        print("\n--- Main Application Placeholder ---")
        print("Application logic for space travel calculations will be orchestrated from other modules.")
        print("This message confirms successful authentication and entry into the main application.")
        # Here, you would typically call a function to start the main calculator interface
        # For example: start_space_calculator_interface()
    else:
        print("Authentication not completed. Application will exit.")


if __name__ == "__main__":
    main()