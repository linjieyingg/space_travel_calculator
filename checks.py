"""
Course Number: ENGR 13300
Semester: Fall 2024

Description:
    This module provides utility functions for validating date inputs,
    primarily for a space travel calculator application.

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

from datetime import datetime

def is_valid_date(date_str, format_str='%Y-%m-%d'):
    """
    Checks if a given string is a valid date according to a specified format.

    Inputs:
        date_str (str): The string representing the date to validate.
        format_str (str): The expected format of the date string (default: '%Y-%m-%d').

    Outputs:
        bool: True if date_str can be successfully parsed into a datetime object
              using format_str, False otherwise.
    """
    try:
        datetime.strptime(date_str, format_str)
        return True
    except ValueError:
        return False