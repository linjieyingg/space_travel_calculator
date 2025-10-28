"""
Course Number: ENGR 13300
Semester: Fall 2024

Description:
    Space travel calculator that calculates time it takes to 
    reach a destination chosen by the user. 

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

# Checks if string is a valid date in ISO format
def is_valid_date(date_str, format_str='%Y-%m-%d'):
    try:
        datetime.strptime(date_str, format_str)
        return True
    except ValueError:
        return False
    
# Checks if speed is valid speed (equal to or less than speed of light)
def is_valid_speed(speed):
    # convert speed to m/s
    speed = speed * 10**6
    c = 299792458
    if speed <= c:
        if speed > 0:
            return True
    else:
        return False
    
# Checks if age is valid age    
def is_valid_age(age):
    # not allowed to go to space is user is a minor or age is < 0
    if age < 18:
        return False
    # not allowed to go to space is user is over 75
    elif age >= 75:
        return False
    else: 
        return True