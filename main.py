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

from datetime import datetime, timedelta
import pandas as pd
import checks as c
import math

## return max distance user can travel in their lifetime
def calc_max_dis(speed, age):
    years = 100 - age # given one lifetime is 100 years
    # convert Mm/s to m/s
    speed = speed * 10**6 
    # convert years to seconds
    secs = years * 365.25 * 24 * 60 * 60
    max_distance = speed * secs # Renamed max to max_distance to avoid shadowing built-in max()
    return max_distance

## calculate how long it would take for user to reach destination given speed
def calc_time(destination, speed, planets):
    # convert Mm/s to m/s
    speed = speed * 10**6 
    # get index of dictionary in list
    # The destination is guaranteed to exist in planets list due to validation in main()
    index = next((index for (index, d) in enumerate(planets) if d["name"] == destination), None)
    
    # Add a check here just in case, though main's validation should prevent None.
    if index is None:
        raise ValueError(f"Destination '{destination}' not found in the list of viable planets.")
        
    secs = (planets[index]['distance']) / speed
    years = secs / 365.25 /24 / 60 / 60
    return years

## calculate how long it would take for user to reach destination in respect to the Earth's reference frame
def calc_time_earth(years,speed):
    # convert Mm/s to m/s
    speed_ms = speed * 10**6 
    secs = years * 365.25 * 24 * 60 * 60
    c = 299792458
    
    user_years = years # Initialize with a fallback value

    try:
        # Calculate the Lorentz factor term (1 - v^2/c^2)
        sqrt_term_arg = 1 - ((speed_ms**2) / (c**2))

        # Handle cases where relativistic calculation is problematic
        if sqrt_term_arg < 0:
            # This implies speed_ms was effectively >= c due to float precision or input bypass.
            # Time dilation would involve sqrt of a negative number, which is undefined.
            # Fallback to non-relativistic time.
            user_years = years
        elif sqrt_term_arg == 0:
            # If speed_ms == c, the denominator becomes zero, leading to division by zero or infinite dilation.
            # Fallback to non-relativistic time for practical purposes.
            user_years = years 
        else:
            user_secs = secs / math.sqrt(sqrt_term_arg)
            user_years = user_secs / 365.25 / 24 / 60 / 60
    except Exception as e:
        # Catch any other unexpected errors during the relativistic calculation.
        # print(f"Warning in calc_time_earth: {e}. Falling back to non-relativistic time.") # For debugging
        user_years = years # Fallback if any unexpected error occurs
        
    return user_years

## returns a list of dictionaries with information of planets that the user can reach within their lifetime
def gen_planets(max_distance):
    df = pd.read_csv('exoplanets.csv')
    df = df[["pl_name",'sy_dist']]
    df = df.reset_index()
    planets = []
    for i in df.index:
        try:
            # convert stellar distance (sy_dist) in parsecs to meters
            distance = df.loc[i,'sy_dist'] * 3.085677581e16
            if distance <= max_distance:
                planets.append({'name': df.loc[i,'pl_name'], 'distance': distance}) 
        except KeyError:
            # Handle cases where 'pl_name' or 'sy_dist' might be missing for a row
            # print(f"Warning: Missing data for a row in exoplanets.csv at index {i}. Skipping.") # For debugging
            continue
        except Exception:
            # Catch any other unexpected error for individual row processing
            # print(f"Warning: An error occurred processing row {i} in exoplanets.csv. Skipping.") # For debugging
            continue
    return planets

##  Calculate age of user when they reach destination
def calc_age(years, age):
    new_age = math.floor(age + years)
    return new_age

## Calculate the date of when user reaches their destination
def calc_arrival(date, years):
    try:
        days = years * 365.25 
        td = timedelta(days=days)
        new_date = date + td
        return new_date
    except OverflowError:
        # This occurs if the number of days results in a date outside the datetime object's valid range.
        # print("Warning: Arrival date calculation resulted in a date out of range due to extremely long travel time.") # For debugging
        return None
    except Exception as e:
        # Catch any other unexpected errors during date calculation.
        # print(f"Warning: An unexpected error occurred calculating arrival date: {e}") # For debugging
        return None

## converts date string into datetime object
def convert_date(date_str):
    # This function is typically called after c.is_valid_date, which should pre-validate the format.
    # However, for defensive programming, a try-except here could catch unforeseen parsing issues.
    # For this refactor, we rely on main's try-except block when calling this function.
    date = datetime.strptime(date_str, '%Y-%m-%d')
    return date

def main():
    # Prompt user to enter speed until entered speed is valid
    while True:
        try:
            speed = float(input("Speed of your spacecraft in megameters per second (Mm/s): "))
            if c.is_valid_speed(speed):
                break
            else:
                print("Speed must be a number greater than 0 and cannot be greater than the speed of light. Re-enter speed.")
        except ValueError:
            print("Invalid input. Speed must be a numeric value. Re-enter speed.")
        except Exception:
            print("An unexpected error occurred with speed input. Please re-enter speed.")
            
    # Prompt user to enter date until input is in valid format
    while True:    
        date_str = input("Date of departure from Earth in ISO format (YYYY-MM-DD): ")
        if c.is_valid_date(date_str):
            try:
                date = convert_date(date_str)
                break
            except ValueError:
                # This catches errors from datetime.strptime if c.is_valid_date was not strict enough.
                print("Error: Date format is superficially correct but could not be parsed. Please re-enter a valid date.")
            except Exception:
                print("An unexpected error occurred while processing the date. Please re-enter a valid date.")
        else:
            print("Error with date format. Re-enter valid date.")
            
    # Prompt user to enter age until age is valid
    while True:
        try:
            age = int(input("Enter your age in years: "))
            if c.is_valid_age(age):
                break
            else:
                print("Enter a valid age between 18 and 74.")
        except ValueError:
            print("Invalid input. Age must be a whole number. Re-enter age.")
        except Exception:
            print("An unexpected error occurred with age input. Please re-enter age.")
            
    max_distance = calc_max_dis(speed, age)    
    
    print("\nChoose a destination from the list below: ")
    # Get list of viable planets and display their names and distances
    try:
        planets = gen_planets(max_distance) 
    except FileNotFoundError:
        print("Error: 'exoplanets.csv' not found. Please ensure the file is in the same directory.")
        return # Exit the program gracefully if the data file is missing
    except pd.errors.EmptyDataError:
        print("Error: 'exoplanets.csv' is empty. No planet data to process.")
        return # Exit if the CSV is empty
    except Exception as e:
        print(f"An unexpected error occurred while reading or processing 'exoplanets.csv': {e}")
        return # Exit for other critical pandas or file-related errors
        
    if not planets:
        print("No viable planets found within your lifetime given the speed and age.")
        print("Try increasing your speed or lowering your age to find destinations.")
        return # Exit gracefully if no planets meet the criteria
        
    for dic in planets:
        print(f"Planet \"{dic['name']}\" is {(dic['distance']):.5g} meters away.")
        
    while True:
        destination = input("Choose a destination from the list above: ")
        if any(d['name'] == destination for d in planets):
            break
        else: 
            print("Please choose a planet from the list above.")
            
    try:
        years = calc_time(destination, speed, planets)
    except ValueError as e:
        print(f"Error calculating travel time: {e}")
        return
    except Exception as e:
        print(f"An unexpected error occurred during travel time calculation: {e}")
        return

    years_ref = calc_time_earth(years, speed) # years_ref is guaranteed to be a float now
    new_age = calc_age(years,age)
    
    arrival_date = calc_arrival(date, years_ref) # arrival_date can be datetime object or None

    if arrival_date is not None:
        print(f"\nYou will reach your destination in {years:.3g} years on {arrival_date:%m/%d/%Y} when you are {new_age} years old while {years_ref:.5g} years would have passed on Earth.")
    else:
        # Fallback if arrival date could not be calculated (e.g., overflow or other error)
        print(f"\nYou will reach your destination in {years:.3g} years when you are {new_age} years old while {years_ref:.5g} years would have passed on Earth. (Arrival date could not be precisely determined due to extreme duration or date limits.)")

if __name__ == "__main__":
    main()