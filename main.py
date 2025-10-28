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
    max = speed * secs
    return max

## calculate how long it would take for user to reach destination given speed
def calc_time(destination, speed, planets):
    # convert Mm/s to m/s
    speed = speed * 10**6 
    # get index of dictionary in list
    index = next((index for (index, d) in enumerate(planets) if d["name"] == destination), None)
    secs = (planets[index]['distance']) / speed
    years = secs / 365.25 /24 / 60 / 60
    return years

## calculate how long it would take for user to reach destination in respect to the Earth's reference frame
def calc_time_earth(years,speed):
    # convert Mm/s to m/s
    speed = speed * 10**6 
    secs = years * 365.25 * 24 * 60 * 60
    c = 299792458
    try:
        user_secs = secs / math.sqrt(1- ((speed**2)/ (c**2) ))
        user_years = user_secs / 365.25 / 24 / 60 / 60
    except:
        pass
    return user_years

## returns a list of dictionaries with information of planets that the user can reach within their lifetime
def gen_planets(max):
    df = pd.read_csv('exoplanets.csv')
    df = df[["pl_name",'sy_dist']]
    df = df.reset_index()
    planets = []
    for i in df.index:
        try:
            # convert distance to meters from pl
            distance = df.loc[i,'sy_dist'] * 3.085677581e16
            if distance <= max:
                planets.append({'name': df.loc[i,'pl_name'], 'distance': distance}) 
            else: 
                continue
        except:
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
    except:
        pass
    return new_date

## converts date string into datetime object
def convert_date(date_str):
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
        except:
            print("Speed must be a number greater than 0 and cannot be greater than the speed of light. Re-enter speed.")
    # Prompt user to enter date until input is in valid format
    while True:    
        date = input("Date of departure from Earth in ISO format (YYYY-MM-DD): ")
        if c.is_valid_date(date):
            date = convert_date(date)
            break
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
        except:
            print("Enter a valid age between 18 and 74.")
    max_distance = calc_max_dis(speed, age)    
    print("Choose a destination from the list below: ")
    # get list of viable planets and display their names and distances
    planets = gen_planets(max_distance) 
    for dic in planets:
        print(f"Planet \"{dic['name']}\" is {(dic['distance']):.5g} meters away.")
    while True:
        destination = input("Choose a destination from the list above: ")
        if any(d['name'] == destination for d in planets):
            break
        else: 
            print("Please choose a planet from the list above.")
    years= calc_time(destination, speed, planets)
    years_ref = calc_time_earth(years, speed)
    new_age = calc_age(years,age)
    try:
        arrival_date = calc_arrival(date,years_ref)
        print(f"You will reach your destination in {years:.3g} years on {arrival_date:%m/%d/%Y} when you are {new_age} years old while {years_ref:.5g} years would have passed on Earth.")
    except:
        print(f"You will reach your destination in {years:.3g} years when you are {new_age} years old while {years_ref:.5g} years would have passed on Earth.")

if __name__ == "__main__":
    main()
