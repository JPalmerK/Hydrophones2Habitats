# -*- coding: utf-8 -*-
"""
Created on Fri Jan  3 15:01:01 2025

@author: kaity
"""

import numpy as np


############################################################################
# Variables in Python
############################################################################

a = np.log10(100)

m = 10

# Check to the Variable Explorer to make sure that the variables you created
# esist 


############################################################################
# Functions in Python, called 'definitions'
############################################################################

def print_hello_world():
   print("Hello, World!")

# Call the function
print_hello_world()

def calculateTL(r, m=15):
    """
    Explain what the function does here. 

    Parameters
    ----------
    r : float
        range from source in km.
    m : TYPE, optional
        TL coefficient unitless. The default is 15.

    Returns
    -------
    None.

    """
    TL = m*np.log10(r)
    return(TL)

calculateTL(1000, m=20)


####################################################################
# Powers in Python
####################################################################

6**2 # six squared
np.sqrt(36)
36**(.5)

# using numpy
5**3 # 5 to the 3rd = 125 
np.power(5,3)
np.power(125, (1/3)) # cube root of 125 should be 2

##################################################################
# For loops in python
##################################################################


# Calculate the squares of numbers from 1 to 10
for number in range(1, 11):
    square = number ** 2  # Compute the square of the number
    print(f"The square of {number} is {square}")
        
    
# List of cetacens
cetaceans = ["strap toothed dolphin", "vaquita", "pilot whale", "puffin pig"]

# For loop to iterate over the list
for whale in cetaceans:
    print(f"I love {whale}!")

