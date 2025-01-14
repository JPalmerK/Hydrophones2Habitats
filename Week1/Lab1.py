# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 09:29:36 2024

@author: Kait Palmer

Introduction to Python
"""

# Load libraryies or parts of libraries
import numpy as np
import matplotlib.pyplot as plt




##################################################################
# Harmonic Oscilator

# Use python to simulate a mass on a spring from the equations in class
##################################################################

def location(s,m,t, A =1, B=1):
    """
    Dispacement (x) of a mass on a string given Amplitude A, B and stiffness 
    and angular frequency ()

    Parameters
    ----------
    s : TYPE
        DESCRIPTION.
    m : TYPE
        DESCRIPTION.
    t : TYPE
        DESCRIPTION.
    A : TYPE, optional
        DESCRIPTION. The default is 1.
    B : TYPE, optional
        DESCRIPTION. The default is 1.

    Returns
    -------
    None.

    """
    w = np.sqrt(s/m)
    x = A* np.cos(w*t)+ B* np.sin(w*t)
    
    return x
    
def velocity(s,m,t, A =1, B=1):
    """
    Velocity (x) of a mass on a string given Amplitude A, B and stiffness 
    and angular frequency ()


    Parameters
    ----------
    s : TYPE
        DESCRIPTION.
    m : TYPE
        DESCRIPTION.
    t : TYPE
        DESCRIPTION.
    A : TYPE, optional
        DESCRIPTION. The default is 1.
    B : TYPE, optional
        DESCRIPTION. The default is 1.

    Returns
    -------
    vel : TYPE
        DESCRIPTION.

    """
    w = np.sqrt(s/m)
    vel = -A*w*np.sin(w*t)+B*w*np.cos(w*t)
    
    return vel
    
    
def acceleration(s,m,t, A =1, B=1):
    """
    Velocity (x) of a mass on a string given Amplitude A, B and stiffness 
    and angular frequency ()


    Parameters
    ----------
    s : TYPE
        DESCRIPTION.
    m : TYPE
        DESCRIPTION.
    t : TYPE
        DESCRIPTION.
    A : TYPE, optional
        DESCRIPTION. The default is 1.
    B : TYPE, optional
        DESCRIPTION. The default is 1.

    Returns
    -------
    accel : TYPE
        DESCRIPTION.

    """
    w = np.sqrt(s/m)
    accel = -A*w**2*np.cos(w*t)+B*w**2*np.sin(w*t)
    
    return accel


# Declair variables (same as R, matlab)
# Test it out. Pick some random starting variables
m =200 # Mass (kg)
s =50 # spring stiffness
t =0 # time 

# We've declared defalutl varibales in the functions so we don't need to worry
# about setting A and B

accel =  acceleration(s,m,t, A =1, B=1)
displacement = location(s,m,t)


# Does the acceleartion equal the stiffness over the mass times the displacment?
accel == -(s/m)*displacement



# Usematlplot lib to plot the velocity displacement and acceleation as a 
# function of time


t = np.linspace(0, 13, num =200)
accel =  acceleration(s,m,t, A = 0)
vel = velocity(s,m,t, A = 0)
displacement = location(s,m,t, A = 0)

plt.figure()
plt.plot(t, accel, label='Acceleration')
plt.plot(t, vel, label='Velocity')
plt.plot(t, displacement, label='Displacement')
plt.legend()
plt.show()

#############################################################################
# Complex numbers
#############################################################################


##############################################################################
#
#############################################################################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Constants
c = 340  # Speed of sound (m/s)
frequency = 2  # Frequency of the wave (Hz)
wavelength = c / frequency
k = 2 * np.pi / wavelength  # Wave number
omega = 2 * np.pi * frequency  # Angular frequency
amplitude = 1  # Amplitude of the wave

# Space and time grid
x = np.linspace(0, 10, 500)  # 10 meters, 500 points
t = np.linspace(0, 1, 100)  # 1 second, 100 time points

# Initialize figure
fig, ax = plt.subplots()
line, = ax.plot(x, np.zeros_like(x), lw=2)
ax.set_ylim(-1.5, 1.5)
ax.set_xlim(0, 10)
ax.set_xlabel("Position (m)")
ax.set_ylabel("Pressure")
ax.set_title("1D Wave Propagation")

# Animation function
def update(frame):
    y = amplitude * np.sin(k * x - omega * frame)
    line.set_ydata(y)
    return line,

ani = FuncAnimation(fig, update, frames=t, interval=50, blit=True)
plt.show()



#############################################################################


###########



# Constants
frequency = 1000  # Frequency in Hz
c = 343  # Speed of sound in m/s
wavelength = c / frequency
k = 2 * np.pi / wavelength  # Wave number
rho0 = 1.21  # Air density in kg/m³
A = 1 # Amplitude, unitless for now but generally Pascals


###################################################################
# Boundary conditions
###################################################################

# Grid setup to plot the instatnaneous pressure on a square grid
x = np.linspace(-2, 2, 500)  # x-axis (meters)
y = np.linspace(-2, 2, 500)  # y-axis (meters)
X, Y = np.meshgrid(x, y)
r = np.sqrt(X**2 + Y**2)

# Source position
source_x, source_y = 0, 0.5  # Source at (0, 0.5 m)

# Compute distance between source and all grid positions
source_r = np.sqrt((X - source_x)**2 + (Y - source_y)**2)

# Radiation pattern for monopole source wiht phase and amplitude
pressure_field = np.exp(1j * k * source_r) / (4 * np.pi * source_r)

# Add boundary effects
boundary_y = 0  # Boundary at y=0

# Pressure-release boundary (anti-symmetric)
image_r = np.sqrt((X - source_x)**2 + (Y + source_y)**2)
pressure_release_field = pressure_field - np.exp(1j * k * image_r) / (4 * np.pi * image_r)

# Rigid boundary (symmetric)
rigid_field = pressure_field + np.exp(1j * k * image_r) / (4 * np.pi * image_r)

# Plotting function
def plot_field(field, title):
    plt.figure(figsize=(8, 6))
    plt.pcolormesh(x, y, np.real(field), shading='auto', cmap='RdBu')
    plt.colorbar(label="Pressure Amplitude (Real Part)")
    plt.title(title)
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    plt.axhline(boundary_y, color="black", linestyle="--", label="Boundary")
    plt.legend()
    plt.tight_layout()
    plt.show()

# Visualizations
plot_field(pressure_field, "Free Field Radiation")
plot_field(pressure_release_field, "Pressure-Release Boundary")
plot_field(rigid_field, "Rigid Boundary")
