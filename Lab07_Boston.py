# -*- coding: utf-8 -*-
"""
Created on Wed Feb 27 11:27:18 2019

@author: btown
"""

"""
Some Notes for myself

omega=angular frequency
omega=sqrt(k/m)

eigenvalue=omega^2
or
omega=+/- sqrt(eigenvalue)
"""
import numpy as np
from numpy.linalg import eigh
import matplotlib.pyplot as plt


def calc_frequencies_and_modes(matrix, k_over_m):
    """
    This function calculates the frequency and eigenvectors of the spring setup
    
    Parameters
    ----------
    matrix: Array
        A 2x2 matrix corresponding to the the setup of springs we are simulating
    k_over_m: float
        The constant defined by the spring constant over the mass of the each of the masses
    
    Returns
    -------
    freq0: float
        the first frequency of the spring setup
    freq1: float
        the second frequency of the spring setup
    eigenvector: Array
        An array containing the 2 eigenvectors associated with this spring setup
    """
    eigenvalue,eigenvector=eigh(k_over_m*matrix)
    freq0=np.sqrt(np.abs(eigenvalue[0]))
    freq1=np.sqrt(np.abs(eigenvalue[1]))
    return freq0, freq1, eigenvector

N_times=1000
t_init=0
t_end=10
matrix1=np.array([[-2,1],[1,-2]])
time=np.linspace(t_init,t_end,num=N_times)
time=time.reshape(N_times,1)

"""
print(calc_frequencies_and_modes(matrix1,1)[0])
print(calc_frequencies_and_modes(matrix1,1)[1])
print(calc_frequencies_and_modes(matrix1,1)[2])
"""

def calc_components_from_initial_conditions(x_init,modes): #modes are eigenvectors
    """
    Function to calculate the components of the equation of motion
    
    Parameters
    ----------
    x_init: array of initial starting points
        This is the array of the starting points for the masses on the spring.
        Positive values mean the mass is moved to the right of the equilibrium position.
    modes: array 
        an array that contains the eigenvectors for the spring setup.
        
    Returns
    -------
    a: float
        one of the components to the equation of motion
    b: float
        one of the other components to the equation of motion
    """
    a=x_init@modes[0]
    b=x_init@modes[1]
    return a, b


def position_of_masses(x_init, time, matrix, k_over_m):
    """
    Function to calculate the position's of the masses on the springs
    over time. 
    
    Parameters
    ----------
    x_init: array of initial starting points
        This is the array of the starting points for the masses on the spring.
        Positive values mean the mass is moved to the right of the equilibrium position.
    time : array of times
        Times at which the positions have been calculated.
    matrix: Array
        A 2x2 matrix corresponding to the the setup of springs we are simulating
    k_over_m: float
        The constant defined by the spring constant over the mass of the each of the masses
        
    Returns:
    --------
    x: Array
        an array that contains the positions of the two masses over time
    """
    freq0,freq1,eigenvector=calc_frequencies_and_modes(matrix, k_over_m)
    a,b=calc_components_from_initial_conditions(x_init,eigenvector)
    x=a*np.cos(freq0*time)*eigenvector[0]+b*np.cos(freq1*time)*eigenvector[1]
    return x
    

def plot_motion_of_masses(x, time, title='Motion of Springs'):
    """
    Function to make a plot of motion of masses as a function of time. The time
    should be on the vertical axis and the position on the horizontal axis.
    Parameters
    ----------
    x : array of position, N_times by 2 elements
        The array of positions, set up so that x[:, 0] is the position of mass
        1 relative to equilibrium and x[:, 1] is the position of mass 2.
    time : array of times
        Times at which the positions have been calculated.
    title : str
        A descriptive title for the plot to make grading easier.
    """
    # Nothing special about these, but they look nice
    x1_equilibrium_pos = 3
    x2_equilibrium_pos = 6

    x1 = x[:, 0] + x1_equilibrium_pos
    x2 = x[:, 1] + x2_equilibrium_pos

    plt.plot(x1, time, label='Mass 1')
    plt.plot(x2, time, label='Mass 2')
    plt.xlim(0, 9)
    plt.xlabel('Position')
    plt.ylabel('Time')
    plt.legend()
    plt.title(title)
    plt.show()

position0=position_of_masses([1,-1], time, matrix1, 1)
position1=position_of_masses([1,1], time, matrix1, 1)
position2=position_of_masses([1,0], time, matrix1, 1)
position3=position_of_masses([0,1], time, matrix1, 1)
position4=position_of_masses([-.82,.3], time, matrix1, 1)

plot_motion_of_masses(position0,time)
plot_motion_of_masses(position1,time)
plot_motion_of_masses(position2,time)
"""
This plot doesn't oscillate as a pure sine wave, its harder to predict
"""
plot_motion_of_masses(position3,time)
"""
This also doesn't oscillate as a pure sine wave, though its not the same
as the previous case since the initial conditions aren't equivalent
"""
plot_motion_of_masses(position4,time)

"""
Answer for part 5:
The  matrix should remain the same, but the coefficient in front of the matrix should change
The springs should oscillate just like they did before, except with a higher frequency. 
"""
position5=position_of_masses([1,1], time, matrix1, 10)
position6=position_of_masses([1,-1], time, matrix1, 10)
position7=position_of_masses([1,0], time, matrix1, 10)
position8=position_of_masses([0,1], time, matrix1, 10)

plot_motion_of_masses(position5,time, 'Heavy Masses') 
plot_motion_of_masses(position6,time, 'Heavy Masses')
"""
The motion when they are set to their eigenvectors is either in phase or antiphase, 
just like what you would expect
"""
plot_motion_of_masses(position7,time, 'Heavy Masses')
plot_motion_of_masses(position8,time, 'Heavy Masses')
"""
These are different because they are not eigenvectors and therefore wont have a motion 
that ends up in phase or antiphase, its something that oscillates in a way thats not a pure
sine wave
"""




