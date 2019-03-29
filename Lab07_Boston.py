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
    eigenvalue,eigenvector=eigh(k_over_m*matrix)
    freq0=np.sqrt(np.abs(eigenvalue[0]))
    freq1=np.sqrt(np.abs(eigenvalue[1]))
    return freq0, freq1, eigenvector


matrix1=np.array([[-2,1],[1,-2]])
time=np.arange(0,10,.01)

"""
print(calc_frequencies_and_modes(matrix1,1)[0])
print(calc_frequencies_and_modes(matrix1,1)[1])
print(calc_frequencies_and_modes(matrix1,1)[2])
"""

def calc_components_from_initial_conditions(x_init,modes): #modes are eigenvectors
    a=x_init@modes[0]
    b=x_init@modes[1]
    return a, b


def position_of_masses(x_init, time, matrix, k_over_m):
    freq0,freq1,eigenvector=calc_frequencies_and_modes(matrix, k_over_m)
    a,b=calc_components_from_initial_conditions(x_init,eigenvector)
    x=a*np.cos(freq0*time)*eigenvector[0]+b*np.cos(freq1*time)*eigenvector[1]
    return x
    
position=position_of_masses([1,1], time, matrix1, 1)


def plot_motion_of_masses(x, time, title='bad title'):
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
    plt.legend()
    plt.title(title)

print(calc_components_from_initial_conditions(np.array([.1,0]),np.array([.70710678,.70710678])))







