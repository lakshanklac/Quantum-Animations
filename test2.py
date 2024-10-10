# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 11:12:07 2024

@author: K.L.A.C. LAKSHAN
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy.special import hermite
import math

# Constants
hbar = 1.054571817e-34  # Reduced Planck's constant, J·s
me = 9.10938356e-31  # Mass of an electron, kg
w = 1.0  # Frequency (arbitrary for now)

# Spatial grid parameters
L = 1e-8  # Spatial limit in meters (10 nm)
N = 1000  # Number of spatial points
x = np.linspace(-L, L, N)  # Position array in meters

# Function to calculate the wavefunction for quantum harmonic oscillator
def wavefunction(n, x):
    epsilon = np.sqrt(me * w / hbar) * x  # Dimensionless variable ε
    Hn = hermite(n)(epsilon)  # Hermite polynomial H_n(ε)
    norm_factor = (me * w / (np.pi * hbar)) ** 0.25 / np.sqrt(2 ** n * math.factorial(n))  # Normalization factor
    psi_n = norm_factor * Hn * np.exp((-epsilon ** 2) / 2)  # Wavefunction Ψ_n(x)
    return psi_n

# Function to plot the wavefunction for a given quantum number n
def plot_wavefunction(n):
    psi_n = wavefunction(n, x)
    
    plt.plot(x * 1e9, psi_n, label=f"n={n}")
    plt.xlabel("Position x (nm)")
    plt.ylabel(r'$\psi_n(x)$')
    plt.title(f"Wavefunction for Quantum Harmonic Oscillator (n={n})")
    plt.grid(True)
    plt.legend()
    plt.show()

# Plot the wavefunction for different quantum numbers
n_values = [0, 1, 2, 3]  # You can change these values as needed

for n in n_values:
    plot_wavefunction(n)
