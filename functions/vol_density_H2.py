# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 19:34:08 2024

@author: fism
"""
# Sources: 
#     https://demaco-cryogenics.com/blog/energy-density-of-hydrogen/
#     https://h2tools.org/bestpractices/liquid-storage-vessels
    
import matplotlib.pyplot as plt
import numpy as np

# Data
pressures = [1, 350, 700]  # Atmospheric, 350 barg, 700 barg (in bar)
density_values = [0.083, 23.715, 39.75]  # Corresponding volumetric mass density values for gaseous hydrogen (in kg/m³)
liquid_pressure = 8.5  # Pressure for liquid hydrogen (in bar) 
liquid_density = 70.9  # Volumetric mass density for liquid hydrogen (in kg/m³)

# Plotting
plt.plot(pressures, density_values, marker='o', linestyle='-', color='b', label='Gaseous Hydrogen @ 20°C')
plt.scatter(liquid_pressure, liquid_density, color='r', label='Liquid Hydrogen @ -252.9°C')
plt.title('Volumetric Mass Density of Hydrogen')
plt.xlabel('Pressure (bar)')
plt.ylabel('Volumetric Mass Density (kg/m³)')
plt.legend()
plt.grid(True)
plt.show()
