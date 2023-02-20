# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt

# Define the density function and its parameters
def density_func(x, mu, sigma):
    return np.exp(-(x - mu)**2 / (2 * sigma**2)) / (np.sqrt(2 * np.pi) * sigma)

mu = 0.5
sigma = 0.1

# Generate a set of x values for the plot
x = np.linspace(0, 1, 1000)

# Compute the density function at each x value
y = density_func(x, mu, sigma)

# Compute the cumulative density function
cdf = np.cumsum(y) / np.sum(y)

# Compute the derivative of the cumulative density function
dydx = np.gradient(cdf) / np.gradient(x)

# Plot the density function and its cumulative distribution
fig, ax = plt.subplots(nrows=2, sharex=True)

ax[0].plot(x, y)
ax[0].set_ylabel('Density')

ax[1].plot(x, cdf)
ax[1].set_ylabel('Cumulative Distribution')

# Add a shaded area to the plot to illustrate the density flux conservation
ax[0].fill_between(x, 0, y, where=(dydx >= 0), alpha=0.2, color='green')
ax[0].fill_between(x, 0, y, where=(dydx < 0), alpha=0.2, color='red')
ax[1].fill_between(x, 0, 1, where=(dydx >= 0), alpha=0.2, color='green')
ax[1].fill_between(x, 0, 1, where=(dydx < 0), alpha=0.2, color='red')

ax[1].set_xlabel('x')

plt.show()