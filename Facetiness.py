# -*- coding: utf-8 -*-
"""
@author: Thor Parmentier
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Data
x = np.array([13.993575, 20.997532, 28.0100945, 35.0341342, 42.0725260000002, 22.2590748, 5])
y = np.array([0.976000000000001, 1.6, 1.928, 2.188, 2.564, 1.647, 0])

# Logarithmic function that passes through (5, 0)
def log_fit(x, a):
    return a * (np.log(x) - np.log(5))

# Exclude the point (5, 0)
x_fit = x[x != 5]
y_fit = y[x != 5]

# Fit the data
a, _ = curve_fit(log_fit, x_fit, y_fit)

# Result
print(f"a = {a[0]}")

# Plotting
plt.scatter(x[:6], y[:6], label="Fukuzawa and Akitaya (1993)", color = "blue", zorder = 3)
plt.scatter(x[6], y[6], label="LaChapelle and Armstrong (1977)", color = "green", zorder = 4)
x_line = np.linspace(min(x), max(x) + 10, 500)
y_line = log_fit(x_line, a[0])
plt.plot(x_line, y_line, label = f"y = {a[0]:.4f} * (ln(x) - ln(5))", color = "red", zorder = 0)
# plt.title("Facetedness Constrained to (5, 0)")
plt.xlabel("Vapor pressure gradient [Pa/cm]")
plt.ylabel("'Facetedness' growth rate [nm/s]")
plt.xlim(0, 45)
plt.ylim(-0.05, 3)
plt.grid(True, linewidth=0.5, alpha=0.7, zorder = 0)
plt.legend()
plt.show()
