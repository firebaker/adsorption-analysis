# Ed Colosky
# example for testing adsorption-analysis

# system modules
import sys
import py_compile

# change working directory
sys.path.append('/home/firebaker/git/adsorption-analysis')
py_compile.compile('/home/firebaker/git/adsorption-analysis/adsorption.py')
py_compile.compile('/home/firebaker/git/adsorption-analysis/isotherm.py')

# third party modules
import numpy as np
import matplotlib.pyplot as plt

# import adsorption analysis modules
import adsorption
import isotherm

# Set general constants
np.random.seed(0)
n = 101
xmin = 0.0
xmax = 100.0

# Create Freundlich isotherm; w noise
noise = np.random.normal(scale=0.04, size=n)
x = np.linspace(xmin, xmax, n)
y = isotherm.Freundlich.IsothermFunc(x, Kf=.2, n=5) + noise

# plot isotherm
plt.scatter(x, y)

# example fitting using individual isotherms
example_Linear = isotherm.Linear([x, y])
example_Freundlich = isotherm.Freundlich([x, y])
example_Langmuir = isotherm.Langmuir([x, y])

# example fitting using Adsorption Analysis module
example_AdsorptionAnalysis = adsorption.AdsorptionAnalysis([x, y])