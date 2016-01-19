# Ed Colosky
# example for testing adsorption-analysis

# system modules
import sys
import py_compile

# change working directory
sys.path.append('/home/firebaker/git/adsorption-analysis')
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

# Create isotherm; w/wo noise
noise = np.random.normal(scale=0.04, size=n)
x = np.linspace(xmin, xmax, n)
y = isotherm.Freundlich.IsothermFunc(x, Kf=.2, n=5) + noise

#x = [0.611, 0.714, 0.895, 0.636, 0.618, 0.713, 0.697, 1.102, 0.939, 10.453,\
#10.564, 20.699, 29.438, 23.872, 28.801]
#y = [72.8644863643, 66.4529783924, 64.5926417041, 183.9867974246,\
#179.4697105899, 190.6825726419, 406.0358170324, 372.6665359867,\
#373.8061894016, 464.2506629708, 456.8038525485, 188.7308634721,\
#351.3011293429, 507.9787582089, 340.1717770768]

output = adsorption.AdsorptionAnalysis(data=[x, y])

plt.scatter(x, y)

x_iso = np.linspace(0, xmax, n)
y_Lin = isotherm.Linear.IsothermFunc(x_iso, **output.Linear.minimizedFit.params.valuesdict())
plt.plot(x_iso, y_Lin, 'b-')
y_Fre = isotherm.Freundlich.IsothermFunc(x_iso, **output.Freundlich.nelder.params.valuesdict())
plt.plot(x_iso, y_Fre, 'r-')
y_Lan = isotherm.Langmuir.IsothermFunc(x_iso, **output.Langmuir.nelder.params.valuesdict())
plt.plot(x_iso, y_Lan, 'g-')