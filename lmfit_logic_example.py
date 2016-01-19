# Ed Colosky
# lmfit logic example

# third party modules
import numpy as np
import matplotlib.pyplot as plt
from lmfit import Parameters, minimize, fit_report

# Isotherm function to test
def Freundlich(x, Kf, n):
    return Kf * x ** (1/n)

# lmfit 
def residual(pars, x, data=None):
    vals = pars.valuesdict()
    Kf = vals['Kf']
    n = vals['n']
    model = Freundlich(x, Kf=Kf, n=n)
    if data is None:
        return model
    return (model - data)

# Set general constants
np.random.seed(0)
n = 101
xmin = 0.0
xmax = 100.0

# Set example isotherm constants 
p_true = Parameters()
p_true.add('Kf', value=.2)
p_true.add('n', value=5)

# Create isotherm; w/ noise
noise = np.random.normal(scale=0.04, size=n)
x = np.linspace(xmin, xmax, n)
data = residual(p_true, x) + noise

# Set initial lmfit params for freundlich function
fit_params = Parameters()
fit_params.add('Kf', value=1.0)
fit_params.add('n', value=1.0)

# fit model against data first with nelder
nelder = minimize(residual, fit_params, args=(x, data), method='nelder')
# fit model using leastsqr method but with nelder fitted parameters
leastsqr = minimize(residual, nelder.params, args=(x, data), method='leastsqr')

print('nelder fit')
print(fit_report(nelder))
print('leastsqr fit with nelder fit params')
print(fit_report(leastsqr))

# plot 
plt.scatter(x, data)
plt.plot(x, Freundlich(x,
                       Kf=leastsqr.params.valuesdict()['Kf'],
                       n=leastsqr.params.valuesdict()['n']))

