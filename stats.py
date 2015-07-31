""" statstical functions """
from scipy.stats.distributions import t
import numpy as np


# calculate linear/non-linear regression SSR (sum of squared residuals)
# requires 'popt' variable from curve_fit()
def _regSSR(func, xdata, ydata, popt):
    residuals = ydata - func(xdata, *popt)
    SSR = sum(residuals**2)
    return SSR


# calculate the confidence interval of the regression line
# requires 'popt' and 'pcov' variables from curve_fit()
def _confSSR(func, xdata, ydata, alpha, popt, pcov):
    # my calculation of general statistics
    n = len(xdata)  # number of data points
    p = len(popt)  # number of parameters
    dof = max(0, n - p)  # degrees of freedom
    tval = t.ppf(1.0-alpha/2., dof)  # calculate student-t value
    # return popt for lower and upper conf intervals
    lower = []
    upper = []
    for p, var in zip(popt, np.diag(pcov)):
        sigma = var**0.5
        lower.append(p-sigma * tval)
        upper.append(p+sigma * tval)
    return lower, upper

""" scratch code from kitchin group

http://kitchingroup.cheme.cmu.edu/blog/2013/02/12/Nonlinear-curve-fitting-with-parameter-confidence-intervals/

# Nonlinear curve fit with confidence interval
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats.distributions import  t

x = np.array([0.5, 0.387, 0.24, 0.136, 0.04, 0.011])
y = np.array([1.255, 1.25, 1.189, 1.124, 0.783, 0.402])

# this is the function we want to fit to our data
def func(x, a, b):
    'nonlinear function in a and b to fit to data'
    return a * x / (b + x)

initial_guess = [1.2, 0.03]
pars, pcov = curve_fit(func, x, y, p0=initial_guess)

alpha = 0.05 # 95% confidence interval = 100*(1-alpha)

n = len(y)    # number of data points
p = len(pars) # number of parameters

dof = max(0, n - p) # number of degrees of freedom

# student-t value for the dof and confidence level
tval = t.ppf(1.0-alpha/2., dof)

lower = []
upper = []
for i, p,var in zip(range(n), pars, np.diag(pcov)):
    sigma = var**0.5
    print 'p{0}: {1} [{2}  {3}]'.format(i, p,
                                  p - sigma*tval,
                                  p + sigma*tval)
    lower.append(p-sigma*tval)
    upper.append(p+sigma*tval)

import matplotlib.pyplot as plt
plt.plot(x,y,'bo ')
xfit = np.linspace(0,1)
yfit = func(xfit, pars[0], pars[1])
plt.plot(xfit,yfit,'b-')
yfit_lower = (func(xfit, *lower))
plt.plot(xfit, yfit_lower, '--', color='r')
yfit_upper = (func(xfit, *upper))
plt.plot(xfit, yfit_upper, '--', color='r')
plt.legend(['data','fit'],loc='best')
plt.show()
"""
