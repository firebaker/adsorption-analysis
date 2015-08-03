""" statstical functions """
from scipy.stats.distributions import t
import numpy as np


# calculate linear/non-linear regression SSR (sum of squared residuals)
# requires 'popt' variable from curve_fit()
def _regSSR(func, xdata, ydata, popt):
    residuals = ydata - func(xdata, *popt)
    SSR = sum(residuals**2)
    return SSR


# calculate AIC (Akaike information criterion)
# requires 'popt' variable from curve_fit()
def _regAIC(func, xdata, ydata, popt, k=2):
    n_obs = len(ydata)
    n_params = len(ydata)
    SSR = _regSSR(func, xdata, ydata, popt)
    MLE = SSR / len(ydata)  # biased maximum likelihood estimation
    loglik = (- n_obs / 2 * np.log(2 * np.pi)
              - n_obs / 2 * np.log(MLE) - 1 /
              (2 * MLE) * SSR)
    AIC = - 2 * loglik + k * (n_params + 1)
    return AIC


# calculate the confidence interval of the regression line
# requires 'popt' and 'pcov' variables from curve_fit()
def _reg_conf(func, xdata, ydata, alpha, popt, pcov):
    # my calculation of general statistics
    n_obs = len(ydata)  # number of data points
    n_params = len(popt)  # number of parameters
    dof = max(0, n_obs - n_params)  # degrees of freedom
    tval = t.ppf(1.0-alpha/2., dof)  # calculate student-t value
    # return popt for lower and upper conf intervals
    lower = []
    upper = []
    for param, sigma in zip(popt, np.sqrt(np.diag(pcov))):
        lower.append(param - sigma * tval)
        upper.append(param + sigma * tval)
    return lower, upper
