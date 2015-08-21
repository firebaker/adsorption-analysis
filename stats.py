""" statstical functions """

# third party modules
import numpy as np
from scipy.stats.distributions import t


# linearRegresion
def _linReg(x, m=1, b=1):
    return m * x + b


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


# confidence interval calculation of the regression line
# with case specific modifications (e.g. 'freundlich')
# requires 'popt' and 'pcov' variables from curve_fit()
def _reg_conf(ydata, alpha, isoName, popt, pcov):
    # init variables
    n_obs = len(ydata)  # number of data points
    n_params = len(popt)  # number of parameters
    dof = max(0, n_obs - n_params)  # degrees of freedom
    tval = t.ppf(1.0-alpha/2., dof)  # calculate student-t value
    # return popt for lower and upper conf intervals
    upper = []
    lower = []
    # case specific: freundlich
    if isoName == 'freundlich':
        # the upper and lower freundlich confidence interval
        # requires the lowest and highest, respectively,
        # possible 'n', which is the second index of popt
        # (i indicates the confidence intrval parameter
        # currently being calculated, i == 1 is n)
        i = 0
        for param, sigma in zip(popt, np.sqrt(np.diag(pcov))):
            if i == 1:
                if popt[0] >= 0:
                    upper.append(param - sigma * tval)
                    lower.append(param + sigma * tval)
                if popt[0] < 0:
                    upper.append(param + sigma * tval)
                    lower.append(param - sigma * tval)
            else:
                upper.append(param + sigma * tval)
                lower.append(param - sigma * tval)
                i = i + 1
        return upper, lower
    # case general
    else:
        for param, sigma in zip(popt, np.sqrt(np.diag(pcov))):
            upper.append(param + sigma * tval)
            lower.append(param - sigma * tval)
        return upper, lower


# check if regression is statistically different from zero
# H0: regression conf intrvl contains zero for all values of x
# Ha: regression conf intrvl does not contain zero for all values of x
def _chk_reg_diff_zero(func, x, upper, lower, chk_pts=50):
    # 'number of chk_pts' from min to max of (x) x-axis values to test
    x_LReg = np.linspace(np.min(x), np.max(x), chk_pts)
    # initialize conf intervals lists
    y_LREG_upper = []
    y_LREG_lower = []
    # input conf intervals based upon upper and lower popt values
    for val in x_LReg:
        y_LREG_upper.append(func(val, *upper))
        y_LREG_lower.append(func(val, *lower))
    # Test H0 rejection validity; reject null hypothesis?
    # False -> unable to reject; True -> able to reject
    for up, low, in zip(y_LREG_upper, y_LREG_lower):
        if not up >= 0 >= low:
            return True
    return False


# check if confidenc interval at point x
# or list of points in iterable x
# contains zero
def _conf_cont_zero(func, x, upper, lower):
    if hasattr(x, '__iter__'):
        for val in x:
            if func(x, *upper) >= 0 >= func(x, *lower):
                return True
        return False
    else:
        if func(x, *upper) >= 0 >= func(x, *lower):
            return True
        else:
            return False
