""" statstical functions """


# Python standard library modules and functions
from operator import itemgetter

# third party modules
import random
import time
import numpy as np
import math
from scipy.stats.distributions import t
from scipy.optimize import curve_fit


# linearRegresion
def _linReg(x, m=1, b=1):
    return m * x + b


# calculate linear/non-linear regression SSR (sum of squared residuals)
# requires 'popt' variable from curve_fit()
def _regSSR(func, xdata, ydata, popt):
    residuals = ydata - func(xdata, *popt)
    SSR = sum(residuals ** 2)
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


# Asymptotic confidence interval calculation (good for linear regression)
def _reg_conf_asym(xdata, alpha, popt, pcov):
    n_obs = len(xdata)
    n_params = len(popt)
    dof = max(0, n_obs - n_params)
    tval = t.ppf(1.0 - alpha / 2., dof)
    upper = []
    lower = []
    for param, sigma in zip(popt, np.sqrt(np.diag(pcov))):
        upper.append(param + sigma * tval)
        lower.append(param - sigma * tval)
    return upper, lower


# Monte Carlo simulation confidence interval calculation
def _reg_conf_monte(xdata, alpha, popt, pcov, func,
                    seed=time.gmtime(), sim=1000):
    random.seed(seed)
    popt_sigma = np.sqrt(np.diag(pcov))
    sim_list = []
    for i in range(1, sim):
        ydata_sim = []
        for x in xdata:
            sim_popt = []
            for mu, sigma in zip(popt, popt_sigma):
                sim_popt.append(random.gauss(mu, sigma))
            ydata_sim.append(func(x, *sim_popt))
        sim_popt, sim_pcov = curve_fit(func, xdata, ydata_sim, popt)
        sim_list.append((sim_popt, func(max(xdata), *sim_popt)))
    sorted_sim_list = sorted(sim_list, key=lambda x: x[1])
    half_sim_alpha = int(math.ceil(alpha * sim / 2))
    del sorted_sim_list[0:half_sim_alpha]
    del sorted_sim_list[-1:-half_sim_alpha]
    return sorted_sim_list[0][0], sorted_sim_list[-1][0]


# Monte Carlo simulation prediction interval calculation
def _reg_pred_monte(xdata, alpha, popt, pcov, func,
                    seed=time.gmtime(), sim=1000):
    random.seed(seed)
    popt_sigma = np.sqrt(np.diag(pcov))
    sim_data = []
    for simulated in range(1, sim):
        for x in xdata:
            sim_popt = []
            for mu, sigma in zip(popt, popt_sigma):
                sim_popt.append(random.gauss(mu, sigma))
            y_sim = func(x, *sim_popt)
            y_sim_mu = y_sim / func(x, *popt)
            sim_data.append([x, y_sim, y_sim_mu])
    sim_data.sort(key=itemgetter(2))


# check if regression is statistically different from zero
# H0: regression conf intrvl contains zero for all values of x
# Ha: regression conf intrvl does not contain zero for all values of x
# Test H0 rejection validity; reject null hypothesis?
# False -> unable to reject; True -> able to reject
def _chk_reg_diff_zero(func, x, upper, lower, chk_pts=50):
    # 'number of chk_pts' from min to max of (x) values to test
    x_LReg = np.linspace(np.min(x), np.max(x), chk_pts)
    y_LREG_upper = []
    y_LREG_lower = []
    for val in x_LReg:
        y_LREG_upper.append(func(val, *upper))
        y_LREG_lower.append(func(val, *lower))
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
