""" statstical functions """

# third party modules
import numpy as np
from scipy.stats.distributions import t


# Asymptotic variable initialization
def InitAsymVars(xdata, params, alpha=0.05):
    n_obs = len(xdata)
    n_params = len(params)
    dof = max(0, n_obs - n_params)
    tval = t.ppf(1 - alpha / 2, dof)
    upper = []
    lower = []
    return n_obs, n_params, dof, tval, upper, lower


# Asymptotic confidence interval calculation
def RegConfAsym(xdata, params, covar, alpha=0.05):
    n_obs, n_params, dof, tval, upper, lower =\
        InitAsymVars(xdata, params, alpha)
    for param, sigma in zip(params, np.sqrt(np.diag(covar))):
        upper.append(param + sigma * tval)
        lower.append(param - sigma * tval)
    return {'upper': upper,
            'lower': lower}
