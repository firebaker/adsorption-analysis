""" background functions strictly related to Adsorption_Analysis.py """

# third party modules
from scipy.optimize import curve_fit

# adsorption-analysis modules
import checkPass as cP
import stats as sts


""" bare adsorption isotherm equations """


# Linear isotherm equation
def linearIsotherm(x, Kd):
    return Kd * x


# Freundlich isotherm equation
def freundlichIsotherm(x, Kf, n):
    return Kf * x**(1/n)


# Langmuir isotherm equation
def langmuirIsotherm(x, Qmax, Kl):
    return Qmax * Kl * x / (1 + Kl * x)


""" Check for isotherm specific errors """


# check for isotherm specific statistical or theory errors
# return dictionary warning: True/Flase;
# if warning: True return plot_capability: True/Flase, and message
# if plot_cpability: True, return modified upper or lower conf
# values
def isothermSpecificCheck(isoName, popt, upper, lower):

    # if user-defined, return 'unable to check warning'
    if isoName == 'user-defined':
        return {'warning': True,
                'plot_capable': True,
                'message': """unable to check or correct for
                isotherm specific errors on user-defined isotherm
                algorithms"""}

    # linear specific check
    elif isoName == 'linear':
        # check if Kd > 0 (indicates adsorption)
        if popt[0] <= 0:
            return {'warning': True,
                    'plot_capable': False,
                    'message': """"unable to plot: The Kd < 0
                    indicates non-adsorption"""}
        # check if confidence interval contains zero
        if sts._conf_cont_zero(linearIsotherm, 0.001, upper, lower):
            return {'warning': True,
                    'plot_capable': False,
                    'message': """unable to plot: Although the data is
                    statistically different from zero, the fitted linear
                    isotherm is not"""}
        return {'warning': False}

    # freundlich specific checks
    elif isoName == 'freundlich':
        # check if 'n' >= 1
        if not popt[1] >= 1.0:
            return {'warning': True,
                    'plot_capable': False,
                    'message': """unable to plot: Freundlich theory
                    requires that 'n' >= 1"""}
        # check if Kf > 0 (indicates adsorption)
        if popt[0] <= 0:
            return {'warning': True,
                    'plot_capable': False,
                    'message': """"unable to plot:
                    'Kf' (={0}) < 0
                    indicates non-adsorption"""
                    .format(popt[0])}
        # check lower Kf > 0
        # if not then the confidence interval
        if lower[0] <= 0:
            return {'warning': True,
                    'plot_capable': False,
                    'message': """unable to plot:
                    lower confidence interval Kf (={0}) < 0;
                    indicates adsorption is
                    non-statistically significant"""}
        # check upper n >= 1
        # if less, indicates very poor fit
        if upper[1] < 1:
            return {'warning': True,
                    'plot_capable': False,
                    'message': """ unable to plot:
                    upper n (={0}) < 1 indicates very poor fit"""}
        return {'warning': False}

    # langmuir specific checks
    elif isoName == 'langmuir':
        # check if Qmax, Kl > 0 (indicates adsorption)
        if popt[0] <= 0 or popt[1] <= 0:
            return {'warning': True,
                    'plot_capable': False,
                    'message': """unable to plot:
                    Qmax (={0}) and/or Kf (={1}) < 0
                    indicates non-adsorption"""
                    .format(popt[0], popt[1])}
        # check lower Qmax and Kl are both > 0
        if lower[0] <= 0 or lower[1] <= 0:
            return {'warning': True,
                    'plot_capable': False,
                    'message': """unable to plot:
                    lower Qmax (={0}) and/or Kf (={1}) < 0
                    indicates adsorption is
                    non-statistically significant"""
                    .format(lower[0], lower[1])}
        return {'warning': False}


# check sorption input
def checkSorptionInput(x, y):
    # set variables for input
    errorcheck = 0
    errorlist = []
    # verify values are numeric
    for array in (x, y):
        try:
            cP._verify_numeric(array)
        except TypeError:
            errorcheck = errorcheck + 1
            errorlist.append(
                'TypeError: not all input values are numeric')
    # verify arrays of equal length
    try:
        cP._verify_similar_length(x, y)
    except:
        errorcheck = errorcheck + 1
        errorlist.append(
            'ArrayError: input arrays are not of equal length')
    # stop function if errorcheck > 0
    if errorcheck > 0:
        print 'Error in -> adsorptionAnalysis:'
        for error in errorlist:
            print error
        return None


# statistically check for adsorption occurrence
# H0: data not different from zero
# Ha: data > zero
def checkStatAdsorption(x, y, alpha=0.05):
    # fit linear regression to data
    s_popt, s_pcov = curve_fit(sts._linReg, x, y)
    if s_popt[0] <= 0:
        return False
    s_upper, s_lower = sts._reg_conf(y, alpha, 'linear', s_popt, s_pcov)
    # test rejection of H0: LinReg not different from zero line
    # False -> unable to reject; True -> able to reject
    if sts._chk_reg_diff_zero(
            sts._linReg, x, s_upper, s_lower):
        return True
    else:
        return False
