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
        # Kd must be > 0 (indicates adsorption)
        if popt[0] <= 0:
            return {'warning': True,
                    'plot_capable': False,
                    'message': """UNABLE TO PLOT:
                    Kd ({0:.3}) < 0; indicates non-adsorption"""
                    .format(popt[0])}
        # confidence interval must not contain zero
        if sts._conf_cont_zero(linearIsotherm, 0.001, upper, lower):
            return {'warning': True,
                    'plot_capable': False,
                    'message': """UNABLE TO PLOT:
                    Although the data is statistically different from
                    zero, the fitted linear isotherm is not"""}
        return {'warning': False}

    # freundlich theory specific checks
    elif isoName == 'freundlich':
        # 'n' must be >= 1
        if not popt[1] >= 1.0:
            return {'warning': True,
                    'plot_capable': False,
                    'message': """UNABLE TO PLOT:
                    Freundlich theory requires that
                    'n' ({0:.3}) >= 1"""
                    .format(popt[1])}
        # Kf must be > 0
        if popt[0] <= 0:
            return {'warning': True,
                    'plot_capable': False,
                    'message': """"UNABLE TO PLOT:
                    'Kf' ({0:.3}) < 0
                    indicates non-adsorption"""
                    .format(popt[0])}
        # lower Kf must be > 0
        if lower[0] <= 0:
            return {'warning': True,
                    'plot_capable': False,
                    'message': """UNABLE TO PLOT:
                    lower confidence interval Kf ({0:.3}) < 0;
                    indicates adsorption is
                    non-statistically significant"""
                    .format(lower[0])}
        # upper n must be >= 1
        if upper[1] < 1:
            return {'warning': True,
                    'plot_capable': False,
                    'message': """UNABLE TO PLOT:
                    upper n ({0:.3}) < 1;
                    indicates poor fit"""
                    .format(upper[1])}
        return {'warning': False}

    # langmuir theory specific checks
    elif isoName == 'langmuir':
        # Qmax & Kl must be > 0
        if popt[0] <= 0 or popt[1] <= 0:
            return {'warning': True,
                    'plot_capable': False,
                    'message': """UNABLE TO PLOT:
                    Qmax ({0:.3}) and/or Kf ({1:.3}) < 0
                    indicates non-adsorption"""
                    .format(popt[0], popt[1])}
        # lower Qmax & Kl must be > 0
        if lower[0] <= 0 or lower[1] <= 0:
            return {'warning': True,
                    'plot_capable': False,
                    'message': """UNABLE TO PLOT:
                    lower Qmax ({0:.3}) and/or Kf ({1:.3}) < 0
                    indicates adsorption is
                    non-statistically significant"""
                    .format(lower[0], lower[1])}
        return {'warning': False}


# check sorption input for logical requirements
def checkSorptionInput(x, y):
    errorcheck = 0
    errorlist = []
    for array in (x, y):
        try:
            cP._verify_numeric(array)
        except TypeError:
            errorcheck += 1
            errorlist.append(
                'TypeError: not all input values are numeric')
    try:
        cP._verify_similar_length(x, y)
    except:
        errorcheck += 1
        errorlist.append(
            'ArrayError: input arrays are not of equal length')
    if errorcheck > 0:
        print 'Error in -> adsorptionAnalysis:'
        for error in errorlist:
            print error
        return None


# statistically check for removal occurrence
# H0: data not different from zero
# Ha: data > zero
# test rejection of H0: LinReg not different from zero line
# False -> unable to reject; True -> able to reject
def checkStatRemoval(x, y, alpha=0.05):
    s_popt, s_pcov = curve_fit(sts._linReg, x, y)
    s_upper, s_lower = sts._reg_conf(y, alpha, 'linear', s_popt, s_pcov)
    if sts._chk_reg_diff_zero(
            sts._linReg, x, s_upper, s_lower):
        if s_popt[0] <= 0:
            return 'ADDITION'
        return 'REMOVAL'
    else:
        return 'NOCHANGE'
