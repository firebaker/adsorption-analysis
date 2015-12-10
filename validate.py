""" validation and check functions for various objects in
adsorption-analysis.py"""

# Python standard library modules
import math

# third party modules
from scipy.optimize import curve_fit

# adsorption-analysis modules
import stats as sts


# validate values are numeric
def validateNumeric(item):
    if hasattr(item, '__iter__'):
        for val in item:
            if not isinstance(val, (int, float, complex)):
                return False
    else:
        if not isinstance(item, (int, float, complex)):
            return False
    return True


def validateData(data):
    errorList = []
    for array in data:
        if not validateNumeric(array):
            errorList.append(
                'DataObsError: input values are not all numeric')
    if len(data) != 2:
        errorList.append(
            'DataInputError: data needs to have exactly two arrays')
    if len(data[0]) != len(data[1]):
        errorList.append(
            'DataInpuError: input arrays are not of equal length')
    return errorList


def validateAlpha(alpha):
    errorList = []
    if not isinstance(alpha, (int, float, complex)):
        return errorList.append(
            'AlphaInputError: alpha must be a single numeric value')
    if alpha < 0 or alpha > 1:
        errorList.append(
            'AlphaInputError: 0 < alpha < 1 is not true')
    return errorList


def validateInput(data, alpha):
    errorList = validateData(data) + validateAlpha(alpha)
    return errorList


def validate_popt0User(popt0User, isoVars, isoDefault):
    if len(popt0User) != len(isoVars):
        return """There are a different number of
        variables in popt0User than required;
        popt0 was set to the default values: {0}""".format(
            isoDefault), isoDefault
    elif not validateNumeric(popt0User):
        return """popt0User values are not all numeric;
        popt0 was set to the default values: {0}""".format(
            isoDefault), isoDefault
    return False, popt0User


def checkBehavior(data, alpha=0.05):
    x, y = data[0], data[1]
    b_popt, b_pcov = curve_fit(sts._linReg, x, y)
    conf = sts._reg_conf_asym(y, alpha, b_popt, b_pcov)
    conflict_phrase = """analyte behavior does not show adsorption occurance;
    further isotherm computation is terminated"""
    if sts._chk_reg_diff_zero(
            sts._linReg, x, conf['upper'], conf['lower']):
        if b_popt[0] <= 0:
            return 'ADDITION', conflict_phrase
        return 'REMOVAL', False
    return 'NOCHANGE', conflict_phrase


# check if nan exists in item
def check_for_NaN(item):
    if hasattr(item, '__iter__'):
        for val in item:
            if math.isnan(val):
                raise TypeError
    else:
        if math.isnan(item):
            raise TypeError
