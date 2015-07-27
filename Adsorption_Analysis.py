import numpy as np

from checkPass import _verify_numeric
from dataStructure import _chk_key_dict
from Isotherms import fitLinear
from Isotherms import fitFreundlich
from Isotherms import fitLangmuir


def adsorptionAnalysis(x, y, Kd=1, Kf=1, n=1, Smax=1, Kl=1):
    # set variables for input
    errorcheck = 0
    errorlist = []
    # verify values are numeric
    for array in (x, y):
        try:
            _verify_numeric(array)
        except TypeError:
            errorcheck = errorcheck + 1
            errorlist.append(
                'TypeError: not all input values are numeric')
    # verify arrays of equal length
    if len(x) != len(y):
        errorcheck = errorcheck + 1
        errorlist.append(
            'ArrayError: input arrays are not of equal length')
    # stop function if errorcheck > 1
    if errorcheck >= 1:
        print 'Error in -> adsorptionAnalysis:'
        for error in errorlist:
            print error
        return None
    # convert input values to float64
    x_iso = np.float64(x)
    y_iso = np.float64(y)
    # create output dictionary
    output = {}
    _chk_key_dict((
        'linear',
        'freundlich',
        'langmuir'), output)
    # fit isotherms to data; return to output dictionary
    output['linear'].update(fitLinear(x_iso, y_iso, Kd))
    output['freundlich'].update(fitFreundlich(x_iso, y_iso, Kf, n))
    output['langmuir'].update(fitLangmuir(x_iso, y_iso, Smax, Kl))
    for isotherm in output:
        if 'warning' in output[isotherm]:
            print 'warning %s' % output[isotherm]['warning']
    return output
