""" main group of functions for adsorption analysis module"""


import numpy as np
from scipy.optimize import curve_fit
from dataStructure import _chk_key_dict
from checkPass import _check_nan
import background_Adsorption_Analysis as bAA


""""
algorithms to fit isotherms to experimental results using
scipy.optimize curve_fit; return popt and pcov variables
for each isotherm in an output dictionary
"""


# fit linear isotherm
def fitLinear(x, y, Kd=1):
    bAA.checkAdsorptionInput(x, y)
    try:
        popt, pcov = curve_fit(bAA.linearIsotherm, x, y, p0=(Kd))
        _check_nan(popt)
        return {'popt': popt, 'pcov': pcov}
    except:
        return {'warning': 'unable to fit linear isotherm'}


# fit freundlich isotherm
def fitFreundlich(x, y, Kf=1, n=1):
    bAA.checkAdsorptionInput(x, y)
    try:
        popt, pcov = curve_fit(bAA.freundlichIsotherm, x, y, p0=(Kf, n))
        _check_nan(popt)
        return {'popt': popt, 'pcov': pcov}
    except:
        return {'warning': 'unable to fit freundlich isotherm'}


# fit langmuir isotherm
def fitLangmuir(x, y, Qmax=1, Kl=1):
    bAA.checkAdsorptionInput(x, y)
    try:
        popt, pcov = curve_fit(bAA.langmuirIsotherm, x, y, p0=(Qmax, Kl))
        _check_nan(popt)
        return {'popt': popt, 'pcov': pcov}
    except:
        return {'warning': 'unable to fit langmuir isotherm'}


def adsorptionAnalysis(x, y, Kd=1, Kf=1, n=1, Qmax=1, Kl=1):
    bAA.checkAdsorptionInput(x, y)
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
    output['langmuir'].update(fitLangmuir(x_iso, y_iso, Qmax, Kl))
    for isotherm in output:
        if 'warning' in output[isotherm]:
            print 'warning: %s' % output[isotherm]['warning']
    return output


# guess Qmax, optional for fitLangmuir()
def guessQmax(x, y, xndpts=3):
    xlist = sorted(range(len(x)), key=lambda i: x[i])[-xndpts:]
    ylist = []
    for index in xlist:
        ylist.append(y[index])
    return np.mean(ylist)
