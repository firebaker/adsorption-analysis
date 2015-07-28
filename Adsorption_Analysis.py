""" main group of functions for adsorption analysis module"""


import numpy as np
from scipy.optimize import curve_fit
from dataStructure import _chk_key_dict
from checkPass import _check_nan
from stats import _regSSR
import background_Adsorption_Analysis as bAA


""""
algorithms to fit isotherms to experimental results using
scipy.optimize curve_fit; return popt and pcov variables
for each isotherm in an output dictionary
"""


# fit linear isotherm; return popt, pcov, SSR
def fitLinear(x, y, Kd=1, checkInput=True):
    if checkInput:
        bAA.checkAdsorptionInput(x, y)
    try:
        popt, pcov = curve_fit(bAA.linearIsotherm, x, y, p0=(Kd))
        _check_nan(popt)
        SSR = _regSSR(bAA.linearIsotherm, x, y, popt)
        return {'popt': popt, 'pcov': pcov, 'SSR': SSR}
    except:
        return {'warning': 'unable to fit linear isotherm'}


# fit freundlich isotherm, return popt, pcov, SSR
def fitFreundlich(x, y, Kf=1, n=1, checkInput=True):
    if checkInput:
        bAA.checkAdsorptionInput(x, y)
    try:
        popt, pcov = curve_fit(bAA.freundlichIsotherm, x, y, p0=(Kf, n))
        _check_nan(popt)
        return {'popt': popt, 'pcov': pcov}
    except:
        return {'warning': 'unable to fit freundlich isotherm'}


# fit langmuir isotherm, return popt, pcov, SSE
def fitLangmuir(x, y, Qmax=1, Kl=1, checkInput=True):
    if checkInput:
        bAA.checkAdsorptionInput(x, y)
    try:
        popt, pcov = curve_fit(bAA.langmuirIsotherm, x, y, p0=(Qmax, Kl))
        _check_nan(popt)
        return {'popt': popt, 'pcov': pcov}
    except:
        return {'warning': 'unable to fit langmuir isotherm'}


def adsorptionAnalysis(
        x, y, Kd=1, Kf=1, n=1, Qmax=1, Kl=1, checkInput=True):
    if checkInput:
        bAA.checkAdsorptionInput(x, y)
        checkInput = False
    # convert input values to float64
    x_iso = np.float64(x)
    y_iso = np.float64(y)
    # create working dictionary
    working = {}
    _chk_key_dict((
        'linear',
        'freundlich',
        'langmuir'), working)
    # fit isotherms to data; return to working dictionary
    working['linear'].update(fitLinear(
        x_iso, y_iso, Kd, checkInput))
    working['freundlich'].update(fitFreundlich(
        x_iso, y_iso, Kf, n, checkInput))
    working['langmuir'].update(fitLangmuir(
        x_iso, y_iso, Qmax, Kl, checkInput))
    for isotherm in working:
        if 'warning' in working[isotherm]:
            print 'warning: %s' % working[isotherm]['warning']
    output = working
    return output


# guess Qmax, optional for fitLangmuir(); average y value of x-axis
# endpoints
def guessQmax(x, y, xndpts=3):
    xlist = sorted(range(len(x)), key=lambda i: x[i])[-xndpts:]
    ylist = []
    for index in xlist:
        ylist.append(y[index])
    return np.mean(ylist)
