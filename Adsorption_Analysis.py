""" main group of functions for adsorption analysis module"""


import numpy as np
# note: curve_fit returns popt and pcov as arrays not list
from scipy.optimize import curve_fit
from dataStructure import _chk_key_dict
from checkPass import _check_nan
from stats import _regSSR
from stats import _confSSR
import background_Adsorption_Analysis as bAA


""""
algorithms to fit isotherms to experimental results using
scipy.optimize curve_fit; return popt (list) and pcov (2d list)
variables for each isotherm in an output dictionary
"""


# generic fit isotherm function (optionally user defined);
# return popt, pcov, SSR
def fitIsotherm(
        isotherm,
        x, y,
        p0,
        alpha=0.05,
        checkInput=True, isoName='user-defined'):
    if checkInput:
        bAA.checkAdsorptionInput(x, y)
        x = np.float64(x)
        y = np.float64(y)
    try:
        popt, pcov = curve_fit(isotherm, x, y, p0)
        _check_nan(popt)
        SSR = _regSSR(isotherm, x, y, popt)
        try:
            upper, lower = _confSSR(isotherm, x, y, alpha, popt, pcov)
        except:
            upper, lower = 'unable to fit conf interval'
        # convert popt and pcov from np.array to list
        popt_return = np.array(popt).tolist()
        pcov_return = np.array(pcov).tolist()
        return {'popt': popt_return, 'pcov': pcov_return,
                'SSR': SSR, 'upper': upper, 'lower': lower}
    except:
        return {'warning': 'unable to fit %s isotherm' % isoName}


# fit linear isotherm; return popt, pcov, SSR
def fitLinear(x, y, Kd=1, alpha=0.05, checkInput=True):
    p0 = (Kd)
    return fitIsotherm(
        bAA.linearIsotherm, x, y, p0, alpha, checkInput, isoName='linear')


# fit freundlich isotherm, return popt, pcov, SSR
def fitFreundlich(x, y, Kf=1, n=1, alpha=0.05, checkInput=True):
    p0 = (Kf, n)
    return fitIsotherm(
        bAA.freundlichIsotherm, x, y, p0, alpha, checkInput, 'freundlich')


# fit langmuir isotherm, return popt, pcov, SSE
def fitLangmuir(x, y, Qmax=1, Kl=1, alpha=0.05, checkInput=True):
    p0 = (Qmax, Kl)
    return fitIsotherm(
        bAA.langmuirIsotherm, x, y, p0, alpha, checkInput, isoName='langmuir')


# fit linear, freundlich, and langmuir, return dict of values
def adsorptionAnalysis(
        x, y, Kd=1, Kf=1, n=1, Qmax=1, Kl=1, alpha=0.05, checkInput=True):
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
        x_iso, y_iso, Kd, alpha, checkInput))
    working['freundlich'].update(fitFreundlich(
        x_iso, y_iso, Kf, n, alpha, checkInput))
    working['langmuir'].update(fitLangmuir(
        x_iso, y_iso, Qmax, Kl, alpha, checkInput))
    for isotherm in working:
        if 'warning' in working[isotherm]:
            print 'warning: %s' % working[isotherm]['warning']
    output = working
    return output


# guess Qmax, optional for fitLangmuir(); average y value of x-axis
# endpoints
def guessQmax(x, y, xndpts=3, checkInput=True):
    if checkInput:
        bAA.checkAdsorptionInput(x, y)
    xlist = sorted(range(len(x)), key=lambda i: x[i])[-xndpts:]
    ylist = []
    for index in xlist:
        ylist.append(y[index])
    return np.mean(ylist)
