""" main functions for adsorption-analysis package"""

# third party modules and functions
import numpy as np
from scipy.optimize import curve_fit

# adsorption-analysis modules
import dataStructure as dS
import checkPass as cP
import stats as sts
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
        checkInput=True,
        adsorptionOccurrence=True,
        isoName='user-defined'):
    # init parameter
    if checkInput:
        bAA.checkSorptionInput(x, y)
        x = np.float64(x)
        y = np.float64(y)
        # statistically check for sorption occurrence
        adsorptionOccurrence = bAA.checkStatSorption(x, y, alpha)
    if adsorptionOccurrence:
        # firebaker Note - 8/12/15: I trust curve_fit,
        # but I ultimately do not know what is happening
        # (or why it would break), therefore the try/catch statement
        try:
            popt, pcov = curve_fit(isotherm, x, y, p0)
            # some situations arise where Inf (NaN in a numpy array)
            # is returned, therefore _check_nan is in the
            # try/catch statement
            cP._check_nan(popt)
        except:
            return {'warning': 'unable to fit %s isotherm' % isoName}
        # calculate confidence interval
        upper, lower = sts._reg_conf(y, alpha, isoName, popt, pcov)
        # check for isotherm specific errors
        iso_spec_warn = bAA.isothermSpecificCheck(
            isoName, popt, upper, lower)
        return_dict = {}
        if iso_spec_warn['warning']:
            return_dict.update({'warning': iso_spec_warn['message']})
            if not iso_spec_warn['plot_capable']:
                return return_dict
            else:
                if 'upper' in set(iso_spec_warn):
                    upper = iso_spec_warn['upper']
                elif 'lower' in set(iso_spec_warn):
                    lower = iso_spec_warn['lower']
        # calculate SSR and AIC
        SSR = sts._regSSR(isotherm, x, y, popt)
        AIC = sts._regAIC(isotherm, x, y, popt)
        return_dict.update({'popt': popt, 'pcov': pcov,
                            'SSR': SSR, 'AIC': AIC,
                            'upper': upper, 'lower': lower})
        return return_dict
    else:
        return {'warning': 'sorptionCheck False for alpha = %s' % alpha}


# fit linear isotherm; return popt, pcov, SSR
def fitLinear(
        x, y, Kd=1, alpha=0.05,
        checkInput=True, adsorptionOccurrence=True):
    p0 = (Kd)
    return fitIsotherm(
        bAA.linearIsotherm, x, y,
        p0,
        alpha, checkInput, adsorptionOccurrence, isoName='linear')


# fit freundlich isotherm, return popt, pcov, SSR
def fitFreundlich(
        x, y,
        Kf=1, n=1, alpha=0.05,
        checkInput=True, adsorptionOccurrence=True):
    p0 = (Kf, n)
    return fitIsotherm(
        bAA.freundlichIsotherm, x, y,
        p0,
        alpha, checkInput, adsorptionOccurrence, isoName='freundlich')


# fit langmuir isotherm, return popt, pcov, SSE
def fitLangmuir(
        x, y,
        Qmax=1, Kl=1,
        alpha=0.05, checkInput=True, adsorptionOccurrence=True):
    p0 = (Qmax, Kl)
    return fitIsotherm(
        bAA.langmuirIsotherm, x, y,
        p0,
        alpha, checkInput, adsorptionOccurrence, isoName='langmuir')


# attempt to fit all isotherms available in module;
# return dict of values
def AdsorptionAnalysis(
        x, y, K_,
        Kd=1, Kf=1, n=1, Qmax=1, Kl=1,
        alpha=0.05, checkInput=True, adsorptionOccurrence=True):
    if checkInput:
        bAA.checkSorptionInput(x, y)
        checkInput = False
        # statistically check for sorption occurrence
        adsorptionOccurrence = bAA.checkStatAdsorption(x, y, alpha)
    # convert input values to float64
    x_iso = np.float64(x)
    y_iso = np.float64(y)
    # check sorptionOccurrence for t/f
    if adsorptionOccurrence:
        # create isotherm fitting output dictionary
        if K_:
            Kd = Kf = Kl = K_
        output = {}
        dS._chk_key_dict((
            'linear',
            'freundlich',
            'langmuir'), output)
        # fit isotherms to data; return to output dictionary
        output['linear'].update(fitLinear(
            x_iso, y_iso, Kd, alpha,
            checkInput, adsorptionOccurrence))
        output['freundlich'].update(fitFreundlich(
            x_iso, y_iso, Kf, n, alpha,
            checkInput, adsorptionOccurrence))
        output['langmuir'].update(fitLangmuir(
            x_iso, y_iso, Qmax, Kl, alpha,
            checkInput, adsorptionOccurrence))
        for isotherm in output:
            if 'warning' in output[isotherm]:
                print 'warning: %s' % output[isotherm]['warning']
    else:
        output = {'warning':
                  'checkStatAdsorption fail for alpha = %s'
                  .format(alpha)}
        print 'warning: %s' % output['warning']
    return output


# guess K_ for isotherm fitting
# return 'm'_ from fitting  y = m * x + b
def guessK_(x, y, checkInput=True):
    if checkInput:
        bAA.checkSorptionInput(x, y)
        x_iso = np.float64(x)
        y_iso = np.float64(y)
    K_popt, K_pcov = curve_fit(sts._linReg, x_iso, y_iso)
    return K_popt[0]


# guess Qmax, optional for fitLangmuir(); average y value of x-axis
# endpoints
def guessQmax(x, y, xndpts=3, checkInput=True):
    if checkInput:
        bAA.checkSorptionInput(x, y)
        x_iso = np.float64(x)
        y_iso = np.float64(y)
    xlist = sorted(range(len(x_iso)), key=lambda i: x_iso[i])[-xndpts:]
    ylist = []
    for index in xlist:
        ylist.append(y_iso[index])
    return np.mean(ylist)
