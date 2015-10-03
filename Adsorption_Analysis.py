""" main functions for adsorption-analysis package"""

# third party modules and functions
import numpy as np
from scipy.optimize import curve_fit

# adsorption-analysis modules
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
        analyte_behavior="REMOVAL",
        isoName='user-defined'):
    if checkInput:
        bAA.checkSorptionInput(x, y)
        x = np.float64(x)
        y = np.float64(y)
        analyte_behavior = bAA.checkStatRemoval(x, y, alpha)
    if analyte_behavior == "REMOVAL":
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
        # confidence interval -> upper, lower
        upper, lower = sts._reg_conf(y, alpha, isoName, popt, pcov)
        # check for isotherm specific errors
        iso_spec_warn = bAA.isothermSpecificCheck(
            isoName, popt, upper, lower)
        return_dict = {}
        if iso_spec_warn['warning']:
            return_dict.update({'warning': iso_spec_warn['message']})
            if not iso_spec_warn['plot_capable']:
                return return_dict
        # Sum of Squared Residuals and Asike Index Criterion
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
        checkInput=True, analyte_behavior="REMOVAL"):
    p0 = (Kd)
    return fitIsotherm(
        bAA.linearIsotherm, x, y,
        p0,
        alpha, checkInput, analyte_behavior, isoName='linear')


# fit freundlich isotherm, return popt, pcov, SSR
def fitFreundlich(
        x, y,
        Kf=1, n=1, alpha=0.05,
        checkInput=True, analyte_behavior="REMOVAL"):
    p0 = (Kf, n)
    return fitIsotherm(
        bAA.freundlichIsotherm, x, y,
        p0,
        alpha, checkInput, analyte_behavior, isoName='freundlich')


# fit langmuir isotherm, return popt, pcov, SSE
def fitLangmuir(
        x, y,
        Qmax=1, Kl=1,
        alpha=0.05, checkInput=True, analyte_behavior="REMOVAL"):
    p0 = (Qmax, Kl)
    return fitIsotherm(
        bAA.langmuirIsotherm, x, y,
        p0,
        alpha, checkInput, analyte_behavior, isoName='langmuir')


# attempt to fit all isotherms available in module;
# return dict of values and general analyte behavior
def AdsorptionAnalysis(
        x, y, K_,
        Kd=1, Kf=1, n=1, Qmax=1, Kl=1,
        alpha=0.05, checkInput=True, analyte_behavior="REMOVAL"):
    if checkInput:
        bAA.checkSorptionInput(x, y)
        checkInput = False
        analyte_behavior = bAA.checkStatRemoval(x, y, alpha)
    x_iso = np.float64(x)
    y_iso = np.float64(y)
    # print whether addition, removal, or no diff from zero
    if analyte_behavior == "REMOVAL":
        if K_:
            Kd = Kf = Kl = K_
        output = {'linear': fitLinear(
                            x_iso, y_iso, Kd, alpha,
                            checkInput, analyte_behavior),
                  'freundlich': fitFreundlich(
                                x_iso, y_iso, Kf, n, alpha,
                                checkInput, analyte_behavior),
                  'langmuir': fitLangmuir(
                              x_iso, y_iso, Qmax, Kl, alpha,
                              checkInput, analyte_behavior),
                  'analyte_behavior': analyte_behavior}
        for isotherm in output:
            if 'warning' in output[isotherm]:
                print "{0} warning: {1}".format(
                isotherm,
                output[isotherm]['warning'])
        if 'warning' in output['linear']:
            if output['linear']['warning'] == """UNABLE TO PLOT:
                    Although the data is statistically different from
                    zero, the fitted linear isotherm is not""":
                        analyte_behavior = "NOCHANGE"
            else:
                analyte_behavior = "ADDITTION"
            output['analyte_behavior'] = analyte_behavior
    else:
        output = {'warning':
                  """checkStatAdsorption fail for alpha = {0:.2}
                  analyte behavior = {1}"""
                  .format(alpha, analyte_behavior), 
                  'analyte_behavior': analyte_behavior}
        print 'warning: {0:s}'.format(output['warning'])
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
