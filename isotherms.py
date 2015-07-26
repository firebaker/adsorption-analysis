from scipy.optimize import curve_fit
from .checkpass import _check_nan


""" bare isotherm equations """


# Linear isotherm equation
def linearIsotherm(x, Kd):
    return Kd * x


# Freundlich isotherm equation
def freundlichIsotherm(x, Kf, n):
    return Kf * x**(1/n)


# Langmuir isotherm equation and related functions
def langmuirIsotherm(x, Smax, Kl):
    return Smax * Kl * x / (1 + Kl * x)


"""" algorithms to fit isotherms to experimental results"""


# fit linear isotherm
def fitLinear(x, y, Kd):
    try:
        popt, pcov = curve_fit(linearIsotherm, x, y, p0=(Kd))
        _check_nan(popt)
        return popt, pcov
    except:
        return None, 'warning: unable to fit linear isotherm'


# fit freundlich isotherm
def fitFreundlich(x, y, Kf, n):
    try:
        popt, pcov = curve_fit(freundlichIsotherm, x, y, p0=(Kf, n))
        _check_nan(popt)
        return popt, pcov
    except:
        return None, 'warning: unable to fit linear isotherm'


# fit langmuir isotherm
def fitLangmuir(x, y, Smax, Kl):
    try:
        popt, pcov = curve_fit(langmuirIsotherm, x, y, p0=(Smax, Kl))
        _check_nan(popt)
        return popt, pcov
    except:
        return None, 'warning: unable to fit linear isotherm'
