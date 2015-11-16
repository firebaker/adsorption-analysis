"""isotherm.py"""

# Python standard library modules and functions
from abc import ABCMeta, abstractmethod

# Third party modules and functions
import numpy as np
from scipy.optimize import curve_fit

# adsorption-analysis modules
from validate import *
import stats as sts


class Isotherm(object):
    """An abstract base class of functions to fit isotherms to
    data and return pertinent results of the fitted isotherm -
    Isotherm objects.

    Private class attributes:
    _isoName: string; isotherm name
    _isoVars: list of string; isotherm variables names
    _isoDefault: list of int; isotherm default input variables
    _Isotherm(): float64; returns output of isotherm equation
    _IsoConf(): returns upper and lower confidence interval for
                 the best fit calculated isotherm to the data
    _IsoPred(): NotImplemented

    Public isotherm object attributes:
    ***Note: Attribute existence is conditional - see code
    data: A 2-array list of user data input; should be numeric;
          data[0] represents mobile phase analyte concentration;
          data[1] represents imobile phase analyte concentration
    alpha: A float (0.0-1.0), represents alpha level for subsequent
           calculations
    error: Either False, a string, or a list of strings. If it evaluates
           to True, further calculations cease to occur
    popt0: The initial values for scipy.optimize.curve_fit functions
    userWarning:  A non-interrupting user warning. Either False or a
                  string. Currently used only to convey to the user
                  warnings concerning the user's designated popt0. If it
                  evaluates to TRUE, default popt0 values are used.
    popt:  Best fit popt values determined by scipy.optimize.curve_fit
    pcov:  popt covariance matrix determined by
           scipy.optimize.curve_fit
    conf_upper:  Upper confidence interval popt values
    conf_lower:  Lower confidence interval popt values
    pred_upper:  Upper prediction interval popt values
    pred_lower:  Lower prediction interval popt values
    SSR:  Sum of squared Residuals for best fit curve
    AIC:  Akiake information criterion for best fit curve
    BIC:  Bayesian information criterion for best fit curve
    """

    __metaclass__ = ABCMeta

    # __isoName should be the name of the isotherm
    _isoName = NotImplementedError

    # __isoVars should be the names of the unkown isotherm variables
    _isoVars = NotImplementedError

    # __isoDefault should be the default initial isotherm values (popt0)
    _isoDefault = NotImplementedError

    def __init__(
            self, data, alpha=0.05, popt0User=None,
            requireInputValidation=True):
        if requireInputValidation:
            self.data = data
            self.alpha = alpha
            self.error = validateInput(self.data, self.alpha)
            if self.error:
                return None
            self.data = [np.float64(data[0]), np.float64(data[1])]
            self.behavior, self.error = checkBehavior(data, alpha)
            if self.behavior != 'REMOVAL':
                return None
        self.popt0User = popt0User
        if self.popt0User:
            self.userWarning, self.popt0 = validate_popt0User(
                self.popt0User, self._isoVars, self._isoDefault)
        else:
            self.userWarning = """popt0User is set to None;
            popt0 was set to the default values: {0}""".format(
                self._isoDefault)
            self.popt0 = self._isoDefault
        """firebaker - 8/12/15: I like curve_fit, but I ultimately do
        not know what is happening (or why it would break), therefore
        the try/catch statement"""
        try:
            self.popt, self.pcov = curve_fit(
                self._Isotherm, data[0], data[1], self.popt0)
            """firebaker - 8/14/15: some situations arise where Inf
            (NaN in a numpy array) is returned, therefore
            check_for_NaN"""
            check_for_NaN(self.popt)
        except:
            self.error = 'unable to fit %s isotherm'.format(
                self._isoName)
            return None
        self.conf_upper, self.conf_lower = _isoConf
        # self.pred_upper, self.pred_lower =
        self.SSR = sts._regSSR(
            self._Isotherm(), self.data[0], self.data[1], self.popt)
        self.AIC = sts._regAIC(
            self.Isotherm, self.data[0], self.data[1], self.popt)
        # self.BIC

    @abstractmethod
    def _Isotherm():
        """Should return output of isotherm equation
        (e.g. for Langmuir - return Qmax * Kl * x / (1 + Kl * x))
        (ex. see clases Linear, Freundlich, and Langmuir)
        """
        raise NotImplementedError

    @abstractmethod
    def _IsoConf(self):
        """Should return upper and lower confidence interval for
        the best fit calculated isotherm to the data
        (e.g. for Linear - return sts._reg_conf_asym(
        self.data[0], self.alpha, self.popt, self.pcov)
        (e.g. for Freundlich & Langmuir - return sts._reg_conf_monte(
        self.data[0], self.alpha, self.popt, self.pcov, self._Isotherm)
        """
        raise NotImplementedError


#    @abstractmethod
#    def _IsoPred(self):
#        """Should return upper and lower prediction interval for
#        the best fit calculated isotherm to the data
#        """
#        raise NotImplementedError


class Linear(Isotherm):

    _isoName = "Linear"
    _isoVars = ["Kd"]
    _isoDefault = [1.0]

    def _Isotherm(x, Kd):
        return Kd * x

    def _IsoConf(self):
        return sts._reg_conf_asym(
            self.data[0], self.alpha, self.popt, self.pcov)


class Freundlich(Isotherm):

    _isoName = "Freundlich"
    _isoVars = ["Kf, n"]
    _isoDefault = [1, 1]

    def _Isotherm(x, Kf, n):
        return Kf * x ** (1 / n)

    def _IsoConf(self):
        return sts._reg_conf_monte(
            self.data[0], self.alpha, self.popt, self.pcov)


class Langmuir(Isotherm):

    _isoName = "Langmuir"
    _isoVars = ["Qmax, Kl"]
    _isoDefault = [1, 1]

    def _Isotherm(x, Qmax, Kl):
        return Qmax * Kl * x / (1 + Kl * x)

    def _IsoConf(self):
        return sts._reg_conf_monte(
            self.data[0], self.alpha, self.popt, self.pcov)
