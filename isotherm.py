"""isotherm.py"""

# Python standard library modules and functions
from abc import ABCMeta, abstractmethod, abstractproperty

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

    Isotherm objects have the following attributes:
        "Attribute existence is conditional - see code"
    data: A 2-array list of user data input; should be numeric;
          data[0] represents mobile phase analyte concentration;
          data[1] represents imobile phase analyte concentration
    alpha:  # TODO
    error:  # TODO
    userWarning:  # TODO
    popt0: # TODO
    popt:  # TODO
    pcov:  # TODO
    conf_upper:  # TODO
    conf_lower:  # TODO
    pred_upper  # TODO
    pred_lower:  # TODO
    SSR:  # TODO
    AIC:  # TODO
    BIC:  # TODO
    """

    __metaclass__ = ABCMeta

    # __isoName should be the name of the isotherm
    __isoName = NotImplementedError

    # __isoVars should be the names of the unkown isotherm variables
    __isoVars = NotImplementedError

    # __isoDefault should be the default initial isotherm values (popt0)
    __isoDefault = NotImplementedError

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
                self.popt0User, self.isoVars, self.__isoDefault)
        else:
            self.userWarning = """popt0User is set to None;
            popt0 was set to the default values: {0}""".format(__isoDefault)
            self.popt0 = self.__isoDefault
        """firebaker - 8/12/15: I like curve_fit, but I ultimately do
        not know what is happening (or why it would break), therefore
        the try/catch statement"""
        try:
            self.popt, self.pcov = curve_fit(
                self.__Isotherm, data[0], data[1], self.popt0)
            """firebaker - 8/14/15: some situations arise where Inf
            (NaN in a numpy array) is returned, therefore
            check_for_NaN"""
            check_for_NaN(self.popt)
        except:
            self.error = 'unable to fit %s isotherm'.format(
                self.__isoName)
            return None
        self.conf_upper, self.conf_lower = __isoConf
        # self.pred_upper, self.pred_lower =
        self.SSR = sts._regSSR(
            self.__Isotherm(), self.data[0], self.data[1], self.popt)
        self.AIC = sts._regAIC(
            self.Isotherm, self.data[0], self.data[1], self.popt)
        # self.BIC

    @abstractmethod
    def __Isotherm(self):
        """Should return output of isotherm equation
        (e.g. for Langmuir - return Qmax * Kl * x / (1 + Kl * x))
        (ex. see clases Linear, Freundlich, and Langmuir)
        """
        raise NotImplementedError

    @abstractmethod
    def __IsoConf(self):
        """Should return upper and lower confidence interval for
        the best fit calculated isotherm to the data
        (e.g. for Linear - return sts._reg_conf_asym(
        self.data[0], self.alpha, self.popt, self.pcov)
        (e.g. for Freundlich & Langmuir - return sts._reg_conf_monte(
        self.data[0], self.alpha, self.popt, self.pcov, self.__Isotherm)
        """
        raise NotImplementedError

    @abstractmethod
    def __IsoPred(self):
        """Should return upper and lower prediction interval for
        the best fit calculated isotherm to the data
        """
        raise NotImplementedError


class Linear(Isotherm):

    __isoName = "Linear"
    __isoVars = ["Kd"]
    __isoDefault = [1.0]

    def __Isotherm(self):
        return Kd * x

    def __IsoConf(self):
        return sts._reg_conf_asym(
            self.data[0], self.alpha, self.popt, self.pcov)
