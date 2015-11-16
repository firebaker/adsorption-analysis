"""adsorption.py"""

# third party modules
import numpy as np

# adsorption-analysis modules
import isotherm
from validate import validateInput, checkBehavior


class AdsorptionAnalysis(object):
    """A collection of isotherm class methods to analyze the adsorption
    properties of sample materials (samples).

    Arguments:
    data: A 2-array list of float64 (or convertable numeric),
          data[0] represents non-adsobed phase analyte concentration
          data[1] represents adsorbed phase analyte concentration
    alpha: A float (0.0-1.0), represents alpha level for subsequent
           calculations
    linear: Initial popt0 values designated by the User; "[Kd]"
    freundlich: see linear; "[Kf, n]"
    langmuir: see langmuir; "[Qmax, Kl]"

    Sample attributes:
        data: see input arguments
        alpha: see input arguments
        error: Either FALSE, a string, or a list of strings. If it
               evaluates to TRUE further calculations cease to occur.
        ***isotherm: This is an abstract base class for isotherm objects.
                     In other words, this attribute doesn't exist as an
                     explicit object, instead it acts as a parent class
                     for the real isotherm objects (e.g. linear,
                     freundlich, langmuir). Isotherm objects contain
                     pertinent attributes of the fitted isotherm object:
        *** see isotherm.py for detailed explanation of isotherm
            attributes
        *** Isotherm attributes:
                isotherm.pop0
                isotherm.popt
                isotherm.pcov
                isotherm.conf_upper
                isotherm.conf_lower
                isotherm.pred_upper
                isotherm.pred_lower
                isotherm.SSR
                isotherm.AIC
                isotherm.BIC
        linear: see *isotherm
        freundlich: see *isotherm
        langmuir: see *isotherm
    """

    def __init__(
            self, data, alpha=0.05,
            linear=None, freundlich=None, langmuir=None):
        self.data = data
        self.alpha = alpha
        self.error = validateInput(self.data, self.alpha)
        if self.error:
            return None
        self.data = [np.float64(data[0]), np.float64(data[1])]
        self.behavior, self.error = checkBehavior(data, alpha)
        if self.behavior != 'REMOVAL':
            return None
        self.linear = isotherm.Linear(
            self.data, self.alpha, linear, requireInputValidation=False)
        self.freundlich = isotherm.Freundlich(
            self.data, self.alpha, freundlich, requireInputValidation=False)
        self.langmuir = isotherm.Langmuir(
            self.data, self.alpha, langmuir, requireInputValidation=False)

    def bestfit(self, selection="AIC"):
        """Return the isotherm with the lowest selection criteria values.
        Selection criteria may be:
        AIC: Akaike information criterion (default)
        BIC: Bayesian information criterion
        SSR: Sum of squared residuals
        """
        # TODO
