"""adsorption.py"""

# third party modules
import numpy as np

# adsorption-analysis modules
import isotherm
import AAvalidate


class AdsorptionAnalysis(object):
    """A collection of isotherm class methods to analyze the adsorption
    properties of sample materials (samples).

    Arguments:
    data: A 2-array list of float64 (or convertable numeric),
          data[0] represents non-adsobed phase analyte concentration
          data[1] represents adsorbed phase analyte concentration
    linear: Initial param values designated by the user; "[Kd]"
    freundlich: see linear; "[Kf, n]"
    langmuir: see langmuir; "[Qmax, Kl]"

    Sample attributes:
        data: see input arguments
        error: Either FALSE, a string, or a list of strings. If it
               evaluates to TRUE further calculations cease to occur.
        Isotherm: abstract base calss for child Isotherm classes
            Isotherm attributes:
                isotherm.params (see isotherm.py)
                isotherm.userParams (see isotherm.py)
                isotherm.userWarning (see isotherm.py)
                isotherm.minimizedFit (see isotherm.py)
                isotherm.modelValidity (see isotherm.py)
                isotherm.modelValidtyMsg (see isotherm.py)
        Linear: see Isotherm
        Freundlich: see Isotherm
        Langmuir: see Isotherm

    Sample Method:
        bestfit(): method to determine the best fit isotherm;
                   can use chisqr, aic, or bic selection criteria
    """

    def __init__(self, data, alpha=0.05,
            linear=None, freundlich=None, langmuir=None):
        self.data = data
        self.error = AAvalidate.validateInput(self.data)
        if self.error:
            return None
        self.data = [np.float64(data[0]), np.float64(data[1])]
        self.Linear = isotherm.Linear(
            self.data, alpha=alpha, userParams=linear, validateInput=False)
        if not self.Linear.modelValidity:
            msg = """invalid linear fit;
            To force a {0} model fit,
            call isotherm.{0}() explicitly"""
            self.Freundlich = msg.format("Freundlich")
            self.Langmuir = msg.format("Langmuir")
            return None
        self.Freundlich = isotherm.Freundlich(
            self.data, alpha=alpha, userParams=freundlich, validateInput=False)
        self.Langmuir = isotherm.Langmuir(
            self.data, alpha=alpha, userParams=langmuir, validateInput=False)

    def bestfit(self, selection="aic"):
        """Return the isotherm with the lowest selection criteria value.
        Selection criteria may be:
        chisqr: sum((residuals array) ** 2)
        aic: Akaike information criterion - (default)
        bic: Bayesian information criterion
        """

        # check for sample isotherm fitting errors
        error_msg = "unable to compute bestfit isotherm"
        if self.error:
            return error_msg.append("\
                check self.error for more information")
        if not self.Linear.modelValidity:
            return error_msg.append("\
                check self.Linear.modelValidtyMsg for more information")

        # initialize selection with linear isotherm
        isotherm = 'Linear'
        slctn_val = getattr(self.Linear.minimizedFit, selection)

        # compare linear fit to other isotherm fits
        isotherms = [self.Freundlich, self.Langmuir]
#         isotherms = [self.Langmuir]
        for iso in isotherms:
            if not iso.modelValidity:
                continue
            iso_slctn_val = getattr(iso.minimizedFit, selection)
            if iso_slctn_val < slctn_val:
                isotherm = iso.__name__
                slctn_val = iso_slctn_val
        return isotherm
