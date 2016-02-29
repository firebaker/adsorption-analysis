"""adsorption.py"""

# adsorption-analysis modules
import isotherm
import AAvalidate as val
import AAstats


class AdsorptionAnalysis(object):
    """A collection of isotherm class methods to analyze the adsorption
    properties of sample materials (samples).

    fit methods is nelder-mead method, and then a second fit using leastsqr
    according to lmfit. see isotherm.py for more information

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
        Isotherm: abstract base class for child Isotherm classes
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

    def __init__(
            self, data, alpha=0.05,
            linear=None, freundlich=None, langmuir=None):

        # validate input
        try:
            val.validateData(data)
            val.validateAlpha(alpha)
        except val.InputError as ie:
            print(ie)
            return None

        # intialize fitting with Linear isotherm
        self.Linear = isotherm.Linear(
            data, userParams=linear, validateInput=False)

        # test for linear model validity
        if not self.Linear.modelValidity:
            self.LinearFailCode()
            return None

        # test for statistical significance (Asymptotic Confidence Interval)
        Kd = [self.Linear.isoModelResult.params.valuesdict()['Kd']]
        covar = self.Linear.isoModelResult.covar
        confIntrvl = AAstats.RegConfAsym(data[0], Kd, covar, alpha)
        if confIntrvl['lower'][0] <= 0:
            self.LinearFailCode
            return None

        # fit the other isotherms
        self.Freundlich = isotherm.Freundlich(
            data, userParams=freundlich, validateInput=False)
        self.Langmuir = isotherm.Langmuir(
            data, userParams=langmuir, validateInput=False)

    def LinearFailCode(self):
        """In the event that the best fit Linear isotherm fails either
        theory or statistical tests, run this"""

        message = """invalid linear fit;
        To force a {0} model fit, call isotherm.{0}() explicitly"""
        self.Freundlich = message.format("Freundlich")
        self.Langmuir = message.format("Langmuir")

    def bestfit(self, selection="aic"):
        """Return the isotherm with the lowest selection criteria value.
        Selection criteria may be:
        chisqr: sum((residuals array) ** 2)
        aic: Akaike information criterion (default)
        bic: Bayesian information criterion
        """

        # check Linear isotherm
        error_message = "unable to compute bestfit isotherm; \n"
        if not self.Linear.modelValidity:
            return error_message.append("self.Linear.modelValidty = False")

        # initialize selection with linear isotherm
        isotherm = 'Linear'
        slctn_val = getattr(self.Linear.isoModelResult, selection)

        # compare other isotherm fits
        isotherms = [self.Freundlich, self.Langmuir]
        for iso in isotherms:
            if not iso.modelValidity:
                continue
            iso_slctn_val = getattr(iso.isoModelResult, selection)
            if iso_slctn_val < slctn_val:
                isotherm = iso.__name__
                slctn_val = iso_slctn_val
        return isotherm
