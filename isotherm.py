"""isotherm.py"""

# Python standard library modules and functions
from abc import ABCMeta, abstractmethod
import collections

# Third party modules and functions
import lmfit

# adsorption-analysis modules
import AAvalidate as val


class Isotherm(metaclass=ABCMeta):
    """An abstract base class of isotherm fitting functions.

    Arguments:
    data: A 2-array list of float64 (or convertable numeric),
        * data[0] represents non-adsobed phase analyte concentration
        * data[1] represents adsorbed phase analyte concentration
    userParams: Initial param values designated by the user;
        * see child isotherm classes for specific params
    fitMethod: Either a str or a list of strs, that are used in sequential
        order to fit and refit data (see note below for an explanation of
        the default values and why sequential refitting may be useful);
        * Values must enumerate with what is avaliable in the lmfit pacakge
        algorithms for optimizing curve fits;
        * Default value is ['nelder', 'leastsqr'];
        *** Note: Nelder-Mead minimiztion is generally considered more robust
        at fitting for special model cases than Levenburg-Marquardt
        Least Squares (useful for Freundlich isotherm fitting).
        However, lmfit's nelder method does not provide uncertainties'
        measurements, whereas lmfit's leastsqr method does, therefore we
        implement Nelder-Mead first to obtain robust fits, then a subsequent
        Levenburg-Marquardt fit is implemented to obtain the uncertainty
        values of the Nelder-Mead fit parameters.
    validateInput: A boolean used to check input values;
        * set to false if being called from AdsorptionAnalysis
        module since the input should already have been checked

    Attributes:
    error: conditional upon "validateInput"
        * Either False, a string, or a list of strings;
        * If evaluating to True, further evaluation ceases
    params: initial isotherm parameters for the minimizing function
    userParams: Initial isotherm parameters given by the user
    userWarning: boolean or string; conveys to the user warnings
                 concerning the user's input params "userParams".
    nelder: lmfit's fitting output using the Nelder-Mead method
    leastsqr: lmfit's fitting output using the Levenburg-Marquardt method
    minimizedFit: pointer to leastsqr
    modelValidity: A boolean that tests theoretical impossibilities;
                   If True, tests were passed;
                   If False, tests did not pass
    modelValidtyMsg: A message relaying the cause of model validity
                     failure
    """

    default_fitMethods = collections.OrderedDict(
        [('nelder', {}), ('leastsqr', {})])

    def __init__(
            self, data, userParams=None,
            fitMethods=default_fitMethods,
            validateInput=True):

        # validate input
        if validateInput:
            try:
                val.validateData(data)
                if userParams:
                    val.validate_userParams(userParams, self.IsothermFunc)
            except val.InputError as ie:
                print(ie)
                return None
            except Exception as inst:
                print(type(inst))
                print(inst)
                return None

        # initiate lmfit's wrapper around the isotherm function
        isoModel = lmfit.Model(self.IsothermFunc)

        # replace inititial params with userParams
        if userParams:
            for key in userParams:
                isoModel.__dict__['def_vals'][key] = userParams[key]

        # fit models using given fitMethods
        isoModelResult = None
        for fit_method in fitMethods:
            if isoModelResult:
                isoModelResult = isoModel.fit(
                    data=data[1], x=data[0],
                    params=isoModelResult.params,  # this is what's different
                    method=fit_method,
                    fit_kws=fitMethods[fit_method])
            else:
                isoModelResult = isoModel.fit(
                    data=data[1], x=data[0],
                    method=fit_method,
                    fit_kws=fitMethods[fit_method])
        self.isoModelResult = isoModelResult

        # validate fit model against isotherm theory
        self.modelValidity = self.ValidateFit()

    @abstractmethod
    def IsothermFunc():
        """ isotherm equation
        """
        raise NotImplementedError

    @abstractmethod
    def ValidateFit():
        """ validate the fitted parameters against specific
        isotherm theory
        """
        raise NotImplementedError


class Linear(Isotherm):

    __name__ = "Linear"

    @staticmethod
    def IsothermFunc(x, Kd=1):
        return Kd * x

    def ValidateFit(self):
        params_dict = self.isoModelResult.params.valuesdict()
        Kd = params_dict['Kd']
        # check model theory
        if Kd <= 0:
            return [False,
                    """best fit Kd = {0};
                    Linear adsorption theory requires Kd > 0"""
                    .format(Kd)]
        return [True, None]


class Freundlich(Isotherm):

    __name__ = "Freundlich"

    @staticmethod
    def IsothermFunc(x, Kf=1, n=1):
        return Kf * x ** (1 / n)

    def ValidateFit(self):
        params_dict = self.isoModelResult.params.valuesdict()
        Kf = params_dict['Kf']
        n = params_dict['n']
        if Kf <= 0:
            return [False,
                    """best fit Kf = {0};
                    Freundlich adsorption theory requires Kf > 0"""
                    .format(Kf)]
        if n < 1:
            return [False,
                    """best fit n = {0};
                    Freundlich adsorption theory requires n >= 1"""
                    .format(n)]
        return [True, None]


class Langmuir(Isotherm):

    __name__ = "Langmuir"

    @staticmethod
    def IsothermFunc(x, Qmax=1, Kl=1):
        return Qmax * Kl * x / (1 + Kl * x)

    def ValidateFit(self):
        params_dict = self.isoModelResult.params.valuesdict()
        Qmax = params_dict['Qmax']
        Kl = params_dict['Kl']
        if Qmax <= 0:
            return [False,
                    """best fit Qmax = {0};
                    Langmuir adsorption theory requires Qmax > 0"""
                    .format(Qmax)]
        if Kl <= 0:
            return [False,
                    """best fit Kl = {0};
                    Langmuir adsorption theory requires Kl > 0"""
                    .format(Kl)]
        return [True, None]
