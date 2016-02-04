"""isotherm.py"""

# Python standard library modules and functions
from abc import ABCMeta, abstractmethod

# Third party modules and functions
import numpy as np
import lmfit

# adsorption-analysis modules
import AAvalidate as val
import AAstats


class Isotherm(metaclass=ABCMeta):
    """An abstract base class of isotherm fitting functions.

    Arguments:
    data: A 2-array list of float64 (or convertable numeric),
        * data[0] represents non-adsobed phase analyte concentration
        * data[1] represents adsorbed phase analyte concentration
    alpha: # TODO
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
    modelValidtyMsg: A message relaying the cause of model validity failure
    """

    def __init__(self, data, alpha=0.05, userParams=None,
                 fitMethod=['nelder', 'leastsqr'], validateInput=True, 
                 confMethod='default'):

        # validate input
        # TODO check alpha
        # TODO check fitMethod
        if validateInput:
            self.error = val.validateInput(data)
            #self.error = val.validateInput(data, alpha, fitMethod)
            if self.error:
                return None

        # initialize parameters
        self.params = self.InitParams()
        self.userParams = userParams
        self.userWarning = "parameters initialized with default values"
        if self.userParams:
            self.userWarning = val.validate_userParams(
                self.userParams, self.params)
            if not self.userWarning:
                for userParam, param in zip(
                        self.userParams, self.params):
                    self.params[param].set(value=userParam)

        # check if fitMethod is composed of multiple methods
        # and fit isotherm
        if hasattr(fitMethod, '__iter__'):
            workingfit = None
            for fit in fitMethod:
                if workingfit:
                    params = workingfit.params
                else:
                    params = self.params
                workingfit = lmfit.minimize(self.MinimizeFunc, params,
                                            kws={'x': data[0], 'y': data[1]},
                                            method=fit)
            self.minimizedFit = workingfit
        else:
            self.minimizedFit = lmfit.minimize(self.MinimizeFunc, self.params,
                                            kws={'x': data[0], 'y': data[1]},
                                            method=fitMethod)

        # validate fit model against isotherm theory
#        self.modelValidity, self.modelValidtyMsg = \
#            self.ValidateFit(self.minimizedFit.params)
        self.modelValidity, self.modelValidtyMsg, self.confInterval = \
            self.ValidateFit(data[0], alpha, confMethod)

    """All the methods must be made static in inherited instance
    to allow for users to call these functions in their own
    applications"""

    @abstractmethod
    def IsothermFunc():
        """ isotherm equation
        """
        raise NotImplementedError

    @abstractmethod
    def InitParams():
        """ initialize isotherm parameters
        """
        raise NotImplementedError

    @abstractmethod
    def MinimizeFunc():
        """ isotherm equation
        """
        raise NotImplementedError

    @abstractmethod
    def ValidateFit():
        """ validate the fitted parameters against specific
        isotherm theory
        """
        raise NotImplementedError

    _Asym = ['Asymptotic', 'asymptotic', 'Asym', 'asym']
    @staticmethod
    def IsoConfAsym(xdata, popt, pcov, alpha=0.05):
        """Default asymptotic confidence interval calculation;
        This may need to be modified for specific isotherms,
        e.g. Freundlich Isotherm.
        """
        return AAstats.RegConfAsym(xdata, popt, pcov, alpha)

    
    _Monte = ['Monte Carlo', 'monte carlo', 'Monte', 'monte']
    @staticmethod
    def IsoConfMonte(xdata, popt, pcov, func=NotImplementedError, alpha=0.05,
                     seed=123456789, sim=1000):
        """ Default Monte Carlo confidence interval calculation;
        This may need to be modified for specific isotherms,
        e.g. Freundlich Isotherm.
        """
        return AAstats.RegConfMonte(func=func, xdata=xdata, popt=popt, pcov=pcov, alpha=alpha,
                    seed=seed, sim=sim)


class Linear(Isotherm):

    __name__ = "Linear"

    @staticmethod
    def IsothermFunc(x, Kd):
        return Kd * x

    @staticmethod
    def MinimizeFunc(params, x, y):
        paramvals = params.valuesdict()
        Kd = paramvals['Kd']
        model = Linear.IsothermFunc(x=x, Kd=Kd)
        return model - y

    @staticmethod
    def InitParams():
        params = lmfit.Parameters()
        params.add('Kd', value=1.0)
        return params
    
#    @staticmethod
#    def ValidateFit(params):
#        paramvals = params.valuesdict()
#        Kd = paramvals['Kd']
#        if Kd <= 0:
#            return False, """ best fit Kd = {0};
#                Linear adsorption theory requires Kd > 0
#                """.format(Kd)
#        return True, None

    def ValidateFit(self, xdata, alpha=0.05, confMethod='default'):
        paramvals = self.minimizedFit.params.valuesdict()
        Kd = paramvals['Kd']
        # check model theory
        if Kd <= 0:
            return False, """best fit Kd = {0};
                Linear adsorption theory requires Kd > 0
                """.format(Kd), None
        # check Confidence Interval theory
        popt = [Kd]
        if confMethod in ['default', *self._Asym]:
            confIntrvl = Linear.IsoConfAsym(xdata, popt, self.minimizedFit.covar, alpha)
        elif confMethod in [*self._Monte]:
            confIntrvl = Linear.IsoConfMonte(xdata=xdata, popt=popt, pcov=self.minimizedFit.covar, func=Linear.IsothermFunc)
        else:
            return False, """confMethod: {0}, is not an option""", None
        if confIntrvl['lower'][0] < 0:
            return False, """lower conf interval Kd = {0};
                For the given alpha ({1}), the best fit Linear Isotherm
                is not statistically different from zero;
                Linear adsorption theory requires the best fit Isotherm
                to be statistically greater than 0.
                """.format(confIntrvl['lower'][0], alpha), confIntrvl
        return True, None, confIntrvl


class Freundlich(Isotherm):

    __name__ = "Freundlich"

    @staticmethod
    def IsothermFunc(x, Kf, n):
        return Kf * x ** (1 / n)

    @staticmethod
    def MinimizeFunc(params, x, y):
        paramvals = params.valuesdict()
        Kf = paramvals['Kf']
        n = paramvals['n']
        model = Freundlich.IsothermFunc(x=x, Kf=Kf, n=n)
        return model - y

    # n must be > 0
    @staticmethod
    def InitParams():
        params = lmfit.Parameters()
        params.add('Kf', value=1)
        params.add('n', value=1)
        return params

#    @staticmethod
#    def ValidateFit(params, cov, alpha=0.05):
#        validModel = True
#        message = None
#        message_list = []
#        paramvals = params.valuesdict()
#        Kf = paramvals['Kf']
#        n = paramvals['n']
#        if Kf <= 0:
#            validModel = False
#            message_list.append(""" best fit Kf = {0};
#            Freundlich adsorption theory requires Kf > 0
#            """.format(Kf))
#        if n < 1:
#            validModel = False
#            message_list.append(""" best fit n = {0};
#            Freundlich adsorption theory requires n >= 1
#            """.format(n))
#        if message_list:
#            if len(message_list) > 1:
#                message = message_list
#            else:
#                message = message_list[0]
#        return validModel, message

    def ValidateFit(self, xdata, alpha=0.05, confMethod='default'):
        validModel = True
        message = None
        message_list = []
        paramvals = self.minimizedFit.params.valuesdict()
        Kf = paramvals['Kf']
        n = paramvals['n']
        confIntrvl = None
        if Kf <= 0:
            validModel = False
            message_list.append("""best fit Kf = {0};
            Freundlich adsorption theory requires Kf > 0
            """.format(Kf))
        if n < 1:
            validModel = False
            message_list.append("""best fit n = {0};
            Freundlich adsorption theory requires n >= 1
            """.format(n))
        # check Confidence Interval theory
        if validModel:
            popt = [Kf, n]
            confIntrvl = Freundlich.IsoConfAsym(xdata, popt,
                                                self.minimizedFit.covar,
                                                alpha)
            if confIntrvl['upper'][1] < 1:
                validModel = False
                message_list.append("""upper confidence interval n = {0};
                Freundlich adsorption theory requires n >= 1
                """.format(confIntrvl['upper'][1]))
            if confIntrvl['lower'][0] <= 0:
                validModel = False
                message_list.append("""lower confidence interval Kf = {0};
                Freundlich adsorption theory requires Kf > 0
                """.format(confIntrvl['lower'][0]))
        if message_list:
            if len(message_list) > 1:
                message = message_list
            else:
                message = message_list[0]
        return validModel, message, confIntrvl
            
    @staticmethod
    def IsoConfAsym(xdata, popt, pcov, alpha=0.05):
        """Freundlich asymptotic confidence interval calculation;
        Modified from default"""
        confIntrvl = AAstats.RegConfAsym(xdata, popt, pcov, alpha)
        confIntrvl['upper'][1], confIntrvl['lower'][1] =\
            confIntrvl['lower'][1], confIntrvl['upper'][1]
        return confIntrvl


class Langmuir(Isotherm):

    __name__ = "Langmuir"

    @staticmethod
    def IsothermFunc(x, Qmax, Kl):
        return Qmax * Kl * x / (1 + Kl * x)

    @staticmethod
    def MinimizeFunc(params, x, y):
        paramvals = params.valuesdict()
        Qmax = paramvals['Qmax']
        Kl = paramvals['Kl']
        model = Langmuir.IsothermFunc(x=x, Qmax=Qmax, Kl=Kl)
        return model - y

    @staticmethod
    def InitParams():
        params = lmfit.Parameters()
        params.add('Qmax', value=1.0)
        params.add('Kl', value=1.0)
        return params

#    @staticmethod
#    def ValidateFit(params, cov, alpha=0.05):
#        validModel = True
#        message = None
#        message_list = []
#        paramvals = params.valuesdict()
#        Qmax = paramvals['Qmax']
#        Kl = paramvals['Kl']
#        if Qmax <= 0:
#            validModel = False
#            message_list.append(""" best fit Qmax = {0};
#            Langmuir adsorption theory requires Qmax > 0
#            """.format(Qmax))
#        if Kl <= 0:
#            validModel = False
#            message_list.append(""" best fit Kl = {0};
#            Langmuir adsorption theory requires Kl > 0
#            """.format(Kl))
#        if message_list:
#            if len(message_list) > 1:
#                message = message_list
#            else:
#                message = message_list[0]
#        return validModel, message

    def ValidateFit(self, xdata, alpha=0.05, confMethod='default'):
        validModel = True
        message = None
        message_list = []
        paramvals = self.minimizedFit.params.valuesdict()
        Qmax = paramvals['Qmax']
        Kl = paramvals['Kl']
        confIntrvl = None
        if Qmax <= 0:
            validModel = False
            message_list.append(""" best fit Qmax = {0};
            Langmuir adsorption theory requires Qmax > 0
            """.format(Qmax))
        if Kl <= 0:
            validModel = False
            message_list.append(""" best fit Kl = {0};
            Langmuir adsorption theory requires Kl > 0
            """.format(Kl))
        # check Confidence Interval theory
        if validModel:
            popt = [Qmax, Kl]
            confIntrvl = Langmuir.IsoConfAsym(xdata, popt,
                                              self.minimizedFit.covar,
                                              alpha=0.05)
            if confIntrvl['lower'][0] <= 0:
                validModel = False
                message_list.append("""upper confidence interval Qmax = {0};
                Langmuir adsorption theory requires Qmax > 0
                """.format(confIntrvl['lower'][0]))
            if confIntrvl['lower'][1] <= 0:
                validModel = False
                message_list.append("""lower confidence interval Kl = {0};
                Langmuir adsorption theory requires Kl > 0
                """.format(confIntrvl['lower'][1]))
        if message_list:
            if len(message_list) > 1:
                message = message_list
            else:
                message = message_list[0]
        return validModel, message, confIntrvl
