"""isotherm.py"""

# Python standard library modules and functions
from abc import ABCMeta, abstractmethod

# Third party modules and functions
import numpy as np
import lmfit
from lmfit import minimize

# adsorption-analysis modules
import AAvalidate as val

#Error Handler Class
class Error_Handler:
    errorcount = 0
    def __init__(self, param, error_message):
        self.paramter = param
        self.error = true
        Error_Handler.errorcount += 1
        self.message = error_message
        self.description = "This class handles the errors in the lmfit routines. "
        
        # getNumberofErrors returns the number of errors collected.
        
    def getNumberofErrors(self):
        return Error_Handler.errorcount
    

class Isotherm(metaclass=ABCMeta):
    """An abstract base class of functions to fit isotherms to
    data and return fitted isotherm results.

    Arguments:
    data: A 2-array list of float64 (or convertable numeric),
          data[0] represents non-adsobed phase analyte concentration
          data[1] represents adsorbed phase analyte concentration
    userParams: Initial param values designated by the user;
                see individual isotherm classes for specific params

    Attributes:
    ***Note: some attributes are conditional upon "validateInput"
    data: A 2-array list of user data input; should be numeric;
          data[0] represents mobile phase analyte concentration;
          data[1] represents imobile phase analyte concentration
    error: Either False, a string, or a list of strings. If it evaluates
           to True, further calculations cease to occur.
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

    def __init__(self, data, userParams=None,
                 validateInput=True):
        # validate input
        if validateInput:
            self.data = data
            self.error = val.validateInput(self.data)
            if self.error:
                return None
            self.data = [np.float_(data[0]), np.float_(data[1])]
            data = self.data
        # initialize parameters
        self.params = self.InitParams()
        #
        #
        self.NumberofParameters = 0
        self.NumberofErrors = 0
        self.Error_List = []
        
        # validate/initialize user parameters
        self.userParams = userParams
        self.userWarning = "parameters initialized with default values"
        if self.userParams:
            self.userWarning = val.validate_userParams(
                self.userParams, self.params)
            if not self.userWarning:
                for userParam, param in zip(
                        self.userParams, self.params):
                    self.params[param].set(value=userParam)
        """ ***Note***
        Nelder-Mead minimiztion is generally considered more robust
        at fitting for special model cases than Levenburg-Marquardt
        Least Squares (useful for Freundlich isotherm fitting).
        However, where lmfit's leastsqr method provides uncertainties'
        measurements, lmfit's nelder method does not, therefore we use
        Nelder-Mead to obtain robust fits, we then use the subsequent
        Levenburg-Marquardt fit to obtain the uncertainty values of
        the Nelder-Mead fit parameters.
        """
        # fit model using Nelder-Mead minimization
        self.nelder = minimize(
            self.MinimizeFunc,
            self.params,
            kws={'x': data[0], 'y': data[1]},
            method='nelder')
        # fit model using Levenburg-Marquart minimization with
        # Nelder-Mead minimization fit parameters
        self.leastsqr = lmfit.minimize(
            self.MinimizeFunc,
            self.nelder.params,
            kws={'x': data[0], 'y': data[1]},
            method='leastsqr')
        self.minimizedFit = self.leastsqr
        # validate fit model against isotherm theory
        self.modelValidity, self.modelValidtyMsg = \
            self.ValidateFit(self.leastsqr.params)

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


class Linear(Isotherm):

    __name__ = "Linear"

    @staticmethod
    def IsothermFunc(x, Kd):
        return Kd * x

    @staticmethod
    def MinimizeFunc(params, x, y):
        paramvals = params.valuesdict()
        Kd = paramvals['Kd']
        NumberofParameters = 1
        model = Linear.IsothermFunc(x=x, Kd=Kd)
        return model - y

    @staticmethod
    def InitParams():
        params = lmfit.Parameters()
        params.add('Kd', value=1.0)
        return params

    @staticmethod    
    def ValidateFit(params):
        Error_List=[]
        NumberofErrors = 0
        paramvals = params.valuesdict()
        Kd = paramvals['Kd']
        
        if Kd <= 0:
           Error_List.add(Error_Handler(1, """ best fit Kd = {0};
                Linear adsorption theory requires Kd > 0
                """.format(Kd)))
           NumberofErrors +=1
       


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
        NumberofParameters=2
        params.add('Kf', value=1)
        params.add('n', value=1)
        return params

    @staticmethod
    def ValidateFit(params):
        # clear error list
        ErrorList=[]
        NumberofErrors=0
        message = None
        message_list = []
        paramvals = params.valuesdict()
        Kf = paramvals['Kf']
        n = paramvals['n']
        if Kf <= 0:
           
             Error_List.add(Error_Handler(1,""" best fit Kf = {0};
            Freundlich adsorption theory requires Kf > 0
            """.format(Kf)))
            NumberofErrors +=1
        if n < 1:
            Error_List.add(Error_Handler(2,""" best fit n = {0};
            Freundlich adsorption theory requires n >= 1
            """.format(n)))
            NumberofErrors +=1
        


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
        NumberofParameters = 2
        params.add('Qmax', value=1.0)
        params.add('Kl', value=1.0)
        return params

    @staticmethod
    def ValidateFit(params):
       # clear error list
        ErrorList=[]
        NumberofErrors=0
        paramvals = params.valuesdict()
        Qmax = paramvals['Qmax']
        Kl = paramvals['Kl']
        if Qmax <= 0:
            Error_List.add(Error_Handler(1,""" best fit Qmax = {0};
            Langmuir adsorption theory requires Qmax > 0
            """.format(Qmax)))
            NumberofErrors +=1
        if Kl <= 0:
            Error_List.add(Error_Handler(2,""" best fit Kl = {0};
            Langmuir adsorption theory requires Kl > 0
            """.format(Kl)))
            NumberofErrors +=1
       
       
