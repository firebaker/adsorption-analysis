# adsorption-analysis
adsorption-anlysis provides tools to analyze experimental adsorption results.

11/16/15 - Things are working with outside module calls in an OOP fashion. 

11/13/15 - We're switching from a script based approach to an OOP based approach. Things have become buggy, and the feature list is out of date.

Current Features List:
- fit linear, freundlich, and/or langmuir isotherm to an experimental dataset using scipy.optimize.curve_fit;

- adsorption.AdsorptionAnlysis(data, alpha) returns a sample object of the currently implemented best fit isotherms. (linear, freundlich, langmuir)
 isotherm attributes include:
 isotherm.pop0 - initial parameters for scipy.curvefit()
 isotherm.popt - output parameters from scipy.curvefit()
 isotherm.pcov - output covariance matrix from scipy.curvefit()
 isotherm.conf_upper - upper confidence interval parameters
 isotherm.conf_lower - lower confidence interval parameters
 # isotherm.pred_upper (not implemented yet)
 # isotherm.pred_lower (not implemented yet)
 isotherm.SSR - Sum of Squared Residuals
 isotherm.AIC - Akiake Information Criterion
 isotherm.BIC - Baysian Information Criterion

- isotherm.py contains the base class from which to build isotherm calculation functions, incase the currently included isotherms are not suffience for you, i.e. you can build your own (see code)