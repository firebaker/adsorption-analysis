# adsorption-analysis
Adsorption-anlysis provides tools to analyze experimental adsorption results

Current Features List:
- fit linear, freundlich, or langmuir isotherm to an experimental dataset using scipy.optimize.curve_fit; returns dictionary with popt, pcov, SSR (sum of squared residuals), and confidence interval (alpha dependent), or warning if unable to fit isotherm; initial parameters -> "1" if not declared, alpha -> 0.05 (conf intrvl = 95%)
- use adsorptionAnalysis() to fit linear, freundlich, and langmuir in one function 
- (optional) guessQmax() -> use to guess initial Qmax for subsequent langmuir fitting; declare "xndpts" as needed, see function
- (optional) guessK_() -> use to guess initial K for isotherm fitting, see function
- fit user-defined isotherm; see function

CAUTIONARY DISCLAIMER: Currently confidence intervals are produced through regression-covariance analysis. Any errors in the confidence intervals results in loss of that particular isotherm.  Outliers may be the cause of errrors, or poor fit. It is recommended to do outlier analysis before fitting.

TODO List:
- modifty the creating of confidence intervals from regression-covariance analysis to a bootstrapping analysis to help with outlier analysis
- create function to plot isotherms and confidence intervals
