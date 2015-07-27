# adsorption-analysis
Adsorption-anlysis provides tools to analyze experimental adsorption results.

Current Features List:
- individually or all-at-once fit linear, freundlich, or langmuir isotherms to an experimental dataset using scipy.optimize.curve_fit; returns dictionary with popt and pcov or warning if unable to fit isotherm

TODO Feature List:
- statistically test for best fit isotherm using Sum of Squared Errors (SSE)
- create confidence interval for isotherms
- plot isotherms with pertinent data (optionally including confidence intervals)
- examine multiple adsorption datasets and compare isotherms
- check isotherms results for impossible outputs