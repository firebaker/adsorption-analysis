# adsorption-analysis
Adsorption-anlysis provides tools to analyze experimental adsorption results.

Current Features List:
- fit linear, freundlich, or langmuir isotherm to an experimental dataset using scipy.optimize.curve_fit; returns dictionary with popt, pcov, SSR (sum of squared residuals), and confidence interval (alpha dependent), or warning if unable to fit isotherm; initial parameters for fitting are set to "1" if not declared, alpha = 0.05 (conf intrvl = 95%)
- fit linear, freundlich, and langmuir isotherms at once using adsorptionAnalysis()
- (optional) guessQmax() -> use to guess initial Qmax for subsequent langmuir fitting; declare "xndpts" as needed, see function
- fit user-defined isotherm; see function

TODO List:
- create function to plot isotherms
- add pertinent data to plots
- optionally, including confidence intervals to plots
- examine multiple adsorption datasets and compare isotherms
- check isotherms output for impossible 