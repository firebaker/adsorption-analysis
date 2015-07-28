# adsorption-analysis
Adsorption-anlysis provides tools to analyze experimental adsorption results.

Current Features List:
- fit linear, freundlich, or langmuir isotherm to an experimental dataset using scipy.optimize.curve_fit; returns dictionary with popt, pcov, SSR (sum of squared residuals), or warning if unable to fit isotherm; initial parameters for fitting are set to "1" if not declared
- fit all isotherms to dataset; returns dictionary of each isotherm with popt, pcov, and SSR or warning if unable to fit isotherm; initial parameters for fitting are set to "1" if not declared
- (optional) guessQmax() -> use to guess initial Qmax for subsequent langmuir fitting; declare "xndpts" as needed, see function

TODO List:
- create function to plot isotherms
- add pertinent data to plots
- optionally, including confidence intervals to plots
- examine multiple adsorption datasets and compare isotherms
- check isotherms output for impossible 