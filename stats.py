""" statstical functions """


# calculate linear/non-linear regression SSR (sum of squared residuals)
# requires 'popt' variable from from curve_fit()
def _regSSR(func, x, y, popt):
    residuals = y - func(x, y, popt)
    SSR = sum(residuals**2)
    return SSR
