""" background functions strictly related to Adsorption_Analysis.py """


from checkPass import _verify_numeric
from checkPass import _verify_similar_length


""" bare isotherm equations """


# Linear isotherm equation
def linearIsotherm(x, Kd):
    return Kd * x


# Freundlich isotherm equation
def freundlichIsotherm(x, Kf, n):
    return Kf * x**(1/n)


# Langmuir isotherm equation
def langmuirIsotherm(x, Qmax, Kl):
    return Qmax * Kl * x / (1 + Kl * x)


""" checkPass """


# check adsorption input
def checkAdsorptionInput(x, y):
    # set variables for input
    errorcheck = 0
    errorlist = []
    # verify values are numeric
    for array in (x, y):
        try:
            _verify_numeric(array)
        except TypeError:
            errorcheck = errorcheck + 1
            errorlist.append(
                'TypeError: not all input values are numeric')
    # verify arrays of equal length
    try:
        _verify_similar_length(x, y)
    except:
        errorcheck = errorcheck + 1
        errorlist.append(
            'ArrayError: input arrays are not of equal length')
    # stop function if errorcheck > 0
    if errorcheck > 0:
        print 'Error in -> adsorptionAnalysis:'
        for error in errorlist:
            print error
        return None
