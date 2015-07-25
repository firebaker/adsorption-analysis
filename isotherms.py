# Linear isotherm equation
def linearIsotherm(x, Kd):
    return Kd * x


# Freundlich isotherm equation
def freundlichIsotherm(x, Kf, n):
    return Kf * x**(1/n)


# Langmuir isotherm equation and related functions
def langmuirIsotherm(x, Smax, Kl):
    return Smax * Kl * x / (1 + Kl * x)
