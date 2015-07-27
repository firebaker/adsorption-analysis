""" furnctions to check data input and output"""


import math


# check if nan exists in item
def _check_nan(item):
    if hasattr(item, '__iter__'):
        for val in item:
            if math.isnan(val):
                raise TypeError
    else:
        if math.isnan(item):
            raise TypeError


# check if values are numeric
def _verify_numeric(item):
    if hasattr(item, '__iter__'):
        for val in item:
            if not isinstance(val, (int, long, float, complex)):
                raise TypeError
    else:
        if not isinstance(item, (int, long, float, complex)):
            raise TypeError


# check if arrays are of same length
def _verify_similar_length(x, y):
    if len(x) != len(y):
        raise Exception
