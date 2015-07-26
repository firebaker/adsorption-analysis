import math


def _check_nan(item):
    for val in item:
        if math.isnan(val):
            raise TypeError


def _verify_numeric(item):
    if len(item) > 1:
        for val in item:
            if not isinstance(val, (int, long, float, complex)):
                raise TypeError
    else:
        if not isinstance(item, (int, long, float, complex)):
            raise TypeError
