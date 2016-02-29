""" validation and check functions for various objects in
adsorption-analysis.py"""

import inspect


class NumericError(Exception):
    pass


class InputError(Exception):
    pass


# validate values are numeric
def validateNumeric(value):
    if not isinstance(value, (int, float, complex)):
        raise NumericError("value is non-numeric")


def validateData(data):
    data_error_message = "Data Error: "
    if len(data) != 2:
        raise InputError(
            data_error_message +
            'input data needs to be composed of exactly two arrays')
    if len(data[0]) != len(data[1]):
        raise InputError(
            data_error_message +
            'input data arrays are not of equal length')
    for array in data:
        array_pos = 0
        for value in array:
            value_pos = 0
            try:
                validateNumeric(value)
            except NumericError:
                raise InputError(
                    data_error_message +
                    "value @ data[{0}][{1}] is non-numeric"
                    .format(array_pos, value_pos))
            value_pos += 1
        array_pos += 1


def validateAlpha(alpha):
    alpha_error_message = "Alpha Error: "
    try:
        validateNumeric(alpha)
    except NumericError:
        raise InputError(
            alpha_error_message +
            'alpha must be a single numeric value')
    if alpha < 0 or alpha > 1:
        raise InputError(
            alpha_error_message +
            'alpha must be > 0 and < 1')


def validate_userParams(userParams, IsothermFunc):
    userParams_error_message = "UserParams Error: "
    if not isinstance(userParams, 'dict'):
        raise InputError(
            userParams_error_message +
            "userParams must be a dictionary type")
    param_names = inspect.getargspec(IsothermFunc)[0]
    del param_names[0]
    if not len(param_names) == len(userParams):
        raise InputError(
            userParams_error_message +
            'len(userParams) is not the expected length')
    # check that keys in the param names are in dict
    for param in param_names:
        if param not in set(userParams):
            raise InputError(
                userParams_error_message +
                'expected vars not in userParams')
