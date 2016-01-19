""" validation and check functions for various objects in
adsorption-analysis.py"""

# validate values are numeric
def validateNumeric(item):
    if hasattr(item, '__iter__'):
        for val in item:
            if not isinstance(val, (int, float, complex)):
                return False
    else:
        if not isinstance(item, (int, float, complex)):
            return False
    return True


def validateData(data):
    errorList = []
    for array in data:
        if not validateNumeric(array):
            errorList.append(
                'DataInputError: input values are not all numeric')
    if len(data) != 2:
        errorList.append(
            'DataInputError: data needs to have exactly two arrays')
    if len(data[0]) != len(data[1]):
        errorList.append(
            'DataInputError: input arrays are not of equal length')
    return errorList


def validateInput(data):
    errorList = validateData(data)
    return errorList


def validate_userParams(userParams, params):
    if len(userParams) != len(params):
        return """len(userParams) != len(params);
        isotherm parameters initialized with default values"""
    elif not validateNumeric(userParams):
        return """userParams values are not all numeric;
        isotherm parameters initialized with default values"""
    return False
