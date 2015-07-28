""" functions to aid in building data structures"""


# check if key-of-dict exists in dict, if not create
def _chk_key_dict(keydict, dict):
    if hasattr(keydict, '__iter__'):
        for key in keydict:
            if not (key in set(dict)):
                dict.update({key: {}})
    else:
        if not (keydict in set(dict)):
            dict.update({keydict: {}})


# check if key-of-list exists in dict, if not create
def _chk_key_list(keylist, dict):
    if hasattr(keylist, '__iter__'):
        for key in keylist:
            if not (key in set(dict)):
                dict.update({key: []})
    else:
        if not (keylist in set(dict)):
            dict.update({keylist: []})
