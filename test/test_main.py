import sys
sys.path.append("/Volumes/Storage/LEDA/interfits")

from interfits import *

def compare_dicts(dict_a, dict_b):
    """ Compare two dictionaries to confirm they contain same data.

    Second dict may have extra keys, but must have all keys of first dict.
    """
    all_ok = True
    ok_exceptions = ['STAXOF', 'POLCALA', 'POLCALB']

    for k in dict_a:
        if dict_b.has_key(k):
            try:
                if type(dict_a[k]) is float:
                    # Check that values are the same within tolerance for floats
                    assert np.allclose(dict_a[k], dict_b[k])
                elif type(dict_a[k]) is np.ndarray:
                    if type(dict_a[k][0]) is str:
                        assert all(dict_a[k] == dict_b[k])
                    else:
                        assert np.allclose(dict_a[k], dict_b[k])
                elif type(dict_a[k]) is np.core.defchararray.chararray:
                    assert all(dict_a[k] == dict_b[k])
                else:
                    assert dict_a[k] == dict_b[k]
            except:
                if k not in ok_exceptions:
                    print "Error:", k, dict_a[k], dict_b[k]
                    print type(dict_a[k])
                    all_ok = False
                else:
                    print "INFO: Known exception: %s"%k
        else:
            if k not in ok_exceptions:
                print "ERROR: %s not in both dictionaries"%k
                all_ok = False
            else:
                print "INFO: Known exception: %s"%k

    if all_ok:
        print "PASSED"
    else:
        print "ERROR"

    return all_ok