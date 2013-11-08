from test_main import *
import numpy as np

def test_json():    
    npDict_orig = {
        'ARRNAM'  : np.chararray([2]),
        'CHICKEN' : np.array([ii for ii in range(128)]),
        'STRING'  : 'qwertyuiop',
        'LIST'    : [1,7,2,9]
        }

    # Dump to file
    npDict_filename = 'data/jsontest.json'
    dump_json(npDict_orig, npDict_filename)

    # Load from file
    npDict_ff = load_json(npDict_filename)

    # Check the two dictionaries match
    ok_count = 0
    ok_count += compare_dicts(npDict_orig, npDict_ff)
    assert ok_count == 1
    

if __name__ == '__main__':

    test_json()