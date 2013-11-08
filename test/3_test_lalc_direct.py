from test_main import *
import pylab as plt
import time


def test_lalc():
    """ Check that LA->interfits->LA->interfits creates identical data. """

    lalc = LedaFits('data/test_lalc.LA')
    lalc.exportFitsidi('data/test_lalc_direct.fitsidi', '../config/config.xml')
    uvf  = LedaFits('data/test_lalc_direct.fitsidi')

    ok_count = 0
    # Check all header values
    h1("Testing header keywords")
    h2("Common")
    ok_count += compare_dicts(uvf.h_common, lalc.h_common)
    h2("Parameters")
    ok_count += compare_dicts(uvf.h_params, lalc.h_params)
    h2("Antenna")
    ok_count += compare_dicts(uvf.h_antenna, lalc.h_antenna)
    h2("Array Geometry")
    ok_count += compare_dicts(uvf.h_array_geometry, lalc.h_array_geometry)
    h2("Frequency")
    ok_count += compare_dicts(uvf.h_frequency, lalc.h_frequency)
    h2("Source")
    ok_count += compare_dicts(uvf.h_source, lalc.h_source)
    h2("UV DATA")
    ok_count += compare_dicts(uvf.h_uv_data, lalc.h_uv_data)

    h1("Testing data tables")
    h2("Antenna")
    ok_count += compare_dicts(uvf.d_antenna, lalc.d_antenna)
    h2("Array Geometry")
    ok_count += compare_dicts(uvf.d_array_geometry, lalc.d_array_geometry)
    h2("Frequency")
    ok_count += compare_dicts(uvf.d_frequency, lalc.d_frequency)
    h2("Source")
    ok_count += compare_dicts(uvf.d_source, lalc.d_source)

    assert ok_count == 11

    try:
        assert repr(lalc) == repr(uvf)
    except AssertionError:
        print "ERROR: __repr__ outputs do not match"
        print repr(uvf)
        print repr(lalc)

    os.remove('data/test_lalc_direct.fitsidi')

if __name__ == '__main__':
    
    test_lalc()