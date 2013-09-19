from test_main import *

import pprint
pp = pprint.PrettyPrinter(indent=4)
pp = pp.pprint

def test_dada():

    uvf = InterFits('test_lalc.fitsidi')
    ddd = InterFits('test.dada')

    h2("HEADER: COMMON")
    pp(ddd.h_common)

    h2("HEADER: PARAMS")
    pp( ddd.h_params)

    h2("HEADER: UV DATA")
    pp( ddd.h_uv_data)

    h2("HEADER: ARRAY GEOMETRY")
    pp( ddd.h_array_geometry)

    h2("HEADER: ANTENNA")
    pp( ddd.h_antenna)

    h2("HEADER: FREQUENCY")
    pp( ddd.h_frequency)

    h2("HEADER: SOURCE")
    pp( ddd.h_frequency)

    h2("DATA: FREQUENCY")
    pp(ddd.d_frequency)


if __name__ == '__main__':

    test_dada()