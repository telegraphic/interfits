from test_main import *
import pylab as plt
import time


def generate_lalc():
    """ Check that creating a fits-idi->interfits->fits-idi creates identical data. """

    lalc = InterFits('test_lalc.LA')
    lalc.export_fitsidi('test_lalc.fitsidi', 'test_lalc.xml')

    ok_count = 0
    assert ok_count == 0

def test_compare_lalc():

    uvf  = InterFits('test_lalc.fitsidi')
    lalc = InterFits('test_lalc.LA')


    lalc_flux = lalc.d_uv_data['FLUX']
    uvf_flux  = uvf.d_uv_data['FLUX']

    lalc_bls = lalc.d_uv_data['BASELINE']
    uvf_bls  = uvf.d_uv_data['BASELINE']

    #print lalc_bls
    #print uvf_bls

    assert np.allclose(lalc_bls, uvf_bls)
    print "PASS: BASELINE IDS MATCH"

    assert lalc_flux.shape == uvf_flux.shape
    print "PASS: FLUX SHAPE MATCH"

    print "Testing flux data..."
    for row in range(0, lalc_flux.shape[0]):
        if not row % 1000:
            print "%i of %i"%(row, lalc_flux.shape[0])
        try:
            xxl = lalc_flux[row][::8]**2 + lalc_flux[row][1::8]**2
            yyl = lalc_flux[row][2::8]**2 + lalc_flux[row][3::8]**2
            xyl = lalc_flux[row][4::8]**2 + lalc_flux[row][5::8]**2
            yxl = lalc_flux[row][6::8]**2 + lalc_flux[row][7::8]**2

            xxu = uvf_flux[row][::8]**2  + uvf_flux[row][1::8]**2
            yyu = uvf_flux[row][2::8]**2 + uvf_flux[row][3::8]**2
            xyu = uvf_flux[row][4::8]**2 + uvf_flux[row][5::8]**2
            yxu = uvf_flux[row][6::8]**2 + uvf_flux[row][7::8]**2

            assert np.allclose(xxl, xxu)
            assert np.allclose(yyl, yyu)
            assert np.allclose(xyl, xyu)
            assert np.allclose(yxl, yxu)

        except AssertionError:
            print "ERROR: Flux values do not agree"
            print uvf_flux[row, 0:10]
            print lalc_flux[row, 0:10]
            raise

    print "PASS: FLUX DATA MATCH"

if __name__ == '__main__':
    
    #generate_lalc()

    test_compare_lalc()