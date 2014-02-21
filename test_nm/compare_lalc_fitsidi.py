from test_main import *

lalc     = LedaFits('2014-01-29-22h44m13_00.0.dada_0.LA')
uvw_lalc = lalc.d_uv_data["FLUX"]

dada     = LedaFits('2014-01-29-22h44m13_00.0.dada')
uvw_dada = lalc.d_uv_data["FLUX"]

print "\nChecking data array shapes match..."
print uvw_dada.shape
print uvw_lalc.shape
assert uvw_dada.shape == uvw_lalc.shape
print "OK"

print "\nChecking FLUX data match..."
print uvw_dada[0]
print uvw_lalc[0]
assert np.allclose(uvw_dada, uvw_lalc)
print "OK"