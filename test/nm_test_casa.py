fitsidifile = 'nm.fitsidi'
vis = 'nm.vis'
importfitsidi()

caltable = 'nm.B'
bandpass()

gaintable = 'nm.B'
applycal()

imagename = 'nm.img'
niter = 0
cell = '15.0arcmin'
clean()
