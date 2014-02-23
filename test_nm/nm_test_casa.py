# Delete old versions of files
_ip.system("rm -i -rf nm-cyg.img.*")
_ip.system("rm -i -rf nm-cyg.vis")
_ip.system("rm -i -rf nm-cyg.B")

# Import data
#fitsidifile = 'nm-cyg.fitsidi'
fitsidifile = 'Zenith/Zenith_b1_d20140129_utc224413.uvfits'
vis = 'nm-zen.vis'
importfitsidi()

#plotms()

# Apply bandpass correction
caltable = 'nm-cyg.B'
bandpass()
gaintable = 'nm-cyg.B'
applycal()

# Make dirty image
imagename = 'nm-cyg.img'
niter = 0
imsize = [512, 512]
cell = '15.0arcmin'
clean()