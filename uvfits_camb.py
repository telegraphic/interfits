# Bojan Nikolic <bn204@mrao.cam.ac.uk>, <bojan@bnikolic.co.uk>
#
# Initial version October 2007
#
# Some useful routines for handling UV-FITS data

from itertools import izip

import numpy.numarray as numarray
import numpy.numarray.linear_algebra as la
import numpy
import re

#from setup import *

import pyfits
import iofits4

import interptools



def AutoFITSOpen(fn):

    "If passed in string assume file to be opened"

    def _autoopen(x, *args, **kwargs):
        if type(x)== str:
            fin=pyfits.open(x)
            return fn(fin, *args, **kwargs)
        else:
            return fn(x, *args, **kwargs)
    _autoopen.__doc__=fn.__doc__
    _autoopen.func_name=fn.func_name+"wAO"
    return _autoopen

def DefMask(fn):

    "If a parameter called mask is passed and is ==None, fill default"

    def _defmask(x, *args, **kwargs):
        if "mask" in kwargs and kwargs["mask"]==None:
            if type(x) == pyfits.HDUList:
                l = len(x[0].data)
            elif type(x) == str :
                l = len(pyfits.open(x)[0].data)
            else:
                l = len(x)
            tmask=numpy.ones( l , numpy.bool)
            kwargs["mask"]=tmask
        return fn(x, *args, **kwargs)

    _defmask.__doc__=fn.__doc__
    return _defmask

def ParamDict(hin):

    "Return the dictionary of the random parameters"

    """
    The keys of the dictionary are the parameter names uppercased for
    consistency. The values are the column numbers.

    If multiple parameters have the same name (e.g., DATE) their
    columns are entered as a list.
    """

    pre=re.compile(r"PTYPE(?P<i>\d+)")
    res={}
    for k,v in hin.header.items():
        m=pre.match(k)
        if m :
            vu=v.upper()
            if vu in res:
                res[ vu ] = [ res[vu], int(m.group("i")) ]
            else:
                res[ vu ] = int(m.group("i"))
    return res
            

@AutoFITSOpen
def BaselineCoords(fin):

    "Return uu, vv , ww coordinates"

    d=ParamDict(fin[0])

    if "UU" in d:
        cu,cv,cw = d["UU"], d["VV"], d["WW"]
    elif "UU---SIN" in d:
        cu,cv,cw = d["UU---SIN"], d["VV---SIN"], d["WW---SIN"]
    else:
        raise "File %s appears not to have baseline coordinates" % str(fin)

    uu, vv, ww = [ numpy.array(fin[0].data.field("C%i" %x)) for x in [cu,cv,cw] ]
    
    return uu, vv, ww
    

@AutoFITSOpen
def BaselineAntennas(fin):

    "Return the antennas on each baseline"

    d=ParamDict(fin[0])

    cb=d["BASELINE"]

    return BaselineFileToBaseLine( fin[0].data.field("c%i" % cb) )
    

@AutoFITSOpen
def NVis(fin):

    "Return the number of visibilities in the file"

    return len(fin[0].data)

    
def BaselineFileToBaseLine(bf):
    "Turn the compressed baseline format to antena pair"

    """
    The field BASELINE is the baseline number (256ant1 + ant2 +
    subarray/100.)  See
    http://www.aoc.nrao.edu/aips/CookHTML/CookBookse142.html
    """

    ant2= numpy.mod(bf, 256)
    ant1= (bf - ant2)/256

    return ( numpy.array(ant1, numpy.int) ,
             numpy.array(ant2, numpy.int) )

@DefMask
def UVDataToBaselineMatrix(fin,
                           mask=None):

    """
    Transform data from UV-fits data to a matrix indexed on baseline
    pair, in which each element contains the data element of the
    uvfits file. 
    """

    din=fin[0].data
    mask=numarray.array(mask)
    a1v, a2v =BaselineFileToBaseLine(din.field("c6")[mask])

    msize= a2v.max() +1

    res= numpy.zeros( (msize,msize) + din.field("data")[0].shape,
                      numpy.float32)

    for i, row in enumerate(din.field("data")[mask]):

        res[a1v[i], a2v[i] ] = row
        

    return res

def UVFilePhase(fnamein):

    return UVDataPhase(pyfits.open(fnamein)[0].data)

def UVDataPhase(din,
                mask=None):

    "Return phase of each measurement"


    if mask is None:
        c=din.field("data")
        ph = numpy.arctan2( c[:,:,:,:,:,:,1] , c[:,:,:,:,:,:,0] )
    else:
        c=din.field("data")[numarray.array(mask)]
        ph = numpy.arctan2( c[:,:,:,:,:,:,1] , c[:,:,:,:,:,:,0] )

    return ph

def UVDataAmp(din):

    "Return phase of each measurement"

    c=din.field("data")

    a = numpy.sqrt( c[:,:,:,:,:,:,1]**2 + c[:,:,:,:,:,:,0]**2 )

    return a
    

def UVDataKAbs(din):

    "Return absolute value of the  baseline"

    """Note that units of u,v,w are conventionally seconds
    see
    ftp://ftp.aoc.nrao.edu/pub/software/aips/TEXT/PUBL/FITS-IDI.pdf
    """
    

    uu, vv, ww = [din.field("c%i" % x) for x in [1,2,3]]

    return numpy.sqrt( uu**2 + vv**2 + ww**2 )


@AutoFITSOpen
def AntennaTable(fin):

    antennahduname="AIPS AN"
    try:
        ahdu  = fin.index_of(antennahduname)
    except KeyError, e:
        ahdu=2

    return fin[ahdu]
    
    
def AntennaPositions(fnamein):

    """
    Note that as per the FITS-IDI standard, antennas are stored in
    geocentric (centre-of-mass) coordinates. Also note returned are
    offsets from the array centre.
    """

    return pyfits.open(fnamein)[2].data.field("STABXYZ")


def AntennaPosProjectF(fnamein,
                       **kwargs):

    return AntennaPosProject(pyfits.open(fnamein),
                             **kwargs)
    
def AntennaPosProject(fin,
                      localmeridian=False):

    "Positions of antennas in the array projected on a tangent plane to array centre"

    """

    localmeridian: FITS-IDI defines the coordinate sytem as: z from
    centre of earth to north pole. x from centre to intersection of
    grenwich meridian and equator. If localmeridian then assume x goes
    from centre of eart to centre of array.
    
    """

    hin= AntennaTable(fin)

    X,Y,Z=  [hin.header["ARRAY"+x] for x in ["X", "Y", "Z"]]
    if localmeridian:
        X= (X**2 + Y**2)**0.5
        Y= 0
    dx,dy,dz=  [ hin.data.field("STABXYZ")[:,i] for i in range(3)]

    return AntennaPosProjCalc( X, Y, Z,
                               dx, dy, dz )


def AntennaPosProjCalc( X, Y, Z,
                        dx, dy, dz ):

    "Does the actual calculation"
    
    rho = (X**2+Y**2)**0.5
    r   = (X**2+Y**2+Z**2)**0.5

    drho = ( X*dx+Y*dy) / rho

    xprime = (X*dy - Y*dx) / rho
    yprime = (rho * dz - Z* drho ) / r

    res  = numpy.array( (xprime, yprime) )
    
    return  res.transpose()

def AntennaNames(fin):

    "Return dictionary conncecting antenna names with seequency numbers"

    hin= AntennaTable(fin)

    res= {}

    for i,row in enumerate(hin.data):
        res[row.field("ANNAME") ]  = i +1

    return res
    
    
    

    

def BaselinePosList(fin   ):

    "Return the positions of the two antennas involved in each visibility measurement"

    """
    fin is the opened pyfits object
    """
    
    apos=       AntennaPosProject(fin)
    a1l, a2l  = BaselineFileToBaseLine(fin[0].data.field("c6"))

    # a1l, a2l are 1-based, convert to zero based :
    a1l -= 1
    a2l -= 1

    a1pos = numpy.take(apos, a1l, axis=0 )
    a2pos = numpy.take(apos, a2l, axis=0 )

    return a1pos, a2pos


def IntegrationDayTime(fin):

    "Return the days and times of integrations in file"

    dc=ParamDict(fin[0])["DATE"]

    return ( fin[0].data.field("c%i" % dc[0]) ,
             fin[0].data.field("c%i" % dc[1]) )

def IntegrationTimeSeries(fin):

    "Return time stamp for each data point, starting w zero"

    """
    Return vector in units of seconds
    """

    dv, tv = IntegrationDayTime(fin)
    ts =  ( (dv-dv[0])  + tv-tv[0]  ) * 24 * 3600
    
    return ts
    

@DefMask
def IntegrationList(fin,
                    mask=None):

    "Return the time stapms of all integrations in a file"

    dayv = fin[0].data.field("c4")[numarray.array(mask)]
    timev= fin[0].data.field("c5")[numarray.array(mask)]

    daylist = numpy.unique( dayv )
    if len(daylist) > 1 :
        raise "Don't know how to handle obs split across days yet"

    timelist = numpy.unique( timev)
    daylist  = numpy.take(daylist ,  [0 for x in timelist]  )

    return daylist, timelist

@AutoFITSOpen    
def SourcesColumn(fin):

    "Return the column containing the sources id for the uv data"
    hin=fin[0]
    pd =ParamDict(hin)

    return hin.data.field("c%i" % pd["SOURCE"])

@AutoFITSOpen    
def SourceVisMask(fin, src):
    "Return mask for visibilities corresponding to source src"

    return (SourcesColumn(fin) == src)

def SourcesList(fin):

    "Return a list of source id's present in uv set"

    sc=SourcesColumn(fin)
    return numpy.unique(sc)


def StokesCol(fin):

    "Findout which dimension is stokes"

    sre=re.compile(r"CTYPE(?P<i>\d+)")
    
    for k in fin[0].header.items():
        if k[1] == "STOKES":
            stokescol= int(sre.match(k[0]).group("i"))

    return stokescol


@AutoFITSOpen    
def StokesDictionary(fin):

    "Find out what stokes parameters sit in which column"

    ntos = { 1 : "I",
             2 : "Q",
             3 : "U",
             4 : "V",
             -1 : "RR",             
             -2 : "LL",
             -3 : "RL",
             -4 : "LR",
             -5 : "XX",
             -6 : "YY",
             -7 : "XY",
             -8 : "YX" }

    stokescol = StokesCol(fin)

    delt=int( fin[0].header["CDELT%i"%stokescol] )
    if delt != -1 :
        raise "Don't know how to handle stokes specifications"
    refval = int( fin[0].header["CRVAL%i"%stokescol] )

    res = {}

    for i in range( int(fin[0].header["NAXIS%i"%stokescol])):
        res[ ntos[refval + i *delt] ] = i

    return res

def DePhaseDataAcross(fnamein,
                      fnameout,
                      phase):

    "Like DePhaseData but assume same phase across all IFs/winds/etc"

    fin=pyfits.open(fnamein)
    for l in fin[0].data.field("data").shape[1:-1]:
        phase= numpy.array( [phase for x in range(l) ] )
    phase=phase.transpose()

    DePhaseData(fnamein,
                fnameout,
                phase)
    
def DePhaseData(fnamein,
                fnameout,
                phase):

    "Dephase data"

    """
    Rotates the phase of visibilities by the amount specified in the
    phase vector.
    """

    f= pyfits.open(fnamein)
    c=f[0].data.field("data")
    
    dr = c[:,:,:,:,:,:,0]
    di = c[:,:,:,:,:,:,1]

    cosph = numpy.cos(phase)
    sinph = numpy.sin(phase)

    res_real = dr * cosph - di *sinph
    res_imag = di * cosph + dr *sinph

    
    c[:,:,:,:,:,:,1] = res_imag
    c[:,:,:,:,:,:,0] = res_real

    iofits4.Write( f,
                   fnameout,
                   overwrite=1)    

    
    
def GaussianDephase(fnamein,
                    fnameout,
                    rms=0.1):

    "Simple procedure to add gaussian noise to phases"


    f= pyfits.open(fnamein)
    DePhaseData(fnamein,
                fnameout,
                numpy.random.normal(scale=rms,
                                    size=f[0].data.field("data").shape[:-1] ) )

class CalibrationData:

    "A class to contain data for calibration"

    """
    ... to avoid unnecessary reads from disk
    """

    def __init__(self,
                 fnamein):
        self.fin=pyfits.open(fnamein)

        self.a1v, self.a2v = BaselineFileToBaseLine(self.fin[0].data.field("c6"))
        
        # rebase to 0
        self.a1v -= 1; self.a2v -=1;

        # Visibility phases
        self.visphase= UVDataPhase(self.fin[0].data)        


def PhaseMatrix(din,
                mask=None):

    "Returns the matrix representing the antenna phase differences formed by correlator"

    """
    Note that this representation is underconstrained when solving for
    phases since the rank of the matrix is one less then number of
    antennas (only measuring differences!!)
    """


    if mask ==None:
        mask = numarray.ones( len(din) , numarray.Bool)
    a1v, a2v =BaselineFileToBaseLine(din.field("c6")[numarray.array(mask)])
    # rebase to 0
    a1v -= 1; a2v -=1;
    amax=max( max(a1v), max(a2v))+1
    res=numarray.zeros( (len(a1v), amax ))

    for i, (a1,a2) in enumerate( izip(a1v, a2v)) :
        res[i, a1] =1
        res[i, a2] =-1
        

    return res

def PhaseMatrixZeroMean(din,
                        mask):

    "As Phase Matrix Assume that the mean of antennas is zero"

    """
    This additional constraint allows us to fully solve for the
    phases. Below simply assume that the phase of the last atntenna is
    sum of phases of all of the other atnennas.
    """

    if isinstance(din, CalibrationData):
        a1v, a2v = din.a1v[mask], din.a2v[mask]
    else:
        a1v, a2v =BaselineFileToBaseLine(din.field("c6")[numarray.array(mask)])
        # rebase to 0
        a1v -= 1; a2v -=1;

    amax=max( max(a1v), max(a2v))+1    
    res=numarray.zeros( (len(a1v), amax-1 ))

    for i, (a1,a2) in enumerate( izip(a1v, a2v)) :
        
        if a1 == (amax-1):
            res[i,:] -= 1
        else:
            res[i, a1] =1
        if a2 == (amax-1):
            res[i,:]   += 1
        else:
            res[i, a2] =-1
        

    return res
    
    

def SolveForAntennaPhases(din,
                          mask=None,
                          ref="ZeroMean",
                          pout=None,
                          pol="mean"):

    "Solve for antenna phases, erm,..."

    """
    pout allows pickled output of the arrays to check for solver
    instability

    pol: how to deal with polarisations. "mean" -> take mean of all polarisations. integer: consider only this
    polarisation
    """

    if mask ==None:
        mask = numpy.ones( len(din) , numpy.bool)

    if ref== "ZeroMean" :
        A=PhaseMatrixZeroMean(din,mask)
    else:
        A=PhaseMatrix(din,mask)

    if isinstance(din, CalibrationData):
        pd=din.visphase[mask]
    else:
        pd=UVDataPhase(din,mask)
        
    if pol == "mean":
        p=pd.mean(5).flatten()
    elif type(pol) == int:
        p=pd[:,:,:,:,:,pol].flatten()
    else:
        raise "Don't know that polarisation handling mode"

    # Note that the numarray version of the least-squares is unstable
    # on opterons! Should be phasing out numarray.
    solution=numpy.linalg.lstsq(A,p)
    x= solution[0]
    if pout :
        fout=open(pout, "w")
        pickle.dump([A,p,x], fout)
        
    if ref== "ZeroMean" :
        return  numarray.array( x.tolist() + [ -1 * x.sum() ] )
    else:
        return solution[0]

@AutoFITSOpen    
def AntennaPhaseModel(fin,
                      refids):

    "Model antenna-based phase drifts"

    """
    Returns a model for each antenna phase for each integration time
    in the supplied UV file
    """

    sc=SourcesColumn(fin)
    ts=IntegrationTimeSeries(fin)

    # The observed data:
    aphases, atimes =[] , []
    for refid in refids:
        mask = ( numpy.array(sc) == refid)
        aphases.append( SolveForAntennaPhases( fin[0].data, mask))
        atimes.append( ts[numarray.array(mask)].mean() )

    # Times to interpolate onto
    itimes = numpy.unique(ts)

    aphases= numpy.array(aphases)
    atimes = numpy.array(atimes)

    r = interptools.LinearInterpolate(itimes, atimes, aphases)

    return r

def AntennaPhaseCalibrate(fnamein,
                          fnameout,
                          refids):

    """
    Assumes sources in list refids are calibrators and calibrate aways
    the antenna-based phase  using these sources.
    """

    fin = pyfits.open(fnamein)
    am = AntennaPhaseModel(fin , refids)

    ts=numpy.array(IntegrationTimeSeries(fin))
    itimes = numpy.unique(ts)

    phases=numpy.zeros( len(fin[0].data))
    for i,itime in enumerate(itimes):
        mask= ( ts == itime)
        A=PhaseMatrix(fin[0].data,mask)
        phases [ ts== itime] = -1* numarray.dot(A, numarray.array(am[i]) )

    DePhaseDataAcross(fnamein,
                      fnameout,
                      phases)    


@DefMask
def AntennaPhaseTimeSeries(fnamein,
                           mask=None,
                           **kwargs):

    "Return a matrix showing how phase of each antenna evolves with time"

    din=CalibrationData(fnamein)
    
    ild, ilt= IntegrationList(din.fin, mask)
    td, tt =  IntegrationDayTime(din.fin)

    for  i,t in enumerate(ilt):
        r=SolveForAntennaPhases(din,
                                mask= numpy.array((numpy.array(tt)==t)),
                                **kwargs)
        if  i==0:
            res= numpy.zeros( (len(ilt), len(r)))

        res[i]=r
    return res


#AntennaPhaseModel("temp/base/test12.ms.fits", [1])
