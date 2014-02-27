# Generated by dada2uvfits.py

FIELDNAME Zenith
N_SCANS   $N_SCANS$
N_INPUTS  64
N_CHANS   $N_CHANS$     # number of channels in spectrum
CORRTYPE  B             # correlation type to use. 'C'(cross), 'B'(both), or 'A'(auto)
INT_TIME  $INT_TIME$    # integration time of scan in seconds
FREQCENT  $FREQCENT$    # observing center freq in MHz
BANDWIDTH $BANDWIDTH$   # total bandwidth in MHz

# To phase to the zenith, these mush be the HA, RA and Dec of the zenith.
HA_HRS    0.000000      # the RA of the desired phase centre (hours)
RA_HRS    $RA_HRS$      # the RA of the desired phase centre (hours)
DEC_DEGS  $DEC_DEGS$    # the DEC of the desired phase centre (degs)
DATE      $DATE$        # YYYYMMDD
TIME      $TIME$        # HHMMSS
INVERT_FREQ 0           # 1 if the freq decreases with channel number
CONJUGATE   1           # conjugate the raw data to fix sign convention problem if necessary

