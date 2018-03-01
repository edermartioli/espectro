# -*- coding: iso-8859-1 -*-
"""
    Shebang options:
            #!/usr/bin/python
            #!/opt/anaconda/bin/python

    Created on Mar 01 2018
    
    Description: Combine a series of cross-correlation plots
    
    @author: Eder Martioli
    
    INPE / Laboratorio Nacional de Astrofisica, Brazil.
    
    Simple usage example:
    python $PATH/App_combinexcorr.py --inputdir=./ --output=outmap
    """

__version__ = "1.0"

__copyright__ = """
    Copyright (c) ...  All rights reserved.
    """

from optparse import OptionParser
import os,sys

import glob
from astropy.io import fits

import matplotlib.pyplot as plt
import numpy as np

parser = OptionParser()
parser.add_option("-i", "--inputdir", dest="inputdir", help="Input directory with x-corr data",type='string', default="")
parser.add_option("-o", "--output", dest="output", help="Output xcorr map",type='string', default="")
parser.add_option("-v", action="store_true", dest="verbose", help="verbose", default=False)

try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print "Error: check usage with App_combinexcorr.py -h "; sys.exit(1);

if options.verbose:
    print 'Input directory: ', options.inputdir
    print 'Output x-corr map: ', options.output

pattern = options.inputdir + "/*m.rv.fits"

files = glob.glob(pattern)

phases = []
xcorrvec = []

period = 4.2308
t0 = 2450203.947
rv = None

for file in files :
    hdu = fits.open(file)
    header = hdu[0].header
    
    bjd = hdu[0].header['BJD']
    rv_measured = hdu[0].header['RV']
    rverr_measured = hdu[0].header['RVERR']
    
    rv = hdu[0].data[0]
    xcorr = hdu[0].data[1]
    xcorr_err = hdu[0].data[2]
    
    xcorrvec.append(xcorr)
    
    phase = ((bjd-t0)/period - float(np.floor((bjd-t0)/period)))
    
    phases.append(phase)
    
    print file, bjd, phase, rv_measured, rverr_measured

#zeros = np.zeros_like(xcorrvec[0])
ones = np.ones_like(xcorrvec[0])

nbins = 100

ph = np.linspace(0.0, 1.0, num=nbins)

dph = 1.0/(2.0*nbins)

imgdata = []
for i in range(len(ph)) :
    
    xcorr_tmp = []
    hasxcorrinrange = False
    for j in range(len(phases)) :
        if ph[i] - dph <= phases[j] < ph[i] + dph :
            xcorr_tmp.append(xcorrvec[j])
            hasxcorrinrange = True

    if hasxcorrinrange :
        mean = np.mean(xcorr_tmp, axis=0)
        imgdata.append(mean)
    else :
        imgdata.append(ones)

if os.path.exists(options.output) :
    os.remove(options.output)

hdu = fits.PrimaryHDU(imgdata)
hdulist = fits.HDUList([hdu])
hdulist.writeto(options.output)


