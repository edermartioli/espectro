# -*- coding: iso-8859-1 -*-
"""
    Shebang options:
            #!/usr/bin/python
            #!/opt/anaconda/bin/python
    Created on Mar 29 2017
    
    Description: Normalize spectrum from OPERA fits products and save to .s.fits
    
    @author: Eder Martioli
    
    Laboratorio Nacional de Astrofisica, Brazil.
    
    Simple usage example:
    
    python $PATH/App_normalize.py --input=N20150811G0074.m.fits.gz --linesmask=LineMaskForContinuumDetection.txt --wavemask=OrderValidWavelengthRanges.txt --spectype=raw --rvsampling=2.0 --output=N20150811G0074.s.fits -tr
    """

__version__ = "1.0"

__copyright__ = """
    Copyright (c) ...  All rights reserved.
    """

from optparse import OptionParser
import os,sys
from spectralclass import Spectrum
import espectrolib

import matplotlib.pyplot as plt

parser = OptionParser()
parser.add_option("-i", "--input", dest="input", help="Input spectrum file",type='string', default="")
parser.add_option("-o", "--output", dest="output", help="Output spectrum file",type='string', default="")
parser.add_option("-l", "--linesmask", dest="linesmask", help="Wavelength ranges file containing lines to mask",type='string', default="")
parser.add_option("-w", "--wavemask", dest="wavemask", help="File with each order wavelength range",type='string', default="")
parser.add_option("-s", "--spectype", dest="spectype", help="Spectrum type: raw, norm, or fcal",type='string', default="raw")
parser.add_option("-m", "--rvsampling", dest="rvsampling", help="RV sampling for output spectrum in km/s",type='string', default="10.0")
parser.add_option("-p", action="store_true", dest="plot", help="plot spectrum", default=False)
parser.add_option("-v", action="store_true", dest="verbose", help="verbose", default=False)
parser.add_option("-t", action="store_true", dest="telluric", help="telluric correction", default=False)
parser.add_option("-r", action="store_true", dest="helio", help="heliocentric correction", default=False)


try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print "Error: check usage with App_normalize.py -h "; sys.exit(1);

if options.verbose:
    print 'Input spectrum: ', options.input
    print 'Output spectrum file: ', options.output
    print 'Lines mask: ', options.linesmask
    print 'File with wavelength ranges: ', options.wavemask
    print 'Spectrum type: ', options.spectype
    print 'RV sampling for output spectrum in km/s: ', options.rvsampling
    print 'Plot output spectrum: ', options.plot
    print 'Telluric correction: ', options.telluric
    print 'Heliocentric correction: ', options.helio

if options.verbose :
    print "Loading spectrum"
spc = Spectrum(options.input, options.spectype, False, options.telluric, options.helio, Sort=False)

if options.verbose :
    print "Calculating continuum"
spc.continuumOrderByOrder(binsize = 30, overlap = 10, linesmask=options.linesmask)

if options.verbose :
    print "Normalizing by the continuum"
spc.normalizeByContinuum(constant=10.0)

if options.verbose :
    print "Binning data"
spc.binningOrderByOrder(rvsampling_kps=float(options.rvsampling), median=False)

if options.verbose :
    print "Applying mask"
spc.applyMask(options.wavemask)

if options.verbose :
    print "Sorting data by wavelength"
spc.sortByWavelength()

if options.plot :
    plt.plot(spc.wl, spc.flux)
    plt.show()

if options.output == 'print' or options.output == '' :
    spc.info()
    spc.printdata(printerrors=False)
else :
    if options.verbose :
        print "Saving data to file"
    spc.saveToFile(options.output)




