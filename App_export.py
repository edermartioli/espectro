# -*- coding: iso-8859-1 -*-
"""
    Shebang options:
            #!/usr/bin/python
            #!/opt/anaconda/bin/python
    Created on Mar 29 2017
    
    Description: Estract spectrum from OPERA fits products, apply mask, resample and bin data, and save to file: *.s.fits
    
    @author: Eder Martioli
    
    Laboratorio Nacional de Astrofisica, Brazil.
    
    Simple usage example:
    
    python $PATH/App_export.py --input=spectrum.m.fits.gz --output=spectrum.s.fits --rvsampling=2.0 --spectype=fcal -tr
    """

__version__ = "1.0"

__copyright__ = """
    Copyright (c) ...  All rights reserved.
    """

from optparse import OptionParser
import os,sys
from spectralclass import Spectrum
from scipy import constants
import matplotlib.pyplot as plt

parser = OptionParser()
parser.add_option("-i", "--input", dest="input", help="Input spectrum file",type='string', default="")
parser.add_option("-o", "--output", dest="output", help="Output spectrum file",type='string', default="")
parser.add_option("-m", "--rvsampling", dest="rvsampling", help="RV sampling for output spectrum in km/s",type='string', default="")
parser.add_option("-s", "--spectype", dest="spectype", help="Spectrum type: raw, norm, or fcal",type='string', default="raw")
parser.add_option("-w", "--wavemask", dest="wavemask", help="File with each order wavelength range",type='string', default="")
parser.add_option("-p", action="store_true", dest="polar", help="polar spectrum", default=False)
parser.add_option("-v", action="store_true", dest="verbose", help="verbose", default=False)
parser.add_option("-t", action="store_true", dest="telluric", help="telluric correction", default=False)
parser.add_option("-r", action="store_true", dest="helio", help="heliocentric correction", default=False)

try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print "Error: check usage with App_extract.py -h "; sys.exit(1);

if options.verbose:
    print 'Input spectrum: ', options.input
    print 'Output spectrum file: ', options.output
    print 'RV sampling for output spectrum in km/s: ', options.rvsampling
    print 'Spectrum type: ', options.spectype
    print 'Polar option: ', options.polar
    print 'Telluric correction: ', options.telluric
    print 'Heliocentric correction: ', options.helio
    print 'File with wavelength ranges: ', options.wavemask

spc = Spectrum(options.input, options.spectype, options.polar, options.telluric, options.helio, Sort=False)

#spc.binningOrderByOrder(rvsampling_kps=float(options.rvsampling), median=False)

if options.wavemask :
    spc.applyMask(options.wavemask)

spc.sortByWavelength()

if options.output.endswith(".fits") :
    
    wlsampling = (float(options.rvsampling)*1000.0/constants.c)*(spc.wl[-1] + spc.wl[0])/2.0
    spc.resampling(wlsampling)
    spc.saveToFile(options.output, format='fits')

elif options.output.endswith(".txt") :
    
    spc.binning(float(options.rvsampling), median=False)
    spc.saveToFile(options.output, format='ascii')

else :
    spc.printdata(printerrors=False)
