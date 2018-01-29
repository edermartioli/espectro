# -*- coding: iso-8859-1 -*-
"""
    Shebang options:
    #!/usr/bin/python
    #!/opt/anaconda/bin/python
    Created on Mar 29 2017
    
    Description: Measure continuum and save to .s.fits
    
    @author: Eder Martioli
    
    Laboratorio Nacional de Astrofisica, Brazil.
    
    Simple usage example:
    
    python $PATH/App_continuum.py --input=spectrum.m.fits.gz --output=output-spectrum.s.fits --rvsampling=2.0 -tr
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
parser.add_option("-o", "--output", dest="output", help="Output continuum spectrum file",type='string', default="")
parser.add_option("-m", "--rvsampling", dest="rvsampling", help="RV sampling for output spectrum in km/s",type='string', default="2.0")
parser.add_option("-p", action="store_true", dest="polar", help="polar spectrum", default=False)
parser.add_option("-v", action="store_true", dest="verbose", help="verbose", default=False)
parser.add_option("-t", action="store_true", dest="telluric", help="telluric correction", default=False)
parser.add_option("-r", action="store_true", dest="helio", help="heliocentric correction", default=False)

try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print "Error: check usage with App_continuum.py -h "; sys.exit(1);

if options.verbose:
    print 'Input spectrum: ', options.input
    print 'Output continuum spectrum file: ', options.output
    print 'RV sampling for output spectrum in km/s: ', options.rvsampling
    print 'Polar option: ', options.polar
    print 'Telluric correction: ', options.telluric
    print 'Heliocentric correction: ', options.helio

spc = Spectrum(options.input, 'fcal', options.polar, options.telluric, options.helio)

spc.applyMask()

spc.binning(float(options.rvsampling), median=True)

spc.continuum(binsize=200, overlap=100, sigmaclip=3.0, window=3)

spc.saveContinuumToFile(options.output)
