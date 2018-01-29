# -*- coding: iso-8859-1 -*-
"""
    Shebang options:
            #!/usr/bin/python
            #!/opt/anaconda/bin/python

    Created on Mar 29 2017
    
    Description: Time series of spectral quantities
    
    @author: Eder Martioli
    
    INPE / Laboratorio Nacional de Astrofisica, Brazil.
    
    Simple usage example:
    python $PATH/App_timeseries.py --inputdir=./espectrosflux/ --wlrange="650 665" --spectype=norm --object="AM Her" -tr
    """

__version__ = "1.0"

__copyright__ = """
    Copyright (c) ...  All rights reserved.
    """

from optparse import OptionParser
import os,sys
from spectralclass import Spectrum
from spectralclass import SpectrumChunk
import espectrolib

parser = OptionParser()
parser.add_option("-i", "--inputdir", dest="inputdir", help="Input directory with spectral data",type='string', default="")
parser.add_option("-o", "--object", dest="object", help="Object name",type='string', default="")
parser.add_option("-w", "--wlrange", dest="wlrange", help="Output wavelength range (nm)",type='string', default="")
parser.add_option("-s", "--spectype", dest="spectype", help="Spectrum type: raw, norm, or fcal",type='string', default="raw")
parser.add_option("-p", action="store_true", dest="polar", help="polar spectrum", default=False)
parser.add_option("-v", action="store_true", dest="verbose", help="verbose", default=False)
parser.add_option("-t", action="store_true", dest="telluric", help="telluric correction", default=False)
parser.add_option("-r", action="store_true", dest="helio", help="heliocentric correction", default=False)

try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print "Error: check usage with App_timeseries.py -h "; sys.exit(1);

if options.verbose:
    print 'Input directory: ', options.inputdir
    print 'Object name: ', options.object
    print 'Output wavelength range: ', options.wlrange
    print 'Spectrum type: ', options.spectype
    print 'Polar option: ', options.polar
    print 'Telluric correction: ', options.telluric
    print 'Heliocentric correction: ', options.helio

wl0,wlf = espectrolib.wlrange(options.wlrange, spc)

filelist = espectrolib.generateList(options.inputdir, options.object)

for filepath in filelist :
    spc = Spectrum(filepath, options.spectype, options.polar, options.telluric, options.helio)
    wl,flux,fluxerr = spc.extractChunk(float(wl0), float(wlf))
    chunk = SpectrumChunk(wl,flux,fluxerr)
    chunk.removeBackground(2.0)
    init_guess = [656.270,2500, 0.5]
    wlcen, amp, sigma = chunk.fitgaussian(init_guess)
    print spc.getTimeHJDTT(), wlcen[0], amp[0], sigma[0]

