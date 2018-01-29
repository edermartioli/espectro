# -*- coding: iso-8859-1 -*-
"""
    Shebang options:
            #!/usr/bin/python
            #!/opt/anaconda/bin/python
    Created on Mar 29 2017
    
    Description: Estract spectrum from OPERA fits products
    
    @author: Eder Martioli
    
    Laboratorio Nacional de Astrofisica, Brazil.
    
    Simple usage example:
    
    python $PATH/App_extract.py --input=spectrum.m.fits.gz --wlrange="650 665" --spectype=norm -tr
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

import matplotlib.pyplot as plt

parser = OptionParser()
parser.add_option("-i", "--input", dest="input", help="Input spectrum file",type='string', default="")
parser.add_option("-w", "--wlrange", dest="wlrange", help="Output wavelength range (nm)",type='string', default="")
parser.add_option("-s", "--spectype", dest="spectype", help="Spectrum type: raw, norm, or fcal",type='string', default="raw")
parser.add_option("-m", "--wavemask", dest="wavemask", help="File with each order wavelength range",type='string', default="")
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
    print 'Output wavelength range: ', options.wlrange
    print 'Spectrum type: ', options.spectype
    print 'Polar option: ', options.polar
    print 'Telluric correction: ', options.telluric
    print 'Heliocentric correction: ', options.helio
    print 'File with wavelength ranges: ', options.wavemask


spc = Spectrum(options.input, options.spectype, options.polar, options.telluric, options.helio)

if options.wavemask :
    spc.applyMask(options.wavemask)

minorder, maxorder = spc.getMinMaxOders()
for o in range(minorder,maxorder) :
    wl,flux,fluxerr = spc.extractOrder(o)
    print o, wl[0], wl[-1]
    plt.plot(wl,flux)


#wl0,wlf = espectrolib.wlrange(options.wlrange, spc)
#wl,flux,fluxerr = spc.extractChunk(wl0,wlf)
#plt.plot(wl,flux)

#chunk = SpectrumChunk(wl,flux,fluxerr)

#wlsampling=0.0045
#chunk.binning(wlsampling, median=True)
#wl,flux,yerr = chunk.getSpectrum()
#plt.plot(wl,flux)

#plt.errorbar(wl,flux,yerr=fluxerr,marker='o')

plt.show()
