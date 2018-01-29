# -*- coding: iso-8859-1 -*-
"""
    Shebang options:
            #!/usr/bin/python
            #!/opt/anaconda/bin/python

    Created on Mar 29 2017
    
    Description: Stack spectra
    
    @author: Eder Martioli
    
    INPE / Laboratorio Nacional de Astrofisica, Brazil.
    
    Simple usage example:
    python $PATH/App_stack.py --inputdir=./espectrosflux/ --spectype=norm --object="AM Her" -tr
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
import numpy as np

import copy

parser = OptionParser()
parser.add_option("-i", "--inputdir", dest="inputdir", help="Input directory with spectral data",type='string', default="")
parser.add_option("-u", "--output", dest="output", help="Output stacked spectrum",type='string', default="")
parser.add_option("-o", "--object", dest="object", help="Object name",type='string', default="")
parser.add_option("-l", "--linesmask", dest="linesmask", help="Wavelength ranges file containing lines to mask",type='string', default="")
parser.add_option("-w", "--wavemask", dest="wavemask", help="File with each order wavelength range",type='string', default="")
parser.add_option("-s", "--spectype", dest="spectype", help="Spectrum type: raw, norm, or fcal",type='string', default="raw")
parser.add_option("-p", action="store_true", dest="polar", help="polar spectrum", default=False)
parser.add_option("-v", action="store_true", dest="verbose", help="verbose", default=False)
parser.add_option("-t", action="store_true", dest="telluric", help="telluric correction", default=False)
parser.add_option("-r", action="store_true", dest="helio", help="heliocentric correction", default=False)

try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print "Error: check usage with App_stack.py -h "; sys.exit(1);

if options.verbose:
    print 'Input directory: ', options.inputdir
    print 'Output stacked spectrum: ', options.output
    print 'Object name: ', options.object
    print 'Lines mask: ', options.linesmask
    print 'File with wavelength ranges: ', options.wavemask
    print 'Spectrum type: ', options.spectype
    print 'Polar option: ', options.polar
    print 'Telluric correction: ', options.telluric
    print 'Heliocentric correction: ', options.helio

useMedian = True      # method to bin data (if False -> use Mean)
outputRVSampling = 10.0 # in km/s

filelist, rvfilelist = espectrolib.generateListWithRVs(options.inputdir, options.object)

sourceRVs = espectrolib.loadSourceRVs(rvfilelist)
stackspc = Spectrum(filelist[0], "raw", options.polar, options.telluric, options.helio, SourceRV=sourceRVs[0]/1000.0)

stackspc.setAllFluxesToZero()

spc = []

for i in range(len(filelist)) :
    print filelist[i]
    spc.append(Spectrum(filelist[i], options.spectype, options.polar, options.telluric, options.helio, SourceRV=sourceRVs[i]/1000.0))
    spc[i].normalizeOrderByOrder()
    spc[i].resample(stackspc.wl)

for j in range(len(stackspc.flux)) :
    x = []
    for i in range(len(filelist)) :
        x.append((spc[i].flux)[j])

    if useMedian :
        stackspc.flux[j] = np.median(x)
        stackspc.fluxerr[j] = np.std(x)
    else :
        stackspc.flux[j] = np.mean(x)
        stackspc.fluxerr[j] = np.std(x)

stackspc.continuumOrderByOrder(binsize = 30, overlap = 10, linesmask=options.linesmask)
#stackspc_copy1 = copy.deepcopy(stackspc)
stackspc.normalizeByContinuum(constant=10.0)
stackspc.binningOrderByOrder(rvsampling_kps=outputRVSampling, median=True)
stackspc.applyMask(options.wavemask)

#stackspc.binning(rvsampling_kps=outputRVSampling, median=True)

"""
#for o in range(22,55) :
for o in range(45,52) :
    print o

    #continuum = stackspc.extractOrderContinuum(o)
    #plt.plot(orderwl, continuum)
    
    #orderwlc, orderfluxc, orderfluxerrc = stackspc_copy1.extractOrder(o)
    #continuum = stackspc_copy1.extractOrderContinuum(o)
    #plt.plot(orderwlc, orderfluxc)
    #plt.plot(orderwlc, continuum)
    
    orderwl, orderflux, orderfluxerr = stackspc.extractOrder(o)
    plt.plot(orderwl, orderflux, marker='.')

plt.show()
"""

stackspc.sortByWavelength()
stackspc.saveToFile(options.output)
