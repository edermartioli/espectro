# -*- coding: iso-8859-1 -*-
"""
    Shebang options:
            #!/usr/bin/python
            #!/opt/anaconda/bin/python
    Created on Mar 29 2017
    
    Description: Radial velocity time series using a list of isolated spectral lines
    
    @author: Eder Martioli <emartioli@lna.br>
    
    Laboratorio Nacional de Astrofisica, Brazil.
    
    Simple usage example:
    python $PATH/App_rvhotstar.py --inputdir=/Users/edermartioli/Reductions/GRACES/20150807/ --spectype=norm --object="KIC 09472174" --wlrange="667.3 668.3" -tr
    """

__version__ = "1.0"

__copyright__ = """
    Copyright (c) ...  All rights reserved.
    """

from optparse import OptionParser
import os,sys
from spectralclass import Spectrum
from spectralclass import SpectrumChunk
from scipy import constants
import numpy as np
import espectrolib

parser = OptionParser()
parser.add_option("-i", "--inputdir", dest="inputdir", help="Input directory with spectral data",type='string', default="")
parser.add_option("-o", "--object", dest="object", help="Object name",type='string', default="")
parser.add_option("-s", "--spectype", dest="spectype", help="Spectrum type: raw, norm, or fcal",type='string', default="raw")
parser.add_option("-p", action="store_true", dest="polar", help="polar spectrum", default=False)
parser.add_option("-v", action="store_true", dest="verbose", help="verbose", default=False)
parser.add_option("-t", action="store_true", dest="telluric", help="telluric correction", default=False)
parser.add_option("-r", action="store_true", dest="helio", help="heliocentric correction", default=False)

try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print "Error: check usage with App_rvhotstar.py -h "; sys.exit(1);

if options.verbose:
    print 'Input directory: ', options.inputdir
    print 'Object name: ', options.object
    print 'Spectrum type: ', options.spectype
    print 'Polar option: ', options.polar
    print 'Telluric correction: ', options.telluric
    print 'Heliocentric correction: ', options.helio


heliumlines = [447.148, 492.193, 501.568, 587.561, 667.815, 706.518]

delta_wl = 0.5

filelist = espectrolib.generateList(options.inputdir, options.object)

for filepath in filelist :
    spc = Spectrum(filepath, options.spectype, options.polar, options.telluric, options.helio)

    rvvalues = []
    for line in heliumlines :
        wl_line = line
        wl0, wlf = wl_line-delta_wl, wl_line+delta_wl
    
        wl,flux,fluxerr = spc.extractChunk(float(wl0), float(wlf))
        chunk = SpectrumChunk(wl,flux,fluxerr)
        chunk.removeBackground(0.1)
        init_guess = [wl_line,-0.1, 0.3]
        wlcen, amp, sigma = chunk.fitgaussian(init_guess)
        rv = (wlcen[0] - wl_line)*constants.c/wl_line
        rvvalues.append(rv)
    
    median_rv = np.median(rvvalues)
    sig_rv = np.std(rvvalues)
    
    print spc.getTimeHJDTT(), median_rv, sig_rv

